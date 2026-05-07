/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundRigid;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNBatch;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNS7;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNS8;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNS9;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;

/**
 * Implementation of the K* algorithm to predict protein sequence mutations that improve
 * binding affinity by computing provably accurate Boltzmann-weighted ensembles
 * {@cite Lilien2008 Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson, and Bruce R. Donald, 2005.
 * A Novel Ensemble-Based Scoring and Search Algorithm for Protein Redesign and Its Application
 * to Modify the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine Adenylation Enzyme
 * In Journal of Computational Biology (vol 12. num. 6 pp. 740–761).}.
 */
public class MARKStar {

	public interface ConfEnergyCalculatorFactory {
		ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	public interface ConfSearchFactory {
		ConfSearch make(EnergyMatrix emat, RCs rcs);
	}

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			/**
			 * Value of epsilon in (0,1] for the epsilon-approximation to a partition function.
			 *
			 * Smaller values for epsilon yield more accurate predictions, but can take
			 * longer to run.
			 */
			private double epsilon = 0.683;

			/**
			 * Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence.
			 * Defined in units of kcal/mol.
			 *
			 * More precisely, a sequence is pruned when the following expression is true:
			 *
			 * U(Z_s) < L(W_s) * B(t)
			 *
			 * where:
			 *   - s represents the unbound protein strand, or unbound ligand strand
			 *   - U(Z_s) is the upper bound on the partition function for strand s
			 *   - L(W_s) is the lower bound on the partition function for strand s in the wild type
			 *   - t is the stability threshold
			 *   - B() is the Boltzmann weighting function
			 *
			 * Set to null to disable the filter entirely.
			 */
			private Double stabilityThreshold = 5.0;

			/** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
			private int maxSimultaneousMutations = 1;

			private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();

			/**
			 * If true, prints out information to the console for each minimized conformation during
			 * partition function approximation
			 */
			private boolean showPfuncProgress = false;

			/**
			 * Pattern of the filename to cache energy matrices.
			 *
			 * K*-type algorithms must calculate multiple energy matrices.
			 * By default, these energy matrices are not cached between runs.
			 * To cache energy matrices between runs, supply a pattern such as:
			 *
			 * "theFolder/emat.*.dat"
			 *
			 * The * in the pattern is a wildcard character that will be replaced with
			 * each type of energy matrix used by the K*-type algorithm.
			 */
			private String energyMatrixCachePattern = null;

			private Parallelism parallelism = null;
			private int maxNumConfs = -1;
			private boolean reduceMinimizations = true;
			private boolean useBranchDecomposition = false;
			
			private boolean useGridDP = false;

			public Builder setEpsilon(double val) {
				epsilon = val;
				return this;
			}

			public Builder setStabilityThreshold(Double val) {
				if (val != null && val.isInfinite()) {
					throw new IllegalArgumentException("only finite values allowed. To turn off the filter, pass null");
				}
				stabilityThreshold = val;
				return this;
			}

			public Builder setMaxSimultaneousMutations(int val) {
				maxSimultaneousMutations = val;
				return this;
			}

			public Builder addScoreWriter(KStarScoreWriter val) {
				scoreWriters.add(val);
				return this;
			}

			public Builder addScoreConsoleWriter(KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToConsole(val));
			}

			public Builder addScoreConsoleWriter() {
				return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequenceKStarPfuncs());
			}

			public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
			}

			public Builder addScoreFileWriter(File file) {
				return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
			}

			public Builder setShowPfuncProgress(boolean val) {
				showPfuncProgress = val;
				return this;
			}

			public Builder setEnergyMatrixCachePattern(String val) {
				energyMatrixCachePattern = val;
				return this;
			}

			public Builder setParallelism(Parallelism p) {
				parallelism = p;
				return this;
			}

			public Builder setMaxNumConfs(int maxNumConfs) {
				this.maxNumConfs = maxNumConfs;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, stabilityThreshold, maxSimultaneousMutations, scoreWriters,
						showPfuncProgress, energyMatrixCachePattern, parallelism, maxNumConfs, reduceMinimizations,
						useBranchDecomposition, useGridDP);
			}

			public Builder setReduceMinimizations(boolean reudceMinimizations) {
			    this.reduceMinimizations = reudceMinimizations;
			    return this;
			}

			public Builder setUseBranchDecomposition(boolean val) {
			    this.useBranchDecomposition = val;
			    return this;
			}


			public Builder setUseGridDP(boolean val) {
			    this.useGridDP = val;
			    return this;
			}
		}

		public final double epsilon;
		public final Double stabilityThreshold;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;
		public final boolean showPfuncProgress;
		public final String energyMatrixCachePattern;
		public final Parallelism parallelism;
		public final int maxNumConfs;
		public final boolean reduceMinimizations;
		public final boolean useBranchDecomposition;
		
		public final boolean useGridDP;

		public Settings(double epsilon, Double stabilityThreshold, int maxSimultaneousMutations,
						KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, String energyMatrixCachePattern,
						Parallelism parallelism, int maxNumConfs, boolean reduceMinimizations,
						boolean useBranchDecomposition, boolean useGridDP) {
			this.epsilon = epsilon;
			this.stabilityThreshold = stabilityThreshold;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
			this.showPfuncProgress = dumpPfuncConfs;
			this.energyMatrixCachePattern = energyMatrixCachePattern;
			this.parallelism = parallelism;
			this.maxNumConfs = maxNumConfs;
			this.reduceMinimizations = reduceMinimizations;
			this.useBranchDecomposition = useBranchDecomposition;
			this.useGridDP = useGridDP;
		}

		public String applyEnergyMatrixCachePattern(String type) {

			// the pattern has a * right?
			if (energyMatrixCachePattern.indexOf('*') < 0) {
				throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + energyMatrixCachePattern + "') has no wildcard character (which is *)");
			}

			return energyMatrixCachePattern.replace("*", type);
		}
	}

	public static class ScoredSequence {

		public final Sequence sequence;
		public final KStarScore score;

		public ScoredSequence(Sequence sequence, KStarScore score) {
			this.sequence = sequence;
			this.score = score;
		}

		@Override
		public String toString() {
			return "sequence: " + sequence + "   K*(log10): " + score;
		}

		public String toString(Sequence wildtype) {
			return "sequence: " + sequence.toString(Sequence.Renderer.AssignmentMutations) + "   K*(log10): " + score;
		}
	}

	public static class InitException extends RuntimeException {

		public InitException(ConfSpaceType type, String name) {
			super(String.format("set %s for the %s conf space info before running", name, type.name()));
		}
	}
	public enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public class ConfSpaceInfo {

		public final ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator rigidConfEcalc;
		public final ConfEnergyCalculator minimizingConfEcalc;
		public UpdatingEnergyMatrix correctionEmat;
		public ConfSearchFactory confSearchFactory = null;
		public File confDBFile = null;

		public EnergyMatrix rigidEmat = null;
		public EnergyMatrix minimizingEmat = null;
		public edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator gnnCalc = null;
		public edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator luteCalc = null; // optional, for accuracy comparison
		public MARKStarBoundGNNBatch.GNNStrategy gnnStrategy = null; // null = don't use batch GNN
		// Strategy 6 (GNN_CP_POOL) Conformal Prediction parameters
		// Calibrated from val set: protein q(α=0.001)=0.055, complex q(α=0.001)=0.058
		public double cpAlpha = 0.001;    // per-prediction miscoverage rate
		public double cpDelta = 0.1;      // total pfunc failure probability
		public double cpQ = 0.06;         // CP quantile bound (kcal/mol) — max(protein,complex) rounded up
		public int gnnMiniBatch = -1;     // GNN pool mini batch size (-1 = use default in MARKStarBoundGNNBatch)
		public boolean useStrategy7 = false;  // Strategy 7: decoupled GNN pool
		public int s7GPUBatchSize = 1000;     // GPU batch size for Strategy 7
		public boolean useStrategy8 = false;  // Strategy 8: S7 + subtree GNN for internal nodes
		public boolean useStrategy9 = false;  // Strategy 9: subtree GNN as search router
		public GNNSubtreeEnergyCalculator subtreeGnnCalc = null;
		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator rigidConfEcalc, ConfEnergyCalculator minimizingConfEcalc) {
			this.type = type;
			this.confSpace = confSpace;
			this.rigidConfEcalc = rigidConfEcalc;
			this.minimizingConfEcalc = minimizingConfEcalc;
		}

		private void check() {
			if (rigidConfEcalc == null) {
				throw new InitException(type, "rigidConfEcalc");
			}
			if (minimizingConfEcalc == null) {
				throw new InitException(type, "minimizingConfEcalc");
			}
			if (confSearchFactory == null) {
				throw new InitException(type, "confSearchFactory");
			}
		}

		public void clear() {
			pfuncResults.clear();
		}

		public void calcEmats() {
			SimplerEnergyMatrixCalculator.Builder rigidBuilder = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc);
			if (settings.energyMatrixCachePattern != null) {
				rigidBuilder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase()+".rigid")));
			}
			SimplerEnergyMatrixCalculator.Builder minimizingBuilder = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc);
			if (settings.energyMatrixCachePattern != null) {
				minimizingBuilder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase()+".minimizing")));
			}
			rigidEmat = rigidBuilder.build().calcEnergyMatrix();
			minimizingEmat = minimizingBuilder.build().calcEnergyMatrix();
			correctionEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat);
		}

		public PartitionFunction.Result calcPfunc(int sequenceIndex, BigDecimal stabilityThreshold) {

			Sequence sequence = sequences.get(sequenceIndex);

			// check the cache first
			PartitionFunction.Result result = pfuncResults.get(sequence);
			if (result != null) {
				return result;
			}

			// cache miss, need to compute the partition function

			// make the partition function
			MARKStarBound pfunc;
			if (settings.useBranchDecomposition && useStrategy8 && gnnCalc != null && subtreeGnnCalc != null) {
				// BranchMARK* + Strategy 8 (leaf GNN + subtree GNN)
				BranchMARKStarBound branchS8Pfunc = new BranchMARKStarBound(
						confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
				branchS8Pfunc.setGNNBatchCalculator(gnnCalc);
				branchS8Pfunc.setSubtreeGNN(subtreeGnnCalc);
				branchS8Pfunc.setCPParams(cpAlpha, cpDelta, cpQ);
				branchS8Pfunc.setGPUBatchSize(s7GPUBatchSize);
				pfunc = branchS8Pfunc;
			} else if (useStrategy8 && gnnCalc != null && subtreeGnnCalc != null) {
				MARKStarBoundGNNS8 s8Pfunc = new MARKStarBoundGNNS8(
						confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
				s8Pfunc.setGNNBatchCalculator(gnnCalc);
				s8Pfunc.setSubtreeGNN(subtreeGnnCalc);
				s8Pfunc.setCPParams(cpAlpha, cpDelta, cpQ);
				s8Pfunc.setGPUBatchSize(s7GPUBatchSize);
				pfunc = s8Pfunc;
				} else if (useStrategy9 && gnnCalc != null && subtreeGnnCalc != null) {
					MARKStarBoundGNNS9 s9Pfunc = new MARKStarBoundGNNS9(
							confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
							sequence.makeRCs(confSpace), settings.parallelism);
					s9Pfunc.setGNNBatchCalculator(gnnCalc);
					s9Pfunc.setSubtreeGNN(subtreeGnnCalc);
					s9Pfunc.setCPParams(cpAlpha, cpDelta, cpQ);
					s9Pfunc.setGPUBatchSize(s7GPUBatchSize);
					pfunc = s9Pfunc;
			} else if (settings.useBranchDecomposition && useStrategy7 && gnnCalc != null) {
				// BranchMARK* + GNN (Strategy 7-style decoupled pool)
				BranchMARKStarBound branchPfunc = new BranchMARKStarBound(
						confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
				branchPfunc.setGNNBatchCalculator(gnnCalc);
				branchPfunc.setCPParams(cpAlpha, cpDelta, cpQ);
				branchPfunc.setGPUBatchSize(s7GPUBatchSize);
				pfunc = branchPfunc;
			} else if (useStrategy7 && gnnCalc != null) {
				MARKStarBoundGNNS7 s7Pfunc = new MARKStarBoundGNNS7(
						confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
				s7Pfunc.setGNNBatchCalculator(gnnCalc);
				s7Pfunc.setCPParams(cpAlpha, cpDelta, cpQ);
				s7Pfunc.setGPUBatchSize(s7GPUBatchSize);
				pfunc = s7Pfunc;
			} else if (gnnStrategy != null && gnnCalc != null) {
				MARKStarBoundGNNBatch batchPfunc = new MARKStarBoundGNNBatch(
						confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
				batchPfunc.setGNNStrategy(gnnStrategy);
				batchPfunc.setGNNBatchCalculator(gnnCalc);
				if (gnnStrategy == MARKStarBoundGNNBatch.GNNStrategy.GNN_CP_POOL) {
					batchPfunc.setCPParams(cpAlpha, cpDelta, cpQ);
					if (gnnMiniBatch > 0) {
						batchPfunc.setGNNMiniBatch(gnnMiniBatch);
					}
				}
				if (luteCalc != null) {
					batchPfunc.setLUTECalculator(luteCalc);
				}
				pfunc = batchPfunc;
			} else if (settings.useBranchDecomposition) {
				pfunc = new BranchMARKStarBound(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
			} else {
				pfunc = new MARKStarBoundFastQueues(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
						sequence.makeRCs(confSpace), settings.parallelism);
			}
			confSearchFactory = (emat, rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
			RCs rcs = sequence.makeRCs(confSpace);
			//GradientDescentMARKStarPfunc pfunc = new GradientDescentMARKStarPfunc(confSpace, rigidEmat, minimizingEmat,
			//		rcs, minimizingConfEcalc);
			pfunc.reduceMinimizations = settings.reduceMinimizations;
			pfunc.stateName = type.name();
			if(settings.maxNumConfs > 0)
				pfunc.setMaxNumConfs(settings.maxNumConfs);
			pfunc.setReportProgress(settings.showPfuncProgress);

			pfunc.setCorrections(correctionEmat);

			// GNN energy surrogate (optional) — only for non-batch strategies
			// Strategy 6 uses gnnBatchCalc for CP bounds; gnnCalc single-conf predictions
			// are inaccurate and corrupt bounds in processFullConfNode
			if (gnnCalc != null && gnnStrategy == null) {
				pfunc.setGNNCalculator(gnnCalc);
			}

			// Phase 7: Grid DP upper bound
			if (settings.useGridDP) {
				pfunc.setUseGridDP(true);
			}

			// compute it
			pfunc.init(settings.epsilon);
			Stopwatch computeTimer = new Stopwatch().start();
			pfunc.compute();
			computeTimer.stop();
			// System.out.println("Computation for "+sequence.toString("+":"+computeTimer.getTime(2));

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);
			return result;
		}
	}

	private interface Scorer {
		KStarScore score(int sequenceNumber, PartitionFunction.Result proteinResult, PartitionFunction.Result ligandResult, PartitionFunction.Result complexResult);
	}

	/** A configuration space containing just the protein strand */
	public final ConfSpaceInfo protein;

	/** A configuration space containing just the ligand strand */
	public final ConfSpaceInfo ligand;

	/** A configuration space containing both the protein and ligand strands */
	public final ConfSpaceInfo complex;

	/** Calculates the rigid energy for a molecule */
	public final EnergyCalculator rigidEcalc;

	/** Calculates the minimized energy for a molecule */
	public final EnergyCalculator minimizingEcalc;

	/** A function that makes a ConfEnergyCalculator with the desired options */
	public final ConfEnergyCalculatorFactory confEcalcFactory;

	/** Optional and overridable settings for K* */
	public final Settings settings;

	private List<Sequence> sequences;

	public MARKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex,
					EnergyCalculator rigidEcalc, EnergyCalculator minimizingEcalc,
					ConfEnergyCalculatorFactory confEcalcFactory, Settings settings) {
		this.protein = new ConfSpaceInfo(ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, rigidEcalc), confEcalcFactory.make(protein, minimizingEcalc));
		this.ligand = new ConfSpaceInfo(ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, rigidEcalc), confEcalcFactory.make(ligand, minimizingEcalc));
		this.complex = new ConfSpaceInfo(ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, rigidEcalc), confEcalcFactory.make(complex, minimizingEcalc));
		this.rigidEcalc = rigidEcalc;
		this.minimizingEcalc = minimizingEcalc;
		this.confEcalcFactory = confEcalcFactory;
		this.settings = settings;
		this.sequences = new ArrayList();
	}

	public void precalcEmats() {
		// compute energy matrices
		protein.calcEmats();
		ligand.calcEmats();
		complex.calcEmats();
	}

	public List<ScoredSequence> run() {

		List<ScoredSequence> scores = new ArrayList<>();

		// compute energy matrices
		protein.calcEmats();
		ligand.calcEmats();
		complex.calcEmats();


		// collect all the seque// collect all the sequences explicitly
		if (complex.confSpace.seqSpace.containsWildTypeSequence()) {
			sequences.add(complex.confSpace.seqSpace.makeWildTypeSequence());
		}
		sequences.addAll(complex.confSpace.seqSpace.getMutants(settings.maxSimultaneousMutations, true));

		// TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

		// now we know how many sequences there are in total
		int n = sequences.size();

		// make the sequence scorer and reporter
		Scorer scorer = (sequenceNumber, proteinResult, ligandResult, complexResult) -> {

			// compute the K* score
			KStarScore kstarScore = new KStarScore(proteinResult, ligandResult, complexResult);
			Sequence sequence = sequences.get(sequenceNumber);
			scores.add(new ScoredSequence(sequence, kstarScore));

			// report scores
			settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
				sequenceNumber,
				n,
				sequence,
				kstarScore
			));

			return kstarScore;
		};

		System.out.println("computing K* scores for " + sequences.size() + " sequences to epsilon = " + settings.epsilon + " ...");
		settings.scoreWriters.writeHeader();
		// TODO: progress bar?

		// compute wild type partition functions first (always at pos 0)
		PartitionFunction.Result wtProtein = protein.calcPfunc(0, BigDecimal.ZERO);
		PartitionFunction.Result wtLigand = ligand.calcPfunc(0, BigDecimal.ZERO);
		PartitionFunction.Result wtComplex = complex.calcPfunc(0, BigDecimal.ZERO);
		System.out.println("[PFUNC seq 0 WT] protein: qstar=" + String.format("%.6e", wtProtein.values.qstar)
			+ " qprime=" + String.format("%.6e", wtProtein.values.qprime)
			+ " eps=" + String.format("%.6f", wtProtein.values.getEffectiveEpsilon()));
		System.out.println("[PFUNC seq 0 WT] ligand:  qstar=" + String.format("%.6e", wtLigand.values.qstar)
			+ " qprime=" + String.format("%.6e", wtLigand.values.qprime)
			+ " eps=" + String.format("%.6f", wtLigand.values.getEffectiveEpsilon()));
		System.out.println("[PFUNC seq 0 WT] complex: qstar=" + String.format("%.6e", wtComplex.values.qstar)
			+ " qprime=" + String.format("%.6e", wtComplex.values.qprime)
			+ " eps=" + String.format("%.6f", wtComplex.values.getEffectiveEpsilon()));
		KStarScore wildTypeScore = scorer.score(0, wtProtein, wtLigand, wtComplex);
		BigDecimal proteinStabilityThreshold = null;
		BigDecimal ligandStabilityThreshold = null;
		if (settings.stabilityThreshold != null) {
			BigDecimal stabilityThresholdFactor = new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(settings.stabilityThreshold);
			proteinStabilityThreshold = wildTypeScore.protein.values.calcLowerBound().multiply(stabilityThresholdFactor);
			ligandStabilityThreshold = wildTypeScore.ligand.values.calcLowerBound().multiply(stabilityThresholdFactor);
		}

		// Sequence index filter: only compute specified sequences (for debugging)
		java.util.Set<Integer> seqFilter = null;
		String seqFilterProp = System.getProperty("osprey.seqFilter");
		if (seqFilterProp != null && !seqFilterProp.isEmpty()) {
			seqFilter = new java.util.HashSet<>();
			for (String s : seqFilterProp.split(",")) {
				seqFilter.add(Integer.parseInt(s.trim()));
			}
			System.out.println("[SEQ_FILTER] Only computing sequences: " + seqFilter);
		}

		// compute all the partition functions and K* scores for the rest of the sequences
		for (int i=1; i<n; i++) {

			// Skip sequences not in filter
			if (seqFilter != null && !seqFilter.contains(i)) {
				continue;
			}

			// get the pfuncs, with short circuits as needed
			long pfuncT0 = System.currentTimeMillis();
			final PartitionFunction.Result proteinResult = protein.calcPfunc(i, proteinStabilityThreshold);
			long proteinMs = System.currentTimeMillis() - pfuncT0;
			final PartitionFunction.Result ligandResult;
			final PartitionFunction.Result complexResult;
			if (!KStarScore.isLigandComplexUseful(proteinResult)) {
				ligandResult = PartitionFunction.Result.makeAborted();
				complexResult = PartitionFunction.Result.makeAborted();
			} else {
				long ligT0 = System.currentTimeMillis();
				ligandResult = ligand.calcPfunc(i, ligandStabilityThreshold);
				long ligandMs = System.currentTimeMillis() - ligT0;
				if (!KStarScore.isComplexUseful(proteinResult, ligandResult)) {
					complexResult = PartitionFunction.Result.makeAborted();
				} else {
					long cplxT0 = System.currentTimeMillis();
					complexResult = complex.calcPfunc(i, BigDecimal.ZERO);
					long complexMs = System.currentTimeMillis() - cplxT0;
					System.out.println("[PFUNC seq " + i + "] complex: time=" + complexMs + "ms"
						+ " qstar=" + String.format("%.6e", complexResult.values.qstar)
						+ " qprime=" + String.format("%.6e", complexResult.values.qprime)
						+ " eps=" + String.format("%.6f", complexResult.values.getEffectiveEpsilon()));
				}
				System.out.println("[PFUNC seq " + i + "] ligand: time=" + ligandMs + "ms"
					+ " qstar=" + String.format("%.6e", ligandResult.values.qstar)
					+ " eps=" + String.format("%.6f", ligandResult.values.getEffectiveEpsilon()));
			}
			System.out.println("[PFUNC seq " + i + "] protein: time=" + proteinMs + "ms"
				+ " qstar=" + String.format("%.6e", proteinResult.values.qstar)
				+ " qprime=" + String.format("%.6e", proteinResult.values.qprime)
				+ " eps=" + String.format("%.6f", proteinResult.values.getEffectiveEpsilon()));

			// DEBUG: print pfunc details
			System.out.println("[DEBUG seq " + i + "] seq=" + sequences.get(i).toString(Sequence.Renderer.ResType)
				+ " protein: status=" + proteinResult.status + " qstar=" + proteinResult.values.qstar
				+ " | ligand: status=" + ligandResult.status + " qstar=" + ligandResult.values.qstar
				+ " | complex: status=" + complexResult.status + " qstar=" + complexResult.values.qstar);

			scorer.score(i, proteinResult, ligandResult, complexResult);
		}

		return scores;
	}
	public Iterable<ConfSpaceInfo> confSpaceInfos() {
		return Arrays.asList(protein, ligand, complex);
	}
}
