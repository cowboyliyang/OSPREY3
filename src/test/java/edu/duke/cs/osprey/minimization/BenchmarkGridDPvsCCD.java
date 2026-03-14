package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.ConfSpaceSpecificSurrogateFactory;
import edu.duke.cs.osprey.energy.approximation.MLPSurrogateMatrix;
import edu.duke.cs.osprey.energy.approximation.TaskGlobalMLPSurrogate;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.lute.BranchResidualConfEnergyCalculator;
import edu.duke.cs.osprey.markstar.framework.ApproximatedGridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.GridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.MLPApproximatedGridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.branch.BranchDecomposition;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeNode;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.ematrix.NewEPICMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Test;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.SimpleBounds;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.*;

/**
 * Benchmark: Grid DP Minimizer vs CCD Minimizer.
 *
 * Compares speed and energy quality of the Grid DP approach
 * (discrete grid + branch decomposition DP) against standard CCD
 * (cyclic coordinate descent) on the same set of conformations.
 */
public class BenchmarkGridDPvsCCD {

    private static final double DIST_CUTOFF = 8.0;
    private static final double ENERGY_CUTOFF = 0.1;
    private static final String SURROGATE_PROP_PREFIX = "osprey.bench.surrogate.";

    @Test
    public void benchmark1CC8() {
        long wallStart = System.nanoTime();

        // 1. Build conf space: 1CC8, 3 flexible residues
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("Conf space: " + confSpace.positions.size() + " positions");
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }

        Parallelism parallelism = Parallelism.makeCpu(1);  // single thread for fair timing

        // 2. Compute energy matrices (rigid + minimizing)
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build()
                .calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build()
                .calcEnergyMatrix();
        System.out.println("Energy matrices computed.");

        // 3. Build interaction graph with dual cutoff
        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + ig.getNumPositions() + " positions, "
                + countEdges(ig) + " edges (vs " + ig.getNumPositions() * (ig.getNumPositions() - 1) / 2 + " complete)");

        // 4. Compute branch decomposition and root the tree
        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        // 5. Get N conformations via A*
        int numConfs = 15;
        System.out.println("Finding " + numConfs + " conformations via A*...");
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional()
                .build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        System.out.println("Got " + scoredConfs.size() + " conformations.\n");

        SurrogateRunConfig surrogateCfg = SurrogateRunConfig.fromSystemProperties();
        ConfEnergyCalculator surrogateConfEcalc = null;
        MLPSurrogateMatrix mlpSurrogate = null;
        TaskGlobalMLPSurrogate globalMLPSurrogate = null;
        BranchResidualConfEnergyCalculator branchResidualConfEcalc = null;
        if (surrogateCfg.enabled) {
            if (surrogateCfg.modelType == SurrogateModelType.Quadratic) {
                surrogateConfEcalc = ConfSpaceSpecificSurrogateFactory.withTaskSpecificApproximator(
                        confEcalcMinimized,
                        surrogateCfg.cacheRoot,
                        surrogateCfg.taskTag,
                        surrogateCfg.samplesPerParam,
                        surrogateCfg.errorBudget
                );
            } else if (surrogateCfg.modelType == SurrogateModelType.MLP) {
                if (surrogateCfg.mlpSingleTaskModel) {
                    globalMLPSurrogate = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificGlobalMLPSurrogate(
                            confEcalcMinimized,
                            surrogateCfg.cacheRoot,
                            surrogateCfg.taskTag,
                            surrogateCfg.mlpGlobalNumTuples,
                            surrogateCfg.mlpGlobalSamplesPerTuple,
                            surrogateCfg.mlpHidden1,
                            surrogateCfg.mlpHidden2,
                            surrogateCfg.mlpEpochs,
                            surrogateCfg.mlpBatchSize,
                            surrogateCfg.mlpLearningRate
                    );
                } else {
                    mlpSurrogate = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificMLPSurrogate(
                            confEcalcMinimized,
                            surrogateCfg.cacheRoot,
                            surrogateCfg.taskTag,
                            surrogateCfg.samplesPerParam,
                            surrogateCfg.mlpMinSamplesPerModel,
                            surrogateCfg.mlpMaxSamplesPerModel,
                            surrogateCfg.mlpHidden1,
                            surrogateCfg.mlpHidden2,
                            surrogateCfg.mlpEpochs,
                            surrogateCfg.mlpBatchSize,
                            surrogateCfg.mlpLearningRate
                    );
                }
            } else if (surrogateCfg.modelType == SurrogateModelType.BranchResidual) {
                ConfSpaceSpecificSurrogateFactory.BranchResidualRidgeConfig branchCfg =
                        new ConfSpaceSpecificSurrogateFactory.BranchResidualRidgeConfig();
                branchCfg.cacheRoot = surrogateCfg.cacheRoot;
                branchCfg.taskTag = surrogateCfg.taskTag;
                branchCfg.forceRetrain = surrogateCfg.forceRetrain;
                branchCfg.featureHashDim = surrogateCfg.branchResidualFeatureHashDim;
                branchCfg.samplerMode = surrogateCfg.branchResidualSamplerMode;
                branchCfg.pretrainGridDPSamples = surrogateCfg.branchResidualPretrainGridDPSamples;
                branchCfg.finetuneCCDSamples = surrogateCfg.branchResidualFineTuneCCDSamples;
                branchCfg.lambdaPre = surrogateCfg.branchResidualLambdaPre;
                branchCfg.lambdaFine = surrogateCfg.branchResidualLambdaFine;
                branchCfg.pretrainUseForcefieldFinalEval = surrogateCfg.branchResidualPretrainUseFinalEval;
                branchCfg.gridSize = surrogateCfg.branchResidualGridSize;
                branchCfg.gridThreads = surrogateCfg.branchResidualGridThreads;
                branchCfg.lowEnergyPoolLimit = surrogateCfg.branchResidualLowEnergyPoolLimit;
                branchCfg.randomSeed = surrogateCfg.branchResidualRandomSeed;
                branchCfg.ridgeMaxIterations = surrogateCfg.branchResidualRidgeMaxIter;
                branchCfg.ridgeTolerance = surrogateCfg.branchResidualRidgeTolerance;
                branchResidualConfEcalc = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificBranchResidualRidge(
                        confEcalcMinimized,
                        ematMinimized,
                        ig,
                        rootEdge,
                        ffparams,
                        eref,
                        branchCfg
                );
            }
            System.out.println(String.format(
                    "Surrogate enabled: type=%s, taskTag=%s, cache=%s, samplesPerParam=%d, mlpMode=%s, mlpMinSamples=%d, mlpMaxSamples=%d, mlpGlobalTuples=%d, mlpGlobalSamplesPerTuple=%d, errorBudget=%g, requireFull=%s, fallback=%s, mlpHidden=%dx%d, mlpEpochs=%d, branchEnabled=%s, branchSampler=%s, branchPre=%d, branchFine=%d, branchHash=%d",
                    surrogateCfg.modelType,
                    surrogateCfg.taskTag,
                    surrogateCfg.cacheRoot.getAbsolutePath(),
                    surrogateCfg.samplesPerParam,
                    surrogateCfg.mlpSingleTaskModel ? "single-task-global" : "decomposed",
                    surrogateCfg.mlpMinSamplesPerModel,
                    surrogateCfg.mlpMaxSamplesPerModel,
                    surrogateCfg.mlpGlobalNumTuples,
                    surrogateCfg.mlpGlobalSamplesPerTuple,
                    surrogateCfg.errorBudget,
                    surrogateCfg.requireFullApproximation,
                    surrogateCfg.gridDPFallbackToForcefield,
                    surrogateCfg.mlpHidden1,
                    surrogateCfg.mlpHidden2,
                    surrogateCfg.mlpEpochs,
                    surrogateCfg.branchResidualEnabled,
                    surrogateCfg.branchResidualSamplerMode,
                    surrogateCfg.branchResidualPretrainGridDPSamples,
                    surrogateCfg.branchResidualFineTuneCCDSamples,
                    surrogateCfg.branchResidualFeatureHashDim
            ));
        }

        // 6. Run benchmark — Hybrid GridDP and Multi Warm Start only
        int n = scoredConfs.size();
        int numCores = Runtime.getRuntime().availableProcessors();
        double[] rigidEnergies = new double[n];
        double[] ccdEnergies = new double[n];
        double[] gridDPEnergies = new double[n];
        double[] hybridEnergies = new double[n];
        double[] multiWarmEnergies = new double[n];
        double[] gridDPForcefieldEnergies = new double[n];
        double[] hybridForcefieldEnergies = new double[n];
        double[] multiWarmForcefieldEnergies = new double[n];
        for (int i = 0; i < n; i++) {
            gridDPForcefieldEnergies[i] = Double.NaN;
            hybridForcefieldEnergies[i] = Double.NaN;
            multiWarmForcefieldEnergies[i] = Double.NaN;
        }
        boolean compareMLPtoForcefield = surrogateCfg.enabled && surrogateCfg.modelType == SurrogateModelType.MLP;
        boolean compareBranchResidualToForcefield = surrogateCfg.enabled
                && surrogateCfg.modelType == SurrogateModelType.BranchResidual
                && branchResidualConfEcalc != null;
        double[] branchResidualPredEnergies = new double[n];
        Arrays.fill(branchResidualPredEnergies, Double.NaN);

        // Compute rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;

        // Warm up JIT
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            if (surrogateCfg.enabled && surrogateCfg.modelType == SurrogateModelType.Quadratic) {
                surrogateConfEcalc.calcEnergy(new RCTuple(warmConf));
            } else if (surrogateCfg.enabled && surrogateCfg.modelType == SurrogateModelType.MLP) {
                RCTuple tuple = new RCTuple(warmConf);
                var obj = makeObjectiveFunction(
                        confSpace, tuple, confEcalcMinimized, surrogateConfEcalc, mlpSurrogate,
                        surrogateCfg, ecalcContext, globalMLPSurrogate
                );
                obj.objective.getValue(obj.objective.getDOFsCenter());
                obj.close();
            } else if (compareBranchResidualToForcefield) {
                branchResidualConfEcalc.calcEnergy(warmConf);
            } else {
                confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            }
            makeGridDPMinimizer(confSpace, ig, rootEdge, ffparams, eref, surrogateCfg, surrogateConfEcalc, mlpSurrogate).minimize(warmConf);
        }

        ConfEnergyCalculator ccdReferenceEcalc = (surrogateCfg.enabled && surrogateCfg.modelType == SurrogateModelType.Quadratic)
                ? surrogateConfEcalc : confEcalcMinimized;

        // CCD reference energies (original SurfingLineSearcher, maxIterations=30)
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            RCTuple tuple = new RCTuple(scoredConfs.get(i).getAssignments());
            if (surrogateCfg.enabled && surrogateCfg.modelType == SurrogateModelType.MLP) {
                ObjectiveWithCleanup obj = makeObjectiveFunction(
                        confSpace, tuple, confEcalcMinimized, surrogateConfEcalc, mlpSurrogate,
                        surrogateCfg, ecalcContext, globalMLPSurrogate
                );
                SimpleCCDMinimizer ccd = new SimpleCCDMinimizer();
                ccd.init(obj.objective);
                ccd.setMaxIterations(SimpleCCDMinimizer.DefaultMaxIterations);
                ccdEnergies[i] = ccd.minimizeFromCenter().energy;
                ccd.close();
                obj.close();
            } else {
                ccdEnergies[i] = ccdReferenceEcalc.calcEnergy(tuple).energy;
            }
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        if (compareBranchResidualToForcefield) {
            for (int i = 0; i < n; i++) {
                int[] conf = scoredConfs.get(i).getAssignments();
                branchResidualPredEnergies[i] = branchResidualConfEcalc.calcEnergy(conf);
            }
        }

        // ===== Method 1: Hybrid GridDP + CCD =====
        long hybridTotalStart = System.nanoTime();
        GridDPMinimizer gridDPForHybrid = makeGridDPMinimizer(
                confSpace, ig, rootEdge, ffparams, eref, surrogateCfg, surrogateConfEcalc, mlpSurrogate
        );
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);
            GridDPMinimizer.Result gridResult = gridDPForHybrid.minimize(conf);
            gridDPEnergies[i] = gridResult.energy;
            if (compareMLPtoForcefield && gridResult.bestDOFValues != null) {
                gridDPForcefieldEnergies[i] = evaluateForcefieldEnergyAt(
                        confSpace, confEcalcMinimized, ecalcContext, tuple, gridResult.bestDOFValues
                );
            }
            ObjectiveWithCleanup obj = makeObjectiveFunction(
                    confSpace, tuple, confEcalcMinimized, surrogateConfEcalc, mlpSurrogate,
                    surrogateCfg, ecalcContext, globalMLPSurrogate
            );
            edu.duke.cs.osprey.minimization.ObjectiveFunction objFunc = obj.objective;
            SimpleCCDMinimizer ccd = new SimpleCCDMinimizer((context) -> new QuadraticLineSearcher());
            ccd.init(objFunc);
            ccd.setMaxIterations(2);
            double[] gridDOFs = gridResult.bestDOFValues;
            cern.colt.matrix.DoubleMatrix1D startx = cern.colt.matrix.DoubleFactory1D.dense.make(objFunc.getNumDOFs());
            if (gridDOFs != null) {
                if (gridDOFs.length != objFunc.getNumDOFs()) {
                    throw new IllegalStateException(String.format(
                            "GridDP DOF count (%d) does not match CCD objective DOFs (%d) for conf %s",
                            gridDOFs.length, objFunc.getNumDOFs(), tuple.stringListing()
                    ));
                }
                for (int d = 0; d < objFunc.getNumDOFs(); d++) {
                    startx.set(d, gridDOFs[d]);
                }
            } else {
                startx.assign(objFunc.getDOFsCenter());
            }
            Minimizer.Result ccdResult = ccd.minimizeFrom(startx);
            hybridEnergies[i] = ccdResult.energy;
            if (compareMLPtoForcefield) {
                hybridForcefieldEnergies[i] = evaluateForcefieldEnergyAt(
                        confSpace, confEcalcMinimized, ecalcContext, tuple, ccdResult.dofValues
                );
            }
            ccd.close();
            obj.close();
        }
        long hybridTotalNs = System.nanoTime() - hybridTotalStart;

        // ===== Method 2: Multi Warm Start Hybrid with phase timing =====
        int multiN = 20;
        System.out.println("Running Multi Warm Start Hybrid with top-" + multiN + " starts, " + numCores + " threads...");
        ExecutorService multiExecutor = Executors.newFixedThreadPool(numCores);
        long multiWarmTotalStart = System.nanoTime();
        long totalMWGridDPNs = 0, totalMWInitNs = 0, totalMWPopNs = 0, totalMWCCDCreateNs = 0, totalMWCCDMinNs = 0;
        GridDPMinimizer gridDPForMulti = makeGridDPMinimizer(
                confSpace, ig, rootEdge, ffparams, eref, surrogateCfg, surrogateConfEcalc, mlpSurrogate
        );
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            // Phase 1: GridDP minimize
            long t0 = System.nanoTime();
            gridDPForMulti.minimize(conf);
            totalMWGridDPNs += System.nanoTime() - t0;

            // Phase 2: Init enumeration (pre-sort lambdas, seed PQ)
            t0 = System.nanoTime();
            gridDPForMulti.initTopNEnumeration();
            totalMWInitNs += System.nanoTime() - t0;

            // Phase 3: Pop top-N starting points
            t0 = System.nanoTime();
            List<double[]> topNDOFs = new ArrayList<>();
            for (int k = 0; k < multiN; k++) {
                double[] dofs = gridDPForMulti.popNextDOFValues();
                if (dofs == null) break;
                topNDOFs.add(dofs);
            }
            totalMWPopNs += System.nanoTime() - t0;

            if (i == 0) {
                System.out.println("  Conf 0: got " + topNDOFs.size() + " starting points (requested " + multiN + ")");
            }

            // Phase 4: Parallel CCD — pre-create objects in main thread to avoid contention
            t0 = System.nanoTime();
            int nStarts = topNDOFs.size();
            edu.duke.cs.osprey.minimization.ObjectiveFunction[] objFuncs = new edu.duke.cs.osprey.minimization.ObjectiveFunction[nStarts];
            ObjectiveWithCleanup[] objHolders = new ObjectiveWithCleanup[nStarts];
            SimpleCCDMinimizer[] ccds = new SimpleCCDMinimizer[nStarts];
            cern.colt.matrix.DoubleMatrix1D[] startxs = new cern.colt.matrix.DoubleMatrix1D[nStarts];
            for (int k = 0; k < nStarts; k++) {
                objHolders[k] = makeObjectiveFunction(
                        confSpace, tuple, confEcalcMinimized, surrogateConfEcalc, mlpSurrogate,
                        surrogateCfg, ecalcContext, globalMLPSurrogate
                );
                objFuncs[k] = objHolders[k].objective;
                ccds[k] = new SimpleCCDMinimizer((context) -> new QuadraticLineSearcher());
                ccds[k].init(objFuncs[k]);
                ccds[k].setMaxIterations(2);
                startxs[k] = cern.colt.matrix.DoubleFactory1D.dense.make(objFuncs[k].getNumDOFs());
                double[] startDOFs = topNDOFs.get(k);
                if (startDOFs.length != objFuncs[k].getNumDOFs()) {
                    throw new IllegalStateException(String.format(
                            "Warm-start DOF count (%d) does not match CCD objective DOFs (%d) for conf %s",
                            startDOFs.length, objFuncs[k].getNumDOFs(), tuple.stringListing()
                    ));
                }
                for (int d = 0; d < objFuncs[k].getNumDOFs(); d++) {
                    startxs[k].set(d, startDOFs[d]);
                }
            }
            long objCreateNs = System.nanoTime() - t0;
            totalMWCCDCreateNs += objCreateNs;

            // Submit only minimizeFrom() to thread pool
            long t1 = System.nanoTime();
            List<Future<Minimizer.Result>> ccdFutures = new ArrayList<>();
            for (int k = 0; k < nStarts; k++) {
                final int kk = k;
                ccdFutures.add(multiExecutor.submit(() -> {
                    return ccds[kk].minimizeFrom(startxs[kk]);
                }));
            }
            double bestE = Double.POSITIVE_INFINITY;
            Minimizer.Result bestResult = null;
            for (Future<Minimizer.Result> f : ccdFutures) {
                try {
                    Minimizer.Result r = f.get();
                    double e = r.energy;
                    if (e < bestE) {
                        bestE = e;
                        bestResult = r;
                    }
                } catch (Exception ex) {
                    throw new RuntimeException("Multi warm start CCD failed", ex);
                }
            }
            multiWarmEnergies[i] = bestE;
            if (compareMLPtoForcefield && bestResult != null) {
                multiWarmForcefieldEnergies[i] = evaluateForcefieldEnergyAt(
                        confSpace, confEcalcMinimized, ecalcContext, tuple, bestResult.dofValues
                );
            }
            totalMWCCDMinNs += System.nanoTime() - t1;
            // Cleanup
            for (int k = 0; k < nStarts; k++) {
                objHolders[k].close();
                ccds[k].close();
            }
        }
        multiExecutor.shutdown();
        long multiWarmTotalNs = System.nanoTime() - multiWarmTotalStart;

        // 7. Print comparison table
        if (compareMLPtoForcefield) {
            System.out.println(String.format("%-5s| %11s | %11s | %11s | %11s | %8s | %8s | %8s | %8s | %8s | %8s | %8s",
                    "Conf", "CCD_E", "GridDP_E", "Hybrid_E", "MWarm_E", "CCD_drop", "GDP_rec%", "Hyb_rec%", "MW_rec%",
                    "GDP_dFF", "Hyb_dFF", "MW_dFF"));
        } else if (compareBranchResidualToForcefield) {
            System.out.println(String.format("%-5s| %11s | %11s | %11s | %11s | %11s | %8s | %8s | %8s | %8s",
                    "Conf", "CCD_E", "GridDP_E", "Hybrid_E", "MWarm_E", "Branch_E",
                    "CCD_drop", "GDP_rec%", "Hyb_rec%", "MW_rec%"));
        } else {
            System.out.println(String.format("%-5s| %11s | %11s | %11s | %11s | %8s | %8s | %8s | %8s",
                    "Conf", "CCD_E", "GridDP_E", "Hybrid_E", "MWarm_E", "CCD_drop", "GDP_rec%", "Hyb_rec%", "MW_rec%"));
        }
        System.out.println(new String(new char[105]).replace('\0', '-'));

        double sumGapGridDP = 0, sumGapHybrid = 0, sumGapMultiWarm = 0;
        double sumRecGridDP = 0, sumRecHybrid = 0, sumRecMultiWarm = 0;
        double sumAbsGDPdFF = 0, sumAbsHybdFF = 0, sumAbsMWdFF = 0;
        double sumAbsBranchdFF = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double gdpDrop = rigidEnergies[i] - gridDPEnergies[i];
            double hybDrop = rigidEnergies[i] - hybridEnergies[i];
            double mwDrop = rigidEnergies[i] - multiWarmEnergies[i];
            double gdpRec = (ccdDrop > 0) ? (gdpDrop / ccdDrop * 100.0) : 0;
            double hybRec = (ccdDrop > 0) ? (hybDrop / ccdDrop * 100.0) : 0;
            double mwRec = (ccdDrop > 0) ? (mwDrop / ccdDrop * 100.0) : 0;
            sumGapGridDP += gridDPEnergies[i] - ccdEnergies[i];
            sumGapHybrid += hybridEnergies[i] - ccdEnergies[i];
            sumGapMultiWarm += multiWarmEnergies[i] - ccdEnergies[i];
            sumRecGridDP += gdpRec; sumRecHybrid += hybRec; sumRecMultiWarm += mwRec;

            if (compareMLPtoForcefield) {
                double gdpdFF = gridDPEnergies[i] - gridDPForcefieldEnergies[i];
                double hybdFF = hybridEnergies[i] - hybridForcefieldEnergies[i];
                double mwdFF = multiWarmEnergies[i] - multiWarmForcefieldEnergies[i];
                sumAbsGDPdFF += Math.abs(gdpdFF);
                sumAbsHybdFF += Math.abs(hybdFF);
                sumAbsMWdFF += Math.abs(mwdFF);
                System.out.println(String.format(
                        "%-5d| %11.4f | %11.4f | %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%% | %7.1f%% | %8.4f | %8.4f | %8.4f",
                        i, ccdEnergies[i], gridDPEnergies[i], hybridEnergies[i], multiWarmEnergies[i],
                        ccdDrop, gdpRec, hybRec, mwRec, gdpdFF, hybdFF, mwdFF
                ));
            } else if (compareBranchResidualToForcefield) {
                double branchE = branchResidualPredEnergies[i];
                sumAbsBranchdFF += Math.abs(branchE - ccdEnergies[i]);
                System.out.println(String.format(
                        "%-5d| %11.4f | %11.4f | %11.4f | %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%% | %7.1f%%",
                        i, ccdEnergies[i], gridDPEnergies[i], hybridEnergies[i], multiWarmEnergies[i], branchE,
                        ccdDrop, gdpRec, hybRec, mwRec
                ));
            } else {
                System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%% | %7.1f%%",
                        i, ccdEnergies[i], gridDPEnergies[i], hybridEnergies[i], multiWarmEnergies[i],
                        ccdDrop, gdpRec, hybRec, mwRec));
            }
        }
        System.out.println(new String(new char[105]).replace('\0', '-'));

        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap GridDP:       %.4f kcal/mol", sumGapGridDP / n));
        System.out.println(String.format("  Avg gap Hybrid:       %.4f kcal/mol", sumGapHybrid / n));
        System.out.println(String.format("  Avg gap MultiWarm:    %.4f kcal/mol", sumGapMultiWarm / n));
        System.out.println(String.format("  Avg recovery GridDP:  %.1f%%", sumRecGridDP / n));
        System.out.println(String.format("  Avg recovery Hybrid:  %.1f%%", sumRecHybrid / n));
        System.out.println(String.format("  Avg recovery MWarm:   %.1f%%", sumRecMultiWarm / n));
        if (compareMLPtoForcefield) {
            System.out.println(String.format("  Avg |GridDP - FF|:    %.4f kcal/mol", sumAbsGDPdFF / n));
            System.out.println(String.format("  Avg |Hybrid - FF|:    %.4f kcal/mol", sumAbsHybdFF / n));
            System.out.println(String.format("  Avg |MWarm - FF|:     %.4f kcal/mol", sumAbsMWdFF / n));
        } else if (compareBranchResidualToForcefield) {
            System.out.println(String.format("  Avg |BranchPred - FF|: %.4f kcal/mol", sumAbsBranchdFF / n));
        }

        System.out.println(String.format("\nTiming Summary (%d confs, %d cores, top-%d starts):", n, numCores, multiN));
        System.out.println(String.format("  CCD total:          %10.1f ms", ccdTotalNs / 1e6));
        System.out.println(String.format("  Hybrid total:       %10.1f ms  (speedup vs CCD: %.2fx)", hybridTotalNs / 1e6, (double) ccdTotalNs / hybridTotalNs));
        System.out.println(String.format("  MWarm total:        %10.1f ms  (speedup vs CCD: %.2fx)", multiWarmTotalNs / 1e6, (double) ccdTotalNs / multiWarmTotalNs));
        System.out.println(String.format("\n  MWarm Phase Breakdown:"));
        System.out.println(String.format("    GridDP minimize:  %10.1f ms  (%.1f%%)", totalMWGridDPNs / 1e6, 100.0 * totalMWGridDPNs / multiWarmTotalNs));
        System.out.println(String.format("    Top-N init:       %10.1f ms  (%.1f%%)", totalMWInitNs / 1e6, 100.0 * totalMWInitNs / multiWarmTotalNs));
        System.out.println(String.format("    Top-N pop (%dx): %10.1f ms  (%.1f%%)", multiN, totalMWPopNs / 1e6, 100.0 * totalMWPopNs / multiWarmTotalNs));
        System.out.println(String.format("    CCD obj create:   %10.1f ms  (%.1f%%)", totalMWCCDCreateNs / 1e6, 100.0 * totalMWCCDCreateNs / multiWarmTotalNs));
        System.out.println(String.format("    CCD minimize:     %10.1f ms  (%.1f%%)", totalMWCCDMinNs / 1e6, 100.0 * totalMWCCDMinNs / multiWarmTotalNs));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    @Test
    public void benchmark1CC8CCDOnlyMLPAccuracy() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        SurrogateRunConfig surrogateCfg = SurrogateRunConfig.fromSystemProperties();
        boolean supportsThisBenchmark = surrogateCfg.modelType == SurrogateModelType.MLP
                || surrogateCfg.modelType == SurrogateModelType.BranchResidual;
        if (!surrogateCfg.enabled || !supportsThisBenchmark) {
            throw new IllegalStateException(
                    "benchmark1CC8CCDOnlyMLPAccuracy requires -Dosprey.bench.surrogate.enabled=true and -Dosprey.bench.surrogate.type=mlp|branchResidual"
            );
        }

        int numConfs = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "ccdOnly.numConfs", "15"));
        int ccdMaxIter = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "ccdOnly.maxIter", "30"));
        int topK = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "ccdOnly.topK", "5"));

        System.out.println("=== CCD-only surrogate accuracy benchmark ===");
        System.out.println("Model type: " + surrogateCfg.modelType);
        System.out.println("Conf space positions: " + confSpace.positions.size());
        System.out.println("numConfs=" + numConfs + ", ccdMaxIter=" + ccdMaxIter + ", topK=" + topK);

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();

        System.out.println("Computing energy matrix for A*...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build()
                .calcEnergyMatrix();

        EnergyCalculator rigidEcalc = null;
        EnergyMatrix ematRigid = null;
        if (surrogateCfg.modelType == SurrogateModelType.BranchResidual) {
            rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                    .setParallelism(parallelism)
                    .setIsMinimizing(false)
                    .build();
            ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);
            ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                    .build()
                    .calcEnergyMatrix();
        }

        RCs rcs = new RCs(confSpace);
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional()
                .build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        System.out.println("Got " + scoredConfs.size() + " conformations.");

        MLPSurrogateMatrix mlpSurrogate = null;
        TaskGlobalMLPSurrogate globalMLPSurrogate = null;
        BranchResidualConfEnergyCalculator branchResidualConfEcalc = null;
        if (surrogateCfg.modelType == SurrogateModelType.MLP) {
            if (surrogateCfg.mlpSingleTaskModel) {
                globalMLPSurrogate = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificGlobalMLPSurrogate(
                        confEcalcMinimized,
                        surrogateCfg.cacheRoot,
                        surrogateCfg.taskTag,
                        surrogateCfg.mlpGlobalNumTuples,
                        surrogateCfg.mlpGlobalSamplesPerTuple,
                        surrogateCfg.mlpHidden1,
                        surrogateCfg.mlpHidden2,
                        surrogateCfg.mlpEpochs,
                        surrogateCfg.mlpBatchSize,
                        surrogateCfg.mlpLearningRate
                );
                System.out.println(String.format(
                        "Global MLP surrogate loaded/trained: taskTag=%s cache=%s tuples=%d samplesPerTuple=%d hidden=%dx%d epochs=%d",
                        surrogateCfg.taskTag,
                        surrogateCfg.cacheRoot.getAbsolutePath(),
                        surrogateCfg.mlpGlobalNumTuples,
                        surrogateCfg.mlpGlobalSamplesPerTuple,
                        surrogateCfg.mlpHidden1,
                        surrogateCfg.mlpHidden2,
                        surrogateCfg.mlpEpochs
                ));
            } else {
                mlpSurrogate = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificMLPSurrogate(
                        confEcalcMinimized,
                        surrogateCfg.cacheRoot,
                        surrogateCfg.taskTag,
                        surrogateCfg.samplesPerParam,
                        surrogateCfg.mlpMinSamplesPerModel,
                        surrogateCfg.mlpMaxSamplesPerModel,
                        surrogateCfg.mlpHidden1,
                        surrogateCfg.mlpHidden2,
                        surrogateCfg.mlpEpochs,
                        surrogateCfg.mlpBatchSize,
                        surrogateCfg.mlpLearningRate
                );
                System.out.println(String.format(
                        "MLP surrogate loaded/trained: taskTag=%s cache=%s samplesPerParam=%d minSamples=%d maxSamples=%d hidden=%dx%d epochs=%d",
                        surrogateCfg.taskTag,
                        surrogateCfg.cacheRoot.getAbsolutePath(),
                        surrogateCfg.samplesPerParam,
                        surrogateCfg.mlpMinSamplesPerModel,
                        surrogateCfg.mlpMaxSamplesPerModel,
                        surrogateCfg.mlpHidden1,
                        surrogateCfg.mlpHidden2,
                        surrogateCfg.mlpEpochs
                ));
            }
        } else {
            InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                    confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF
            );
            BranchDecomposition bd = new BranchDecomposition(ig);
            bd.compute();
            RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
            RootedTreeEdge.postOrderCompLlambda(rootedRoot);
            RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
            rootEdge.compactTree();

            ConfSpaceSpecificSurrogateFactory.BranchResidualRidgeConfig branchCfg =
                    new ConfSpaceSpecificSurrogateFactory.BranchResidualRidgeConfig();
            branchCfg.cacheRoot = surrogateCfg.cacheRoot;
            branchCfg.taskTag = surrogateCfg.taskTag;
            branchCfg.forceRetrain = surrogateCfg.forceRetrain;
            branchCfg.featureHashDim = surrogateCfg.branchResidualFeatureHashDim;
            branchCfg.samplerMode = surrogateCfg.branchResidualSamplerMode;
            branchCfg.pretrainGridDPSamples = surrogateCfg.branchResidualPretrainGridDPSamples;
            branchCfg.finetuneCCDSamples = surrogateCfg.branchResidualFineTuneCCDSamples;
            branchCfg.lambdaPre = surrogateCfg.branchResidualLambdaPre;
            branchCfg.lambdaFine = surrogateCfg.branchResidualLambdaFine;
            branchCfg.pretrainUseForcefieldFinalEval = surrogateCfg.branchResidualPretrainUseFinalEval;
            branchCfg.gridSize = surrogateCfg.branchResidualGridSize;
            branchCfg.gridThreads = surrogateCfg.branchResidualGridThreads;
            branchCfg.lowEnergyPoolLimit = surrogateCfg.branchResidualLowEnergyPoolLimit;
            branchCfg.randomSeed = surrogateCfg.branchResidualRandomSeed;
            branchCfg.ridgeMaxIterations = surrogateCfg.branchResidualRidgeMaxIter;
            branchCfg.ridgeTolerance = surrogateCfg.branchResidualRidgeTolerance;

            branchResidualConfEcalc = ConfSpaceSpecificSurrogateFactory.loadOrTrainTaskSpecificBranchResidualRidge(
                    confEcalcMinimized,
                    ematMinimized,
                    ig,
                    rootEdge,
                    ffparams,
                    eref,
                    branchCfg
            );
            System.out.println(String.format(
                    "BranchResidual surrogate loaded/trained: taskTag=%s cache=%s hashDim=%d pre=%d fine=%d",
                    surrogateCfg.taskTag,
                    surrogateCfg.cacheRoot.getAbsolutePath(),
                    surrogateCfg.branchResidualFeatureHashDim,
                    surrogateCfg.branchResidualPretrainGridDPSamples,
                    surrogateCfg.branchResidualFineTuneCCDSamples
            ));
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;

        double[] predEnergies = new double[scoredConfs.size()];
        double[] ffEnergies = new double[scoredConfs.size()];
        double sumAbs = 0.0;
        double sumSq = 0.0;
        double sumBias = 0.0;
        double maxAbs = 0.0;

        String predHeader = surrogateCfg.modelType == SurrogateModelType.MLP ? "MLP_CCD_E" : "Branch_E";
        System.out.println(String.format("%-5s| %11s | %11s | %11s | %8s",
                "Conf", predHeader, "FF_final_E", "Delta", "AbsDelta"));
        System.out.println("-".repeat(66));

        for (int i = 0; i < scoredConfs.size(); i++) {
            RCTuple tuple = new RCTuple(scoredConfs.get(i).getAssignments());
            double predEnergy;
            double ffEnergy;
            if (surrogateCfg.modelType == SurrogateModelType.MLP) {
                ObjectiveWithCleanup obj = makeObjectiveFunction(
                        confSpace, tuple, confEcalcMinimized, null, mlpSurrogate,
                        surrogateCfg, ecalcContext, globalMLPSurrogate
                );
                SimpleCCDMinimizer ccd = new SimpleCCDMinimizer();
                ccd.init(obj.objective);
                ccd.setMaxIterations(ccdMaxIter);
                Minimizer.Result result = ccd.minimizeFromCenter();
                predEnergy = result.energy;
                ffEnergy = evaluateForcefieldEnergyAt(
                        confSpace, confEcalcMinimized, ecalcContext, tuple, result.dofValues
                );
                ccd.close();
                obj.close();
            } else {
                predEnergy = branchResidualConfEcalc.calcEnergy(scoredConfs.get(i).getAssignments());
                ffEnergy = confEcalcMinimized.calcEnergy(tuple).energy;
            }

            double delta = predEnergy - ffEnergy;
            double absDelta = Math.abs(delta);

            predEnergies[i] = predEnergy;
            ffEnergies[i] = ffEnergy;
            sumAbs += absDelta;
            sumSq += delta * delta;
            sumBias += delta;
            maxAbs = Math.max(maxAbs, absDelta);

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %8.4f",
                    i, predEnergy, ffEnergy, delta, absDelta));
        }

        int n = scoredConfs.size();
        double mae = sumAbs / n;
        double rmse = Math.sqrt(sumSq / n);
        double bias = sumBias / n;
        int k = Math.max(1, Math.min(topK, n));
        double recallAtK = topKRecall(predEnergies, ffEnergies, k);
        double rankCorr = spearman(predEnergies, ffEnergies);

        System.out.println("-".repeat(66));
        System.out.println("Surrogate vs Forcefield final-eval error summary:");
        System.out.println(String.format("  MAE:   %.6f kcal/mol", mae));
        System.out.println(String.format("  RMSE:  %.6f kcal/mol", rmse));
        System.out.println(String.format("  Bias:  %.6f kcal/mol", bias));
        System.out.println(String.format("  Max|d| %.6f kcal/mol", maxAbs));
        System.out.println(String.format("  Top-%d recall: %.2f%%", k, recallAtK * 100.0));
        System.out.println(String.format("  Spearman rho: %.4f", rankCorr));
        System.out.println(String.format("  Total wall-clock time: %.2f s", (System.nanoTime() - wallStart) / 1e9));

        minimizingEcalc.close();
        if (rigidEcalc != null) {
            rigidEcalc.close();
        }
    }

    @Test
    public void benchmarkGridSizes() {
        // Same setup but compare g=3,5,7
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                        .build().calcReferenceEnergies()).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized).build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid).build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        // Get 20 confs
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs).setTraditional().build();
        List<ScoredConf> confs = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            confs.add(sc);
        }

        // CCD reference
        double[] ccdEnergies = new double[confs.size()];
        long[] ccdTimesNs = new long[confs.size()];
        for (int i = 0; i < confs.size(); i++) {
            int[] conf = confs.get(i).getAssignments();
            long t0 = System.nanoTime();
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(conf)).energy;
            ccdTimesNs[i] = System.nanoTime() - t0;
        }

        // Test each grid size
        int[] gridSizes = {3, 5, 7};
        for (int g : gridSizes) {
            System.out.println("\n" + "=".repeat(80));
            System.out.println("GRID SIZE g=" + g);
            System.out.println("=".repeat(80));

            GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, g, ffparams);

            double[] gridEnergies = new double[confs.size()];
            long[] precompNs = new long[confs.size()];
            long[] dpNs = new long[confs.size()];

            for (int i = 0; i < confs.size(); i++) {
                GridDPMinimizer.Result r = gridDP.minimize(confs.get(i).getAssignments());
                gridEnergies[i] = r.energy;
                precompNs[i] = r.precomputeTimeNs;
                dpNs[i] = r.dpTimeNs;
            }

            printResults(confs, ccdEnergies, ccdEnergies, gridEnergies, ccdTimesNs, precompNs, dpNs);
        }

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    private void printResults(List<ScoredConf> confs, double[] rigidE, double[] ccdE, double[] gridE,
                              long[] ccdNs, long[] precompNs, long[] dpNs) {
        System.out.println(String.format("%-6s | %12s | %12s | %12s | %10s | %10s | %10s | %12s | %10s | %8s",
                "Conf#", "Rigid_E", "CCD_E", "GridDP_E", "Gap", "Reduction%", "CCD_ms", "Precomp_ms", "DP_ms", "Speedup"));
        System.out.println("-".repeat(135));

        double sumGap = 0;
        double sumCcdMs = 0;
        double sumGridMs = 0;
        double sumReduction = 0;

        for (int i = 0; i < confs.size(); i++) {
            double gap = gridE[i] - ccdE[i];
            double totalDrop = rigidE[i] - ccdE[i];  // CCD improvement over rigid
            double gridDrop = rigidE[i] - gridE[i];   // Grid DP improvement over rigid
            double reductionPct = (totalDrop > 0) ? (gridDrop / totalDrop * 100.0) : 0;
            double ccdMs = ccdNs[i] / 1e6;
            double preMs = precompNs[i] / 1e6;
            double dpMs = dpNs[i] / 1e6;
            double gridTotalMs = preMs + dpMs;
            double speedup = (gridTotalMs > 0) ? ccdMs / gridTotalMs : 0;

            sumGap += gap;
            sumCcdMs += ccdMs;
            sumGridMs += gridTotalMs;
            sumReduction += reductionPct;

            System.out.println(String.format("%-6d | %12.4f | %12.4f | %12.4f | %10.4f | %9.1f%% | %10.2f | %12.2f | %10.2f | %8.2f",
                    i, rigidE[i], ccdE[i], gridE[i], gap, reductionPct, ccdMs, preMs, dpMs, speedup));
        }

        System.out.println("-".repeat(135));
        int n = confs.size();
        System.out.println(String.format("%-6s | %12s | %12s | %12s | %10.4f | %9.1f%% | %10.2f | %12s | %10s | %8.2f",
                "AVG", "", "", "", sumGap / n, sumReduction / n, sumCcdMs / n, "", "", sumCcdMs / sumGridMs));
        System.out.println();
        System.out.println("Summary:");
        System.out.println("  Avg gap (GridDP - CCD):      " + String.format("%.4f", sumGap / n) + " kcal/mol");
        System.out.println("  Avg reduction (GridDP/CCD):   " + String.format("%.1f", sumReduction / n) + "%");
        System.out.println("  Total CCD time:   " + String.format("%.1f", sumCcdMs) + " ms");
        System.out.println("  Total GridDP time: " + String.format("%.1f", sumGridMs) + " ms");
        System.out.println("  Overall speedup:   " + String.format("%.2f", sumCcdMs / sumGridMs) + "x");
    }

    private void printResultsWithCache(List<ScoredConf> confs,
                                         double[] rigidE, double[] ccdE, double[] gridE, double[] cachedE,
                                         long[] ccdNs, long[] precompNs, long[] dpNs,
                                         long[] cachedPrecompNs, long[] cachedDPNs,
                                         int[] cacheHits, int[] cacheLookups,
                                         int[] dpHits, int[] dpLookups) {
        System.out.println(String.format("%-5s | %11s | %11s | %11s | %8s | %8s | %10s | %10s | %8s | %10s | %10s | %8s | %7s | %7s",
                "Conf", "CCD_E", "GridDP_E", "Cached_E", "Gap", "CCD_ms",
                "Precomp_ms", "DP_ms", "Speedup", "CPrecomp_ms", "CDP_ms", "CSpeedup", "EHits", "DPHits"));
        System.out.println("-".repeat(170));

        double sumGap = 0, sumCcdMs = 0, sumGridMs = 0, sumCachedMs = 0, sumReduction = 0;

        for (int i = 0; i < confs.size(); i++) {
            double gap = gridE[i] - ccdE[i];
            double cachedGap = cachedE[i] - gridE[i];  // should be 0
            double totalDrop = rigidE[i] - ccdE[i];
            double gridDrop = rigidE[i] - gridE[i];
            double reductionPct = (totalDrop > 0) ? (gridDrop / totalDrop * 100.0) : 0;
            double ccdMs = ccdNs[i] / 1e6;
            double preMs = precompNs[i] / 1e6;
            double dpMs = dpNs[i] / 1e6;
            double gridTotalMs = preMs + dpMs;
            double speedup = (gridTotalMs > 0) ? ccdMs / gridTotalMs : 0;
            double cPreMs = cachedPrecompNs[i] / 1e6;
            double cDpMs = cachedDPNs[i] / 1e6;
            double cachedTotalMs = cPreMs + cDpMs;
            double cSpeedup = (cachedTotalMs > 0) ? ccdMs / cachedTotalMs : 0;

            sumGap += gap;
            sumCcdMs += ccdMs;
            sumGridMs += gridTotalMs;
            sumCachedMs += cachedTotalMs;
            sumReduction += reductionPct;

            System.out.println(String.format(
                    "%-5d | %11.4f | %11.4f | %11.4f | %8.4f | %8.2f | %10.2f | %10.2f | %8.2f | %10.2f | %10.2f | %8.2f | %3d/%3d | %3d/%3d",
                    i, ccdE[i], gridE[i], cachedE[i], gap, ccdMs,
                    preMs, dpMs, speedup, cPreMs, cDpMs, cSpeedup,
                    cacheHits[i], cacheLookups[i], dpHits[i], dpLookups[i]));

            // Sanity check: cached and uncached should produce identical energies
            if (Math.abs(cachedGap) > 1e-10) {
                System.out.println("  *** WARNING: cached energy differs from uncached by " + cachedGap + " ***");
            }
        }

        int n = confs.size();
        System.out.println("-".repeat(170));
        System.out.println(String.format(
                "%-5s | %11s | %11s | %11s | %8.4f | %8.2f | %10s | %10s | %8.2f | %10s | %10s | %8.2f |",
                "AVG", "", "", "", sumGap / n, sumCcdMs / n,
                "", "", sumCcdMs / sumGridMs, "", "", sumCcdMs / sumCachedMs));
        System.out.println();
        System.out.println("Summary:");
        System.out.println("  Avg gap (GridDP - CCD):       " + String.format("%.4f", sumGap / n) + " kcal/mol");
        System.out.println("  Avg reduction (GridDP/CCD):    " + String.format("%.1f", sumReduction / n) + "%");
        System.out.println("  Total CCD time:      " + String.format("%.1f", sumCcdMs) + " ms");
        System.out.println("  Total GridDP time:    " + String.format("%.1f", sumGridMs) + " ms  (speedup: " + String.format("%.2f", sumCcdMs / sumGridMs) + "x)");
        System.out.println("  Total Cached time:    " + String.format("%.1f", sumCachedMs) + " ms  (speedup: " + String.format("%.2f", sumCcdMs / sumCachedMs) + "x)");
        System.out.println("  Cache benefit:        " + String.format("%.2f", sumGridMs / sumCachedMs) + "x faster than uncached");
    }

    /**
     * Benchmark: CCD (maxIter=30) vs Continuous DP (quadratic relaxation on branch tree).
     * Continuous DP uses gridSize=3 (minimum for quadratic fitting).
     */
    @Test
    public void benchmarkContinuousDP() {
        long wallStart = System.nanoTime();

        // 1. Build conf space: 1CC8, 20 flexible residues
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== Continuous DP Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        // Get conformations
        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            GridDPMinimizer warmGrid = new GridDPMinimizer(confSpace, ig, rootEdge, 3, ffparams, eref, true);
            warmGrid.minimize(warmConf);
            warmGrid.continuousRelax();
        }

        double[] ccdEnergies = new double[n];
        double[] contDPEnergies = new double[n];
        double[] gridDPEnergies = new double[n];
        double[] rigidEnergies = new double[n];
        long[] ccdTimesNs = new long[n];
        long[] contDPTimesNs = new long[n];  // grid precompute + discrete DP + continuous DP + eval
        long[] gridPrecompNs = new long[n];
        long[] contDPOnlyNs = new long[n];
        long[] evalNs = new long[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD (maxIter=30, default)
        for (int i = 0; i < n; i++) {
            long t0 = System.nanoTime();
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
            ccdTimesNs[i] = System.nanoTime() - t0;
        }

        // Continuous DP (gridSize=3)
        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 3, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            long t0 = System.nanoTime();
            GridDPMinimizer.Result gridResult = gridDP.minimize(conf);
            gridDPEnergies[i] = gridResult.energy;
            gridPrecompNs[i] = gridResult.precomputeTimeNs + gridResult.dpTimeNs;

            GridDPMinimizer.ContinuousResult contResult = gridDP.continuousRelax();
            contDPEnergies[i] = contResult.energy;
            contDPOnlyNs[i] = contResult.continuousDPTimeNs;
            evalNs[i] = contResult.evalTimeNs;
            contDPTimesNs[i] = System.nanoTime() - t0;
        }

        // Print results
        System.out.println(String.format("%-5s| %11s | %11s | %11s | %11s | %8s | %8s | %8s | %8s | %8s | %8s",
                "Conf", "CCD_E", "GridDP_E", "ContDP_E", "CCD_drop",
                "GDP_rec%", "CDP_rec%", "CCD_ms", "Grid_ms", "CDP_ms", "Total_ms"));
        System.out.println("-".repeat(130));

        double sumCCDMs = 0, sumContDPMs = 0;
        double sumGDPRec = 0, sumCDPRec = 0;
        double sumGapGDP = 0, sumGapCDP = 0;

        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double gdpDrop = rigidEnergies[i] - gridDPEnergies[i];
            double cdpDrop = rigidEnergies[i] - contDPEnergies[i];
            double gdpRec = (ccdDrop > 0) ? (gdpDrop / ccdDrop * 100.0) : 0;
            double cdpRec = (ccdDrop > 0) ? (cdpDrop / ccdDrop * 100.0) : 0;
            double ccdMs = ccdTimesNs[i] / 1e6;
            double gridMs = gridPrecompNs[i] / 1e6;
            double cdpMs = contDPOnlyNs[i] / 1e6;
            double totalMs = contDPTimesNs[i] / 1e6;

            sumCCDMs += ccdMs;
            sumContDPMs += totalMs;
            sumGDPRec += gdpRec;
            sumCDPRec += cdpRec;
            sumGapGDP += gridDPEnergies[i] - ccdEnergies[i];
            sumGapCDP += contDPEnergies[i] - ccdEnergies[i];

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%% | %8.2f | %8.2f | %8.2f | %8.2f",
                    i, ccdEnergies[i], gridDPEnergies[i], contDPEnergies[i], ccdDrop,
                    gdpRec, cdpRec, ccdMs, gridMs, cdpMs, totalMs));
        }

        System.out.println("-".repeat(130));
        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap GridDP-CCD:     %.4f kcal/mol", sumGapGDP / n));
        System.out.println(String.format("  Avg gap ContDP-CCD:     %.4f kcal/mol", sumGapCDP / n));
        System.out.println(String.format("  Avg recovery GridDP:    %.1f%%", sumGDPRec / n));
        System.out.println(String.format("  Avg recovery ContDP:    %.1f%%", sumCDPRec / n));

        System.out.println(String.format("\nTiming Summary (%d confs):", n));
        System.out.println(String.format("  CCD total:          %10.1f ms", sumCCDMs));
        System.out.println(String.format("  ContDP total:       %10.1f ms  (speedup vs CCD: %.2fx)",
                sumContDPMs, sumCCDMs / sumContDPMs));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    /**
     * Benchmark: CCD (maxIter=30) vs Multi Warm Start + Powell's Conjugate Direction.
     */
    @Test
    public void benchmarkPowellHybrid() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== Powell Hybrid Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;
        int multiN = 20;
        int numCores = Runtime.getRuntime().availableProcessors();

        double[] ccdEnergies = new double[n];
        double[] powellEnergies = new double[n];
        double[] rigidEnergies = new double[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD baseline (maxIter=30)
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        // Multi Warm Start + Powell
        System.out.println("Running Multi Warm Start + Powell (top-" + multiN + " starts, " + numCores + " threads)...");
        ExecutorService executor = Executors.newFixedThreadPool(numCores);
        long powellTotalStart = System.nanoTime();
        long totalGridDPNs = 0, totalInitNs = 0, totalPopNs = 0, totalPowellNs = 0;

        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            // Phase 1: GridDP
            long t0 = System.nanoTime();
            gridDP.minimize(conf);
            totalGridDPNs += System.nanoTime() - t0;

            // Phase 2: Init enumeration
            t0 = System.nanoTime();
            gridDP.initTopNEnumeration();
            totalInitNs += System.nanoTime() - t0;

            // Phase 3: Pop top-N
            t0 = System.nanoTime();
            List<double[]> topNDOFs = new ArrayList<>();
            for (int k = 0; k < multiN; k++) {
                double[] dofs = gridDP.popNextDOFValues();
                if (dofs == null) break;
                topNDOFs.add(dofs);
            }
            totalPopNs += System.nanoTime() - t0;

            // Phase 4: Parallel Powell minimization
            t0 = System.nanoTime();
            int nStarts = topNDOFs.size();
            List<Future<Double>> futures = new ArrayList<>();

            for (int k = 0; k < nStarts; k++) {
                final double[] startDOFs = topNDOFs.get(k);
                futures.add(executor.submit(() -> {
                    ParametricMolecule pmolK = confSpace.makeMolecule(tuple);
                    ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
                    EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmolK.mol);
                    MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmolK, efunc);

                    int nDofs = pmolK.dofs.size();
                    cern.colt.matrix.DoubleMatrix1D[] constraints = objFunc.getConstraints();
                    double[] lo = new double[nDofs];
                    double[] hi = new double[nDofs];
                    for (int d = 0; d < nDofs; d++) {
                        lo[d] = constraints[0].get(d);
                        hi[d] = constraints[1].get(d);
                    }

                    // Clamp start point to bounds
                    double[] start = startDOFs.clone();
                    for (int d = 0; d < nDofs; d++) {
                        start[d] = Math.max(lo[d], Math.min(hi[d], start[d]));
                    }

                    // Wrap with bounds clamping
                    MultivariateFunction powellFunc = (double[] point) -> {
                        cern.colt.matrix.DoubleMatrix1D x = cern.colt.matrix.DoubleFactory1D.dense.make(nDofs);
                        for (int d = 0; d < nDofs; d++) {
                            x.set(d, Math.max(lo[d], Math.min(hi[d], point[d])));
                        }
                        return objFunc.getValue(x);
                    };

                    PowellOptimizer powell = new PowellOptimizer(1e-4, 1e-6);
                    try {
                        PointValuePair result = powell.optimize(
                                new MaxEval(500),
                                new ObjectiveFunction(powellFunc),
                                GoalType.MINIMIZE,
                                new InitialGuess(start));
                        double energy = result.getValue();
                        efunc.close();
                        return energy;
                    } catch (Exception ex) {
                        // Fallback: evaluate at start point
                        double energy = powellFunc.value(start);
                        efunc.close();
                        return energy;
                    }
                }));
            }

            double bestE = Double.POSITIVE_INFINITY;
            for (Future<Double> f : futures) {
                try {
                    double e = f.get();
                    if (e < bestE) bestE = e;
                } catch (Exception ex) {
                    throw new RuntimeException("Powell failed", ex);
                }
            }
            powellEnergies[i] = bestE;
            totalPowellNs += System.nanoTime() - t0;
        }
        executor.shutdown();
        long powellTotalNs = System.nanoTime() - powellTotalStart;

        // Print results
        System.out.println(String.format("\n%-5s| %11s | %11s | %8s | %8s | %8s",
                "Conf", "CCD_E", "Powell_E", "CCD_drop", "CCD_rec%", "Pow_rec%"));
        System.out.println("-".repeat(70));

        double sumCCDRec = 0, sumPowRec = 0, sumGapPow = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double powDrop = rigidEnergies[i] - powellEnergies[i];
            double powRec = (ccdDrop > 0) ? (powDrop / ccdDrop * 100.0) : 0;
            sumCCDRec += 100.0;
            sumPowRec += powRec;
            sumGapPow += powellEnergies[i] - ccdEnergies[i];

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%%",
                    i, ccdEnergies[i], powellEnergies[i], ccdDrop, 100.0, powRec));
        }
        System.out.println("-".repeat(70));

        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap Powell-CCD:     %.4f kcal/mol", sumGapPow / n));
        System.out.println(String.format("  Avg recovery CCD:      100.0%%"));
        System.out.println(String.format("  Avg recovery Powell:    %.1f%%", sumPowRec / n));

        System.out.println(String.format("\nTiming Summary (%d confs, %d cores, top-%d starts):", n, numCores, multiN));
        System.out.println(String.format("  CCD total:          %10.1f ms", ccdTotalNs / 1e6));
        System.out.println(String.format("  Powell total:       %10.1f ms  (speedup vs CCD: %.2fx)",
                powellTotalNs / 1e6, (double) ccdTotalNs / powellTotalNs));
        System.out.println(String.format("\n  Powell Phase Breakdown:"));
        System.out.println(String.format("    GridDP minimize:  %10.1f ms", totalGridDPNs / 1e6));
        System.out.println(String.format("    Top-N init:       %10.1f ms", totalInitNs / 1e6));
        System.out.println(String.format("    Top-N pop:        %10.1f ms", totalPopNs / 1e6));
        System.out.println(String.format("    Powell minimize:  %10.1f ms", totalPowellNs / 1e6));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    @Test
    public void benchmarkHybridIterations() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== Hybrid vs CCD Iteration Comparison ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;

        double[] rigidEnergies = new double[n];
        double[] ccdEnergies = new double[n];
        double[] hybridEnergies = new double[n];
        int[] ccdIters = new int[n];
        int[] hybridIters = new int[n];
        long[] ccdTimesNs = new long[n];
        long[] hybridTimesNs = new long[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD from center (maxIter=30, SurfingLineSearcher) — track iterations
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);
            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmol, efunc);
            SimpleCCDMinimizer ccd = new SimpleCCDMinimizer();
            ccd.init(objFunc);
            // maxIter=30 (default)

            long t0 = System.nanoTime();
            Minimizer.Result result = ccd.minimizeFromCenter();
            ccdTimesNs[i] = System.nanoTime() - t0;
            ccdEnergies[i] = result.energy;
            ccdIters[i] = ccd.getLastIterations();

            efunc.close();
            ccd.clean();
        }

        // Hybrid: GridDP start + CCD (maxIter=30, SurfingLineSearcher) — track iterations
        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            long t0 = System.nanoTime();
            GridDPMinimizer.Result gridResult = gridDP.minimize(conf);

            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmol, efunc);
            SimpleCCDMinimizer ccd = new SimpleCCDMinimizer();
            ccd.init(objFunc);
            // maxIter=30 (default), same as CCD baseline

            cern.colt.matrix.DoubleMatrix1D startx = cern.colt.matrix.DoubleFactory1D.dense.make(pmol.dofs.size());
            double[] gridDOFs = gridResult.bestDOFValues;
            if (gridDOFs != null) {
                for (int d = 0; d < pmol.dofs.size(); d++) {
                    startx.set(d, gridDOFs[d]);
                }
            } else {
                pmol.dofBounds.getCenter(startx);
            }

            Minimizer.Result result = ccd.minimizeFrom(startx);
            hybridTimesNs[i] = System.nanoTime() - t0;
            hybridEnergies[i] = result.energy;
            hybridIters[i] = ccd.getLastIterations();

            efunc.close();
            ccd.clean();
        }

        // Print results
        System.out.println(String.format("%-5s| %11s | %11s | %8s | %8s | %8s | %8s | %10s | %10s",
                "Conf", "CCD_E", "Hybrid_E", "CCD_drop", "Hyb_rec%", "CCD_iter", "Hyb_iter", "CCD_ms", "Hyb_ms"));
        System.out.println("-".repeat(105));

        double sumRec = 0;
        int sumCCDIter = 0, sumHybIter = 0;
        double sumCCDMs = 0, sumHybMs = 0;
        double sumGap = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double hybDrop = rigidEnergies[i] - hybridEnergies[i];
            double hybRec = (ccdDrop > 0) ? (hybDrop / ccdDrop * 100.0) : 0;
            sumRec += hybRec;
            sumCCDIter += ccdIters[i];
            sumHybIter += hybridIters[i];
            sumCCDMs += ccdTimesNs[i] / 1e6;
            sumHybMs += hybridTimesNs[i] / 1e6;
            sumGap += hybridEnergies[i] - ccdEnergies[i];

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %7.1f%% | %8d | %8d | %10.1f | %10.1f",
                    i, ccdEnergies[i], hybridEnergies[i], ccdDrop, hybRec,
                    ccdIters[i], hybridIters[i], ccdTimesNs[i] / 1e6, hybridTimesNs[i] / 1e6));
        }
        System.out.println("-".repeat(105));

        System.out.println(String.format("\nSummary (%d confs):", n));
        System.out.println(String.format("  Avg gap Hybrid-CCD:    %.4f kcal/mol", sumGap / n));
        System.out.println(String.format("  Avg recovery Hybrid:   %.1f%%", sumRec / n));
        System.out.println(String.format("  Avg CCD iterations:    %.1f", (double) sumCCDIter / n));
        System.out.println(String.format("  Avg Hybrid iterations: %.1f", (double) sumHybIter / n));
        System.out.println(String.format("  Iteration reduction:   %.1f%%", 100.0 * (1.0 - (double) sumHybIter / sumCCDIter)));
        System.out.println(String.format("  Total CCD time:        %.1f ms", sumCCDMs));
        System.out.println(String.format("  Total Hybrid time:     %.1f ms  (speedup: %.2fx)", sumHybMs, sumCCDMs / sumHybMs));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    @Test
    public void benchmarkBOBYQAHybrid() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== BOBYQA Hybrid Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;
        int multiN = 20;
        int numCores = Runtime.getRuntime().availableProcessors();

        double[] ccdEnergies = new double[n];
        double[] bobyqaEnergies = new double[n];
        double[] rigidEnergies = new double[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD baseline (maxIter=30)
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        // Multi Warm Start + BOBYQA
        System.out.println("Running Multi Warm Start + BOBYQA (top-" + multiN + " starts, " + numCores + " threads)...");
        ExecutorService executor = Executors.newFixedThreadPool(numCores);
        long bobyqaTotalStart = System.nanoTime();
        long totalGridDPNs = 0, totalInitNs = 0, totalPopNs = 0, totalBOBYQANs = 0;

        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            // Phase 1: GridDP
            long t0 = System.nanoTime();
            gridDP.minimize(conf);
            totalGridDPNs += System.nanoTime() - t0;

            // Phase 2: Init enumeration
            t0 = System.nanoTime();
            gridDP.initTopNEnumeration();
            totalInitNs += System.nanoTime() - t0;

            // Phase 3: Pop top-N
            t0 = System.nanoTime();
            List<double[]> topNDOFs = new ArrayList<>();
            for (int k = 0; k < multiN; k++) {
                double[] dofs = gridDP.popNextDOFValues();
                if (dofs == null) break;
                topNDOFs.add(dofs);
            }
            totalPopNs += System.nanoTime() - t0;

            // Phase 4: Parallel BOBYQA minimization with native bounds
            t0 = System.nanoTime();
            int nStarts = topNDOFs.size();
            List<Future<Double>> futures = new ArrayList<>();

            for (int k = 0; k < nStarts; k++) {
                final double[] startDOFs = topNDOFs.get(k);
                futures.add(executor.submit(() -> {
                    ParametricMolecule pmolK = confSpace.makeMolecule(tuple);
                    ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
                    EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmolK.mol);
                    MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmolK, efunc);

                    int nDofs = pmolK.dofs.size();
                    cern.colt.matrix.DoubleMatrix1D[] constraints = objFunc.getConstraints();
                    double[] lo = new double[nDofs];
                    double[] hi = new double[nDofs];
                    for (int d = 0; d < nDofs; d++) {
                        lo[d] = constraints[0].get(d);
                        hi[d] = constraints[1].get(d);
                    }

                    // Clamp start point to bounds (strictly inside for BOBYQA)
                    double[] start = startDOFs.clone();
                    for (int d = 0; d < nDofs; d++) {
                        double eps = 1e-10 * (hi[d] - lo[d]);
                        start[d] = Math.max(lo[d] + eps, Math.min(hi[d] - eps, start[d]));
                    }

                    // BOBYQA objective: evaluate forcefield with DOFs
                    MultivariateFunction bobyqaFunc = (double[] point) -> {
                        cern.colt.matrix.DoubleMatrix1D x = cern.colt.matrix.DoubleFactory1D.dense.make(nDofs);
                        for (int d = 0; d < nDofs; d++) {
                            x.set(d, point[d]);
                        }
                        return objFunc.getValue(x);
                    };

                    // BOBYQA: interpolation points = 2*n+1 (recommended default)
                    int interpPoints = 2 * nDofs + 1;
                    BOBYQAOptimizer bobyqa = new BOBYQAOptimizer(interpPoints);
                    try {
                        PointValuePair result = bobyqa.optimize(
                                new MaxEval(500),
                                new ObjectiveFunction(bobyqaFunc),
                                GoalType.MINIMIZE,
                                new InitialGuess(start),
                                new SimpleBounds(lo, hi));
                        double energy = result.getValue();
                        efunc.close();
                        return energy;
                    } catch (Exception ex) {
                        // Fallback: evaluate at start point
                        double energy = bobyqaFunc.value(start);
                        efunc.close();
                        return energy;
                    }
                }));
            }

            double bestE = Double.POSITIVE_INFINITY;
            for (Future<Double> f : futures) {
                try {
                    double e = f.get();
                    if (e < bestE) bestE = e;
                } catch (Exception ex) {
                    throw new RuntimeException("BOBYQA failed", ex);
                }
            }
            bobyqaEnergies[i] = bestE;
            totalBOBYQANs += System.nanoTime() - t0;
        }
        executor.shutdown();
        long bobyqaTotalNs = System.nanoTime() - bobyqaTotalStart;

        // Print results
        System.out.println(String.format("\n%-5s| %11s | %11s | %8s | %8s | %8s",
                "Conf", "CCD_E", "BOBYQA_E", "CCD_drop", "CCD_rec%", "BOB_rec%"));
        System.out.println("-".repeat(70));

        double sumCCDRec = 0, sumBobRec = 0, sumGapBob = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double bobDrop = rigidEnergies[i] - bobyqaEnergies[i];
            double bobRec = (ccdDrop > 0) ? (bobDrop / ccdDrop * 100.0) : 0;
            sumCCDRec += 100.0;
            sumBobRec += bobRec;
            sumGapBob += bobyqaEnergies[i] - ccdEnergies[i];

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%%",
                    i, ccdEnergies[i], bobyqaEnergies[i], ccdDrop, 100.0, bobRec));
        }
        System.out.println("-".repeat(70));

        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap BOBYQA-CCD:     %.4f kcal/mol", sumGapBob / n));
        System.out.println(String.format("  Avg recovery CCD:      100.0%%"));
        System.out.println(String.format("  Avg recovery BOBYQA:    %.1f%%", sumBobRec / n));

        System.out.println(String.format("\nTiming Summary (%d confs, %d cores, top-%d starts):", n, numCores, multiN));
        System.out.println(String.format("  CCD total:          %10.1f ms", ccdTotalNs / 1e6));
        System.out.println(String.format("  BOBYQA total:       %10.1f ms  (speedup vs CCD: %.2fx)",
                bobyqaTotalNs / 1e6, (double) ccdTotalNs / bobyqaTotalNs));
        System.out.println(String.format("\n  BOBYQA Phase Breakdown:"));
        System.out.println(String.format("    GridDP minimize:  %10.1f ms", totalGridDPNs / 1e6));
        System.out.println(String.format("    Top-N init:       %10.1f ms", totalInitNs / 1e6));
        System.out.println(String.format("    Top-N pop:        %10.1f ms", totalPopNs / 1e6));
        System.out.println(String.format("    BOBYQA minimize:  %10.1f ms", totalBOBYQANs / 1e6));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    @Test
    public void benchmarkWallJumpHybrid() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== Wall Jump Hybrid Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;
        int multiN = 20;
        int numCores = Runtime.getRuntime().availableProcessors();

        double[] ccdEnergies = new double[n];
        double[] wjEnergies = new double[n];
        double[] rigidEnergies = new double[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD baseline (maxIter=30, SurfingLineSearcher)
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        // Multi Warm Start + WallJumping
        System.out.println("Running Multi Warm Start + WallJumping (top-" + multiN + " starts, " + numCores + " threads)...");
        ExecutorService executor = Executors.newFixedThreadPool(numCores);
        long wjTotalStart = System.nanoTime();
        long totalGridDPNs = 0, totalInitNs = 0, totalPopNs = 0, totalWJNs = 0;

        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            // Phase 1: GridDP
            long t0 = System.nanoTime();
            gridDP.minimize(conf);
            totalGridDPNs += System.nanoTime() - t0;

            // Phase 2: Init enumeration
            t0 = System.nanoTime();
            gridDP.initTopNEnumeration();
            totalInitNs += System.nanoTime() - t0;

            // Phase 3: Pop top-N
            t0 = System.nanoTime();
            List<double[]> topNDOFs = new ArrayList<>();
            for (int k = 0; k < multiN; k++) {
                double[] dofs = gridDP.popNextDOFValues();
                if (dofs == null) break;
                topNDOFs.add(dofs);
            }
            totalPopNs += System.nanoTime() - t0;

            // Phase 4: Parallel WallJumping CCD
            t0 = System.nanoTime();
            int nStarts = topNDOFs.size();
            MoleculeObjectiveFunction[] objFuncs = new MoleculeObjectiveFunction[nStarts];
            SimpleCCDMinimizer[] ccds = new SimpleCCDMinimizer[nStarts];
            EnergyFunction[] efuncs = new EnergyFunction[nStarts];
            cern.colt.matrix.DoubleMatrix1D[] startxs = new cern.colt.matrix.DoubleMatrix1D[nStarts];
            for (int k = 0; k < nStarts; k++) {
                ParametricMolecule pmolK = confSpace.makeMolecule(tuple);
                ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
                efuncs[k] = ecalcContext.efuncs.make(inters, pmolK.mol);
                objFuncs[k] = new MoleculeObjectiveFunction(pmolK, efuncs[k]);
                ccds[k] = new SimpleCCDMinimizer((context) -> new WallJumpingLineSearcher());
                ccds[k].init(objFuncs[k]);
                ccds[k].setMaxIterations(2);
                startxs[k] = cern.colt.matrix.DoubleFactory1D.dense.make(pmolK.dofs.size());
                double[] startDOFs = topNDOFs.get(k);
                for (int d = 0; d < pmolK.dofs.size(); d++) {
                    startxs[k].set(d, startDOFs[d]);
                }
            }

            // Submit minimizeFrom() to thread pool
            List<Future<Double>> futures = new ArrayList<>();
            for (int k = 0; k < nStarts; k++) {
                final int kk = k;
                futures.add(executor.submit(() -> {
                    Minimizer.Result r = ccds[kk].minimizeFrom(startxs[kk]);
                    return r.energy;
                }));
            }
            double bestE = Double.POSITIVE_INFINITY;
            for (Future<Double> f : futures) {
                try {
                    double e = f.get();
                    if (e < bestE) bestE = e;
                } catch (Exception ex) {
                    throw new RuntimeException("WallJump CCD failed", ex);
                }
            }
            wjEnergies[i] = bestE;
            totalWJNs += System.nanoTime() - t0;

            // Cleanup
            for (int k = 0; k < nStarts; k++) {
                efuncs[k].close();
                ccds[k].close();
            }
        }
        executor.shutdown();
        long wjTotalNs = System.nanoTime() - wjTotalStart;

        // Print results
        System.out.println(String.format("\n%-5s| %11s | %11s | %8s | %8s | %8s",
                "Conf", "CCD_E", "WJ_E", "CCD_drop", "CCD_rec%", "WJ_rec%"));
        System.out.println("-".repeat(70));

        double sumWJRec = 0, sumGapWJ = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double wjDrop = rigidEnergies[i] - wjEnergies[i];
            double wjRec = (ccdDrop > 0) ? (wjDrop / ccdDrop * 100.0) : 0;
            sumWJRec += wjRec;
            sumGapWJ += wjEnergies[i] - ccdEnergies[i];

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%%",
                    i, ccdEnergies[i], wjEnergies[i], ccdDrop, 100.0, wjRec));
        }
        System.out.println("-".repeat(70));

        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap WJ-CCD:        %.4f kcal/mol", sumGapWJ / n));
        System.out.println(String.format("  Avg recovery CCD:      100.0%%"));
        System.out.println(String.format("  Avg recovery WJ:       %.1f%%", sumWJRec / n));

        System.out.println(String.format("\nTiming Summary (%d confs, %d cores, top-%d starts):", n, numCores, multiN));
        System.out.println(String.format("  CCD total:          %10.1f ms", ccdTotalNs / 1e6));
        System.out.println(String.format("  WJ total:           %10.1f ms  (speedup vs CCD: %.2fx)",
                wjTotalNs / 1e6, (double) ccdTotalNs / wjTotalNs));
        System.out.println(String.format("\n  WJ Phase Breakdown:"));
        System.out.println(String.format("    GridDP minimize:  %10.1f ms", totalGridDPNs / 1e6));
        System.out.println(String.format("    Top-N init:       %10.1f ms", totalInitNs / 1e6));
        System.out.println(String.format("    Top-N pop:        %10.1f ms", totalPopNs / 1e6));
        System.out.println(String.format("    WJ minimize:      %10.1f ms", totalWJNs / 1e6));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    /**
     * Benchmark: CCD vs GridDP+L-BFGS vs L-BFGS from center.
     * Tests the SimpleLBFGSMinimizer with sparse numerical gradients.
     */
    @Test
    public void benchmarkLBFGSHybrid() {
        long wallStart = System.nanoTime();

        // 1. Build conf space: 1CC8, 20 flexible residues
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A33").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A47").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A50").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A51").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== L-BFGS Hybrid Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");

        Parallelism parallelism = Parallelism.makeCpu(1);  // single thread for fair timing
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        int numConfs = 15;
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;

        double[] rigidEnergies = new double[n];
        double[] ccdEnergies = new double[n];
        double[] lbfgsCenterEnergies = new double[n];
        double[] lbfgsHybridEnergies = new double[n];
        int[] ccdIters = new int[n];
        int[] lbfgsCenterIters = new int[n];
        int[] lbfgsHybridIters = new int[n];
        long[] ccdTimesNs = new long[n];
        long[] lbfgsCenterTimesNs = new long[n];
        long[] lbfgsHybridTimesNs = new long[n]; // includes GridDP time

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD baseline (maxIter=30, SurfingLineSearcher) — track iterations
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);
            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmol, efunc);
            SimpleCCDMinimizer ccd = new SimpleCCDMinimizer();
            ccd.init(objFunc);

            long t0 = System.nanoTime();
            Minimizer.Result result = ccd.minimizeFromCenter();
            ccdTimesNs[i] = System.nanoTime() - t0;
            ccdEnergies[i] = result.energy;
            ccdIters[i] = ccd.getLastIterations();

            efunc.close();
            ccd.clean();
        }

        // L-BFGS from center (maxIter=30) — no GridDP warm start
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);
            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmol, efunc);
            SimpleLBFGSMinimizer lbfgs = new SimpleLBFGSMinimizer();
            lbfgs.init(objFunc);

            long t0 = System.nanoTime();
            Minimizer.Result result = lbfgs.minimizeFromCenter();
            lbfgsCenterTimesNs[i] = System.nanoTime() - t0;
            lbfgsCenterEnergies[i] = result.energy;
            lbfgsCenterIters[i] = lbfgs.getLastIterations();

            efunc.close();
            lbfgs.clean();
        }

        // Hybrid: GridDP(g=2) start → L-BFGS (maxIter=30)
        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            long t0 = System.nanoTime();
            GridDPMinimizer.Result gridResult = gridDP.minimize(conf);

            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            MoleculeObjectiveFunction objFunc = new MoleculeObjectiveFunction(pmol, efunc);
            SimpleLBFGSMinimizer lbfgs = new SimpleLBFGSMinimizer();
            lbfgs.init(objFunc);

            cern.colt.matrix.DoubleMatrix1D startx = cern.colt.matrix.DoubleFactory1D.dense.make(pmol.dofs.size());
            double[] gridDOFs = gridResult.bestDOFValues;
            if (gridDOFs != null) {
                for (int d = 0; d < pmol.dofs.size(); d++) {
                    startx.set(d, gridDOFs[d]);
                }
            } else {
                pmol.dofBounds.getCenter(startx);
            }

            Minimizer.Result result = lbfgs.minimizeFrom(startx);
            lbfgsHybridTimesNs[i] = System.nanoTime() - t0;
            lbfgsHybridEnergies[i] = result.energy;
            lbfgsHybridIters[i] = lbfgs.getLastIterations();

            efunc.close();
            lbfgs.clean();
        }

        // Print comparison table
        System.out.println(String.format("%-5s| %11s | %11s | %11s | %8s | %8s | %8s | %6s | %6s | %6s | %10s | %10s | %10s",
                "Conf", "CCD_E", "LBFGS_ctr_E", "LBFGS_hyb_E", "CCD_drop",
                "Lctr_rec%", "Lhyb_rec%", "CCD_it", "Lc_it", "Lh_it",
                "CCD_ms", "Lctr_ms", "Lhyb_ms"));
        System.out.println("-".repeat(145));

        double sumCCDRec = 0, sumLctrRec = 0, sumLhybRec = 0;
        double sumGapLctr = 0, sumGapLhyb = 0;
        int sumCCDIter = 0, sumLctrIter = 0, sumLhybIter = 0;
        double sumCCDMs = 0, sumLctrMs = 0, sumLhybMs = 0;

        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double lctrDrop = rigidEnergies[i] - lbfgsCenterEnergies[i];
            double lhybDrop = rigidEnergies[i] - lbfgsHybridEnergies[i];
            double lctrRec = (ccdDrop > 0) ? (lctrDrop / ccdDrop * 100.0) : 0;
            double lhybRec = (ccdDrop > 0) ? (lhybDrop / ccdDrop * 100.0) : 0;

            sumCCDRec += 100.0;
            sumLctrRec += lctrRec;
            sumLhybRec += lhybRec;
            sumGapLctr += lbfgsCenterEnergies[i] - ccdEnergies[i];
            sumGapLhyb += lbfgsHybridEnergies[i] - ccdEnergies[i];
            sumCCDIter += ccdIters[i];
            sumLctrIter += lbfgsCenterIters[i];
            sumLhybIter += lbfgsHybridIters[i];
            sumCCDMs += ccdTimesNs[i] / 1e6;
            sumLctrMs += lbfgsCenterTimesNs[i] / 1e6;
            sumLhybMs += lbfgsHybridTimesNs[i] / 1e6;

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %8.4f | %7.1f%% | %7.1f%% | %6d | %6d | %6d | %10.1f | %10.1f | %10.1f",
                    i, ccdEnergies[i], lbfgsCenterEnergies[i], lbfgsHybridEnergies[i], ccdDrop,
                    lctrRec, lhybRec, ccdIters[i], lbfgsCenterIters[i], lbfgsHybridIters[i],
                    ccdTimesNs[i] / 1e6, lbfgsCenterTimesNs[i] / 1e6, lbfgsHybridTimesNs[i] / 1e6));
        }
        System.out.println("-".repeat(145));

        System.out.println(String.format("\nEnergy Summary (%d confs):", n));
        System.out.println(String.format("  Avg gap LBFGS_center - CCD:  %.4f kcal/mol", sumGapLctr / n));
        System.out.println(String.format("  Avg gap LBFGS_hybrid - CCD:  %.4f kcal/mol", sumGapLhyb / n));
        System.out.println(String.format("  Avg recovery CCD:           100.0%%"));
        System.out.println(String.format("  Avg recovery LBFGS_center:   %.1f%%", sumLctrRec / n));
        System.out.println(String.format("  Avg recovery LBFGS_hybrid:   %.1f%%", sumLhybRec / n));

        System.out.println(String.format("\nIteration Summary:"));
        System.out.println(String.format("  Avg CCD iterations:          %.1f", (double) sumCCDIter / n));
        System.out.println(String.format("  Avg LBFGS_center iterations: %.1f", (double) sumLctrIter / n));
        System.out.println(String.format("  Avg LBFGS_hybrid iterations: %.1f", (double) sumLhybIter / n));

        System.out.println(String.format("\nTiming Summary:"));
        System.out.println(String.format("  CCD total:           %10.1f ms", sumCCDMs));
        System.out.println(String.format("  LBFGS_center total:  %10.1f ms  (speedup vs CCD: %.2fx)", sumLctrMs, sumCCDMs / sumLctrMs));
        System.out.println(String.format("  LBFGS_hybrid total:  %10.1f ms  (speedup vs CCD: %.2fx)", sumLhybMs, sumCCDMs / sumLhybMs));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    /**
     * Compare three minimization approaches on the same conformations:
     *   1. CCD original (maxIter=30, default)
     *   2. EPIC (polynomial surrogate + CCD)
     *   3. Multi Warm Start Hybrid GridDP (top-20 starts → CCD maxIter=3)
     *
     * Reports energy quality (gap vs CCD), timing, and EPIC precomputation cost.
     */
    @Test
    public void benchmarkCCDvsEPICvsMultiWarm() {
        long wallStart = System.nanoTime();

        // 1. Build conf space: 1CC8, 8 flexible residues (reduced for EPIC)
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== CCD vs EPIC vs Multi Warm Start Benchmark ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }

        Parallelism parallelism = Parallelism.makeCpu(1);

        // 2. Compute energy matrices
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();
        System.out.println("Energy matrices computed.");

        // 3. Build interaction graph + branch decomposition (for GridDP)
        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);
        System.out.println("Interaction graph: " + ig.getNumPositions() + " positions, "
                + countEdges(ig) + " edges");

        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        System.out.println("Branch decomposition: branchwidth=" + bd.getBranchwidth());

        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        // 4. Compute EPIC matrix (timed separately as one-time precomputation)
        // System.out.println("\nComputing EPIC matrix...");
        // EPICSettings epicSettings = EPICSettings.defaultEPIC();
        // // Disable SAPE — its deep-copy fails with SimpleConfSpace molecules
        // try {
        //     java.lang.reflect.Field f = EPICSettings.class.getDeclaredField("useSAPE");
        //     f.setAccessible(true);
        //     f.setBoolean(epicSettings, false);
        // } catch (Exception e) {
        //     throw new RuntimeException("Failed to disable SAPE", e);
        // }
        // PruningMatrix pruneMat = new PruningMatrix(confSpace);
        // long epicPrecompStart = System.nanoTime();
        // NewEPICMatrixCalculator epicCalc = new NewEPICMatrixCalculator(
        //         confSpace, confEcalcMinimized, pruneMat, epicSettings);
        // epicCalc.calcPEM();
        // NewEPICMatrix epicMat = epicCalc.getEPICMatrix();
        // long epicPrecompNs = System.nanoTime() - epicPrecompStart;
        // System.out.println("EPIC precomputation: " + String.format("%.1f", epicPrecompNs / 1e6) + " ms");

        // 5. Get N conformations via A*
        int numConfs = 1;
        System.out.println("\nFinding " + numConfs + " conformations via A*...");
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            // epicMat.minimizeEnergy(warmConf);
            new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true).minimize(warmConf);
        }

        EnergyCalculator.Type.Context ecalcContext = minimizingEcalc.context;
        int multiN = 20;
        int ccdIterForHybrid = 3;
        int numCores = Runtime.getRuntime().availableProcessors();

        double[] rigidEnergies = new double[n];
        double[] ccdEnergies = new double[n];
        // double[] epicEnergies = new double[n];
        double[] mwarmEnergies = new double[n];
        long[] ccdTimesNs = new long[n];
        // long[] epicTimesNs = new long[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // ===== Method 1: CCD original (maxIter=30, default) =====
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            long t0 = System.nanoTime();
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
            ccdTimesNs[i] = System.nanoTime() - t0;
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        // ===== Method 2: EPIC (polynomial surrogate + CCD minimization) =====
        // long epicTotalStart = System.nanoTime();
        // for (int i = 0; i < n; i++) {
        //     long t0 = System.nanoTime();
        //     epicEnergies[i] = epicMat.minimizeEnergy(scoredConfs.get(i).getAssignments());
        //     epicTimesNs[i] = System.nanoTime() - t0;
        // }
        // long epicTotalNs = System.nanoTime() - epicTotalStart;

        // ===== Method 3: Multi Warm Start Hybrid GridDP (top-20 starts, CCD maxIter=3) =====
        System.out.println("Running Multi Warm Start Hybrid (top-" + multiN + " starts, CCD iter="
                + ccdIterForHybrid + ", " + numCores + " threads)...");
        ExecutorService executor = Executors.newFixedThreadPool(numCores);
        long mwarmTotalStart = System.nanoTime();
        long totalMWGridDPNs = 0, totalMWInitNs = 0, totalMWPopNs = 0;
        long totalMWCCDCreateNs = 0, totalMWCCDMinNs = 0;

        GridDPMinimizer gridDP = new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            RCTuple tuple = new RCTuple(conf);

            // Phase 1: GridDP minimize
            long t0 = System.nanoTime();
            gridDP.minimize(conf);
            totalMWGridDPNs += System.nanoTime() - t0;

            // Phase 2: Init enumeration
            t0 = System.nanoTime();
            gridDP.initTopNEnumeration();
            totalMWInitNs += System.nanoTime() - t0;

            // Phase 3: Pop top-N starting points
            t0 = System.nanoTime();
            List<double[]> topNDOFs = new ArrayList<>();
            for (int k = 0; k < multiN; k++) {
                double[] dofs = gridDP.popNextDOFValues();
                if (dofs == null) break;
                topNDOFs.add(dofs);
            }
            totalMWPopNs += System.nanoTime() - t0;

            // Phase 4: Pre-create CCD objects
            t0 = System.nanoTime();
            int nStarts = topNDOFs.size();
            MoleculeObjectiveFunction[] objFuncs = new MoleculeObjectiveFunction[nStarts];
            SimpleCCDMinimizer[] ccds = new SimpleCCDMinimizer[nStarts];
            EnergyFunction[] efuncs = new EnergyFunction[nStarts];
            cern.colt.matrix.DoubleMatrix1D[] startxs = new cern.colt.matrix.DoubleMatrix1D[nStarts];
            for (int k = 0; k < nStarts; k++) {
                ParametricMolecule pmolK = confSpace.makeMolecule(tuple);
                ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
                efuncs[k] = ecalcContext.efuncs.make(inters, pmolK.mol);
                objFuncs[k] = new MoleculeObjectiveFunction(pmolK, efuncs[k]);
                ccds[k] = new SimpleCCDMinimizer((context) -> new QuadraticLineSearcher());
                ccds[k].init(objFuncs[k]);
                ccds[k].setMaxIterations(ccdIterForHybrid);
                startxs[k] = cern.colt.matrix.DoubleFactory1D.dense.make(pmolK.dofs.size());
                double[] startDOFs = topNDOFs.get(k);
                for (int d = 0; d < pmolK.dofs.size(); d++) {
                    startxs[k].set(d, startDOFs[d]);
                }
            }
            totalMWCCDCreateNs += System.nanoTime() - t0;

            // Phase 5: Parallel CCD minimization
            long t1 = System.nanoTime();
            List<Future<Double>> futures = new ArrayList<>();
            for (int k = 0; k < nStarts; k++) {
                final int kk = k;
                futures.add(executor.submit(() -> {
                    Minimizer.Result result = ccds[kk].minimizeFrom(startxs[kk]);
                    return result.energy;
                }));
            }
            double bestE = Double.POSITIVE_INFINITY;
            for (Future<Double> f : futures) {
                try {
                    double e = f.get();
                    if (e < bestE) bestE = e;
                } catch (Exception ex) {
                    throw new RuntimeException("Multi warm start CCD failed", ex);
                }
            }
            mwarmEnergies[i] = bestE;
            totalMWCCDMinNs += System.nanoTime() - t1;

            // Cleanup
            for (int k = 0; k < nStarts; k++) {
                efuncs[k].close();
                ccds[k].close();
            }
        }
        executor.shutdown();
        long mwarmTotalNs = System.nanoTime() - mwarmTotalStart;

        // 6. Print comparison table
        System.out.println(String.format("\n%-5s| %11s | %11s | %8s | %8s | %8s",
                "Conf", "CCD_E", "MWarm_E", "CCD_drop", "MW_rec%", "CCD_ms"));
        System.out.println("-".repeat(70));

        double sumGapMW = 0;
        double sumRecMW = 0;
        double sumCCDMs = 0;

        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double mwDrop = rigidEnergies[i] - mwarmEnergies[i];
            double mwRec = (ccdDrop > 0) ? (mwDrop / ccdDrop * 100.0) : 0;
            double ccdMs = ccdTimesNs[i] / 1e6;

            sumGapMW += mwarmEnergies[i] - ccdEnergies[i];
            sumRecMW += mwRec;
            sumCCDMs += ccdMs;

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %7.1f%% | %8.2f",
                    i, ccdEnergies[i], mwarmEnergies[i],
                    ccdDrop, mwRec, ccdMs));
        }
        System.out.println("-".repeat(70));

        System.out.println(String.format("\nEnergy Summary:"));
        System.out.println(String.format("  Avg gap MultiWarm-CCD:    %.4f kcal/mol", sumGapMW / n));
        System.out.println(String.format("  Avg recovery MultiWarm:   %.1f%%", sumRecMW / n));

        System.out.println(String.format("\nTiming Summary (%d confs, %d cores, top-%d starts, CCD iter=%d):",
                n, numCores, multiN, ccdIterForHybrid));
        System.out.println(String.format("  CCD total:          %10.1f ms", ccdTotalNs / 1e6));
        System.out.println(String.format("  MultiWarm total:    %10.1f ms  (speedup vs CCD: %.2fx)",
                mwarmTotalNs / 1e6, (double) ccdTotalNs / mwarmTotalNs));
        System.out.println(String.format("\n  MultiWarm Phase Breakdown:"));
        System.out.println(String.format("    GridDP minimize:  %10.1f ms  (%.1f%%)",
                totalMWGridDPNs / 1e6, 100.0 * totalMWGridDPNs / mwarmTotalNs));
        System.out.println(String.format("    Top-N init:       %10.1f ms  (%.1f%%)",
                totalMWInitNs / 1e6, 100.0 * totalMWInitNs / mwarmTotalNs));
        System.out.println(String.format("    Top-N pop (%dx): %10.1f ms  (%.1f%%)",
                multiN, totalMWPopNs / 1e6, 100.0 * totalMWPopNs / mwarmTotalNs));
        System.out.println(String.format("    CCD obj create:   %10.1f ms  (%.1f%%)",
                totalMWCCDCreateNs / 1e6, 100.0 * totalMWCCDCreateNs / mwarmTotalNs));
        System.out.println(String.format("    CCD minimize:     %10.1f ms  (%.1f%%)",
                totalMWCCDMinNs / 1e6, 100.0 * totalMWCCDMinNs / mwarmTotalNs));

        long wallEnd = System.nanoTime();
        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (wallEnd - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    private GridDPMinimizer makeGridDPMinimizer(
            SimpleConfSpace confSpace,
            InteractionGraph ig,
            RootedTreeEdge rootEdge,
            ForcefieldParams ffparams,
            SimpleReferenceEnergies eref,
            SurrogateRunConfig surrogateCfg,
            ConfEnergyCalculator surrogateConfEcalc,
            MLPSurrogateMatrix mlpSurrogate
    ) {
        if (!surrogateCfg.enabled || surrogateCfg.modelType == SurrogateModelType.BranchResidual) {
            return new GridDPMinimizer(confSpace, ig, rootEdge, 2, ffparams, eref, true);
        }
        if (surrogateCfg.modelType == SurrogateModelType.Quadratic) {
            return new ApproximatedGridDPMinimizer(
                    confSpace, ig, rootEdge, 2, ffparams, eref, true,
                    surrogateConfEcalc.amat, surrogateCfg.gridDPFallbackToForcefield
            );
        }
        if (surrogateCfg.modelType == SurrogateModelType.MLP) {
            if (surrogateCfg.mlpSingleTaskModel) {
                throw new IllegalStateException(
                        "Global single-task MLP is not supported for GridDP decomposition yet; use CCD-only benchmark"
                );
            }
            return new MLPApproximatedGridDPMinimizer(
                    confSpace, ig, rootEdge, 2, ffparams, eref, true,
                    mlpSurrogate, surrogateCfg.gridDPFallbackToForcefield
            );
        }
        throw new IllegalStateException("Unsupported surrogate model type: " + surrogateCfg.modelType);
    }

    private ObjectiveWithCleanup makeObjectiveFunction(
            SimpleConfSpace confSpace,
            RCTuple tuple,
            ConfEnergyCalculator confEcalcMinimized,
            ConfEnergyCalculator surrogateConfEcalc,
            MLPSurrogateMatrix mlpSurrogate,
            SurrogateRunConfig surrogateCfg,
            EnergyCalculator.Type.Context ecalcContext,
            TaskGlobalMLPSurrogate globalMLPSurrogate
    ) {
        if (!surrogateCfg.enabled || surrogateCfg.modelType == SurrogateModelType.BranchResidual) {
            ParametricMolecule pmol = confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
            EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol);
            return new ObjectiveWithCleanup(new MoleculeObjectiveFunction(pmol, efunc), efunc);
        }
        if (surrogateCfg.modelType == SurrogateModelType.Quadratic) {
            return new ObjectiveWithCleanup(
                    ConfSpaceSpecificSurrogateFactory.makeApproximationObjective(
                            surrogateConfEcalc, tuple, surrogateCfg.requireFullApproximation
                    ),
                    null
            );
        }
        if (surrogateCfg.modelType == SurrogateModelType.MLP) {
            if (surrogateCfg.mlpSingleTaskModel) {
                if (globalMLPSurrogate == null) {
                    throw new IllegalStateException("Global MLP surrogate is not available");
                }
                return new ObjectiveWithCleanup(
                        ConfSpaceSpecificSurrogateFactory.makeGlobalMLPSurrogateObjective(
                                confEcalcMinimized, globalMLPSurrogate, tuple
                        ),
                        null
                );
            }
            return new ObjectiveWithCleanup(
                    ConfSpaceSpecificSurrogateFactory.makeMLPSurrogateObjective(
                            confEcalcMinimized, mlpSurrogate, tuple, surrogateCfg.requireFullApproximation
                    ),
                    null
            );
        }
        throw new IllegalStateException("Unsupported surrogate model type: " + surrogateCfg.modelType);
    }

    private double evaluateForcefieldEnergyAt(
            SimpleConfSpace confSpace,
            ConfEnergyCalculator confEcalcMinimized,
            EnergyCalculator.Type.Context ecalcContext,
            RCTuple tuple,
            DoubleMatrix1D x
    ) {
        ParametricMolecule pmol = confSpace.makeMolecule(tuple);
        for (int d = 0; d < x.size(); d++) {
            pmol.dofs.get(d).apply(x.get(d));
        }
        ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
        try (EnergyFunction efunc = ecalcContext.efuncs.make(inters, pmol.mol)) {
            return efunc.getEnergy();
        }
    }

    private double evaluateForcefieldEnergyAt(
            SimpleConfSpace confSpace,
            ConfEnergyCalculator confEcalcMinimized,
            EnergyCalculator.Type.Context ecalcContext,
            RCTuple tuple,
            double[] x
    ) {
        return evaluateForcefieldEnergyAt(
                confSpace, confEcalcMinimized, ecalcContext, tuple, DoubleFactory1D.dense.make(x)
        );
    }

    private double topKRecall(double[] pred, double[] truth, int k) {
        int[] predOrder = argsortAscending(pred);
        int[] truthOrder = argsortAscending(truth);

        Set<Integer> truthTopK = new HashSet<>();
        for (int i = 0; i < k; i++) {
            truthTopK.add(truthOrder[i]);
        }

        int hit = 0;
        for (int i = 0; i < k; i++) {
            if (truthTopK.contains(predOrder[i])) {
                hit++;
            }
        }
        return (double) hit / k;
    }

    private double spearman(double[] x, double[] y) {
        if (x.length != y.length || x.length == 0) {
            return Double.NaN;
        }

        int n = x.length;
        int[] rankX = rankAscending(x);
        int[] rankY = rankAscending(y);

        double mean = (n - 1) / 2.0;
        double cov = 0.0;
        double varX = 0.0;
        double varY = 0.0;
        for (int i = 0; i < n; i++) {
            double dx = rankX[i] - mean;
            double dy = rankY[i] - mean;
            cov += dx * dy;
            varX += dx * dx;
            varY += dy * dy;
        }
        if (varX <= 0 || varY <= 0) {
            return 0.0;
        }
        return cov / Math.sqrt(varX * varY);
    }

    private int[] rankAscending(double[] arr) {
        int[] order = argsortAscending(arr);
        int[] rank = new int[arr.length];
        for (int i = 0; i < order.length; i++) {
            rank[order[i]] = i;
        }
        return rank;
    }

    private int[] argsortAscending(double[] arr) {
        Integer[] idx = new Integer[arr.length];
        for (int i = 0; i < arr.length; i++) {
            idx[i] = i;
        }
        Arrays.sort(idx, Comparator.comparingDouble(i -> arr[i]));
        int[] out = new int[arr.length];
        for (int i = 0; i < arr.length; i++) {
            out[i] = idx[i];
        }
        return out;
    }

    private static class ObjectiveWithCleanup {
        final edu.duke.cs.osprey.minimization.ObjectiveFunction objective;
        final EnergyFunction efunc;

        ObjectiveWithCleanup(edu.duke.cs.osprey.minimization.ObjectiveFunction objective, EnergyFunction efunc) {
            this.objective = objective;
            this.efunc = efunc;
        }

        void close() {
            if (efunc != null) {
                efunc.close();
            }
        }
    }

    private enum SurrogateModelType {
        Quadratic,
        MLP,
        BranchResidual;

        static SurrogateModelType fromString(String s) {
            String normalized = s == null ? "" : s.trim().toLowerCase();
            if (normalized.equals("mlp")) {
                return MLP;
            }
            if (normalized.equals("branchresidual") || normalized.equals("branch-residual")
                    || normalized.equals("branch_ridge") || normalized.equals("branch")) {
                return BranchResidual;
            }
            return Quadratic;
        }
    }

    private static class SurrogateRunConfig {
        final boolean enabled;
        final boolean forceRetrain;
        final File cacheRoot;
        final String taskTag;
        final SurrogateModelType modelType;
        final int samplesPerParam;
        final double errorBudget;
        final boolean requireFullApproximation;
        final boolean gridDPFallbackToForcefield;
        final int mlpHidden1;
        final int mlpHidden2;
        final boolean mlpSingleTaskModel;
        final int mlpMinSamplesPerModel;
        final int mlpMaxSamplesPerModel;
        final int mlpGlobalNumTuples;
        final int mlpGlobalSamplesPerTuple;
        final int mlpEpochs;
        final int mlpBatchSize;
        final double mlpLearningRate;
        final boolean branchResidualEnabled;
        final String branchResidualSamplerMode;
        final int branchResidualPretrainGridDPSamples;
        final int branchResidualFineTuneCCDSamples;
        final double branchResidualLambdaPre;
        final double branchResidualLambdaFine;
        final int branchResidualFeatureHashDim;
        final boolean branchResidualPretrainUseFinalEval;
        final int branchResidualGridSize;
        final int branchResidualGridThreads;
        final int branchResidualLowEnergyPoolLimit;
        final int branchResidualRandomSeed;
        final int branchResidualRidgeMaxIter;
        final double branchResidualRidgeTolerance;

        private SurrogateRunConfig(boolean enabled, boolean forceRetrain, File cacheRoot, String taskTag, SurrogateModelType modelType, int samplesPerParam,
                                   double errorBudget, boolean requireFullApproximation,
                                   boolean gridDPFallbackToForcefield,
                                   int mlpHidden1, int mlpHidden2, boolean mlpSingleTaskModel,
                                   int mlpMinSamplesPerModel, int mlpMaxSamplesPerModel,
                                   int mlpGlobalNumTuples, int mlpGlobalSamplesPerTuple,
                                   int mlpEpochs, int mlpBatchSize, double mlpLearningRate,
                                   boolean branchResidualEnabled, String branchResidualSamplerMode,
                                   int branchResidualPretrainGridDPSamples, int branchResidualFineTuneCCDSamples,
                                   double branchResidualLambdaPre, double branchResidualLambdaFine,
                                   int branchResidualFeatureHashDim, boolean branchResidualPretrainUseFinalEval,
                                   int branchResidualGridSize, int branchResidualGridThreads,
                                   int branchResidualLowEnergyPoolLimit, int branchResidualRandomSeed,
                                   int branchResidualRidgeMaxIter, double branchResidualRidgeTolerance) {
            this.enabled = enabled;
            this.forceRetrain = forceRetrain;
            this.cacheRoot = cacheRoot;
            this.taskTag = taskTag;
            this.modelType = modelType;
            this.samplesPerParam = samplesPerParam;
            this.errorBudget = errorBudget;
            this.requireFullApproximation = requireFullApproximation;
            this.gridDPFallbackToForcefield = gridDPFallbackToForcefield;
            this.mlpHidden1 = mlpHidden1;
            this.mlpHidden2 = mlpHidden2;
            this.mlpSingleTaskModel = mlpSingleTaskModel;
            this.mlpMinSamplesPerModel = mlpMinSamplesPerModel;
            this.mlpMaxSamplesPerModel = mlpMaxSamplesPerModel;
            this.mlpGlobalNumTuples = mlpGlobalNumTuples;
            this.mlpGlobalSamplesPerTuple = mlpGlobalSamplesPerTuple;
            this.mlpEpochs = mlpEpochs;
            this.mlpBatchSize = mlpBatchSize;
            this.mlpLearningRate = mlpLearningRate;
            this.branchResidualEnabled = branchResidualEnabled;
            this.branchResidualSamplerMode = branchResidualSamplerMode;
            this.branchResidualPretrainGridDPSamples = branchResidualPretrainGridDPSamples;
            this.branchResidualFineTuneCCDSamples = branchResidualFineTuneCCDSamples;
            this.branchResidualLambdaPre = branchResidualLambdaPre;
            this.branchResidualLambdaFine = branchResidualLambdaFine;
            this.branchResidualFeatureHashDim = branchResidualFeatureHashDim;
            this.branchResidualPretrainUseFinalEval = branchResidualPretrainUseFinalEval;
            this.branchResidualGridSize = branchResidualGridSize;
            this.branchResidualGridThreads = branchResidualGridThreads;
            this.branchResidualLowEnergyPoolLimit = branchResidualLowEnergyPoolLimit;
            this.branchResidualRandomSeed = branchResidualRandomSeed;
            this.branchResidualRidgeMaxIter = branchResidualRidgeMaxIter;
            this.branchResidualRidgeTolerance = branchResidualRidgeTolerance;
        }

        static SurrogateRunConfig fromSystemProperties() {
            boolean enabled = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "enabled", "false"));
            boolean forceRetrain = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "forceRetrain", "false"));
            String cacheRootPath = System.getProperty(SURROGATE_PROP_PREFIX + "cacheRoot", "/tmp/osprey-confspace-surrogate-cache");
            String taskTag = System.getProperty(SURROGATE_PROP_PREFIX + "taskTag", "benchmark-griddp-vs-ccd");
            SurrogateModelType modelType = SurrogateModelType.fromString(System.getProperty(SURROGATE_PROP_PREFIX + "type", "quadratic"));
            int samplesPerParam = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "samplesPerParam", "10"));
            double errorBudget = Double.parseDouble(System.getProperty(SURROGATE_PROP_PREFIX + "errorBudget", "1000000.0"));
            boolean requireFullApproximation = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "requireFullApproximation", "false"));
            boolean fallback = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "gridDPFallbackToForcefield", "true"));
            int mlpHidden1 = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.hidden1", "32"));
            int mlpHidden2 = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.hidden2", "32"));
            boolean mlpSingleTaskModel = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.singleTaskModel", "false"));
            int mlpMinSamplesPerModel = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.minSamplesPerModel", "64"));
            int mlpMaxSamplesPerModel = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.maxSamplesPerModel", "4096"));
            int mlpGlobalNumTuples = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.globalNumTuples", "256"));
            int mlpGlobalSamplesPerTuple = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.globalSamplesPerTuple", "32"));
            int mlpEpochs = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.epochs", "300"));
            int mlpBatchSize = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.batchSize", "64"));
            double mlpLearningRate = Double.parseDouble(System.getProperty(SURROGATE_PROP_PREFIX + "mlp.learningRate", "0.001"));
            boolean branchResidualEnabled = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.enabled", "false"));
            String branchResidualSamplerMode = System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.samplerMode", "sparse");
            int branchResidualPretrainGridDPSamples = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.pretrain.griddpSamples", "20000"));
            int branchResidualFineTuneCCDSamples = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.finetune.ccdSamples", "2000"));
            double branchResidualLambdaPre = Double.parseDouble(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.ridge.lambdaPre", "0.01"));
            double branchResidualLambdaFine = Double.parseDouble(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.ridge.lambdaFine", "0.1"));
            int branchResidualFeatureHashDim = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.feature.hashDim", "65536"));
            boolean branchResidualPretrainUseFinalEval = Boolean.parseBoolean(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.pretrain.useFinalEval", "true"));
            int branchResidualGridSize = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.gridSize", "2"));
            int branchResidualGridThreads = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.gridThreads", "1"));
            int branchResidualLowEnergyPoolLimit = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.lowEnergyPoolLimit", "4096"));
            int branchResidualRandomSeed = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.randomSeed", "12345"));
            int branchResidualRidgeMaxIter = Integer.parseInt(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.ridge.maxIter", "256"));
            double branchResidualRidgeTolerance = Double.parseDouble(System.getProperty(SURROGATE_PROP_PREFIX + "branchResidual.ridge.tol", "1.0e-6"));

            if (branchResidualEnabled && modelType != SurrogateModelType.BranchResidual) {
                modelType = SurrogateModelType.BranchResidual;
                enabled = true;
            }
            return new SurrogateRunConfig(
                    enabled,
                    forceRetrain,
                    new File(cacheRootPath),
                    taskTag,
                    modelType,
                    samplesPerParam,
                    errorBudget,
                    requireFullApproximation,
                    fallback,
                    mlpHidden1,
                    mlpHidden2,
                    mlpSingleTaskModel,
                    mlpMinSamplesPerModel,
                    mlpMaxSamplesPerModel,
                    mlpGlobalNumTuples,
                    mlpGlobalSamplesPerTuple,
                    mlpEpochs,
                    mlpBatchSize,
                    mlpLearningRate,
                    branchResidualEnabled,
                    branchResidualSamplerMode,
                    branchResidualPretrainGridDPSamples,
                    branchResidualFineTuneCCDSamples,
                    branchResidualLambdaPre,
                    branchResidualLambdaFine,
                    branchResidualFeatureHashDim,
                    branchResidualPretrainUseFinalEval,
                    branchResidualGridSize,
                    branchResidualGridThreads,
                    branchResidualLowEnergyPoolLimit,
                    branchResidualRandomSeed,
                    branchResidualRidgeMaxIter,
                    branchResidualRidgeTolerance
            );
        }
    }

    /**
     * Benchmark: CCD vs GNN on 1CC8 with 8 flexible residues.
     *
     * Loads a pre-trained GNN ONNX model and compares its energy predictions
     * against CCD ground truth. Reports MAE, RMSE, bias, max error, Spearman
     * correlation, top-K recall, and timing.
     *
     * Requires: gnn_data/1CC8_8res/model/gnn_model.onnx (run exportGNNData1CC8_8res + train.py first)
     *
     * System properties:
     *   -Dosprey.gnn.eval.numConfs=100   number of conformations to evaluate
     *   -Dosprey.gnn.eval.topK=10        top-K recall parameter
     *   -Dosprey.gnn.eval.modelPath=...  path to ONNX model (default: gnn_data/1CC8_8res/model/gnn_model.onnx)
     */
    @Test
    public void benchmarkCCDvsGNN1CC8_8res() {
        long wallStart = System.nanoTime();

        // 1. Build conf space: 1CC8, 8 flexible residues (A38-A45)
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        int numConfs = Integer.parseInt(System.getProperty("osprey.gnn.eval.numConfs", "100"));
        int topK = Integer.parseInt(System.getProperty("osprey.gnn.eval.topK", "10"));
        String modelPath = System.getProperty("osprey.gnn.eval.modelPath",
                "gnn_data/1CC8_8res/model/gnn_model.onnx");

        System.out.println("=== CCD vs GNN Benchmark: 1CC8, 8 flexible residues ===");
        System.out.println("Positions: " + confSpace.positions.size());
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }
        System.out.println("numConfs=" + numConfs + ", topK=" + topK);
        System.out.println("Model: " + modelPath);

        // 2. Compute energy matrices
        Parallelism parallelism = Parallelism.makeCpu(1);  // single thread for fair CCD timing
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();

        System.out.println("Computing energy matrix...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();

        // 3. Enumerate conformations via A*
        RCs rcs = new RCs(confSpace);
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional()
                .build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        // 4. Load GNN model
        File onnxFile = new File(modelPath);
        if (!onnxFile.exists()) {
            throw new IllegalStateException("ONNX model not found: " + onnxFile.getAbsolutePath()
                    + "\nRun exportGNNData1CC8_8res + gnn/train.py first.");
        }
        GNNConfEnergyCalculator gnnCalc = new GNNConfEnergyCalculator(
                onnxFile, ematMinimized, confSpace.positions.size());

        // 5. Warm up: one CCD + one GNN call
        if (n > 0) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            gnnCalc.calcEnergy(warmConf);
        }

        // 6. CCD energies
        double[] ccdEnergies = new double[n];
        long ccdStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            RCTuple tuple = new RCTuple(scoredConfs.get(i).getAssignments());
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(tuple).energy;
        }
        long ccdNs = System.nanoTime() - ccdStart;

        // 7. GNN energies
        double[] gnnEnergies = new double[n];
        long gnnStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            gnnEnergies[i] = gnnCalc.calcEnergy(conf);
        }
        long gnnNs = System.nanoTime() - gnnStart;

        // 8. emat baseline energies
        double[] ematEnergies = new double[n];
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            ematEnergies[i] = ematMinimized.confE(conf);
        }

        // 9. Per-conf comparison table
        System.out.println(String.format("%-5s| %11s | %11s | %11s | %8s | %8s",
                "Conf", "CCD_E", "GNN_E", "emat_E", "GNN_err", "emat_err"));
        System.out.println("-".repeat(70));

        double sumAbsGNN = 0, sumSqGNN = 0, sumBiasGNN = 0, maxAbsGNN = 0;
        double sumAbsEmat = 0, sumSqEmat = 0, sumBiasEmat = 0, maxAbsEmat = 0;

        for (int i = 0; i < n; i++) {
            double gnnErr = gnnEnergies[i] - ccdEnergies[i];
            double ematErr = ematEnergies[i] - ccdEnergies[i];
            sumAbsGNN += Math.abs(gnnErr);
            sumSqGNN += gnnErr * gnnErr;
            sumBiasGNN += gnnErr;
            maxAbsGNN = Math.max(maxAbsGNN, Math.abs(gnnErr));
            sumAbsEmat += Math.abs(ematErr);
            sumSqEmat += ematErr * ematErr;
            sumBiasEmat += ematErr;
            maxAbsEmat = Math.max(maxAbsEmat, Math.abs(ematErr));

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %8.4f | %8.4f",
                    i, ccdEnergies[i], gnnEnergies[i], ematEnergies[i], gnnErr, ematErr));
        }
        System.out.println("-".repeat(70));

        // 10. Summary statistics
        double maeGNN = sumAbsGNN / n;
        double rmseGNN = Math.sqrt(sumSqGNN / n);
        double biasGNN = sumBiasGNN / n;
        double maeEmat = sumAbsEmat / n;
        double rmseEmat = Math.sqrt(sumSqEmat / n);
        double biasEmat = sumBiasEmat / n;

        int k = Math.max(1, Math.min(topK, n));
        double gnnRecall5 = topKRecall(gnnEnergies, ccdEnergies, Math.min(5, n));
        double gnnRecall10 = topKRecall(gnnEnergies, ccdEnergies, Math.min(10, n));
        double gnnRecall20 = topKRecall(gnnEnergies, ccdEnergies, Math.min(20, n));
        double gnnRecallK = topKRecall(gnnEnergies, ccdEnergies, k);
        double ematRecallK = topKRecall(ematEnergies, ccdEnergies, k);
        double gnnSpearman = spearman(gnnEnergies, ccdEnergies);
        double ematSpearman = spearman(ematEnergies, ccdEnergies);

        System.out.println("\n======================================================================");
        System.out.println("GNN vs CCD  (N=" + n + ")");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE:      %.6f kcal/mol", maeGNN));
        System.out.println(String.format("  RMSE:     %.6f kcal/mol", rmseGNN));
        System.out.println(String.format("  Bias:     %.6f kcal/mol", biasGNN));
        System.out.println(String.format("  Max|d|:   %.6f kcal/mol", maxAbsGNN));
        System.out.println(String.format("  Spearman: %.4f", gnnSpearman));
        System.out.println(String.format("  Top-5  recall: %.1f%%", gnnRecall5 * 100.0));
        System.out.println(String.format("  Top-10 recall: %.1f%%", gnnRecall10 * 100.0));
        System.out.println(String.format("  Top-20 recall: %.1f%%", gnnRecall20 * 100.0));
        System.out.println(String.format("  Top-%d recall: %.1f%%", k, gnnRecallK * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("emat vs CCD  (baseline, N=" + n + ")");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE:      %.6f kcal/mol", maeEmat));
        System.out.println(String.format("  RMSE:     %.6f kcal/mol", rmseEmat));
        System.out.println(String.format("  Bias:     %.6f kcal/mol", biasEmat));
        System.out.println(String.format("  Max|d|:   %.6f kcal/mol", maxAbsEmat));
        System.out.println(String.format("  Spearman: %.4f", ematSpearman));
        System.out.println(String.format("  Top-%d recall: %.1f%%", k, ematRecallK * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("Improvement: GNN over emat baseline");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE  reduction: %.1f%%", (1.0 - maeGNN / maeEmat) * 100.0));
        System.out.println(String.format("  RMSE reduction: %.1f%%", (1.0 - rmseGNN / rmseEmat) * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("Timing");
        System.out.println("======================================================================");
        System.out.println(String.format("  CCD total:  %10.1f ms  (%.2f ms/conf)", ccdNs / 1e6, ccdNs / 1e6 / n));
        System.out.println(String.format("  GNN total:  %10.1f ms  (%.2f ms/conf)", gnnNs / 1e6, gnnNs / 1e6 / n));
        System.out.println(String.format("  Speedup:    %.1fx", (double) ccdNs / gnnNs));

        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (System.nanoTime() - wallStart) / 1e9));

        gnnCalc.close();
        minimizingEcalc.close();
    }

    /**
     * Export GNN training data for 1CC8 with 8 flexible residues.
     * Parallel CCD generates ground-truth energies; also reports emat vs CCD accuracy baseline.
     */
    @Test
    public void exportGNNData1CC8_8res() {
        long wallStart = System.nanoTime();

        // 1CC8 with 8 flexible residues (A38-A45)
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== GNN Data Export: 1CC8, 8 flexible residues ===");
        System.out.println("Positions: " + confSpace.positions.size());
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }

        int numCPUs = Runtime.getRuntime().availableProcessors();
        System.out.println("Available CPUs: " + numCPUs);
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);

        int numSamples = Integer.parseInt(System.getProperty(
                "osprey.gnn.numSamples", "10000"));
        File outputDir = new File(System.getProperty(
                "osprey.gnn.outputDir", "gnn_data/1CC8_8res"));

        System.out.println("Enumerating " + numSamples + " conformations via A*...");
        System.out.println("Output dir: " + outputDir.getAbsolutePath());

        try {
            GNNDataExporter exporter = new GNNDataExporter(
                    confEcalcMinimized, ematMinimized, ig, rcs);
            exporter.export(numSamples, outputDir);
        } catch (Exception e) {
            throw new RuntimeException("GNN data export failed", e);
        }

        System.out.println(String.format("Total wall-clock time: %.2f s",
                (System.nanoTime() - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    /**
     * Export GNN training data for 1CC8 with 25 flexible residues (A26-A50).
     */
    @Test
    public void exportGNNData1CC8_25res() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        for (int r = 26; r <= 50; r++) {
            strand.flexibility.get("A" + r).setLibraryRotamers(Strand.WildType).setContinuous();
        }
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== GNN Data Export: 1CC8, 25 flexible residues (A26-A50) ===");
        System.out.println("Positions: " + confSpace.positions.size());
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }

        int numCPUs = Runtime.getRuntime().availableProcessors();
        System.out.println("Available CPUs: " + numCPUs);
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("Computing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, DIST_CUTOFF, ENERGY_CUTOFF);

        int numSamples = Integer.parseInt(System.getProperty(
                "osprey.gnn.numSamples", "50000"));
        File outputDir = new File(System.getProperty(
                "osprey.gnn.outputDir", "gnn_data/1CC8_25res"));

        System.out.println("Enumerating " + numSamples + " conformations via A*...");
        System.out.println("Output dir: " + outputDir.getAbsolutePath());

        try {
            GNNDataExporter exporter = new GNNDataExporter(
                    confEcalcMinimized, ematMinimized, ig, rcs);
            exporter.export(numSamples, outputDir);
        } catch (Exception e) {
            throw new RuntimeException("GNN data export failed", e);
        }

        System.out.println(String.format("Total wall-clock time: %.2f s",
                (System.nanoTime() - wallStart) / 1e9));

        minimizingEcalc.close();
        rigidEcalc.close();
    }

    /**
     * Benchmark: CCD vs GNN on 1CC8 with 25 flexible residues (A26-A50).
     */
    @Test
    public void benchmarkCCDvsGNN1CC8_25res() {
        long wallStart = System.nanoTime();

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        for (int r = 26; r <= 50; r++) {
            strand.flexibility.get("A" + r).setLibraryRotamers(Strand.WildType).setContinuous();
        }
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        int numConfs = Integer.parseInt(System.getProperty("osprey.gnn.eval.numConfs", "200"));
        int topK = Integer.parseInt(System.getProperty("osprey.gnn.eval.topK", "10"));
        String modelPath = System.getProperty("osprey.gnn.eval.modelPath",
                "gnn_data/1CC8_25res/model/gnn_model.onnx");

        System.out.println("=== CCD vs GNN Benchmark: 1CC8, 25 flexible residues (A26-A50) ===");
        System.out.println("Positions: " + confSpace.positions.size());
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }
        System.out.println("numConfs=" + numConfs + ", topK=" + topK);
        System.out.println("Model: " + modelPath);

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build()
                .calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref)
                .build();

        System.out.println("Computing energy matrix...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional()
                .build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("Got " + n + " conformations.\n");

        File onnxFile = new File(modelPath);
        if (!onnxFile.exists()) {
            throw new IllegalStateException("ONNX model not found: " + onnxFile.getAbsolutePath()
                    + "\nRun exportGNNData1CC8_25res + gnn/train.py first.");
        }
        GNNConfEnergyCalculator gnnCalc = new GNNConfEnergyCalculator(
                onnxFile, ematMinimized, confSpace.positions.size());

        if (n > 0) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            gnnCalc.calcEnergy(warmConf);
        }

        double[] ccdEnergies = new double[n];
        long ccdStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            RCTuple tuple = new RCTuple(scoredConfs.get(i).getAssignments());
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(tuple).energy;
        }
        long ccdNs = System.nanoTime() - ccdStart;

        double[] gnnEnergies = new double[n];
        long gnnStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            gnnEnergies[i] = gnnCalc.calcEnergy(conf);
        }
        long gnnNs = System.nanoTime() - gnnStart;

        double[] ematEnergies = new double[n];
        for (int i = 0; i < n; i++) {
            int[] conf = scoredConfs.get(i).getAssignments();
            ematEnergies[i] = ematMinimized.confE(conf);
        }

        System.out.println(String.format("%-5s| %11s | %11s | %11s | %8s | %8s",
                "Conf", "CCD_E", "GNN_E", "emat_E", "GNN_err", "emat_err"));
        System.out.println("-".repeat(70));

        double sumAbsGNN = 0, sumSqGNN = 0, sumBiasGNN = 0, maxAbsGNN = 0;
        double sumAbsEmat = 0, sumSqEmat = 0, sumBiasEmat = 0, maxAbsEmat = 0;

        for (int i = 0; i < n; i++) {
            double gnnErr = gnnEnergies[i] - ccdEnergies[i];
            double ematErr = ematEnergies[i] - ccdEnergies[i];
            sumAbsGNN += Math.abs(gnnErr);
            sumSqGNN += gnnErr * gnnErr;
            sumBiasGNN += gnnErr;
            maxAbsGNN = Math.max(maxAbsGNN, Math.abs(gnnErr));
            sumAbsEmat += Math.abs(ematErr);
            sumSqEmat += ematErr * ematErr;
            sumBiasEmat += ematErr;
            maxAbsEmat = Math.max(maxAbsEmat, Math.abs(ematErr));

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %11.4f | %8.4f | %8.4f",
                    i, ccdEnergies[i], gnnEnergies[i], ematEnergies[i], gnnErr, ematErr));
        }
        System.out.println("-".repeat(70));

        double maeGNN = sumAbsGNN / n;
        double rmseGNN = Math.sqrt(sumSqGNN / n);
        double biasGNN = sumBiasGNN / n;
        double maeEmat = sumAbsEmat / n;
        double rmseEmat = Math.sqrt(sumSqEmat / n);
        double biasEmat = sumBiasEmat / n;

        int k = Math.max(1, Math.min(topK, n));
        double gnnRecall5 = topKRecall(gnnEnergies, ccdEnergies, Math.min(5, n));
        double gnnRecall10 = topKRecall(gnnEnergies, ccdEnergies, Math.min(10, n));
        double gnnRecall20 = topKRecall(gnnEnergies, ccdEnergies, Math.min(20, n));
        double gnnRecallK = topKRecall(gnnEnergies, ccdEnergies, k);
        double ematRecallK = topKRecall(ematEnergies, ccdEnergies, k);
        double gnnSpearman = spearman(gnnEnergies, ccdEnergies);
        double ematSpearman = spearman(ematEnergies, ccdEnergies);

        System.out.println("\n======================================================================");
        System.out.println("GNN vs CCD  (N=" + n + ", 25 residues)");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE:      %.6f kcal/mol", maeGNN));
        System.out.println(String.format("  RMSE:     %.6f kcal/mol", rmseGNN));
        System.out.println(String.format("  Bias:     %.6f kcal/mol", biasGNN));
        System.out.println(String.format("  Max|d|:   %.6f kcal/mol", maxAbsGNN));
        System.out.println(String.format("  Spearman: %.4f", gnnSpearman));
        System.out.println(String.format("  Top-5  recall: %.1f%%", gnnRecall5 * 100.0));
        System.out.println(String.format("  Top-10 recall: %.1f%%", gnnRecall10 * 100.0));
        System.out.println(String.format("  Top-20 recall: %.1f%%", gnnRecall20 * 100.0));
        System.out.println(String.format("  Top-%d recall: %.1f%%", k, gnnRecallK * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("emat vs CCD  (baseline, N=" + n + ", 25 residues)");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE:      %.6f kcal/mol", maeEmat));
        System.out.println(String.format("  RMSE:     %.6f kcal/mol", rmseEmat));
        System.out.println(String.format("  Bias:     %.6f kcal/mol", biasEmat));
        System.out.println(String.format("  Max|d|:   %.6f kcal/mol", maxAbsEmat));
        System.out.println(String.format("  Spearman: %.4f", ematSpearman));
        System.out.println(String.format("  Top-%d recall: %.1f%%", k, ematRecallK * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("Improvement: GNN over emat baseline");
        System.out.println("======================================================================");
        System.out.println(String.format("  MAE  reduction: %.1f%%", (1.0 - maeGNN / maeEmat) * 100.0));
        System.out.println(String.format("  RMSE reduction: %.1f%%", (1.0 - rmseGNN / rmseEmat) * 100.0));

        System.out.println("\n======================================================================");
        System.out.println("Timing");
        System.out.println("======================================================================");
        System.out.println(String.format("  CCD total:  %10.1f ms  (%.2f ms/conf)", ccdNs / 1e6, ccdNs / 1e6 / n));
        System.out.println(String.format("  GNN total:  %10.1f ms  (%.2f ms/conf)", gnnNs / 1e6, gnnNs / 1e6 / n));
        System.out.println(String.format("  Speedup:    %.1fx", (double) ccdNs / gnnNs));

        System.out.println(String.format("\n  Total wall-clock time: %.2f s", (System.nanoTime() - wallStart) / 1e9));

        gnnCalc.close();
        minimizingEcalc.close();
    }

    private int countEdges(InteractionGraph ig) {
        int count = 0;
        for (int i = 0; i < ig.getNumPositions(); i++) {
            for (int j = i + 1; j < ig.getNumPositions(); j++) {
                if (ig.hasEdge(i, j)) count++;
            }
        }
        return count;
    }
}
