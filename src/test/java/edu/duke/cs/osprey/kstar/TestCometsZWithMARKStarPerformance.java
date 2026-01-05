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

package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.tools.Log.log;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.GradientDescentMARKStarPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Performance comparison test between CometsZ using GradientDescentPfunc
 * and CometsZ using MARK* (GradientDescentMARKStarPfunc).
 *
 * This test suite compares the runtime performance to evaluate the speedup
 * achieved by using MARK* for partition function calculations.
 */
public class TestCometsZWithMARKStarPerformance {

	private static ForcefieldParams ffparams = new ForcefieldParams();

	/**
	 * Helper class to store test results and timing information
	 */
	private static class PerformanceResult {
		String algorithm;
		long setupTimeMs;
		long executionTimeMs;
		long totalTimeMs;
		int numSequencesFound;
		List<String> topSequences;

		PerformanceResult(String algorithm) {
			this.algorithm = algorithm;
			this.topSequences = new ArrayList<>();
		}

		@Override
		public String toString() {
			return String.format(
				"%s Results:\n" +
				"  Setup Time:     %,d ms\n" +
				"  Execution Time: %,d ms\n" +
				"  Total Time:     %,d ms\n" +
				"  Sequences:      %d\n" +
				"  Top Sequences:  %s",
				algorithm, setupTimeMs, executionTimeMs, totalTimeMs,
				numSequencesFound, topSequences
			);
		}
	}

	/**
	 * Initialize CometsZ states with GradientDescentPfunc
	 */
	private static void initCometsZStatesWithGradientDescent(
		List<CometsZ.State> states,
		ForcefieldParams ffparams,
		boolean boundedMemory,
		Map<CometsZ.State, EnergyMatrix> rigidEmats,
		Map<CometsZ.State, EnergyMatrix> minimizingEmats
	) {
		List<SimpleConfSpace> confSpaceList = states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {
			for (CometsZ.State state : states) {
				// how should we define energies of conformations?
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// calc minimizing energy matrix
				EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = minimizingEmat;
				minimizingEmats.put(state, minimizingEmat);

				// calc rigid energy matrix (without minimization)
				// FIX: Use confEcalc with reference energies, then disable minimization
				ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(state.confEcalc.eref)  // Reuse same reference energies
					.build();
				EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
					.build()
					.calcEnergyMatrix();
				rigidEmats.put(state, rigidEmat);

				// how should confs be ordered and searched?
				state.confTreeFactory = (rcs) -> new edu.duke.cs.osprey.astar.conf.ConfAStarTree.Builder(minimizingEmat, rcs)
					.setMaxNumNodes(boundedMemory ? 100000L : null)
					.setTraditional()
					.build();

				// Use GradientDescentPfunc for partition function calculation
				state.pfuncFactory = (rcs) -> new GradientDescentPfunc(
					state.confEcalc,
					state.confTreeFactory.apply(rcs),
					state.confTreeFactory.apply(rcs),
					rcs.getNumConformations()
				);
			}
		}
	}

	/**
	 * Initialize CometsZ states with MARK* (MARKStarBound)
	 * Following the pattern from TestMSKStar.initStates()
	 */
	private static void initCometsZStatesWithMARKStar(
		List<CometsZ.State> states,
		ForcefieldParams ffparams,
		boolean boundedMemory,
		Map<CometsZ.State, EnergyMatrix> rigidEmats,
		Map<CometsZ.State, EnergyMatrix> minimizingEmats
	) {
		List<SimpleConfSpace> confSpaceList = states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());

		// Create minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.setIsMinimizing(true)
			.build();

		// Create rigid energy calculator using SharedBuilder
		EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
			.setIsMinimizing(false)
			.build();

		try {
			for (CometsZ.State state : states) {
				// Setup minimizing conf energy calculator
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, minimizingEcalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// Calculate minimizing energy matrix
				EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = minimizingEmat;
				minimizingEmats.put(state, minimizingEmat);

				// Calculate rigid energy matrix with same reference energies
				ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(state.confSpace, rigidEcalc)
					.setReferenceEnergies(state.confEcalc.eref)  // Reuse same reference energies
					.build();
				EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
					.build()
					.calcEnergyMatrix();
				rigidEmats.put(state, rigidEmat);

				// how should confs be ordered and searched?
				state.confTreeFactory = (rcs) -> new edu.duke.cs.osprey.astar.conf.ConfAStarTree.Builder(minimizingEmat, rcs)
					.setMaxNumNodes(boundedMemory ? 100000L : null)
					.setTraditional()
					.build();

				// Use MARK* (MARKStarBound) for partition function calculation
				// Following TestMSKStar pattern: create MARKStarBound in the factory
				state.pfuncFactory = (rcs) -> {
					edu.duke.cs.osprey.markstar.framework.MARKStarBound markstarPfunc =
						new edu.duke.cs.osprey.markstar.framework.MARKStarBound(
							state.confSpace,
							rigidEmats.get(state),
							minimizingEmats.get(state),
							state.confEcalc,
							rcs,
							Parallelism.makeCpu(8)
						);
					// CRITICAL: Set the correction matrix to avoid NullPointerException
					markstarPfunc.setCorrections(new edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix(
						state.confSpace,
						minimizingEmats.get(state),
						state.confEcalc
					));
					return markstarPfunc;
				};
			}
		} finally {
			// Clean up energy calculators
			minimizingEcalc.close();
			rigidEcalc.close();
		}
	}

	/**
	 * Create CometsZ instance for full 2RL0 test case (8 positions) with GradientDescentPfunc
	 */
	private static CometsZ makeCometsZ2RL0WithGradientDescent(boolean boundedMemory) {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithGradientDescent(cometsZ.states, ffparams, boundedMemory, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Create CometsZ instance for full 2RL0 test case (8 positions) with MARK*
	 */
	private static CometsZ makeCometsZ2RL0WithMARKStar(boolean boundedMemory) {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithMARKStar(cometsZ.states, ffparams, boundedMemory, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Create CometsZ instance for 2RL0 with only one mutant - GradientDescent version
	 */
	private static CometsZ makeCometsZ2RL0OnlyOneMutantWithGradientDescent() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0OnlyOneMutant();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithGradientDescent(cometsZ.states, ffparams, false, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Create CometsZ instance for 2RL0 with only one mutant - MARK* version
	 */
	private static CometsZ makeCometsZ2RL0OnlyOneMutantWithMARKStar() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0OnlyOneMutant();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithMARKStar(cometsZ.states, ffparams, false, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Create CometsZ instance for 2RL0 space without wild type - GradientDescent version
	 */
	private static CometsZ makeCometsZ2RL0SpaceWithoutWildTypeWithGradientDescent() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0SpaceWithoutWildType();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(2)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithGradientDescent(cometsZ.states, ffparams, false, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Create CometsZ instance for 2RL0 space without wild type - MARK* version
	 */
	private static CometsZ makeCometsZ2RL0SpaceWithoutWildTypeWithMARKStar() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0SpaceWithoutWildType();
		final double epsilon = 0.95;

		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(2)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.build();

		Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
		Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
		initCometsZStatesWithMARKStar(cometsZ.states, ffparams, false, rigidEmats, minimizingEmats);

		return cometsZ;
	}

	/**
	 * Prepare CometsZ states and run block while EnergyCalculator is active
	 */
	private static void prepCometsZStates(CometsZ cometsZ, ForcefieldParams ffparams, Runnable block) {
		List<SimpleConfSpace> confSpaces = new ArrayList<>();
		for (CometsZ.State state : cometsZ.states) {
			confSpaces.add(state.confSpace);
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			// refresh the conf ecalcs
			for (CometsZ.State state : cometsZ.states) {
				state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
			}

			block.run();
		}
	}

	/**
	 * Run performance comparison on full 2RL0 test (8 positions, 25 sequences)
	 */
	@Test
	public void compare2RL0FullGradientDescentVsMARKStar() {
		log("\n========================================");
		log("Performance Comparison: CometsZ Full 2RL0 Test");
		log("GradientDescentPfunc vs MARK* (GradientDescentMARKStarPfunc)");
		log("========================================\n");

		int numSequences = 25;
		boolean boundedMemory = false;

		// Test CometsZ with GradientDescentPfunc
		PerformanceResult gradientDescentResult = new PerformanceResult("CometsZ (GradientDescent)");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZGradient = makeCometsZ2RL0WithGradientDescent(boundedMemory);
		gradientDescentResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZGradient, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZGradient.findBestSequences(numSequences);
			gradientDescentResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			gradientDescentResult.totalTimeMs = gradientDescentResult.setupTimeMs + gradientDescentResult.executionTimeMs;
			gradientDescentResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				gradientDescentResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test CometsZ with MARK*
		PerformanceResult markstarResult = new PerformanceResult("CometsZ (MARK*)");
		long startSetup2 = System.currentTimeMillis();
		CometsZ cometsZMarkstar = makeCometsZ2RL0WithMARKStar(boundedMemory);
		markstarResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsZStates(cometsZMarkstar, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZMarkstar.findBestSequences(numSequences);
			markstarResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			markstarResult.totalTimeMs = markstarResult.setupTimeMs + markstarResult.executionTimeMs;
			markstarResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				markstarResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + gradientDescentResult.toString());
		log("\n" + markstarResult.toString());
		log("\n========================================");
		log(String.format("MARK* Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) gradientDescentResult.totalTimeMs / markstarResult.totalTimeMs,
			(double) gradientDescentResult.setupTimeMs / markstarResult.setupTimeMs,
			(double) gradientDescentResult.executionTimeMs / markstarResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on 2RL0 with only one mutant
	 */
	@Test
	public void compare2RL0OnlyOneMutantGradientDescentVsMARKStar() {
		log("\n========================================");
		log("Performance Comparison: CometsZ 2RL0 Only One Mutant");
		log("GradientDescentPfunc vs MARK* (GradientDescentMARKStarPfunc)");
		log("========================================\n");

		int numSequences = 1;

		// Test CometsZ with GradientDescentPfunc
		PerformanceResult gradientDescentResult = new PerformanceResult("CometsZ (GradientDescent)");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZGradient = makeCometsZ2RL0OnlyOneMutantWithGradientDescent();
		gradientDescentResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZGradient, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZGradient.findBestSequences(numSequences);
			gradientDescentResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			gradientDescentResult.totalTimeMs = gradientDescentResult.setupTimeMs + gradientDescentResult.executionTimeMs;
			gradientDescentResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				gradientDescentResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test CometsZ with MARK*
		PerformanceResult markstarResult = new PerformanceResult("CometsZ (MARK*)");
		long startSetup2 = System.currentTimeMillis();
		CometsZ cometsZMarkstar = makeCometsZ2RL0OnlyOneMutantWithMARKStar();
		markstarResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsZStates(cometsZMarkstar, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZMarkstar.findBestSequences(numSequences);
			markstarResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			markstarResult.totalTimeMs = markstarResult.setupTimeMs + markstarResult.executionTimeMs;
			markstarResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				markstarResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + gradientDescentResult.toString());
		log("\n" + markstarResult.toString());
		log("\n========================================");
		log(String.format("MARK* Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) gradientDescentResult.totalTimeMs / markstarResult.totalTimeMs,
			(double) gradientDescentResult.setupTimeMs / markstarResult.setupTimeMs,
			(double) gradientDescentResult.executionTimeMs / markstarResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on 2RL0 space without wild type
	 */
	@Test
	public void compare2RL0SpaceWithoutWildTypeGradientDescentVsMARKStar() {
		log("\n========================================");
		log("Performance Comparison: CometsZ 2RL0 Space Without Wild Type");
		log("GradientDescentPfunc vs MARK* (GradientDescentMARKStarPfunc)");
		log("========================================\n");

		int numSequences = 2;

		// Test CometsZ with GradientDescentPfunc
		PerformanceResult gradientDescentResult = new PerformanceResult("CometsZ (GradientDescent)");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZGradient = makeCometsZ2RL0SpaceWithoutWildTypeWithGradientDescent();
		gradientDescentResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZGradient, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZGradient.findBestSequences(numSequences);
			gradientDescentResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			gradientDescentResult.totalTimeMs = gradientDescentResult.setupTimeMs + gradientDescentResult.executionTimeMs;
			gradientDescentResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				gradientDescentResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test CometsZ with MARK*
		PerformanceResult markstarResult = new PerformanceResult("CometsZ (MARK*)");
		long startSetup2 = System.currentTimeMillis();
		CometsZ cometsZMarkstar = makeCometsZ2RL0SpaceWithoutWildTypeWithMARKStar();
		markstarResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsZStates(cometsZMarkstar, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> seqs = cometsZMarkstar.findBestSequences(numSequences);
			markstarResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			markstarResult.totalTimeMs = markstarResult.setupTimeMs + markstarResult.executionTimeMs;
			markstarResult.numSequencesFound = seqs.size();
			for (int i = 0; i < Math.min(3, seqs.size()); i++) {
				markstarResult.topSequences.add(seqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + gradientDescentResult.toString());
		log("\n" + markstarResult.toString());
		log("\n========================================");
		log(String.format("MARK* Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) gradientDescentResult.totalTimeMs / markstarResult.totalTimeMs,
			(double) gradientDescentResult.setupTimeMs / markstarResult.setupTimeMs,
			(double) gradientDescentResult.executionTimeMs / markstarResult.executionTimeMs
		));
		log("========================================\n");
	}
}
