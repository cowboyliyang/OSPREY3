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

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.Comets;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Test;

import java.util.*;

/**
 * Performance comparison test between CometsZ and COMETS algorithms.
 *
 * This test suite compares the runtime performance of CometsZ and COMETS
 * on identical protein design problems to evaluate their relative efficiency.
 */
public class TestCometsVsCometsZPerformance {

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
	 * Initialize states for CometsZ with energy matrices and partition functions
	 */
	private static void initCometsZStates(List<CometsZ.State> states, ForcefieldParams ffparams, boolean boundedMemory) {

		// Create a single shared EnergyCalculator for all states
		List<SimpleConfSpace> confSpaceList = new ArrayList<>();
		for (CometsZ.State state : states) {
			confSpaceList.add(state.confSpace);
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			for (CometsZ.State state : states) {

				// how should we define energies of conformations?
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// calc energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				// how should confs be ordered and searched?
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setMaxNumNodes(boundedMemory ? 100000L : null)
					.setTraditional()
					.build();

				// how should partition functions be calculated?
				state.pfuncFactory = (rcs) -> new GradientDescentPfunc(
					state.confEcalc,
					new ConfAStarTree.Builder(emat, rcs)
						.setMaxNumNodes(boundedMemory ? 100000L : null)
						.setTraditional()
						.build(),
					new ConfAStarTree.Builder(emat, rcs)
						.setMaxNumNodes(boundedMemory ? 100000L : null)
						.setTraditional()
						.build(),
					rcs.getNumConformations()
				);
			}
		}
	}

	/**
	 * Initialize states for COMETS with energy matrices and pruning
	 */
	private static void initCometsStates(List<Comets.State> states, boolean boundedMemory) {

		// Create a single shared EnergyCalculator for all states
		List<SimpleConfSpace> confSpaceList = new ArrayList<>();
		for (Comets.State state : states) {
			confSpaceList.add(state.confSpace);
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			for (Comets.State state : states) {

				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();

				state.fragmentEnergies = emat;

				// configure A*
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();
			}
		}
	}

	/**
	 * Create CometsZ instance for 2RL0 test case (small)
	 */
	private static CometsZ makeCometsZ2RL0Small(boolean boundedMemory) {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the CometsZ states
		CometsZ.State complex = new CometsZ.State(
			"Complex",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		CometsZ.State unbound = new CometsZ.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);

		CometsZ.State ligandOnly = new CometsZ.State(
			"Ligand",
			new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build()
		);

		// configure CometsZ
		CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(unbound, -1.0)
			.addState(ligandOnly, -1.0)
			.build();

		CometsZ cometsZ = new CometsZ.Builder(objective)
			.setEpsilon(0.95)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0)
			.setObjectiveWindowMax(100.0)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initCometsZStates(cometsZ.states, ffparams, boundedMemory);

		return cometsZ;
	}

	/**
	 * Create COMETS instance for 2RL0 test case (small)
	 */
	private static Comets makeComets2RL0Small(boolean boundedMemory) {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states
		Comets.State bound = new Comets.State(
			"Bound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		Comets.State unbound = new Comets.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);

		// configure COMETS
		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();

		Comets comets = new Comets.Builder(objective)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initCometsStates(comets.states, boundedMemory);

		return comets;
	}

	/**
	 * Create CometsZ instance for full 2RL0 test case (8 positions)
	 */
	private static CometsZ makeCometsZ2RL0Full(boolean boundedMemory) {
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

		initCometsZStates(cometsZ.states, ffparams, boundedMemory);

		return cometsZ;
	}

	/**
	 * Create COMETS instance for full 2RL0 test case (8 positions)
	 */
	private static Comets makeComets2RL0Full(boolean boundedMemory) {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		Comets.State bound = new Comets.State("Bound", confSpaces.complex);
		Comets.State unbound = new Comets.State("Unbound", confSpaces.protein);

		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();

		Comets comets = new Comets.Builder(objective)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initCometsStates(comets.states, boundedMemory);

		return comets;
	}

	/**
	 * Create CometsZ instance for 2RL0 with only one mutant
	 */
	private static CometsZ makeCometsZ2RL0OnlyOneMutant() {
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

		initCometsZStates(cometsZ.states, ffparams, false);

		return cometsZ;
	}

	/**
	 * Create COMETS instance for 2RL0 with only one mutant
	 */
	private static Comets makeComets2RL0OnlyOneMutant() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0OnlyOneMutant();

		Comets.State bound = new Comets.State("Bound", confSpaces.complex);
		Comets.State unbound = new Comets.State("Unbound", confSpaces.protein);

		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();

		Comets comets = new Comets.Builder(objective)
			.build();

		initCometsStates(comets.states, false);

		return comets;
	}

	/**
	 * Create CometsZ instance for 2RL0 space without wild type
	 */
	private static CometsZ makeCometsZ2RL0SpaceWithoutWildType() {
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

		initCometsZStates(cometsZ.states, ffparams, false);

		return cometsZ;
	}

	/**
	 * Create COMETS instance for 2RL0 space without wild type
	 */
	private static Comets makeComets2RL0SpaceWithoutWildType() {
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0SpaceWithoutWildType();

		Comets.State bound = new Comets.State("Bound", confSpaces.complex);
		Comets.State unbound = new Comets.State("Unbound", confSpaces.protein);

		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();

		Comets comets = new Comets.Builder(objective)
			.build();

		initCometsStates(comets.states, false);

		return comets;
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
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			// refresh the conf ecalcs
			for (CometsZ.State state : cometsZ.states) {
				state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
			}

			block.run();
		}
	}

	/**
	 * Prepare COMETS states and run block while EnergyCalculator is active
	 */
	private static void prepCometsStates(Comets comets, Runnable block) {
		List<SimpleConfSpace> confSpaces = new ArrayList<>();
		for (Comets.State state : comets.states) {
			confSpaces.add(state.confSpace);
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			// refresh the conf ecalcs
			for (Comets.State state : comets.states) {
				state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
			}

			block.run();
		}
	}

	/**
	 * Run performance comparison on small test case
	 */
	@Test
	public void compareSmall2RL0() {
		log("\n========================================");
		log("Performance Comparison: Small 2RL0 Test");
		log("========================================\n");

		int numSequences = 10;
		boolean boundedMemory = false;

		// Test CometsZ
		PerformanceResult cometsZResult = new PerformanceResult("CometsZ");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZ = makeCometsZ2RL0Small(boundedMemory);
		cometsZResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZ, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> cometsZSeqs = cometsZ.findBestSequences(numSequences);
			cometsZResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsZResult.totalTimeMs = cometsZResult.setupTimeMs + cometsZResult.executionTimeMs;
			cometsZResult.numSequencesFound = cometsZSeqs.size();
			for (int i = 0; i < Math.min(3, cometsZSeqs.size()); i++) {
				cometsZResult.topSequences.add(cometsZSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test COMETS
		PerformanceResult cometsResult = new PerformanceResult("COMETS");
		long startSetup2 = System.currentTimeMillis();
		Comets comets = makeComets2RL0Small(boundedMemory);
		cometsResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsStates(comets, () -> {
			long startExecution = System.currentTimeMillis();
			List<Comets.SequenceInfo> cometsSeqs = comets.findBestSequences(numSequences);
			cometsResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsResult.totalTimeMs = cometsResult.setupTimeMs + cometsResult.executionTimeMs;
			cometsResult.numSequencesFound = cometsSeqs.size();
			for (int i = 0; i < Math.min(3, cometsSeqs.size()); i++) {
				cometsResult.topSequences.add(cometsSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + cometsZResult.toString());
		log("\n" + cometsResult.toString());
		log("\n========================================");
		log(String.format("Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) cometsResult.totalTimeMs / cometsZResult.totalTimeMs,
			(double) cometsResult.setupTimeMs / cometsZResult.setupTimeMs,
			(double) cometsResult.executionTimeMs / cometsZResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on small test case with bounded memory
	 */
	@Test
	public void compareSmall2RL0BoundedMemory() {
		log("\n========================================");
		log("Performance Comparison: Small 2RL0 Test (Bounded Memory)");
		log("========================================\n");

		int numSequences = 10;
		boolean boundedMemory = true;

		// Test CometsZ
		PerformanceResult cometsZResult = new PerformanceResult("CometsZ");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZ = makeCometsZ2RL0Small(boundedMemory);
		cometsZResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZ, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> cometsZSeqs = cometsZ.findBestSequences(numSequences);
			cometsZResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsZResult.totalTimeMs = cometsZResult.setupTimeMs + cometsZResult.executionTimeMs;
			cometsZResult.numSequencesFound = cometsZSeqs.size();
			for (int i = 0; i < Math.min(3, cometsZSeqs.size()); i++) {
				cometsZResult.topSequences.add(cometsZSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test COMETS
		PerformanceResult cometsResult = new PerformanceResult("COMETS");
		long startSetup2 = System.currentTimeMillis();
		Comets comets = makeComets2RL0Small(boundedMemory);
		cometsResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsStates(comets, () -> {
			long startExecution = System.currentTimeMillis();
			List<Comets.SequenceInfo> cometsSeqs = comets.findBestSequences(numSequences);
			cometsResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsResult.totalTimeMs = cometsResult.setupTimeMs + cometsResult.executionTimeMs;
			cometsResult.numSequencesFound = cometsSeqs.size();
			for (int i = 0; i < Math.min(3, cometsSeqs.size()); i++) {
				cometsResult.topSequences.add(cometsSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + cometsZResult.toString());
		log("\n" + cometsResult.toString());
		log("\n========================================");
		log(String.format("Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) cometsResult.totalTimeMs / cometsZResult.totalTimeMs,
			(double) cometsResult.setupTimeMs / cometsZResult.setupTimeMs,
			(double) cometsResult.executionTimeMs / cometsZResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on full 2RL0 test (8 positions, 25 sequences)
	 */
	@Test
	public void compare2RL0Full() {
		log("\n========================================");
		log("Performance Comparison: Full 2RL0 Test (8 positions, 25 sequences)");
		log("========================================\n");

		int numSequences = 25;
		boolean boundedMemory = false;

		// Test CometsZ
		PerformanceResult cometsZResult = new PerformanceResult("CometsZ");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZ = makeCometsZ2RL0Full(boundedMemory);
		cometsZResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZ, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> cometsZSeqs = cometsZ.findBestSequences(numSequences);
			cometsZResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsZResult.totalTimeMs = cometsZResult.setupTimeMs + cometsZResult.executionTimeMs;
			cometsZResult.numSequencesFound = cometsZSeqs.size();
			for (int i = 0; i < Math.min(3, cometsZSeqs.size()); i++) {
				cometsZResult.topSequences.add(cometsZSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test COMETS
		PerformanceResult cometsResult = new PerformanceResult("COMETS");
		long startSetup2 = System.currentTimeMillis();
		Comets comets = makeComets2RL0Full(boundedMemory);
		cometsResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsStates(comets, () -> {
			long startExecution = System.currentTimeMillis();
			List<Comets.SequenceInfo> cometsSeqs = comets.findBestSequences(numSequences);
			cometsResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsResult.totalTimeMs = cometsResult.setupTimeMs + cometsResult.executionTimeMs;
			cometsResult.numSequencesFound = cometsSeqs.size();
			for (int i = 0; i < Math.min(3, cometsSeqs.size()); i++) {
				cometsResult.topSequences.add(cometsSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + cometsZResult.toString());
		log("\n" + cometsResult.toString());
		log("\n========================================");
		log(String.format("Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) cometsResult.totalTimeMs / cometsZResult.totalTimeMs,
			(double) cometsResult.setupTimeMs / cometsZResult.setupTimeMs,
			(double) cometsResult.executionTimeMs / cometsZResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on 2RL0 with only one mutant
	 */
	@Test
	public void compare2RL0OnlyOneMutant() {
		log("\n========================================");
		log("Performance Comparison: 2RL0 Only One Mutant");
		log("========================================\n");

		int numSequences = 1;

		// Test CometsZ
		PerformanceResult cometsZResult = new PerformanceResult("CometsZ");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZ = makeCometsZ2RL0OnlyOneMutant();
		cometsZResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZ, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> cometsZSeqs = cometsZ.findBestSequences(numSequences);
			cometsZResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsZResult.totalTimeMs = cometsZResult.setupTimeMs + cometsZResult.executionTimeMs;
			cometsZResult.numSequencesFound = cometsZSeqs.size();
			for (int i = 0; i < Math.min(3, cometsZSeqs.size()); i++) {
				cometsZResult.topSequences.add(cometsZSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test COMETS
		PerformanceResult cometsResult = new PerformanceResult("COMETS");
		long startSetup2 = System.currentTimeMillis();
		Comets comets = makeComets2RL0OnlyOneMutant();
		cometsResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsStates(comets, () -> {
			long startExecution = System.currentTimeMillis();
			List<Comets.SequenceInfo> cometsSeqs = comets.findBestSequences(numSequences);
			cometsResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsResult.totalTimeMs = cometsResult.setupTimeMs + cometsResult.executionTimeMs;
			cometsResult.numSequencesFound = cometsSeqs.size();
			for (int i = 0; i < Math.min(3, cometsSeqs.size()); i++) {
				cometsResult.topSequences.add(cometsSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + cometsZResult.toString());
		log("\n" + cometsResult.toString());
		log("\n========================================");
		log(String.format("Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) cometsResult.totalTimeMs / cometsZResult.totalTimeMs,
			(double) cometsResult.setupTimeMs / cometsZResult.setupTimeMs,
			(double) cometsResult.executionTimeMs / cometsZResult.executionTimeMs
		));
		log("========================================\n");
	}

	/**
	 * Run performance comparison on 2RL0 space without wild type
	 */
	@Test
	public void compare2RL0SpaceWithoutWildType() {
		log("\n========================================");
		log("Performance Comparison: 2RL0 Space Without Wild Type");
		log("========================================\n");

		int numSequences = 2;

		// Test CometsZ
		PerformanceResult cometsZResult = new PerformanceResult("CometsZ");
		long startSetup = System.currentTimeMillis();
		CometsZ cometsZ = makeCometsZ2RL0SpaceWithoutWildType();
		cometsZResult.setupTimeMs = System.currentTimeMillis() - startSetup;

		prepCometsZStates(cometsZ, ffparams, () -> {
			long startExecution = System.currentTimeMillis();
			List<CometsZ.SequenceInfo> cometsZSeqs = cometsZ.findBestSequences(numSequences);
			cometsZResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsZResult.totalTimeMs = cometsZResult.setupTimeMs + cometsZResult.executionTimeMs;
			cometsZResult.numSequencesFound = cometsZSeqs.size();
			for (int i = 0; i < Math.min(3, cometsZSeqs.size()); i++) {
				cometsZResult.topSequences.add(cometsZSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Test COMETS
		PerformanceResult cometsResult = new PerformanceResult("COMETS");
		long startSetup2 = System.currentTimeMillis();
		Comets comets = makeComets2RL0SpaceWithoutWildType();
		cometsResult.setupTimeMs = System.currentTimeMillis() - startSetup2;

		prepCometsStates(comets, () -> {
			long startExecution = System.currentTimeMillis();
			List<Comets.SequenceInfo> cometsSeqs = comets.findBestSequences(numSequences);
			cometsResult.executionTimeMs = System.currentTimeMillis() - startExecution;
			cometsResult.totalTimeMs = cometsResult.setupTimeMs + cometsResult.executionTimeMs;
			cometsResult.numSequencesFound = cometsSeqs.size();
			for (int i = 0; i < Math.min(3, cometsSeqs.size()); i++) {
				cometsResult.topSequences.add(cometsSeqs.get(i).sequence.toString(Sequence.Renderer.ResTypeMutations));
			}
		});

		// Print results
		log("\n" + cometsZResult.toString());
		log("\n" + cometsResult.toString());
		log("\n========================================");
		log(String.format("Speedup: %.2fx (Setup: %.2fx, Execution: %.2fx)",
			(double) cometsResult.totalTimeMs / cometsZResult.totalTimeMs,
			(double) cometsResult.setupTimeMs / cometsZResult.setupTimeMs,
			(double) cometsResult.executionTimeMs / cometsZResult.executionTimeMs
		));
		log("========================================\n");
	}
}
