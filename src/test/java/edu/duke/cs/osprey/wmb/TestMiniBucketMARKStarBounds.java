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
** Contact Info:
**    Bruce Donald, Duke University, Department of Computer Science
**    e-mail: www.cs.duke.edu/brd/
*/

package edu.duke.cs.osprey.wmb;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

/**
 * Direct checks for WMB/MF partition-function bounds.
 */
public class TestMiniBucketMARKStarBounds {

	private static SimpleConfSpace confSpace;
	private static EnergyMatrix minimizingEmat;
	private static EnergyMatrix rigidEmat;
	private static RCs rcs;
	private static final double RT = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;

	@BeforeAll
	static void setup() {
		ForcefieldParams ffparams = new ForcefieldParams();
		Molecule mol = PDBIO.readFile("examples/python.KStar/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		confSpace = new SimpleConfSpace.Builder().addStrand(ligand).build();

		Parallelism parallelism = Parallelism.makeCpu(2);
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(parallelism)
			.build();
		EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
			.setIsMinimizing(false)
			.build();

		ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
			.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
				.build()
				.calcReferenceEnergies())
			.setEnergyPartition(EnergyPartition.Traditional)
			.build();
		minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
			.build()
			.calcEnergyMatrix();

		ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(confSpace, rigidEcalc)
			.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, rigidEcalc)
				.build()
				.calcReferenceEnergies())
			.setEnergyPartition(EnergyPartition.Traditional)
			.build();
		rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
			.build()
			.calcEnergyMatrix();

		rcs = new RCs(confSpace);
		rigidEcalc.clean();
		minimizingEcalc.clean();
	}

	@Test
	@Timeout(value = 4, unit = TimeUnit.MINUTES)
	public void miniBucketUpperBoundsExactLogZ() {
		int[] assignments = unassigned();
		double exact = exactLogZ(minimizingEmat, rcs);
		double upper = WeightedMiniBucket.upperLogZ(minimizingEmat, rcs, assignments, 2, RT);

		assertThat(upper, greaterThanOrEqualTo(exact - 1e-8));
	}

	@Test
	@Timeout(value = 4, unit = TimeUnit.MINUTES)
	public void meanFieldLowerBoundsExactLogZ() {
		int[] assignments = unassigned();
		double exact = exactLogZ(rigidEmat, rcs);
		double lower = MeanFieldBound.lowerLogZ(rigidEmat, rcs, assignments,
			100, MeanFieldBound.DEFAULT_TOLERANCE, RT);

		assertThat(lower, lessThanOrEqualTo(exact + 1e-8));
		assertThat(lower, not(Double.NEGATIVE_INFINITY));
	}

	private static int[] unassigned() {
		int[] assignments = new int[rcs.getNumPos()];
		Arrays.fill(assignments, -1);
		return assignments;
	}

	private static double exactLogZ(EnergyMatrix emat, RCs rcs) {
		return exactLogZ(emat, rcs, new int[rcs.getNumPos()], 0, Double.NEGATIVE_INFINITY);
	}

	private static double exactLogZ(EnergyMatrix emat, RCs rcs, int[] conf, int pos, double logZ) {
		if (pos == rcs.getNumPos()) {
			return logAdd(logZ, -emat.confE(conf) / RT);
		}
		for (int i = 0; i < rcs.getNum(pos); i++) {
			conf[pos] = rcs.get(pos, i);
			logZ = exactLogZ(emat, rcs, conf, pos + 1, logZ);
		}
		return logZ;
	}

	private static double logAdd(double a, double b) {
		if (a == Double.NEGATIVE_INFINITY) {
			return b;
		}
		if (b == Double.NEGATIVE_INFINITY) {
			return a;
		}
		double max = Math.max(a, b);
		return max + Math.log(Math.exp(a - max) + Math.exp(b - max));
	}
}
