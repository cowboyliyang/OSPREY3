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
import static org.junit.jupiter.api.Assertions.assertThrows;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import org.junit.jupiter.api.Test;

import java.util.Random;

/**
 * The mean-field bound is checked against the exact partition function (by
 * independent brute-force enumeration): it must never exceed it, must reach it
 * when the model factorizes, and must beat the worst-case extreme-energy lower
 * bound that MARK* uses for an un-minimized subtree.
 */
public class TestMeanFieldBound {

	private static final double RT = 1.9891e-3 * 298.15;

	// ---- fixtures ----

	private static EnergyMatrix randomEmat(int[] numRCs, long seed, boolean coupled) {
		Random rand = new Random(seed);
		EnergyMatrix emat = new EnergyMatrix(numRCs.length, numRCs, 0.0);
		for (int pos = 0; pos < numRCs.length; pos++) {
			for (int rc = 0; rc < numRCs[pos]; rc++) {
				emat.setOneBody(pos, rc, rand.nextDouble() * 4.0 - 2.0);
			}
		}
		for (int p1 = 0; p1 < numRCs.length; p1++) {
			for (int p2 = p1 + 1; p2 < numRCs.length; p2++) {
				for (int r1 = 0; r1 < numRCs[p1]; r1++) {
					for (int r2 = 0; r2 < numRCs[p2]; r2++) {
						emat.setPairwise(p1, r1, p2, r2, coupled ? rand.nextDouble() * 2.0 - 1.0 : 0.0);
					}
				}
			}
		}
		return emat;
	}

	private static RCs allRcs(int[] numRCs) {
		int[][] rcsAtPos = new int[numRCs.length][];
		for (int pos = 0; pos < numRCs.length; pos++) {
			rcsAtPos[pos] = new int[numRCs[pos]];
			for (int rc = 0; rc < numRCs[pos]; rc++) {
				rcsAtPos[pos][rc] = rc;
			}
		}
		return new RCs(rcsAtPos);
	}

	private static int[] unassigned(int numPos) {
		int[] a = new int[numPos];
		java.util.Arrays.fill(a, -1);
		return a;
	}

	private static double[] bruteForce(EnergyMatrix emat, RCs rcs, int[] assignments) {
		// returns {logZ, minEnergy, maxEnergy} over all completions
		int numPos = rcs.getNumPos();
		java.util.List<Integer> free = new java.util.ArrayList<>();
		for (int pos = 0; pos < numPos; pos++) {
			if (assignments[pos] < 0) {
				free.add(pos);
			}
		}
		long n = 1;
		for (int pos : free) {
			n *= rcs.getNum(pos);
		}
		double[] logTerms = new double[(int) n];
		double minE = Double.POSITIVE_INFINITY;
		double maxE = Double.NEGATIVE_INFINITY;
		int[] full = assignments.clone();
		int[] choice = new int[free.size()];
		for (int idx = 0; idx < n; idx++) {
			for (int f = 0; f < free.size(); f++) {
				full[free.get(f)] = rcs.get(free.get(f), choice[f]);
			}
			double energy = 0.0;
			for (int i = 0; i < numPos; i++) {
				energy += emat.getEnergy(i, full[i]);
				for (int j = i + 1; j < numPos; j++) {
					energy += emat.getEnergy(i, full[i], j, full[j]);
				}
			}
			logTerms[idx] = -energy / RT;
			minE = Math.min(minE, energy);
			maxE = Math.max(maxE, energy);
			for (int f = free.size() - 1; f >= 0; f--) {
				choice[f]++;
				if (choice[f] < rcs.getNum(free.get(f))) {
					break;
				}
				choice[f] = 0;
			}
		}
		return new double[]{logSumExp(logTerms), minE, maxE, Math.log(n)};
	}

	private static double logSumExp(double[] v) {
		double max = Double.NEGATIVE_INFINITY;
		for (double x : v) {
			max = Math.max(max, x);
		}
		if (max == Double.NEGATIVE_INFINITY) {
			return Double.NEGATIVE_INFINITY;
		}
		double sum = 0.0;
		for (double x : v) {
			sum += Math.exp(x - max);
		}
		return max + Math.log(sum);
	}

	private static WmbModel model(EnergyMatrix emat, RCs rcs, int[] assign) {
		return new WmbModel(emat, rcs, assign, RT);
	}

	// ---- validity ----

	@Test
	public void lowerBoundBelowExactForCoupledModels() {
		for (long seed = 0; seed < 6; seed++) {
			int[] numRCs = {3, 2, 3, 2};
			EnergyMatrix emat = randomEmat(numRCs, seed, true);
			RCs rcs = allRcs(numRCs);
			double exact = bruteForce(emat, rcs, unassigned(numRCs.length))[0];
			double mf = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)));
			assertThat("seed " + seed, mf, lessThanOrEqualTo(exact + 1e-7));
		}
	}

	@Test
	public void exactWhenModelFactorizes() {
		int[] numRCs = {3, 2, 4, 2};
		EnergyMatrix emat = randomEmat(numRCs, 11L, false); // no pairwise coupling
		RCs rcs = allRcs(numRCs);
		double exact = bruteForce(emat, rcs, unassigned(numRCs.length))[0];
		double mf = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)));
		assertThat(mf, closeTo(exact, 1e-6));
	}

	@Test
	public void sweepsMonotonicallyTighten() {
		int[] numRCs = {3, 3, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 21L, true);
		RCs rcs = allRcs(numRCs);
		double afterOne = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)), 1, 1e-12);
		double afterMany = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)), 50, 1e-12);
		double exact = bruteForce(emat, rcs, unassigned(numRCs.length))[0];
		assertThat(afterMany, greaterThanOrEqualTo(afterOne - 1e-9));
		assertThat(afterMany, lessThanOrEqualTo(exact + 1e-7));
	}

	@Test
	public void beatsWorstCaseExtremeLowerBound() {
		// MARK*'s un-minimized subtree lower bound is N*exp(-E_max/RT): log N - E_max/RT.
		for (long seed = 30; seed < 35; seed++) {
			int[] numRCs = {3, 3, 2, 2};
			EnergyMatrix emat = randomEmat(numRCs, seed, true);
			RCs rcs = allRcs(numRCs);
			double[] bf = bruteForce(emat, rcs, unassigned(numRCs.length));
			double exact = bf[0];
			double maxE = bf[2];
			double logN = bf[3];
			double crude = logN - maxE / RT;
			double mf = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)));
			assertThat("seed " + seed + " vs crude", mf, greaterThanOrEqualTo(crude - 1e-9));
			assertThat("seed " + seed + " vs exact", mf, lessThanOrEqualTo(exact + 1e-7));
		}
	}

	@Test
	public void constTermIsIncluded() {
		int[] numRCs = {3, 2, 3};
		EnergyMatrix emat = randomEmat(numRCs, 4L, true);
		double constTerm = 7.5;
		emat.setConstTerm(constTerm);
		RCs rcs = allRcs(numRCs);
		double exactWithConst = bruteForce(emat, rcs, unassigned(numRCs.length))[0] - constTerm / RT;
		double mf = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)));
		assertThat(mf, lessThanOrEqualTo(exactWithConst + 1e-7));
	}

	@Test
	public void handlesHugeClashEnergyWithoutNaN() {
		int[] numRCs = {3, 3, 3};
		EnergyMatrix emat = randomEmat(numRCs, 5L, true);
		emat.setOneBody(1, 0, 1.0e6);            // a rigid-style clash rotamer
		emat.setPairwise(0, 1, 2, 1, 1.0e6);     // and a clashing pair
		RCs rcs = allRcs(numRCs);
		double exact = bruteForce(emat, rcs, unassigned(numRCs.length))[0];
		double mf = MeanFieldBound.logZLower(model(emat, rcs, unassigned(numRCs.length)));
		assertThat(Double.isNaN(mf), is(false));
		assertThat(Double.isFinite(mf), is(true));
		assertThat(mf, lessThanOrEqualTo(exact + 1e-6));
	}

	@Test
	public void conditioningOnPartialAssignmentStaysBelowExact() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 7L, true);
		RCs rcs = allRcs(numRCs);
		int[] assign = unassigned(numRCs.length);
		assign[0] = 1;
		assign[2] = 0;
		double exact = bruteForce(emat, rcs, assign)[0];
		double mf = MeanFieldBound.logZLower(model(emat, rcs, assign));
		assertThat(mf, lessThanOrEqualTo(exact + 1e-7));
	}

	@Test
	public void fullyAssignedReturnsTheSingleTerm() {
		int[] numRCs = {2, 2, 2};
		EnergyMatrix emat = randomEmat(numRCs, 3L, true);
		RCs rcs = allRcs(numRCs);
		int[] assign = {1, 0, 1};
		double exact = bruteForce(emat, rcs, assign)[0];
		double mf = MeanFieldBound.logZLower(model(emat, rcs, assign));
		assertThat(mf, closeTo(exact, 1e-9));
	}

	@Test
	public void maxSweepsBelowOneIsRejected() {
		int[] numRCs = {2, 2};
		EnergyMatrix emat = randomEmat(numRCs, 1L, true);
		RCs rcs = allRcs(numRCs);
		WmbModel m = model(emat, rcs, unassigned(2));
		assertThrows(IllegalArgumentException.class, () -> MeanFieldBound.logZLower(m, 0, 1e-9));
	}
}
