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
 * The exact partition function over a small energy-matrix subtree is the oracle:
 * it is computed by independent brute-force enumeration of completions, and the
 * mini-bucket bracket is checked against it.
 */
public class TestWeightedMiniBucket {

	private static final double RT = 1.9891e-3 * 298.15; // OSPREY Classic/Room

	// ---- fixtures ----

	private static EnergyMatrix randomEmat(int[] numRCs, long seed) {
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
						emat.setPairwise(p1, r1, p2, r2, rand.nextDouble() * 2.0 - 1.0);
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

	/** Independent oracle: exact log Z over all completions of the partial assignment. */
	private static double bruteForceLogZ(EnergyMatrix emat, RCs rcs, int[] assignments) {
		int numPos = rcs.getNumPos();
		java.util.List<Integer> free = new java.util.ArrayList<>();
		for (int pos = 0; pos < numPos; pos++) {
			if (assignments[pos] < 0) {
				free.add(pos);
			}
		}
		double[] logTerms = new double[(int) countCompletions(rcs, free)];
		int[] full = assignments.clone();
		int[] choice = new int[free.size()];
		for (int idx = 0; idx < logTerms.length; idx++) {
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
			for (int f = free.size() - 1; f >= 0; f--) {
				choice[f]++;
				if (choice[f] < rcs.getNum(free.get(f))) {
					break;
				}
				choice[f] = 0;
			}
		}
		return logSumExp(logTerms);
	}

	private static long countCompletions(RCs rcs, java.util.List<Integer> free) {
		long n = 1;
		for (int pos : free) {
			n *= rcs.getNum(pos);
		}
		return n;
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

	private static double logBaseline(EnergyMatrix emat, RCs rcs, int[] assignments) {
		// MARK*-style N * exp(-Emin/RT): log N - Emin/RT
		int numPos = rcs.getNumPos();
		java.util.List<Integer> free = new java.util.ArrayList<>();
		for (int pos = 0; pos < numPos; pos++) {
			if (assignments[pos] < 0) {
				free.add(pos);
			}
		}
		long n = countCompletions(rcs, free);
		double minEnergy = Double.POSITIVE_INFINITY;
		int[] full = assignments.clone();
		int[] choice = new int[free.size()];
		for (long idx = 0; idx < n; idx++) {
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
			minEnergy = Math.min(minEnergy, energy);
			for (int f = free.size() - 1; f >= 0; f--) {
				choice[f]++;
				if (choice[f] < rcs.getNum(free.get(f))) {
					break;
				}
				choice[f] = 0;
			}
		}
		return Math.log(n) - minEnergy / RT;
	}

	private static double confLogWeight(EnergyMatrix emat, int[] conf) {
		double energy = 0.0;
		for (int i = 0; i < conf.length; i++) {
			energy += emat.getEnergy(i, conf[i]);
			for (int j = i + 1; j < conf.length; j++) {
				energy += emat.getEnergy(i, conf[i], j, conf[j]);
			}
		}
		return -energy / RT;
	}

	// ---- validity / exactness ----

	@Test
	public void exactAtIBoundEqualInducedWidth() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 1L);
		RCs rcs = allRcs(numRCs);
		double exact = bruteForceLogZ(emat, rcs, unassigned(numRCs.length));

		int width = numRCs.length - 1; // dense graph -> complete
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(
				new WmbModel(emat, rcs, unassigned(numRCs.length), RT), width);
		assertThat(b.isExact(), is(true));
		assertThat(b.logZUpper, closeTo(exact, 1e-9));
		assertThat(b.logZLower, closeTo(exact, 1e-9));
	}

	@Test
	public void upperGeExactAndLowerLeExactForAllIBound() {
		for (long seed = 0; seed < 5; seed++) {
			int[] numRCs = {3, 2, 3, 2};
			EnergyMatrix emat = randomEmat(numRCs, seed);
			RCs rcs = allRcs(numRCs);
			WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);
			double exact = bruteForceLogZ(emat, rcs, unassigned(numRCs.length));
			for (int ib = 1; ib <= numRCs.length - 1; ib++) {
				MiniBucketBound b = WeightedMiniBucket.boundsForModel(model, ib);
				assertThat("seed " + seed + " ib " + ib + " upper",
						b.logZUpper, greaterThanOrEqualTo(exact - 1e-7));
				assertThat("seed " + seed + " ib " + ib + " lower",
						b.logZLower, lessThanOrEqualTo(exact + 1e-7));
			}
		}
	}

	@Test
	public void constTermIsIncludedInTheBound() {
		// the energy-matrix constant offset shifts every conformation energy, so it
		// must shift log Z by -constTerm/RT; the brute-force oracle here ignores the
		// constant, so the bound should sit exactly that far below it.
		int[] numRCs = {3, 2, 3};
		EnergyMatrix emat = randomEmat(numRCs, 4L);
		double constTerm = 7.5;
		emat.setConstTerm(constTerm);
		RCs rcs = allRcs(numRCs);
		int width = numRCs.length - 1;
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(
				new WmbModel(emat, rcs, unassigned(numRCs.length), RT), width);
		double exactWithConst = bruteForceLogZ(emat, rcs, unassigned(numRCs.length)) - constTerm / RT;
		assertThat(b.logZUpper, closeTo(exactWithConst, 1e-9));
		assertThat(b.logZLower, closeTo(exactWithConst, 1e-9));
	}

	@Test
	public void conditioningOnPartialAssignmentMatchesBruteForce() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 7L);
		RCs rcs = allRcs(numRCs);
		int[] assign = unassigned(numRCs.length);
		assign[0] = 1; // fix the first position to rotamer 1
		assign[2] = 0; // fix the third position to rotamer 0
		double exact = bruteForceLogZ(emat, rcs, assign);

		int width = WeightedMiniBucket.inducedWidth(new WmbModel(emat, rcs, assign, RT),
				WeightedMiniBucket.naturalOrder(2));
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(
				new WmbModel(emat, rcs, assign, RT), Math.max(width, 1));
		assertThat(b.logZUpper, closeTo(exact, 1e-9));
		assertThat(b.logZLower, closeTo(exact, 1e-9));
	}

	@Test
	public void exactProposalLogProbabilityMatchesBoltzmannDistribution() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 11L);
		RCs rcs = allRcs(numRCs);
		int[] assignments = unassigned(numRCs.length);
		WmbModel model = new WmbModel(emat, rcs, assignments, RT);
		double logZ = bruteForceLogZ(emat, rcs, assignments);
		WeightedMiniBucket.Proposal proposal = WeightedMiniBucket.proposalForModel(model, numRCs.length - 1);

		int[] domainValues = new int[numRCs.length];
		int n = (int) countCompletions(rcs, java.util.Arrays.asList(0, 1, 2, 3));
		for (int idx = 0; idx < n; idx++) {
			int[] conf = new int[numRCs.length];
			for (int v = 0; v < domainValues.length; v++) {
				conf[model.position(v)] = model.rotamer(v, domainValues[v]);
			}
			assertThat(proposal.logProbability(domainValues),
					closeTo(confLogWeight(emat, conf) - logZ, 1e-9));
			for (int v = domainValues.length - 1; v >= 0; v--) {
				domainValues[v]++;
				if (domainValues[v] < numRCs[v]) {
					break;
				}
				domainValues[v] = 0;
			}
		}
	}

	@Test
	public void exactLocalImportanceWeightCapCoversAllAssignments() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 13L);
		RCs rcs = allRcs(numRCs);
		int[] assignments = unassigned(numRCs.length);
		WmbModel model = new WmbModel(emat, rcs, assignments, RT);
		MiniBucketBound bound = WeightedMiniBucket.boundsForModel(model, 1);
		WeightedMiniBucket.Proposal proposal = WeightedMiniBucket.proposalForModel(model, 1);
		WeightedMiniBucket.Proposal.LogWeightCap cap =
			proposal.logWeightUpperBound(bound.logZUpper, 1000L);
		WeightedMiniBucket.Proposal.LogWeightCap miniBucketCap =
			proposal.logWeightUpperBound(bound.logZUpper, 0L);

		assertThat(cap.exact, is(true));
		assertThat(miniBucketCap.exact, is(false));
		assertThat(cap.assignments, is(countCompletions(rcs, java.util.Arrays.asList(0, 1, 2, 3))));
		assertThat(cap.logWeightUpper, lessThanOrEqualTo(cap.fallbackLogWeightUpper + 1e-12));
		assertThat(miniBucketCap.logWeightUpper, lessThanOrEqualTo(miniBucketCap.fallbackLogWeightUpper + 1e-12));

		int[] domainValues = new int[numRCs.length];
		int n = (int) cap.assignments;
		for (int idx = 0; idx < n; idx++) {
			double logWeight = model.logValue(domainValues) - proposal.logProbability(domainValues);
			assertThat(logWeight, lessThanOrEqualTo(cap.logWeightUpper + 1e-12));
			assertThat(logWeight, lessThanOrEqualTo(miniBucketCap.logWeightUpper + 1e-12));
			for (int v = domainValues.length - 1; v >= 0; v--) {
				domainValues[v]++;
				if (domainValues[v] < numRCs[v]) {
					break;
				}
				domainValues[v] = 0;
			}
		}
	}

	@Test
	public void exactProposalLocalImportanceWeightCapEqualsLogZ() {
		int[] numRCs = {3, 2, 3};
		EnergyMatrix emat = randomEmat(numRCs, 15L);
		RCs rcs = allRcs(numRCs);
		int[] assignments = unassigned(numRCs.length);
		WmbModel model = new WmbModel(emat, rcs, assignments, RT);
		double logZ = bruteForceLogZ(emat, rcs, assignments);
		WeightedMiniBucket.Proposal proposal = WeightedMiniBucket.proposalForModel(model, numRCs.length - 1);
		WeightedMiniBucket.Proposal.LogWeightCap cap =
			proposal.logWeightUpperBound(logZ, 1000L);

		assertThat(cap.exact, is(true));
		assertThat(cap.logWeightUpper, closeTo(logZ, 1e-9));
	}

	@Test
	public void proposalLogProbabilityLowerBoundCoversAllAssignments() {
		int[] numRCs = {3, 2, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 17L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);
		WeightedMiniBucket.Proposal proposal = WeightedMiniBucket.proposalForModel(model, 1);
		double lower = proposal.logProbabilityLowerBound();

		int[] domainValues = new int[numRCs.length];
		int n = (int) countCompletions(rcs, java.util.Arrays.asList(0, 1, 2, 3));
		for (int idx = 0; idx < n; idx++) {
			assertThat(proposal.logProbability(domainValues), greaterThanOrEqualTo(lower - 1e-12));
			for (int v = domainValues.length - 1; v >= 0; v--) {
				domainValues[v]++;
				if (domainValues[v] < numRCs[v]) {
					break;
				}
				domainValues[v] = 0;
			}
		}
	}

	@Test
	public void looseProposalImportanceSamplingMatchesEnumeratedPartitionFunction() {
		int[] numRCs = {3, 2, 3, 2, 2};
		EnergyMatrix emat = randomEmat(numRCs, 19L);
		RCs rcs = allRcs(numRCs);
		int[] assignments = unassigned(numRCs.length);
		WmbModel model = new WmbModel(emat, rcs, assignments, RT);
		WeightedMiniBucket.Proposal proposal = WeightedMiniBucket.proposalForModel(model, 1);

		int samples = 20000;
		Random rng = new Random(23L);
		double[] logWeights = new double[samples];
		for (int i = 0; i < samples; i++) {
			WeightedMiniBucket.Sample sample = proposal.sample(rng);
			int[] conf = new int[numRCs.length];
			for (int v = 0; v < sample.domainValues.length; v++) {
				conf[model.position(v)] = model.rotamer(v, sample.domainValues[v]);
			}
			logWeights[i] = confLogWeight(emat, conf) - sample.logQ;
		}

		double exact = bruteForceLogZ(emat, rcs, assignments);
		double estimate = logSumExp(logWeights) - Math.log(samples);
		System.out.printf("WMB-IS small validation logZ exact=%.10f estimate=%.10f absErr=%.10f samples=%d%n",
				exact, estimate, Math.abs(estimate - exact), samples);
		assertThat(estimate, closeTo(exact, 0.08));
	}

	@Test
	public void fullyAssignedSubtreeIsASingleTerm() {
		int[] numRCs = {2, 2, 2};
		EnergyMatrix emat = randomEmat(numRCs, 3L);
		RCs rcs = allRcs(numRCs);
		int[] assign = {1, 0, 1};
		double exact = bruteForceLogZ(emat, rcs, assign); // single completion
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(
				new WmbModel(emat, rcs, assign, RT), 1);
		assertThat(b.logZUpper, closeTo(exact, 1e-9));
		assertThat(b.logZLower, closeTo(exact, 1e-9));
	}

	// ---- tighter than the MARK* extreme baseline ----

	@Test
	public void certifiedUpperIsValidAndNeverLooserThanBaseline() {
		for (long seed = 10; seed < 15; seed++) {
			int[] numRCs = {3, 3, 2, 2};
			EnergyMatrix emat = randomEmat(numRCs, seed);
			RCs rcs = allRcs(numRCs);
			WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);
			double exact = bruteForceLogZ(emat, rcs, unassigned(numRCs.length));
			double baseline = logBaseline(emat, rcs, unassigned(numRCs.length));
			for (int ib = 1; ib <= numRCs.length - 1; ib++) {
				double wmbUpper = WeightedMiniBucket.boundsForModel(model, ib).logZUpper;
				double certified = Math.min(wmbUpper, baseline);
				assertThat(certified, greaterThanOrEqualTo(exact - 1e-7));
				assertThat(certified, lessThanOrEqualTo(baseline + 1e-9));
			}
		}
	}

	@Test
	public void gapTightensMonotonicallyWithIBound() {
		int[] numRCs = {3, 3, 3, 2};
		EnergyMatrix emat = randomEmat(numRCs, 21L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);
		double prev = Double.POSITIVE_INFINITY;
		for (int ib = 1; ib <= numRCs.length - 1; ib++) {
			double gap = WeightedMiniBucket.boundsForModel(model, ib).logGap();
			assertThat("ib " + ib, gap, lessThanOrEqualTo(prev + 1e-6));
			prev = gap;
		}
		assertThat(prev, closeTo(0.0, 1e-6));
	}

	// ---- cost shape, decoupled from conf count ----

	@Test
	public void tableCellsBoundedByDomainPowIBoundPlusOne() {
		int[] numRCs = {3, 3, 3, 3, 3};
		EnergyMatrix emat = randomEmat(numRCs, 33L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);
		int q = 3;
		for (int ib = 1; ib <= numRCs.length - 1; ib++) {
			MiniBucketBound b = WeightedMiniBucket.boundsForModel(model, ib);
			assertThat(b.maxMiniBucketVars, lessThanOrEqualTo(ib + 1));
			assertThat(b.maxTableCells, lessThanOrEqualTo((long) Math.pow(q, ib + 1)));
		}
	}

	@Test
	public void maxTableCellsLowersEffectiveIBound() {
		int[] numRCs = {3, 3, 3, 3};
		EnergyMatrix emat = randomEmat(numRCs, 34L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);

		MiniBucketBound b = WeightedMiniBucket.boundsForModel(model, 3, 9L);
		assertThat(b.iBound, is(1));
		assertThat(b.inducedWidth, is(3));
		assertThat(b.isExact(), is(false));
		assertThat(b.maxMiniBucketVars, lessThanOrEqualTo(2));
		assertThat(b.maxTableCells, lessThanOrEqualTo(9L));

		MiniBucketBound stateBound = WeightedMiniBucket.bounds(
			emat, emat, rcs, unassigned(numRCs.length), 3, 9L, RT);
		assertThat(stateBound.iBound, is(1));
		assertThat(stateBound.maxTableCells, lessThanOrEqualTo(9L));
		assertThat(WeightedMiniBucket.upperLogZ(
			emat, rcs, unassigned(numRCs.length), 3, 9L, RT), closeTo(b.logZUpper, 1e-9));
	}

	@Test
	public void maxTableCellsTooSmallForPairFactorsIsRejected() {
		int[] numRCs = {4, 4, 4};
		EnergyMatrix emat = randomEmat(numRCs, 35L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numRCs.length), RT);

		assertThrows(IllegalArgumentException.class,
				() -> WeightedMiniBucket.boundsForModel(model, 2, 15L));
		assertThrows(IllegalArgumentException.class,
				() -> WeightedMiniBucket.proposalForModel(model, 2, 15L));
	}

	@Test
	public void decoupledFromConformationCount() {
		// 30 binary positions: 2^30 completions, but iBound=2 stays cheap and finite.
		int numPos = 30;
		int[] numRCs = new int[numPos];
		java.util.Arrays.fill(numRCs, 2);
		EnergyMatrix emat = randomEmat(numRCs, 99L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(numPos), RT);
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(model, 2);
		assertThat(Double.isFinite(b.logZUpper), is(true));
		assertThat(Double.isFinite(b.logZLower), is(true));
		assertThat(b.logZLower, lessThanOrEqualTo(b.logZUpper + 1e-9));
		assertThat(b.maxMiniBucketVars, lessThanOrEqualTo(3));
		assertThat(b.maxTableCells, lessThanOrEqualTo(8L));
	}

	// ---- negative tests ----

	@Test
	public void iBoundBelowOneIsRejected() {
		int[] numRCs = {2, 2};
		EnergyMatrix emat = randomEmat(numRCs, 1L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(2), RT);
		assertThrows(IllegalArgumentException.class,
				() -> WeightedMiniBucket.boundsForModel(model, 0));
	}

	@Test
	public void treatingTheLowerBoundAsAnUpperBoundWouldFailValidity() {
		// Guards against a min/sum mixup: the lower bound must not exceed exact Z.
		int[] numRCs = {3, 3, 3};
		EnergyMatrix emat = randomEmat(numRCs, 5L);
		RCs rcs = allRcs(numRCs);
		WmbModel model = new WmbModel(emat, rcs, unassigned(3), RT);
		double exact = bruteForceLogZ(emat, rcs, unassigned(3));
		MiniBucketBound b = WeightedMiniBucket.boundsForModel(model, 1);
		assertThat(b.logZLower, lessThanOrEqualTo(exact + 1e-7));
		boolean wouldBeValidUpper = b.logZLower >= exact - 1e-7;
		assertThat("lower bound must be strictly below exact here", wouldBeValidUpper, is(false));
	}
}
