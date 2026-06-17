package edu.duke.cs.osprey.branchdp;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.util.Arrays;

/**
 * 0-edge direct independent-position DP (design 1.2) vs brute-force joint
 * enumeration, fully local (no conf space / cluster).
 *
 * When the interaction graph has no pairwise edges, Z factorizes per position:
 *   logZ = sum_pos logsumexp_rc( -E(pos,rc) / RT )
 * which RootedTreeEdge.independentPositionLogZ computes in O(sum cards). These
 * tests check that value equals the brute-force log-sum-exp over the FULL
 * Cartesian product of RCs (O(prod cards)), for both the lower bound (rigid
 * emat) and the upper bound (min emat) -- the regression the design doc asks for
 * and that the method's javadoc claims ("Validated against brute-force ... in
 * TestZeroEdgeDP").
 */
public class TestZeroEdgeDP {

    private static final double RT = 1.9; // arbitrary positive RT

    /** RCs with explicit RC indices 0..card-1 at each position (so get(pos,i)==i). */
    private static RCs rcs(int[] cards) {
        int[][] arr = new int[cards.length][];
        for (int p = 0; p < cards.length; p++) {
            arr[p] = new int[cards[p]];
            for (int i = 0; i < cards[p]; i++) arr[p][i] = i;
        }
        return new RCs(arr);
    }

    /** One-body-only emat with deterministic, distinct, sign-varying energies. */
    private static EnergyMatrix emat(int[] cards, double seed) {
        EnergyMatrix e = new EnergyMatrix(cards.length, cards, 0.0);
        for (int p = 0; p < cards.length; p++) {
            for (int rc = 0; rc < cards[p]; rc++) {
                e.setOneBody(p, rc, seed + 0.37 * p - 0.61 * rc + 0.05 * p * rc);
            }
        }
        return e;
    }

    /** Brute force: logsumexp over the full Cartesian product of RCs of -Etotal/RT. */
    private static double bruteForceLogZ(EnergyMatrix e, RCs r, int[] positions, double RT) {
        int n = positions.length;
        int[] card = new int[n];
        long total = 1;
        for (int i = 0; i < n; i++) { card[i] = r.getNum(positions[i]); total *= card[i]; }

        double max = Double.NEGATIVE_INFINITY; // online (max-shift) log-sum-exp
        double sum = 0.0;
        int[] idx = new int[n];
        for (long c = 0; c < total; c++) {
            double energy = 0.0;
            for (int i = 0; i < n; i++) {
                energy += e.getOneBody(positions[i], r.get(positions[i], idx[i]));
            }
            double x = -energy / RT;
            if (x > max) { sum = sum * Math.exp(max - x) + 1.0; max = x; }
            else { sum += Math.exp(x - max); }
            for (int i = n - 1; i >= 0; i--) { if (++idx[i] < card[i]) break; idx[i] = 0; }
        }
        return max + Math.log(sum);
    }

    private void checkAllPositions(int[] cards) {
        RCs r = rcs(cards);
        EnergyMatrix rigid = emat(cards, 1.0);
        EnergyMatrix min = emat(cards, -0.5); // distinct so lower != upper
        int[] positions = new int[cards.length];
        for (int i = 0; i < cards.length; i++) positions[i] = i;

        double[] z = RootedTreeEdge.independentPositionLogZ(rigid, min, r, positions, RT);
        double expLower = bruteForceLogZ(rigid, r, positions, RT);
        double expUpper = bruteForceLogZ(min, r, positions, RT);

        String where = " cards=" + Arrays.toString(cards);
        assertEquals(expLower, z[0], 1e-9 * (1.0 + Math.abs(expLower)), "lower (rigid) vs brute force," + where);
        assertEquals(expUpper, z[1], 1e-9 * (1.0 + Math.abs(expUpper)), "upper (min) vs brute force," + where);
    }

    @Test
    public void matchesBruteForce_variousShapes() {
        checkAllPositions(new int[]{1});
        checkAllPositions(new int[]{2});
        checkAllPositions(new int[]{2, 3});
        checkAllPositions(new int[]{3, 1, 4});
        checkAllPositions(new int[]{2, 2, 2, 2});
        checkAllPositions(new int[]{5, 4, 3, 2});
    }

    @Test
    public void positionSubset_usesOnlyGivenPositions() {
        // 4 positions exist but only a non-contiguous subset participates
        int[] cards = {3, 2, 4, 2};
        RCs r = rcs(cards);
        EnergyMatrix rigid = emat(cards, 0.8);
        EnergyMatrix min = emat(cards, 0.2);
        int[] positions = {0, 2}; // skip 1 and 3

        double[] z = RootedTreeEdge.independentPositionLogZ(rigid, min, r, positions, RT);
        assertEquals(bruteForceLogZ(rigid, r, positions, RT), z[0], 1e-9, "lower over subset");
        assertEquals(bruteForceLogZ(min, r, positions, RT), z[1], 1e-9, "upper over subset");
    }

    @Test
    public void singleRcPerPosition_isJustNegEnergySum() {
        // every position has exactly one RC -> logsumexp is a no-op, logZ = sum(-E/RT)
        int[] cards = {1, 1, 1};
        RCs r = rcs(cards);
        EnergyMatrix rigid = emat(cards, 1.3);
        int[] positions = {0, 1, 2};
        double expected = 0.0;
        for (int p : positions) expected += -rigid.getOneBody(p, 0) / RT;

        double[] z = RootedTreeEdge.independentPositionLogZ(rigid, rigid, r, positions, RT);
        assertEquals(expected, z[0], 1e-12);
        assertEquals(expected, z[1], 1e-12);
    }
}
