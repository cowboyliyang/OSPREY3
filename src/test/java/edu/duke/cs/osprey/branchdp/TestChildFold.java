package edu.duke.cs.osprey.branchdp;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.LinkedHashSet;

/**
 * Option C (non-leaf child-table folding) correctness, fully local (no cluster).
 *
 * Runs the real RootedTreeEdge.computeFullDP() on a small synthetic non-leaf
 * edge with three F-set children (one lambda-dependent, two lambda-invariant)
 * under three configs and compares the resulting DP table:
 *   - foldChildren=false                 (legacy getMstateForFullState path)
 *   - foldChildren=true, hoist=false     (default; must be BIT-IDENTICAL)
 *   - foldChildren=true, hoist=true      (lambda-invariant hoist; equal within FP tol)
 *
 * The energy expressions are identical across all three paths, so any difference
 * can only come from the child-index folding / summation, which is exactly what
 * Option C changes. Complements TestMStateIndexing (computeIndexInA long /
 * decodeMState round-trip) and TestDPTable (dense vs sharded).
 */
public class TestChildFold {

    // ---- reflection / construction helpers --------------------------------

    private static void set(Object o, String field, Object val) {
        try {
            Field f = RootedTreeEdge.class.getDeclaredField(field);
            f.setAccessible(true);
            f.set(o, val);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("set " + field, e);
        }
    }

    /** RCs whose getNum(positions[i]) == cardinalities[i]; other positions have 1 RC. */
    private static RCs rcsFor(int[] positions, int[] cardinalities) {
        int maxPos = 0;
        for (int p : positions) maxPos = Math.max(maxPos, p);
        int[][] arr = new int[maxPos + 1][];
        for (int i = 0; i < arr.length; i++) arr[i] = new int[1];
        for (int i = 0; i < positions.length; i++) arr[positions[i]] = new int[cardinalities[i]];
        return new RCs(arr);
    }

    private static EnergyMatrix ematFilled(int numPos, int[] cards, double seed) {
        EnergyMatrix e = new EnergyMatrix(numPos, cards, 0.0);
        for (int p = 0; p < numPos; p++) {
            for (int rc = 0; rc < cards[p]; rc++) {
                e.setOneBody(p, rc, seed + 0.13 * p - 0.07 * rc);
            }
        }
        for (int pi = 0; pi < numPos; pi++) {
            for (int pj = pi + 1; pj < numPos; pj++) {
                for (int ri = 0; ri < cards[pi]; ri++) {
                    for (int rj = 0; rj < cards[pj]; rj++) {
                        e.setPairwise(pi, ri, pj, rj, 0.31 * seed + 0.05 * (ri + 1) * (rj + 2) - 0.02 * (pi + pj));
                    }
                }
            }
        }
        return e;
    }

    private static InteractionGraph fullyConnected(int numPos) {
        boolean[][] adj = new boolean[numPos][numPos];
        for (int i = 0; i < numPos; i++) {
            for (int j = 0; j < numPos; j++) {
                if (i != j) adj[i][j] = true;
            }
        }
        try {
            Constructor<InteractionGraph> c =
                    InteractionGraph.class.getDeclaredConstructor(int.class, boolean[][].class);
            c.setAccessible(true);
            return c.newInstance(numPos, adj);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("build InteractionGraph", e);
        }
    }

    private static long stateCount(RCs rcs, int[] pos) {
        long c = 1;
        for (int p : pos) c *= rcs.getNum(p);
        return c;
    }

    /** A child lambda-edge with given M positions and a deterministic, distinct DP table. */
    private static RootedTreeEdge child(RCs rcs, int[] mPos, double base) {
        RootedTreeEdge c = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(c, "mPositionsSorted", mPos.clone());
        set(c, "isLambdaEdge", Boolean.TRUE);
        long mc = stateCount(rcs, mPos);
        set(c, "mStateCount", Long.valueOf(mc));
        set(c, "mArraySize", Integer.valueOf((int) mc));
        DenseDPTable t = new DenseDPTable(mc);
        for (long i = 0; i < mc; i++) t.set(i, base + i, base + 100.0 + i);
        set(c, "dpTable", t);
        return c;
    }

    /** Build a fresh parent edge + 3 children and run computeFullDP under the given flags. */
    private double[][] runDP(boolean fold, boolean hoist) {
        System.setProperty("branchmarkstar.dp.foldChildren", String.valueOf(fold));
        System.setProperty("branchmarkstar.dp.foldChildren.hoistInvariant", String.valueOf(hoist));
        System.setProperty("branchmarkstar.dp.progress", "false");
        System.setProperty("branchmarkstar.dp.parallel", "false");

        int numPos = 4;
        int[] cards = {2, 2, 3, 2};                 // pos0..3 cardinalities
        int[] mPos = {0, 1};
        int[] lambdaPos = {2, 3};
        RCs rcs = rcsFor(new int[]{0, 1, 2, 3}, cards);

        // Fset order [dependent, invariant, invariant] so the hoist genuinely re-associates.
        LinkedHashSet<RootedTreeEdge> fset = new LinkedHashSet<>();
        fset.add(child(rcs, new int[]{1, 2}, 1000.0)); // depends on lambda pos 2
        fset.add(child(rcs, new int[]{0, 1}, 2000.0)); // lambda-invariant (subset of M)
        fset.add(child(rcs, new int[]{0},    3000.0)); // lambda-invariant

        RootedTreeEdge parent = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(parent, "isLambdaEdge", Boolean.TRUE);
        set(parent, "mPositionsSorted", mPos.clone());
        set(parent, "lambdaPositionsSorted", lambdaPos.clone());
        long mc = stateCount(rcs, mPos);
        set(parent, "mStateCount", Long.valueOf(mc));
        set(parent, "mArraySize", Integer.valueOf((int) mc));
        set(parent, "totalLambdaStates", Integer.valueOf((int) stateCount(rcs, lambdaPos)));
        set(parent, "dpTable", new DenseDPTable(mc));
        set(parent, "Fset", fset);
        set(parent, "cachedRigidEmat", ematFilled(numPos, cards, 1.0));
        set(parent, "cachedMinEmat", ematFilled(numPos, cards, 2.0));
        set(parent, "cachedG", fullyConnected(numPos));
        set(parent, "cachedRT", Double.valueOf(1.9));

        parent.computeFullDP();

        double[] lo = new double[(int) mc];
        double[] up = new double[(int) mc];
        for (int i = 0; i < mc; i++) {
            lo[i] = parent.getLogZLower(i);
            up[i] = parent.getLogZUpper(i);
        }
        return new double[][]{lo, up};
    }

    @Test
    public void foldedDefault_isBitIdenticalToLegacy() {
        double[][] legacy = runDP(false, false);
        double[][] folded = runDP(true, false);
        for (int i = 0; i < legacy[0].length; i++) {
            assertEquals(Double.doubleToLongBits(legacy[0][i]), Double.doubleToLongBits(folded[0][i]),
                    "lower bit-identity at mIdx=" + i);
            assertEquals(Double.doubleToLongBits(legacy[1][i]), Double.doubleToLongBits(folded[1][i]),
                    "upper bit-identity at mIdx=" + i);
        }
    }

    @Test
    public void hoistInvariant_isEquivalentWithinTolerance() {
        double[][] legacy = runDP(false, false);
        double[][] hoisted = runDP(true, true);
        for (int i = 0; i < legacy[0].length; i++) {
            assertEquals(legacy[0][i], hoisted[0][i], 1e-9 * (1.0 + Math.abs(legacy[0][i])),
                    "lower hoist-equivalence at mIdx=" + i);
            assertEquals(legacy[1][i], hoisted[1][i], 1e-9 * (1.0 + Math.abs(legacy[1][i])),
                    "upper hoist-equivalence at mIdx=" + i);
        }
        // sanity: the values are real, finite, and not all equal across M-states
        assertTrue(Double.isFinite(legacy[0][0]));
        assertNotEquals(legacy[0][0], legacy[0][legacy[0].length - 1]);
    }
}
