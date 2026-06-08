package edu.duke.cs.osprey.markstar.framework.branch;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.LinkedHashSet;
import java.util.Locale;

/**
 * Numeric validation of the CUDA GPU full-DP fast path
 * (dp.cu full_dp_n_children + DPGpuFullDP, gated by branchmarkstar.dp.gpu=true)
 * against the Java DP path, on the same synthetic non-leaf edge with N children.
 *
 * The GPU run goes through the real gate (tryComputeFullDPGpu via reflection) so
 * the test asserts the GPU path actually fired; if no CUDA GPU is present the gate
 * returns false and the test is skipped (assumeTrue) rather than vacuously passing.
 * The Java reference uses the LEGACY non-fold path (foldChildren=false,
 * getMstateForFullState), an INDEPENDENT child-index computation from the GPU's
 * ChildFoldPlan, so agreement cross-validates the fold/CSR mapping too.
 *
 * A fully-connected graph guarantees real lambda-M pair terms (lmPairs > 0).
 * Multi-child shapes mix lambda-dependent and lambda-invariant children (as in
 * TestChildFold) to exercise the per-child CSR offsets + table bases.
 *
 * Not bit-identical by design (two-pass log-sum-exp + block-parallel summation +
 * GPU FP), so the comparison is to a relative tolerance; the measured max abs/rel
 * difference is printed for the record.
 */
public class TestGpuFullDP {

    private static final double REL_TOL = 1e-6;

    // ---- reflection helpers -----------------------------------------------

    private static void set(Object o, String field, Object val) {
        try {
            Field f = RootedTreeEdge.class.getDeclaredField(field);
            f.setAccessible(true);
            f.set(o, val);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("set " + field, e);
        }
    }

    private static Object invoke(Object o, String method) {
        try {
            Method m = RootedTreeEdge.class.getDeclaredMethod(method);
            m.setAccessible(true);
            return m.invoke(o);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("invoke " + method, e);
        }
    }

    // ---- synthetic edge construction (mirrors TestChildFold) ---------------

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
        for (int i = 0; i < numPos; i++)
            for (int j = 0; j < numPos; j++)
                if (i != j) adj[i][j] = true;
        try {
            Constructor<InteractionGraph> c =
                    InteractionGraph.class.getDeclaredConstructor(int.class, boolean[][].class);
            c.setAccessible(true);
            return c.newInstance(numPos, adj);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("build InteractionGraph", e);
        }
    }

    private static long product(int[] cards, int[] pos) {
        long c = 1;
        for (int p : pos) c *= cards[p];
        return c;
    }

    /** A child lambda-edge with the given M positions and a deterministic dense DP table. */
    private static RootedTreeEdge child(RCs rcs, int[] mPos, int[] cards, double base) {
        RootedTreeEdge c = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(c, "mPositionsSorted", mPos.clone());
        set(c, "isLambdaEdge", Boolean.TRUE);
        long mc = product(cards, mPos);
        set(c, "mStateCount", Long.valueOf(mc));
        set(c, "mArraySize", Integer.valueOf((int) mc));
        DenseDPTable t = new DenseDPTable(mc);
        for (long i = 0; i < mc; i++) t.set(i, base + 0.5 * i, base + 50.0 + 0.3 * i);
        set(c, "dpTable", t);
        return c;
    }

    /** Fresh, fully-initialized N-child parent edge ready for computeFullDP. */
    private RootedTreeEdge buildParent(int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions) {
        return buildParent(cards, mPos, lambdaPos, childMPositions, false, 1);
    }

    private RootedTreeEdge buildParent(int[] cards, int[] mPos, int[] lambdaPos,
                                      int[][] childMPositions, boolean shardedParent, int shardSize) {
        int numPos = cards.length;
        int[] allPos = new int[numPos];
        for (int i = 0; i < numPos; i++) allPos[i] = i;
        RCs rcs = rcsFor(allPos, cards);

        LinkedHashSet<RootedTreeEdge> fset = new LinkedHashSet<>();
        for (int c = 0; c < childMPositions.length; c++) {
            fset.add(child(rcs, childMPositions[c], cards, 1000.0 * (c + 1)));
        }

        RootedTreeEdge parent = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(parent, "isLambdaEdge", Boolean.TRUE);
        set(parent, "mPositionsSorted", mPos.clone());
        set(parent, "lambdaPositionsSorted", lambdaPos.clone());
        long mc = product(cards, mPos);
        set(parent, "mStateCount", Long.valueOf(mc));
        set(parent, "mArraySize", Integer.valueOf((int) mc));
        set(parent, "totalLambdaStates", Integer.valueOf((int) product(cards, lambdaPos)));
        set(parent, "dpTable", shardedParent
                ? new ShardedDPTable(mc, shardSize)
                : new DenseDPTable(mc));
        set(parent, "Fset", fset);

        // Populates lambdaOnlyRigid/Min + cached emat/G/RT (required by the GPU gate).
        parent.initIncrementalEnumeration(
                ematFilled(numPos, cards, 1.0), ematFilled(numPos, cards, 2.0),
                fullyConnected(numPos), 1.9);
        return parent;
    }

    private static double[][] read(RootedTreeEdge p, long mc) {
        double[] lo = new double[(int) mc];
        double[] up = new double[(int) mc];
        for (int i = 0; i < mc; i++) {
            lo[i] = p.getLogZLower(i);
            up[i] = p.getLogZUpper(i);
        }
        return new double[][]{lo, up};
    }

    private double[][] runJava(int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions) {
        System.setProperty("branchmarkstar.dp.gpu", "false");
        System.setProperty("branchmarkstar.dp.nativeKernel", "false");
        System.setProperty("branchmarkstar.dp.parallel", "false");
        System.setProperty("branchmarkstar.dp.progress", "false");
        // Legacy (non-fold) index path: INDEPENDENT child-index from the GPU's
        // ChildFoldPlan, so agreement cross-validates the fold/CSR mapping.
        System.setProperty("branchmarkstar.dp.foldChildren", "false");
        RootedTreeEdge p = buildParent(cards, mPos, lambdaPos, childMPositions);
        p.computeFullDP();
        return read(p, product(cards, mPos));
    }

    private double[][] runGpu(int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions, boolean multiGpu) {
        return runGpu(cards, mPos, lambdaPos, childMPositions, multiGpu, false, 1, Integer.MAX_VALUE);
    }

    private double[][] runGpu(int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions,
                              boolean multiGpu, boolean shardedParent, int shardSize, long outputTileMStates) {
        System.setProperty("branchmarkstar.dp.gpu", "true");
        System.setProperty("branchmarkstar.dp.gpu.minWork", "1");
        System.setProperty("branchmarkstar.dp.gpu.maxBytes", String.valueOf(6L * 1024 * 1024 * 1024));
        System.setProperty("branchmarkstar.dp.gpu.trace", "true");
        System.setProperty("branchmarkstar.dp.gpu.outputTileMStates", String.valueOf(outputTileMStates));
        System.setProperty("branchmarkstar.dp.gpu.persistentContext", "true");
        System.setProperty("branchmarkstar.dp.nativeKernel", "false");
        System.setProperty("branchmarkstar.dp.parallel", "false");
        System.setProperty("branchmarkstar.dp.foldChildren", "true");
        if (multiGpu) {
            // Force the M-state split across all visible GPUs + log the multi-gpu line.
            System.setProperty("branchmarkstar.dp.progress", "true");
            System.setProperty("branchmarkstar.dp.gpu.minMStatesPerGpu", "1");
            System.setProperty("branchmarkstar.dp.gpu.maxGpus", "0");
        } else {
            System.setProperty("branchmarkstar.dp.progress", "false");
            System.setProperty("branchmarkstar.dp.gpu.minMStatesPerGpu", String.valueOf(Long.MAX_VALUE));
            System.setProperty("branchmarkstar.dp.gpu.maxGpus", "1");
        }
        RootedTreeEdge p = buildParent(cards, mPos, lambdaPos, childMPositions, shardedParent, shardSize);
        invoke(p, "ensureChildFoldPlans");
        boolean fired = (Boolean) invoke(p, "tryComputeFullDPGpu");
        Assumptions.assumeTrue(fired,
                "GPU full-DP path did not fire (no CUDA GPU on this node?) - skipping numeric comparison");
        return read(p, product(cards, mPos));
    }

    private void compareShape(String name, int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions) {
        compareShape(name, cards, mPos, lambdaPos, childMPositions, false);
    }

    private void compareShape(String name, int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions, boolean multiGpu) {
        compareShape(name, cards, mPos, lambdaPos, childMPositions, multiGpu, false, 1, Integer.MAX_VALUE);
    }

    private void compareShape(String name, int[] cards, int[] mPos, int[] lambdaPos, int[][] childMPositions,
                              boolean multiGpu, boolean shardedParent, int shardSize, long outputTileMStates) {
        double[][] gpu = runGpu(cards, mPos, lambdaPos, childMPositions, multiGpu,
                shardedParent, shardSize, outputTileMStates);   // may skip if no GPU
        double[][] jav = runJava(cards, mPos, lambdaPos, childMPositions);
        int n = jav[0].length;

        double maxAbsLo = 0, maxRelLo = 0, maxAbsUp = 0, maxRelUp = 0;
        for (int i = 0; i < n; i++) {
            assertTrue(Double.isFinite(gpu[0][i]) && Double.isFinite(gpu[1][i]),
                    name + ": GPU produced non-finite bound at mIdx=" + i);
            double aLo = Math.abs(gpu[0][i] - jav[0][i]);
            double aUp = Math.abs(gpu[1][i] - jav[1][i]);
            maxAbsLo = Math.max(maxAbsLo, aLo);
            maxAbsUp = Math.max(maxAbsUp, aUp);
            maxRelLo = Math.max(maxRelLo, aLo / (1.0 + Math.abs(jav[0][i])));
            maxRelUp = Math.max(maxRelUp, aUp / (1.0 + Math.abs(jav[1][i])));
        }

        System.out.println(String.format(Locale.ROOT,
                "[GPU-DP-VALIDATE] %s children=%d mStates=%d lambdaStates=%d | lower maxAbs=%.3e maxRel=%.3e | upper maxAbs=%.3e maxRel=%.3e | java lo[0]=%.10f up[0]=%.10f gpu lo[0]=%.10f up[0]=%.10f",
                name, childMPositions.length, n, product(cards, lambdaPos),
                maxAbsLo, maxRelLo, maxAbsUp, maxRelUp,
                jav[0][0], jav[1][0], gpu[0][0], gpu[1][0]));

        assertTrue(n > 1 && jav[0][0] != jav[0][n - 1], name + ": Java bounds do not vary across M-states");
        assertTrue(maxRelLo <= REL_TOL, name + ": lower bound rel diff " + maxRelLo + " > " + REL_TOL);
        assertTrue(maxRelUp <= REL_TOL, name + ": upper bound rel diff " + maxRelUp + " > " + REL_TOL);
    }

    @Test
    public void gpuMatchesJava_singleChild_small() {
        compareShape("1child-small", new int[]{3, 3, 4, 4}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}});
    }

    @Test
    public void gpuMatchesJava_singleChild_larger() {
        compareShape("1child-larger", new int[]{8, 8, 16, 16}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}});
    }

    @Test
    public void gpuMatchesJava_threeChildren_mixed() {
        // Mirrors TestChildFold: one lambda-dependent child + two lambda-invariant.
        // {1,2}=M+lambda (dep), {0,1}=subset of parent M (invariant), {0}=invariant.
        compareShape("3child-mixed", new int[]{3, 3, 4, 4}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}, {0, 1}, {0}});
    }

    @Test
    public void gpuMatchesJava_twoChildren_bothLambdaDependent() {
        // {1,2} and {0,3}: both mix one M and one lambda position -> both lambda-dependent.
        compareShape("2child-dep", new int[]{4, 4, 5, 5}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}, {0, 3}});
    }

    @Test
    public void gpuMatchesJava_shardedParent_outputTiling() {
        // Parent output is shard-backed and forced through multiple tiny GPU
        // output tiles; child tables remain dense, which is the current GPU gate.
        compareShape("sharded-parent-tiles", new int[]{5, 6, 7, 8}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}, {0, 3}}, false, true, 11, 7);
    }

    // ---- single-node multi-GPU M-state split (requires >=2 visible GPUs) ----

    @Test
    public void gpuMultiGpu_matchesJava_singleChild() {
        // 100 M-states split across all visible GPUs; M-states are independent so
        // the gathered table must equal the single-threaded Java path exactly.
        compareShape("multigpu-1child", new int[]{10, 10, 8, 8}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}}, true);
    }

    @Test
    public void gpuMultiGpu_matchesJava_threeChildren() {
        compareShape("multigpu-3child", new int[]{10, 10, 8, 8}, new int[]{0, 1}, new int[]{2, 3},
                new int[][]{{1, 2}, {0, 1}, {0}}, true);
    }
}
