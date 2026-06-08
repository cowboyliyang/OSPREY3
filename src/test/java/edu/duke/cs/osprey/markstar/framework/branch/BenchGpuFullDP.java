package edu.duke.cs.osprey.markstar.framework.branch;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.LinkedHashSet;
import java.util.Locale;

/**
 * Throughput micro-benchmark for the GPU full-DP fast path vs the Java DP path,
 * on one large synthetic non-leaf edge whose per-element compute (energy gather
 * over lm-pairs + multi-child fold + two-pass log-sum-exp) matches the real
 * kernel, so the measured speed ratio transfers to real cases. Answers the open
 * question: is the fp64 GPU kernel actually faster than the 32-thread Java path
 * on A5000 (weak fp64), and how much do multi-GPU / Titan V add.
 *
 * Not a unit test (no behavioural assert beyond a correctness cross-check at
 * scale); run explicitly:
 *   ./gradlew test --tests ...BenchGpuFullDP -DtestMaxHeap=8g
 */
public class BenchGpuFullDP {

    // ---- helpers (mirror TestGpuFullDP) -----------------------------------

    private static void set(Object o, String field, Object val) {
        try { Field f = RootedTreeEdge.class.getDeclaredField(field); f.setAccessible(true); f.set(o, val); }
        catch (ReflectiveOperationException e) { throw new RuntimeException("set " + field, e); }
    }
    private static Object invoke(Object o, String method) {
        try { Method m = RootedTreeEdge.class.getDeclaredMethod(method); m.setAccessible(true); return m.invoke(o); }
        catch (ReflectiveOperationException e) { throw new RuntimeException("invoke " + method, e); }
    }
    private static RCs rcsFor(int[] positions, int[] cards) {
        int maxPos = 0; for (int p : positions) maxPos = Math.max(maxPos, p);
        int[][] arr = new int[maxPos + 1][];
        for (int i = 0; i < arr.length; i++) arr[i] = new int[1];
        for (int i = 0; i < positions.length; i++) arr[positions[i]] = new int[cards[i]];
        return new RCs(arr);
    }
    private static EnergyMatrix ematFilled(int numPos, int[] cards, double seed) {
        EnergyMatrix e = new EnergyMatrix(numPos, cards, 0.0);
        for (int p = 0; p < numPos; p++)
            for (int rc = 0; rc < cards[p]; rc++)
                e.setOneBody(p, rc, seed + 0.13 * p - 0.07 * rc);
        for (int pi = 0; pi < numPos; pi++)
            for (int pj = pi + 1; pj < numPos; pj++)
                for (int ri = 0; ri < cards[pi]; ri++)
                    for (int rj = 0; rj < cards[pj]; rj++)
                        e.setPairwise(pi, ri, pj, rj, 0.31 * seed + 0.05 * (ri + 1) * (rj + 2) - 0.02 * (pi + pj));
        return e;
    }
    private static InteractionGraph fullyConnected(int numPos) {
        boolean[][] adj = new boolean[numPos][numPos];
        for (int i = 0; i < numPos; i++) for (int j = 0; j < numPos; j++) if (i != j) adj[i][j] = true;
        try {
            Constructor<InteractionGraph> c = InteractionGraph.class.getDeclaredConstructor(int.class, boolean[][].class);
            c.setAccessible(true); return c.newInstance(numPos, adj);
        } catch (ReflectiveOperationException e) { throw new RuntimeException("InteractionGraph", e); }
    }
    private static long product(int[] cards, int[] pos) { long c = 1; for (int p : pos) c *= cards[p]; return c; }

    private static RootedTreeEdge child(RCs rcs, int[] mPos, int[] cards, double base) {
        RootedTreeEdge c = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(c, "mPositionsSorted", mPos.clone());
        set(c, "isLambdaEdge", Boolean.TRUE);
        long mc = product(cards, mPos);
        set(c, "mStateCount", Long.valueOf(mc));
        set(c, "mArraySize", Integer.valueOf((int) mc));
        DenseDPTable t = new DenseDPTable(mc);
        for (long i = 0; i < mc; i++) t.set(i, base + 0.5 * (i % 97), base + 50.0 + 0.3 * (i % 89));
        set(c, "dpTable", t);
        return c;
    }

    // bench edge: mStates=27^4=531441, lambdaStates=45^2=2025, 2 lambda-dependent children
    private static final int[] CARDS = {27, 27, 27, 27, 45, 45};
    private static final int[] MPOS = {0, 1, 2, 3};
    private static final int[] LPOS = {4, 5};
    private static final int[][] CHILDREN = {{0, 1, 4}, {2, 3, 5}};

    private RootedTreeEdge buildParent() {
        int numPos = CARDS.length;
        int[] allPos = new int[numPos];
        for (int i = 0; i < numPos; i++) allPos[i] = i;
        RCs rcs = rcsFor(allPos, CARDS);
        LinkedHashSet<RootedTreeEdge> fset = new LinkedHashSet<>();
        for (int c = 0; c < CHILDREN.length; c++) fset.add(child(rcs, CHILDREN[c], CARDS, 1000.0 * (c + 1)));
        RootedTreeEdge p = new RootedTreeEdge(null, null, new LinkedHashSet<>(), false, rcs);
        set(p, "isLambdaEdge", Boolean.TRUE);
        set(p, "mPositionsSorted", MPOS.clone());
        set(p, "lambdaPositionsSorted", LPOS.clone());
        long mc = product(CARDS, MPOS);
        set(p, "mStateCount", Long.valueOf(mc));
        set(p, "mArraySize", Integer.valueOf((int) mc));
        set(p, "totalLambdaStates", Integer.valueOf((int) product(CARDS, LPOS)));
        set(p, "dpTable", new DenseDPTable(mc));
        set(p, "Fset", fset);
        p.initIncrementalEnumeration(ematFilled(numPos, CARDS, 1.0), ematFilled(numPos, CARDS, 2.0),
                fullyConnected(numPos), 1.9);
        return p;
    }

    private static void baseGpuProps() {
        System.setProperty("branchmarkstar.dp.gpu.minWork", "1");
        System.setProperty("branchmarkstar.dp.gpu.maxBytes", String.valueOf(20L * 1024 * 1024 * 1024));
        System.setProperty("branchmarkstar.dp.gpu.trace", "true");
        System.setProperty("branchmarkstar.dp.nativeKernel", "false");
        System.setProperty("branchmarkstar.dp.foldChildren", "true");
        System.setProperty("branchmarkstar.dp.progress", "true");
    }

    /** Returns {minElapsedMs, lo0, up0, lastLo, lastUp} for correctness cross-check. */
    private double[] timeJava(int threads, int iters) {
        System.setProperty("branchmarkstar.dp.gpu", "false");
        System.setProperty("branchmarkstar.dp.parallel", threads > 1 ? "true" : "false");
        System.setProperty("branchmarkstar.dp.parallel.threads", String.valueOf(Math.max(1, threads)));
        System.setProperty("branchmarkstar.dp.parallel.minMStates", "1");
        System.setProperty("branchmarkstar.dp.progress", "false");
        double min = Double.MAX_VALUE; double lo0 = 0, up0 = 0;
        for (int it = 0; it < iters; it++) {
            RootedTreeEdge p = buildParent();
            long t0 = System.nanoTime();
            p.computeFullDP();
            double ms = (System.nanoTime() - t0) / 1e6;
            min = Math.min(min, ms);
            lo0 = p.getLogZLower(0); up0 = p.getLogZUpper(0);
        }
        return new double[]{min, lo0, up0};
    }

    /** maxGpus<=0 => all visible. Returns {minElapsedMs, lo0, up0, firedFlag}. */
    private double[] timeGpu(int maxGpus, int iters) {
        baseGpuProps();
        System.setProperty("branchmarkstar.dp.gpu", "true");
        System.setProperty("branchmarkstar.dp.parallel", "false");
        if (maxGpus <= 0) {
            System.setProperty("branchmarkstar.dp.gpu.maxGpus", "0");
            System.setProperty("branchmarkstar.dp.gpu.minMStatesPerGpu", "1");
        } else {
            System.setProperty("branchmarkstar.dp.gpu.maxGpus", String.valueOf(maxGpus));
            System.setProperty("branchmarkstar.dp.gpu.minMStatesPerGpu", String.valueOf(Long.MAX_VALUE));
        }
        double min = Double.MAX_VALUE; double lo0 = 0, up0 = 0; boolean fired = false;
        for (int it = 0; it < iters; it++) {
            RootedTreeEdge p = buildParent();
            invoke(p, "ensureChildFoldPlans");
            long t0 = System.nanoTime();
            boolean ok = (Boolean) invoke(p, "tryComputeFullDPGpu");
            double ms = (System.nanoTime() - t0) / 1e6;
            fired = ok;
            if (!ok) break;
            min = Math.min(min, ms);
            lo0 = p.getLogZLower(0); up0 = p.getLogZUpper(0);
        }
        return new double[]{min, lo0, up0, fired ? 1 : 0};
    }

    @Test
    public void benchmark() {
        long work = product(CARDS, MPOS) * product(CARDS, LPOS);
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println(String.format(Locale.ROOT,
                "[GPU-DP-BENCH] edge mStates=%d lambdaStates=%d children=%d work=%.3e cores=%d",
                product(CARDS, MPOS), product(CARDS, LPOS), CHILDREN.length, (double) work, cores));

        // warmup (JIT) + timed
        timeJava(1, 1);
        double[] j1 = timeJava(1, 3);
        double[] jN = timeJava(cores, 3);
        double[] g1 = timeGpu(1, 4);
        double[] gA = timeGpu(0, 4);

        System.out.println(String.format(Locale.ROOT, "[GPU-DP-BENCH] java-1thread    minMs=%.1f", j1[0]));
        System.out.println(String.format(Locale.ROOT, "[GPU-DP-BENCH] java-%dthread   minMs=%.1f  vsJava1=%.2fx", cores, jN[0], j1[0] / jN[0]));
        if (g1[3] == 1) {
            System.out.println(String.format(Locale.ROOT, "[GPU-DP-BENCH] gpu-1          minMs=%.1f  vsJava1=%.2fx  vsJava%d=%.2fx", g1[0], j1[0] / g1[0], cores, jN[0] / g1[0]));
        } else {
            System.out.println("[GPU-DP-BENCH] gpu-1          NOT AVAILABLE (no CUDA GPU)");
        }
        if (gA[3] == 1) {
            System.out.println(String.format(Locale.ROOT, "[GPU-DP-BENCH] gpu-all        minMs=%.1f  vsJava1=%.2fx  vsJava%d=%.2fx", gA[0], j1[0] / gA[0], cores, jN[0] / gA[0]));
        } else {
            System.out.println("[GPU-DP-BENCH] gpu-all        NOT AVAILABLE (no CUDA GPU)");
        }

        // correctness cross-check at scale (GPU vs Java), if GPU ran
        if (g1[3] == 1) {
            assertEquals(j1[1], g1[1], 1e-6 * (1 + Math.abs(j1[1])), "gpu-1 lower[0] vs java");
            assertEquals(j1[2], g1[2], 1e-6 * (1 + Math.abs(j1[2])), "gpu-1 upper[0] vs java");
        }
        assertTrue(Double.isFinite(j1[1]));
    }
}
