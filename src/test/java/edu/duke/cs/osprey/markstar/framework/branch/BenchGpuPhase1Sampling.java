package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.gpu.cuda.Gpus;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Throughput micro-benchmark for the PAC Phase-1 GPU sampler, answering the open
 * question "which sampling method performs better": GUMBEL (full per-sample
 * argmax-Gumbel scan) vs CDF_DEDUP (one softmax CDF per distinct parent mIdx,
 * reused by binary search across all samples sharing that mIdx).
 *
 * The decisive variable is reuse = samples / distinctMIdx. cdfDedup amortizes
 * its O(distinctMIdx x lambdaStates) CDF build over O(samples x log lambda)
 * draws, so it should win as reuse grows; gumbel does O(samples x lambdaStates)
 * regardless. Real PAC Phase-1 has distinctMIdx << samples, so we sweep reuse.
 *
 * CPU references (mirrors of the kernels, same as TestGpuPhase1Sampling) anchor
 * whether the GPU is worth it at all on this fp64-ish workload.
 *
 * Not a unit test (no assert beyond GPU-availability skip); run explicitly:
 *   ./gradlew test --tests ...BenchGpuPhase1Sampling -DtestMaxHeap=8g
 */
public class BenchGpuPhase1Sampling {

    // synthetic edge shape (kept self-consistent across all tables)
    private static final int LA = 64, LB = 64;        // lambda factors
    private static final int TOTAL_LAMBDA = LA * LB;  // 4096 lambda states
    private static final int MA = 64, MB = 64;        // M factors
    private static final long M_PRODUCT = (long) MA * MB; // distinct parent M states
    private static final int BLOCK_THREADS = 256;
    private static final double INV_RT = 1.25;

    private static SamplingGpuPhase1.Request buildRequest(int samples, int distinctMIdx) {
        SamplingGpuPhase1.Request req = new SamplingGpuPhase1.Request();
        req.totalLambdaStates = TOTAL_LAMBDA;
        req.blockThreads = BLOCK_THREADS;
        req.numChildren = 2;
        req.invRT = INV_RT;
        req.mCounts = new int[] { MA, MB };
        req.lambdaCounts = new int[] { LA, LB };

        req.lambdaOnlyMin = new double[TOTAL_LAMBDA];
        for (int i = 0; i < TOTAL_LAMBDA; i++) {
            req.lambdaOnlyMin[i] = -0.013 * i + 0.07 * (i % 11) - 0.002 * (i % 97);
        }

        // two lm-pairs: (lamSlot0,mSlot0) and (lamSlot1,mSlot1), each block [lambda][m]
        req.lmLamSlots = new int[] { 0, 1 };
        req.lmMSlots   = new int[] { 0, 1 };
        req.lmMCounts  = new int[] { MA, MB };
        int block0 = LA * MA, block1 = LB * MB;
        req.lmOffsets  = new long[] { 0, block0 };
        req.lmMin = new double[block0 + block1];
        for (int i = 0; i < req.lmMin.length; i++) {
            req.lmMin[i] = 0.11 * Math.sin(i + 0.25) - 0.07 * Math.cos(0.5 * i);
        }

        // two children: child0 indexes [mSlot0][lamSlot0], child1 [mSlot1][lamSlot1]
        req.childMSrcAll    = new int[]  { 0, 1 };
        req.childMStrideAll = new long[] { LA, LB };
        req.childMTermOff   = new int[]  { 0, 1 };
        req.childMTermCnt   = new int[]  { 1, 1 };
        req.childLSrcAll    = new int[]  { 0, 1 };
        req.childLStrideAll = new long[] { 1, 1 };
        req.childLTermOff   = new int[]  { 0, 1 };
        req.childLTermCnt   = new int[]  { 1, 1 };
        int childSize0 = MA * LA, childSize1 = MB * LB;
        req.childTableBase  = new long[] { 0, childSize0 };
        req.childUpperAll = new double[childSize0 + childSize1];
        for (int i = 0; i < req.childUpperAll.length; i++) {
            req.childUpperAll[i] = 0.019 * i - 0.04 * (i % 7) + 0.5 * Math.sin(0.001 * i);
        }
        req.childUpperCacheKey = 0x6142d97f2a5b6c31L ^ ((long) samples << 20) ^ distinctMIdx;

        // samples: distinctMIdx distinct parent M states, round-robin -> reuse = samples/distinctMIdx
        req.mIdxPerSample = new long[samples];
        int d = Math.max(1, Math.min(distinctMIdx, (int) Math.min(M_PRODUCT, Integer.MAX_VALUE)));
        for (int s = 0; s < samples; s++) {
            req.mIdxPerSample[s] = Math.floorMod(splitmix64(0x5eed1234L + (s % d)), M_PRODUCT);
        }
        req.baseSeed = 0x1234abcd55aa7711L;
        req.progress = false;
        req.multiGpu = false;          // single GPU primary; multi-GPU overhead dominates small work
        req.maxGpus = 1;
        req.minGroupsPerGpu = 1;
        req.persistentContext = true;
        req.residentChildTables = true;
        req.method = SamplingGpuPhase1.Method.GUMBEL;
        return req;
    }

    private static double timeGpu(SamplingGpuPhase1.Method method, int samples, int distinctMIdx,
                                  int iters, int[] outFired) {
        double min = Double.MAX_VALUE;
        boolean fired = false;
        for (int it = 0; it < iters; it++) {
            SamplingGpuPhase1.Request req = buildRequest(samples, distinctMIdx);
            req.method = method;
            long t0 = System.nanoTime();
            int[] out = SamplingGpuPhase1.sample(req);
            double ms = (System.nanoTime() - t0) / 1e6;
            if (out == null) { outFired[0] = 0; return Double.NaN; }
            fired = true;
            min = Math.min(min, ms);
        }
        outFired[0] = fired ? 1 : 0;
        return min;
    }

    private static double timeCpu(boolean cdf, int samples, int distinctMIdx, int iters) {
        double min = Double.MAX_VALUE;
        for (int it = 0; it < iters; it++) {
            SamplingGpuPhase1.Request req = buildRequest(samples, distinctMIdx);
            long t0 = System.nanoTime();
            int[] out = cdf ? cpuCdfReference(req) : cpuGumbelReference(req);
            double ms = (System.nanoTime() - t0) / 1e6;
            if (out.length != samples) throw new IllegalStateException("bad cpu out");
            min = Math.min(min, ms);
        }
        return min;
    }

    @Test
    public void benchmark() {
        Assumptions.assumeTrue(hasCudaGpu(), "no CUDA GPU available");
        SamplingGpuPhase1.resetForTesting();

        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println(String.format(Locale.ROOT,
                "[GPU-SAMPLE-BENCH] edge totalLambda=%d mProduct=%d children=2 cores=%d",
                TOTAL_LAMBDA, M_PRODUCT, cores));
        System.out.println("[GPU-SAMPLE-BENCH] timing = min over iters of full sample() wall (ms), single GPU");

        int samples = 4000;
        int[] reuseDistinct = { 1, 16, 256, 4000 }; // distinctMIdx; reuse = samples/distinctMIdx
        int iters = 4;

        // warmup both GPU paths + CPU
        int[] fired = new int[1];
        timeGpu(SamplingGpuPhase1.Method.GUMBEL, 256, 16, 1, fired);
        timeGpu(SamplingGpuPhase1.Method.CDF_DEDUP, 256, 16, 1, fired);
        timeCpu(false, 256, 16, 1);
        timeCpu(true, 256, 16, 1);

        System.out.println(String.format(Locale.ROOT,
                "%-10s %-9s %-12s %-12s %-12s %-12s %-12s",
                "samples", "distinct", "reuse", "gpuGumbel", "gpuCdfDedup", "cpuGumbel", "cpuCdfDedup"));
        for (int distinct : reuseDistinct) {
            double reuse = (double) samples / distinct;
            double gG = timeGpu(SamplingGpuPhase1.Method.GUMBEL, samples, distinct, iters, fired);
            double gC = timeGpu(SamplingGpuPhase1.Method.CDF_DEDUP, samples, distinct, iters, fired);
            double cG = timeCpu(false, samples, distinct, iters);
            double cC = timeCpu(true, samples, distinct, iters);
            System.out.println(String.format(Locale.ROOT,
                    "%-10d %-9d %-12.1f %-12.2f %-12.2f %-12.2f %-12.2f",
                    samples, distinct, reuse, gG, gC, cG, cC));
        }
        System.out.println("[GPU-SAMPLE-BENCH] interpretation: cdfDedup should pull ahead as reuse (samples/distinct) grows.");
    }

    private static boolean hasCudaGpu() {
        try { return !Gpus.get().getGpus().isEmpty(); }
        catch (Throwable t) { return false; }
    }

    // ---- CPU references (mirror TestGpuPhase1Sampling kernels) -------------

    private static int[] cpuGumbelReference(SamplingGpuPhase1.Request req) {
        int[] out = new int[req.mIdxPerSample.length];
        int[] mLocal = new int[req.mCounts.length];
        int[] lambdaLocal = new int[req.lambdaCounts.length];
        for (int s = 0; s < out.length; s++) {
            decodeState(req.mIdxPerSample[s], req.mCounts, mLocal);
            double bestScore = Double.NEGATIVE_INFINITY;
            int bestIdx = req.totalLambdaStates - 1;
            for (int lIdx = 0; lIdx < req.totalLambdaStates; lIdx++) {
                double eMin = localEnergy(req, lIdx, mLocal, lambdaLocal);
                double fUpper = childSum(req, mLocal, lambdaLocal);
                double logWeight = -eMin * req.invRT + fUpper;
                if (Double.isNaN(logWeight)) continue;
                double score = logWeight + gumbelNoise(req.baseSeed, s, lIdx);
                if (score > bestScore) { bestScore = score; bestIdx = lIdx; }
            }
            out[s] = bestIdx;
        }
        return out;
    }

    private static int[] cpuCdfReference(SamplingGpuPhase1.Request req) {
        int[] out = new int[req.mIdxPerSample.length];
        Map<Long, double[]> cdfs = new HashMap<>();
        for (int s = 0; s < out.length; s++) {
            long mIdx = req.mIdxPerSample[s];
            double[] cdf = cdfs.computeIfAbsent(mIdx, key -> buildCdf(req, key));
            out[s] = sampleFromCdf(cdf, uniform01(req.baseSeed, s));
        }
        return out;
    }

    private static double[] buildCdf(SamplingGpuPhase1.Request req, long mIdx) {
        double[] cdf = new double[req.totalLambdaStates];
        int[] mLocal = new int[req.mCounts.length];
        int[] lambdaLocal = new int[req.lambdaCounts.length];
        decodeState(mIdx, req.mCounts, mLocal);
        double maxLog = Double.NEGATIVE_INFINITY;
        for (int lIdx = 0; lIdx < req.totalLambdaStates; lIdx++) {
            double eMin = localEnergy(req, lIdx, mLocal, lambdaLocal);
            double fUpper = childSum(req, mLocal, lambdaLocal);
            double logWeight = -eMin * req.invRT + fUpper;
            cdf[lIdx] = logWeight;
            if (!Double.isNaN(logWeight) && logWeight > maxLog) maxLog = logWeight;
        }
        if (!Double.isFinite(maxLog)) { Arrays.fill(cdf, 0.0); return cdf; }
        double running = 0.0;
        for (int lIdx = 0; lIdx < req.totalLambdaStates; lIdx++) {
            double w = Math.exp(cdf[lIdx] - maxLog);
            if (!Double.isFinite(w)) w = 0.0;
            running += w;
            cdf[lIdx] = running;
        }
        return cdf;
    }

    private static int sampleFromCdf(double[] cdf, double u) {
        double total = cdf[cdf.length - 1];
        if (!(total > 0.0) || !Double.isFinite(total)) return cdf.length - 1;
        double target = u * total;
        int lo = 0, hi = cdf.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (target < cdf[mid]) hi = mid; else lo = mid + 1;
        }
        return lo;
    }

    private static double localEnergy(SamplingGpuPhase1.Request req, int lIdx, int[] mLocal, int[] lambdaLocal) {
        decodeState(lIdx, req.lambdaCounts, lambdaLocal);
        double e = req.lambdaOnlyMin[lIdx];
        for (int p = 0; p < req.lmLamSlots.length; p++) {
            int lamSlot = req.lmLamSlots[p];
            int mSlot = req.lmMSlots[p];
            long off = req.lmOffsets[p]
                    + (long) lambdaLocal[lamSlot] * (long) req.lmMCounts[p]
                    + (long) mLocal[mSlot];
            e += req.lmMin[(int) off];
        }
        return e;
    }

    private static double childSum(SamplingGpuPhase1.Request req, int[] mLocal, int[] lambdaLocal) {
        double sum = 0.0;
        for (int c = 0; c < req.numChildren; c++) {
            long idx = 0L;
            for (int t = 0; t < req.childMTermCnt[c]; t++) {
                int off = req.childMTermOff[c] + t;
                idx += (long) mLocal[req.childMSrcAll[off]] * req.childMStrideAll[off];
            }
            for (int t = 0; t < req.childLTermCnt[c]; t++) {
                int off = req.childLTermOff[c] + t;
                idx += (long) lambdaLocal[req.childLSrcAll[off]] * req.childLStrideAll[off];
            }
            sum += req.childUpperAll[(int) (req.childTableBase[c] + idx)];
        }
        return sum;
    }

    private static void decodeState(long idx, int[] counts, int[] out) {
        for (int i = counts.length - 1; i >= 0; i--) { out[i] = (int) (idx % counts[i]); idx /= counts[i]; }
    }
    private static double uniform01(long baseSeed, int sampleSlot) {
        long h = splitmix64(baseSeed + 0x9E3779B97F4A7C15L * (long) sampleSlot);
        return ((double) (h >>> 11) + 0.5) * (1.0 / 9007199254740992.0);
    }
    private static double gumbelNoise(long baseSeed, int sampleSlot, int lIdx) {
        long h = splitmix64(baseSeed + 0x9E3779B97F4A7C15L * (long) sampleSlot);
        h = splitmix64(h + (long) lIdx);
        double u = ((double) (h >>> 11) + 0.5) * (1.0 / 9007199254740992.0);
        return -Math.log(-Math.log(u));
    }
    private static long splitmix64(long x) {
        x += 0x9E3779B97F4A7C15L;
        x = (x ^ (x >>> 30)) * 0xBF58476D1CE4E5B9L;
        x = (x ^ (x >>> 27)) * 0x94D049BB133111EBL;
        return x ^ (x >>> 31);
    }
}
