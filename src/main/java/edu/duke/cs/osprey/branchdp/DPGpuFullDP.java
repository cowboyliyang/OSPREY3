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
*/

package edu.duke.cs.osprey.branchdp;

import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Context;
import edu.duke.cs.osprey.gpu.cuda.Gpu;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Gpus;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import jcuda.Pointer;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUresult;
import jcuda.driver.JCudaDriver;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.util.ArrayList;
import java.util.List;

/**
 * CUDA fast path for branch-DP full-DP: one non-leaf edge with one or more
 * dense child DP tables and full lambda enumeration. Parent output is written
 * tile-by-tile into DPTable, so it may be dense or shard-backed.
 * The caller owns shape checks and falls back to the Java DP path otherwise.
 */
final class DPGpuFullDP {

    static final int MAX_EDGE_POSITIONS = 32;

    // When an edge's GPU footprint exceeds device memory, fail the whole run instead
    // of silently falling back to the (single-threaded, hours-long) Java DP path.
    private static final String DP_GPU_FAIL_IF_EXCEEDS_VRAM = "branchdp.dp.gpu.failIfExceedsVram";
    private static final String DP_GPU_VRAM_HEADROOM_PERCENT = "branchdp.dp.gpu.vramHeadroomPercent";
    private static final int DEFAULT_DP_GPU_VRAM_HEADROOM_PERCENT = 92;

    /** Edge GPU footprint exceeds device memory; propagated as fatal (no CPU fallback). */
    static final class GpuMemoryExceededException extends RuntimeException {
        GpuMemoryExceededException(String message) {
            super(message);
        }
    }

    private static volatile boolean unavailable = false;
    private static volatile boolean unavailableLogged = false;
    // Per-edge (non-sticky) fallback: a single oversized/failed edge falls back to
    // Java without disabling GPU for the whole run; only a genuinely broken GPU
    // (MAX_CONSECUTIVE_EDGE_FAILURES in a row) flips the sticky unavailable flag.
    private static volatile int consecutiveEdgeFailures = 0;
    private static final int MAX_CONSECUTIVE_EDGE_FAILURES = 5;
    private static final Object persistentLock = new Object();
    private static PersistentPool persistentPool = null;
    private static boolean persistentShutdownHook = false;

    private DPGpuFullDP() {
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        return BranchDpConfig.getBackendBoolean(key, defaultValue);
    }

    private static int getConfigInteger(String key, int defaultValue) {
        return BranchDpConfig.getBackendInteger(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    /** True if t (or any cause) is a device-OOM / footprint-too-big condition. */
    private static boolean isGpuMemoryFailure(Throwable t) {
        for (Throwable c = t; c != null; c = c.getCause()) {
            if (c instanceof GpuMemoryExceededException) {
                return true;
            }
            String m = c.getMessage();
            if (m != null && m.contains("OUT_OF_MEMORY")) {
                return true;
            }
        }
        return false;
    }

    private static String rootMessage(Throwable t) {
        Throwable c = t;
        while (c.getCause() != null) {
            c = c.getCause();
        }
        return c.getMessage();
    }

    static final class Request {
        long mStateCount;     // global parent M-state count
        long mStateChunk;     // max M-states per output tile/kernel launch (= numBlocks)
        int totalLambdaStates;
        int blockThreads;
        int numChildren;
        double invRT;

        int[] mCounts;
        int[] lambdaCounts;

        double[] lambdaOnlyRigid;
        double[] lambdaOnlyMin;
        double[] lmRigid;
        double[] lmMin;
        int[] lmLamSlots;
        int[] lmMSlots;
        int[] lmMCounts;
        long[] lmOffsets;

        // Per-child fold plans, concatenated CSR-style across numChildren.
        int[] childMSrcAll;
        long[] childMStrideAll;
        int[] childMTermOff;
        int[] childMTermCnt;
        int[] childLSrcAll;
        long[] childLStrideAll;
        int[] childLTermOff;
        int[] childLTermCnt;
        long[] childTableBase;     // per-child element offset into the device buffer (long: may exceed int)
        long childTableTotal;      // total element count across all child tables (may exceed Integer.MAX_VALUE)
        double[][] childLowerParts; // per-child dense lower tables; never concatenated on the Java heap
        double[][] childUpperParts;

        DPTable outTable;     // parent table written tile-by-tile

        boolean progress;
        long work;
        long estimatedDeviceBytes;
        boolean persistentContext;

        // Multi-GPU (single-node) M-state split controls.
        boolean multiGpu;
        int maxGpus;          // 0 = all available
        long minMStatesPerGpu;
    }

    static synchronized boolean compute(Request req) {
        if (unavailable) {
            return false;
        }

        List<Gpu> gpus;
        try {
            gpus = Gpus.get().getGpus();
        } catch (Throwable t) {
            // GPU discovery itself failed => genuinely no usable CUDA; sticky.
            markUnavailable("GPU discovery failed: " + t.getMessage(), t);
            return false;
        }
        if (gpus.isEmpty()) {
            markUnavailable("no CUDA GPU with double support", null);
            return false;
        }

        int nGpus = chooseGpuCount(req, gpus.size());
        GpuExecutor[] executors = null;
        try {
            executors = makeExecutors(req, gpus, nGpus);
            if (nGpus <= 1) {
                runTiled(req, executors[0], 0L, req.mStateCount);
            } else {
                runMultiGpu(req, executors, nGpus);
            }
            consecutiveEdgeFailures = 0;
            return true;
        } catch (Throwable t) {
            // Device out-of-memory: the only CPU fallback for these huge edges is the
            // single-threaded Java DP path, which can take hours. Fail the run fast so
            // the job can be resubmitted on a larger-VRAM GPU (a6000/rtx_pro_6000),
            // unless the operator explicitly opts back into fallback.
            if (getConfigBoolean(DP_GPU_FAIL_IF_EXCEEDS_VRAM, true) && isGpuMemoryFailure(t)) {
                System.out.println(BranchDpConfig.getBackendLogPrefix()
                        + " GPU DP out of device memory; failing fast (no CPU fallback): "
                        + rootMessage(t));
                if (t instanceof RuntimeException) {
                    throw (RuntimeException)t;
                }
                throw new RuntimeException(t);
            }
            // Per-edge failure (e.g. a buffer too big even for chunked upload, or a
            // device OOM on one oversized edge): fall back to Java for THIS edge and
            // keep GPU enabled. Disable globally only after a run of failures.
            consecutiveEdgeFailures++;
            if (consecutiveEdgeFailures >= MAX_CONSECUTIVE_EDGE_FAILURES) {
                markUnavailable("disabled after " + consecutiveEdgeFailures
                        + " consecutive edge failures; last="
                        + t.getClass().getSimpleName() + ": " + t.getMessage(), t);
            } else {
                edgeFallback(t);
            }
            return false;
        } finally {
            cleanupExecutors(executors);
        }
    }

    /** Number of GPUs to split the M-state range across (>=1). */
    private static int chooseGpuCount(Request req, int available) {
        if (!req.multiGpu) {
            return 1;
        }
        int cap = available;
        if (req.maxGpus > 0) {
            cap = Math.min(cap, req.maxGpus);
        }
        if (req.mStateCount > 0 && req.mStateCount < cap) {
            cap = (int)req.mStateCount;
        }
        if (req.minMStatesPerGpu > 0) {
            long byWork = req.mStateCount/req.minMStatesPerGpu;
            cap = (int)Math.min((long)cap, Math.max(1L, byWork));
        }
        return Math.max(1, cap);
    }

    /** Split [0, mStateCount) into nGpus contiguous ranges; each GPU then tiles its range. */
    private static void runMultiGpu(Request req, GpuExecutor[] executors, int nGpus) throws InterruptedException {
        long total = req.mStateCount;
        long base = total/nGpus;
        long rem = total%nGpus;
        Thread[] threads = new Thread[nGpus];
        Throwable[] errors = new Throwable[nGpus];
        long offset = 0;
        long startNanos = System.nanoTime();
        for (int g=0; g<nGpus; g++) {
            long chunk = base + (g < rem ? 1 : 0);
            final long off = offset;
            final long sz = chunk;
            final GpuExecutor executor = executors[g];
            final int gi = g;
            threads[g] = new Thread(() -> {
                try {
                    runTiled(req, executor, off, sz);
                } catch (Throwable t) {
                    errors[gi] = t;
                }
            }, "bms-gpu-dp-" + g);
            threads[g].start();
            offset += chunk;
        }
        for (Thread th : threads) {
            th.join();
        }
        for (Throwable e : errors) {
            if (e != null) {
                throw new RuntimeException("multi-GPU DP failed: " + e.getMessage(), e);
            }
        }
        if (req.progress) {
            double ms = (System.nanoTime() - startNanos)/1e6;
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP multi-gpu done, gpus=" + nGpus
                    + ", mStates=" + total
                    + ", lambdaStates=" + req.totalLambdaStates
                    + ", children=" + req.numChildren
                    + ", outputTileMStates=" + req.mStateChunk
                    + ", persistentContext=" + req.persistentContext
                    + ", elapsedMs=" + String.format(java.util.Locale.ROOT, "%.1f", ms));
        }
    }

    private static void runTiled(Request req, GpuExecutor executor, long rangeStart, long rangeSize) {
        if (rangeSize <= 0) {
            return;
        }
        if (req.mStateChunk <= 0 || req.mStateChunk > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("GPU DP output tile exceeds CUDA grid int limit: " + req.mStateChunk);
        }
        executor.runRange(req, rangeStart, rangeSize);
    }

    private static GpuExecutor[] makeExecutors(Request req, List<Gpu> gpus, int nGpus) {
        if (req.persistentContext) {
            return persistentExecutors(gpus, nGpus);
        }
        GpuExecutor[] execs = new GpuExecutor[nGpus];
        for (int i = 0; i < nGpus; i++) {
            execs[i] = new GpuExecutor(gpus.get(i), true);
        }
        return execs;
    }

    private static GpuExecutor[] persistentExecutors(List<Gpu> gpus, int nGpus) {
        synchronized (persistentLock) {
            if (persistentPool == null || !persistentPool.matches(gpus, nGpus)) {
                if (persistentPool != null) {
                    persistentPool.cleanup();
                }
                persistentPool = new PersistentPool(gpus, nGpus);
                if (!persistentShutdownHook) {
                    Runtime.getRuntime().addShutdownHook(new Thread(DPGpuFullDP::cleanupPersistentPool,
                            "bms-gpu-dp-persistent-cleanup"));
                    persistentShutdownHook = true;
                }
            }
            return persistentPool.executors(nGpus);
        }
    }

    private static void cleanupExecutors(GpuExecutor[] executors) {
        if (executors == null) {
            return;
        }
        for (GpuExecutor executor : executors) {
            if (executor != null) {
                executor.cleanupIfOwned();
            }
        }
    }

    private static void cleanupPersistentPool() {
        synchronized (persistentLock) {
            if (persistentPool != null) {
                persistentPool.cleanup();
                persistentPool = null;
            }
        }
    }

    private static final class PersistentPool {
        private final Gpu[] gpus;
        private final GpuExecutor[] executors;

        PersistentPool(List<Gpu> available, int nGpus) {
            this.gpus = new Gpu[nGpus];
            this.executors = new GpuExecutor[nGpus];
            for (int i = 0; i < nGpus; i++) {
                this.gpus[i] = available.get(i);
                this.executors[i] = new GpuExecutor(this.gpus[i], false);
            }
        }

        boolean matches(List<Gpu> available, int nGpus) {
            if (executors.length < nGpus) {
                return false;
            }
            for (int i = 0; i < nGpus; i++) {
                if (gpus[i] != available.get(i)) {
                    return false;
                }
            }
            return true;
        }

        GpuExecutor[] executors(int nGpus) {
            GpuExecutor[] out = new GpuExecutor[nGpus];
            System.arraycopy(executors, 0, out, 0, nGpus);
            return out;
        }

        void cleanup() {
            for (int i = executors.length - 1; i >= 0; i--) {
                executors[i].cleanupPersistent();
            }
        }
    }

    private static final class GpuExecutor {
        private final Gpu gpu;
        private final Context context;
        private final GpuStream stream;
        private final boolean owned;

        GpuExecutor(Gpu gpu, boolean owned) {
            this.gpu = gpu;
            this.owned = owned;
            this.context = new Context(gpu);
            this.context.attachCurrentThread();
            this.stream = new GpuStream(context);
        }

        synchronized void runRange(Request req, long rangeStart, long rangeSize) {
            context.attachCurrentThread();
            runOnGpuRange(req, stream, rangeStart, rangeSize);
        }

        void cleanupIfOwned() {
            if (owned) {
                cleanup();
            }
        }

        void cleanupPersistent() {
            cleanup();
        }

        private void cleanup() {
            try {
                stream.cleanup();
            } catch (Throwable t) {
                t.printStackTrace(System.err);
            }
            try {
                context.cleanup();
            } catch (Throwable t) {
                t.printStackTrace(System.err);
            }
        }
    }

    private static void runOnGpuRange(Request req, GpuStream stream, long rangeStart, long rangeSize) {
        List<CUBuffer<?>> buffers = new ArrayList<>();
        // Large double inputs (child tables, lm/lambda energy tables) can exceed the
        // 2 GiB hard cap of a single java.nio direct buffer, so they bypass CUBuffer
        // and are uploaded straight into raw cuMemAlloc'd device memory (size is long;
        // the device side / cuMemcpyHtoD ByteCount are long, only NIO was int-capped).
        List<CUdeviceptr> rawBufs = new ArrayList<>();

        try {
            stream.getContext().attachCurrentThread();

            // Proactive VRAM gate: child tables must be fully resident (their indexing
            // is scattered, so they cannot be streamed in slices). If the footprint
            // does not fit this GPU, stop here rather than allocate-then-OOM.
            long[] freeMem = {0L};
            long[] totalMem = {0L};
            if (JCudaDriver.cuMemGetInfo(freeMem, totalMem) == CUresult.CUDA_SUCCESS) {
                long need = estimateDeviceBytes(req);
                int headroomPct = Math.max(1, Math.min(100,
                        getConfigInteger(DP_GPU_VRAM_HEADROOM_PERCENT, DEFAULT_DP_GPU_VRAM_HEADROOM_PERCENT)));
                long usable = (long)((double)freeMem[0]*headroomPct/100.0);
                if (need > usable) {
                    throw new GpuMemoryExceededException("GPU DP footprint ~" + need
                            + " B exceeds usable VRAM " + usable + " B (free=" + freeMem[0]
                            + " B, total=" + totalMem[0] + " B, headroom=" + headroomPct + "%); "
                            + "lambdaStates=" + req.totalLambdaStates
                            + ", childTableElems=" + req.childTableTotal
                            + " -- resubmit on a larger-VRAM GPU (a6000 48G / rtx_pro_6000 96G)");
                }
            }

            int blockThreads = resolveBlockThreads(req, stream);
            Kernel kernel = new Kernel(stream, "dp");
            Kernel.Function func = kernel.makeFunction("full_dp_n_children");
            func.blockThreads = blockThreads;
            func.sharedMemCalc = bt -> 4*bt*Double.BYTES;

            CUBuffer<IntBuffer> mCounts = uploadInts(stream, req.mCounts, buffers);
            CUBuffer<IntBuffer> lambdaCounts = uploadInts(stream, req.lambdaCounts, buffers);
            CUdeviceptr lambdaOnlyRigid = uploadDoublesBig(req.lambdaOnlyRigid, rawBufs);
            CUdeviceptr lambdaOnlyMin = uploadDoublesBig(req.lambdaOnlyMin, rawBufs);
            CUdeviceptr lmRigid = uploadDoublesBig(req.lmRigid, rawBufs);
            CUdeviceptr lmMin = uploadDoublesBig(req.lmMin, rawBufs);
            CUBuffer<IntBuffer> lmLamSlots = uploadInts(stream, req.lmLamSlots, buffers);
            CUBuffer<IntBuffer> lmMSlots = uploadInts(stream, req.lmMSlots, buffers);
            CUBuffer<IntBuffer> lmMCounts = uploadInts(stream, req.lmMCounts, buffers);
            CUBuffer<LongBuffer> lmOffsets = uploadLongs(stream, req.lmOffsets, buffers);
            CUBuffer<IntBuffer> childMSrcAll = uploadInts(stream, req.childMSrcAll, buffers);
            CUBuffer<LongBuffer> childMStrideAll = uploadLongs(stream, req.childMStrideAll, buffers);
            CUBuffer<IntBuffer> childMTermOff = uploadInts(stream, req.childMTermOff, buffers);
            CUBuffer<IntBuffer> childMTermCnt = uploadInts(stream, req.childMTermCnt, buffers);
            CUBuffer<IntBuffer> childLSrcAll = uploadInts(stream, req.childLSrcAll, buffers);
            CUBuffer<LongBuffer> childLStrideAll = uploadLongs(stream, req.childLStrideAll, buffers);
            CUBuffer<IntBuffer> childLTermOff = uploadInts(stream, req.childLTermOff, buffers);
            CUBuffer<IntBuffer> childLTermCnt = uploadInts(stream, req.childLTermCnt, buffers);
            CUBuffer<LongBuffer> childTableBase = uploadLongs(stream, req.childTableBase, buffers);
            CUdeviceptr childLowerAll = uploadChildTables(req.childLowerParts, req.childTableBase,
                    req.childTableTotal, rawBufs);
            CUdeviceptr childUpperAll = uploadChildTables(req.childUpperParts, req.childTableBase,
                    req.childTableTotal, rawBufs);
            int maxTile = (int)Math.min(req.mStateChunk, rangeSize);
            CUBuffer<DoubleBuffer> outLower = makeDoubles(stream, maxTile, buffers);
            CUBuffer<DoubleBuffer> outUpper = makeDoubles(stream, maxTile, buffers);

            long[] mStateCountArg = { req.mStateCount };
            long[] mStateOffsetArg = { rangeStart };
            int[] totalLambdaStatesArg = { req.totalLambdaStates };
            int[] mCountArg = { req.mCounts.length };
            int[] lambdaCountArg = { req.lambdaCounts.length };
            int[] lmPairCountArg = { req.lmLamSlots.length };
            int[] numChildrenArg = { req.numChildren };
            double[] invRTArg = { req.invRT };

            long remaining = rangeSize;
            long offset = rangeStart;
            while (remaining > 0) {
                int chunk = (int)Math.min(remaining, req.mStateChunk);
                func.numBlocks = chunk;
                mStateOffsetArg[0] = offset;

                func.setArgs(Pointer.to(
                    mCounts.getDevicePointer(),
                    lambdaCounts.getDevicePointer(),
                    Pointer.to(lambdaOnlyRigid),
                    Pointer.to(lambdaOnlyMin),
                    Pointer.to(lmRigid),
                    Pointer.to(lmMin),
                    lmLamSlots.getDevicePointer(),
                    lmMSlots.getDevicePointer(),
                    lmMCounts.getDevicePointer(),
                    lmOffsets.getDevicePointer(),
                    childMSrcAll.getDevicePointer(),
                    childMStrideAll.getDevicePointer(),
                    childMTermOff.getDevicePointer(),
                    childMTermCnt.getDevicePointer(),
                    childLSrcAll.getDevicePointer(),
                    childLStrideAll.getDevicePointer(),
                    childLTermOff.getDevicePointer(),
                    childLTermCnt.getDevicePointer(),
                    childTableBase.getDevicePointer(),
                    Pointer.to(childLowerAll),
                    Pointer.to(childUpperAll),
                    outLower.getDevicePointer(),
                    outUpper.getDevicePointer(),
                    Pointer.to(mStateCountArg),
                    Pointer.to(mStateOffsetArg),
                    Pointer.to(totalLambdaStatesArg),
                    Pointer.to(mCountArg),
                    Pointer.to(lambdaCountArg),
                    Pointer.to(lmPairCountArg),
                    Pointer.to(numChildrenArg),
                    Pointer.to(invRTArg)
                ));

                long startNanos = System.nanoTime();
                func.runAsync();
                outLower.downloadAsync();
                outUpper.downloadAsync();
                stream.waitForGpu();

                req.outTable.copyFrom(offset, outLower.getHostBuffer(),
                        outUpper.getHostBuffer(), chunk);

                if (req.progress) {
                    double ms = (System.nanoTime() - startNanos)/1e6;
                    System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP done, mStates=" + req.mStateCount
                            + " (chunk=" + chunk + "@" + offset + ")"
                            + ", lambdaStates=" + req.totalLambdaStates
                            + ", children=" + req.numChildren
                            + ", work=" + req.work
                            + ", blockThreads=" + blockThreads
                            + ", persistentContext=" + req.persistentContext
                            + ", elapsedMs=" + String.format(java.util.Locale.ROOT, "%.1f", ms));
                }

                offset += chunk;
                remaining -= chunk;
            }
        } finally {
            for (int i=buffers.size() - 1; i>=0; i--) {
                try {
                    buffers.get(i).cleanup();
                } catch (Throwable t) {
                    t.printStackTrace(System.err);
                }
            }
            for (int i=rawBufs.size() - 1; i>=0; i--) {
                try {
                    JCudaDriver.cuMemFree(rawBufs.get(i));
                } catch (Throwable t) {
                    t.printStackTrace(System.err);
                }
            }
        }
    }

    private static int resolveBlockThreads(Request req, GpuStream stream) {
        int requested = Math.max(1, req.blockThreads);
        int capped = Math.min(requested, stream.getContext().getGpu().getMaxBlockThreads());
        int powerOfTwo = Integer.highestOneBit(capped);
        return Math.max(1, powerOfTwo);
    }

    private static CUBuffer<IntBuffer> uploadInts(GpuStream stream, int[] values, List<CUBuffer<?>> buffers) {
        CUBuffer<IntBuffer> buf = stream.makeIntBuffer(Math.max(1, values.length));
        IntBuffer host = buf.getHostBuffer();
        host.clear();
        if (values.length > 0) {
            host.put(values);
        }
        while (host.position() < host.capacity()) {
            host.put(0);
        }
        buf.uploadAsync();
        buffers.add(buf);
        return buf;
    }

    private static CUBuffer<LongBuffer> uploadLongs(GpuStream stream, long[] values, List<CUBuffer<?>> buffers) {
        CUBuffer<LongBuffer> buf = stream.makeLongBuffer(Math.max(1, values.length));
        LongBuffer host = buf.getHostBuffer();
        host.clear();
        if (values.length > 0) {
            host.put(values);
        }
        while (host.position() < host.capacity()) {
            host.put(0L);
        }
        buf.uploadAsync();
        buffers.add(buf);
        return buf;
    }

    private static CUBuffer<DoubleBuffer> makeDoubles(GpuStream stream, int size, List<CUBuffer<?>> buffers) {
        CUBuffer<DoubleBuffer> buf = stream.makeDoubleBuffer(Math.max(1, size));
        buffers.add(buf);
        return buf;
    }

    /**
     * Upload a double[] straight into raw device memory, bypassing the 2 GiB cap of a
     * single java.nio direct buffer. cuMemAlloc takes a long byte count and a single
     * cuMemcpyHtoD from the Java array copies the whole (possibly > 2 GiB) buffer.
     * The caller must cuMemFree every returned pointer (tracked in rawBufs).
     */
    private static CUdeviceptr uploadDoublesBig(double[] values, List<CUdeviceptr> rawBufs) {
        long count = Math.max(1L, values.length);
        long bytes = count*(long)Double.BYTES;
        CUdeviceptr dptr = new CUdeviceptr();
        cudaCheck(JCudaDriver.cuMemAlloc(dptr, bytes), "cuMemAlloc(" + bytes + " bytes)");
        rawBufs.add(dptr);
        if (values.length > 0) {
            cudaCheck(JCudaDriver.cuMemcpyHtoD(dptr, Pointer.to(values),
                    (long)values.length*(long)Double.BYTES), "cuMemcpyHtoD");
        } else {
            cudaCheck(JCudaDriver.cuMemsetD8(dptr, (byte)0, bytes), "cuMemsetD8");
        }
        return dptr;
    }

    /**
     * Allocate one device buffer for all concatenated child tables and copy each
     * child's dense table straight into its slot at {@code base[c]} (element offset).
     * The child tables are NEVER concatenated on the Java heap, so their combined
     * length may exceed Integer.MAX_VALUE; only device memory (long byte count via
     * cuMemAlloc) bounds it. Device layout is byte-identical to the old single
     * concatenated upload, so the kernel and childTableBase semantics are unchanged.
     */
    private static CUdeviceptr uploadChildTables(double[][] parts, long[] base, long total,
            List<CUdeviceptr> rawBufs) {
        long count = Math.max(1L, total);
        long bytes = count*(long)Double.BYTES;
        CUdeviceptr dptr = new CUdeviceptr();
        cudaCheck(JCudaDriver.cuMemAlloc(dptr, bytes), "cuMemAlloc(" + bytes + " bytes)");
        rawBufs.add(dptr);
        for (int c = 0; c < parts.length; c++) {
            double[] part = parts[c];
            if (part == null || part.length == 0) {
                continue;
            }
            CUdeviceptr slot = dptr.withByteOffset(base[c]*(long)Double.BYTES);
            cudaCheck(JCudaDriver.cuMemcpyHtoD(slot, Pointer.to(part),
                    (long)part.length*(long)Double.BYTES), "cuMemcpyHtoD(child " + c + ")");
        }
        return dptr;
    }

    private static void cudaCheck(int result, String what) {
        if (result != CUresult.CUDA_SUCCESS) {
            throw new RuntimeException("CUDA " + what + " failed: " + CUresult.stringFor(result));
        }
    }

    private static void edgeFallback(Throwable t) {
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP edge fallback ("
                + t.getClass().getSimpleName() + ": " + t.getMessage()
                + "); this edge on Java DP path");
        if (getConfigBoolean("branchdp.dp.gpu.trace", false)) {
            t.printStackTrace(System.err);
        }
    }

    static long estimateDeviceBytes(Request req) {
        long bytes = 0L;
        bytes = addBytes(bytes, req.mCounts.length, Integer.BYTES);
        bytes = addBytes(bytes, req.lambdaCounts.length, Integer.BYTES);
        bytes = addBytes(bytes, req.lambdaOnlyRigid.length, Double.BYTES);
        bytes = addBytes(bytes, req.lambdaOnlyMin.length, Double.BYTES);
        bytes = addBytes(bytes, req.lmRigid.length, Double.BYTES);
        bytes = addBytes(bytes, req.lmMin.length, Double.BYTES);
        bytes = addBytes(bytes, req.lmLamSlots.length, Integer.BYTES);
        bytes = addBytes(bytes, req.lmMSlots.length, Integer.BYTES);
        bytes = addBytes(bytes, req.lmMCounts.length, Integer.BYTES);
        bytes = addBytes(bytes, req.lmOffsets.length, Long.BYTES);
        bytes = addBytes(bytes, req.childMSrcAll.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childMStrideAll.length, Long.BYTES);
        bytes = addBytes(bytes, req.childMTermOff.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childMTermCnt.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childLSrcAll.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childLStrideAll.length, Long.BYTES);
        bytes = addBytes(bytes, req.childLTermOff.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childLTermCnt.length, Integer.BYTES);
        bytes = addBytes(bytes, req.childTableBase.length, Long.BYTES);
        bytes = addBytes(bytes, req.childTableTotal, Double.BYTES);
        bytes = addBytes(bytes, req.childTableTotal, Double.BYTES);
        bytes = addBytes(bytes, req.mStateChunk, Double.BYTES);
        bytes = addBytes(bytes, req.mStateChunk, Double.BYTES);
        return bytes;
    }

    private static long addBytes(long current, int count, int elemBytes) {
        return addBytes(current, (long)count, elemBytes);
    }

    private static long addBytes(long current, long count, int elemBytes) {
        long c = Math.max(1L, count);
        if (c > Long.MAX_VALUE/(long)elemBytes) {
            return Long.MAX_VALUE;
        }
        long add = c*(long)elemBytes;
        if (current > Long.MAX_VALUE - add) {
            return Long.MAX_VALUE;
        }
        return current + add;
    }

    private static void markUnavailable(String reason, Throwable t) {
        unavailable = true;
        if (!unavailableLogged) {
            unavailableLogged = true;
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP unavailable (" + reason
                    + "); falling back to Java DP path");
            if (t != null && getConfigBoolean("branchdp.dp.gpu.trace", false)) {
                t.printStackTrace(System.err);
            }
        }
    }
}
