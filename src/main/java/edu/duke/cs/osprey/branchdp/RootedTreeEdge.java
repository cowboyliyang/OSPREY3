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

import com.sun.jna.Native;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.LongStream;

/**
 * An edge in the rooted branch decomposition tree with incremental partition function DP.
 *
 * Key design: lambda-states are pre-sorted by energy per M-state. Enumeration is incremental:
 *   - logZ_lower[m] = logSumExp of enumerated rigid-energy terms (valid: omitting positive terms)
 *   - logZ_upper[m] = logSumExp of enumerated min-energy terms + tail bound for remaining
 *
 * "Expanding a node" = enumerating the next lambda-state for a given (edge, mIdx).
 * This is analogous to original MARK*'s A* node expansion with H-score bounds.
 *
 * For leaf edges (no F-set): lambda-state energies come purely from emat, pre-sorted once.
 * For non-leaf edges (has F-set): re-enumeration needed when children change (uses looser tail bound).
 */
public class RootedTreeEdge {

    // ========== Tree structure ==========

    private RootedTreeNode parent;
    private RootedTreeNode child;
    private boolean isRootEdge;
    private boolean isLambdaEdge;

    // ========== Set structure ==========

    private LinkedHashSet<Integer> M;
    private LinkedHashSet<Integer> lambda;
    private LinkedHashSet<Integer> L;

    // F-set: nearest lambda-edge descendants
    private LinkedHashSet<RootedTreeEdge> Fset;

    // Compact tree pointers
    private RootedTreeNode compactLeftChild;
    private RootedTreeNode compactRightChild;
    private RootedTreeEdge compactParent;

    // ========== M-state indexing ==========

    private int[] mPositionsSorted;
    private int[] lambdaPositionsSorted;
    private RCs rcs;
    private long mStateCount;
    private int mArraySize;

    // ========== Log-space partition function bounds per M-state ==========

    private DPTable dpTable;      // log(Z_lower/upper) indexed by M-state

    // Option C (non-leaf DP kernel): precomputed child-fold plans. Built once,
    // before the (possibly parallel) DP loop, then read-only across DP workers.
    // null => fold disabled for this edge (legacy getMstateForFullState path).
    private ChildFoldPlan[] childFoldPlans;
    private boolean childFoldHoistInvariant;

    // Option A: native kernel toggle (resolved once before the DP loop) + per-thread
    // scratch row buffers reused across edges/M-states to avoid per-call allocation.
    private boolean nativeKernelEnabled;
    private static final ThreadLocal<double[]> NATIVE_LOWER_BUF = ThreadLocal.withInitial(() -> new double[0]);
    private static final ThreadLocal<double[]> NATIVE_UPPER_BUF = ThreadLocal.withInitial(() -> new double[0]);
    private static final ThreadLocal<double[]> NATIVE_OUT_BUF = ThreadLocal.withInitial(() -> new double[2]);

    private static final double NEG_INF = Double.NEGATIVE_INFINITY;
    private static final String DP_PARALLEL_PROPERTY = "branchdp.dp.parallel";
    private static final String DP_PARALLEL_MIN_M_STATES_PROPERTY = "branchdp.dp.parallel.minMStates";
    private static final String DP_PARALLEL_THREADS_PROPERTY = "branchdp.dp.parallel.threads";
    private static final String DP_COMPUTE_SHARD_SIZE_PROPERTY = "branchdp.dp.computeShardSize";
    private static final String DP_PROGRESS_PROPERTY = "branchdp.dp.progress";
    private static final int DEFAULT_DP_COMPUTE_SHARD_SIZE = 65_536;
    private static final String DP_MAX_M_STATES_PROPERTY = "branchdp.dp.maxMStates";
    private static final String DP_MAX_MATERIALIZED_PAIRS_PROPERTY = "branchdp.dp.maxMaterializedPairs";
    private static final String DP_TABLE_MODE_PROPERTY = "branchdp.dp.tableMode";
    private static final String DP_SHARD_SIZE_PROPERTY = "branchdp.dp.shardSize";
    private static final int DEFAULT_DP_PARALLEL_MIN_M_STATES = 1024;
    private static final long DEFAULT_DP_MAX_MATERIALIZED_PAIRS = 50_000_000L;
    private static final int DEFAULT_DP_SHARD_SIZE = 1_048_576;
    // Option C: non-leaf child-table folding (pure Java, no SIMD/JIT-vector risk).
    private static final String DP_FOLD_CHILDREN_PROPERTY = "branchdp.dp.foldChildren";
    private static final String DP_FOLD_HOIST_PROPERTY = "branchdp.dp.foldChildren.hoistInvariant";
    // Option A: native SIMD log-sum-exp kernel. Default OFF (wired but dormant:
    // libOspreyLogSumExp.so is built/validated separately, see src/main/c/).
    private static final String DP_NATIVE_KERNEL_PROPERTY = "branchdp.dp.nativeKernel";
    // CUDA full-DP fast path. Default OFF and intentionally narrow: one non-leaf
    // child, dense parent/child DP tables, and enough m*lambda work to amortize
    // GPU packing/transfer overhead.
    private static final String DP_GPU_PROPERTY = "branchdp.dp.gpu";
    private static final String DP_GPU_MIN_WORK_PROPERTY = "branchdp.dp.gpu.minWork";
    private static final String DP_GPU_MAX_BYTES_PROPERTY = "branchdp.dp.gpu.maxBytes";
    private static final String DP_GPU_FAIL_IF_EXCEEDS_VRAM_PROPERTY = "branchdp.dp.gpu.failIfExceedsVram";
    private static final String DP_GPU_BLOCK_THREADS_PROPERTY = "branchdp.dp.gpu.blockThreads";
    private static final String DP_GPU_OUTPUT_TILE_MSTATES_PROPERTY = "branchdp.dp.gpu.outputTileMStates";
    private static final String DP_GPU_PERSISTENT_CONTEXT_PROPERTY = "branchdp.dp.gpu.persistentContext";
    private static final long DEFAULT_DP_GPU_MIN_WORK = 50_000_000L;
    // Single device buffers can now exceed 2 GiB (chunked raw upload in DPGpuFullDP),
    // so 0 = no static cap: the real budget is detected free VRAM (cuMemGetInfo, with
    // headroom) inside DPGpuFullDP. A positive value re-enables an explicit cap.
    private static final long DEFAULT_DP_GPU_MAX_BYTES = 0L;
    private static final int DEFAULT_DP_GPU_BLOCK_THREADS = 256;
    private static final long DEFAULT_DP_GPU_OUTPUT_TILE_MSTATES = 1_048_576L;
    private static final boolean DEFAULT_DP_GPU_PERSISTENT_CONTEXT = true;
    private static final int DP_GPU_MAX_CHILDREN = 64;
    private static final String DP_GPU_MULTI_PROPERTY = "branchdp.dp.gpu.multiGpu";
    private static final String DP_GPU_MAX_GPUS_PROPERTY = "branchdp.dp.gpu.maxGpus";
    private static final String DP_GPU_MIN_MSTATES_PER_GPU_PROPERTY = "branchdp.dp.gpu.minMStatesPerGpu";
    private static final long DEFAULT_DP_GPU_MIN_MSTATES_PER_GPU = 4096L;
    private static final String PAC_SAMPLING_GPU_MULTI_PROPERTY = "branchdp.pac.sampling.gpu.multiGpu";
    private static final String PAC_SAMPLING_GPU_MAX_GPUS_PROPERTY = "branchdp.pac.sampling.gpu.maxGpus";
    private static final String PAC_SAMPLING_GPU_MIN_GROUPS_PER_GPU_PROPERTY = "branchdp.pac.sampling.gpu.minGroupsPerGpu";
    private static final String PAC_SAMPLING_GPU_PERSISTENT_CONTEXT_PROPERTY = "branchdp.pac.sampling.gpu.persistentContext";
    private static final String PAC_SAMPLING_GPU_RESIDENT_CHILD_TABLES_PROPERTY = "branchdp.pac.sampling.gpu.residentChildTables";
    private static final String PAC_SAMPLING_GPU_METHOD_PROPERTY = "branchdp.pac.sampling.gpu.method";
    private static final int DEFAULT_PAC_SAMPLING_GPU_MIN_GROUPS_PER_GPU = 1;
    private static final boolean DEFAULT_PAC_SAMPLING_GPU_PERSISTENT_CONTEXT = true;
    private static final boolean DEFAULT_PAC_SAMPLING_GPU_RESIDENT_CHILD_TABLES = true;
    private static final SamplingGpuPhase1.Method DEFAULT_PAC_SAMPLING_GPU_METHOD =
            SamplingGpuPhase1.Method.GUMBEL;

    // ========== Incremental lambda enumeration state ==========

    /**
     * Per-M-state sorted lambda-states and enumeration progress.
     * sortedLambdaIndices[mIdx] = array of lambda-state indices sorted by E_min ascending.
     * enumeratedCount[mIdx] = how many lambda-states have been enumerated so far (k out of K).
     */
    private int[][] sortedLambdaIndices;  // [mIdx][rank] -> lambdaIdx
    private int[] enumeratedCount;         // [mIdx] -> k
    private int totalLambdaStates;         // K = total number of lambda-states

    // Pre-computed per-lambda-state energies (M-independent part: lambda-only)
    // Full energy per (m, λ) = lambdaOnlyRigid[λ] + lambdaMRigid(m,λ) for rigid
    private double[] lambdaOnlyRigid;      // [lambdaIdx]
    private double[] lambdaOnlyMin;        // [lambdaIdx]

    // Cached energy matrices and parameters for on-demand computation
    private EnergyMatrix cachedRigidEmat;
    private EnergyMatrix cachedMinEmat;
    private InteractionGraph cachedG;
    private double cachedRT;

    // Per (mIdx, lambdaIdx) full energies, lazily computed and cached
    // fullEnergyRigid[mIdx][lambdaIdx], fullEnergyMin[mIdx][lambdaIdx]
    private double[][] fullEnergyRigid;
    private double[][] fullEnergyMin;
    private boolean fullEnergyTablesMaterialized;

    // ========== Constructor ==========

    public RootedTreeEdge(RootedTreeNode parent, RootedTreeNode child,
                           LinkedHashSet<Integer> mSet, boolean isRootEdge, RCs rcs) {
        this.parent = parent;
        this.child = child;
        this.M = new LinkedHashSet<>(mSet);
        this.isRootEdge = isRootEdge;
        this.rcs = rcs;
        this.isLambdaEdge = false;
    }

    // ========== Getters ==========

    public RootedTreeNode getParent() { return parent; }
    public RootedTreeNode getChild() { return child; }
    public boolean getIsRootEdge() { return isRootEdge; }
    public boolean getIsLambdaEdge() { return isLambdaEdge; }
    public LinkedHashSet<Integer> getM() { return M; }
    public LinkedHashSet<Integer> getLambda() { return lambda; }
    public LinkedHashSet<Integer> getL() { return L; }
    public double[] getLogZLower() { return dpTable == null ? null : dpTable.lowerArrayUnsafe(); }
    public double[] getLogZUpper() { return dpTable == null ? null : dpTable.upperArrayUnsafe(); }
    public double getLogZLower(long mIdx) { return dpTable.lower(mIdx); }
    public double getLogZUpper(long mIdx) { return dpTable.upper(mIdx); }
    public void setLogZ(long mIdx, double lower, double upper) { dpTable.set(mIdx, lower, upper); }
    public boolean hasDPTable() { return dpTable != null; }
    public boolean hasDenseDPTable() { return dpTable != null && dpTable.isDenseArrayBacked(); }
    public long getDPTableBytes() { return dpTable == null ? estimateDPTableBytes(mStateCount) : dpTable.estimatedBytes(); }
    public RootedTreeNode getCompactLeftChild() { return compactLeftChild; }
    public RootedTreeNode getCompactRightChild() { return compactRightChild; }
    public void setCompactParent(RootedTreeEdge p) { this.compactParent = p; }
    public RootedTreeEdge getCompactParent() { return compactParent; }
    public int getMArraySize() { return requireDenseMStateCount(); }
    public long getMStateCount() { return mStateCount; }
    public int[] getMPositionsSorted() { return mPositionsSorted; }
    public int[] getLambdaPositionsSorted() { return lambdaPositionsSorted; }
    public LinkedHashSet<RootedTreeEdge> getFset() { return Fset; }
    public int getTotalLambdaStates() { return totalLambdaStates; }
    public RCs getRCs() { return rcs; }
    public int getEnumeratedCount(int mIdx) { return enumeratedCount != null ? enumeratedCount[mIdx] : 0; }
    public double[][] getFullEnergyMin() { return fullEnergyMin; }
    public double[][] getFullEnergyRigid() { return fullEnergyRigid; }

    // ========== Log-space utility ==========

    public static double logSumExp(double a, double b) {
        if (a == NEG_INF) return b;
        if (b == NEG_INF) return a;
        double max = Math.max(a, b);
        return max + Math.log1p(Math.exp(-Math.abs(a - b)));
    }

    private static String getConfigProperty(String key, String defaultValue) {
        return BranchDpConfig.getBackendProperty(key, defaultValue);
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        return BranchDpConfig.getBackendBoolean(key, defaultValue);
    }

    private static int getConfigInteger(String key, int defaultValue) {
        return BranchDpConfig.getBackendInteger(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    private static long getConfigLong(String key, long defaultValue) {
        return BranchDpConfig.getBackendLong(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    // ========== compLlambda: compute L and lambda sets ==========

    public void compLlambda() {
        compLlambda(true);
    }

    public void compLlambda(boolean initEnumerationArrays) {
        RootedTreeNode clc = child.getLeftChild();

        if (clc == null) {
            L = new LinkedHashSet<>();
            if (child.getV1() >= 0) L.add(child.getV1());
            if (child.getV2() >= 0) L.add(child.getV2());
            L.removeAll(M);
            lambda = new LinkedHashSet<>(L);
        } else {
            RootedTreeEdge leftEdge = clc.getChildOfEdge();
            RootedTreeNode crc = child.getRightChild();
            RootedTreeEdge rightEdge = crc == null ? null : crc.getChildOfEdge();

            LinkedHashSet<Integer> uMc = new LinkedHashSet<>(leftEdge.getM());
            if (rightEdge != null) {
                uMc.addAll(rightEdge.getM());
            }
            lambda = new LinkedHashSet<>(uMc);
            lambda.removeAll(M);

            L = new LinkedHashSet<>(lambda);
            L.addAll(leftEdge.getL());
            if (rightEdge != null) {
                L.addAll(rightEdge.getL());
            }
        }

        if (!lambda.isEmpty()) {
            isLambdaEdge = true;
            Fset = new LinkedHashSet<>();
            computeFset(child);
            initializeIndexing();
            if (initEnumerationArrays) {
                initializeArrays();
            }
        } else {
            isLambdaEdge = false;
        }
    }

    // ========== F-set computation ==========

    private void computeFset(RootedTreeNode tn) {
        if (tn.getLeftChild() != null) {
            RootedTreeEdge leftEdge = tn.getLeftChild().getChildOfEdge();
            if (leftEdge.isLambdaEdge) {
                Fset.add(leftEdge);
            } else {
                computeFset(tn.getLeftChild());
            }
        }
        if (tn.getRightChild() != null) {
            RootedTreeEdge rightEdge = tn.getRightChild().getChildOfEdge();
            if (rightEdge.isLambdaEdge) {
                Fset.add(rightEdge);
            } else {
                computeFset(tn.getRightChild());
            }
        }
    }

    // ========== Array initialization ==========

    private void initializeIndexing() {
        mPositionsSorted = M.stream().mapToInt(Integer::intValue).sorted().toArray();
        lambdaPositionsSorted = lambda.stream().mapToInt(Integer::intValue).sorted().toArray();

        mStateCount = checkedStateCountLong("M", mPositionsSorted, true);
        mArraySize = mStateCount <= Integer.MAX_VALUE ? (int) mStateCount : -1;
        totalLambdaStates = checkedStateCount("lambda", lambdaPositionsSorted, false);
    }

    private int checkedStateCount(String label, int[] positions, boolean isMStateCount) {
        long count = checkedStateCountLong(label, positions, isMStateCount);
        if (count > Integer.MAX_VALUE) {
            throw stateCountException(label, positions, count);
        }
        return (int) count;
    }

    private long checkedStateCountLong(String label, int[] positions, boolean isMStateCount) {
        long count = 1L;
        for (int pos : positions) {
            int n = rcs.getNum(pos);
            if (n < 0 || count > Long.MAX_VALUE / Math.max(1, n)) {
                throw stateCountException(label, positions, Long.MAX_VALUE);
            }
            count *= n;
        }

        long maxStates = isMStateCount
                ? getConfigLong(DP_MAX_M_STATES_PROPERTY, Long.MAX_VALUE)
                : Integer.MAX_VALUE;
        if (count > maxStates) {
            throw new DPTableTooLargeException(label, count, positions,
                    "DP " + label + "-state count exceeds "
                            + DP_MAX_M_STATES_PROPERTY + "=" + maxStates
                            + ". Reduce branchwidth/rooting, enable a shard-backed DP table, "
                            + "or raise the limit only if memory allows.");
        }
        return count;
    }

    private DPTableTooLargeException stateCountException(String label, int[] positions, long count) {
        return new DPTableTooLargeException(label, count, positions,
                "DP " + label + "-state count exceeds Java array indexing limits. "
                        + "This used to overflow into NegativeArraySizeException; "
                        + "use a smaller separator or a shard-backed DP table.");
    }

    private void initializeArrays() {
        dpTable = makeDPTable();
        dpTable.fill(NEG_INF, NEG_INF);

        enumeratedCount = mStateCount <= Integer.MAX_VALUE
                ? new int[(int) mStateCount]
                : null;
        // sortedLambdaIndices allocated lazily per M-state in initIncrementalEnumeration
    }

    private DPTable makeDPTable() {
        String mode = getConfigProperty(DP_TABLE_MODE_PROPERTY, "auto").trim().toLowerCase(Locale.ROOT);
        int shardSize = Math.max(1, getConfigInteger(DP_SHARD_SIZE_PROPERTY, DEFAULT_DP_SHARD_SIZE));
        switch (mode) {
            case "":
            case "auto":
                return mStateCount <= Integer.MAX_VALUE
                        ? new DenseDPTable(mStateCount)
                        : new ShardedDPTable(mStateCount, shardSize);
            case "dense":
                if (mStateCount > Integer.MAX_VALUE) {
                    throw stateCountException("M", mPositionsSorted, mStateCount);
                }
                return new DenseDPTable(mStateCount);
            case "sharded":
            case "shard":
                return new ShardedDPTable(mStateCount, shardSize);
            default:
                System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown DP table mode '" + mode
                        + "', using auto.");
                return mStateCount <= Integer.MAX_VALUE
                        ? new DenseDPTable(mStateCount)
                        : new ShardedDPTable(mStateCount, shardSize);
        }
    }

    private int requireDenseMStateCount() {
        if (mArraySize < 0) {
            throw new IllegalStateException("DP M-state count " + mStateCount
                    + " exceeds dense int indexing; use getMStateCount() and DPTable access.");
        }
        return mArraySize;
    }

    private static long estimateDPTableBytes(long mStates) {
        if (mStates > Long.MAX_VALUE / (2L * (long) Double.BYTES)) {
            return Long.MAX_VALUE;
        }
        return 2L * mStates * (long) Double.BYTES;
    }

    // ========== M-state indexing ==========

    public long computeIndexInA(int[] mRCIndices) {
        if (isRootEdge || mRCIndices == null || mPositionsSorted.length == 0) {
            return 0L;
        }

        long index = mRCIndices[mPositionsSorted.length - 1];
        long stride = 1L;
        for (int i = mPositionsSorted.length - 2; i >= 0; i--) {
            stride *= rcs.getNum(mPositionsSorted[i + 1]);
            index += (long) mRCIndices[i] * stride;
        }
        return index;
    }

    public int[] getMstateForFullState(int[] mRCs, int[] lambdaRCs, RootedTreeEdge childEdge) {
        int[] childMPositions = childEdge.mPositionsSorted;
        int[] result = new int[childMPositions.length];

        for (int i = 0; i < childMPositions.length; i++) {
            int targetPos = childMPositions[i];
            result[i] = -1;

            for (int j = 0; j < mPositionsSorted.length; j++) {
                if (mPositionsSorted[j] == targetPos) {
                    result[i] = mRCs[j];
                    break;
                }
            }
            if (result[i] >= 0) continue;

            for (int j = 0; j < lambdaPositionsSorted.length; j++) {
                if (lambdaPositionsSorted[j] == targetPos) {
                    result[i] = lambdaRCs[j];
                    break;
                }
            }
        }
        return result;
    }

    public int[] decodeMStatePublic(long mIndex) {
        return decodeMState(mIndex);
    }

    public int[] decodeLambdaStatePublic(int lambdaIndex) {
        return decodeLambdaState(lambdaIndex);
    }

    private int[] decodeMState(long mIndex) {
        int[] mRCs = new int[mPositionsSorted.length];
        long remaining = mIndex;
        for (int i = mPositionsSorted.length - 1; i >= 0; i--) {
            int numRCs = rcs.getNum(mPositionsSorted[i]);
            mRCs[i] = (int) (remaining % numRCs);
            remaining /= numRCs;
        }
        return mRCs;
    }

    private int[] decodeLambdaState(int lambdaIndex) {
        int[] lambdaRCs = new int[lambdaPositionsSorted.length];
        int remaining = lambdaIndex;
        for (int i = lambdaPositionsSorted.length - 1; i >= 0; i--) {
            int numRCs = rcs.getNum(lambdaPositionsSorted[i]);
            lambdaRCs[i] = remaining % numRCs;
            remaining /= numRCs;
        }
        return lambdaRCs;
    }

    // ========== Compact tree ==========

    public void compactTree() {
        compactLeftChild = searchSubtree(child.getLeftChild());
        compactRightChild = searchSubtree(child.getRightChild());

        if (compactLeftChild != null) {
            compactLeftChild.getChildOfEdge().compactTree();
            compactLeftChild.getChildOfEdge().setCompactParent(this);
        }
        if (compactRightChild != null) {
            compactRightChild.getChildOfEdge().compactTree();
            compactRightChild.getChildOfEdge().setCompactParent(this);
        }

        if (compactLeftChild == null && compactRightChild != null) {
            compactLeftChild = compactRightChild;
            compactRightChild = null;
        }
    }

    private RootedTreeNode searchSubtree(RootedTreeNode tn) {
        if (tn == null) return null;
        if (tn.getChildOfEdge() != null && tn.getChildOfEdge().isLambdaEdge) {
            return tn;
        }
        if (tn.getIsLeaf()) return null;

        RootedTreeNode left = searchSubtree(tn.getLeftChild());
        RootedTreeNode right = searchSubtree(tn.getRightChild());

        if (left != null && right != null) return tn;
        if (left != null) return left;
        return right;
    }

    // ========== TESS computation ==========

    public long computeTESS() {
        long leftTESS = 0;
        if (compactLeftChild != null) {
            leftTESS = compactLeftChild.getChildOfEdge().computeTESS();
        }
        long rightTESS = 0;
        if (compactRightChild != null) {
            rightTESS = compactRightChild.getChildOfEdge().computeTESS();
        }

        long selfTESS = 1;
        if (mPositionsSorted != null) {
            for (int pos : mPositionsSorted) selfTESS *= rcs.getNum(pos);
        }
        if (lambdaPositionsSorted != null) {
            for (int pos : lambdaPositionsSorted) selfTESS *= rcs.getNum(pos);
        }

        return leftTESS + rightTESS + selfTESS;
    }

    public double computeLogTESS() {
        List<Double> logTerms = new ArrayList<>();
        collectLogTESSTerms(logTerms);

        if (logTerms.isEmpty()) return 0.0;

        double maxLog = Collections.max(logTerms);
        double sumExp = 0.0;
        for (double logTerm : logTerms) {
            sumExp += Math.exp(logTerm - maxLog);
        }
        return maxLog + Math.log(sumExp);
    }

    private void collectLogTESSTerms(List<Double> logTerms) {
        if (compactLeftChild != null) {
            compactLeftChild.getChildOfEdge().collectLogTESSTerms(logTerms);
        }
        if (compactRightChild != null) {
            compactRightChild.getChildOfEdge().collectLogTESSTerms(logTerms);
        }

        double logSelf = 0.0;
        if (mPositionsSorted != null) {
            for (int pos : mPositionsSorted) logSelf += Math.log(rcs.getNum(pos));
        }
        if (lambdaPositionsSorted != null) {
            for (int pos : lambdaPositionsSorted) logSelf += Math.log(rcs.getNum(pos));
        }
        logTerms.add(logSelf);
    }

    // ========== Energy computation helpers ==========

    private double computeLambdaOnlyEnergy(int[] lambdaRCs, EnergyMatrix emat, InteractionGraph G) {
        double energy = 0.0;
        for (int i = 0; i < lambdaPositionsSorted.length; i++) {
            int pos = lambdaPositionsSorted[i];
            int rc = rcs.get(pos, lambdaRCs[i]);
            energy += emat.getOneBody(pos, rc);
        }
        for (int i = 0; i < lambdaPositionsSorted.length; i++) {
            int posI = lambdaPositionsSorted[i];
            int rcI = rcs.get(posI, lambdaRCs[i]);
            for (int j = i + 1; j < lambdaPositionsSorted.length; j++) {
                int posJ = lambdaPositionsSorted[j];
                if (G.hasEdge(posI, posJ)) {
                    int rcJ = rcs.get(posJ, lambdaRCs[j]);
                    energy += emat.getPairwise(posI, rcI, posJ, rcJ);
                }
            }
        }
        return energy;
    }

    private double computeLambdaMEnergy(int[] mRCs, int[] lambdaRCs,
                                         EnergyMatrix emat, InteractionGraph G) {
        double energy = 0.0;
        for (int i = 0; i < lambdaPositionsSorted.length; i++) {
            int posI = lambdaPositionsSorted[i];
            int rcI = rcs.get(posI, lambdaRCs[i]);
            for (int j = 0; j < mPositionsSorted.length; j++) {
                int posJ = mPositionsSorted[j];
                if (G.hasEdge(posI, posJ)) {
                    int rcJ = rcs.get(posJ, mRCs[j]);
                    energy += emat.getPairwise(posI, rcI, posJ, rcJ);
                }
            }
        }
        return energy;
    }

    private double computeLocalEnergy(int[] mRCs, int[] lambdaRCs,
                                       EnergyMatrix emat, InteractionGraph G) {
        return computeLambdaOnlyEnergy(lambdaRCs, emat, G)
                + computeLambdaMEnergy(mRCs, lambdaRCs, emat, G);
    }

    // ========== Incremental lambda enumeration ==========

    /**
     * Initialize incremental enumeration for all edges in the tree (post-order).
     * Pre-computes lambda-only energies and per-M-state sorted lambda indices.
     * Initial logZ bounds are loose (k=0 for all M-states):
     *   logZ_lower = -inf (nothing enumerated)
     *   logZ_upper = tail bound (all K lambda-states in tail)
     */
    public void initIncrementalEnumeration(EnergyMatrix rigidEmat, EnergyMatrix minEmat,
                                            InteractionGraph G, double RT) {
        initIncrementalEnumeration(rigidEmat, minEmat, G, RT, false);
    }

    public void initIncrementalEnumeration(EnergyMatrix rigidEmat, EnergyMatrix minEmat,
                                            InteractionGraph G, double RT,
                                            boolean materializeFullEnergyTables) {
        if (!isLambdaEdge) return;

        this.cachedRigidEmat = rigidEmat;
        this.cachedMinEmat = minEmat;
        this.cachedG = G;
        this.cachedRT = RT;
        this.fullEnergyTablesMaterialized = materializeFullEnergyTables;

        // Pre-compute lambda-only energies (M-independent)
        lambdaOnlyRigid = new double[totalLambdaStates];
        lambdaOnlyMin = new double[totalLambdaStates];
        for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
            int[] lambdaRCs = decodeLambdaState(lIdx);
            lambdaOnlyRigid[lIdx] = computeLambdaOnlyEnergy(lambdaRCs, rigidEmat, G);
            lambdaOnlyMin[lIdx] = computeLambdaOnlyEnergy(lambdaRCs, minEmat, G);
        }

        fullEnergyRigid = null;
        fullEnergyMin = null;
        sortedLambdaIndices = null;

        if (materializeFullEnergyTables) {
            requireDenseMStateCount();
            checkMaterializedPairCount();

            // Allocate per-(mIdx, lambdaIdx) energy cache only for legacy paths
            // that mutate or sort per-entry energies (region-atom DP integration,
            // incremental enumeration). Full DP uses streaming and does not need it.
            fullEnergyRigid = new double[mArraySize][totalLambdaStates];
            fullEnergyMin = new double[mArraySize][totalLambdaStates];
            for (double[] row : fullEnergyRigid) Arrays.fill(row, Double.NaN);
            for (double[] row : fullEnergyMin) Arrays.fill(row, Double.NaN);

            // For leaf edges (no F-set): pre-compute full energies and sort per M-state.
            // For non-leaf edges: sorting depends on child logZ, handled differently.
            sortedLambdaIndices = new int[mArraySize][];

            if (!hasFsetChildren()) {
                // Leaf edge: full energy = lambdaOnly + lambdaM (no child logZ)
                for (int mIdx = 0; mIdx < mArraySize; mIdx++) {
                    int[] mRCs = decodeMState(mIdx);
                    double[] minEnergies = new double[totalLambdaStates];
                    for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                        int[] lambdaRCs = decodeLambdaState(lIdx);
                        double lmRigid = computeLambdaMEnergy(mRCs, lambdaRCs, rigidEmat, G);
                        double lmMin = computeLambdaMEnergy(mRCs, lambdaRCs, minEmat, G);
                        fullEnergyRigid[mIdx][lIdx] = lambdaOnlyRigid[lIdx] + lmRigid;
                        fullEnergyMin[mIdx][lIdx] = lambdaOnlyMin[lIdx] + lmMin;
                        minEnergies[lIdx] = fullEnergyMin[mIdx][lIdx];
                    }
                    // Sort lambda indices by E_min ascending (highest Boltzmann first)
                    sortedLambdaIndices[mIdx] = sortIndicesBy(minEnergies);
                }
            } else {
                // Non-leaf edge: use identity order initially (will recompute on demand)
                for (int mIdx = 0; mIdx < mArraySize; mIdx++) {
                    sortedLambdaIndices[mIdx] = new int[totalLambdaStates];
                    for (int i = 0; i < totalLambdaStates; i++) {
                        sortedLambdaIndices[mIdx][i] = i;
                    }
                }
            }
        }

        // Set initial logZ bounds: k=0 for all M-states
        dpTable.fill(NEG_INF, NEG_INF);
        if (materializeFullEnergyTables) {
            // logZ_upper = tail bound for all K lambda-states in incremental mode.
            for (int mIdx = 0; mIdx < mArraySize; mIdx++) {
                dpTable.set(mIdx, NEG_INF, computeTailBound(mIdx, 0));
            }
        }
    }

    private void checkMaterializedPairCount() {
        long pairs = mStateCount * (long) totalLambdaStates;
        long maxPairs = getConfigLong(DP_MAX_MATERIALIZED_PAIRS_PROPERTY,
                DEFAULT_DP_MAX_MATERIALIZED_PAIRS);
        if (pairs > maxPairs) {
            throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix() + " Refusing to materialize DP full-energy table with "
                    + pairs + " (M,lambda) entries for M=" + Arrays.toString(mPositionsSorted)
                    + ", lambda=" + Arrays.toString(lambdaPositionsSorted)
                    + ". Full DP now streams this table; region-atom/incremental paths can raise "
                    + DP_MAX_MATERIALIZED_PAIRS_PROPERTY + " only if memory allows.");
        }
    }

    /**
     * Sort indices [0..n-1] by values[i] ascending. Returns sorted index array.
     */
    private static int[] sortIndicesBy(double[] values) {
        Integer[] indices = new Integer[values.length];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, Comparator.comparingDouble(i -> values[i]));
        int[] result = new int[indices.length];
        for (int i = 0; i < indices.length; i++) result[i] = indices[i];
        return result;
    }

    /**
     * Compute tail bound for un-enumerated lambda-states.
     *
     * For leaf edges: remaining (K-k) states each have E_min >= E_min of the k-th sorted state.
     *   tail = log(K-k) + (-E_min(sorted[k]) / RT)
     *
     * For non-leaf edges: uses a looser bound based on best possible local energy
     *   + best child logZ_upper values.
     *
     * @param mIdx  M-state index
     * @param k     number already enumerated
     * @return log-space tail bound (upper)
     */
    private double computeTailBound(int mIdx, int k) {
        int remaining = totalLambdaStates - k;
        if (remaining <= 0) return NEG_INF;

        if (!hasFsetChildren()) {
            // Leaf edge: sorted by E_min, so sorted[k] is the best remaining
            int bestRemainingLambdaIdx = sortedLambdaIndices[mIdx][k];
            double bestRemainingEMin = fullEnergyMin[mIdx][bestRemainingLambdaIdx];
            return Math.log(remaining) + (-bestRemainingEMin / cachedRT);
        } else {
            // Non-leaf edge: use conservative bound
            // Best local energy among all lambda-states
            double bestLocalEMin = Double.MAX_VALUE;
            for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                int[] mRCs = decodeMState(mIdx);
                int[] lambdaRCs = decodeLambdaState(lIdx);
                double localE = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);
                bestLocalEMin = Math.min(bestLocalEMin, localE);
            }
            // Best child logZ_upper for any M-state mapping
            double bestFSumUpper = 0.0;
            if (Fset != null) {
                for (RootedTreeEdge fEdge : Fset) {
                    double bestChildUpper = NEG_INF;
                    for (long fi = 0; fi < fEdge.getMStateCount(); fi++) {
                        bestChildUpper = Math.max(bestChildUpper, fEdge.getLogZUpper(fi));
                    }
                    bestFSumUpper += bestChildUpper;
                }
            }
            return Math.log(remaining) + (-bestLocalEMin / cachedRT) + bestFSumUpper;
        }
    }

    /**
     * Enumerate the next lambda-state for a given M-state.
     * Returns true if a state was enumerated, false if fully enumerated already.
     *
     * After calling this, logZ_lower and logZ_upper for this mIdx are updated.
     */
    public boolean enumerateNextLambda(int mIdx) {
        if (sortedLambdaIndices == null) {
            throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix() + " Incremental lambda enumeration requires materialized DP energy tables. "
                    + "Initialize with materializeFullEnergyTables=true or use computeFullDP().");
        }
        int k = enumeratedCount[mIdx];
        if (k >= totalLambdaStates) return false;

        int lambdaIdx = sortedLambdaIndices[mIdx][k];

        if (!hasFsetChildren()) {
            // Leaf edge: use pre-computed energies
            double eRigid = fullEnergyRigid[mIdx][lambdaIdx];
            double eMin = fullEnergyMin[mIdx][lambdaIdx];

            double logTermLower = -eRigid / cachedRT;
            double logTermUpper = -eMin / cachedRT;

            dpTable.set(mIdx, logSumExp(dpTable.lower(mIdx), logTermLower), dpTable.upper(mIdx));

            // Upper = enumerated part + new tail
            // Recompute upper from scratch: enumerated terms + tail
            // (simpler and avoids accumulation errors)
            double enumeratedUpper = NEG_INF;
            for (int i = 0; i <= k; i++) {
                int li = sortedLambdaIndices[mIdx][i];
                enumeratedUpper = logSumExp(enumeratedUpper, -fullEnergyMin[mIdx][li] / cachedRT);
            }
            double tailUpper = computeTailBound(mIdx, k + 1);
            dpTable.set(mIdx, dpTable.lower(mIdx), logSumExp(enumeratedUpper, tailUpper));

        } else {
            // Non-leaf edge: need child logZ values
            int[] mRCs = decodeMState(mIdx);
            int[] lambdaRCs = decodeLambdaState(lambdaIdx);

            double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
            double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

            double fSumLower = 0.0;
            double fSumUpper = 0.0;
            if (Fset != null) {
                for (RootedTreeEdge fEdge : Fset) {
                    int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                    long fIndex = fEdge.computeIndexInA(fM);
                    fSumLower += fEdge.getLogZLower(fIndex);
                    fSumUpper += fEdge.getLogZUpper(fIndex);
                }
            }

            double logTermLower = -eRigid / cachedRT + fSumLower;
            double logTermUpper = -eMin / cachedRT + fSumUpper;

            dpTable.set(mIdx, logSumExp(dpTable.lower(mIdx), logTermLower), dpTable.upper(mIdx));

            // Recompute upper: all enumerated terms + tail
            double enumeratedUpper = NEG_INF;
            for (int i = 0; i <= k; i++) {
                int li = sortedLambdaIndices[mIdx][i];
                int[] liRCs = decodeLambdaState(li);
                double liEMin = computeLocalEnergy(mRCs, liRCs, cachedMinEmat, cachedG);
                double liFSumUpper = 0.0;
                if (Fset != null) {
                    for (RootedTreeEdge fEdge : Fset) {
                        int[] fM = getMstateForFullState(mRCs, liRCs, fEdge);
                        long fIndex = fEdge.computeIndexInA(fM);
                        liFSumUpper += fEdge.getLogZUpper(fIndex);
                    }
                }
                enumeratedUpper = logSumExp(enumeratedUpper, -liEMin / cachedRT + liFSumUpper);
            }
            double tailUpper = computeTailBound(mIdx, k + 1);
            dpTable.set(mIdx, dpTable.lower(mIdx), logSumExp(enumeratedUpper, tailUpper));
        }

        enumeratedCount[mIdx] = k + 1;
        return true;
    }

    /**
     * Check if this M-state can be further tightened by lambda enumeration.
     */
    public boolean canEnumerate(int mIdx) {
        return enumeratedCount != null && enumeratedCount[mIdx] < totalLambdaStates;
    }

    /**
     * Check if all M-states are fully enumerated.
     */
    public boolean isFullyEnumerated() {
        if (enumeratedCount == null) return true;
        for (int k : enumeratedCount) {
            if (k < totalLambdaStates) return false;
        }
        return true;
    }

    /** Whether this edge has F-set children (non-leaf in compact tree). */
    public boolean hasFsetChildren() {
        return Fset != null && !Fset.isEmpty();
    }

    // ========== Error and bounds ==========

    /**
     * Get the error (logUpper - logLower) for a specific M-state.
     */
    public double getErrorForMState(int mIdx) {
        if (dpTable.upper(mIdx) == NEG_INF && dpTable.lower(mIdx) == NEG_INF) return 0.0;
        if (dpTable.lower(mIdx) == NEG_INF) return Double.MAX_VALUE;
        return dpTable.upper(mIdx) - dpTable.lower(mIdx);
    }

    // ========== Recompute non-leaf edge bounds from children ==========

    /**
     * Recompute logZ bounds for a non-leaf edge after children have been tightened.
     * Re-enumerates all previously enumerated lambda-states with current child bounds,
     * plus recomputes tail bound.
     */
    public void recomputeFromChildren() {
        if (!hasFsetChildren()) return;

        for (int mIdx = 0; mIdx < mArraySize; mIdx++) {
            int k = enumeratedCount[mIdx];
            int[] mRCs = decodeMState(mIdx);

            double newLower = NEG_INF;
            double newUpper = NEG_INF;

            // Re-evaluate all enumerated lambda-states with current child bounds
            for (int i = 0; i < k; i++) {
                int lambdaIdx = sortedLambdaIndices[mIdx][i];
                int[] lambdaRCs = decodeLambdaState(lambdaIdx);

                double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

                double fSumLower = 0.0;
                double fSumUpper = 0.0;
                if (Fset != null) {
                    for (RootedTreeEdge fEdge : Fset) {
                        int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                        long fIndex = fEdge.computeIndexInA(fM);
                        fSumLower += fEdge.getLogZLower(fIndex);
                        fSumUpper += fEdge.getLogZUpper(fIndex);
                    }
                }

                newLower = logSumExp(newLower, -eRigid / cachedRT + fSumLower);
                newUpper = logSumExp(newUpper, -eMin / cachedRT + fSumUpper);
            }

            // Add tail bound
            double tail = computeTailBound(mIdx, k);
            newUpper = logSumExp(newUpper, tail);

            dpTable.set(mIdx, newLower, newUpper);
        }
    }

    /**
     * Propagate updated bounds upward through the compact tree.
     * After tightening this edge's bounds, all ancestor edges
     * that reference this edge via their F-set need recomputation.
     */
    public void propagateBoundsUp() {
        RootedTreeEdge ancestor = compactParent;
        while (ancestor != null) {
            ancestor.recomputeFromChildren();
            ancestor = ancestor.compactParent;
        }
    }

    // ========== Conformation generation (for Phase 2 full minimization) ==========

    /**
     * Get the lambda-state index with the largest error contribution for a given M-state.
     * Used for top-down conformation generation in Phase 2.
     * Returns the lambda-state index (flat), or -1 if no error.
     */
    public int getBestErrorLambdaState(int mIdx) {
        if (!hasFsetChildren()) {
            // Leaf edge: best error = largest gap between E_rigid and E_min
            int k = enumeratedCount[mIdx];
            double bestError = 0.0;
            int bestLambdaIdx = -1;
            // Among enumerated states, find the one with largest gap
            for (int i = 0; i < k; i++) {
                int lIdx = sortedLambdaIndices != null ? sortedLambdaIndices[mIdx][i] : i;
                double eRigid = getFullEnergyRigid(mIdx, lIdx);
                double eMin = getFullEnergyMin(mIdx, lIdx);
                // Error in Boltzmann-weighted sense: exp(-eMin/RT) - exp(-eRigid/RT)
                double error = (-eMin / cachedRT) - (-eRigid / cachedRT);
                if (error > bestError) {
                    bestError = error;
                    bestLambdaIdx = lIdx;
                }
            }
            return bestLambdaIdx;
        } else {
            // Non-leaf edge: find lambda-state with largest contribution to error
            int k = enumeratedCount[mIdx];
            int[] mRCs = decodeMState(mIdx);
            double bestError = 0.0;
            int bestLambdaIdx = -1;
            for (int i = 0; i < k; i++) {
                int lIdx = sortedLambdaIndices != null ? sortedLambdaIndices[mIdx][i] : i;
                int[] lambdaRCs = decodeLambdaState(lIdx);

                double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

                double fSumLower = 0.0;
                double fSumUpper = 0.0;
                if (Fset != null) {
                    for (RootedTreeEdge fEdge : Fset) {
                        int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                        long fIndex = fEdge.computeIndexInA(fM);
                        fSumLower += fEdge.getLogZLower(fIndex);
                        fSumUpper += fEdge.getLogZUpper(fIndex);
                    }
                }

                double termUpper = -eMin / cachedRT + fSumUpper;
                double termLower = -eRigid / cachedRT + fSumLower;
                double error = termUpper - termLower;
                if (error > bestError) {
                    bestError = error;
                    bestLambdaIdx = lIdx;
                }
            }
            return bestLambdaIdx;
        }
    }

    /**
     * Decode a lambda-state index to global RC assignments for lambda positions.
     * Returns array where result[i] = global RC at lambdaPositionsSorted[i].
     */
    public int[] getLambdaGlobalRCs(int lambdaIdx) {
        int[] localRCs = decodeLambdaState(lambdaIdx);
        int[] globalRCs = new int[lambdaPositionsSorted.length];
        for (int i = 0; i < lambdaPositionsSorted.length; i++) {
            globalRCs[i] = rcs.get(lambdaPositionsSorted[i], localRCs[i]);
        }
        return globalRCs;
    }

    /**
     * Decode a M-state index to global RC assignments for M positions.
     */
    public int[] getMGlobalRCs(long mIdx) {
        int[] localRCs = decodeMState(mIdx);
        int[] globalRCs = new int[mPositionsSorted.length];
        for (int i = 0; i < mPositionsSorted.length; i++) {
            globalRCs[i] = rcs.get(mPositionsSorted[i], localRCs[i]);
        }
        return globalRCs;
    }

    /**
     * Get the rigid and min energies for a specific (mIdx, lambdaIdx) on a leaf edge.
     * Used by BranchDpBackend for correction computation.
     */
    public double getFullEnergyRigid(long mIdx, int lambdaIdx) {
        if (fullEnergyRigid != null) {
            return fullEnergyRigid[(int) mIdx][lambdaIdx];
        }
        return computeFullEnergy(mIdx, lambdaIdx, cachedRigidEmat, lambdaOnlyRigid);
    }

    public double getFullEnergyMin(long mIdx, int lambdaIdx) {
        if (fullEnergyMin != null) {
            return fullEnergyMin[(int) mIdx][lambdaIdx];
        }
        return computeFullEnergy(mIdx, lambdaIdx, cachedMinEmat, lambdaOnlyMin);
    }

    private double computeFullEnergy(long mIdx, int lambdaIdx,
                                     EnergyMatrix emat, double[] lambdaOnlyEnergies) {
        int[] mRCs = decodeMState(mIdx);
        int[] lambdaRCs = decodeLambdaState(lambdaIdx);
        double lambdaOnly = lambdaOnlyEnergies != null
                ? lambdaOnlyEnergies[lambdaIdx]
                : computeLambdaOnlyEnergy(lambdaRCs, emat, cachedG);
        return lambdaOnly + computeLambdaMEnergy(mRCs, lambdaRCs, emat, cachedG);
    }

    // ========== One-shot bottom-up DP ==========

    /**
     * Compute logZ bounds for all M-states by enumerating ALL lambda-states at once.
     * This replaces the incremental enumeration scheme: instead of enumerating 1 λ-state
     * at a time with propagateBoundsUp() after each step, we enumerate everything in one pass.
     *
     * Must be called in post-order: all F-set children must have their logZ computed first.
     *
     * After this call:
     *   logZLower[mIdx] = logSumExp over all λ-states of (-eRigid/RT + fSumLower)
     *   logZUpper[mIdx] = logSumExp over all λ-states of (-eMin/RT + fSumUpper)
     * No tail bound needed — all λ-states are enumerated.
     */
    public void computeFullDP() {
        if (!isLambdaEdge) return;

        // Build the read-only child-fold plan before any parallel dispatch.
        ensureChildFoldPlans();

        // Option A: resolve the native kernel toggle once. If requested but the
        // .so is missing/unloadable, log and fall back to the Java path so a
        // flag flip without a built library degrades gracefully (compile-only).
        nativeKernelEnabled = getConfigBoolean(DP_NATIVE_KERNEL_PROPERTY, false);
        if (nativeKernelEnabled && !ensureNativeKernelLoaded()) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " native DP kernel requested but libOspreyLogSumExp "
                    + "unavailable; falling back to Java DP path");
            nativeKernelEnabled = false;
        }

        if (tryComputeFullDPGpu()) {
            if (enumeratedCount != null) {
                Arrays.fill(enumeratedCount, totalLambdaStates);
            }
            return;
        }

        if (shouldParallelizeFullDP()) {
            computeFullDPSharded();
        } else {
            LongStream.range(0, mStateCount).forEach(this::computeFullDPForMState);
        }

        // Mark all M-states as fully enumerated (for Phase 2 conformation generation)
        if (enumeratedCount != null) {
            Arrays.fill(enumeratedCount, totalLambdaStates);
        }
    }

    private boolean shouldParallelizeFullDP() {
        return getConfigBoolean(DP_PARALLEL_PROPERTY, true)
                && mStateCount >= getConfigInteger(DP_PARALLEL_MIN_M_STATES_PROPERTY,
                        DEFAULT_DP_PARALLEL_MIN_M_STATES);
    }

    /**
     * Shard-scheduled parallel full DP (design 2.5).
     *
     * The M-state range [0, mStateCount) is split into contiguous shards of
     * {@code computeShardSize} M-states. A bounded, fixed-size thread pool pulls
     * shards off a shared atomic cursor (dynamic load balancing, so uneven shard
     * cost does not let one shard monopolize a worker). Contiguous shards keep
     * each worker's reads/writes local to a few backing arrays, which is friendly
     * to both cache and the ShardedDPTable layout.
     *
     * Safety: each M-state writes a distinct dpTable index, child edges are
     * already fully computed (post-order) and read-only here, and the log-sum-exp
     * accumulators are per-M-state locals, so contiguous-range partitioning has no
     * data races.
     */
    private void computeFullDPSharded() {
        int threads = resolveDPThreads();
        long shardSize = resolveComputeShardSize();
        long numShards = (mStateCount + shardSize - 1L) / shardSize;
        boolean progress = getConfigBoolean(DP_PROGRESS_PROPERTY, true);

        if (numShards <= 1 || threads <= 1) {
            LongStream.range(0, mStateCount).forEach(this::computeFullDPForMState);
            return;
        }

        int workers = (int) Math.min(threads, numShards);
        ExecutorService pool = Executors.newFixedThreadPool(workers,
                daemonThreadFactory(BranchDpConfig.getBackendThreadNamePrefix() + "-dp"));
        AtomicLong nextShard = new AtomicLong(0);
        AtomicLong completedShards = new AtomicLong(0);
        long startNanos = System.nanoTime();
        long logEvery = Math.max(1L, numShards / 20L);

        if (progress) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP shard schedule start, mStates=" + mStateCount
                    + ", shards=" + numShards + ", shardSize=" + shardSize
                    + ", workers=" + workers);
        }

        List<Future<?>> futures = new ArrayList<>(workers);
        for (int w = 0; w < workers; w++) {
            futures.add(pool.submit(() -> {
                long s;
                while ((s = nextShard.getAndIncrement()) < numShards) {
                    long shardStartNanos = System.nanoTime();
                    long mStart = s * shardSize;
                    long mEnd = Math.min(mStart + shardSize, mStateCount);
                    for (long mIdx = mStart; mIdx < mEnd; mIdx++) {
                        computeFullDPForMState(mIdx);
                    }
                    long done = completedShards.incrementAndGet();
                    if (progress && (done % logEvery == 0 || done == numShards)) {
                        double shardMs = (System.nanoTime() - shardStartNanos) / 1e6;
                        double elapsedS = (System.nanoTime() - startNanos) / 1e9;
                        System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP shards " + done + "/" + numShards
                                + " (" + String.format(Locale.ROOT, "%.1f", 100.0 * done / numShards)
                                + "%), lastShard=" + String.format(Locale.ROOT, "%.1f", shardMs)
                                + "ms, elapsed=" + String.format(Locale.ROOT, "%.1f", elapsedS) + "s");
                    }
                }
            }));
        }

        try {
            for (Future<?> f : futures) {
                f.get();
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("DP shard scheduling interrupted", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new RuntimeException("DP shard computation failed", cause);
        } finally {
            pool.shutdownNow();
        }

        if (progress) {
            double totalS = (System.nanoTime() - startNanos) / 1e9;
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP shard schedule done, shards=" + numShards
                    + ", elapsed=" + String.format(Locale.ROOT, "%.1f", totalS) + "s");
        }
    }

    private int resolveDPThreads() {
        int configured = getConfigInteger(DP_PARALLEL_THREADS_PROPERTY, 0);
        if (configured > 0) {
            return configured;
        }
        return Math.max(1, Runtime.getRuntime().availableProcessors());
    }

    private long resolveComputeShardSize() {
        long shardSize = getConfigInteger(DP_COMPUTE_SHARD_SIZE_PROPERTY, DEFAULT_DP_COMPUTE_SHARD_SIZE);
        return Math.max(1L, shardSize);
    }

    private boolean tryComputeFullDPGpu() {
        if (!getConfigBoolean(DP_GPU_PROPERTY, false)) {
            return false;
        }
        if (nativeKernelEnabled) {
            return false;
        }
        if (!hasFsetChildren() || childFoldPlans == null
                || childFoldPlans.length < 1
                || childFoldPlans.length > DP_GPU_MAX_CHILDREN) {
            return false;
        }
        if (dpTable == null) {
            return false;
        }
        if (mStateCount <= 0 || totalLambdaStates <= 0) {
            return false;
        }
        if (mPositionsSorted.length > DPGpuFullDP.MAX_EDGE_POSITIONS
                || lambdaPositionsSorted.length > DPGpuFullDP.MAX_EDGE_POSITIONS) {
            return false;
        }
        for (ChildFoldPlan plan : childFoldPlans) {
            if (plan == null || plan.child == null || !plan.child.hasDenseDPTable()
                    || plan.mSrcIdx.length > DPGpuFullDP.MAX_EDGE_POSITIONS
                    || plan.lamSrcIdx.length > DPGpuFullDP.MAX_EDGE_POSITIONS) {
                return false;
            }
        }
        if (lambdaOnlyRigid == null || lambdaOnlyMin == null
                || cachedRigidEmat == null || cachedMinEmat == null
                || cachedG == null || cachedRT == 0.0 || !Double.isFinite(cachedRT)) {
            return false;
        }
        if (mStateCount > Long.MAX_VALUE/(long)totalLambdaStates) {
            return false;
        }

        long work = mStateCount*(long)totalLambdaStates;
        long minWork = getConfigLong(DP_GPU_MIN_WORK_PROPERTY, DEFAULT_DP_GPU_MIN_WORK);
        if (work < minWork) {
            return false;
        }

        DPGpuFullDP.Request req = buildGpuFullDPRequest(work);
        if (req == null) {
            return false;
        }

        long maxBytes = getConfigLong(DP_GPU_MAX_BYTES_PROPERTY, DEFAULT_DP_GPU_MAX_BYTES);
        req.estimatedDeviceBytes = DPGpuFullDP.estimateDeviceBytes(req);
        boolean progress = getConfigBoolean(DP_PROGRESS_PROPERTY, true);
        req.progress = progress;
        // maxBytes <= 0 means "no static cap": defer to the detected-VRAM gate in
        // DPGpuFullDP. A positive cap is an explicit operator limit.
        if (maxBytes > 0 && req.estimatedDeviceBytes > maxBytes) {
            String msg = "GPU DP footprint " + req.estimatedDeviceBytes + " B exceeds "
                    + DP_GPU_MAX_BYTES_PROPERTY + "=" + maxBytes + " B (lambdaStates="
                    + totalLambdaStates + ", childTableElems=" + req.childTableTotal + ")";
            if (getConfigBoolean(DP_GPU_FAIL_IF_EXCEEDS_VRAM_PROPERTY, true)) {
                throw new DPGpuFullDP.GpuMemoryExceededException(msg
                        + " -- resubmit on a larger-VRAM GPU or raise " + DP_GPU_MAX_BYTES_PROPERTY);
            }
            if (progress) {
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP skipped, " + msg);
            }
            return false;
        }

        if (progress) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " GPU DP start, mStates=" + mStateCount
                    + ", lambdaStates=" + totalLambdaStates
                    + ", children=" + req.numChildren
                    + ", work=" + work
                    + ", lmPairs=" + req.lmLamSlots.length
                    + ", lmTerms=" + req.lmRigid.length
                    + ", outputTileMStates=" + req.mStateChunk
                    + ", outputTableDense=" + dpTable.isDenseArrayBacked()
                    + ", persistentContext=" + req.persistentContext
                    + ", estimatedDeviceBytes=" + req.estimatedDeviceBytes);
        }
        return DPGpuFullDP.compute(req);
    }

    private DPGpuFullDP.Request buildGpuFullDPRequest(long work) {
        long lmTermCount = 0L;
        List<int[]> lmPairs = new ArrayList<>();
        for (int li = 0; li < lambdaPositionsSorted.length; li++) {
            int lambdaPos = lambdaPositionsSorted[li];
            int lambdaCount = rcs.getNum(lambdaPos);
            for (int mi = 0; mi < mPositionsSorted.length; mi++) {
                int mPos = mPositionsSorted[mi];
                if (!cachedG.hasEdge(lambdaPos, mPos)) {
                    continue;
                }
                int mCount = rcs.getNum(mPos);
                long terms = (long)lambdaCount*(long)mCount;
                if (terms < 0 || lmTermCount > Integer.MAX_VALUE - terms) {
                    return null;
                }
                lmPairs.add(new int[]{li, mi, lambdaPos, mPos, lambdaCount, mCount});
                lmTermCount += terms;
            }
        }

        DPGpuFullDP.Request req = new DPGpuFullDP.Request();
        req.mStateCount = mStateCount;
        req.totalLambdaStates = totalLambdaStates;
        req.blockThreads = Math.max(1, getConfigInteger(DP_GPU_BLOCK_THREADS_PROPERTY,
                DEFAULT_DP_GPU_BLOCK_THREADS));
        req.invRT = 1.0/cachedRT;
        req.work = work;
        req.mCounts = countsForPositions(mPositionsSorted);
        req.lambdaCounts = countsForPositions(lambdaPositionsSorted);
        req.lambdaOnlyRigid = lambdaOnlyRigid;
        req.lambdaOnlyMin = lambdaOnlyMin;
        // Concatenate the per-child fold plans CSR-style for full_dp_n_children.
        int nChildren = childFoldPlans.length;
        req.numChildren = nChildren;
        req.mStateChunk = resolveGpuOutputTileMStates();
        req.outTable = dpTable;
        req.persistentContext = getConfigBoolean(DP_GPU_PERSISTENT_CONTEXT_PROPERTY,
                DEFAULT_DP_GPU_PERSISTENT_CONTEXT);
        req.multiGpu = getConfigBoolean(DP_GPU_MULTI_PROPERTY, true);
        req.maxGpus = getConfigInteger(DP_GPU_MAX_GPUS_PROPERTY, 0);
        req.minMStatesPerGpu = getConfigLong(DP_GPU_MIN_MSTATES_PER_GPU_PROPERTY,
                DEFAULT_DP_GPU_MIN_MSTATES_PER_GPU);

        long mTermTotal = 0L, lTermTotal = 0L, tableTotal = 0L;
        for (ChildFoldPlan plan : childFoldPlans) {
            mTermTotal += plan.mSrcIdx.length;
            lTermTotal += plan.lamSrcIdx.length;
            tableTotal += plan.child.getLogZLower().length;
        }
        // The CSR index arrays ARE concatenated into single Java int[] buffers, so they
        // must each fit in an int (in practice tiny: ~#positions). The child *tables* are
        // no longer concatenated on the Java heap -- each child keeps its own double[] and
        // is copied straight into its device-buffer slot -- so their combined length
        // (tableTotal) may exceed Integer.MAX_VALUE; only device VRAM bounds it (enforced
        // at upload, where an over-budget edge fails the run instead of crawling on CPU).
        if (mTermTotal > Integer.MAX_VALUE || lTermTotal > Integer.MAX_VALUE) {
            return null;
        }

        req.childMSrcAll = new int[(int)mTermTotal];
        req.childMStrideAll = new long[(int)mTermTotal];
        req.childMTermOff = new int[nChildren];
        req.childMTermCnt = new int[nChildren];
        req.childLSrcAll = new int[(int)lTermTotal];
        req.childLStrideAll = new long[(int)lTermTotal];
        req.childLTermOff = new int[nChildren];
        req.childLTermCnt = new int[nChildren];
        req.childTableBase = new long[nChildren];
        req.childTableTotal = tableTotal;
        req.childLowerParts = new double[nChildren][];
        req.childUpperParts = new double[nChildren][];

        int mOff = 0, lOff = 0;
        long tBase = 0L;
        for (int c = 0; c < nChildren; c++) {
            ChildFoldPlan plan = childFoldPlans[c];
            req.childMTermOff[c] = mOff;
            req.childMTermCnt[c] = plan.mSrcIdx.length;
            System.arraycopy(plan.mSrcIdx, 0, req.childMSrcAll, mOff, plan.mSrcIdx.length);
            System.arraycopy(plan.mStride, 0, req.childMStrideAll, mOff, plan.mStride.length);
            mOff += plan.mSrcIdx.length;

            req.childLTermOff[c] = lOff;
            req.childLTermCnt[c] = plan.lamSrcIdx.length;
            System.arraycopy(plan.lamSrcIdx, 0, req.childLSrcAll, lOff, plan.lamSrcIdx.length);
            System.arraycopy(plan.lamStride, 0, req.childLStrideAll, lOff, plan.lamStride.length);
            lOff += plan.lamSrcIdx.length;

            // Reference each child's dense table directly; copied per-child into the
            // device buffer at childTableBase[c] (a long element offset).
            req.childTableBase[c] = tBase;
            req.childLowerParts[c] = plan.child.getLogZLower();
            req.childUpperParts[c] = plan.child.getLogZUpper();
            tBase += req.childLowerParts[c].length;
        }

        int nPairs = lmPairs.size();
        req.lmLamSlots = new int[nPairs];
        req.lmMSlots = new int[nPairs];
        req.lmMCounts = new int[nPairs];
        req.lmOffsets = new long[nPairs];
        req.lmRigid = new double[(int)lmTermCount];
        req.lmMin = new double[(int)lmTermCount];

        int offset = 0;
        for (int p = 0; p < nPairs; p++) {
            int[] pair = lmPairs.get(p);
            int lambdaSlot = pair[0];
            int mSlot = pair[1];
            int lambdaPos = pair[2];
            int mPos = pair[3];
            int lambdaCount = pair[4];
            int mCount = pair[5];

            req.lmLamSlots[p] = lambdaSlot;
            req.lmMSlots[p] = mSlot;
            req.lmMCounts[p] = mCount;
            req.lmOffsets[p] = offset;

            for (int lrc = 0; lrc < lambdaCount; lrc++) {
                int globalLrc = rcs.get(lambdaPos, lrc);
                for (int mrc = 0; mrc < mCount; mrc++) {
                    int globalMrc = rcs.get(mPos, mrc);
                    int idx = offset + lrc*mCount + mrc;
                    req.lmRigid[idx] = cachedRigidEmat.getPairwise(lambdaPos, globalLrc, mPos, globalMrc);
                    req.lmMin[idx] = cachedMinEmat.getPairwise(lambdaPos, globalLrc, mPos, globalMrc);
                }
            }
            offset += (int)((long)lambdaCount*(long)mCount);
        }

        return req;
    }

    private long resolveGpuOutputTileMStates() {
        long requested = getConfigLong(DP_GPU_OUTPUT_TILE_MSTATES_PROPERTY,
                DEFAULT_DP_GPU_OUTPUT_TILE_MSTATES);
        long tile = Math.max(1L, requested);
        tile = Math.min(tile, (long)Integer.MAX_VALUE);
        return Math.min(tile, mStateCount);
    }

    // ===== PAC Phase-1 GPU sampling (Gumbel-max over the DP upper-bound weight) =====

    /** Structural gate: this lambda edge can be sampled by the GPU Gumbel kernel. */
    public boolean canUseGpuSampling() {
        if (!isLambdaEdge) {
            return false;
        }
        ensureChildFoldPlans();
        if (!hasFsetChildren() || childFoldPlans == null
                || childFoldPlans.length < 1
                || childFoldPlans.length > DP_GPU_MAX_CHILDREN) {
            return false;
        }
        if (totalLambdaStates <= 0) {
            return false;
        }
        if (mPositionsSorted.length > SamplingGpuPhase1.MAX_EDGE_POSITIONS
                || lambdaPositionsSorted.length > SamplingGpuPhase1.MAX_EDGE_POSITIONS) {
            return false;
        }
        for (ChildFoldPlan plan : childFoldPlans) {
            if (plan == null || plan.child == null || !plan.child.hasDenseDPTable()
                    || plan.mSrcIdx.length > SamplingGpuPhase1.MAX_EDGE_POSITIONS
                    || plan.lamSrcIdx.length > SamplingGpuPhase1.MAX_EDGE_POSITIONS) {
                return false;
            }
        }
        return lambdaOnlyMin != null && cachedMinEmat != null && cachedG != null
                && cachedRT != 0.0 && Double.isFinite(cachedRT);
    }

    /**
     * Draw one lambda-state index per sample on the GPU via Gumbel-max, for the
     * parent M-states in mIdxPerSample. Returns null on any failure so the caller
     * falls back to the Java sampler. Statistically (not bit-) equivalent.
     */
    public int[] sampleLambdaStatesGpu(long[] mIdxPerSample, long baseSeed, boolean progress) {
        if (mIdxPerSample == null || mIdxPerSample.length == 0) {
            return new int[0];
        }
        SamplingGpuPhase1.Request req = buildGpuSamplingRequest(mIdxPerSample, baseSeed, progress);
        if (req == null) {
            return null;
        }
        return SamplingGpuPhase1.sample(req);
    }

    /** Build the sampling request: the DP upper/min half + per-sample mIdx (no rigid/lower/dpTable). */
    private SamplingGpuPhase1.Request buildGpuSamplingRequest(long[] mIdxPerSample, long baseSeed, boolean progress) {
        long lmTermCount = 0L;
        List<int[]> lmPairs = new ArrayList<>();
        for (int li = 0; li < lambdaPositionsSorted.length; li++) {
            int lambdaPos = lambdaPositionsSorted[li];
            int lambdaCount = rcs.getNum(lambdaPos);
            for (int mi = 0; mi < mPositionsSorted.length; mi++) {
                int mPos = mPositionsSorted[mi];
                if (!cachedG.hasEdge(lambdaPos, mPos)) {
                    continue;
                }
                int mCount = rcs.getNum(mPos);
                long terms = (long)lambdaCount*(long)mCount;
                if (terms < 0 || lmTermCount > Integer.MAX_VALUE - terms) {
                    return null;
                }
                lmPairs.add(new int[]{li, mi, lambdaPos, mPos, lambdaCount, mCount});
                lmTermCount += terms;
            }
        }

        SamplingGpuPhase1.Request req = new SamplingGpuPhase1.Request();
        req.totalLambdaStates = totalLambdaStates;
        req.blockThreads = Math.max(1, getConfigInteger(DP_GPU_BLOCK_THREADS_PROPERTY,
                DEFAULT_DP_GPU_BLOCK_THREADS));
        req.invRT = 1.0/cachedRT;
        req.mCounts = countsForPositions(mPositionsSorted);
        req.lambdaCounts = countsForPositions(lambdaPositionsSorted);
        req.lambdaOnlyMin = lambdaOnlyMin;
        req.mIdxPerSample = mIdxPerSample;
        req.baseSeed = baseSeed;
        req.progress = progress;
        req.multiGpu = getConfigBoolean(PAC_SAMPLING_GPU_MULTI_PROPERTY, true);
        req.maxGpus = getConfigInteger(PAC_SAMPLING_GPU_MAX_GPUS_PROPERTY, 0);
        req.minGroupsPerGpu = Math.max(1, getConfigInteger(PAC_SAMPLING_GPU_MIN_GROUPS_PER_GPU_PROPERTY,
                DEFAULT_PAC_SAMPLING_GPU_MIN_GROUPS_PER_GPU));
        req.persistentContext = getConfigBoolean(PAC_SAMPLING_GPU_PERSISTENT_CONTEXT_PROPERTY,
                DEFAULT_PAC_SAMPLING_GPU_PERSISTENT_CONTEXT);
        req.residentChildTables = getConfigBoolean(PAC_SAMPLING_GPU_RESIDENT_CHILD_TABLES_PROPERTY,
                DEFAULT_PAC_SAMPLING_GPU_RESIDENT_CHILD_TABLES);
        req.method = SamplingGpuPhase1.Method.fromProperty(
                getConfigProperty(PAC_SAMPLING_GPU_METHOD_PROPERTY,
                        DEFAULT_PAC_SAMPLING_GPU_METHOD.propertyValue),
                DEFAULT_PAC_SAMPLING_GPU_METHOD);

        int nChildren = childFoldPlans.length;
        req.numChildren = nChildren;

        long mTermTotal = 0L, lTermTotal = 0L, tableTotal = 0L;
        for (ChildFoldPlan plan : childFoldPlans) {
            mTermTotal += plan.mSrcIdx.length;
            lTermTotal += plan.lamSrcIdx.length;
            tableTotal += plan.child.getLogZUpper().length;
        }
        if (mTermTotal > Integer.MAX_VALUE || lTermTotal > Integer.MAX_VALUE
                || tableTotal > Integer.MAX_VALUE) {
            return null;
        }

        req.childMSrcAll = new int[(int)mTermTotal];
        req.childMStrideAll = new long[(int)mTermTotal];
        req.childMTermOff = new int[nChildren];
        req.childMTermCnt = new int[nChildren];
        req.childLSrcAll = new int[(int)lTermTotal];
        req.childLStrideAll = new long[(int)lTermTotal];
        req.childLTermOff = new int[nChildren];
        req.childLTermCnt = new int[nChildren];
        req.childTableBase = new long[nChildren];
        req.childUpperAll = new double[(int)tableTotal];

        int mOff = 0, lOff = 0, tBase = 0;
        long childUpperKey = 0x6a09e667f3bcc909L;
        childUpperKey = mixGpuSamplingKey(childUpperKey, nChildren);
        childUpperKey = mixGpuSamplingKey(childUpperKey, tableTotal);
        for (int c = 0; c < nChildren; c++) {
            ChildFoldPlan plan = childFoldPlans[c];
            req.childMTermOff[c] = mOff;
            req.childMTermCnt[c] = plan.mSrcIdx.length;
            System.arraycopy(plan.mSrcIdx, 0, req.childMSrcAll, mOff, plan.mSrcIdx.length);
            System.arraycopy(plan.mStride, 0, req.childMStrideAll, mOff, plan.mStride.length);
            mOff += plan.mSrcIdx.length;

            req.childLTermOff[c] = lOff;
            req.childLTermCnt[c] = plan.lamSrcIdx.length;
            System.arraycopy(plan.lamSrcIdx, 0, req.childLSrcAll, lOff, plan.lamSrcIdx.length);
            System.arraycopy(plan.lamStride, 0, req.childLStrideAll, lOff, plan.lamStride.length);
            lOff += plan.lamSrcIdx.length;

            double[] cu = plan.child.getLogZUpper();
            req.childTableBase[c] = tBase;
            childUpperKey = mixGpuSamplingKey(childUpperKey, c);
            childUpperKey = mixGpuSamplingKey(childUpperKey, cu.length);
            for (int i = 0; i < cu.length; i++) {
                double value = cu[i];
                req.childUpperAll[tBase + i] = value;
                childUpperKey = mixGpuSamplingKey(childUpperKey, Double.doubleToRawLongBits(value));
            }
            tBase += cu.length;
        }
        req.childUpperCacheKey = childUpperKey;

        int nPairs = lmPairs.size();
        req.lmLamSlots = new int[nPairs];
        req.lmMSlots = new int[nPairs];
        req.lmMCounts = new int[nPairs];
        req.lmOffsets = new long[nPairs];
        req.lmMin = new double[(int)lmTermCount];

        int offset = 0;
        for (int p = 0; p < nPairs; p++) {
            int[] pair = lmPairs.get(p);
            int lambdaSlot = pair[0];
            int mSlot = pair[1];
            int lambdaPos = pair[2];
            int mPos = pair[3];
            int lambdaCount = pair[4];
            int mCount = pair[5];

            req.lmLamSlots[p] = lambdaSlot;
            req.lmMSlots[p] = mSlot;
            req.lmMCounts[p] = mCount;
            req.lmOffsets[p] = offset;

            for (int lrc = 0; lrc < lambdaCount; lrc++) {
                int globalLrc = rcs.get(lambdaPos, lrc);
                for (int mrc = 0; mrc < mCount; mrc++) {
                    int globalMrc = rcs.get(mPos, mrc);
                    int idx = offset + lrc*mCount + mrc;
                    req.lmMin[idx] = cachedMinEmat.getPairwise(lambdaPos, globalLrc, mPos, globalMrc);
                }
            }
            offset += (int)((long)lambdaCount*(long)mCount);
        }

        return req;
    }

    private static long mixGpuSamplingKey(long h, long x) {
        long z = x + 0x9E3779B97F4A7C15L + (h << 6) + (h >>> 2);
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return h ^ (z ^ (z >>> 31));
    }

    private int[] countsForPositions(int[] positions) {
        int[] counts = new int[positions.length];
        for (int i = 0; i < positions.length; i++) {
            counts[i] = rcs.getNum(positions[i]);
        }
        return counts;
    }

    private static ThreadFactory daemonThreadFactory(String namePrefix) {
        AtomicLong counter = new AtomicLong(0);
        return runnable -> {
            Thread t = new Thread(runnable, namePrefix + "-" + counter.getAndIncrement());
            t.setDaemon(true);
            return t;
        };
    }

    private void computeFullDPForMState(long mIdx) {
        if (nativeKernelEnabled) {
            computeFullDPForMStateNative(mIdx);
            return;
        }

        LogSumExpAccumulator lower = new LogSumExpAccumulator();
        LogSumExpAccumulator upper = new LogSumExpAccumulator();

        if (!hasFsetChildren()) {
            // Leaf edge: stream lambda states without materializing [M x lambda].
            if (fullEnergyRigid != null && fullEnergyMin != null) {
                // Materialized: read the contiguous per-M rows directly.
                double[] rigidRow = fullEnergyRigid[(int) mIdx];
                double[] minRow = fullEnergyMin[(int) mIdx];
                for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                    lower.add(-rigidRow[lIdx] / cachedRT);
                    upper.add(-minRow[lIdx] / cachedRT);
                }
            } else {
                // Streamed: hoist the M-state decode (constant across the lambda loop) and
                // decode each lambda-state ONCE for both rigid and min energy. The previous
                // getFullEnergyRigid/getFullEnergyMin path re-decoded the M-state twice and
                // the lambda-state twice per element. This is bit-for-bit identical (same
                // operations, same order) — it only removes redundant decode work/allocations.
                int[] mRCs = decodeMState(mIdx);
                for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                    int[] lambdaRCs = decodeLambdaState(lIdx);
                    double eRigid = (lambdaOnlyRigid != null
                            ? lambdaOnlyRigid[lIdx]
                            : computeLambdaOnlyEnergy(lambdaRCs, cachedRigidEmat, cachedG))
                            + computeLambdaMEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                    double eMin = (lambdaOnlyMin != null
                            ? lambdaOnlyMin[lIdx]
                            : computeLambdaOnlyEnergy(lambdaRCs, cachedMinEmat, cachedG))
                            + computeLambdaMEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);
                    lower.add(-eRigid / cachedRT);
                    upper.add(-eMin / cachedRT);
                }
            }
        } else if (childFoldPlans != null) {
            // Non-leaf edge (Option C): fold in child bounds with a precomputed plan.
            // The plan replaces the per-element getMstateForFullState() int[]
            // allocation + linear search with a direct mixed-radix index, and
            // hoists each child's M-only ("base") index out of the lambda loop.
            int[] mRCs = decodeMState(mIdx);
            int nChildren = childFoldPlans.length;

            long[] baseIdx = new long[nChildren];
            double constLower = 0.0;
            double constUpper = 0.0;
            for (int c = 0; c < nChildren; c++) {
                ChildFoldPlan plan = childFoldPlans[c];
                baseIdx[c] = plan.baseIndex(mRCs);
                if (childFoldHoistInvariant && !plan.lambdaDependent) {
                    // λ-invariant child: constant fIdx across the whole lambda
                    // loop, so gather its bounds once. NOTE: this changes the
                    // FP add order vs the legacy/unfolded path, so it is not
                    // strictly bit-identical (validate like Step 1).
                    constLower += plan.child.getLogZLower(baseIdx[c]);
                    constUpper += plan.child.getLogZUpper(baseIdx[c]);
                }
            }

            for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                int[] lambdaRCs = decodeLambdaState(lIdx);
                double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

                double fSumLower = constLower;
                double fSumUpper = constUpper;
                for (int c = 0; c < nChildren; c++) {
                    ChildFoldPlan plan = childFoldPlans[c];
                    if (childFoldHoistInvariant && !plan.lambdaDependent) {
                        continue; // already folded into constLower/constUpper
                    }
                    long fIdx = baseIdx[c] + plan.lambdaIndex(lambdaRCs);
                    fSumLower += plan.child.getLogZLower(fIdx);
                    fSumUpper += plan.child.getLogZUpper(fIdx);
                }

                lower.add(-eRigid / cachedRT + fSumLower);
                upper.add(-eMin / cachedRT + fSumUpper);
            }
        } else {
            // Non-leaf edge (legacy path): fold disabled or unmapped child slot.
            int[] mRCs = decodeMState(mIdx);
            for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                int[] lambdaRCs = decodeLambdaState(lIdx);
                double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

                double fSumLower = 0.0;
                double fSumUpper = 0.0;
                for (RootedTreeEdge fEdge : Fset) {
                    int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                    long fIdx = fEdge.computeIndexInA(fM);
                    fSumLower += fEdge.getLogZLower(fIdx);
                    fSumUpper += fEdge.getLogZUpper(fIdx);
                }

                lower.add(-eRigid / cachedRT + fSumLower);
                upper.add(-eMin / cachedRT + fSumUpper);
            }
        }

        dpTable.set(mIdx, lower.value(), upper.value());
    }

    /**
     * Option A: native-kernel variant of {@link #computeFullDPForMState}. Reuses
     * the exact same per-lambda value expressions (and the Option C fold plan for
     * non-leaf edges), but materializes the per-row lower/upper value vectors into
     * reusable scratch buffers and offloads the log-sum-exp reduction to the
     * native SIMD kernel. One fused native call per M-state amortizes JNA overhead.
     *
     * NOT bit-identical to the streaming-accumulator path (different summation
     * order in the reduction) — this is the documented Option A caveat; validate
     * bounds before trusting. Default OFF.
     */
    private void computeFullDPForMStateNative(long mIdx) {
        int n = totalLambdaStates;
        double[] vLower = NATIVE_LOWER_BUF.get();
        double[] vUpper = NATIVE_UPPER_BUF.get();
        if (vLower.length < n) { vLower = new double[n]; NATIVE_LOWER_BUF.set(vLower); }
        if (vUpper.length < n) { vUpper = new double[n]; NATIVE_UPPER_BUF.set(vUpper); }

        if (!hasFsetChildren()) {
            if (fullEnergyRigid != null && fullEnergyMin != null) {
                double[] rigidRow = fullEnergyRigid[(int) mIdx];
                double[] minRow = fullEnergyMin[(int) mIdx];
                for (int lIdx = 0; lIdx < n; lIdx++) {
                    vLower[lIdx] = -rigidRow[lIdx] / cachedRT;
                    vUpper[lIdx] = -minRow[lIdx] / cachedRT;
                }
            } else {
                int[] mRCs = decodeMState(mIdx);
                for (int lIdx = 0; lIdx < n; lIdx++) {
                    int[] lambdaRCs = decodeLambdaState(lIdx);
                    double eRigid = (lambdaOnlyRigid != null
                            ? lambdaOnlyRigid[lIdx]
                            : computeLambdaOnlyEnergy(lambdaRCs, cachedRigidEmat, cachedG))
                            + computeLambdaMEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                    double eMin = (lambdaOnlyMin != null
                            ? lambdaOnlyMin[lIdx]
                            : computeLambdaOnlyEnergy(lambdaRCs, cachedMinEmat, cachedG))
                            + computeLambdaMEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);
                    vLower[lIdx] = -eRigid / cachedRT;
                    vUpper[lIdx] = -eMin / cachedRT;
                }
            }
        } else {
            int[] mRCs = decodeMState(mIdx);
            if (childFoldPlans != null) {
                int nChildren = childFoldPlans.length;
                long[] baseIdx = new long[nChildren];
                for (int c = 0; c < nChildren; c++) {
                    baseIdx[c] = childFoldPlans[c].baseIndex(mRCs);
                }
                for (int lIdx = 0; lIdx < n; lIdx++) {
                    int[] lambdaRCs = decodeLambdaState(lIdx);
                    double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                    double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);
                    double fSumLower = 0.0;
                    double fSumUpper = 0.0;
                    for (int c = 0; c < nChildren; c++) {
                        ChildFoldPlan plan = childFoldPlans[c];
                        long fIdx = baseIdx[c] + plan.lambdaIndex(lambdaRCs);
                        fSumLower += plan.child.getLogZLower(fIdx);
                        fSumUpper += plan.child.getLogZUpper(fIdx);
                    }
                    vLower[lIdx] = -eRigid / cachedRT + fSumLower;
                    vUpper[lIdx] = -eMin / cachedRT + fSumUpper;
                }
            } else {
                for (int lIdx = 0; lIdx < n; lIdx++) {
                    int[] lambdaRCs = decodeLambdaState(lIdx);
                    double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                    double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);
                    double fSumLower = 0.0;
                    double fSumUpper = 0.0;
                    for (RootedTreeEdge fEdge : Fset) {
                        int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                        long fIdx = fEdge.computeIndexInA(fM);
                        fSumLower += fEdge.getLogZLower(fIdx);
                        fSumUpper += fEdge.getLogZUpper(fIdx);
                    }
                    vLower[lIdx] = -eRigid / cachedRT + fSumLower;
                    vUpper[lIdx] = -eMin / cachedRT + fSumUpper;
                }
            }
        }

        double[] out = NATIVE_OUT_BUF.get();
        NativeLogSumExp.osprey_logsumexp2_f64(vLower, vUpper, n, out);
        dpTable.set(mIdx, out[0], out[1]);
    }

    /**
     * JNA binding for the Option A native kernel (libOspreyLogSumExp). Registered
     * lazily and defensively via {@link #ensureNativeKernelLoaded()} so a missing
     * .so cannot poison class init — it just disables the native path.
     */
    private static final class NativeLogSumExp {
        static native int osprey_logsumexp_version();
        static native double osprey_logsumexp_f64(double[] vals, long n);
        static native void osprey_logsumexp2_f64(double[] lower, double[] upper, long n, double[] out2);
    }

    private static volatile boolean nativeKernelRegistered = false;
    private static volatile boolean nativeKernelUnavailable = false;

    private static synchronized boolean ensureNativeKernelLoaded() {
        if (nativeKernelRegistered) return true;
        if (nativeKernelUnavailable) return false;
        try {
            Native.register(NativeLogSumExp.class, "OspreyLogSumExp");
            NativeLogSumExp.osprey_logsumexp_version(); // probe the symbol
            nativeKernelRegistered = true;
            return true;
        } catch (Throwable t) {
            nativeKernelUnavailable = true;
            return false;
        }
    }

    private static class LogSumExpAccumulator {
        private double max = NEG_INF;
        private double scaledSum = 0.0;

        void add(double value) {
            if (Double.isNaN(value) || Double.isNaN(max)) {
                max = Double.NaN;
                scaledSum = Double.NaN;
                return;
            }
            if (value == NEG_INF) return;
            if (max == NEG_INF) {
                max = value;
                scaledSum = 1.0;
            } else if (value > max) {
                scaledSum = scaledSum * Math.exp(max - value) + 1.0;
                max = value;
            } else {
                scaledSum += Math.exp(value - max);
            }
        }

        double value() {
            if (max == NEG_INF) return NEG_INF;
            return max + Math.log(scaledSum);
        }
    }

    /**
     * Direct independent-position DP for a 0-edge interaction graph (no pairwise
     * terms). Z then factorizes per position, so
     *   logZ = sum_pos logsumexp_rc( -E(pos, rc) / RT ),
     * computed exactly in O(sum cardinalities) instead of the O(prod cardinalities)
     * generic joint leaf DP. Returns {logZLower (rigid emat), logZUpper (min emat)},
     * matching the leaf-DP convention (lower uses rigid, upper uses min).
     *
     * Validated against brute-force joint enumeration in TestZeroEdgeDP. Wired
     * behind {@code branchdp.dp.zeroEdgeDirect} (default off) at the 0-edge
     * fallback in BranchDpBackend.
     */
    public static double[] independentPositionLogZ(EnergyMatrix rigidEmat, EnergyMatrix minEmat,
                                                   RCs rcs, int[] positions, double RT) {
        double logZLower = 0.0;
        double logZUpper = 0.0;
        for (int pos : positions) {
            LogSumExpAccumulator lo = new LogSumExpAccumulator();
            LogSumExpAccumulator hi = new LogSumExpAccumulator();
            int n = rcs.getNum(pos);
            for (int i = 0; i < n; i++) {
                int rc = rcs.get(pos, i);
                lo.add(-rigidEmat.getOneBody(pos, rc) / RT);
                hi.add(-minEmat.getOneBody(pos, rc) / RT);
            }
            logZLower += lo.value();
            logZUpper += hi.value();
        }
        return new double[]{logZLower, logZUpper};
    }

    /**
     * Precomputed plan to fold one child edge's logZ bounds into this (non-leaf)
     * edge's DP. Replaces the per-element getMstateForFullState() linear search +
     * int[] allocation with a direct mixed-radix index built from the parent's
     * already-decoded M/lambda RC arrays.
     *
     * For child M-slot i, computeIndexInA() weights it by
     *   stride[i] = product_{k>i} rcs.getNum(child.mPositionsSorted[k]).
     * By the branch-decomposition invariant each child M-position is one of the
     * parent's M-positions or lambda-positions, so the index splits into an
     * M-only "base" (constant across the lambda loop) plus a lambda-dependent
     * part. A child with no lambda-sourced slots has a constant index across the
     * whole lambda loop ("λ-invariant") and can be hoisted out entirely.
     */
    private static final class ChildFoldPlan {
        final RootedTreeEdge child;
        final int[] mSrcIdx;    // parent mRCs index for each M-sourced child slot
        final long[] mStride;   // mixed-radix stride for each M-sourced child slot
        final int[] lamSrcIdx;  // parent lambdaRCs index for each λ-sourced child slot
        final long[] lamStride; // mixed-radix stride for each λ-sourced child slot
        final boolean lambdaDependent;

        ChildFoldPlan(RootedTreeEdge child, int[] mSrcIdx, long[] mStride,
                      int[] lamSrcIdx, long[] lamStride) {
            this.child = child;
            this.mSrcIdx = mSrcIdx;
            this.mStride = mStride;
            this.lamSrcIdx = lamSrcIdx;
            this.lamStride = lamStride;
            this.lambdaDependent = lamStride.length > 0;
        }

        long baseIndex(int[] mRCs) {
            long idx = 0L;
            for (int t = 0; t < mStride.length; t++) {
                idx += (long) mRCs[mSrcIdx[t]] * mStride[t];
            }
            return idx;
        }

        long lambdaIndex(int[] lambdaRCs) {
            long idx = 0L;
            for (int t = 0; t < lamStride.length; t++) {
                idx += (long) lambdaRCs[lamSrcIdx[t]] * lamStride[t];
            }
            return idx;
        }
    }

    /**
     * Build the child-fold plans once, before the (possibly parallel) DP loop.
     * The result is read-only afterwards, so it is safe to share across DP
     * worker threads. Disabled (null plans => legacy path) when
     * {@code branchdp.dp.foldChildren=false}, when there are no F-set
     * children, or when any child M-position is not covered by this edge's
     * M ∪ lambda (in which case the legacy getMstateForFullState() path is used
     * unchanged, preserving behavior on malformed/edge-case decompositions).
     */
    private void ensureChildFoldPlans() {
        childFoldHoistInvariant = getConfigBoolean(DP_FOLD_HOIST_PROPERTY, false);
        if (!getConfigBoolean(DP_FOLD_CHILDREN_PROPERTY, true) || !hasFsetChildren()) {
            childFoldPlans = null;
            return;
        }

        ChildFoldPlan[] plans = new ChildFoldPlan[Fset.size()];
        int ci = 0;
        for (RootedTreeEdge fEdge : Fset) {
            ChildFoldPlan plan = buildChildFoldPlan(fEdge);
            if (plan == null) {
                childFoldPlans = null; // unmapped slot -> fall back to legacy path
                return;
            }
            plans[ci++] = plan;
        }
        childFoldPlans = plans;

        if (getConfigBoolean(DP_PROGRESS_PROPERTY, true)) {
            int invariant = 0;
            for (ChildFoldPlan p : plans) {
                if (!p.lambdaDependent) invariant++;
            }
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " childFold lambdaStates=" + totalLambdaStates
                    + ", mStates=" + mStateCount + ", children=" + plans.length
                    + ", lambdaInvariant=" + invariant
                    + ", lambdaDependent=" + (plans.length - invariant)
                    + ", hoist=" + childFoldHoistInvariant);
        }
    }

    private ChildFoldPlan buildChildFoldPlan(RootedTreeEdge child) {
        int[] childM = child.mPositionsSorted;
        int len = childM.length;

        // Mixed-radix stride per child M-slot (matches child.computeIndexInA):
        // stride[len-1] = 1; stride[i] = product_{k>i} rcs.getNum(childM[k]).
        long[] stride = new long[len];
        if (len > 0) {
            stride[len - 1] = 1L;
            for (int i = len - 2; i >= 0; i--) {
                stride[i] = stride[i + 1] * child.rcs.getNum(childM[i + 1]);
            }
        }

        int[] mSrcTmp = new int[len];
        long[] mStrideTmp = new long[len];
        int mCount = 0;
        int[] lamSrcTmp = new int[len];
        long[] lamStrideTmp = new long[len];
        int lamCount = 0;

        for (int i = 0; i < len; i++) {
            int targetPos = childM[i];
            // Search parent M first, then parent lambda (matches getMstateForFullState).
            int j = indexOf(mPositionsSorted, targetPos);
            if (j >= 0) {
                mSrcTmp[mCount] = j;
                mStrideTmp[mCount] = stride[i];
                mCount++;
                continue;
            }
            j = indexOf(lambdaPositionsSorted, targetPos);
            if (j >= 0) {
                lamSrcTmp[lamCount] = j;
                lamStrideTmp[lamCount] = stride[i];
                lamCount++;
                continue;
            }
            return null; // child M-position outside parent M ∪ lambda -> bail to legacy
        }

        return new ChildFoldPlan(child,
                Arrays.copyOf(mSrcTmp, mCount), Arrays.copyOf(mStrideTmp, mCount),
                Arrays.copyOf(lamSrcTmp, lamCount), Arrays.copyOf(lamStrideTmp, lamCount));
    }

    private static int indexOf(int[] arr, int value) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == value) return i;
        }
        return -1;
    }

    /**
     * Post-order traversal: compute full DP for all lambda-edges.
     * Leaves are computed first, then their parents fold in the leaf bounds.
     */
    public static void postOrderComputeFullDP(RootedTreeNode node) {
        if (node == null) return;
        postOrderComputeFullDP(node.getLeftChild());
        postOrderComputeFullDP(node.getRightChild());
        if (node.getChildOfEdge() != null && node.getChildOfEdge().isLambdaEdge) {
            node.getChildOfEdge().computeFullDP();
        }
    }

    // ========== Post-order traversal helpers ==========

    public static void postOrderCompLlambda(RootedTreeNode node) {
        postOrderCompLlambda(node, true);
    }

    public static void postOrderCompLlambda(RootedTreeNode node, boolean initEnumerationArrays) {
        if (node == null) return;
        postOrderCompLlambda(node.getLeftChild(), initEnumerationArrays);
        postOrderCompLlambda(node.getRightChild(), initEnumerationArrays);
        if (node.getChildOfEdge() != null) {
            node.getChildOfEdge().compLlambda(initEnumerationArrays);
        }
    }

    /**
     * Post-order traversal: initialize incremental enumeration for all lambda-edges.
     * After this call, all edges have k=0 with loose initial bounds.
     */
    public static void postOrderInitIncremental(RootedTreeNode node,
                                                  EnergyMatrix rigidEmat, EnergyMatrix minEmat,
                                                  InteractionGraph G, double RT) {
        postOrderInitIncremental(node, rigidEmat, minEmat, G, RT, false);
    }

    public static void postOrderInitIncremental(RootedTreeNode node,
                                                  EnergyMatrix rigidEmat, EnergyMatrix minEmat,
                                                  InteractionGraph G, double RT,
                                                  boolean materializeFullEnergyTables) {
        if (node == null) return;
        postOrderInitIncremental(node.getLeftChild(), rigidEmat, minEmat, G, RT,
                materializeFullEnergyTables);
        postOrderInitIncremental(node.getRightChild(), rigidEmat, minEmat, G, RT,
                materializeFullEnergyTables);

        if (node.getChildOfEdge() != null && node.getChildOfEdge().isLambdaEdge) {
            node.getChildOfEdge().initIncrementalEnumeration(rigidEmat, minEmat, G, RT,
                    materializeFullEnergyTables);
        }
    }

    /**
     * Collect all lambda-edges in the tree (post-order).
     */
    public static void collectLambdaEdges(RootedTreeNode node, List<RootedTreeEdge> result) {
        if (node == null) return;
        collectLambdaEdges(node.getLeftChild(), result);
        collectLambdaEdges(node.getRightChild(), result);
        if (node.getChildOfEdge() != null && node.getChildOfEdge().isLambdaEdge) {
            result.add(node.getChildOfEdge());
        }
    }
}
