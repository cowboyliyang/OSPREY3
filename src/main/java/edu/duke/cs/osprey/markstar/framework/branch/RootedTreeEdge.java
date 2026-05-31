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

package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import java.util.*;
import java.util.stream.IntStream;

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
    private int mArraySize;

    // ========== Log-space partition function bounds per M-state ==========

    private double[] logZLower;   // log(Z_lower) indexed by M-state
    private double[] logZUpper;   // log(Z_upper) indexed by M-state

    private static final double NEG_INF = Double.NEGATIVE_INFINITY;
    private static final String DP_PARALLEL_PROPERTY = "branchmarkstar.dp.parallel";
    private static final String DP_PARALLEL_MIN_M_STATES_PROPERTY = "branchmarkstar.dp.parallel.minMStates";
    private static final String DP_MAX_M_STATES_PROPERTY = "branchmarkstar.dp.maxMStates";
    private static final String DP_MAX_MATERIALIZED_PAIRS_PROPERTY = "branchmarkstar.dp.maxMaterializedPairs";
    private static final int DEFAULT_DP_PARALLEL_MIN_M_STATES = 1024;
    private static final long DEFAULT_DP_MAX_MATERIALIZED_PAIRS = 50_000_000L;

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
    public double[] getLogZLower() { return logZLower; }
    public double[] getLogZUpper() { return logZUpper; }
    public RootedTreeNode getCompactLeftChild() { return compactLeftChild; }
    public RootedTreeNode getCompactRightChild() { return compactRightChild; }
    public void setCompactParent(RootedTreeEdge p) { this.compactParent = p; }
    public RootedTreeEdge getCompactParent() { return compactParent; }
    public int getMArraySize() { return mArraySize; }
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
        String value = System.getProperty(key);
        if (value == null) {
            value = System.getProperty("osprey." + key);
        }
        return value != null ? value : defaultValue;
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        String value = getConfigProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    private static int getConfigInteger(String key, int defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("BranchMARK*: Invalid integer for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static long getConfigLong(String key, long defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Long.parseLong(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("BranchMARK*: Invalid long for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
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

        mArraySize = checkedStateCount("M", mPositionsSorted, true);
        totalLambdaStates = checkedStateCount("lambda", lambdaPositionsSorted, false);
    }

    private int checkedStateCount(String label, int[] positions, boolean isMStateCount) {
        long count = 1L;
        for (int pos : positions) {
            int n = rcs.getNum(pos);
            if (n < 0 || count > Long.MAX_VALUE / Math.max(1, n)) {
                throw stateCountException(label, positions, Long.MAX_VALUE);
            }
            count *= n;
            if (count > Integer.MAX_VALUE) {
                throw stateCountException(label, positions, count);
            }
        }

        long maxStates = isMStateCount
                ? getConfigLong(DP_MAX_M_STATES_PROPERTY, Integer.MAX_VALUE)
                : Integer.MAX_VALUE;
        if (count > maxStates) {
            throw new DPTableTooLargeException(label, count, positions,
                    "DP " + label + "-state count exceeds "
                            + DP_MAX_M_STATES_PROPERTY + "=" + maxStates
                            + ". Reduce branchwidth/rooting, enable a shard-backed DP table, "
                            + "or raise the limit only if memory allows.");
        }
        return (int) count;
    }

    private DPTableTooLargeException stateCountException(String label, int[] positions, long count) {
        return new DPTableTooLargeException(label, count, positions,
                "DP " + label + "-state count exceeds Java array indexing limits. "
                        + "This used to overflow into NegativeArraySizeException; "
                        + "use a smaller separator or a shard-backed DP table.");
    }

    private void initializeArrays() {
        logZLower = new double[mArraySize];
        logZUpper = new double[mArraySize];
        Arrays.fill(logZLower, NEG_INF);
        Arrays.fill(logZUpper, NEG_INF);

        enumeratedCount = new int[mArraySize];
        // sortedLambdaIndices allocated lazily per M-state in initIncrementalEnumeration
    }

    // ========== M-state indexing ==========

    public int computeIndexInA(int[] mRCIndices) {
        if (isRootEdge || mRCIndices == null || mPositionsSorted.length == 0) {
            return 0;
        }

        int index = mRCIndices[mPositionsSorted.length - 1];
        int stride = 1;
        for (int i = mPositionsSorted.length - 2; i >= 0; i--) {
            stride *= rcs.getNum(mPositionsSorted[i + 1]);
            index += mRCIndices[i] * stride;
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

    public int[] decodeMStatePublic(int mIndex) {
        return decodeMState(mIndex);
    }

    public int[] decodeLambdaStatePublic(int lambdaIndex) {
        return decodeLambdaState(lambdaIndex);
    }

    private int[] decodeMState(int mIndex) {
        int[] mRCs = new int[mPositionsSorted.length];
        int remaining = mIndex;
        for (int i = mPositionsSorted.length - 1; i >= 0; i--) {
            int numRCs = rcs.getNum(mPositionsSorted[i]);
            mRCs[i] = remaining % numRCs;
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
        Arrays.fill(logZLower, NEG_INF);
        Arrays.fill(logZUpper, NEG_INF);
        if (materializeFullEnergyTables) {
            // logZ_upper = tail bound for all K lambda-states in incremental mode.
            for (int mIdx = 0; mIdx < mArraySize; mIdx++) {
                logZUpper[mIdx] = computeTailBound(mIdx, 0);
            }
        }
    }

    private void checkMaterializedPairCount() {
        long pairs = (long) mArraySize * (long) totalLambdaStates;
        long maxPairs = getConfigLong(DP_MAX_MATERIALIZED_PAIRS_PROPERTY,
                DEFAULT_DP_MAX_MATERIALIZED_PAIRS);
        if (pairs > maxPairs) {
            throw new IllegalStateException("BranchMARK*: Refusing to materialize DP full-energy table with "
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
                    for (int fi = 0; fi < fEdge.mArraySize; fi++) {
                        bestChildUpper = Math.max(bestChildUpper, fEdge.logZUpper[fi]);
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
            throw new IllegalStateException("BranchMARK*: Incremental lambda enumeration requires materialized DP energy tables. "
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

            logZLower[mIdx] = logSumExp(logZLower[mIdx], logTermLower);

            // Upper = enumerated part + new tail
            // Recompute upper from scratch: enumerated terms + tail
            // (simpler and avoids accumulation errors)
            double enumeratedUpper = NEG_INF;
            for (int i = 0; i <= k; i++) {
                int li = sortedLambdaIndices[mIdx][i];
                enumeratedUpper = logSumExp(enumeratedUpper, -fullEnergyMin[mIdx][li] / cachedRT);
            }
            double tailUpper = computeTailBound(mIdx, k + 1);
            logZUpper[mIdx] = logSumExp(enumeratedUpper, tailUpper);

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
                    int fIndex = fEdge.computeIndexInA(fM);
                    fSumLower += fEdge.logZLower[fIndex];
                    fSumUpper += fEdge.logZUpper[fIndex];
                }
            }

            double logTermLower = -eRigid / cachedRT + fSumLower;
            double logTermUpper = -eMin / cachedRT + fSumUpper;

            logZLower[mIdx] = logSumExp(logZLower[mIdx], logTermLower);

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
                        int fIndex = fEdge.computeIndexInA(fM);
                        liFSumUpper += fEdge.logZUpper[fIndex];
                    }
                }
                enumeratedUpper = logSumExp(enumeratedUpper, -liEMin / cachedRT + liFSumUpper);
            }
            double tailUpper = computeTailBound(mIdx, k + 1);
            logZUpper[mIdx] = logSumExp(enumeratedUpper, tailUpper);
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
        if (logZUpper[mIdx] == NEG_INF && logZLower[mIdx] == NEG_INF) return 0.0;
        if (logZLower[mIdx] == NEG_INF) return Double.MAX_VALUE;
        return logZUpper[mIdx] - logZLower[mIdx];
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
                        int fIndex = fEdge.computeIndexInA(fM);
                        fSumLower += fEdge.logZLower[fIndex];
                        fSumUpper += fEdge.logZUpper[fIndex];
                    }
                }

                newLower = logSumExp(newLower, -eRigid / cachedRT + fSumLower);
                newUpper = logSumExp(newUpper, -eMin / cachedRT + fSumUpper);
            }

            // Add tail bound
            double tail = computeTailBound(mIdx, k);
            newUpper = logSumExp(newUpper, tail);

            logZLower[mIdx] = newLower;
            logZUpper[mIdx] = newUpper;
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
                        int fIndex = fEdge.computeIndexInA(fM);
                        fSumLower += fEdge.logZLower[fIndex];
                        fSumUpper += fEdge.logZUpper[fIndex];
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
    public int[] getMGlobalRCs(int mIdx) {
        int[] localRCs = decodeMState(mIdx);
        int[] globalRCs = new int[mPositionsSorted.length];
        for (int i = 0; i < mPositionsSorted.length; i++) {
            globalRCs[i] = rcs.get(mPositionsSorted[i], localRCs[i]);
        }
        return globalRCs;
    }

    /**
     * Get the rigid and min energies for a specific (mIdx, lambdaIdx) on a leaf edge.
     * Used by BranchMARKStarBound for correction computation.
     */
    public double getFullEnergyRigid(int mIdx, int lambdaIdx) {
        if (fullEnergyRigid != null) {
            return fullEnergyRigid[mIdx][lambdaIdx];
        }
        return computeFullEnergy(mIdx, lambdaIdx, cachedRigidEmat, lambdaOnlyRigid);
    }

    public double getFullEnergyMin(int mIdx, int lambdaIdx) {
        if (fullEnergyMin != null) {
            return fullEnergyMin[mIdx][lambdaIdx];
        }
        return computeFullEnergy(mIdx, lambdaIdx, cachedMinEmat, lambdaOnlyMin);
    }

    private double computeFullEnergy(int mIdx, int lambdaIdx,
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

        IntStream mStates = IntStream.range(0, mArraySize);
        if (shouldParallelizeFullDP()) {
            mStates = mStates.parallel();
        }
        mStates.forEach(this::computeFullDPForMState);

        // Mark all M-states as fully enumerated (for Phase 2 conformation generation)
        if (enumeratedCount != null) {
            Arrays.fill(enumeratedCount, totalLambdaStates);
        }
    }

    private boolean shouldParallelizeFullDP() {
        return getConfigBoolean(DP_PARALLEL_PROPERTY, true)
                && mArraySize >= getConfigInteger(DP_PARALLEL_MIN_M_STATES_PROPERTY,
                        DEFAULT_DP_PARALLEL_MIN_M_STATES);
    }

    private void computeFullDPForMState(int mIdx) {
        LogSumExpAccumulator lower = new LogSumExpAccumulator();
        LogSumExpAccumulator upper = new LogSumExpAccumulator();

        if (!hasFsetChildren()) {
            // Leaf edge: stream lambda states without materializing [M x lambda].
            for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                lower.add(-getFullEnergyRigid(mIdx, lIdx) / cachedRT);
                upper.add(-getFullEnergyMin(mIdx, lIdx) / cachedRT);
            }
        } else {
            // Non-leaf edge: fold in child bounds (children already computed via post-order).
            int[] mRCs = decodeMState(mIdx);
            for (int lIdx = 0; lIdx < totalLambdaStates; lIdx++) {
                int[] lambdaRCs = decodeLambdaState(lIdx);
                double eRigid = computeLocalEnergy(mRCs, lambdaRCs, cachedRigidEmat, cachedG);
                double eMin = computeLocalEnergy(mRCs, lambdaRCs, cachedMinEmat, cachedG);

                double fSumLower = 0.0;
                double fSumUpper = 0.0;
                for (RootedTreeEdge fEdge : Fset) {
                    int[] fM = getMstateForFullState(mRCs, lambdaRCs, fEdge);
                    int fIdx = fEdge.computeIndexInA(fM);
                    fSumLower += fEdge.logZLower[fIdx];
                    fSumUpper += fEdge.logZUpper[fIdx];
                }

                lower.add(-eRigid / cachedRT + fSumLower);
                upper.add(-eMin / cachedRT + fSumUpper);
            }
        }

        logZLower[mIdx] = lower.value();
        logZUpper[mIdx] = upper.value();
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
