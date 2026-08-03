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

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.markstar.framework.GridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Shared branch-DP partition function backend.
 *
 * Top-down edge-level traversal of the decomposition tree.
 * Each "expand" assigns an entire edge's lambda-set.
 * A priority queue decides between expand (internal) and minimize (leaf).
 * Z bounds maintained incrementally from h-scores (no DP needed).
 * O(1) bound updates per operation.
 */
public abstract class BranchDpBackend extends MARKStarBound {

    // Configuration
    private static final double DEFAULT_DIST_CUTOFF = 8.0;      // Angstroms
    private static final double DEFAULT_ENERGY_CUTOFF = 0.1;    // kcal/mol
    private static final double DEFAULT_RESIDUAL_BUDGET = 0.5;  // kcal/mol
    private static final double TESS_FALLBACK_THRESHOLD = 0.5;  // fallback if TESS/Naive > this
    private static final String CUTOFF_STRATEGY_PROPERTY = "branchdp.cutoff.strategy";
    private static final String DIST_CUTOFF_PROPERTY = "branchdp.cutoff.dist";
    private static final String ENERGY_CUTOFF_PROPERTY = "branchdp.cutoff.energy";
    private static final String RESIDUAL_BUDGET_PROPERTY = "branchdp.cutoff.residualBudget";
    private static final String CUTOFF_KEEP_CONNECTED_PROPERTY = "branchdp.cutoff.keepConnected";
    private static final String CERTIFY_FULL_BOUNDS_PROPERTY = "branchdp.certifyFullBounds";
    private static final String DECOMP_STRATEGY_PROPERTY = "branchdp.decomp.strategy";
    private static final String DECOMP_STRATEGY_SHORT_PROPERTY = "branchdp.decomp";
    private static final String ADAPTIVE_GPU_WORK_THRESHOLD_PROPERTY =
            "branchdp.decomp.adaptive.gpuWorkThreshold";
    private static final BigInteger DEFAULT_ADAPTIVE_GPU_WORK_THRESHOLD =
            new BigInteger("100000000000000");
    private static final String ADMISSION_MAX_GPU_WORK_PROPERTY =
            "branchdp.admission.maxGpuWork";
    private static final String EXACT_TREEWIDTH_PROPERTY = "branchdp.decomp.exactTreewidth";
    private static final String EXACT_TREEWIDTH_MAX_POSITIONS_PROPERTY = "branchdp.decomp.exactTreewidth.maxPositions";
    private static final String EDGE_SELECTION_PROPERTY = "branchdp.edgeSelection";
    private static final String EDGE_LOOKAHEAD_MAX_STATES_PROPERTY = "branchdp.edgeSelection.maxExactStates";
    private static final String EDGE_LOOKAHEAD_MAX_PENDING_EDGES_PROPERTY = "branchdp.edgeSelection.maxPendingEdges";
    private static final String EDGE_LOOKAHEAD_PARALLEL_PROPERTY = "branchdp.edgeSelection.parallelLookahead";
    private static final String ROOT_SPLIT_PROPERTY = "branchdp.rootSplit";
    private static final String ROOT_SPLIT_GPU_BUDGET_PROPERTY = "branchdp.rootSplit.gpuBudgetBytes";
    private static final String ROOT_SPLIT_HOST_BUDGET_PROPERTY = "branchdp.rootSplit.hostBudgetBytes";
    private static final double DEFAULT_ROOT_SPLIT_HOST_HEAP_FRACTION = 0.70;
    // 0-edge (no pairwise) graphs: exact independent-position DP. Default off.
    private static final String ZERO_EDGE_DIRECT_PROPERTY = "branchdp.dp.zeroEdgeDirect";
    private static final String ROOT_SPLIT_MAX_FSET_PROPERTY = "branchdp.rootSplit.maxFset";
    private static final String DRY_RUN_PROPERTY = "branchdp.dp.dryRun";
    private static final String DP_GPU_OUTPUT_TILE_MSTATES_PROPERTY = "branchdp.dp.gpu.outputTileMStates";
    private static final String MUTABLE_POSITIONS_PROPERTY = "branchdp.mutablePositions";
    private static final String PARALLEL_INTERNAL_PROPERTY = "branchdp.parallel.internal";
    private static final String PARALLEL_ENUMERATION_PROPERTY = "branchdp.parallel.enumeration";
    private static final String NUMERIC_AUDIT_PROPERTY = "branchdp.numericAudit";
    private static final String NUMERIC_AUDIT_INTERVAL_PROPERTY = "branchdp.numericAudit.interval";
    private static final String TRACE_PROPERTY = "branchdp.trace";
    private static final String TRACE_ROUNDS_PROPERTY = "branchdp.trace.rounds";
    private static final String TRACE_CORRECTIONS_PROPERTY = "branchdp.trace.corrections";
    private static final String CORRECTION_AUDIT_PROPERTY = "branchdp.correctionAudit";
    private static final String ENERGY_MODE_PROPERTY = "branchdp.energyMode";
    private static final String SPARSE_PRUNE_THRESHOLD_PROPERTY = "branchdp.sparse.pruneThreshold";
    private static final String DP_CACHE_ENABLED_PROPERTY = "branchdp.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "branchdp.dp.cache.maxEntries";
    private static final String DP_CACHE_MAX_TABLE_BYTES_PROPERTY = "branchdp.dp.cache.maxTableBytes";
    private static final String DP_CACHE_MAX_TOTAL_BYTES_PROPERTY = "branchdp.dp.cache.maxTotalBytes";
    private static final String DP_CACHE_SKIP_IF_M_STATES_PROPERTY = "branchdp.dp.cache.skipIfMStates";
    private static final int DEFAULT_DP_CACHE_MAX_ENTRIES = 20000;
    private static final long DEFAULT_DP_CACHE_MAX_TABLE_BYTES = 256L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_MAX_TOTAL_BYTES = 4L * 1024L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_SKIP_IF_M_STATES = 8_000_000L;
    private static final long DEFAULT_DP_GPU_OUTPUT_TILE_MSTATES = 1_048_576L;
    private static final int DEFAULT_EXACT_TREEWIDTH_MAX_POSITIONS = 22;
    private static final BigInteger SLICED_TRAFFIC_UNAVAILABLE = BigInteger.ONE.shiftLeft(1024);

    protected enum EdgeSelectionStrategy {
        LAMBDA_STATES,
        CONTRACTION,
        CONTRACTION_PER_STATE;

        static EdgeSelectionStrategy fromProperty(String value) {
            if (value == null || value.trim().isEmpty()) {
                return LAMBDA_STATES;
            }
            switch (value.trim().toLowerCase(Locale.ROOT)) {
                case "lambdastates":
                case "lambda_states":
                case "states":
                case "legacy":
                    return LAMBDA_STATES;
                case "contraction":
                case "drop":
                case "certificate":
                    return CONTRACTION;
                case "contractionperstate":
                case "contraction_per_state":
                case "dropperstate":
                case "drop_per_state":
                case "certificate_per_state":
                    return CONTRACTION_PER_STATE;
                default:
                    System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown edge selection strategy '" + value
                            + "', using lambdaStates.");
                    return LAMBDA_STATES;
            }
        }
    }

    private enum CutoffStrategy {
        RESIDUAL_BUDGET,
        DUAL,
        COMPLETE;

        static CutoffStrategy fromProperty(String value) {
            if (value == null || value.trim().isEmpty()) {
                return RESIDUAL_BUDGET;
            }
            switch (value.trim().toLowerCase(Locale.ROOT)) {
                case "residual":
                case "residualbudget":
                case "residual_budget":
                case "budget":
                    return RESIDUAL_BUDGET;
                case "dual":
                case "dualcutoff":
                case "dual_cutoff":
                case "legacy":
                    return DUAL;
                case "complete":
                case "full":
                case "none":
                    return COMPLETE;
                default:
                    System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown cutoff strategy '" + value
                            + "', using residualBudget.");
                    return RESIDUAL_BUDGET;
            }
        }
    }

    /**
     * Pfunc energy model:
     *  - SPARSE: strict sparse pfunc. The branch interaction graph defines
     *    BOTH the energy function (cut pairs contribute 0) AND the conformation
     *    set: conformations whose cut-edge pair |E| exceeds the sparse prune
     *    threshold are removed from the pfunc sum (pruned in the emat to +inf,
     *    so A* and leaf minimization skip them). This makes Z a well-defined
     *    pfunc over a clear sub-conformation space C_sparse with a clear
     *    sparse energy E_sparse, NOT a heuristic approximation that lets
     *    "cut-but-large-interaction" full confs leak in. Triple corrections
     *    are disabled because they are not certified against E_sparse /
     *    C_sparse.
     *  - FULL: legacy behavior. Sparse graph is only used to drive branch
     *    decomposition; pfunc sums all original confs with the full emat.
     */
    private enum EnergyMode {
        SPARSE,
        FULL;

        static EnergyMode fromProperty(String value) {
            if (value == null || value.trim().isEmpty()) {
                return SPARSE;
            }
            switch (value.trim().toLowerCase(Locale.ROOT)) {
                case "sparse":
                case "graph":
                case "sparse_graph":
                    return SPARSE;
                case "full":
                case "original":
                    return FULL;
                default:
                    System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown energy mode '" + value
                            + "', using sparse.");
                    return SPARSE;
            }
        }
    }

    private final CutoffStrategy cutoffStrategy;
    private final double distCutoff;
    private final double energyCutoff;
    private final double residualBudget;
    private final boolean cutoffKeepConnected;
    private final boolean certifyFullBounds;
    private final BranchDecomposition.Strategy decompStrategy;
    private final boolean exactTreewidthDiagnostic;
    private final int exactTreewidthMaxPositions;
    private final EdgeSelectionStrategy edgeSelectionStrategy;
    private final int edgeLookaheadMaxStates;
    private final int edgeLookaheadMaxPendingEdges;
    private final boolean edgeLookaheadParallel;
    private final String rootSplitStrategy;
    private final int rootSplitMaxFset;
    private final boolean parallelInternal;
    private final boolean parallelEnumeration;
    private final boolean numericAudit;
    private final int numericAuditInterval;
    private final boolean trace;
    private final int traceRounds;
    private final int traceCorrectionLimit;
    private final boolean correctionAudit;
    private final boolean useHigherOrderCorrections;
    private final EnergyMode energyMode;
    private final boolean useSparsePfunc;
    private final double sparsePruneThreshold;
    private final boolean dpCacheEnabled;
    private final int dpCacheMaxEntries;
    private final long dpCacheMaxTableBytes;
    private final long dpCacheMaxTotalBytes;
    private final long dpCacheSkipIfMStates;
    private final boolean dryRun;
    private final BranchDpAdmission.Hardware admissionHardware;
    private final String admissionStateKey;
    private final int admissionDpSweeps;
    private BranchDpAdmission.Prediction admissionPrediction = null;
    private boolean largeMemoryReleased = false;
    private long selectedRootHostBudgetBytes = Long.MAX_VALUE;
    private long selectedRootGpuBudgetBytes = Long.MAX_VALUE;
    private DPTableTooLargeException dpTooLargeException = null;
    /**
     * Energy matrices used by ALL pfunc-relevant operations under SPARSE mode:
     * scorers, A*, computeFullConfPairwiseEnergy, breakdownScoreByPosition,
     * tuple internal-energy queries, and leaf pairwise bounds. In SPARSE mode
     * these are copies of {@code rigidEmat}/{@code minimizingEmat} with cut
     * edges processed by {@link #makeSparseEnergyMatrix} (zeroed if accepted,
     * +inf if pruned). In FULL mode they alias the originals so behavior is
     * identical to the pre-sparse build.
     */
    protected EnergyMatrix branchRigidEmat;
    protected EnergyMatrix branchMinimizingEmat;
    private int sparsePrunedRcPairs = 0;
    private int sparseZeroedRcPairs = 0;
    private int sparseTotalCutPositionPairs = 0;

    /** True once {@link RootedTreeEdge#postOrderComputeFullDP} has populated logZ tables.
     * Only set in SPARSE mode; controls whether {@link #enhanceWithDPBounds} is applied. */
    protected boolean dpTablesReady = false;
    /** RCs used to build the tree / DP tables; stashed for use by enhanceWithDPBounds. */
    protected RCs searchRCs;

    private static final Object DP_TABLE_CACHE_LOCK = new Object();
    private static final LinkedHashMap<String, CachedDPTable> DP_TABLE_CACHE =
            new LinkedHashMap<>(1024, 0.75f, true);
    private static long dpTableCacheBytes = 0L;

    private static class CachedDPTable {
        final double[] lower;
        final double[] upper;
        final long bytes;

        CachedDPTable(double[] lower, double[] upper, long bytes) {
            this.lower = lower;
            this.upper = upper;
            this.bytes = bytes;
        }
    }

    private static class DPCacheStats {
        int hits = 0;
        int misses = 0;
        int stores = 0;
        int skippedLarge = 0;
        int evictions = 0;
    }

    // Whether branch decomposition is active
    private boolean useBranchDecomposition = true;

    // Conformation space
    protected final SimpleConfSpace confSpace;

    // Branch decomposition state
    protected InteractionGraph interactionGraph;
    private BranchDecomposition branchDecomposition;
    private int branchwidth;

    // Rooted tree
    protected RootedTreeNode rootedRoot;
    protected RootedTreeEdge rootedRootEdge;

    // Statistics
    private int totalEnumerationSteps = 0;
    protected int totalMinimizations = 0;
    private int totalDrillDowns = 0;
    private int totalDrillDownEarlyStops = 0;

    private long edgeSelectionCalls = 0;
    private long edgeSelectionSingleEdgeCalls = 0;
    private long edgeSelectionMultiEdgeCalls = 0;
    private long edgeSelectionLambdaStatesCalls = 0;
    private long edgeSelectionExactLookaheadCalls = 0;
    private long edgeSelectionFallbackCalls = 0;
    private long edgeSelectionFallbackMaxStates = 0;
    private long edgeLookaheadCandidates = 0;
    private long edgeLookaheadEnumeratedStates = 0;
    private long edgeLookaheadCacheHits = 0;
    private long edgeLookaheadCacheMisses = 0;
    private long edgeLookaheadParallelCalls = 0;
    private long edgeLookaheadNanos = 0;
    private long edgeLookaheadSelectedStates = 0;
    private BigDecimal edgeLookaheadSelectedDrop = BigDecimal.ZERO;
    private int selectedRootSplitIndex = 0;
    private int maxPendingEdgesSeen = 0;
    private final Map<Integer, Long> pendingEdgeCountHistogram = new TreeMap<>();
    private long numericAuditChecks = 0;
    private long numericAuditWarnings = 0;
    private long correctionChecks = 0;
    private long correctionNotTighter = 0;
    private long correctionApplications = 0;
    private long correctionSparseGenerationSkips = 0;
    private long correctionGenerationCalls = 0;
    private long correctionLargePairs = 0;
    private long correctionTupleCandidates = 0;
    private long correctionTupleBelowThreshold = 0;
    private long correctionTupleAlreadyKnown = 0;
    private long correctionAsyncSubmitted = 0;
    private long correctionGenerated = 0;
    private long correctionZeroGenerated = 0;

    private final Map<DecompSearchNode, Map<Integer, EdgeLookaheadResult>> edgeLookaheadCache =
            Collections.synchronizedMap(new IdentityHashMap<>());

    // Priority queue for decomposition-guided search
    private PriorityQueue<DecompSearchNode> decompQueue;

    // Running Z bounds
    // Z_lower = sum of all nodes' subtreeLowerBound
    // Z_upper = sum of all nodes' subtreeUpperBound
    protected BigDecimal flatSumZLower = BigDecimal.ZERO;
    protected BigDecimal flatSumZUpper = BigDecimal.ZERO;

    // Scorers (created once, reused)
    private AStarScorer gScorerMin;
    private AStarScorer hScorerMin;
    private AStarScorer gScorerRigid;
    private AStarScorer hScorerNegRigid;

    // Set of already-minimized conformations to avoid duplicate work
    private final Set<String> minimizedConfs = new HashSet<>();

    private int roundCounter = 0;

    private static class EdgeLookaheadResult {
        final BigDecimal childErrorSum;
        final int totalStates;
        final int enumeratedStates;
        final long nanos;
        final boolean cacheHit;

        EdgeLookaheadResult(BigDecimal childErrorSum, int totalStates, int enumeratedStates,
                            long nanos, boolean cacheHit) {
            this.childErrorSum = childErrorSum;
            this.totalStates = totalStates;
            this.enumeratedStates = enumeratedStates;
            this.nanos = nanos;
            this.cacheHit = cacheHit;
        }

        EdgeLookaheadResult asCacheHit() {
            return new EdgeLookaheadResult(childErrorSum, totalStates, enumeratedStates, 0, true);
        }
    }

    private static class ChildEnumerationBatch {
        final List<DecompSearchNode> children;
        final BigDecimal lowerSum;
        final BigDecimal upperSum;
        final DecompSearchNode bestChild;

        ChildEnumerationBatch(List<DecompSearchNode> children, BigDecimal lowerSum,
                              BigDecimal upperSum, DecompSearchNode bestChild) {
            this.children = children;
            this.lowerSum = lowerSum;
            this.upperSum = upperSum;
            this.bestChild = bestChild;
        }
    }

    private static class ChildEnumerationChunk {
        final int chunkIndex;
        final List<DecompSearchNode> children;
        final BigDecimal lowerSum;
        final BigDecimal upperSum;
        final DecompSearchNode bestChild;

        ChildEnumerationChunk(int chunkIndex, List<DecompSearchNode> children,
                              BigDecimal lowerSum, BigDecimal upperSum,
                              DecompSearchNode bestChild) {
            this.chunkIndex = chunkIndex;
            this.children = children;
            this.lowerSum = lowerSum;
            this.upperSum = upperSum;
            this.bestChild = bestChild;
        }
    }

    private static class InternalNodeResult {
        final int nodeIndex;
        final List<DecompSearchNode> generatedNodes;
        final BigDecimal lowerDelta;
        final BigDecimal upperDelta;
        final boolean reachedTarget;
        final int enumerationSteps;
        final int drillDowns;
        final int drillDownEarlyStops;

        InternalNodeResult(int nodeIndex, List<DecompSearchNode> generatedNodes,
                           BigDecimal lowerDelta, BigDecimal upperDelta,
                           boolean reachedTarget, int enumerationSteps,
                           int drillDowns, int drillDownEarlyStops) {
            this.nodeIndex = nodeIndex;
            this.generatedNodes = generatedNodes;
            this.lowerDelta = lowerDelta;
            this.upperDelta = upperDelta;
            this.reachedTarget = reachedTarget;
            this.enumerationSteps = enumerationSteps;
            this.drillDowns = drillDowns;
            this.drillDownEarlyStops = drillDownEarlyStops;
        }
    }

    private static class RootingCandidate {
        final RootedTreeNode root;
        final int splitEdgeIndex;
        final double logTESS;
        final int lambdaEdges;
        final int maxFsetSize;
        final int branchingEdges;
        final int totalFsetEdges;
        final int rootFsetSize;
        final int reusableLambdaEdges;
        final double reusableLogWork;
        final double totalLogWork;
        final long maxMStates;
        final long totalMStates;
        final long maxDPTableBytes;
        final long totalDPTableBytes;
        // Estimated live Java-heap storage for DP tables plus per-edge
        // enumeration/lambda arrays. This is the hard root-feasibility metric;
        // DP-table bytes alone omit up to four bytes per int-addressable M state.
        final long maxHostBytes;
        final long totalHostBytes;
        final long maxFileBackedBytes;
        final long totalFileBackedBytes;
        final double logDPWork;
        // Exact count of dominant CUDA lambda-loop iterations across all edges.
        // BigInteger makes the comparison exact even when the total exceeds a long.
        final BigInteger gpuWork;
        final int gpuUnsupportedEdges;
        // Worst full-resident device footprint across lambda edges. Unlike the final
        // DP table metrics above, this counts child inputs and only a bounded output
        // tile, matching DPGpuFullDP's allocation model.
        final long maxSingleGpuBytes;
        final boolean fitsSingleGpu;
        final BigInteger estimatedSlicedTraffic;

        RootingCandidate(RootedTreeNode root, int splitEdgeIndex, double logTESS,
                         int lambdaEdges, int maxFsetSize, int branchingEdges,
                         int totalFsetEdges, int rootFsetSize,
                         int reusableLambdaEdges, double reusableLogWork,
                         double totalLogWork, long maxMStates, long totalMStates,
                         long maxDPTableBytes, long totalDPTableBytes,
                         long maxHostBytes, long totalHostBytes,
                         long maxFileBackedBytes, long totalFileBackedBytes,
                         double logDPWork, BigInteger gpuWork, int gpuUnsupportedEdges,
                         long maxSingleGpuBytes, boolean fitsSingleGpu,
                         BigInteger estimatedSlicedTraffic) {
            this.root = root;
            this.splitEdgeIndex = splitEdgeIndex;
            this.logTESS = logTESS;
            this.lambdaEdges = lambdaEdges;
            this.maxFsetSize = maxFsetSize;
            this.branchingEdges = branchingEdges;
            this.totalFsetEdges = totalFsetEdges;
            this.rootFsetSize = rootFsetSize;
            this.reusableLambdaEdges = reusableLambdaEdges;
            this.reusableLogWork = reusableLogWork;
            this.totalLogWork = totalLogWork;
            this.maxMStates = maxMStates;
            this.totalMStates = totalMStates;
            this.maxDPTableBytes = maxDPTableBytes;
            this.totalDPTableBytes = totalDPTableBytes;
            this.maxHostBytes = maxHostBytes;
            this.totalHostBytes = totalHostBytes;
            this.maxFileBackedBytes = maxFileBackedBytes;
            this.totalFileBackedBytes = totalFileBackedBytes;
            this.logDPWork = logDPWork;
            this.gpuWork = gpuWork;
            this.gpuUnsupportedEdges = gpuUnsupportedEdges;
            this.maxSingleGpuBytes = maxSingleGpuBytes;
            this.fitsSingleGpu = fitsSingleGpu;
            this.estimatedSlicedTraffic = estimatedSlicedTraffic;
        }

        double reusableEdgeRatio() {
            return lambdaEdges == 0 ? 0.0 : (double) reusableLambdaEdges / (double) lambdaEdges;
        }

        double reusableWorkRatio() {
            if (reusableLogWork == Double.NEGATIVE_INFINITY) return 0.0;
            if (totalLogWork == Double.NEGATIVE_INFINITY) return 0.0;
            return Math.exp(reusableLogWork - totalLogWork);
        }
    }

    // ========== Constructor ==========

    protected BranchDpBackend(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism) {
        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism, null);
    }

    protected BranchDpBackend(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism, String stateNameOverride) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        if (stateNameOverride != null && !stateNameOverride.trim().isEmpty()) {
            this.stateName = stateNameOverride;
        }
        this.confSpace = confSpace;
        this.admissionHardware = BranchDpAdmission.Hardware.fromBackendConfig();
        this.admissionStateKey = BranchDpAdmission.stateKey(this.stateName, rcs);
        this.cutoffStrategy = CutoffStrategy.fromProperty(
                getConfigProperty(CUTOFF_STRATEGY_PROPERTY, null));
        this.distCutoff = Math.max(0.0,
                getConfigDouble(DIST_CUTOFF_PROPERTY, DEFAULT_DIST_CUTOFF));
        this.energyCutoff = Math.max(0.0,
                getConfigDouble(ENERGY_CUTOFF_PROPERTY, DEFAULT_ENERGY_CUTOFF));
        this.residualBudget = Math.max(0.0,
                getConfigDouble(RESIDUAL_BUDGET_PROPERTY, DEFAULT_RESIDUAL_BUDGET));
        this.cutoffKeepConnected = getConfigBoolean(
                CUTOFF_KEEP_CONNECTED_PROPERTY, true);
        this.energyMode = EnergyMode.fromProperty(getConfigProperty(ENERGY_MODE_PROPERTY, null));
        this.useSparsePfunc = energyMode == EnergyMode.SPARSE;
        // Default prune threshold = +inf, i.e. NEVER prune any (rc1, rc2) cut-pair config.
        // This keeps the search's conformation set identical to the DP's: the DP iterates
        // only graph pairs (so it would silently treat any +inf cut entry as 0), whereas
        // gScorer iterates ALL defined pairs (so a +inf cut entry would push that child's
        // Z to 0). Letting them disagree means parent_DP_Z > Σ children_g_Z and flatSumZ
        // drains to 0 under expansion. Setting cut pairs to plain 0 keeps both consistent
        // and yields the clean Z_sparse pfunc. Users can still set a finite threshold
        // explicitly via {@link #SPARSE_PRUNE_THRESHOLD_PROPERTY}, with the caveat that
        // the DP's bound would then no longer match the search's pfunc target.
        this.sparsePruneThreshold = Math.max(0.0,
                getConfigDouble(SPARSE_PRUNE_THRESHOLD_PROPERTY, Double.POSITIVE_INFINITY));
        // certifyFullBounds inflates Z by exp(rho/RT) on the assumption that we
        // are summing over the FULL conformation set with a sparse energy. In SPARSE
        // mode we deliberately changed the conformation set to C_sparse, so that
        // inflation no longer reflects a valid sparse->full bound. Force it off.
        boolean requestedCertify = getConfigBoolean(CERTIFY_FULL_BOUNDS_PROPERTY, false);
        if (useSparsePfunc && requestedCertify) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " certifyFullBounds requested but disabled in SPARSE mode "
                    + "(sparse pfunc is a self-consistent Z over C_sparse, not an approximation of full Z).");
        }
        this.certifyFullBounds = requestedCertify && !useSparsePfunc;
        String decompStrategyValue = getConfigProperty(DECOMP_STRATEGY_PROPERTY, null);
        if (decompStrategyValue == null) {
            decompStrategyValue = getConfigProperty(DECOMP_STRATEGY_SHORT_PROPERTY, null);
        }
        this.decompStrategy = BranchDecomposition.Strategy.fromProperty(decompStrategyValue);
        this.exactTreewidthDiagnostic = getConfigBoolean(EXACT_TREEWIDTH_PROPERTY, false);
        this.exactTreewidthMaxPositions = Math.max(0,
                getConfigInteger(EXACT_TREEWIDTH_MAX_POSITIONS_PROPERTY, DEFAULT_EXACT_TREEWIDTH_MAX_POSITIONS));
        this.edgeSelectionStrategy = getEdgeSelectionStrategy();
        this.edgeLookaheadMaxStates = getEdgeLookaheadMaxStates();
        this.edgeLookaheadMaxPendingEdges = getEdgeLookaheadMaxPendingEdges();
        this.edgeLookaheadParallel = getEdgeLookaheadParallel();
        this.rootSplitStrategy = getConfigProperty(ROOT_SPLIT_PROPERTY, "work").trim();
        this.dryRun = getConfigBoolean(DRY_RUN_PROPERTY, false);
        this.admissionDpSweeps = Math.max(1, getAdmissionDpSweeps());
        this.rootSplitMaxFset = Math.max(1,
                getConfigInteger(ROOT_SPLIT_MAX_FSET_PROPERTY, 2));
        this.parallelInternal = getConfigBoolean(PARALLEL_INTERNAL_PROPERTY, true);
        this.parallelEnumeration = getConfigBoolean(PARALLEL_ENUMERATION_PROPERTY, true);
        this.numericAudit = getConfigBoolean(NUMERIC_AUDIT_PROPERTY, true);
        this.numericAuditInterval = Math.max(1,
                getConfigInteger(NUMERIC_AUDIT_INTERVAL_PROPERTY, 1000));
        this.trace = getConfigBoolean(TRACE_PROPERTY, false);
        this.traceRounds = Math.max(0,
                getConfigInteger(TRACE_ROUNDS_PROPERTY, 20));
        this.traceCorrectionLimit = Math.max(0,
                getConfigInteger(TRACE_CORRECTIONS_PROPERTY, 80));
        this.correctionAudit = getConfigBoolean(CORRECTION_AUDIT_PROPERTY, false);
        this.useHigherOrderCorrections = getUseHigherOrderCorrections();
        this.dpCacheEnabled = getConfigBoolean(DP_CACHE_ENABLED_PROPERTY, true);
        this.dpCacheMaxEntries = Math.max(0,
                getConfigInteger(DP_CACHE_MAX_ENTRIES_PROPERTY, DEFAULT_DP_CACHE_MAX_ENTRIES));
        this.dpCacheMaxTableBytes = Math.max(0L,
                getConfigBytes(DP_CACHE_MAX_TABLE_BYTES_PROPERTY, DEFAULT_DP_CACHE_MAX_TABLE_BYTES));
        this.dpCacheMaxTotalBytes = Math.max(0L,
                getConfigBytes(DP_CACHE_MAX_TOTAL_BYTES_PROPERTY, DEFAULT_DP_CACHE_MAX_TOTAL_BYTES));
        this.dpCacheSkipIfMStates = Math.max(0L,
                getConfigLong(DP_CACHE_SKIP_IF_M_STATES_PROPERTY, DEFAULT_DP_CACHE_SKIP_IF_M_STATES));

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Cutoff strategy=" + cutoffStrategy
                + ", distCutoff=" + distCutoff
                + ", energyCutoff=" + energyCutoff
                + ", residualBudget=" + residualBudget
                + ", keepConnected=" + cutoffKeepConnected
                + ", certifyFullBounds=" + certifyFullBounds
                + ", energyMode=" + energyMode
                + ", sparsePruneThreshold=" + String.format("%.6f", sparsePruneThreshold));
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Edge selection=" + edgeSelectionStrategy
                + ", decompStrategy=" + decompStrategy
                + ", exactTreewidthDiagnostic=" + exactTreewidthDiagnostic
                + ", maxExactLookaheadStates=" + edgeLookaheadMaxStates
                + ", maxPendingLookaheadEdges=" + edgeLookaheadMaxPendingEdges
                + ", parallelLookahead=" + edgeLookaheadParallel
                + ", rootSplit=" + rootSplitStrategy
                + ", rootSplitMaxFset=" + rootSplitMaxFset
                + ", tripleCorrections=" + useHigherOrderCorrections
                + ", correctionAudit=" + correctionAudit
                + ", dpCache=" + dpCacheEnabled
                + ", dpCacheMaxTable=" + formatBytes(dpCacheMaxTableBytes)
                + ", dpCacheMaxTotal=" + formatBytes(dpCacheMaxTotalBytes)
                + ", dpCacheSkipIfMStates=" + dpCacheSkipIfMStates
                + ", parallelInternal=" + parallelInternal + " [DEPRECATED]"
                + ", parallelEnumeration=" + parallelEnumeration + " [DEPRECATED]");
        logBackendControlOverrides();
        if (trace) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " trace enabled"
                    + " (rounds=" + traceRounds
                    + ", corrections=" + traceCorrectionLimit + ")");
        }
        if (parallelInternal) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " parallelInternal is DEPRECATED; prefer pending-edge lookahead parallelism.");
        }
        if (parallelEnumeration) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " parallelEnumeration is DEPRECATED; prefer pending-edge lookahead parallelism.");
        }

        // Step 1: Build sparse interaction graph. The residual-budget rule is
        // preferred because it gives a direct sparse-to-full perturbation bound.
        interactionGraph = buildInteractionGraph(rcs);
        printSparseResidualInfo();
        printExactTreewidthDiagnostic();

        // Step 1b: In SPARSE mode, derive the branch-pfunc energy matrices.
        // Cut-edge (rc1, rc2) pairs with |E| > pruneThreshold get +inf (pruned out
        // of the conf set); the rest get 0 (accepted sparse approximation).
        // FULL mode keeps the originals so behavior matches the legacy build.
        if (useSparsePfunc) {
            branchRigidEmat = makeSparseEnergyMatrix("rigid", rigidEmat, interactionGraph,
                    sparsePruneThreshold, true);
            branchMinimizingEmat = makeSparseEnergyMatrix("minimizing", minimizingEmat, interactionGraph,
                    sparsePruneThreshold, false);
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Strict SPARSE pfunc active. "
                    + "Cut-edge confs with |pair E| > " + String.format("%.6f", sparsePruneThreshold)
                    + " kcal/mol are pruned out of the conformation set. "
                    + "Accepted cut-edge pair energies are zeroed. "
                    + "Triple corrections are disabled.");
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse emat stats: cutPositionPairs="
                    + sparseTotalCutPositionPairs
                    + ", zeroedRcPairs=" + sparseZeroedRcPairs
                    + ", prunedRcPairs=" + sparsePrunedRcPairs);
        } else {
            branchRigidEmat = rigidEmat;
            branchMinimizingEmat = minimizingEmat;
        }

        // Step 2: Compute branch decomposition
        BranchDecomposition.Strategy firstStageDecompStrategy =
                decompStrategy == BranchDecomposition.Strategy.ADAPTIVE
                        ? BranchDecomposition.Strategy.WEIGHTED_HICKS
                        : decompStrategy;
        branchDecomposition = new BranchDecomposition(
                interactionGraph, firstStageDecompStrategy, positionStateCounts(rcs));
        branchDecomposition.compute();
        branchwidth = branchDecomposition.getBranchwidth();

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Branch decomposition computed. Branchwidth=" + branchwidth
                + ", positions=" + interactionGraph.getNumPositions()
                + ", edges=" + branchDecomposition.getTree().getNumEdges());

        // Step 3: Root the tree. Legacy uses split edge 0; optional strategies
        // can choose a split that exposes multiple pending edges for lookahead.
        RootingCandidate rooting;
        RootingCandidate firstStageRooting = null;
        int firstStageBranchwidth = branchwidth;
        boolean adaptiveAttempted = false;
        boolean adaptiveAccepted = false;
        long rootSelectionStartNanos = System.nanoTime();
        try {
            rooting = selectRooting(rcs);
            firstStageRooting = rooting;

            if (shouldRetryAdaptiveDecomposition(rooting)) {
                BranchDecomposition firstStageDecomposition = branchDecomposition;
                adaptiveAttempted = true;
                long adaptiveStartNanos = System.nanoTime();
                BranchDpAdmission.ExactPolicy exactPolicy =
                        BranchDpAdmission.getExactPolicy(admissionStateKey);
                BranchDecomposition.ExactImproveOptions exactOptions =
                        exactPolicy == null
                                ? BranchDecomposition.configuredExactImproveOptions()
                                : new BranchDecomposition.ExactImproveOptions(
                                exactPolicy.minDrop, exactPolicy.maxDrop,
                                exactPolicy.maxMillis);
                double firstStageHours = predictedSeconds(rooting) / 3600.0;
                System.out.println(BranchDpConfig.getBackendLogPrefix()
                        + " Adaptive decomposition retry: predictedHours="
                        + String.format(Locale.ROOT, "%.4f", firstStageHours)
                        + ", gpuWork=" + rooting.gpuWork
                        + ", firstStageBranchwidth=" + firstStageBranchwidth
                        + ", targetDrop=" + exactOptions.minDrop
                        + ".." + exactOptions.maxDrop
                        + ", maxMillis=" + exactOptions.maxMillis
                        + ", source=" + (exactPolicy == null
                        ? "state-soft-budget" : "whole-case-final-pass")
                        + ", strategy=WEIGHTED_EXACT_IMPROVE.");

                BranchDecomposition exactImproved = new BranchDecomposition(
                        interactionGraph,
                        BranchDecomposition.Strategy.WEIGHTED_EXACT_IMPROVE,
                        positionStateCounts(rcs), exactOptions);
                exactImproved.compute();
                branchDecomposition = exactImproved;
                branchwidth = exactImproved.getBranchwidth();
                RootingCandidate exactRooting = selectRooting(rcs);

                if (isBetterAdaptiveDecomposition(exactRooting, firstStageRooting)) {
                    rooting = exactRooting;
                    adaptiveAccepted = true;
                    System.out.println(BranchDpConfig.getBackendLogPrefix()
                            + " Adaptive decomposition accepted: branchwidth="
                            + firstStageBranchwidth + "->" + branchwidth
                            + ", gpuWork=" + firstStageRooting.gpuWork
                            + "->" + exactRooting.gpuWork
                            + ", maxDPTableBytes=" + firstStageRooting.maxDPTableBytes
                            + "->" + exactRooting.maxDPTableBytes
                            + ", elapsedMs=" + String.format(Locale.ROOT, "%.1f",
                            (System.nanoTime() - adaptiveStartNanos) / 1_000_000.0));
                } else {
                    branchDecomposition = firstStageDecomposition;
                    branchwidth = firstStageBranchwidth;
                    rooting = firstStageRooting;
                    System.out.println(BranchDpConfig.getBackendLogPrefix()
                            + " Adaptive decomposition rejected: candidate gpuWork="
                            + (exactRooting == null ? "unavailable" : exactRooting.gpuWork)
                            + ", keeping first-stage gpuWork=" + firstStageRooting.gpuWork
                            + ", elapsedMs=" + String.format(Locale.ROOT, "%.1f",
                            (System.nanoTime() - adaptiveStartNanos) / 1_000_000.0));
                }
            }

            admissionPrediction = makeAdmissionPrediction(rooting,
                    firstStageRooting, firstStageBranchwidth,
                    adaptiveAttempted, adaptiveAccepted);
            System.out.println(BranchDpConfig.getBackendLogPrefix()
                    + " admission preview: " + admissionPrediction.format()
                    + ", allocationPhase=before-materialization");
            enforcePerStateGpuWorkAdmission(rooting);
            BranchDpAdmission.enforceRetainedPredictionCeiling(
                    admissionPrediction, BranchDpConfig.getBackendLogPrefix());
            rooting = materializeSelectedRoot(rcs, rooting);
        } catch (DPTableTooLargeException e) {
            dpTooLargeException = e;
            reportDPTooLarge(e);
            return;
        }
        double rootSelectionMs = (System.nanoTime() - rootSelectionStartNanos) / 1_000_000.0;
        rootedRoot = rooting != null ? rooting.root : null;
        if (rootedRoot == null) {
            // 0-edge graph (no pairwise interactions): Z factorizes per position.
            // Default off; when enabled, report the exact independent-position DP
            // bound as a diagnostic. Building a live rootedRootEdge + sampling
            // for this degenerate case is future work (see doc 1.2), so we still
            // fall back to standard MARK* for the actual run.
            if (branchDecomposition.getTree().getNumEdges() == 0
                    && getConfigBoolean(ZERO_EDGE_DIRECT_PROPERTY, false)) {
                double rtForDP = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
                int numPos = interactionGraph.getNumPositions();
                int[] positions = new int[numPos];
                for (int i = 0; i < numPos; i++) positions[i] = i;
                double[] z = RootedTreeEdge.independentPositionLogZ(
                        branchRigidEmat, branchMinimizingEmat, rcs, positions, rtForDP);
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " 0-edge direct independent-position DP:"
                        + " logZLower=" + String.format(Locale.ROOT, "%.4f", z[0])
                        + ", logZUpper=" + String.format(Locale.ROOT, "%.4f", z[1])
                        + " (diagnostic; falling back to standard MARK* for the run)");
            }
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Empty tree, falling back to standard MARK*.");
            useBranchDecomposition = false;
            return;
        }
        selectedRootSplitIndex = rooting.splitEdgeIndex;
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Root split edge=" + selectedRootSplitIndex
                + " (lambdaEdges=" + rooting.lambdaEdges
                + ", rootFset=" + rooting.rootFsetSize
                + ", maxFset=" + rooting.maxFsetSize
                + ", branchingEdges=" + rooting.branchingEdges
                + ", totalFsetEdges=" + rooting.totalFsetEdges
                + ", reusableLambdaEdges=" + rooting.reusableLambdaEdges
                + ", reusableEdgeRatio=" + String.format("%.4f", rooting.reusableEdgeRatio())
                + ", reusableWorkRatio=" + String.format("%.4f", rooting.reusableWorkRatio())
                + ", maxMStates=" + rooting.maxMStates
                + ", totalMStates=" + rooting.totalMStates
                + ", maxDPTable=" + formatBytes(rooting.maxDPTableBytes)
                + ", totalDPTable=" + formatBytes(rooting.totalDPTableBytes)
                + ", maxHostBytes=" + formatBytes(rooting.maxHostBytes)
                + ", totalHostBytes=" + formatBytes(rooting.totalHostBytes)
                + ", maxFileBackedBytes=" + formatBytes(rooting.maxFileBackedBytes)
                + ", totalFileBackedBytes=" + formatBytes(rooting.totalFileBackedBytes)
                + ", hostBudget=" + formatBytes(selectedRootHostBudgetBytes)
                + ", logDPWork=" + String.format(Locale.ROOT, "%.2f", rooting.logDPWork)
                + ", gpuWork=" + rooting.gpuWork
                + ", gpuUnsupportedEdges=" + rooting.gpuUnsupportedEdges
                + ", maxFullDeviceBytes=" + formatBytes(rooting.maxSingleGpuBytes)
                + ", fitsSingleGpu=" + rooting.fitsSingleGpu
                + ", estimatedSlicedTrafficBytes="
                + formatSlicedTraffic(rooting.estimatedSlicedTraffic)
                + ", rootSelectionMs=" + String.format(Locale.ROOT, "%.3f", rootSelectionMs)
                + ")");

        rootedRootEdge = rootedRoot.getLeftChild().getChildOfEdge();

        // Step 4: TESS check
        double logTESS = rootedRootEdge.computeLogTESS();
        double logNaive = computeLogNaive(rcs);
        double tessRatio = Math.exp(logTESS - logNaive);

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Compact tree TESS: logTESS=" + String.format("%.2f", logTESS)
                + ", logNaive=" + String.format("%.2f", logNaive)
                + ", ratio=" + String.format("%.4f", tessRatio));
        long worstCaseSingleGpuBytes = logDPMemoryPredictions("initial", logTESS, rcs);
        checkWorstCaseSingleGpuFits(worstCaseSingleGpuBytes);

        if (dryRun) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DRYRUN summary"
                    + " rootSplit=" + rootSplitStrategy
                    + " splitEdge=" + selectedRootSplitIndex
                    + " branchwidth=" + branchwidth
                    + " lambdaEdges=" + rooting.lambdaEdges
                    + " maxMStates=" + rooting.maxMStates
                    + " totalMStates=" + rooting.totalMStates
                    + " maxDPTableBytes=" + rooting.maxDPTableBytes
                    + " totalDPTableBytes=" + rooting.totalDPTableBytes
                    + " maxHostBytes=" + rooting.maxHostBytes
                    + " totalHostBytes=" + rooting.totalHostBytes
                    + " maxFileBackedBytes=" + rooting.maxFileBackedBytes
                    + " totalFileBackedBytes=" + rooting.totalFileBackedBytes
                    + " hostBudgetBytes=" + selectedRootHostBudgetBytes
                    + " logDPWork=" + String.format(Locale.ROOT, "%.4f", rooting.logDPWork)
                    + " gpuWork=" + rooting.gpuWork
                    + " gpuUnsupportedEdges=" + rooting.gpuUnsupportedEdges
                    + " logTESS=" + String.format(Locale.ROOT, "%.4f", logTESS)
                    + " logNaive=" + String.format(Locale.ROOT, "%.4f", logNaive)
                    + " maxFullDeviceBytes=" + rooting.maxSingleGpuBytes
                    + " fitsSingleGpu=" + rooting.fitsSingleGpu
                    + " estimatedSlicedTrafficBytes="
                    + formatSlicedTraffic(rooting.estimatedSlicedTraffic)
                    + " rootSelectionMs="
                    + String.format(Locale.ROOT, "%.3f", rootSelectionMs)
                    + " worstCaseSingleGpuBytes=" + worstCaseSingleGpuBytes);
            return;
        }

        if (tessRatio > TESS_FALLBACK_THRESHOLD) {
            if (useSparsePfunc) {
                // SPARSE mode uses branch* emat (pruned + zeroed); standard MARK* would
                // sum over the FULL conformation set with the original emat, which is a
                // different pfunc target. Stay on branch decomposition to keep Z = Z_sparse.
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " TESS ratio " + String.format("%.4f", tessRatio)
                        + " > threshold " + TESS_FALLBACK_THRESHOLD
                        + ", continuing on branch decomposition to keep SPARSE pfunc self-consistent.");
            } else {
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " TESS ratio " + String.format("%.4f", tessRatio)
                        + " > threshold " + TESS_FALLBACK_THRESHOLD + ", falling back to standard MARK*.");
                useBranchDecomposition = false;
                return;
            }
        }

        // Initialize scorers and root search node
        initSearch(rcs);
    }

    @Override
    public void setCorrections(UpdatingEnergyMatrix cachedCorrections) {
        super.setCorrections(cachedCorrections);
        if (useHigherOrderCorrections || correctionAudit) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " correction matrix attached"
                    + " (tripleCorrections=" + useHigherOrderCorrections
                    + ", correctionAudit=" + correctionAudit
                    + ", energyMode=" + energyMode + ")");
        }
    }

    private InteractionGraph buildInteractionGraph(RCs rcs) {
        switch (cutoffStrategy) {
            case COMPLETE:
                System.out.println("InteractionGraph complete: no interaction edges cut.");
                return InteractionGraph.buildComplete(rcs.getNumPos());
            case DUAL:
                return InteractionGraph.buildWithDualCutoff(
                        confSpace, rigidEmat, minimizingEmat, rcs,
                        distCutoff, energyCutoff);
            case RESIDUAL_BUDGET:
            default:
                return InteractionGraph.buildWithResidualBudget(
                        confSpace, rigidEmat, minimizingEmat, rcs,
                        residualBudget, cutoffKeepConnected);
        }
    }

    private void printSparseResidualInfo() {
        double rho = interactionGraph.getCutResidualUpperBound();
        if (rho <= 0.0) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse residual rho=0; branch graph matches the full pair graph.");
            return;
        }

        if (!Double.isFinite(rho)) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse residual rho is infinite; no finite sparse-to-full Z certificate.");
            return;
        }

        double rt = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
        double log10Factor = rho / rt / Math.log(10.0);
        // Note: this residual bound assumes we sum over the FULL conformation set
        // with the sparse energy (FULL mode). Under SPARSE mode we additionally
        // remove cut-edge confs from C, so this is reported for context only.
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse residual certificate: cutEdges="
                + interactionGraph.getNumCutEdges()
                + ", rho<=" + String.format("%.6f", rho) + " kcal/mol"
                + ", per-state Z factor<=10^" + String.format("%.3f", log10Factor)
                + " (exp(rho/RT)); full-bound inflation="
                + (certifyFullBounds ? "enabled" : "disabled (forced off in SPARSE mode)"));
    }

    private void printExactTreewidthDiagnostic() {
        if (!exactTreewidthDiagnostic) {
            return;
        }

        int numPositions = interactionGraph.getNumPositions();
        if (numPositions > exactTreewidthMaxPositions) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Exact treewidth diagnostic skipped: positions="
                    + numPositions + " > maxPositions=" + exactTreewidthMaxPositions + ".");
            return;
        }

        ExactTreewidth.Result result = ExactTreewidth.compute(interactionGraph);
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Exact treewidth diagnostic:"
                + " tw=" + result.treewidth
                + ", branchwidthLowerBound=" + result.branchwidthLowerBound
                + ", branchwidthUpperBound=" + result.branchwidthUpperBound
                + ", positions=" + result.numPositions
                + ", graphEdges=" + result.numEdges
                + ", density=" + String.format(Locale.ROOT, "%.3f", interactionGraph.getDensity())
                + ", elapsedMs=" + String.format(Locale.ROOT, "%.2f", result.elapsedMillis()));
    }

    /**
     * Build a sparse copy of {@code source} for use as branchRigidEmat or
     * branchMinimizingEmat under SPARSE mode.
     *
     * For each cut edge (pos1, pos2):
     *   - per (rc1, rc2), compare max(|rigid|, |minimizing|) against {@code pruneThreshold}.
     *   - if it exceeds the threshold, set sparse[pos1, rc1, pos2, rc2] to +inf
     *     so the full conformation containing this (rc1, rc2) is removed from
     *     C_sparse (A* skips it, leaf minimization never reaches it).
     *   - otherwise set sparse[pos1, rc1, pos2, rc2] to 0 (cut-edge pair energy
     *     is dropped from E_sparse for accepted confs).
     *
     * Pruning uses max(|rigid|, |minimizing|) so both bounds stay consistent —
     * we never accept a conf in the lower-bound emat that the upper-bound emat
     * would have rejected.
     *
     * The first call (rigid) populates the shared counters
     * sparseTotalCutPositionPairs / sparseZeroedRcPairs / sparsePrunedRcPairs;
     * the second call (minimizing) overwrites them but the totals are identical
     * by construction since the pruning decision is symmetric.
     */
    private EnergyMatrix makeSparseEnergyMatrix(String label, EnergyMatrix source,
                                                InteractionGraph graph,
                                                double pruneThreshold, boolean recordStats) {
        EnergyMatrix sparse = new EnergyMatrix(source);
        int cutPositionPairs = 0;
        int zeroedRcPairs = 0;
        int prunedRcPairs = 0;

        for (int pos1 = 0; pos1 < sparse.getNumPos(); pos1++) {
            for (int pos2 = pos1 + 1; pos2 < sparse.getNumPos(); pos2++) {
                if (graph.hasEdge(pos1, pos2)) {
                    continue;
                }
                cutPositionPairs++;
                int n1 = sparse.getNumConfAtPos(pos1);
                int n2 = sparse.getNumConfAtPos(pos2);
                for (int rc1 = 0; rc1 < n1; rc1++) {
                    for (int rc2 = 0; rc2 < n2; rc2++) {
                        double rigidE = rigidEmat.getPairwise(pos1, rc1, pos2, rc2);
                        double minE = minimizingEmat.getPairwise(pos1, rc1, pos2, rc2);
                        double worstAbs;
                        if (Double.isNaN(rigidE) || Double.isNaN(minE)) {
                            worstAbs = Double.POSITIVE_INFINITY;
                        } else {
                            worstAbs = Math.max(
                                    Double.isInfinite(rigidE) ? Double.POSITIVE_INFINITY : Math.abs(rigidE),
                                    Double.isInfinite(minE) ? Double.POSITIVE_INFINITY : Math.abs(minE));
                        }
                        if (worstAbs > pruneThreshold) {
                            // Prune this (rc1, rc2): mark unreachable in both emats.
                            sparse.setPairwise(pos1, rc1, pos2, rc2, Double.POSITIVE_INFINITY);
                            prunedRcPairs++;
                        } else {
                            sparse.setPairwise(pos1, rc1, pos2, rc2, 0.0);
                            zeroedRcPairs++;
                        }
                    }
                }
            }
        }

        if (recordStats) {
            this.sparseTotalCutPositionPairs = cutPositionPairs;
            this.sparseZeroedRcPairs = zeroedRcPairs;
            this.sparsePrunedRcPairs = prunedRcPairs;
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse " + label + " emat: cutPositionPairs="
                + cutPositionPairs + ", zeroedRcPairs=" + zeroedRcPairs
                + ", prunedRcPairs=" + prunedRcPairs
                + " (pruneThreshold=" + String.format("%.6f", pruneThreshold) + " kcal/mol)");
        return sparse;
    }

    /**
     * Residue interactions used for leaf minimization under SPARSE mode:
     * singles + only the pair interactions whose (pos1, pos2) is in the
     * interaction graph + shells (if requested by the ecalc).
     *
     * This makes the minimized energy E_true ALSO live in the sparse energy
     * world — minimized leaves are scored against the same pair-set used by
     * the pfunc bounds, keeping Z bounds tight and Z self-consistent.
     */
    private ResidueInteractions makeSparseFullConfInters(int[] conf) {
        ResidueInteractions inters = new ResidueInteractions();

        for (int pos = 0; pos < conf.length; pos++) {
            if (conf[pos] < 0) continue;
            inters.addAll(minimizingEcalc.makeSingleInters(pos, conf[pos]));
        }

        for (int pos1 = 0; pos1 < conf.length; pos1++) {
            if (conf[pos1] < 0) continue;
            for (int pos2 = pos1 + 1; pos2 < conf.length; pos2++) {
                if (conf[pos2] < 0 || !interactionGraph.hasEdge(pos1, pos2)) continue;
                inters.addAll(minimizingEcalc.makePairInters(pos1, conf[pos1], pos2, conf[pos2]));
            }
        }

        if (minimizingEcalc.addShellInters) {
            inters.addAll(minimizingEcalc.makeShellInters());
        }

        return inters;
    }

    protected EdgeSelectionStrategy getEdgeSelectionStrategy() {
        return EdgeSelectionStrategy.fromProperty(getConfigProperty(EDGE_SELECTION_PROPERTY, null));
    }

    protected int getEdgeLookaheadMaxStates() {
        return Math.max(1, getConfigInteger(EDGE_LOOKAHEAD_MAX_STATES_PROPERTY, 2000));
    }

    protected int getEdgeLookaheadMaxPendingEdges() {
        return Math.max(2, getConfigInteger(EDGE_LOOKAHEAD_MAX_PENDING_EDGES_PROPERTY, 4));
    }

    protected boolean getEdgeLookaheadParallel() {
        return getConfigBoolean(EDGE_LOOKAHEAD_PARALLEL_PROPERTY, true);
    }

    protected boolean getUseHigherOrderCorrections() {
        return getConfigBoolean("branchdp.useHigherOrderCorrections", false);
    }

    private static int[] positionStateCounts(RCs rcs) {
        int[] counts = new int[rcs.getNumPos()];
        for (int pos = 0; pos < counts.length; pos++) {
            counts[pos] = Math.max(1, rcs.getNum(pos));
        }
        return counts;
    }

    protected void logBackendControlOverrides() {
    }

    protected ConfEnergyCalculator getMinimizingEcalc() {
        return minimizingEcalc;
    }

    protected String getConfigProperty(String key, String defaultValue) {
        return BranchDpConfig.getBackendProperty(key, defaultValue);
    }

    protected int getConfigInteger(String key, int defaultValue) {
        return BranchDpConfig.getBackendInteger(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    protected long getConfigLong(String key, long defaultValue) {
        return BranchDpConfig.getBackendLong(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    protected long getConfigBytes(String key, long defaultValue) {
        return BranchDpConfig.getBackendBytes(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    protected double getConfigDouble(String key, double defaultValue) {
        return BranchDpConfig.getBackendDouble(key, defaultValue, BranchDpConfig.getBackendLogPrefix());
    }

    protected boolean getConfigBoolean(String key, boolean defaultValue) {
        return BranchDpConfig.getBackendBoolean(key, defaultValue);
    }

    /** Maximum full-tree DP sweeps used for conservative state admission. */
    protected int getAdmissionDpSweeps() {
        return Math.max(1, getConfigInteger(
                BranchDpAdmission.DP_SWEEPS_PROPERTY, 1));
    }

    private RootingCandidate selectRooting(RCs rcs) {
        int numEdges = branchDecomposition.getTree().getNumEdges();
        if (numEdges == 0) return null;

        selectedRootHostBudgetBytes = resolveRootHostBudgetBytes();

        String strategy = rootSplitStrategy.toLowerCase(Locale.ROOT);
        if (strategy.isEmpty() || strategy.equals("auto")) {
            strategy = "memory";
        }
        if (strategy.equals("legacy") || strategy.equals("edge0")) {
            selectedRootGpuBudgetBytes = Long.MAX_VALUE;
            return evaluateRootSplit(rcs, 0, false);
        }

        try {
            int explicitSplit = Integer.parseInt(strategy);
            selectedRootGpuBudgetBytes = Long.MAX_VALUE;
            return evaluateRootSplit(rcs, explicitSplit, false);
        } catch (NumberFormatException ignored) {
            // fall through to named strategies
        }

        boolean useReuseScoring = strategy.equals("reuse");
        boolean useMemoryScoring = strategy.equals("memory")
                || strategy.equals("mem")
                || strategy.equals("dp")
                || strategy.equals("dpmemory")
                || strategy.equals("dp_memory");
        boolean useWorkScoring = strategy.equals("work")
                || strategy.equals("dpwork")
                || strategy.equals("dp_work");
        // GPU-aware exhaustive root scoring: prefer structurally supported roots that
        // fit the single-device budget, minimize exact DP work when resident, and use
        // estimated child-slice traffic before work when every candidate must slice.
        boolean useGpuBytesScoring = strategy.equals("gpubytes")
                || strategy.equals("gpu")
                || strategy.equals("devicebytes")
                || strategy.equals("vram");
        boolean usePredictedScoring = strategy.equals("predicted")
                || strategy.equals("predictedhours")
                || strategy.equals("predicted_hours")
                || strategy.equals("admission")
                || strategy.equals("sla");

        if (!strategy.equals("branching") && !strategy.equals("maxfset")
                && !strategy.equals("lookahead") && !useReuseScoring
                && !useMemoryScoring && !useWorkScoring && !useGpuBytesScoring
                && !usePredictedScoring) {
            System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown root split strategy '" + rootSplitStrategy
                    + "', using legacy split edge 0.");
            selectedRootGpuBudgetBytes = Long.MAX_VALUE;
            return evaluateRootSplit(rcs, 0, false);
        }

        double logNaive = computeLogNaive(rcs);
        Set<Integer> mutablePositions = identifyMutablePositions(rcs);
        if (strategy.equals("reuse")) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " rootSplit=reuse mutablePositions="
                    + formatPositionsWithResidues(mutablePositions));
            if (mutablePositions.isEmpty()) {
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " rootSplit=reuse found no mutable positions; "
                        + "falling back to TESS/fset root scoring.");
            }
        }
        useReuseScoring = useReuseScoring && !mutablePositions.isEmpty();
        long gpuBudgetBytes = (useGpuBytesScoring || usePredictedScoring)
                ? resolveRootGpuBudgetBytes() : Long.MAX_VALUE;
        selectedRootGpuBudgetBytes = gpuBudgetBytes;
        RootingCandidate best = null;
        RootingCandidate leastHost = null;
        for (int splitIdx = 0; splitIdx < numEdges; splitIdx++) {
            RootingCandidate candidate = evaluateRootSplit(rcs, splitIdx, false, gpuBudgetBytes);
            if (candidate == null) continue;
            if (leastHost == null || candidate.totalHostBytes < leastHost.totalHostBytes) {
                leastHost = candidate;
            }
            boolean fitsHostBudget = candidate.totalHostBytes <= selectedRootHostBudgetBytes;
            if (useGpuBytesScoring || usePredictedScoring) {
                System.out.println(BranchDpConfig.getBackendLogPrefix()
                        + " rootSplit=" + (usePredictedScoring ? "predicted" : "gpubytes")
                        + " candidate=" + splitIdx
                        + ", maxHostBytes=" + candidate.maxHostBytes
                        + ", totalHostBytes=" + candidate.totalHostBytes
                        + ", maxFileBackedBytes=" + candidate.maxFileBackedBytes
                        + ", totalFileBackedBytes=" + candidate.totalFileBackedBytes
                        + ", fitsHostBudget=" + fitsHostBudget
                        + ", fullDeviceBytes=" + candidate.maxSingleGpuBytes
                        + ", fitsSingleGpu=" + candidate.fitsSingleGpu
                        + ", estimatedSlicedTrafficBytes="
                        + formatSlicedTraffic(candidate.estimatedSlicedTraffic)
                        + ", logDPWork=" + String.format(Locale.ROOT, "%.4f", candidate.logDPWork)
                        + ", gpuWork=" + candidate.gpuWork
                        + ", predictedHours=" + String.format(Locale.ROOT, "%.4f",
                        predictedSeconds(candidate) / 3600.0)
                        + ", gpuUnsupportedEdges=" + candidate.gpuUnsupportedEdges);
            }
            if (!fitsHostBudget) {
                continue;
            }
            boolean better;
            if (usePredictedScoring) {
                better = isBetterPredictedRooting(candidate, best, logNaive);
            } else if (useGpuBytesScoring) {
                better = isBetterGpuBytesRooting(candidate, best, logNaive);
            } else if (useMemoryScoring) {
                better = isBetterMemoryRooting(candidate, best, logNaive);
            } else if (useWorkScoring) {
                better = isBetterWorkRooting(candidate, best, logNaive);
            } else if (useReuseScoring) {
                better = isBetterReuseRooting(candidate, best, logNaive);
            } else {
                if (candidate.maxFsetSize > rootSplitMaxFset) continue;
                better = isBetterRooting(candidate, best, logNaive);
            }
            if (better) {
                best = candidate;
            }
        }
        if (best == null && leastHost != null) {
            throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix()
                    + " no root split fits the Java-heap budget "
                    + formatBytes(selectedRootHostBudgetBytes)
                    + "; minimum estimated live host storage is "
                    + formatBytes(leastHost.totalHostBytes)
                    + " at split edge " + leastHost.splitEdgeIndex
                    + ". Raise -D" + ROOT_SPLIT_HOST_BUDGET_PROPERTY
                    + " only if the node and -Xmx have sufficient headroom.");
        }
        if (best == null) {
            return evaluateRootSplit(rcs, 0, false, gpuBudgetBytes);
        }
        return best;
    }

    private RootingCandidate materializeSelectedRoot(RCs rcs, RootingCandidate preview) {
        if (preview == null) {
            return null;
        }
        if (preview.totalHostBytes > selectedRootHostBudgetBytes) {
            throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix()
                    + " requested root split edge " + preview.splitEdgeIndex
                    + " needs an estimated " + formatBytes(preview.totalHostBytes)
                    + " of live Java-heap storage, exceeding root host budget "
                    + formatBytes(selectedRootHostBudgetBytes) + ".");
        }
        return dryRun
                ? preview
                : evaluateRootSplit(rcs, preview.splitEdgeIndex, true,
                selectedRootGpuBudgetBytes);
    }

    private boolean shouldRetryAdaptiveDecomposition(RootingCandidate rooting) {
        if (decompStrategy != BranchDecomposition.Strategy.ADAPTIVE
                || rooting == null) {
            return false;
        }
        if (BranchDpAdmission.getExactPolicy(admissionStateKey) != null) {
            return true;
        }

        double softHours = Math.max(0.0, getConfigDouble(
                BranchDpAdmission.SOFT_STATE_HOURS_PROPERTY, 0.0));
        if (softHours > 0.0) {
            double seconds = predictedSeconds(rooting);
            // A missing rate or unavailable traffic is deliberately treated as
            // over budget. Whole-case admission will later reject it unless the
            // exact retry produces a finite, calibrated estimate.
            return !Double.isFinite(seconds) || seconds > softHours * 3600.0;
        }

        // Backward-compatible fallback for callers that have not configured a
        // hardware time model. Production case admission sets softStateHours.
        return rooting.gpuWork.compareTo(adaptiveGpuWorkThreshold()) > 0;
    }

    private BigInteger adaptiveGpuWorkThreshold() {
        return getConfigBigInteger(ADAPTIVE_GPU_WORK_THRESHOLD_PROPERTY,
                DEFAULT_ADAPTIVE_GPU_WORK_THRESHOLD);
    }

    private BigInteger getConfigBigInteger(String key, BigInteger defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            BigInteger parsed = new BigInteger(value.trim().replace("_", ""));
            return parsed.signum() < 0 ? BigInteger.ZERO : parsed;
        } catch (NumberFormatException e) {
            System.err.println(BranchDpConfig.getBackendLogPrefix()
                    + " Invalid integer for '" + key + "': '" + value
                    + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private boolean isBetterAdaptiveDecomposition(RootingCandidate candidate,
                                                   RootingCandidate firstStage) {
        if (candidate == null) return false;
        if (firstStage == null) return true;

        int unsupportedCmp = Integer.compare(candidate.gpuUnsupportedEdges,
                firstStage.gpuUnsupportedEdges);
        if (unsupportedCmp != 0) return unsupportedCmp < 0;

        int predictedCmp = Double.compare(predictedSeconds(candidate),
                predictedSeconds(firstStage));
        if (predictedCmp != 0) return predictedCmp < 0;

        int workCmp = candidate.gpuWork.compareTo(firstStage.gpuWork);
        if (workCmp != 0) return workCmp < 0;

        int trafficCmp = candidate.estimatedSlicedTraffic.compareTo(
                firstStage.estimatedSlicedTraffic);
        if (trafficCmp != 0) return trafficCmp < 0;

        int maxTableCmp = Long.compare(candidate.maxDPTableBytes,
                firstStage.maxDPTableBytes);
        if (maxTableCmp != 0) return maxTableCmp < 0;

        return candidate.totalDPTableBytes < firstStage.totalDPTableBytes;
    }

    private double predictedSeconds(RootingCandidate rooting) {
        if (rooting == null) return 0.0;
        boolean trafficAvailable = rooting.estimatedSlicedTraffic
                .compareTo(SLICED_TRAFFIC_UNAVAILABLE) < 0;
        BigInteger sweeps = BigInteger.valueOf(admissionDpSweeps);
        return admissionHardware.totalSeconds(
                rooting.gpuWork.multiply(sweeps),
                trafficAvailable
                        ? rooting.estimatedSlicedTraffic.multiply(sweeps)
                        : BigInteger.ZERO,
                trafficAvailable, rooting.gpuUnsupportedEdges);
    }

    private BranchDpAdmission.Prediction makeAdmissionPrediction(
            RootingCandidate rooting,
            RootingCandidate firstStageRooting,
            int firstStageBranchwidth,
            boolean adaptiveAttempted,
            boolean adaptiveAccepted) {
        BigInteger gpuWork = rooting == null ? BigInteger.ZERO : rooting.gpuWork;
        BigInteger traffic = rooting == null
                ? BigInteger.ZERO : rooting.estimatedSlicedTraffic;
        boolean trafficAvailable = rooting == null
                || traffic.compareTo(SLICED_TRAFFIC_UNAVAILABLE) < 0;
        if (!trafficAvailable) traffic = BigInteger.ZERO;
        BigInteger firstWork = firstStageRooting == null
                ? gpuWork : firstStageRooting.gpuWork;
        return new BranchDpAdmission.Prediction(
                stateName, admissionStateKey, branchwidth,
                rooting == null ? -1 : rooting.splitEdgeIndex,
                gpuWork, traffic, trafficAvailable,
                rooting == null ? 0 : rooting.gpuUnsupportedEdges,
                admissionHardware, admissionDpSweeps,
                firstStageBranchwidth, firstWork,
                predictedSeconds(firstStageRooting),
                adaptiveAttempted, adaptiveAccepted);
    }

    /** Allocation-free estimate exported to PACK* whole-case preflight. */
    public BranchDpAdmission.Prediction getAdmissionPrediction() {
        return admissionPrediction;
    }

    private void enforcePerStateGpuWorkAdmission(RootingCandidate rooting) {
        if (rooting == null) return;
        BigInteger maxGpuWork = getConfigBigInteger(
                ADMISSION_MAX_GPU_WORK_PROPERTY, BigInteger.ZERO);
        if (maxGpuWork.signum() <= 0) return;
        if (rooting.gpuWork.compareTo(maxGpuWork) > 0) {
            throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix()
                    + " GPU-work admission rejected state=" + stateName
                    + ": predicted gpuWork=" + rooting.gpuWork
                    + " exceeds configured maximum=" + maxGpuWork
                    + ". This check ran before DP-table materialization. Raise -D"
                    + ADMISSION_MAX_GPU_WORK_PROPERTY
                    + " or set it to 0 only for an intentional override.");
        }
    }

    private long resolveRootHostBudgetBytes() {
        long configured = getConfigBytes(ROOT_SPLIT_HOST_BUDGET_PROPERTY, 0L);
        long maxHeap = Runtime.getRuntime().maxMemory();
        long budget = configured > 0L
                ? configured
                : Math.max(1L, (long) Math.floor(maxHeap * DEFAULT_ROOT_SPLIT_HOST_HEAP_FRACTION));
        System.out.println(BranchDpConfig.getBackendLogPrefix()
                + " rootSplit Java-heap budget=" + formatBytes(budget)
                + (configured > 0L
                ? " (configured)"
                : " (auto=" + String.format(Locale.ROOT, "%.0f%%", 100.0 * DEFAULT_ROOT_SPLIT_HOST_HEAP_FRACTION)
                    + " of maxHeap " + formatBytes(maxHeap) + ")"));
        return budget;
    }

    private long resolveRootGpuBudgetBytes() {
        long configured = getConfigBytes(ROOT_SPLIT_GPU_BUDGET_PROPERTY, 0L);
        long budget = configured > 0L ? configured : DPGpuFullDP.queryMinUsableVramBytes();
        if (budget < 0L) {
            budget = 0L;
            System.out.println(BranchDpConfig.getBackendLogPrefix()
                    + " rootSplit=gpubytes could not query GPU VRAM; treating every candidate as sliced. "
                    + "Set -D" + ROOT_SPLIT_GPU_BUDGET_PROPERTY + "=<bytes> for offline/dry-run scoring.");
        } else {
            System.out.println(BranchDpConfig.getBackendLogPrefix()
                    + " rootSplit=gpubytes single-GPU budget=" + formatBytes(budget)
                    + (configured > 0L ? " (configured)" : " (queried)"));
        }
        return budget;
    }

    private RootingCandidate evaluateRootSplit(RCs rcs, int splitEdgeIndex,
                                               boolean initEnumerationArrays) {
        return evaluateRootSplit(rcs, splitEdgeIndex, initEnumerationArrays, Long.MAX_VALUE);
    }

    private RootingCandidate evaluateRootSplit(RCs rcs, int splitEdgeIndex,
                                               boolean initEnumerationArrays,
                                               long gpuBudgetBytes) {
        RootedTreeNode root = branchDecomposition.rootBranchTree(rcs, splitEdgeIndex);
        if (root == null) return null;

        try {
            RootedTreeEdge.postOrderCompLlambda(root, initEnumerationArrays);
        } catch (IllegalStateException e) {
            if (!initEnumerationArrays) {
                System.err.println(BranchDpConfig.getBackendLogPrefix() + " Skipping root split edge " + splitEdgeIndex
                        + " during root selection: " + e.getMessage());
                return null;
            }
            throw e;
        }
        RootedTreeEdge rootEdge = root.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(root, lambdaEdges);

        int maxFsetSize = 0;
        int branchingEdges = 0;
        int totalFsetEdges = 0;
        long maxMStates = 0L;
        long totalMStates = 0L;
        long maxDPTableBytes = 0L;
        long totalDPTableBytes = 0L;
        long maxHostBytes = 0L;
        long totalHostBytes = 0L;
        long maxFileBackedBytes = 0L;
        long totalFileBackedBytes = 0L;
        double logDPWork = Double.NEGATIVE_INFINITY;
        BigInteger gpuWork = BigInteger.ZERO;
        int gpuUnsupportedEdges = 0;
        long maxSingleGpuBytes = 0L;
        BigInteger estimatedSlicedTraffic = BigInteger.ZERO;
        for (RootedTreeEdge edge : lambdaEdges) {
            int fsetSize = edge.getFset() == null ? 0 : edge.getFset().size();
            maxFsetSize = Math.max(maxFsetSize, fsetSize);
            if (fsetSize > 1) {
                branchingEdges++;
            }
            totalFsetEdges += fsetSize;

            long mStates = Math.max(1L, edge.getMStateCount());
            long lambdaStates = Math.max(1, edge.getTotalLambdaStates());
            long tableBytes = estimateDPTableBytes(edge.getMStateCount());
            boolean fileBacked = RootedTreeEdge.shouldUseFileBackedDPTable(edge.getMStateCount());
            long hostBytes = fileBacked
                    ? estimateDPAuxHostBytes(edge.getMStateCount(), edge.getTotalLambdaStates())
                    : estimateDPHostBytes(edge.getMStateCount(), edge.getTotalLambdaStates());
            maxMStates = Math.max(maxMStates, mStates);
            totalMStates = saturatingAdd(totalMStates, mStates);
            maxDPTableBytes = Math.max(maxDPTableBytes, tableBytes);
            totalDPTableBytes = saturatingAdd(totalDPTableBytes, tableBytes);
            maxHostBytes = Math.max(maxHostBytes, hostBytes);
            totalHostBytes = saturatingAdd(totalHostBytes, hostBytes);
            if (fileBacked) {
                maxFileBackedBytes = Math.max(maxFileBackedBytes, tableBytes);
                totalFileBackedBytes = saturatingAdd(totalFileBackedBytes, tableBytes);
            }
            logDPWork = RootedTreeEdge.logSumExp(logDPWork,
                    Math.log((double) mStates) + Math.log((double) lambdaStates));
            gpuWork = gpuWork.add(BigInteger.valueOf(mStates)
                    .multiply(BigInteger.valueOf(lambdaStates)));
            GpuEdgeShape gpuShape = estimateGpuEdgeShape(edge, rcs);
            if (!gpuShape.supported) {
                gpuUnsupportedEdges++;
            }
            maxSingleGpuBytes = Math.max(maxSingleGpuBytes, gpuShape.fullDeviceBytes);
            if (gpuShape.fullDeviceBytes > gpuBudgetBytes) {
                estimatedSlicedTraffic = estimatedSlicedTraffic.add(
                        gpuShape.estimateSlicedTraffic(gpuBudgetBytes));
            }
        }

        int rootFsetSize = rootEdge.getFset() == null ? 0 : rootEdge.getFset().size();
        Set<Integer> mutablePositions = identifyMutablePositions(rcs);
        ReuseStats reuseStats = computeReuseStats(lambdaEdges, mutablePositions, rcs);
        return new RootingCandidate(root, splitEdgeIndex, rootEdge.computeLogTESS(),
                lambdaEdges.size(), maxFsetSize, branchingEdges, totalFsetEdges, rootFsetSize,
                reuseStats.reusableEdges, reuseStats.reusableLogWork, reuseStats.totalLogWork,
                maxMStates, totalMStates, maxDPTableBytes, totalDPTableBytes,
                maxHostBytes, totalHostBytes, maxFileBackedBytes, totalFileBackedBytes,
                logDPWork,
                gpuWork, gpuUnsupportedEdges, maxSingleGpuBytes,
                maxSingleGpuBytes <= gpuBudgetBytes, estimatedSlicedTraffic);
    }

    private static final class GpuEdgeShape {
        final int mPositionCount;
        final int lambdaPositionCount;
        final long totalLambdaStates;
        final long lmTermCount;
        final int lmPairCount;
        final long childMTermTotal;
        final long childLTermTotal;
        final int numChildren;
        final long childTableTotalStates;
        final long outputTileMStates;
        final long fixedSlicedDeviceBytes;
        final long unionStateCount;
        final long[] childRowStates;
        final long[] childLambdaStates;
        final long childSliceMaxBytes;
        final DPGpuOutOfCore.PlanningInput outOfCoreInput;
        final boolean childFoldMapped;
        final boolean supported;
        final long fullDeviceBytes;

        GpuEdgeShape(int mPositionCount, int lambdaPositionCount, long totalLambdaStates,
                     long lmTermCount, int lmPairCount, long childMTermTotal,
                     long childLTermTotal, int numChildren, long childTableTotalStates,
                     long outputTileMStates, long fixedSlicedDeviceBytes,
                     long unionStateCount, long[] childRowStates, long[] childLambdaStates,
                     long childSliceMaxBytes,
                     DPGpuOutOfCore.PlanningInput outOfCoreInput,
                     boolean childFoldMapped, boolean supported, long fullDeviceBytes) {
            this.mPositionCount = mPositionCount;
            this.lambdaPositionCount = lambdaPositionCount;
            this.totalLambdaStates = totalLambdaStates;
            this.lmTermCount = lmTermCount;
            this.lmPairCount = lmPairCount;
            this.childMTermTotal = childMTermTotal;
            this.childLTermTotal = childLTermTotal;
            this.numChildren = numChildren;
            this.childTableTotalStates = childTableTotalStates;
            this.outputTileMStates = outputTileMStates;
            this.fixedSlicedDeviceBytes = fixedSlicedDeviceBytes;
            this.unionStateCount = unionStateCount;
            this.childRowStates = childRowStates;
            this.childLambdaStates = childLambdaStates;
            this.childSliceMaxBytes = childSliceMaxBytes;
            this.outOfCoreInput = outOfCoreInput;
            this.childFoldMapped = childFoldMapped;
            this.supported = supported;
            this.fullDeviceBytes = fullDeviceBytes;
        }

        BigInteger estimateSlicedTraffic(long gpuBudgetBytes) {
            if (!childFoldMapped) {
                return SLICED_TRAFFIC_UNAVAILABLE;
            }
            BigInteger best = SLICED_TRAFFIC_UNAVAILABLE;

            if (gpuBudgetBytes > fixedSlicedDeviceBytes) {
                long available = gpuBudgetBytes - fixedSlicedDeviceBytes;
                long sliceBytes = Math.min(childSliceMaxBytes, available);
                long lambdaStatesPerUnion = 0L;
                for (long states : childLambdaStates) {
                    lambdaStatesPerUnion = saturatingAdd(lambdaStatesPerUnion, states);
                }
                long bytesPerUnion = saturatingMultiply(2L * Double.BYTES,
                        Math.max(1L, lambdaStatesPerUnion));
                if (bytesPerUnion <= sliceBytes) {
                    long unionStatesPerSlice = Math.max(1L, sliceBytes / bytesPerUnion);
                    unionStatesPerSlice = Math.min(unionStatesPerSlice, unionStateCount);
                    long fullSlices = unionStateCount / unionStatesPerSlice;
                    long remainder = unionStateCount % unionStatesPerSlice;

                    BigInteger packedStates = BigInteger.ZERO;
                    for (int c = 0; c < childRowStates.length; c++) {
                        long rowsPerFullSlice = Math.min(childRowStates[c], unionStatesPerSlice);
                        BigInteger rows = BigInteger.valueOf(fullSlices)
                                .multiply(BigInteger.valueOf(rowsPerFullSlice));
                        if (remainder > 0L) {
                            rows = rows.add(BigInteger.valueOf(
                                    Math.min(childRowStates[c], remainder)));
                        }
                        packedStates = packedStates.add(rows.multiply(
                                BigInteger.valueOf(childLambdaStates[c])));
                    }
                    best = packedStates.multiply(BigInteger.valueOf(2L * Double.BYTES));
                }
            }

            if (numChildren >= 2) {
                for (int streamed = 0; streamed < numChildren; streamed++) {
                    long streamedStates = saturatingMultiply(childRowStates[streamed],
                            childLambdaStates[streamed]);
                    if (streamedStates > childTableTotalStates) {
                        continue;
                    }
                    long residentStates = childTableTotalStates - streamedStates;
                    long fixed = DPGpuFullDP.estimateDeviceBytesForShape(
                            mPositionCount,
                            lambdaPositionCount,
                            totalLambdaStates,
                            lmTermCount,
                            lmPairCount,
                            childMTermTotal,
                            childLTermTotal,
                            numChildren,
                            residentStates,
                            outputTileMStates);
                    fixed = saturatingAdd(fixed,
                            saturatingMultiply(outputTileMStates, Long.BYTES));
                    fixed = saturatingAdd(fixed,
                            saturatingMultiply(2L * numChildren, Long.BYTES));
                    if (fixed >= gpuBudgetBytes) {
                        continue;
                    }
                    long rowBytes = saturatingMultiply(2L * Double.BYTES,
                            childLambdaStates[streamed]);
                    long streamBudget = Math.min(childSliceMaxBytes,
                            (gpuBudgetBytes - fixed) / 2L);
                    if (rowBytes <= 0L || rowBytes > streamBudget
                            || childLambdaStates[streamed]
                            > Integer.MAX_VALUE / (long)Double.BYTES) {
                        continue;
                    }
                    BigInteger hybridTraffic = BigInteger.valueOf(childTableTotalStates)
                            .multiply(BigInteger.valueOf(2L * Double.BYTES));
                    if (hybridTraffic.compareTo(best) < 0) {
                        best = hybridTraffic;
                    }
                }
            }

            DPGpuOutOfCore.Plan outOfCorePlan = DPGpuOutOfCore.choosePlan(
                    outOfCoreInput, gpuBudgetBytes);
            BigInteger outOfCoreTraffic = DPGpuOutOfCore.estimateTrafficBytes(
                    outOfCoreInput, outOfCorePlan);
            if (outOfCoreTraffic != null && outOfCoreTraffic.compareTo(best) < 0) {
                best = outOfCoreTraffic;
            }
            return best;
        }
    }

    private GpuEdgeShape estimateGpuEdgeShape(RootedTreeEdge edge, RCs rcs) {
        int[] mPos = edge.getMPositionsSorted();
        int[] lambdaPos = edge.getLambdaPositionsSorted();
        long childTableTotalStates = 0L;
        long childMTermTotal = 0L;
        long childLTermTotal = 0L;
        int numChildren = 0;
        boolean childFoldMapped = true;
        LinkedHashSet<RootedTreeEdge> children = edge.getFset();
        int childCount = children == null ? 0 : children.size();
        long[] childRowStates = new long[childCount];
        long[] childLambdaStates = new long[childCount];
        int[][] childMParentDims = new int[childCount][];
        int[][] childLambdaDims = new int[childCount][];
        boolean[] unionMPositions = new boolean[mPos.length];
        if (children != null) {
            int childIndex = 0;
            for (RootedTreeEdge child : children) {
                numChildren++;
                childTableTotalStates = saturatingAdd(childTableTotalStates,
                        child.getMStateCount());
                long rowStates = 1L;
                long lambdaStates = 1L;
                int[] childMParentDimsTmp = new int[child.getMPositionsSorted().length];
                int[] childLambdaDimsTmp = new int[child.getMPositionsSorted().length];
                int childMCount = 0;
                int childLambdaCount = 0;
                for (int childMPos : child.getMPositionsSorted()) {
                    int parentMIndex = indexOfPosition(mPos, childMPos);
                    if (parentMIndex >= 0) {
                        childMTermTotal = saturatingAdd(childMTermTotal, 1L);
                        unionMPositions[parentMIndex] = true;
                        childMParentDimsTmp[childMCount++] = parentMIndex;
                        rowStates = saturatingMultiply(rowStates, rcs.getNum(childMPos));
                    } else {
                        int parentLambdaIndex = indexOfPosition(lambdaPos, childMPos);
                        if (parentLambdaIndex >= 0) {
                            childLTermTotal = saturatingAdd(childLTermTotal, 1L);
                            childLambdaDimsTmp[childLambdaCount++] = parentLambdaIndex;
                            lambdaStates = saturatingMultiply(lambdaStates,
                                    rcs.getNum(childMPos));
                        } else {
                            childFoldMapped = false;
                        }
                    }
                }
                childMParentDims[childIndex] = Arrays.copyOf(
                        childMParentDimsTmp, childMCount);
                childLambdaDims[childIndex] = Arrays.copyOf(
                        childLambdaDimsTmp, childLambdaCount);
                childRowStates[childIndex] = rowStates;
                childLambdaStates[childIndex] = lambdaStates;
                childIndex++;
            }
        }
        long unionStateCount = 1L;
        long parentFreeStateCount = 1L;
        int[] parentMToUnion = new int[mPos.length];
        Arrays.fill(parentMToUnion, -1);
        int unionDimensionCount = 0;
        for (int i = 0; i < unionMPositions.length; i++) {
            if (unionMPositions[i]) {
                parentMToUnion[i] = unionDimensionCount++;
                unionStateCount = saturatingMultiply(unionStateCount, rcs.getNum(mPos[i]));
            } else {
                parentFreeStateCount = saturatingMultiply(parentFreeStateCount,
                        rcs.getNum(mPos[i]));
            }
        }
        int[] unionMCounts = new int[unionDimensionCount];
        for (int i = 0; i < mPos.length; i++) {
            int unionDimension = parentMToUnion[i];
            if (unionDimension >= 0) {
                unionMCounts[unionDimension] = rcs.getNum(mPos[i]);
            }
        }
        int[][] childMUnionDims = new int[childCount][];
        for (int c = 0; c < childCount; c++) {
            childMUnionDims[c] = new int[childMParentDims[c].length];
            for (int t = 0; t < childMParentDims[c].length; t++) {
                childMUnionDims[c][t] = parentMToUnion[childMParentDims[c][t]];
            }
        }

        int[] lambdaCounts = new int[lambdaPos.length];
        for (int i = 0; i < lambdaPos.length; i++) {
            lambdaCounts[i] = rcs.getNum(lambdaPos[i]);
        }

        long lmTermCount = 0L;
        int lmPairCount = 0;
        for (int lp : lambdaPos) {
            long lambdaCount = rcs.getNum(lp);
            for (int mp : mPos) {
                if (!interactionGraph.hasEdge(lp, mp)) {
                    continue;
                }
                lmPairCount++;
                lmTermCount = saturatingAdd(lmTermCount,
                        saturatingMultiply(lambdaCount, rcs.getNum(mp)));
            }
        }

        long outputTileMStates = resolveGpuOutputTileMStates(edge.getMStateCount());
        long childSliceMaxBytes = DPGpuFullDP.configuredChildSliceMaxBytes();
        long outOfCoreOutputWorkspaceMaxBytes =
                DPGpuFullDP.configuredOutOfCoreOutputWorkspaceMaxBytes();
        boolean supported = childFoldMapped
                && mPos.length <= DPGpuFullDP.MAX_EDGE_POSITIONS
                && lambdaPos.length <= DPGpuFullDP.MAX_EDGE_POSITIONS
                && numChildren <= 64
                && lmTermCount <= Integer.MAX_VALUE
                && childMTermTotal <= Integer.MAX_VALUE
                && childLTermTotal <= Integer.MAX_VALUE;
        DPGpuOutOfCore.PlanningInput outOfCoreInput = supported
                ? new DPGpuOutOfCore.PlanningInput(mPos.length, lambdaCounts,
                        unionMCounts, childMUnionDims, childLambdaDims,
                        parentFreeStateCount, childSliceMaxBytes,
                        outOfCoreOutputWorkspaceMaxBytes,
                        lmTermCount, lmPairCount,
                        childMTermTotal, childLTermTotal)
                : null;
        long fullDeviceBytes = childFoldMapped
                ? DPGpuFullDP.estimateDeviceBytesForShape(
                        mPos.length,
                        lambdaPos.length,
                        edge.getTotalLambdaStates(),
                        lmTermCount,
                        lmPairCount,
                        childMTermTotal,
                        childLTermTotal,
                        numChildren,
                        childTableTotalStates,
                        outputTileMStates)
                : Long.MAX_VALUE;
        long fixedSlicedDeviceBytes = DPGpuFullDP.estimateDeviceBytesForShape(
                mPos.length,
                lambdaPos.length,
                edge.getTotalLambdaStates(),
                lmTermCount,
                lmPairCount,
                childMTermTotal,
                childLTermTotal,
                numChildren,
                0L,
                outputTileMStates);
        fixedSlicedDeviceBytes = saturatingAdd(fixedSlicedDeviceBytes,
                saturatingMultiply(outputTileMStates, Long.BYTES));
        return new GpuEdgeShape(mPos.length, lambdaPos.length, edge.getTotalLambdaStates(),
                lmTermCount, lmPairCount, childMTermTotal, childLTermTotal, numChildren,
                childTableTotalStates, outputTileMStates, fixedSlicedDeviceBytes,
                unionStateCount, childRowStates, childLambdaStates,
                childSliceMaxBytes, outOfCoreInput, childFoldMapped, supported,
                fullDeviceBytes);
    }

    private static class ReuseStats {
        final int reusableEdges;
        final double reusableLogWork;
        final double totalLogWork;

        ReuseStats(int reusableEdges, double reusableLogWork, double totalLogWork) {
            this.reusableEdges = reusableEdges;
            this.reusableLogWork = reusableLogWork;
            this.totalLogWork = totalLogWork;
        }
    }

    private ReuseStats computeReuseStats(List<RootedTreeEdge> lambdaEdges,
                                         Set<Integer> mutablePositions, RCs rcs) {
        int reusableEdges = 0;
        double reusableLogWork = Double.NEGATIVE_INFINITY;
        double totalLogWork = Double.NEGATIVE_INFINITY;

        for (RootedTreeEdge edge : lambdaEdges) {
            double logWork = computeEdgeLogWork(edge, rcs);
            totalLogWork = RootedTreeEdge.logSumExp(totalLogWork, logWork);
            if (!edgeSubtreeTouchesMutable(edge, mutablePositions)) {
                reusableEdges++;
                reusableLogWork = RootedTreeEdge.logSumExp(reusableLogWork, logWork);
            }
        }

        return new ReuseStats(reusableEdges, reusableLogWork, totalLogWork);
    }

    private double computeEdgeLogWork(RootedTreeEdge edge, RCs rcs) {
        double logWork = 0.0;
        for (int pos : edge.getMPositionsSorted()) {
            logWork += Math.log(Math.max(1, rcs.getNum(pos)));
        }
        for (int pos : edge.getLambdaPositionsSorted()) {
            logWork += Math.log(Math.max(1, rcs.getNum(pos)));
        }
        return logWork;
    }

    private boolean edgeSubtreeTouchesMutable(RootedTreeEdge edge, Set<Integer> mutablePositions) {
        if (mutablePositions.isEmpty()) return false;
        for (int pos : edge.getM()) {
            if (mutablePositions.contains(pos)) return true;
        }
        if (edge.getL() != null) {
            for (int pos : edge.getL()) {
                if (mutablePositions.contains(pos)) return true;
            }
        }
        return false;
    }

    private Set<Integer> identifyMutablePositions(RCs rcs) {
        TreeSet<Integer> mutable = new TreeSet<>();

        String override = getConfigProperty(MUTABLE_POSITIONS_PROPERTY, null);
        if (override != null && !override.trim().isEmpty()) {
            for (String field : override.split(",")) {
                String trimmed = field.trim();
                if (trimmed.isEmpty()) continue;
                try {
                    int pos = Integer.parseInt(trimmed);
                    if (pos >= 0 && pos < rcs.getNumPos()) {
                        mutable.add(pos);
                    } else {
                        System.err.println(BranchDpConfig.getBackendLogPrefix() + " Ignoring mutable position " + pos
                                + " outside [0," + rcs.getNumPos() + ").");
                    }
                } catch (NumberFormatException e) {
                    System.err.println(BranchDpConfig.getBackendLogPrefix() + " Invalid mutable position '" + trimmed
                            + "' in " + MUTABLE_POSITIONS_PROPERTY + ", skipping.");
                }
            }
            return mutable;
        }

        if (confSpace != null && confSpace.positions != null) {
            int n = Math.min(rcs.getNumPos(), confSpace.positions.size());
            for (int pos = 0; pos < n; pos++) {
                if (confSpace.positions.get(pos).hasMutations()) {
                    mutable.add(pos);
                }
            }
        }
        return mutable;
    }

    private String formatPositionsWithResidues(Set<Integer> positions) {
        if (positions.isEmpty()) return "[]";
        StringBuilder sb = new StringBuilder("[");
        boolean first = true;
        for (int pos : positions) {
            if (!first) sb.append(',');
            first = false;
            sb.append(pos);
            if (confSpace != null && pos >= 0 && pos < confSpace.positions.size()) {
                sb.append(':').append(confSpace.positions.get(pos).resNum);
            }
        }
        sb.append(']');
        return sb.toString();
    }

    private boolean isBetterRooting(RootingCandidate candidate, RootingCandidate best,
                                    double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        boolean candidateCapped = candidate.maxFsetSize <= rootSplitMaxFset;
        boolean bestCapped = best.maxFsetSize <= rootSplitMaxFset;
        if (candidateCapped != bestCapped) return candidateCapped;

        if (candidate.maxFsetSize != best.maxFsetSize) {
            return candidate.maxFsetSize > best.maxFsetSize;
        }
        if (candidate.branchingEdges != best.branchingEdges) {
            return candidate.branchingEdges > best.branchingEdges;
        }
        if (candidate.rootFsetSize != best.rootFsetSize) {
            return candidate.rootFsetSize > best.rootFsetSize;
        }
        if (candidate.totalFsetEdges != best.totalFsetEdges) {
            return candidate.totalFsetEdges > best.totalFsetEdges;
        }
        int tessCmp = Double.compare(candidate.logTESS, best.logTESS);
        if (tessCmp != 0) return tessCmp < 0;
        return candidate.splitEdgeIndex < best.splitEdgeIndex;
    }

    private boolean isBetterReuseRooting(RootingCandidate candidate, RootingCandidate best,
                                         double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        boolean candidateCapped = candidate.maxFsetSize <= rootSplitMaxFset;
        boolean bestCapped = best.maxFsetSize <= rootSplitMaxFset;
        if (candidateCapped != bestCapped) return candidateCapped;

        int workRatioCmp = Double.compare(candidate.reusableWorkRatio(), best.reusableWorkRatio());
        if (workRatioCmp != 0) return workRatioCmp > 0;

        int reusableEdgesCmp = Integer.compare(candidate.reusableLambdaEdges, best.reusableLambdaEdges);
        if (reusableEdgesCmp != 0) return reusableEdgesCmp > 0;

        int edgeRatioCmp = Double.compare(candidate.reusableEdgeRatio(), best.reusableEdgeRatio());
        if (edgeRatioCmp != 0) return edgeRatioCmp > 0;

        int cleanWorkCmp = Double.compare(candidate.reusableLogWork, best.reusableLogWork);
        if (cleanWorkCmp != 0) return cleanWorkCmp > 0;

        return isBetterMemoryRooting(candidate, best, logNaive);
    }

    private boolean isBetterMemoryRooting(RootingCandidate candidate, RootingCandidate best,
                                          double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        int maxTableCmp = Long.compare(candidate.maxDPTableBytes, best.maxDPTableBytes);
        if (maxTableCmp != 0) return maxTableCmp < 0;

        int totalTableCmp = Long.compare(candidate.totalDPTableBytes, best.totalDPTableBytes);
        if (totalTableCmp != 0) return totalTableCmp < 0;

        int workCmp = Double.compare(candidate.logDPWork, best.logDPWork);
        if (workCmp != 0) return workCmp < 0;

        boolean candidateCapped = candidate.maxFsetSize <= rootSplitMaxFset;
        boolean bestCapped = best.maxFsetSize <= rootSplitMaxFset;
        if (candidateCapped != bestCapped) return candidateCapped;

        return isBetterRooting(candidate, best, logNaive);
    }

    private boolean isBetterGpuBytesRooting(RootingCandidate candidate, RootingCandidate best,
                                            double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        boolean candidateSupported = candidate.gpuUnsupportedEdges == 0;
        boolean bestSupported = best.gpuUnsupportedEdges == 0;
        if (candidateSupported != bestSupported) return candidateSupported;
        int unsupportedCmp = Integer.compare(candidate.gpuUnsupportedEdges, best.gpuUnsupportedEdges);
        if (unsupportedCmp != 0) return unsupportedCmp < 0;

        // Avoid child slicing when any enumerated root can keep every edge resident.
        // Within the same residency class, optimize the exact dominant CUDA work;
        // minimizing bytes further is only a tie-breaker.
        boolean candidateResident = candidate.fitsSingleGpu;
        boolean bestResident = best.fitsSingleGpu;
        if (candidateResident != bestResident) return candidateResident;

        if (!candidateResident) {
            int trafficCmp = candidate.estimatedSlicedTraffic.compareTo(best.estimatedSlicedTraffic);
            if (trafficCmp != 0) return trafficCmp < 0;
        }

        int gpuWorkCmp = candidate.gpuWork.compareTo(best.gpuWork);
        if (gpuWorkCmp != 0) return gpuWorkCmp < 0;

        int gpuBytesCmp = Long.compare(candidate.maxSingleGpuBytes, best.maxSingleGpuBytes);
        if (gpuBytesCmp != 0) return gpuBytesCmp < 0;

        int maxTableCmp = Long.compare(candidate.maxDPTableBytes, best.maxDPTableBytes);
        if (maxTableCmp != 0) return maxTableCmp < 0;

        boolean candidateCapped = candidate.maxFsetSize <= rootSplitMaxFset;
        boolean bestCapped = best.maxFsetSize <= rootSplitMaxFset;
        if (candidateCapped != bestCapped) return candidateCapped;

        return isBetterRooting(candidate, best, logNaive);
    }

    private boolean isBetterPredictedRooting(RootingCandidate candidate,
                                             RootingCandidate best,
                                             double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive)
                <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive)
                <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        int unsupportedCmp = Integer.compare(candidate.gpuUnsupportedEdges,
                best.gpuUnsupportedEdges);
        if (unsupportedCmp != 0) return unsupportedCmp < 0;

        double candidateSeconds = predictedSeconds(candidate);
        double bestSeconds = predictedSeconds(best);
        int predictedCmp = Double.compare(candidateSeconds, bestSeconds);
        if (predictedCmp != 0) return predictedCmp < 0;

        // If the hardware model is incomplete, or two roots have the same
        // modeled duration, retain the established GPU-aware deterministic
        // ordering as a tie-breaker.
        return isBetterGpuBytesRooting(candidate, best, logNaive);
    }

    private boolean isBetterWorkRooting(RootingCandidate candidate, RootingCandidate best,
                                        double logNaive) {
        if (best == null) return true;

        boolean candidateValid = Math.exp(candidate.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        boolean bestValid = Math.exp(best.logTESS - logNaive) <= TESS_FALLBACK_THRESHOLD;
        if (candidateValid != bestValid) return candidateValid;

        int workCmp = Double.compare(candidate.logDPWork, best.logDPWork);
        if (workCmp != 0) return workCmp < 0;

        int maxTableCmp = Long.compare(candidate.maxDPTableBytes, best.maxDPTableBytes);
        if (maxTableCmp != 0) return maxTableCmp < 0;

        int totalTableCmp = Long.compare(candidate.totalDPTableBytes, best.totalDPTableBytes);
        if (totalTableCmp != 0) return totalTableCmp < 0;

        boolean candidateCapped = candidate.maxFsetSize <= rootSplitMaxFset;
        boolean bestCapped = best.maxFsetSize <= rootSplitMaxFset;
        if (candidateCapped != bestCapped) return candidateCapped;

        return isBetterRooting(candidate, best, logNaive);
    }

    private double computeLogNaive(RCs rcs) {
        double logNaive = 0.0;
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            if (rcs.getNum(pos) > 0) {
                logNaive += Math.log(rcs.getNum(pos));
            }
        }
        return logNaive;
    }

    private void reportDPTooLarge(DPTableTooLargeException e) {
        System.err.println(BranchDpConfig.getBackendLogPrefix() + " status=TooLargeForDenseDP"
                + ", state=" + stateName
                + ", tableState=" + e.stateName
                + ", mStates=" + e.mStates
                + ", positions=" + Arrays.toString(e.mPositions)
                + ", estimatedFinalTableBytes=" + formatBytes(saturatingDPTableBytes(e.mStates))
                + ". Suggested actions: lower branchwidth/rooting pressure, run with "
                + "-Dbranchdp.dp.cache=false for dense-only runs, or enable the shard-backed DP table once available.");
    }

    /** Returns the worst-case single-GPU VRAM requirement (bytes) across GPU full-DP
     *  lambda edges, using the same device-buffer accounting as DPGpuFullDP. The
     *  parent output table is still not counted as a full table because the GPU
     *  path writes it through a bounded output tile. */
    private long logDPMemoryPredictions(String context, double logTESS, RCs rcs) {
        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);

        long totalFinalBytes = 0L;
        long totalCacheBytes = 0L;
        // The single-GPU VRAM budget for the GPU full-DP path is the complete set
        // of device buffers allocated by DPGpuFullDP for one edge on one GPU. The
        // dominant resident inputs are the child lower/upper tables, lambda-only
        // arrays, and lambda-M cross-term arrays; metadata/index buffers and the
        // bounded output tile are included too so this preflight matches the
        // runtime estimateDeviceBytes gate.
        long worstCaseSingleGpuBytes = 0L;
        int worstCaseEdgeId = -1;
        for (int edgeId = 0; edgeId < lambdaEdges.size(); edgeId++) {
            RootedTreeEdge edge = lambdaEdges.get(edgeId);
            long finalBytes = estimateDPTableBytes(edge.getMStateCount());
            long cacheBytes = (dpCacheEnabled && dpCacheMaxEntries > 0) ? finalBytes : 0L;
            totalFinalBytes = saturatingAdd(totalFinalBytes, finalBytes);
            totalCacheBytes = saturatingAdd(totalCacheBytes, cacheBytes);

            GpuEdgeShape gpuShape = estimateGpuEdgeShape(edge, rcs);
            long childTableTotalBytes = saturatingMultiply(2L * Double.BYTES,
                    gpuShape.childTableTotalStates);
            long lambdaOnlyBytes = saturatingMultiply(2L * Double.BYTES,
                    gpuShape.totalLambdaStates);
            long lmCrossBytes = saturatingMultiply(2L * Double.BYTES,
                    gpuShape.lmTermCount);
            long residentInputBytes = saturatingAdd(saturatingAdd(childTableTotalBytes, lambdaOnlyBytes), lmCrossBytes);
            long gpuDeviceBytes = gpuShape.fullDeviceBytes;
            long metadataAndTileBytes = gpuDeviceBytes > residentInputBytes
                    ? gpuDeviceBytes - residentInputBytes
                    : 0L;
            long singleGpuBytes = gpuDeviceBytes;
            if (singleGpuBytes > worstCaseSingleGpuBytes) {
                worstCaseSingleGpuBytes = singleGpuBytes;
                worstCaseEdgeId = edgeId;
            }

            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP memory prediction"
                    + " context=" + context
                    + ", edgeId=" + edgeId
                    + ", state=" + stateName
                    + ", M=" + Arrays.toString(edge.getMPositionsSorted())
                    + ", lambda=" + Arrays.toString(edge.getLambdaPositionsSorted())
                    + ", mStateCount=" + edge.getMStateCount()
                    + ", totalLambdaStates=" + edge.getTotalLambdaStates()
                    + ", finalTableBytes=" + formatBytes(finalBytes)
                    + ", cacheCopyBytes=" + formatBytes(cacheBytes)
                    + ", childTableTotalBytes=" + formatBytes(childTableTotalBytes)
                    + ", lambdaOnlyBytes=" + formatBytes(lambdaOnlyBytes)
                    + ", lmCrossBytes=" + formatBytes(lmCrossBytes)
                    + ", metadataAndTileBytes=" + formatBytes(metadataAndTileBytes)
                    + ", outputTileMStates=" + gpuShape.outputTileMStates
                    + ", childFoldMapped=" + gpuShape.childFoldMapped
                    + ", gpuSupported=" + gpuShape.supported
                    + ", gpuDeviceBytes=" + formatBytes(gpuDeviceBytes)
                    + ", singleGpuBytes=" + formatBytes(singleGpuBytes)
                    + ", branchwidth=" + branchwidth
                    + ", logTESS=" + String.format(Locale.ROOT, "%.2f", logTESS));
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP memory prediction summary"
                + " context=" + context
                + ", state=" + stateName
                + ", lambdaEdges=" + lambdaEdges.size()
                + ", totalFinalTableBytes=" + formatBytes(totalFinalBytes)
                + ", totalCacheCopyBytes=" + formatBytes(totalCacheBytes)
                + ", worstCaseSingleGpuBytes=" + formatBytes(worstCaseSingleGpuBytes)
                + " (edgeId=" + worstCaseEdgeId + ")"
                + ", dpCache=" + (dpCacheEnabled && dpCacheMaxEntries > 0)
                + ", dpCacheMaxTable=" + formatBytes(dpCacheMaxTableBytes)
                + ", dpCacheMaxTotal=" + formatBytes(dpCacheMaxTotalBytes)
                + ", branchwidth=" + branchwidth
                + ", logTESS=" + String.format(Locale.ROOT, "%.2f", logTESS));

        return worstCaseSingleGpuBytes;
    }

    private long resolveGpuOutputTileMStates(long mStateCount) {
        long requested = getConfigLong(DP_GPU_OUTPUT_TILE_MSTATES_PROPERTY,
                DEFAULT_DP_GPU_OUTPUT_TILE_MSTATES);
        long tile = Math.max(1L, requested);
        tile = Math.min(tile, (long) Integer.MAX_VALUE);
        return Math.min(tile, Math.max(1L, mStateCount));
    }

    private static int indexOfPosition(int[] positions, int target) {
        for (int i = 0; i < positions.length; i++) {
            if (positions[i] == target) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Sanity-check the worst-case single-GPU VRAM requirement (see
     * logDPMemoryPredictions) against the actually-available VRAM on this hardware,
     * BEFORE any expensive work (CCD, sampling, DP) begins. If it clearly will not fit,
     * abort now rather than burn GPU-hours only to hit the same wall hours or days later
     * (per-edge GPU OOM already fails fast once reached; this just moves that fast-fail
     * to the front of the run instead of wherever the offending edge happens to sit in
     * the tree). Only runs when GPU DP is actually requested; silently skipped if VRAM
     * could not be queried (e.g. no GPU visible) rather than false-aborting.
     */
    private void checkWorstCaseSingleGpuFits(long worstCaseSingleGpuBytes) {
        if (!getConfigBoolean("branchdp.dp.gpu", false)) {
            return; // GPU DP not requested: this VRAM concern doesn't apply
        }
        long usableVramBytes = DPGpuFullDP.queryMinUsableVramBytes();
        if (usableVramBytes < 0) {
            System.out.println(BranchDpConfig.getBackendLogPrefix()
                    + " worst-case single-GPU VRAM check skipped (could not query GPU memory)");
            return;
        }
        System.out.println(BranchDpConfig.getBackendLogPrefix()
                + " worst-case single-GPU VRAM check: need=" + formatBytes(worstCaseSingleGpuBytes)
                + ", usable=" + formatBytes(usableVramBytes));
        if (worstCaseSingleGpuBytes > usableVramBytes) {
            if (DPGpuFullDP.defaultChildSlicingEnabled()) {
                System.out.println(BranchDpConfig.getBackendLogPrefix()
                        + " worst-case full-resident GPU footprint exceeds single-GPU VRAM; "
                        + "child table slicing is enabled, so oversized child tables will be sliced at runtime");
                return;
            }
            throw new DPGpuFullDP.GpuMemoryExceededException(
                    "Preflight abort: worst-case single-GPU full-DP device-buffer requirement "
                    + worstCaseSingleGpuBytes + " B exceeds usable VRAM " + usableVramBytes + " B"
                    + " on this hardware. This includes child tables, lambda/lambda-M energy buffers,"
                    + " metadata, and the output tile for one GPU. Child tables are fully replicated"
                    + " onto every GPU used for an edge (multi-GPU only shards the OUTPUT range, not"
                    + " child tables), so more GPUs does not help this per-device footprint. Aborting"
                    + " before CCD/sampling/DP work begins instead of burning GPU-hours on a run that"
                    + " would fail (or silently crawl on CPU) once it reaches this edge.");
        }
    }

    private static long saturatingDPTableBytes(long mStates) {
        if (mStates > Long.MAX_VALUE / (2L * (long) Double.BYTES)) {
            return Long.MAX_VALUE;
        }
        return 2L * mStates * (long) Double.BYTES;
    }

    private static long saturatingAdd(long a, long b) {
        if (Long.MAX_VALUE - a < b) {
            return Long.MAX_VALUE;
        }
        return a + b;
    }

    private static long saturatingMultiply(long a, long b) {
        if (a == 0L || b == 0L) {
            return 0L;
        }
        if (a > Long.MAX_VALUE / b) {
            return Long.MAX_VALUE;
        }
        return a * b;
    }

    private DPCacheStats computeFullDPTables(String context, boolean allowCache) {
        if (!allowCache || !dpCacheEnabled || dpCacheMaxEntries == 0) {
            RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
            return null;
        }

        DPCacheStats stats = new DPCacheStats();
        IdentityHashMap<RootedTreeEdge, String> edgeKeys = new IdentityHashMap<>();
        String namespace = makeDPCacheNamespace(context);
        postOrderComputeFullDPWithCache(rootedRoot, namespace, edgeKeys, stats);
        return stats;
    }

    private void postOrderComputeFullDPWithCache(RootedTreeNode node, String namespace,
                                                 IdentityHashMap<RootedTreeEdge, String> edgeKeys,
                                                 DPCacheStats stats) {
        if (node == null) return;
        postOrderComputeFullDPWithCache(node.getLeftChild(), namespace, edgeKeys, stats);
        postOrderComputeFullDPWithCache(node.getRightChild(), namespace, edgeKeys, stats);

        RootedTreeEdge edge = node.getChildOfEdge();
        if (edge == null || !edge.getIsLambdaEdge()) return;

        if (!edge.hasDenseDPTable()) {
            edge.computeFullDP();
            edgeKeys.put(edge, buildDPCacheKey(namespace, edge, edgeKeys));
            stats.skippedLarge++;
            return;
        }

        String key = buildDPCacheKey(namespace, edge, edgeKeys);
        CachedDPTable cached = getCachedDPTable(key, edge.getMArraySize());
        if (cached != null) {
            System.arraycopy(cached.lower, 0, edge.getLogZLower(), 0, cached.lower.length);
            System.arraycopy(cached.upper, 0, edge.getLogZUpper(), 0, cached.upper.length);
            stats.hits++;
        } else {
            edge.computeFullDP();
            if (putCachedDPTable(key, edge.getLogZLower(), edge.getLogZUpper(), stats)) {
                stats.stores++;
            }
            stats.misses++;
        }
        edgeKeys.put(edge, key);
    }

    private CachedDPTable getCachedDPTable(String key, int expectedLength) {
        synchronized (DP_TABLE_CACHE_LOCK) {
            CachedDPTable cached = DP_TABLE_CACHE.get(key);
            if (cached == null) return null;
            if (cached.lower.length != expectedLength || cached.upper.length != expectedLength) {
                removeCachedDPTable(key, cached);
                return null;
            }
            return cached;
        }
    }

    private boolean putCachedDPTable(String key, double[] lower, double[] upper,
                                     DPCacheStats stats) {
        long tableBytes = estimateDPTableBytes(lower.length);
        if (dpCacheShouldSkip(lower.length, tableBytes, dpCacheSkipIfMStates,
                dpCacheMaxTableBytes, dpCacheMaxTotalBytes)) {
            stats.skippedLarge++;
            return false;
        }

        synchronized (DP_TABLE_CACHE_LOCK) {
            CachedDPTable previous = DP_TABLE_CACHE.remove(key);
            if (previous != null) {
                dpTableCacheBytes -= previous.bytes;
            }

            DP_TABLE_CACHE.put(key, new CachedDPTable(
                    Arrays.copyOf(lower, lower.length),
                    Arrays.copyOf(upper, upper.length),
                    tableBytes));
            dpTableCacheBytes += tableBytes;

            while (DP_TABLE_CACHE.size() > dpCacheMaxEntries
                    || dpTableCacheBytes > dpCacheMaxTotalBytes) {
                Iterator<Map.Entry<String, CachedDPTable>> it = DP_TABLE_CACHE.entrySet().iterator();
                if (!it.hasNext()) break;
                // Read the value straight off the entry (not via a second DP_TABLE_CACHE.get()
                // call): DP_TABLE_CACHE is access-order, so a get() here would touch
                // afterNodeAccess() and bump modCount out from under this same iterator,
                // throwing ConcurrentModificationException on the following it.remove().
                CachedDPTable evicted = it.next().getValue();
                it.remove();
                if (evicted != null) {
                    dpTableCacheBytes -= evicted.bytes;
                }
                stats.evictions++;
            }
        }
        return true;
    }

    private static void removeCachedDPTable(String key, CachedDPTable cached) {
        DP_TABLE_CACHE.remove(key);
        dpTableCacheBytes -= cached.bytes;
    }

    public static long estimateDPTableBytes(long mStates) {
        if (mStates > Long.MAX_VALUE / (2L * (long) Double.BYTES)) {
            return Long.MAX_VALUE;
        }
        return 2L * (long) mStates * (long) Double.BYTES;
    }

    /**
     * Estimate live Java-heap storage for one initialized lambda edge.
     * Besides the lower/upper DP table, int-addressable edges allocate one
     * {@code int} enumeration counter per M state and every edge owns two
     * lambda-only energy arrays. Object/shard headers are intentionally left to
     * the root budget's 30% heap headroom.
     */
    public static long estimateDPHostBytes(long mStates, long lambdaStates) {
        long tableBytes = estimateDPTableBytes(mStates);
        return saturatingAdd(tableBytes, estimateDPAuxHostBytes(mStates, lambdaStates));
    }

    /** Heap-only arrays retained beside a heap- or file-backed DP table. */
    public static long estimateDPAuxHostBytes(long mStates, long lambdaStates) {
        long enumerationBytes = mStates >= 0L && mStates <= Integer.MAX_VALUE
                ? saturatingMultiply(mStates, Integer.BYTES)
                : 0L;
        long lambdaBytes = lambdaStates <= 0L
                ? 0L
                : saturatingMultiply(lambdaStates, 2L * Double.BYTES);
        return saturatingAdd(enumerationBytes, lambdaBytes);
    }

    /**
     * Whether a DP table is too large to put in the byte-limited cache (design
     * 1.1): skip if its M-state count reaches skipIfMStates, or its estimated
     * byte size exceeds either the per-table or the total cache budget. Pure
     * decision function (no state); unit-tested in TestDPCacheBudget.
     */
    public static boolean dpCacheShouldSkip(long tableLen, long tableBytes, long skipIfMStates,
                                            long maxTableBytes, long maxTotalBytes) {
        return tableLen >= skipIfMStates
                || tableBytes > maxTableBytes
                || tableBytes > maxTotalBytes;
    }

    private static String formatBytes(long bytes) {
        if (bytes < 1024L) {
            return bytes + " B";
        }
        double value = bytes / 1024.0;
        String[] units = {"KiB", "MiB", "GiB", "TiB"};
        int unit = 0;
        while (value >= 1024.0 && unit < units.length - 1) {
            value /= 1024.0;
            unit++;
        }
        return String.format(Locale.ROOT, "%.1f %s", value, units[unit]);
    }

    private static String formatSlicedTraffic(BigInteger bytes) {
        return bytes.compareTo(SLICED_TRAFFIC_UNAVAILABLE) >= 0
                ? "unavailable"
                : bytes.toString();
    }

    private String formatDPCacheStats(DPCacheStats stats) {
        if (stats == null) return "";
        int size;
        long bytes;
        synchronized (DP_TABLE_CACHE_LOCK) {
            size = DP_TABLE_CACHE.size();
            bytes = dpTableCacheBytes;
        }
        return ", dpCacheHits=" + stats.hits
                + ", dpCacheMisses=" + stats.misses
                + ", dpCacheStores=" + stats.stores
                + ", dpCacheSkippedLarge=" + stats.skippedLarge
                + ", dpCacheEvictions=" + stats.evictions
                + ", dpCacheSize=" + size
                + ", dpCacheBytes=" + formatBytes(bytes);
    }

    private String makeDPCacheNamespace(String context) {
        StringBuilder sb = new StringBuilder();
        sb.append(context)
                .append("|conf=").append(System.identityHashCode(confSpace))
                .append("|rigid=").append(System.identityHashCode(rigidEmat))
                .append("|min=").append(System.identityHashCode(minimizingEmat))
                .append("|mode=").append(energyMode)
                .append("|prune=").append(sparsePruneThreshold)
                .append("|root=").append(selectedRootSplitIndex)
                .append("|graph=");
        for (int[] edge : interactionGraph.getEdgeList()) {
            sb.append(edge[0]).append('-').append(edge[1]).append(';');
        }
        return sb.toString();
    }

    private String buildDPCacheKey(String namespace, RootedTreeEdge edge,
                                   IdentityHashMap<RootedTreeEdge, String> edgeKeys) {
        StringBuilder sb = new StringBuilder(namespace);
        sb.append("|M=");
        appendPositionsWithRCs(sb, edge.getMPositionsSorted(), edge.getRCs());
        sb.append("|L=");
        appendPositionsWithRCs(sb, edge.getLambdaPositionsSorted(), edge.getRCs());
        sb.append("|F=");
        if (edge.getFset() != null) {
            for (RootedTreeEdge child : edge.getFset()) {
                String childKey = edgeKeys.get(child);
                sb.append('{').append(childKey == null ? "missing" : childKey).append("},");
            }
        }
        return sb.toString();
    }

    private void appendPositionsWithRCs(StringBuilder sb, int[] positions, RCs rcs) {
        for (int pos : positions) {
            sb.append(pos).append('[');
            int n = rcs.getNum(pos);
            for (int i = 0; i < n; i++) {
                if (i > 0) sb.append(',');
                sb.append(rcs.get(pos, i));
            }
            sb.append("];");
        }
    }

    /** Initialize search: create scorers and root search node */
    private void initSearch(RCs rcs) {
        int numPos = rcs.getNumPos();
        this.searchRCs = rcs;

        // Create scorers on the branch-pfunc emats (sparse-pruned in SPARSE mode,
        // identical to originals in FULL mode). A* and h-scoring therefore both
        // walk only over C_sparse and account only for E_sparse.
        gScorerMin = new PairwiseGScorer(branchMinimizingEmat);
        hScorerMin = new MPLPPairwiseHScorer(new EdgeUpdater(), branchMinimizingEmat, 1, 0.0001);
        gScorerRigid = new PairwiseGScorer(branchRigidEmat);
        hScorerNegRigid = new TraditionalPairwiseHScorer(
                new NegatedEnergyMatrix(confSpace, branchRigidEmat), rcs);

        // DP table precomputation: enabled in SPARSE mode (DP uses branch* emats with
        // cut pairs zeroed, computing exactly Z_sparse = Σ_sep Π_subtree Z_subtree(sep)).
        // In FULL mode, cut pairs are missing from the DP → invalid bound for full Z;
        // keep disabled there.
        if (useSparsePfunc) {
            long dpInitStart = System.currentTimeMillis();
            double rtForDP = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
            RootedTreeEdge.postOrderInitIncremental(rootedRoot,
                    branchRigidEmat, branchMinimizingEmat, interactionGraph, rtForDP,
                    false);
            DPCacheStats dpCacheStats = computeFullDPTables("initial", true);

            dpTablesReady = true;
            long dpInitTime = System.currentTimeMillis() - dpInitStart;
            // Diagnostic: report root edge DP values so we can spot -inf / 0 issues.
            String lzlStr = !rootedRootEdge.hasDPTable() || rootedRootEdge.getMStateCount() == 0
                    ? "null" : String.format("%.4f", rootedRootEdge.getLogZLower(0));
            String lzuStr = !rootedRootEdge.hasDPTable() || rootedRootEdge.getMStateCount() == 0
                    ? "null" : String.format("%.4f", rootedRootEdge.getLogZUpper(0));
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP tables computed in " + dpInitTime + " ms (SPARSE factored Z)"
                    + ", rootEdge logZ[0]: lower=" + lzlStr + ", upper=" + lzuStr
                    + ", rootEdge isLambdaEdge=" + rootedRootEdge.getIsLambdaEdge()
                    + ", mStateCount=" + rootedRootEdge.getMStateCount()
                    + formatDPCacheStats(dpCacheStats));
        } else {
            dpTablesReady = false;
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP table precomputation skipped (FULL mode, DP would be invalid)");
        }

        // Create root search node
        DecompSearchNode rootSearchNode = DecompSearchNode.makeRoot(
                rootedRootEdge, numPos, rcs,
                gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

        // Enhance root node's Z bounds using DP tables
        enhanceWithDPBounds(rootSearchNode);

        // Initialize priority queue
        decompQueue = new PriorityQueue<>();
        decompQueue.add(rootSearchNode);

        // Start at full parallelism so async loopTasks can saturate all threads
        // from the first leaf-dominant round (previously +1/round warmup wasted
        // CPUs while ramping up to parallelism.numThreads).
        maxMinimizations = Math.max(1, parallelism.numThreads);

        // Initialize Z bounds from root node
        flatSumZLower = rootSearchNode.subtreeLowerBound;
        flatSumZUpper = rootSearchNode.subtreeUpperBound;

        // Compute initial epsilon
        updateBound();

        // Count lambda-edges and print decomposition stats
        java.util.List<edu.duke.cs.osprey.branchdp.RootedTreeEdge> lambdaEdges = new java.util.ArrayList<>();
        edu.duke.cs.osprey.branchdp.RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);
        int leafEdgeCount = 0;
        long maxLambdaStates = 0;
        for (edu.duke.cs.osprey.branchdp.RootedTreeEdge e : lambdaEdges) {
            if (!e.hasFsetChildren()) leafEdgeCount++;
            maxLambdaStates = Math.max(maxLambdaStates, e.getTotalLambdaStates());
            // Debug: print each edge's lambda-set details
            int[] lpos = e.getLambdaPositionsSorted();
            int[] mpos = e.getMPositionsSorted();
            StringBuilder sb = new StringBuilder();
            sb.append("  Edge: lambda={");
            for (int i = 0; i < lpos.length; i++) {
                if (i > 0) sb.append(",");
                sb.append(lpos[i]).append("(").append(RCs.getNum(lpos[i])).append("rc)");
            }
            sb.append("}, m={");
            for (int i = 0; i < mpos.length; i++) {
                if (i > 0) sb.append(",");
                sb.append(mpos[i]).append("(").append(RCs.getNum(mpos[i])).append("rc)");
            }
            sb.append("}, totalLambda=").append(e.getTotalLambdaStates());
            sb.append(", fset=").append(e.getFset() == null ? 0 : e.getFset().size());
            sb.append(", leaf=").append(!e.hasFsetChildren());
            System.out.println(sb.toString());
        }
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " " + lambdaEdges.size() + " lambda-edges ("
                + leafEdgeCount + " leaf), maxLambdaStates=" + maxLambdaStates
                + ", pendingEdges=" + rootSearchNode.pendingEdges.size());

        // [DEBUG] Identify isolated positions (not covered by any lambda-edge)
        java.util.Set<Integer> coveredPos = new java.util.HashSet<>();
        for (edu.duke.cs.osprey.branchdp.RootedTreeEdge e : lambdaEdges) {
            for (int p : e.getLambdaPositionsSorted()) coveredPos.add(p);
            for (int p : e.getMPositionsSorted()) coveredPos.add(p);
        }
        java.util.List<Integer> isolatedPos = new java.util.ArrayList<>();
        for (int p = 0; p < rcs.getNumPos(); p++) {
            if (!coveredPos.contains(p)) isolatedPos.add(p);
        }
        StringBuilder isoSb = new StringBuilder();
        isoSb.append(BranchDpConfig.getBackendLogPrefix() + " [DEBUG] isolated positions=").append(isolatedPos.size())
                .append("/").append(rcs.getNumPos()).append(", list=[");
        for (int i = 0; i < isolatedPos.size(); i++) {
            int p = isolatedPos.get(i);
            if (i > 0) isoSb.append(",");
            isoSb.append(p).append("(").append(rcs.getNum(p)).append("rc)");
        }
        isoSb.append("]");
        System.out.println(isoSb.toString());

        // Also print edge incidence count per position (from sparse interaction graph)
        int[] degree = new int[rcs.getNumPos()];
        for (int[] edge : interactionGraph.getEdgeList()) {
            degree[edge[0]]++;
            degree[edge[1]]++;
        }
        StringBuilder degSb = new StringBuilder();
        degSb.append(BranchDpConfig.getBackendLogPrefix() + " [DEBUG] sparse-graph degree per position: [");
        for (int p = 0; p < rcs.getNumPos(); p++) {
            if (p > 0) degSb.append(",");
            degSb.append("p").append(p).append("=").append(degree[p]);
        }
        degSb.append("]");
        System.out.println(degSb.toString());

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Initialized. "
                + "Root node: upper=" + String.format("%.4e", flatSumZUpper.doubleValue())
                + ", lower=" + String.format("%.4e", flatSumZLower.doubleValue())
                + ", epsilon=" + String.format("%.6f", epsilonBound));
    }

    private double computeEpsilon() {
        return computeEpsilon(flatSumZLower, flatSumZUpper);
    }

    private double computeEpsilon(BigDecimal zLower, BigDecimal zUpper) {
        if (zUpper.signum() <= 0) return Double.POSITIVE_INFINITY;
        if (zLower.signum() <= 0) return 1.0;

        double eps = MathTools.bigDivide(zUpper.subtract(zLower), zUpper,
                PartitionFunction.decimalPrecision).doubleValue();
        return Math.max(0.0, eps);
    }

    private void auditFlatSumNumerics(String context, boolean force) {
        if (!numericAudit && !force) return;

        numericAuditChecks++;
        boolean invalidOrdering = flatSumZLower.compareTo(flatSumZUpper) > 0;
        boolean invalidSign = flatSumZLower.signum() < 0 || flatSumZUpper.signum() < 0;
        boolean invalidEpsilon = Double.isNaN(epsilonBound) || epsilonBound < -1e-12;
        boolean warning = invalidOrdering || invalidSign || invalidEpsilon;

        if (!force && !warning && numericAuditChecks % numericAuditInterval != 0) {
            return;
        }
        if (warning) {
            numericAuditWarnings++;
        }

        BigDecimal gap = flatSumZUpper.subtract(flatSumZLower);
        String doubleStatus = (Double.isFinite(flatSumZLower.doubleValue())
                && Double.isFinite(flatSumZUpper.doubleValue())
                && Double.isFinite(gap.doubleValue())) ? "finite" : "overflow";

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Numeric audit [" + context + "] "
                + "lower=" + formatBigExp(flatSumZLower)
                + ", upper=" + formatBigExp(flatSumZUpper)
                + ", gap=" + formatBigExp(gap)
                + ", epsilon=" + String.format("%.8f", epsilonBound)
                + ", double=" + doubleStatus
                + ", warnings=" + numericAuditWarnings);
    }

    private static String formatBigExp(BigDecimal value) {
        if (value == null) return "null";
        if (value.signum() == 0) return "0";
        return new BigExp(value).toString(4);
    }

    private static int countAssigned(int[] conf) {
        int count = 0;
        for (int rc : conf) {
            if (rc >= 0) count++;
        }
        return count;
    }

    private static BigInteger countStates(int[] positions, RCs rcs) {
        BigInteger count = BigInteger.ONE;
        for (int pos : positions) {
            count = count.multiply(BigInteger.valueOf(rcs.getNum(pos)));
        }
        return count;
    }

    private static String formatDelta(BigDecimal before, BigDecimal after) {
        return formatBigExp(before.subtract(after));
    }

    private void trace(String message) {
        if (trace) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " [TRACE] " + message);
        }
    }

    private boolean shouldPrintCorrectionDetail(long count) {
        return correctionAudit && (traceCorrectionLimit <= 0 || count <= traceCorrectionLimit);
    }

    private void correctionAuditLog(String message) {
        System.out.println(BranchDpConfig.getBackendLogPrefix() + " [CORRECTION_AUDIT] " + message);
    }

    // ========== Bound tightening ==========

    @Override
    protected void tightenBoundInPhases() {
        if (useBranchDecomposition) {
            tightenFlatSum();
        } else {
            super.tightenBoundInPhases();
        }
    }

    // ========== Search: Top-down decomposition-guided expand/minimize ==========

    /**
     * One round of the tighten loop.
     * Follows MARK*'s tightenBoundInPhases:
     *   1. Pull up to maxNodes internals + maxMinimizations leaves
     *   2. Apply correction matrix during pulling (corrected nodes put back)
     *   3. Compare internalZ vs leafZ, process only one type (error-driven selection)
     *   4. High-error internals: drill-down with full enum at each edge (no aggregates)
     *      Low-error internals: expand one edge with full enum (no aggregates)
     */
    private void tightenFlatSum() {
        roundCounter++;
        if (roundCounter % 1000 == 0) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " round=" + roundCounter
                    + ", eps=" + String.format("%.6f", epsilonBound)
                    + ", minimizations=" + totalMinimizations
                    + ", queueSize=" + decompQueue.size()
                    + ", enumSteps=" + totalEnumerationSteps
                    + ", drillDowns=" + totalDrillDowns
                    + ", drillDownEarlyStops=" + totalDrillDownEarlyStops);
        }

        // === Phase 1: Pull nodes (mirrors MARK*'s populateQueues) ===
        // Pull up to maxNodes internals and maxMinimizations leaves.
        // Corrected nodes are put back without classification.
        int maxNodes = 1000;
        int queueSizeBeforePull = decompQueue.size();
        List<DecompSearchNode> internalNodes = new ArrayList<>();
        List<DecompSearchNode> leafNodes = new ArrayList<>();
        List<DecompSearchNode> leftoverLeaves = new ArrayList<>();
        BigDecimal internalZ = BigDecimal.ZERO;
        BigDecimal leafZ = BigDecimal.ZERO;
        int correctionsThisRound = 0;
        BigDecimal correctionUpperDrop = BigDecimal.ZERO;
        BigDecimal correctionLowerDrop = BigDecimal.ZERO;

        while (!decompQueue.isEmpty() && internalNodes.size() < maxNodes) {
            DecompSearchNode node = decompQueue.poll();

            if (node.errorBound.signum() <= 0) continue;

            // In SPARSE mode, node Z bounds are certified by the branch DP tables
            // over E_sparse/C_sparse. The legacy correction matrix wraps the full
            // minimizing emat and would also overwrite DP-enhanced Z bounds, so only
            // allow this correction path in FULL mode.
            if (useHigherOrderCorrections && !useSparsePfunc && !node.isAggregate && !node.minimized) {
                correctionChecks++;
                double confCorrection = correctionMatrix.confE(node.partialConf);
                if (confCorrection > node.confLowerBound) {
                    double oldConfLower = node.confLowerBound;
                    double oldConfUpper = node.confUpperBound;
                    BigDecimal oldLowerZ = node.subtreeLowerBound;
                    BigDecimal oldUpperZ = node.subtreeUpperBound;

                    flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
                    flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

                    node.confLowerBound = Math.min(confCorrection, node.confUpperBound);
                    node.recomputeZBounds(bc);

                    flatSumZLower = flatSumZLower.add(node.subtreeLowerBound);
                    flatSumZUpper = flatSumZUpper.add(node.subtreeUpperBound);

                    correctionsThisRound++;
                    correctionApplications++;
                    correctionLowerDrop = correctionLowerDrop.add(oldLowerZ.subtract(node.subtreeLowerBound));
                    correctionUpperDrop = correctionUpperDrop.add(oldUpperZ.subtract(node.subtreeUpperBound));
                    if (trace && correctionApplications <= traceCorrectionLimit) {
                        trace("correction#" + correctionApplications
                                + " round=" + roundCounter
                                + " state=" + stateName
                                + " leaf=" + node.isLeaf()
                                + " assigned=" + countAssigned(node.partialConf) + "/" + node.partialConf.length
                                + " pending=" + node.pendingEdges.size()
                                + " confLower " + String.format("%.4f", oldConfLower)
                                + " -> " + String.format("%.4f", node.confLowerBound)
                                + " legacyCorrection=" + String.format("%.4f", confCorrection)
                                + " confUpper=" + String.format("%.4f", oldConfUpper)
                                + " zUpperDrop=" + formatDelta(oldUpperZ, node.subtreeUpperBound)
                                + " zLowerDrop=" + formatDelta(oldLowerZ, node.subtreeLowerBound)
                                + " conf=" + SimpleConfSpace.formatConfRCs(node.partialConf));
                    }
                    if (shouldPrintCorrectionDetail(correctionApplications)) {
                        correctionAuditLog("APPLY correction#" + correctionApplications
                                + " round=" + roundCounter
                                + " state=" + stateName
                                + " leaf=" + node.isLeaf()
                                + " assigned=" + countAssigned(node.partialConf) + "/" + node.partialConf.length
                                + " pending=" + node.pendingEdges.size()
                                + " confLower " + String.format("%.4f", oldConfLower)
                                + " -> " + String.format("%.4f", node.confLowerBound)
                                + " legacyCorrection=" + String.format("%.4f", confCorrection)
                                + " confUpper=" + String.format("%.4f", oldConfUpper)
                                + " zUpperDrop=" + formatDelta(oldUpperZ, node.subtreeUpperBound)
                                + " zLowerDrop=" + formatDelta(oldLowerZ, node.subtreeLowerBound)
                                + " conf=" + SimpleConfSpace.formatConfRCs(node.partialConf));
                    }

                    leftoverLeaves.add(node);
                    continue;
                } else {
                    correctionNotTighter++;
                    if (shouldPrintCorrectionDetail(correctionNotTighter)) {
                        correctionAuditLog("SKIP not-tighter check#" + correctionChecks
                                + " round=" + roundCounter
                                + " state=" + stateName
                                + " leaf=" + node.isLeaf()
                                + " assigned=" + countAssigned(node.partialConf) + "/" + node.partialConf.length
                                + " pending=" + node.pendingEdges.size()
                                + " confLower=" + String.format("%.4f", node.confLowerBound)
                                + " correction=" + String.format("%.4f", confCorrection)
                                + " confUpper=" + String.format("%.4f", node.confUpperBound)
                                + " conf=" + SimpleConfSpace.formatConfRCs(node.partialConf));
                    }
                }
            }

            if (node.isLeaf()) {
                if (leafNodes.size() < maxMinimizations) {
                    leafNodes.add(node);
                    leafZ = leafZ.add(node.errorBound);
                } else {
                    leftoverLeaves.add(node);
                }
            } else {
                internalNodes.add(node);
                internalZ = internalZ.add(node.errorBound);
            }
        }
        decompQueue.addAll(leftoverLeaves);

        updateBound();

        if (trace && (roundCounter <= traceRounds || correctionsThisRound > 0)) {
            trace("round=" + roundCounter
                    + " pulled queueBefore=" + queueSizeBeforePull
                    + " queueAfterPull=" + decompQueue.size()
                    + " internals=" + internalNodes.size()
                    + " leaves=" + leafNodes.size()
                    + " leftovers=" + leftoverLeaves.size()
                    + " maxMin=" + maxMinimizations
                    + " corrections=" + correctionsThisRound
                    + " corrUpperDrop=" + formatBigExp(correctionUpperDrop)
                    + " corrLowerDrop=" + formatBigExp(correctionLowerDrop)
                    + " internalZ=" + formatBigExp(internalZ)
                    + " leafZ=" + formatBigExp(leafZ)
                    + " flatLower=" + formatBigExp(flatSumZLower)
                    + " flatUpper=" + formatBigExp(flatSumZUpper)
                    + " eps=" + String.format("%.6f", epsilonBound));
        }

        // === Phase 2: Error-driven selection (mirrors MARK*'s decision logic) ===
        // Compare aggregate leaf error vs internal error; process only the dominant type.
        BigDecimal drillDownThreshold = new BigDecimal(1 - targetEpsilon);

        if (MathTools.isLessThan(internalZ, leafZ)) {
            if (trace && (roundCounter <= traceRounds || correctionsThisRound > 0)) {
                trace("round=" + roundCounter + " decision=leaf"
                        + " leafZ=" + formatBigExp(leafZ)
                        + " internalZ=" + formatBigExp(internalZ));
            }
            // Leaf error dominates: minimize leaves, put internals back
            // Submit all leaves asynchronously to loopTasks (parallel CCD min),
            // then wait for completion before continuing. Mirrors MARK*'s pattern.
            for (DecompSearchNode leaf : leafNodes) {
                minimizeDecompLeaf(leaf);
            }
            loopTasks.waitForFinish();
            decompQueue.addAll(internalNodes);
            // maxMinimizations is already at parallelism.numThreads (set in run()),
            // so no per-round ramp-up is needed.
        } else {
            if (trace && (roundCounter <= traceRounds || correctionsThisRound > 0)) {
                trace("round=" + roundCounter + " decision=internal"
                        + " leafZ=" + formatBigExp(leafZ)
                        + " internalZ=" + formatBigExp(internalZ));
            }
            // Internal error dominates: expand/drill-down internals, put leaves back
            if (parallelInternal && parallelism.numThreads > 1 && internalNodes.size() > 1) {
                processInternalNodesParallel(internalNodes, drillDownThreshold);
            } else {
                List<DecompSearchNode> newNodes = new ArrayList<>();
                boolean reachedTarget = false;
                for (int i = 0; i < internalNodes.size(); i++) {
                    DecompSearchNode internal = internalNodes.get(i);

                    // High error check (same criterion as MARK*'s boundLowestBoundConfUnderNode):
                    // lowerBound <= 1 AND upperBound/rootUpperBound > (1 - targetEpsilon)
                    if (!internal.isAggregate && shouldDrillDown(internal, drillDownThreshold)) {
                        // High error: drill-down through all pending edges with full enumeration.
                        // All siblings become explicit children (no aggregates created).
                        reachedTarget = drillDownFullEnum(internal, newNodes, drillDownThreshold);
                    } else {
                        // Low error: expand one pending edge with full enumeration (no aggregates)
                        expandDecompNode(internal);
                        updateBound();
                        reachedTarget = epsilonBound <= targetEpsilon;
                    }

                    if (reachedTarget) {
                        for (int j = i + 1; j < internalNodes.size(); j++) {
                            decompQueue.add(internalNodes.get(j));
                        }
                        break;
                    }
                }
                for (DecompSearchNode n : newNodes) {
                    if (n.errorBound.signum() > 0) {
                        decompQueue.add(n);
                    }
                }
            }
            decompQueue.addAll(leafNodes);
        }

        updateBound();
    }

    private boolean shouldDrillDown(DecompSearchNode node, BigDecimal drillDownThreshold) {
        return shouldDrillDown(node, flatSumZUpper, drillDownThreshold);
    }

    private boolean shouldDrillDown(DecompSearchNode node, BigDecimal zUpper,
                                    BigDecimal drillDownThreshold) {
        if (zUpper.signum() <= 0) return false;
        return !MathTools.isGreaterThan(node.subtreeLowerBound, BigDecimal.ONE) &&
                MathTools.isGreaterThan(
                        MathTools.bigDivide(node.subtreeUpperBound, zUpper,
                                PartitionFunction.decimalPrecision),
                        drillDownThreshold);
    }

    /**
     * Enhance a DecompSearchNode's Z bounds using pre-computed DP tables.
     *
     * Instead of the loose bound: numCompletions × boltz(g + h_position),
     * uses the factored bound: boltz(g_assigned) × prod_pending(exp(logZ_dp[mIdx])).
     *
     * Branch decomposition guarantees no direct interactions between lambda-sets of
     * sibling pending edges (all cross-subtree interactions go through M-set separators),
     * so the product factorization is exact.
     *
     * The DP bound replaces the h-score entirely: logZUpper/logZLower for each pending
     * edge already account for all lambda-state energies (one-body, lambda-lambda pairwise,
     * and lambda-M pairwise) within the edge's subtree.
     *
     * Applies the DP bound only if it is tighter than the existing g+h bound.
     */
    private void enhanceWithDPBounds(DecompSearchNode node) {
        // Only valid in SPARSE mode: DP tables were built using branch* emats (cut pairs = 0)
        // and the same interaction graph, so they exactly compute Z_sparse over completions.
        // In FULL mode the DP would miss cut-pair contributions and is therefore disabled.
        if (!useSparsePfunc || !dpTablesReady) return;
        if (node == null || node.isAggregate || node.minimized) return;
        if (node.pendingEdges == null || node.pendingEdges.isEmpty()) return;

        // Sum logZ_subtree across all pending edges. Branch decomposition guarantees these
        // subtrees share no direct pairwise interactions (all cross-subtree edges go via
        // M-set separators), so the product factorization is exact:
        //   Z_node = exp(-g_assigned/RT) × Π_{pe} Z_subtree(pe.mIdx).
        double logZUpperSum = 0.0;
        double logZLowerSum = 0.0;
        for (DecompSearchNode.PendingEdge pe : node.pendingEdges) {
            if (!pe.edge.hasDPTable()) return;
            if (pe.mIdx < 0 || pe.mIdx >= pe.edge.getMStateCount()) return;
            double lzU = pe.edge.getLogZUpper(pe.mIdx);
            double lzL = pe.edge.getLogZLower(pe.mIdx);
            if (Double.isNaN(lzU) || Double.isNaN(lzL)) return;
            // -inf logZ means subtree has no feasible completions; this prunes the node.
            if (lzU == Double.NEGATIVE_INFINITY) {
                node.subtreeUpperBound = BigDecimal.ZERO;
                node.subtreeLowerBound = BigDecimal.ZERO;
                node.errorBound = BigDecimal.ZERO;
                return;
            }
            logZUpperSum += lzU;
            // lzL may be -inf if the lower-bound DP has no contributions; treat as 0 weight.
            if (lzL != Double.NEGATIVE_INFINITY) {
                logZLowerSum += lzL;
            } else {
                logZLowerSum = Double.NEGATIVE_INFINITY;
            }
        }

        // g_assigned: energy over positions already pinned in partialConf, using the same
        // branch* emats that the DP tables are built from (so the assigned and DP halves
        // share a single energy model).
        ConfIndex confIndex = DecompSearchNode.buildConfIndex(node.partialConf,
                node.partialConf.length);
        double gMin = gScorerMin.calc(confIndex, searchRCs);
        double gRigid = gScorerRigid.calc(confIndex, searchRCs);

        // Z_upper_node = exp(-g_min/RT) × exp(logZUpperSum) → energy equivalent: g_min - RT*logZU
        // Z_lower_node = exp(-g_rigid/RT) × exp(logZLowerSum) → energy equivalent: g_rigid - RT*logZL
        double rtForDP = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
        double upperEnergy = gMin - rtForDP * logZUpperSum;
        BigDecimal dpZUpper = DecompSearchNode.safeBoltzCalc(bc, upperEnergy);

        BigDecimal dpZLower;
        if (logZLowerSum == Double.NEGATIVE_INFINITY) {
            dpZLower = BigDecimal.ZERO;
        } else {
            double lowerEnergy = gRigid - rtForDP * logZLowerSum;
            dpZLower = DecompSearchNode.safeBoltzCalc(bc, lowerEnergy);
        }

        // DP value is the EXACT factored Σ_{completions} exp(-E_sparse/RT) over the
        // sparse rigid/min emats. Always replace the loose g+h envelope with it; do NOT
        // gate on "only if tighter", since the g+h envelope can spuriously hit 0 when
        // MPLP touches +inf entries left over from cut-pair pruning, which would otherwise
        // veto the (correct) DP bound.
        node.subtreeUpperBound = dpZUpper;
        node.subtreeLowerBound = dpZLower;
        if (node.subtreeLowerBound.compareTo(node.subtreeUpperBound) > 0) {
            node.subtreeLowerBound = node.subtreeUpperBound;
        }
        node.errorBound = node.subtreeUpperBound.subtract(node.subtreeLowerBound);
        if (node.errorBound.signum() < 0) node.errorBound = BigDecimal.ZERO;
    }

    /** Number of lambda-states to produce per expand (top-k by energy). */
    private static final int LAMBDA_EXPAND_K = Integer.MAX_VALUE;

    private synchronized int pickPendingEdge(DecompSearchNode node) {
        edgeSelectionCalls++;
        int pendingCount = node.pendingEdges.size();
        maxPendingEdgesSeen = Math.max(maxPendingEdgesSeen, pendingCount);
        pendingEdgeCountHistogram.merge(pendingCount, 1L, Long::sum);

        if (pendingCount <= 1) {
            edgeSelectionSingleEdgeCalls++;
            return node.pickHighestErrorPendingEdge();
        }
        edgeSelectionMultiEdgeCalls++;

        if (edgeSelectionStrategy == EdgeSelectionStrategy.LAMBDA_STATES) {
            edgeSelectionLambdaStatesCalls++;
            return node.pickHighestErrorPendingEdge();
        }

        if (pendingCount > edgeLookaheadMaxPendingEdges) {
            edgeSelectionFallbackCalls++;
            return node.pickHighestErrorPendingEdge();
        }

        int maxStates = maxPendingLambdaStates(node);
        if (maxStates > edgeLookaheadMaxStates) {
            edgeSelectionFallbackCalls++;
            edgeSelectionFallbackMaxStates = Math.max(edgeSelectionFallbackMaxStates, maxStates);
            return node.pickHighestErrorPendingEdge();
        }

        edgeSelectionExactLookaheadCalls++;
        return pickBestPendingEdgeByContraction(node);
    }

    private int maxPendingLambdaStates(DecompSearchNode node) {
        int maxStates = 0;
        for (DecompSearchNode.PendingEdge pe : node.pendingEdges) {
            maxStates = Math.max(maxStates, pe.edge.getTotalLambdaStates());
        }
        return maxStates;
    }

    private int pickBestPendingEdgeByContraction(DecompSearchNode node) {
        int bestIdx = node.pickHighestErrorPendingEdge();
        BigDecimal bestDrop = null;
        double bestPerState = Double.NEGATIVE_INFINITY;
        int bestStates = 1;
        EdgeLookaheadResult[] lookahead = estimateChildErrorSums(node);

        for (int idx = 0; idx < node.pendingEdges.size(); idx++) {
            EdgeLookaheadResult result = lookahead[idx];
            BigDecimal childErrorSum = result.childErrorSum;
            BigDecimal drop = node.errorBound.subtract(childErrorSum);
            int states = Math.max(1, result.totalStates);

            if (edgeSelectionStrategy == EdgeSelectionStrategy.CONTRACTION_PER_STATE) {
                double perState = drop.doubleValue() / states;
                if (perState > bestPerState) {
                    bestPerState = perState;
                    bestDrop = drop;
                    bestStates = states;
                    bestIdx = idx;
                }
            } else if (bestDrop == null || drop.compareTo(bestDrop) > 0) {
                bestDrop = drop;
                bestStates = states;
                bestIdx = idx;
            }
        }

        if (bestDrop != null) {
            edgeLookaheadSelectedDrop = edgeLookaheadSelectedDrop.add(bestDrop);
            edgeLookaheadSelectedStates += bestStates;
        }

        return bestIdx;
    }

    private EdgeLookaheadResult[] estimateChildErrorSums(DecompSearchNode node) {
        int n = node.pendingEdges.size();
        EdgeLookaheadResult[] results;
        if (edgeLookaheadParallel && parallelism.numThreads > 1 && n > 1) {
            edgeLookaheadParallelCalls++;
            results = IntStream.range(0, n)
                    .parallel()
                    .mapToObj(idx -> estimateChildErrorSum(node, idx, true))
                    .toArray(EdgeLookaheadResult[]::new);
        } else {
            results = new EdgeLookaheadResult[n];
            for (int idx = 0; idx < n; idx++) {
                results[idx] = estimateChildErrorSum(node, idx, false);
            }
        }

        for (EdgeLookaheadResult result : results) {
            edgeLookaheadCandidates++;
            edgeLookaheadEnumeratedStates += result.enumeratedStates;
            edgeLookaheadNanos += result.nanos;
            if (result.cacheHit) {
                edgeLookaheadCacheHits++;
            } else {
                edgeLookaheadCacheMisses++;
            }
        }
        return results;
    }

    private EdgeLookaheadResult estimateChildErrorSum(DecompSearchNode node, int pendingIdx,
                                                      boolean useScorerCopies) {
        EdgeLookaheadResult cached = getCachedLookahead(node, pendingIdx);
        if (cached != null) {
            return cached.asCacheHit();
        }

        DecompSearchNode.PendingEdge pe = node.pendingEdges.get(pendingIdx);
        LambdaAStar lambdaAStar = new LambdaAStar(pe.edge, pe.mIdx, branchMinimizingEmat);
        BigDecimal childErrorSum = BigDecimal.ZERO;
        int enumerated = 0;
        long start = System.nanoTime();

        AStarScorer localGScorerMin = useScorerCopies ? gScorerMin.make() : gScorerMin;
        AStarScorer localHScorerMin = useScorerCopies ? hScorerMin.make() : hScorerMin;
        AStarScorer localGScorerRigid = useScorerCopies ? gScorerRigid.make() : gScorerRigid;
        AStarScorer localHScorerNegRigid = useScorerCopies ? hScorerNegRigid.make() : hScorerNegRigid;

        while (true) {
            int lambdaIdx = lambdaAStar.next();
            if (lambdaIdx < 0) break;
            enumerated++;

            DecompSearchNode child = DecompSearchNode.makeChild(
                    node, pendingIdx, lambdaIdx, RCs,
                    localGScorerMin, localHScorerMin, localGScorerRigid, localHScorerNegRigid, bc);

            enhanceWithDPBounds(child);
            childErrorSum = childErrorSum.add(child.errorBound);
        }

        EdgeLookaheadResult result = new EdgeLookaheadResult(
                childErrorSum, pe.edge.getTotalLambdaStates(), enumerated,
                System.nanoTime() - start, false);
        cacheLookahead(node, pendingIdx, result);
        return result;
    }

    private EdgeLookaheadResult getCachedLookahead(DecompSearchNode node, int pendingIdx) {
        synchronized (edgeLookaheadCache) {
            Map<Integer, EdgeLookaheadResult> nodeCache = edgeLookaheadCache.get(node);
            if (nodeCache == null) return null;
            return nodeCache.get(pendingIdx);
        }
    }

    private void cacheLookahead(DecompSearchNode node, int pendingIdx, EdgeLookaheadResult result) {
        synchronized (edgeLookaheadCache) {
            Map<Integer, EdgeLookaheadResult> nodeCache = edgeLookaheadCache.get(node);
            if (nodeCache == null) {
                nodeCache = new HashMap<>();
                edgeLookaheadCache.put(node, nodeCache);
            }
            nodeCache.put(pendingIdx, result);
        }
    }

    private void clearLookaheadCache(DecompSearchNode node) {
        synchronized (edgeLookaheadCache) {
            edgeLookaheadCache.remove(node);
        }
    }

    private ChildEnumerationBatch enumerateAllLambdaChildren(DecompSearchNode node, int pendingIdx,
                                                             RootedTreeEdge edge) {
        return enumerateAllLambdaChildren(node, pendingIdx, edge, parallelEnumeration);
    }

    private ChildEnumerationBatch enumerateAllLambdaChildren(DecompSearchNode node, int pendingIdx,
                                                             RootedTreeEdge edge, boolean allowParallel) {
        int totalStates = edge.getTotalLambdaStates();
        if (totalStates <= 0) {
            return new ChildEnumerationBatch(Collections.emptyList(), BigDecimal.ZERO, BigDecimal.ZERO, null);
        }

        if (!allowParallel || parallelism.numThreads <= 1 || totalStates == 1 || loopTasks == null) {
            ChildEnumerationChunk chunk = enumerateLambdaChildChunk(
                    0, node, pendingIdx, 0, totalStates, !allowParallel && parallelism.numThreads > 1);
            return new ChildEnumerationBatch(chunk.children, chunk.lowerSum, chunk.upperSum, chunk.bestChild);
        }

        int numChunks = Math.min(parallelism.numThreads, totalStates);
        int chunkSize = (totalStates + numChunks - 1) / numChunks;
        ChildEnumerationChunk[] chunks = new ChildEnumerationChunk[numChunks];

        for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
            int startLambdaIdx = chunkIndex * chunkSize;
            int endLambdaIdx = Math.min(totalStates, startLambdaIdx + chunkSize);
            final int taskChunkIndex = chunkIndex;
            loopTasks.submit(() -> enumerateLambdaChildChunk(
                    taskChunkIndex, node, pendingIdx, startLambdaIdx, endLambdaIdx, true),
                    chunk -> chunks[chunk.chunkIndex] = chunk);
        }
        loopTasks.waitForFinish();

        List<DecompSearchNode> children = new ArrayList<>(totalStates);
        BigDecimal lowerSum = BigDecimal.ZERO;
        BigDecimal upperSum = BigDecimal.ZERO;
        DecompSearchNode bestChild = null;
        double bestChildLower = Double.POSITIVE_INFINITY;

        for (ChildEnumerationChunk chunk : chunks) {
            if (chunk == null) {
                throw new IllegalStateException(BranchDpConfig.getBackendLogPrefix() + " missing lambda enumeration chunk result");
            }
            children.addAll(chunk.children);
            lowerSum = lowerSum.add(chunk.lowerSum);
            upperSum = upperSum.add(chunk.upperSum);
            if (chunk.bestChild != null && chunk.bestChild.confLowerBound < bestChildLower) {
                bestChild = chunk.bestChild;
                bestChildLower = chunk.bestChild.confLowerBound;
            }
        }

        return new ChildEnumerationBatch(children, lowerSum, upperSum, bestChild);
    }

    private ChildEnumerationChunk enumerateLambdaChildChunk(int chunkIndex, DecompSearchNode node,
                                                           int pendingIdx, int startLambdaIdx,
                                                           int endLambdaIdx, boolean useScorerCopies) {
        AStarScorer localGScorerMin = useScorerCopies ? gScorerMin.make() : gScorerMin;
        AStarScorer localHScorerMin = useScorerCopies ? hScorerMin.make() : hScorerMin;
        AStarScorer localGScorerRigid = useScorerCopies ? gScorerRigid.make() : gScorerRigid;
        AStarScorer localHScorerNegRigid = useScorerCopies ? hScorerNegRigid.make() : hScorerNegRigid;

        List<DecompSearchNode> children = new ArrayList<>(Math.max(0, endLambdaIdx - startLambdaIdx));
        BigDecimal lowerSum = BigDecimal.ZERO;
        BigDecimal upperSum = BigDecimal.ZERO;
        DecompSearchNode bestChild = null;
        double bestChildLower = Double.POSITIVE_INFINITY;

        for (int lambdaIdx = startLambdaIdx; lambdaIdx < endLambdaIdx; lambdaIdx++) {
            DecompSearchNode child = DecompSearchNode.makeChild(
                    node, pendingIdx, lambdaIdx, RCs,
                    localGScorerMin, localHScorerMin, localGScorerRigid, localHScorerNegRigid, bc);
            enhanceWithDPBounds(child);

            children.add(child);
            lowerSum = lowerSum.add(child.subtreeLowerBound);
            upperSum = upperSum.add(child.subtreeUpperBound);
            if (child.confLowerBound < bestChildLower) {
                bestChild = child;
                bestChildLower = child.confLowerBound;
            }
        }

        return new ChildEnumerationChunk(chunkIndex, children, lowerSum, upperSum, bestChild);
    }

    private void processInternalNodesParallel(List<DecompSearchNode> internalNodes,
                                              BigDecimal drillDownThreshold) {
        final BigDecimal lowerSnapshot = flatSumZLower;
        final BigDecimal upperSnapshot = flatSumZUpper;

        InternalNodeResult[] results = IntStream.range(0, internalNodes.size())
                .parallel()
                .mapToObj(idx -> processInternalNodeLocal(
                        idx, internalNodes.get(idx), drillDownThreshold, lowerSnapshot, upperSnapshot))
                .toArray(InternalNodeResult[]::new);

        boolean reachedTarget = false;
        for (InternalNodeResult result : results) {
            flatSumZLower = flatSumZLower.add(result.lowerDelta);
            flatSumZUpper = flatSumZUpper.add(result.upperDelta);
            totalEnumerationSteps += result.enumerationSteps;
            totalDrillDowns += result.drillDowns;
            totalDrillDownEarlyStops += result.drillDownEarlyStops;
            reachedTarget |= result.reachedTarget;

            for (DecompSearchNode node : result.generatedNodes) {
                if (node.errorBound.signum() > 0) {
                    decompQueue.add(node);
                }
            }
        }

        if (reachedTarget) {
            epsilonBound = computeEpsilon();
        }
    }

    private InternalNodeResult processInternalNodeLocal(int nodeIndex, DecompSearchNode internal,
                                                       BigDecimal drillDownThreshold,
                                                       BigDecimal lowerSnapshot,
                                                       BigDecimal upperSnapshot) {
        if (!internal.isAggregate && shouldDrillDown(internal, upperSnapshot, drillDownThreshold)) {
            return drillDownFullEnumLocal(
                    nodeIndex, internal, drillDownThreshold, lowerSnapshot, upperSnapshot);
        }
        return expandDecompNodeLocal(nodeIndex, internal);
    }

    private InternalNodeResult expandDecompNodeLocal(int nodeIndex, DecompSearchNode node) {
        if (node.isAggregate) {
            return expandAggregateNodeLocal(nodeIndex, node);
        }

        int pendingIdx = pickPendingEdge(node);
        clearLookaheadCache(node);
        DecompSearchNode.PendingEdge pe = node.pendingEdges.get(pendingIdx);
        RootedTreeEdge edge = pe.edge;

        List<DecompSearchNode> generatedNodes = new ArrayList<>();
        BigDecimal lowerDelta = node.subtreeLowerBound.negate();
        BigDecimal upperDelta = node.subtreeUpperBound.negate();
        BigDecimal childrenZUpper = BigDecimal.ZERO;

        int totalLambdaStates = edge.getTotalLambdaStates();
        int k = Math.min(LAMBDA_EXPAND_K, totalLambdaStates);

        if (k >= totalLambdaStates) {
            ChildEnumerationBatch batch = enumerateAllLambdaChildren(node, pendingIdx, edge, false);
            generatedNodes.addAll(batch.children);
            lowerDelta = lowerDelta.add(batch.lowerSum);
            upperDelta = upperDelta.add(batch.upperSum);
            return new InternalNodeResult(
                    nodeIndex, generatedNodes, lowerDelta, upperDelta,
                    false, 1, 0, 0);
        }

        LambdaAStar lambdaAStar = new LambdaAStar(edge, pe.mIdx, branchMinimizingEmat);
        AStarScorer localGScorerMin = gScorerMin.make();
        AStarScorer localHScorerMin = hScorerMin.make();
        AStarScorer localGScorerRigid = gScorerRigid.make();
        AStarScorer localHScorerNegRigid = hScorerNegRigid.make();

        for (int i = 0; i < k; i++) {
            int lambdaIdx = lambdaAStar.next();
            if (lambdaIdx < 0) break;

            DecompSearchNode child = DecompSearchNode.makeChild(
                    node, pendingIdx, lambdaIdx, RCs,
                    localGScorerMin, localHScorerMin, localGScorerRigid, localHScorerNegRigid, bc);
            enhanceWithDPBounds(child);

            generatedNodes.add(child);
            lowerDelta = lowerDelta.add(child.subtreeLowerBound);
            upperDelta = upperDelta.add(child.subtreeUpperBound);
            childrenZUpper = childrenZUpper.add(child.subtreeUpperBound);
        }

        if (!lambdaAStar.isExhausted()) {
            DecompSearchNode aggregate = DecompSearchNode.makeAggregate(
                    node, pendingIdx, lambdaAStar,
                    node.subtreeUpperBound, childrenZUpper, bc);
            generatedNodes.add(aggregate);
            lowerDelta = lowerDelta.add(aggregate.subtreeLowerBound);
            upperDelta = upperDelta.add(aggregate.subtreeUpperBound);
        }

        return new InternalNodeResult(
                nodeIndex, generatedNodes, lowerDelta, upperDelta,
                false, 1, 0, 0);
    }

    private InternalNodeResult expandAggregateNodeLocal(int nodeIndex, DecompSearchNode aggregate) {
        BigDecimal lowerDelta = aggregate.subtreeLowerBound.negate();
        BigDecimal upperDelta = aggregate.subtreeUpperBound.negate();
        List<DecompSearchNode> generatedNodes = new ArrayList<>();

        DecompSearchNode child = aggregate.popNextFromAggregate(
                RCs, gScorerMin.make(), hScorerMin.make(), gScorerRigid.make(), hScorerNegRigid.make(), bc);

        if (child != null) {
            enhanceWithDPBounds(child);
            generatedNodes.add(child);
            lowerDelta = lowerDelta.add(child.subtreeLowerBound);
            upperDelta = upperDelta.add(child.subtreeUpperBound);

            DecompSearchNode.computeAggregateZBounds(aggregate, aggregate.aggregateParent, bc);
            if (aggregate.errorBound.signum() > 0 && !aggregate.lambdaAStar.isExhausted()) {
                generatedNodes.add(aggregate);
                lowerDelta = lowerDelta.add(aggregate.subtreeLowerBound);
                upperDelta = upperDelta.add(aggregate.subtreeUpperBound);
            }
        }

        return new InternalNodeResult(
                nodeIndex, generatedNodes, lowerDelta, upperDelta,
                false, 0, 0, 0);
    }

    private InternalNodeResult drillDownFullEnumLocal(int nodeIndex, DecompSearchNode node,
                                                     BigDecimal drillDownThreshold,
                                                     BigDecimal lowerSnapshot,
                                                     BigDecimal upperSnapshot) {
        List<DecompSearchNode> generatedNodes = new ArrayList<>();
        BigDecimal lowerDelta = node.subtreeLowerBound.negate();
        BigDecimal upperDelta = node.subtreeUpperBound.negate();
        BigDecimal localLower = lowerSnapshot.add(lowerDelta);
        BigDecimal localUpper = upperSnapshot.add(upperDelta);
        int enumerationSteps = 0;
        int earlyStops = 0;

        DecompSearchNode current = node;

        while (!current.isLeaf()) {
            int pendingIdx = pickPendingEdge(current);
            clearLookaheadCache(current);
            DecompSearchNode.PendingEdge pe = current.pendingEdges.get(pendingIdx);
            RootedTreeEdge edge = pe.edge;

            enumerationSteps++;

            ChildEnumerationBatch batch = enumerateAllLambdaChildren(current, pendingIdx, edge, false);
            DecompSearchNode bestChild = batch.bestChild;

            lowerDelta = lowerDelta.add(batch.lowerSum);
            upperDelta = upperDelta.add(batch.upperSum);
            localLower = localLower.add(batch.lowerSum);
            localUpper = localUpper.add(batch.upperSum);
            generatedNodes.addAll(batch.children);

            if (bestChild == null) {
                boolean reachedTarget = computeEpsilon(localLower, localUpper) <= targetEpsilon;
                return new InternalNodeResult(
                        nodeIndex, generatedNodes, lowerDelta, upperDelta,
                        reachedTarget, enumerationSteps, 1, earlyStops);
            }

            // Early-stop removed: always drill to a leaf (mirrors drillDownFullEnum).

            generatedNodes.remove(bestChild);
            lowerDelta = lowerDelta.subtract(bestChild.subtreeLowerBound);
            upperDelta = upperDelta.subtract(bestChild.subtreeUpperBound);
            localLower = localLower.subtract(bestChild.subtreeLowerBound);
            localUpper = localUpper.subtract(bestChild.subtreeUpperBound);

            current = bestChild;
        }

        lowerDelta = lowerDelta.add(current.subtreeLowerBound);
        upperDelta = upperDelta.add(current.subtreeUpperBound);
        localLower = localLower.add(current.subtreeLowerBound);
        localUpper = localUpper.add(current.subtreeUpperBound);
        generatedNodes.add(current);

        boolean reachedTarget = computeEpsilon(localLower, localUpper) <= targetEpsilon;
        return new InternalNodeResult(
                nodeIndex, generatedNodes, lowerDelta, upperDelta,
                reachedTarget, enumerationSteps, 1, earlyStops);
    }

    /**
     * Expand a decomposition search node: assign one pending edge's lambda-set.
     *
     * Full enumeration is scored in parallel. If LAMBDA_EXPAND_K is lowered later,
     * the old LambdaAStar top-k path still creates an aggregate for the remaining states.
     */
    private void expandDecompNode(DecompSearchNode node) {
        if (node.isAggregate) {
            expandAggregateNode(node);
            return;
        }

        // Pick the pending edge to expand
        int pendingIdx = pickPendingEdge(node);
        clearLookaheadCache(node);
        DecompSearchNode.PendingEdge pe = node.pendingEdges.get(pendingIdx);
        RootedTreeEdge edge = pe.edge;

        // Remove parent's contribution from running Z totals
        flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

        totalEnumerationSteps++;

        int totalLambdaStates = edge.getTotalLambdaStates();
        int k = Math.min(LAMBDA_EXPAND_K, totalLambdaStates);
        BigDecimal childrenZUpper = BigDecimal.ZERO;

        if (k >= totalLambdaStates) {
            BigDecimal parentLower = node.subtreeLowerBound;
            BigDecimal parentUpper = node.subtreeUpperBound;
            ChildEnumerationBatch batch = enumerateAllLambdaChildren(node, pendingIdx, edge);
            flatSumZLower = flatSumZLower.add(batch.lowerSum);
            flatSumZUpper = flatSumZUpper.add(batch.upperSum);
            childrenZUpper = childrenZUpper.add(batch.upperSum);
            if (trace && totalEnumerationSteps <= traceRounds) {
                trace("expand enumStep=" + totalEnumerationSteps
                        + " state=" + stateName
                        + " assigned=" + countAssigned(node.partialConf) + "/" + node.partialConf.length
                        + " lambdaStates=" + totalLambdaStates
                        + " children=" + batch.children.size()
                        + " parentUpper=" + formatBigExp(parentUpper)
                        + " childUpperSum=" + formatBigExp(batch.upperSum)
                        + " upperDrop=" + formatDelta(parentUpper, batch.upperSum)
                        + " parentLower=" + formatBigExp(parentLower)
                        + " childLowerSum=" + formatBigExp(batch.lowerSum)
                        + " bestLower=" + (batch.bestChild == null ? "none"
                                : String.format("%.4f", batch.bestChild.confLowerBound)));
            }
            for (DecompSearchNode child : batch.children) {
                if (child.errorBound.signum() > 0) {
                    decompQueue.add(child);
                }
            }
            return;
        }

        // Lazy A* enumeration: produce top-k, aggregate the rest
        LambdaAStar lambdaAStar = new LambdaAStar(edge, pe.mIdx, branchMinimizingEmat);

        for (int i = 0; i < k; i++) {
            int lambdaIdx = lambdaAStar.next();
            if (lambdaIdx < 0) break;

            DecompSearchNode child = DecompSearchNode.makeChild(
                    node, pendingIdx, lambdaIdx, RCs,
                    gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

            // Tighten Z bounds using DP table lookups for remaining pending edges
            enhanceWithDPBounds(child);

            flatSumZLower = flatSumZLower.add(child.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.add(child.subtreeUpperBound);
            childrenZUpper = childrenZUpper.add(child.subtreeUpperBound);

            if (child.errorBound.signum() > 0) {
                decompQueue.add(child);
            }
        }

        // Create aggregate node for remaining lambda-states
        if (!lambdaAStar.isExhausted()) {
            DecompSearchNode aggregate = DecompSearchNode.makeAggregate(
                    node, pendingIdx, lambdaAStar,
                    node.subtreeUpperBound, childrenZUpper, bc);

            flatSumZLower = flatSumZLower.add(aggregate.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.add(aggregate.subtreeUpperBound);

            if (aggregate.errorBound.signum() > 0) {
                decompQueue.add(aggregate);
            }
        }
    }

    /**
     * Expand an aggregate node: pop 1 more lambda-state from its A* search.
     * The new child is added to the queue, and the aggregate's bounds are updated.
     */
    private void expandAggregateNode(DecompSearchNode aggregate) {
        // Pop one more lambda-state from the A*
        DecompSearchNode child = aggregate.popNextFromAggregate(
                RCs, gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

        if (child == null) {
            // A* exhausted — aggregate contributes nothing more
            flatSumZLower = flatSumZLower.subtract(aggregate.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.subtract(aggregate.subtreeUpperBound);
            return;
        }

        // Update running Z totals:
        // 1. Remove old aggregate contribution
        flatSumZLower = flatSumZLower.subtract(aggregate.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(aggregate.subtreeUpperBound);

        // 2. Tighten child's Z bounds using DP tables, then add contribution
        enhanceWithDPBounds(child);
        flatSumZLower = flatSumZLower.add(child.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.add(child.subtreeUpperBound);

        if (child.errorBound.signum() > 0) {
            decompQueue.add(child);
        }

        // 3. Recompute aggregate bounds from updated A* frontier
        DecompSearchNode.computeAggregateZBounds(aggregate, aggregate.aggregateParent, bc);

        // 4. Add updated aggregate back
        flatSumZLower = flatSumZLower.add(aggregate.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.add(aggregate.subtreeUpperBound);

        if (aggregate.errorBound.signum() > 0 && !aggregate.lambdaAStar.isExhausted()) {
            decompQueue.add(aggregate);
        }
    }

    /**
     * Drill-down with full enumeration at each edge (no aggregates).
     *
     * Like MARK*'s boundLowestBoundConfUnderNode: greedily pick the best child at
     * each edge, but enumerate ALL lambda-states as explicit children (siblings go
     * to newNodes for later queue insertion). The best child continues drilling.
     *
     * This produces one leaf per call, plus all siblings at every level as explicit
     * nodes with tight individual bounds. No aggregate nodes are created.
     *
     * @param node     the internal node to drill down from
     * @param newNodes output list: all generated children (siblings + final leaf)
     * @return true if the current flat-sum certificate has reached target epsilon
     */
    private boolean drillDownFullEnum(DecompSearchNode node, List<DecompSearchNode> newNodes,
                                      BigDecimal drillDownThreshold) {
        totalDrillDowns++;

        // Remove the original node's contribution (replaced by its children)
        flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

        DecompSearchNode current = node;

        while (!current.isLeaf()) {
            int pendingIdx = pickPendingEdge(current);
            clearLookaheadCache(current);
            DecompSearchNode.PendingEdge pe = current.pendingEdges.get(pendingIdx);
            RootedTreeEdge edge = pe.edge;

            totalEnumerationSteps++;

            // Enumerate ALL lambda-states for this edge (no aggregate)
            BigDecimal parentLower = current.subtreeLowerBound;
            BigDecimal parentUpper = current.subtreeUpperBound;
            ChildEnumerationBatch batch = enumerateAllLambdaChildren(current, pendingIdx, edge);
            DecompSearchNode bestChild = batch.bestChild;

            flatSumZLower = flatSumZLower.add(batch.lowerSum);
            flatSumZUpper = flatSumZUpper.add(batch.upperSum);
            newNodes.addAll(batch.children);

            if (trace && totalEnumerationSteps <= traceRounds) {
                trace("drill enumStep=" + totalEnumerationSteps
                        + " state=" + stateName
                        + " assigned=" + countAssigned(current.partialConf) + "/" + current.partialConf.length
                        + " lambdaStates=" + edge.getTotalLambdaStates()
                        + " children=" + batch.children.size()
                        + " parentUpper=" + formatBigExp(parentUpper)
                        + " childUpperSum=" + formatBigExp(batch.upperSum)
                        + " upperDrop=" + formatDelta(parentUpper, batch.upperSum)
                        + " parentLower=" + formatBigExp(parentLower)
                        + " childLowerSum=" + formatBigExp(batch.lowerSum)
                        + " bestLower=" + (bestChild == null ? "none"
                                : String.format("%.4f", bestChild.confLowerBound))
                        + " bestAssigned=" + (bestChild == null ? 0
                                : countAssigned(bestChild.partialConf)));
            }

            if (bestChild == null) {
                epsilonBound = computeEpsilon();
                return epsilonBound <= targetEpsilon;
            }

            // Early-stop removed: always drill to a leaf to enable leaf minimization,
            // which is the dominant gap-closing mechanism (matches pre-regression behavior).

            // Best child continues drilling; remove it from newNodes and flat sum
            // (it will be re-added at the next level or as final leaf)
            newNodes.remove(bestChild);
            flatSumZLower = flatSumZLower.subtract(bestChild.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.subtract(bestChild.subtreeUpperBound);

            current = bestChild;
        }

        // current is now a leaf — add to flat sum and newNodes
        // (will be pulled and minimized in a future round when leafZ dominates)
        flatSumZLower = flatSumZLower.add(current.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.add(current.subtreeUpperBound);
        newNodes.add(current);

        epsilonBound = computeEpsilon();
        return epsilonBound <= targetEpsilon;
    }

    /**
     * Minimize a leaf node (full conformation).
     * calcEnergy + triple corrections.
     */
    private void minimizeDecompLeaf(DecompSearchNode leaf) {
        final int[] fullConf = leaf.partialConf;

        // Check if already minimized (main thread, before async dispatch).
        final String confKey = SimpleConfSpace.formatConfRCs(fullConf);
        if (minimizedConfs.contains(confKey)) {
            // Put back with current bounds
            decompQueue.add(leaf);
            return;
        }
        minimizedConfs.add(confKey);

        // Use branch emats so pairwise scoring stays inside the sparse pfunc world.
        // Cheap emat lookups: keep on main thread; capture pre-min state for the callback.
        final double eRigid = computeFullConfPairwiseEnergy(fullConf, branchRigidEmat);
        final double eMin = computeFullConfPairwiseEnergy(fullConf, branchMinimizingEmat);
        final double oldConfLower = leaf.confLowerBound;
        final double oldConfUpper = leaf.confUpperBound;
        final double epsilonBefore = computeEpsilon();
        final BigDecimal oldLeafLower = leaf.subtreeLowerBound;
        final BigDecimal oldLeafUpper = leaf.subtreeUpperBound;
        final ConfSearch.ScoredConf scoredConf = new ConfSearch.ScoredConf(fullConf, eMin);

        // Submit the CCD minimization (and any triple-correction work that needs the
        // ConfAnalysis object) to loopTasks. Multiple leaves run in parallel up to
        // parallelism.numThreads. All shared-state updates happen inside synchronized
        // blocks so flatSumZ / totalMinimizations / decompQueue stay consistent.
        loopTasks.submit(() -> {
            // Worker thread: run CCD and (in FULL mode) compute triple corrections.
            double eTrueLocal;
            if (useSparsePfunc) {
                eTrueLocal = minimizingEcalc.calcEnergy(new RCTuple(fullConf), makeSparseFullConfInters(fullConf)).energy;
                if (useHigherOrderCorrections) {
                    synchronized (this) {
                        correctionSparseGenerationSkips++;
                        if (shouldPrintCorrectionDetail(correctionSparseGenerationSkips)) {
                            correctionAuditLog("NO-GENERATE sparse leaf min#" + (totalMinimizations + 1)
                                    + " state=" + stateName
                                    + " reason=SPARSE uses sparse interactions and does not run ConfAnalyzer triple corrections"
                                    + " conf=" + SimpleConfSpace.formatConfRCs(fullConf));
                        }
                    }
                }
            } else {
                ConfAnalyzer.ConfAnalysis analysis = new ConfAnalyzer(minimizingEcalc).analyze(scoredConf);
                eTrueLocal = analysis.epmol.energy;
                if (useHigherOrderCorrections) {
                    // computeTripleCorrections mutates shared correction state; serialize.
                    synchronized (this) {
                        computeTripleCorrections(analysis, scoredConf);
                    }
                }
            }
            return eTrueLocal;
        }, (Double eTrueObj) -> {
            // Listener: update shared bookkeeping under a single lock so flatSumZ,
            // totalMinimizations, and per-leaf profile recording stay atomic.
            double eTrue = eTrueObj.doubleValue();
            synchronized (this) {
                totalMinimizations++;

                // Update Z bounds: remove old contribution, add new (exact) contribution.
                flatSumZLower = flatSumZLower.subtract(oldLeafLower);
                flatSumZUpper = flatSumZUpper.subtract(oldLeafUpper);

                BigDecimal boltzTrue = bc.calc(eTrue);
                flatSumZLower = flatSumZLower.add(boltzTrue);
                flatSumZUpper = flatSumZUpper.add(boltzTrue);
                double epsilonAfter = computeEpsilon();

                // Mark as minimized (don't re-add to queue — it's exact now).
                leaf.minimized = true;
                leaf.minimizedEnergy = eTrue;

                recordLeafMinimizationProfile(BranchDpConfig.getBackendName() + "-" + edgeSelectionStrategy.name(),
                        totalMinimizations, fullConf,
                        oldConfLower, oldConfUpper, eTrue,
                        oldLeafLower, oldLeafUpper, boltzTrue,
                        epsilonBefore, epsilonAfter);

                if (totalMinimizations % 100 == 0 || totalMinimizations <= 5) {
                    System.out.println(BranchDpConfig.getBackendLogPrefix() + " " + totalMinimizations
                            + " minimizations, epsilon=" + String.format("%.6f", epsilonBound)
                            + ", flatSumZLower=" + String.format("%.6e", flatSumZLower.doubleValue())
                            + ", flatSumZUpper=" + String.format("%.6e", flatSumZUpper.doubleValue())
                            + ", oldConf=[" + String.format("%.4f", oldConfLower)
                            + "," + String.format("%.4f", oldConfUpper) + "]"
                            + ", eRigid=" + String.format("%.4f", eRigid)
                            + ", eMin=" + String.format("%.4f", eMin)
                            + ", eTrue=" + String.format("%.4f", eTrue)
                            + ", boltzTrue=" + String.format("%.6e", boltzTrue.doubleValue())
                            + ", leafOldLower=" + String.format("%.6e", oldLeafLower.doubleValue())
                            + ", leafOldUpper=" + String.format("%.6e", oldLeafUpper.doubleValue())
                            + (trace ? ", conf=" + SimpleConfSpace.formatConfRCs(fullConf) : ""));
                }
            }
        });
    }

    /**
     * Apply a correction to a leaf node without full minimization.
     * Similar to MARK*'s triple correction skip.
     */
    private void applyCorrection(DecompSearchNode leaf, double correctionEnergy) {
        // Remove old contribution
        flatSumZLower = flatSumZLower.subtract(leaf.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(leaf.subtreeUpperBound);

        // Update lower bound: use correction energy instead of rigid
        BigDecimal boltzCorr = bc.calc(correctionEnergy);
        // Upper bound stays: minimizing emat bound is still valid
        double eMin = computeFullConfPairwiseEnergy(leaf.partialConf, branchMinimizingEmat);
        BigDecimal boltzMin = bc.calc(eMin);

        // New contribution
        leaf.subtreeLowerBound = boltzCorr;
        leaf.subtreeUpperBound = boltzMin;
        leaf.errorBound = leaf.subtreeUpperBound.subtract(leaf.subtreeLowerBound);
        if (leaf.errorBound.signum() < 0) leaf.errorBound = BigDecimal.ZERO;

        flatSumZLower = flatSumZLower.add(leaf.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.add(leaf.subtreeUpperBound);

        // Put back in queue with updated bounds
        if (leaf.errorBound.signum() > 0) {
            decompQueue.add(leaf);
        }
    }

    /**
     * Compute triple corrections from a minimized conformation.
     */
    private void computeTripleCorrections(ConfAnalyzer.ConfAnalysis analysis,
                                            ConfSearch.ScoredConf conf) {
        if (conf.getAssignments().length < 3) return;
        correctionGenerationCalls++;

        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(branchMinimizingEmat);
        EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);

        List<int[]> largePairs = new ArrayList<>();
        double minDifference = 0.9;
        for (int pos1 = 0; pos1 < diff.getNumPos(); pos1++) {
            for (int pos2 = pos1 + 1; pos2 < diff.getNumPos(); pos2++) {
                double sum = diff.getOneBody(pos1, 0)
                        + diff.getPairwise(pos1, 0, pos2, 0)
                        + diff.getOneBody(pos2, 0);
                if (sum >= minDifference) {
                    largePairs.add(new int[]{pos1, pos2});
                }
            }
        }
        correctionLargePairs += largePairs.size();
        if (shouldPrintCorrectionDetail(correctionGenerationCalls)) {
            correctionAuditLog("GENERATE call#" + correctionGenerationCalls
                    + " state=" + stateName
                    + " conf=" + SimpleConfSpace.formatConfRCs(conf.getAssignments())
                    + " largePairs=" + largePairs.size());
        }

        double tripleThreshold = 0.3;
        for (int[] pair : largePairs) {
            int pos1 = pair[0], pos2 = pair[1];
            for (int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos1 || pos3 == pos2) continue;
                correctionTupleCandidates++;
                RCTuple tuple = new RCTuple(
                        pos1, conf.getAssignments()[pos1],
                        pos2, conf.getAssignments()[pos2],
                        pos3, conf.getAssignments()[pos3]);
                double tupleBounds = branchRigidEmat.getInternalEnergy(tuple)
                        - branchMinimizingEmat.getInternalEnergy(tuple);
                if (tupleBounds < tripleThreshold) {
                    correctionTupleBelowThreshold++;
                    continue;
                }
                if (correctionMatrix.hasHigherOrderTermFor(tuple)) {
                    correctionTupleAlreadyKnown++;
                    continue;
                }

                correctionAsyncSubmitted++;
                minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
                    double tripleEnergy = minimizedTuple.energy;
                    double lowerbound = branchMinimizingEmat.getInternalEnergy(tuple);
                    double correction = tripleEnergy - lowerbound;
                    if (correction != 0) {
                        synchronized (BranchDpBackend.this) {
                            correctionMatrix.setHigherOrder(tuple, correction);
                            correctionGenerated++;
                            if (shouldPrintCorrectionDetail(correctionGenerated)) {
                                correctionAuditLog("GENERATED correction#" + correctionGenerated
                                        + " state=" + stateName
                                        + " tuple=" + tuple.stringListing()
                                        + " lower=" + String.format("%.4f", lowerbound)
                                        + " triple=" + String.format("%.4f", tripleEnergy)
                                        + " correction=" + String.format("%.4f", correction));
                            }
                        }
                    } else {
                        synchronized (BranchDpBackend.this) {
                            correctionZeroGenerated++;
                        }
                    }
                });
            }
        }
        minimizingEcalc.tasks.waitForFinish();
    }

    // ========== Shared utilities ==========


    private double computeFullConfPairwiseEnergy(int[] conf, EnergyMatrix emat) {
        double energy = emat.getConstTerm();
        int numPos = conf.length;
        for (int i = 0; i < numPos; i++) {
            if (conf[i] < 0) continue;
            energy += emat.getOneBody(i, conf[i]);
            for (int j = i + 1; j < numPos; j++) {
                if (conf[j] < 0) continue;
                energy += emat.getPairwise(i, conf[i], j, conf[j]);
            }
        }
        return energy;
    }

    private void printEdgeSelectionStats() {
        BigDecimal avgSelectedDropPerState = BigDecimal.ZERO;
        if (edgeLookaheadSelectedStates > 0) {
            avgSelectedDropPerState = MathTools.bigDivide(
                    edgeLookaheadSelectedDrop,
                    BigDecimal.valueOf(edgeLookaheadSelectedStates),
                    PartitionFunction.decimalPrecision);
        }

        int cachedEntries = 0;
        synchronized (edgeLookaheadCache) {
            for (Map<Integer, EdgeLookaheadResult> nodeCache : edgeLookaheadCache.values()) {
                cachedEntries += nodeCache.size();
            }
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Edge selection stats: "
                + "calls=" + edgeSelectionCalls
                + ", rootSplit=" + selectedRootSplitIndex
                + ", lambdaStates=" + edgeSelectionLambdaStatesCalls
                + ", singleEdge=" + edgeSelectionSingleEdgeCalls
                + ", multiEdge=" + edgeSelectionMultiEdgeCalls
                + ", maxPendingEdges=" + maxPendingEdgesSeen
                + ", pendingHistogram=" + pendingEdgeCountHistogram
                + ", exactLookahead=" + edgeSelectionExactLookaheadCalls
                + ", fallback=" + edgeSelectionFallbackCalls
                + ", fallbackMaxStates=" + edgeSelectionFallbackMaxStates
                + ", candidates=" + edgeLookaheadCandidates
                + ", enumeratedStates=" + edgeLookaheadEnumeratedStates
                + ", cacheHits=" + edgeLookaheadCacheHits
                + ", cacheMisses=" + edgeLookaheadCacheMisses
                + ", cacheEntries=" + cachedEntries
                + ", parallelCalls=" + edgeLookaheadParallelCalls
                + ", lookaheadTimeMs=" + String.format("%.1f", edgeLookaheadNanos / 1e6)
                + ", selectedDrop=" + formatBigExp(edgeLookaheadSelectedDrop)
                + ", avgSelectedDropPerState=" + formatBigExp(avgSelectedDropPerState));
    }

    private void printCorrectionAuditStats() {
        if (!useHigherOrderCorrections && !correctionAudit
                && correctionChecks == 0 && correctionGenerationCalls == 0
                && correctionSparseGenerationSkips == 0) {
            return;
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Correction audit summary: "
                + "requested=" + useHigherOrderCorrections
                + ", applyEnabled=" + (useHigherOrderCorrections && !useSparsePfunc)
                + ", mode=" + energyMode
                + ", checks=" + correctionChecks
                + ", applied=" + correctionApplications
                + ", notTighter=" + correctionNotTighter
                + ", sparseGenerationSkips=" + correctionSparseGenerationSkips
                + ", generationCalls=" + correctionGenerationCalls
                + ", largePairs=" + correctionLargePairs
                + ", tupleCandidates=" + correctionTupleCandidates
                + ", belowThreshold=" + correctionTupleBelowThreshold
                + ", alreadyKnown=" + correctionTupleAlreadyKnown
                + ", asyncSubmitted=" + correctionAsyncSubmitted
                + ", generated=" + correctionGenerated
                + ", zeroGenerated=" + correctionZeroGenerated);
        if (useHigherOrderCorrections && useSparsePfunc) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Correction audit note: SPARSE mode skips correctionMatrix applications "
                    + "and triple-correction generation; internal node bounds come from branch DP tables.");
        }
    }

    // ========== Override compute ==========

    @Override
    protected void updateBound() {
        if (useBranchDecomposition) {
            epsilonBound = computeEpsilon();
            auditFlatSumNumerics("update", false);
        } else {
            super.updateBound();
        }
    }

    @Override
    public void compute(int maxNumConfs) {
        if (dryRun) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " DRYRUN abort (no DP/CCD computed); rootSplit="
                    + rootSplitStrategy);
            values = new Values();
            values.qprime = MathTools.BigPositiveInfinity;
            setStatus(Status.Aborted);
            return;
        }
        if (dpTooLargeException != null) {
            reportDPTooLarge(dpTooLargeException);
            values = new Values();
            values.qprime = MathTools.BigPositiveInfinity;
            setStatus(Status.Aborted);
            return;
        }

        if (!useBranchDecomposition) {
            computeWithoutBranchDecomposition(maxNumConfs);
            return;
        }

        roundCounter = 0;

        if (beforeBranchDecompositionSearch(maxNumConfs)) {
            return;
        }

        double lastEps = 1.0;

        while (epsilonBound > targetEpsilon) {
            tightenBoundInPhases();

            if (lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println(BranchDpConfig.getBackendLogPrefix() + " Warning - bounds got looser. eps=" + epsilonBound);
            }
            lastEps = epsilonBound;
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Finished. epsilon=" + String.format("%.6f", epsilonBound)
                + " after " + totalEnumerationSteps + " enum steps, "
                + totalMinimizations + " minimizations, "
                + totalDrillDowns + " drill-downs, "
                + totalDrillDownEarlyStops + " early stops.");
        printLeafMinimizationProfile();
        printEdgeSelectionStats();
        printCorrectionAuditStats();
        auditFlatSumNumerics("final", true);
        edgeLookaheadCache.clear();

        // Set final partition function values
        PartitionFunction.Values vals = getValues();
        double reportedEpsilon = setReportedPartitionFunctionValues(vals);

        if (certifyFullBounds && reportedEpsilon > targetEpsilon) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Sparse-to-full certificate epsilon="
                    + String.format("%.6f", reportedEpsilon)
                    + " > target=" + String.format("%.6f", targetEpsilon)
                    + "; leaving status as " + getStatus() + ".");
        }

        if (reportedEpsilon <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
    }

    protected void computeWithoutBranchDecomposition(int maxNumConfs) {
        super.compute(maxNumConfs);
    }

    protected boolean beforeBranchDecompositionSearch(int maxNumConfs) {
        return false;
    }

    /**
     * Deterministically detach the rooted DP/search structures after a caller
     * has copied out the scalar partition-function result.
     *
     * <p>This is deliberately separate from {@link #compute(int)}: callers may
     * inspect values and status after computation, but none of those values
     * depend on the rooted tree.  Production PACK* tables can occupy hundreds
     * of GiB, so relying on a best-effort {@code System.gc()} while the facade
     * still references this backend can exhaust the heap across K* sequences.</p>
     */
    public synchronized void releaseLargeMemory() {
        if (largeMemoryReleased) {
            return;
        }
        largeMemoryReleased = true;

        edgeLookaheadCache.clear();
        minimizedConfs.clear();
        if (decompQueue != null) {
            decompQueue.clear();
            decompQueue = null;
        }
        queue.clear();
        rootNode = null;

        if (rootedRoot != null) {
            RootedTreeEdge.postOrderReleaseLargeMemory(rootedRoot);
        }
        rootedRootEdge = null;
        rootedRoot = null;
        dpTablesReady = false;

        gScorerMin = null;
        hScorerMin = null;
        gScorerRigid = null;
        hScorerNegRigid = null;
        searchRCs = null;
        branchRigidEmat = null;
        branchMinimizingEmat = null;
    }

    /**
     * Finish a pfunc using the current deterministic DP/search bounds.
     */
    protected void finishWithCurrentDPBounds(String message) {
        System.out.println(message);

        PartitionFunction.Values vals = getValues();
        double reportedEpsilon = setReportedPartitionFunctionValues(vals);
        if (reportedEpsilon <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
    }

    protected void restoreInitialDPTables(String logPrefix) {
        long restoreStart = System.currentTimeMillis();
        double rtForDP = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
        RootedTreeEdge.postOrderInitIncremental(rootedRoot,
                branchRigidEmat, branchMinimizingEmat, interactionGraph, rtForDP);
        DPCacheStats restoreStats = computeFullDPTables("initial", true);
        System.out.println(logPrefix + " restored original DP tables in "
                + (System.currentTimeMillis() - restoreStart) + " ms"
                + formatDPCacheStats(restoreStats));
    }

    private double setReportedPartitionFunctionValues(PartitionFunction.Values vals) {
        if (!certifyFullBounds) {
            vals.qstar = flatSumZLower;
            vals.pstar = flatSumZUpper;
            vals.qprime = vals.pstar.subtract(vals.qstar);
            trace("reported values state=" + stateName
                    + " certifyFullBounds=false"
                    + " qstar=" + formatBigExp(vals.qstar)
                    + " pstar=" + formatBigExp(vals.pstar)
                    + " qprime=" + formatBigExp(vals.qprime)
                    + " eps=" + String.format("%.6f", epsilonBound)
                    + " corrections=" + correctionApplications);
            return epsilonBound;
        }

        double rho = interactionGraph.getCutResidualUpperBound();
        if (rho <= 0.0) {
            vals.qstar = flatSumZLower;
            vals.pstar = flatSumZUpper;
            vals.qprime = vals.pstar.subtract(vals.qstar);
            trace("reported values state=" + stateName
                    + " rho=0"
                    + " qstar=" + formatBigExp(vals.qstar)
                    + " pstar=" + formatBigExp(vals.pstar)
                    + " qprime=" + formatBigExp(vals.qprime)
                    + " eps=" + String.format("%.6f", epsilonBound)
                    + " corrections=" + correctionApplications);
            return epsilonBound;
        }

        if (!Double.isFinite(rho)) {
            vals.qstar = BigDecimal.ZERO;
            vals.pstar = MathTools.BigPositiveInfinity;
            vals.qprime = MathTools.BigPositiveInfinity;
            return 1.0;
        }

        double exponent = rho / (bc.R * bc.T);
        BigDecimal lowerFactor = bc.exp(-exponent);
        BigDecimal upperFactor = bc.exp(exponent);

        vals.qstar = flatSumZLower.multiply(lowerFactor, PartitionFunction.decimalPrecision);
        vals.pstar = flatSumZUpper.multiply(upperFactor, PartitionFunction.decimalPrecision);
        vals.qprime = vals.pstar.subtract(vals.qstar);

        double certifiedEpsilon = MathTools.bigDivide(
                vals.pstar.subtract(vals.qstar),
                vals.pstar,
                PartitionFunction.decimalPrecision).doubleValue();

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Applied sparse-to-full pfunc inflation: rho="
                + String.format("%.6f", rho)
                + ", lowerFactor=" + String.format("%.6e", lowerFactor.doubleValue())
                + ", upperFactor=" + String.format("%.6e", upperFactor.doubleValue())
                + ", sparseEpsilon=" + String.format("%.6f", epsilonBound)
                + ", certifiedEpsilon=" + String.format("%.6f", certifiedEpsilon));

        return certifiedEpsilon;
    }
}
