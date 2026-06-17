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
    private static final String CUTOFF_STRATEGY_PROPERTY = "branchmarkstar.cutoff.strategy";
    private static final String DIST_CUTOFF_PROPERTY = "branchmarkstar.cutoff.dist";
    private static final String ENERGY_CUTOFF_PROPERTY = "branchmarkstar.cutoff.energy";
    private static final String RESIDUAL_BUDGET_PROPERTY = "branchmarkstar.cutoff.residualBudget";
    private static final String CUTOFF_KEEP_CONNECTED_PROPERTY = "branchmarkstar.cutoff.keepConnected";
    private static final String CERTIFY_FULL_BOUNDS_PROPERTY = "branchmarkstar.certifyFullBounds";
    private static final String EDGE_SELECTION_PROPERTY = "branchmarkstar.edgeSelection";
    private static final String EDGE_LOOKAHEAD_MAX_STATES_PROPERTY = "branchmarkstar.edgeSelection.maxExactStates";
    private static final String EDGE_LOOKAHEAD_MAX_PENDING_EDGES_PROPERTY = "branchmarkstar.edgeSelection.maxPendingEdges";
    private static final String EDGE_LOOKAHEAD_PARALLEL_PROPERTY = "branchmarkstar.edgeSelection.parallelLookahead";
    private static final String ROOT_SPLIT_PROPERTY = "branchmarkstar.rootSplit";
    // 0-edge (no pairwise) graphs: exact independent-position DP. Default off.
    private static final String ZERO_EDGE_DIRECT_PROPERTY = "branchmarkstar.dp.zeroEdgeDirect";
    private static final String ROOT_SPLIT_MAX_FSET_PROPERTY = "branchmarkstar.rootSplit.maxFset";
    private static final String DRY_RUN_PROPERTY = "branchmarkstar.dp.dryRun";
    private static final String MUTABLE_POSITIONS_PROPERTY = "branchmarkstar.mutablePositions";
    private static final String PARALLEL_INTERNAL_PROPERTY = "branchmarkstar.parallel.internal";
    private static final String PARALLEL_ENUMERATION_PROPERTY = "branchmarkstar.parallel.enumeration";
    private static final String NUMERIC_AUDIT_PROPERTY = "branchmarkstar.numericAudit";
    private static final String NUMERIC_AUDIT_INTERVAL_PROPERTY = "branchmarkstar.numericAudit.interval";
    private static final String TRACE_PROPERTY = "branchmarkstar.trace";
    private static final String TRACE_ROUNDS_PROPERTY = "branchmarkstar.trace.rounds";
    private static final String TRACE_CORRECTIONS_PROPERTY = "branchmarkstar.trace.corrections";
    private static final String CORRECTION_AUDIT_PROPERTY = "branchmarkstar.correctionAudit";
    private static final String ENERGY_MODE_PROPERTY = "branchmarkstar.energyMode";
    private static final String SPARSE_PRUNE_THRESHOLD_PROPERTY = "branchmarkstar.sparse.pruneThreshold";
    private static final String REGION_ATOM_ENABLED_PROPERTY = "branchmarkstar.regionAtom.enabled";
    private static final String REGION_ATOM_REGIONS_PROPERTY = "branchmarkstar.regionAtom.regions";
    private static final String REGION_ATOM_SCOUT_ONLY_PROPERTY = "branchmarkstar.regionAtom.scoutOnly";
    private static final String REGION_ATOM_SCOUT_MAX_LEAVES_PROPERTY = "branchmarkstar.regionAtom.scout.maxLeaves";
    private static final String REGION_ATOM_SCOUT_TOP_CELLS_PROPERTY = "branchmarkstar.regionAtom.scout.topCells";
    private static final String REGION_ATOM_TABLE_MAX_JOBS_PROPERTY = "branchmarkstar.regionAtom.table.maxJobs";
    private static final String REGION_ATOM_WHATIF_DELTAS_PROPERTY = "branchmarkstar.regionAtom.whatIfDeltas";
    private static final String REGION_ATOM_CERTIFY_PROPERTY = "branchmarkstar.regionAtom.certify";
    private static final String REGION_ATOM_CERTIFY_USE_DP_PROPERTY = "branchmarkstar.regionAtom.certify.useDP";
    /** Top-K boundary cells (by Phase 0 emat gap) to certify. 0 = all cells. */
    private static final String REGION_ATOM_CERTIFY_TOPK_PROPERTY = "branchmarkstar.regionAtom.certify.topK";
    private static final String DP_CACHE_ENABLED_PROPERTY = "branchmarkstar.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "branchmarkstar.dp.cache.maxEntries";
    private static final String DP_CACHE_MAX_TABLE_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTableBytes";
    private static final String DP_CACHE_MAX_TOTAL_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTotalBytes";
    private static final String DP_CACHE_SKIP_IF_M_STATES_PROPERTY = "branchmarkstar.dp.cache.skipIfMStates";
    private static final int DEFAULT_DP_CACHE_MAX_ENTRIES = 20000;
    private static final long DEFAULT_DP_CACHE_MAX_TABLE_BYTES = 256L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_MAX_TOTAL_BYTES = 4L * 1024L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_SKIP_IF_M_STATES = 8_000_000L;

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
    private final boolean regionAtomEnabled;
    private final boolean regionAtomScoutOnly;
    private final int regionAtomScoutMaxLeaves;
    private final int regionAtomScoutTopCells;
    private final long regionAtomTableMaxJobs;
    private final double[] regionAtomWhatIfDeltas;
    private final Set<String> regionAtomCertifyNames;
    private final boolean regionAtomCertifyUseDP;
    private final int regionAtomCertifyTopK;
    private final boolean dpCacheEnabled;
    private final int dpCacheMaxEntries;
    private final long dpCacheMaxTableBytes;
    private final long dpCacheMaxTotalBytes;
    private final long dpCacheSkipIfMStates;
    private final boolean dryRun;
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

    // ========== Region-atom table scout (diagnostic, pre-minimization) ==========

    private final List<RegionAtomTable> regionAtomTables = new ArrayList<>();
    private final Set<String> regionAtomSeenLeafConfs = new HashSet<>();
    private boolean regionAtomTablesInitialized = false;
    private boolean regionAtomScoutComplete = false;
    private boolean regionAtomScoutSummaryPrinted = false;
    private long regionAtomScoutLeavesSeen = 0;
    private BigDecimal regionAtomScoutLower = BigDecimal.ZERO;
    private BigDecimal regionAtomScoutUpper = BigDecimal.ZERO;
    private BigDecimal regionAtomScoutGap = BigDecimal.ZERO;

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
        final double logDPWork;

        RootingCandidate(RootedTreeNode root, int splitEdgeIndex, double logTESS,
                         int lambdaEdges, int maxFsetSize, int branchingEdges,
                         int totalFsetEdges, int rootFsetSize,
                         int reusableLambdaEdges, double reusableLogWork,
                         double totalLogWork, long maxMStates, long totalMStates,
                         long maxDPTableBytes, long totalDPTableBytes,
                         double logDPWork) {
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
            this.logDPWork = logDPWork;
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
        this.edgeSelectionStrategy = getEdgeSelectionStrategy();
        this.edgeLookaheadMaxStates = getEdgeLookaheadMaxStates();
        this.edgeLookaheadMaxPendingEdges = getEdgeLookaheadMaxPendingEdges();
        this.edgeLookaheadParallel = getEdgeLookaheadParallel();
        this.rootSplitStrategy = getConfigProperty(ROOT_SPLIT_PROPERTY, "work").trim();
        this.dryRun = getConfigBoolean(DRY_RUN_PROPERTY, false);
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
        this.useHigherOrderCorrections = getConfigBoolean("branchmarkstar.useHigherOrderCorrections", false);
        this.regionAtomEnabled = getRegionAtomEnabled();
        this.regionAtomScoutOnly = getRegionAtomScoutOnly();
        this.regionAtomScoutMaxLeaves = Math.max(1,
                getConfigInteger(REGION_ATOM_SCOUT_MAX_LEAVES_PROPERTY, 2000));
        this.regionAtomScoutTopCells = Math.max(1,
                getConfigInteger(REGION_ATOM_SCOUT_TOP_CELLS_PROPERTY, 12));
        this.regionAtomTableMaxJobs = Math.max(1L,
                getConfigLong(REGION_ATOM_TABLE_MAX_JOBS_PROPERTY, 50000000L));
        this.regionAtomWhatIfDeltas = getRegionAtomWhatIfDeltas();
        this.regionAtomCertifyNames = getRegionAtomCertifyNames();
        this.regionAtomCertifyUseDP = getRegionAtomCertifyUseDP();
        this.regionAtomCertifyTopK = getRegionAtomCertifyTopK();
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
        if (regionAtomEnabled) {
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " Region-atom table scout enabled"
                    + " (regions=" + getConfigProperty(REGION_ATOM_REGIONS_PROPERTY, "(auto)")
                    + ", scoutOnly=" + regionAtomScoutOnly
                    + ", maxLeaves=" + regionAtomScoutMaxLeaves
                    + ", topCells=" + regionAtomScoutTopCells
                    + ", tableMaxJobs=" + regionAtomTableMaxJobs + ")");
        }
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
        branchDecomposition = new BranchDecomposition(interactionGraph);
        branchDecomposition.compute();
        branchwidth = branchDecomposition.getBranchwidth();

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Branch decomposition computed. Branchwidth=" + branchwidth
                + ", positions=" + interactionGraph.getNumPositions()
                + ", edges=" + branchDecomposition.getTree().getNumEdges());

        // Step 3: Root the tree. Legacy uses split edge 0; optional strategies
        // can choose a split that exposes multiple pending edges for lookahead.
        RootingCandidate rooting;
        try {
            rooting = selectRooting(rcs);
        } catch (DPTableTooLargeException e) {
            dpTooLargeException = e;
            reportDPTooLarge(e);
            return;
        }
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
                + ", logDPWork=" + String.format(Locale.ROOT, "%.2f", rooting.logDPWork)
                + ")");

        rootedRootEdge = rootedRoot.getLeftChild().getChildOfEdge();

        // Step 4: TESS check
        double logTESS = rootedRootEdge.computeLogTESS();
        double logNaive = computeLogNaive(rcs);
        double tessRatio = Math.exp(logTESS - logNaive);

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Compact tree TESS: logTESS=" + String.format("%.2f", logTESS)
                + ", logNaive=" + String.format("%.2f", logNaive)
                + ", ratio=" + String.format("%.4f", tessRatio));
        logDPMemoryPredictions("initial", logTESS);

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
                    + " logDPWork=" + String.format(Locale.ROOT, "%.4f", rooting.logDPWork)
                    + " logTESS=" + String.format(Locale.ROOT, "%.4f", logTESS)
                    + " logNaive=" + String.format(Locale.ROOT, "%.4f", logNaive));
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

        // Initialize region-atom diagnostic tables (Phase 0 emat-only scout).
        // Runs only when branch decomposition is active and the user opted in
        // with branchmarkstar.regionAtom.enabled=true. Diagnostic only — does
        // NOT modify flatSumZ, DP tables, or correctionMatrix.
        if (regionAtomEnabled) {
            initRegionAtoms(rcs);

            // DP-level integration: if any tables were certified, apply corrections
            // to DP table entries and recompute. Must happen AFTER certification
            // (initRegionAtoms) and AFTER initial DP build (initSearch).
            if (regionAtomCertifyUseDP && dpTablesReady
                    && regionAtomTables.stream().anyMatch(t -> t.certified)) {
                double rtForDP = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
                applyRegionAtomToDP(rtForDP);
                RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
                // Rebuild root search node with corrected DP bounds
                DecompSearchNode rootSearchNode = DecompSearchNode.makeRoot(
                        rootedRootEdge, rcs.getNumPos(), rcs,
                        gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);
                enhanceWithDPBounds(rootSearchNode);
                decompQueue.clear();
                decompQueue.add(rootSearchNode);
                flatSumZLower = rootSearchNode.subtreeLowerBound;
                flatSumZUpper = rootSearchNode.subtreeUpperBound;
                updateBound();
                System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP tables recomputed after region-atom certification."
                        + " Root upper=" + String.format("%.4e", flatSumZUpper.doubleValue())
                        + ", lower=" + String.format("%.4e", flatSumZLower.doubleValue())
                        + ", epsilon=" + String.format("%.6f", epsilonBound));
            }
        }
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

    protected boolean getRegionAtomEnabled() {
        return getConfigBoolean(REGION_ATOM_ENABLED_PROPERTY, false);
    }

    protected boolean getRegionAtomScoutOnly() {
        return getConfigBoolean(REGION_ATOM_SCOUT_ONLY_PROPERTY, false);
    }

    protected double[] getRegionAtomWhatIfDeltas() {
        return parseDoubleList(getConfigProperty(REGION_ATOM_WHATIF_DELTAS_PROPERTY, "0.1,0.3,0.5"));
    }

    protected Set<String> getRegionAtomCertifyNames() {
        return parseRegionAtomCertifyNames(getConfigProperty(REGION_ATOM_CERTIFY_PROPERTY, ""));
    }

    protected boolean getRegionAtomCertifyUseDP() {
        return getConfigBoolean(REGION_ATOM_CERTIFY_USE_DP_PROPERTY, false);
    }

    protected int getRegionAtomCertifyTopK() {
        return Math.max(0, getConfigInteger(REGION_ATOM_CERTIFY_TOPK_PROPERTY, 0));
    }

    protected void logBackendControlOverrides() {
    }

    protected ConfEnergyCalculator getMinimizingEcalc() {
        return minimizingEcalc;
    }

    protected String getConfigProperty(String key, String defaultValue) {
        String value = System.getProperty(key);
        if (value == null) {
            value = System.getProperty("osprey." + key);
        }
        return value != null ? value : defaultValue;
    }

    protected int getConfigInteger(String key, int defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value.trim());
        } catch (NumberFormatException e) {
            warnConfig("Invalid integer", key, value, Integer.toString(defaultValue));
            return defaultValue;
        }
    }

    protected long getConfigLong(String key, long defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Long.parseLong(value.trim());
        } catch (NumberFormatException e) {
            warnConfig("Invalid long", key, value, Long.toString(defaultValue));
            return defaultValue;
        }
    }

    protected long getConfigBytes(String key, long defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return parseByteCount(value.trim());
        } catch (NumberFormatException e) {
            warnConfig("Invalid byte count", key, value, Long.toString(defaultValue));
            return defaultValue;
        }
    }

    protected double getConfigDouble(String key, double defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            warnConfig("Invalid double", key, value, Double.toString(defaultValue));
            return defaultValue;
        }
    }

    protected boolean getConfigBoolean(String key, boolean defaultValue) {
        String value = getConfigProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    private static long parseByteCount(String value) {
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace("_", "");
        long multiplier = 1L;
        if (normalized.endsWith("kib")) {
            multiplier = 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("kb") || normalized.endsWith("k")) {
            multiplier = 1024L;
            normalized = normalized.replaceAll("kb?$", "");
        } else if (normalized.endsWith("mib")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("mb") || normalized.endsWith("m")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.replaceAll("mb?$", "");
        } else if (normalized.endsWith("gib")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("gb") || normalized.endsWith("g")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("gb?$", "");
        } else if (normalized.endsWith("tib")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("tb") || normalized.endsWith("t")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("tb?$", "");
        } else if (normalized.endsWith("b")) {
            normalized = normalized.substring(0, normalized.length() - 1);
        }

        double amount = Double.parseDouble(normalized.trim());
        if (amount < 0 || amount > Long.MAX_VALUE / (double) multiplier) {
            throw new NumberFormatException(value);
        }
        return (long) (amount * multiplier);
    }

    private static void warnConfig(String kind, String key, String value, String defaultValue) {
        System.err.println(BranchDpConfig.getBackendLogPrefix() + " " + kind + " for '" + key
                + "': '" + value + "', using " + defaultValue + ".");
    }

    private static double[] parseDoubleList(String value) {
        if (value == null || value.trim().isEmpty()) {
            return new double[0];
        }
        String[] fields = value.split(",");
        List<Double> parsed = new ArrayList<>();
        for (String field : fields) {
            String trimmed = field.trim();
            if (trimmed.isEmpty()) continue;
            try {
                parsed.add(Double.parseDouble(trimmed));
            } catch (NumberFormatException e) {
                System.err.println(BranchDpConfig.getBackendLogPrefix() + " Invalid double in list: '" + trimmed + "', skipping.");
            }
        }
        double[] result = new double[parsed.size()];
        for (int i = 0; i < parsed.size(); i++) {
            result[i] = parsed.get(i);
        }
        return result;
    }

    private RootingCandidate selectRooting(RCs rcs) {
        int numEdges = branchDecomposition.getTree().getNumEdges();
        if (numEdges == 0) return null;

        String strategy = rootSplitStrategy.toLowerCase(Locale.ROOT);
        if (strategy.isEmpty() || strategy.equals("auto")) {
            strategy = "memory";
        }
        if (strategy.equals("legacy") || strategy.equals("edge0")) {
            return evaluateRootSplit(rcs, 0, true);
        }

        try {
            int explicitSplit = Integer.parseInt(strategy);
            return evaluateRootSplit(rcs, explicitSplit, true);
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

        if (!strategy.equals("branching") && !strategy.equals("maxfset")
                && !strategy.equals("lookahead") && !useReuseScoring
                && !useMemoryScoring && !useWorkScoring) {
            System.err.println(BranchDpConfig.getBackendLogPrefix() + " Unknown root split strategy '" + rootSplitStrategy
                    + "', using legacy split edge 0.");
            return evaluateRootSplit(rcs, 0, true);
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
        RootingCandidate best = null;
        for (int splitIdx = 0; splitIdx < numEdges; splitIdx++) {
            RootingCandidate candidate = evaluateRootSplit(rcs, splitIdx, false);
            if (candidate == null) continue;
            boolean better;
            if (useMemoryScoring) {
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
        return best != null ? evaluateRootSplit(rcs, best.splitEdgeIndex, true) : evaluateRootSplit(rcs, 0, true);
    }

    private RootingCandidate evaluateRootSplit(RCs rcs, int splitEdgeIndex,
                                               boolean initEnumerationArrays) {
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
        double logDPWork = Double.NEGATIVE_INFINITY;
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
            maxMStates = Math.max(maxMStates, mStates);
            totalMStates = saturatingAdd(totalMStates, mStates);
            maxDPTableBytes = Math.max(maxDPTableBytes, tableBytes);
            totalDPTableBytes = saturatingAdd(totalDPTableBytes, tableBytes);
            logDPWork = RootedTreeEdge.logSumExp(logDPWork,
                    Math.log((double) mStates) + Math.log((double) lambdaStates));
        }

        int rootFsetSize = rootEdge.getFset() == null ? 0 : rootEdge.getFset().size();
        Set<Integer> mutablePositions = identifyMutablePositions(rcs);
        ReuseStats reuseStats = computeReuseStats(lambdaEdges, mutablePositions, rcs);
        return new RootingCandidate(root, splitEdgeIndex, rootEdge.computeLogTESS(),
                lambdaEdges.size(), maxFsetSize, branchingEdges, totalFsetEdges, rootFsetSize,
                reuseStats.reusableEdges, reuseStats.reusableLogWork, reuseStats.totalLogWork,
                maxMStates, totalMStates, maxDPTableBytes, totalDPTableBytes, logDPWork);
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

    private static class RegionAtomSpec {
        final String name;
        final int[] region;
        final int[] boundary;
        final BigInteger boundaryCells;
        final BigInteger localStatesPerBoundary;
        final BigInteger totalJobs;

        RegionAtomSpec(String name, int[] region, int[] boundary, RCs rcs) {
            this.name = name;
            this.region = region;
            this.boundary = boundary;
            this.boundaryCells = countStates(boundary, rcs);
            this.localStatesPerBoundary = countStates(region, rcs);
            this.totalJobs = boundaryCells.multiply(localStatesPerBoundary);
        }
    }

    private static class RegionAtomCell {
        int count = 0;
        BigDecimal lower = BigDecimal.ZERO;
        BigDecimal upper = BigDecimal.ZERO;
        BigDecimal gap = BigDecimal.ZERO;
        double preGapSum = 0.0;
        double preGapMax = Double.NEGATIVE_INFINITY;

        void add(BigDecimal lowerZ, BigDecimal upperZ, double preGap) {
            count++;
            lower = lower.add(lowerZ);
            upper = upper.add(upperZ);
            BigDecimal cellGap = upperZ.subtract(lowerZ);
            if (cellGap.signum() > 0) {
                gap = gap.add(cellGap);
            }
            preGapSum += preGap;
            preGapMax = Math.max(preGapMax, preGap);
        }

        double meanPreGap() {
            return count == 0 ? 0.0 : preGapSum / count;
        }
    }

    /**
     * Per-boundary certified table cell:
     *   bKey -> (L_R(b), U_R(b)) with L_R = sum_{x_R} exp(-E^+_rigid/RT) (pre-min lower bound)
     *           and U_R = sum_{x_R} exp(-E_R_certified(x_R; b)/RT) (CCD-min upper bound)
     *
     * Filled by {@link #certifyRegionAtomTable}. Consumed by
     * {@link #lookupCertifiedBoundary} when the DP factor replacement is active.
     */
    private static class CertifiedBoundaryCell {
        final String bKey;
        final int[] bAssignment;
        final int totalLocalStates;
        BigDecimal lowerZ = BigDecimal.ZERO;
        BigDecimal upperZ = BigDecimal.ZERO;
        int localStatesDone = 0;
        int localStatesPruned = 0;
        double minCertEnergy = Double.POSITIVE_INFINITY;
        double maxCertEnergy = Double.NEGATIVE_INFINITY;
        /** Per local assignment x_R: CCD-min owned energy E_R_cert(x_R; b). Used by
         *  the Phase 4 per-node bound tightening to apply the certified value to the
         *  EXACT (b, x_R) cell, not a summary. */
        final Map<String, Double> certifiedEnergyPerLocal = new java.util.concurrent.ConcurrentHashMap<>();

        CertifiedBoundaryCell(String bKey, int[] bAssignment, int totalLocalStates) {
            this.bKey = bKey;
            this.bAssignment = bAssignment;
            this.totalLocalStates = totalLocalStates;
        }

        synchronized void addLocal(String localKey, BigDecimal lower, BigDecimal upper, double eCert) {
            lowerZ = lowerZ.add(lower);
            upperZ = upperZ.add(upper);
            localStatesDone++;
            if (eCert < minCertEnergy) minCertEnergy = eCert;
            if (eCert > maxCertEnergy) maxCertEnergy = eCert;
            certifiedEnergyPerLocal.put(localKey, eCert);
        }

        synchronized void markPruned() {
            localStatesPruned++;
        }
    }

    private static class RegionAtomTable {
        final RegionAtomSpec spec;
        final Map<String, RegionAtomCell> boundaryCells = new HashMap<>();
        final Map<String, RegionAtomCell> localCells = new HashMap<>();
        BigDecimal localLower = BigDecimal.ZERO;
        BigDecimal localUpper = BigDecimal.ZERO;
        BigDecimal localGap = BigDecimal.ZERO;
        int localBoundaryCellsComputed = 0;
        boolean localTableComputed = false;
        boolean localTableSkipped = false;

        // ===== Phase 3 certified table =====
        final Map<String, CertifiedBoundaryCell> certifiedCells = new java.util.concurrent.ConcurrentHashMap<>();
        boolean certified = false;
        BigDecimal certifiedSumLower = BigDecimal.ZERO;
        BigDecimal certifiedSumUpper = BigDecimal.ZERO;
        BigDecimal certifiedSumGap = BigDecimal.ZERO;

        RegionAtomTable(RegionAtomSpec spec) {
            this.spec = spec;
        }

        void recordLeaf(String boundaryKey, String localKey,
                        BigDecimal lowerZ, BigDecimal upperZ, double preGap) {
            boundaryCells.computeIfAbsent(boundaryKey, k -> new RegionAtomCell())
                    .add(lowerZ, upperZ, preGap);
            localCells.computeIfAbsent(boundaryKey + "|" + localKey, k -> new RegionAtomCell())
                    .add(lowerZ, upperZ, preGap);
        }
    }

    // ========== Region-atom: parse + boundary + Phase 0 emat-only scout ==========

    /**
     * Initialize region-atom diagnostics: parse atom specs from the user property,
     * compute the boundary B = N_G(R) \ R against the active interaction graph,
     * and (when affordable) compute the emat-only pre-min L_R(b) / U_R(b) tables.
     *
     * This is Phase 0 of the scout described in
     * {@code slurm/analysis/2026-05-20-region-atom-table-plan.md}. It DOES NOT
     * modify flatSumZ, DP tables, or correctionMatrix — bounds are unaffected.
     */
    private void initRegionAtoms(RCs rcs) {
        if (regionAtomTablesInitialized) return;
        regionAtomTablesInitialized = true;

        String regionsStr = getConfigProperty(REGION_ATOM_REGIONS_PROPERTY, "");
        if (regionsStr == null || regionsStr.trim().isEmpty()) {
            System.out.println("[REGION_ATOM] init skipped: property "
                    + REGION_ATOM_REGIONS_PROPERTY + " is empty");
            return;
        }

        List<int[]> regions = parseRegionAtomRegions(regionsStr, rcs);
        if (regions.isEmpty()) {
            System.out.println("[REGION_ATOM] init skipped: no valid regions parsed from '"
                    + regionsStr + "'");
            return;
        }

        System.out.println("[REGION_ATOM] init: parsed " + regions.size() + " region(s); "
                + "graph=" + cutoffStrategy + ", energyMode=" + energyMode
                + ", numPos=" + rcs.getNumPos() + ", graphEdges=" + interactionGraph.getNumEdges()
                + ", tableMaxJobs=" + regionAtomTableMaxJobs);

        for (int[] region : regions) {
            int[] boundary = computeBoundary(region, interactionGraph);
            String name = formatPositions(region);
            RegionAtomSpec spec = new RegionAtomSpec(name, region, boundary, rcs);
            RegionAtomTable table = new RegionAtomTable(spec);
            regionAtomTables.add(table);

            System.out.println("[REGION_ATOM_SPEC] name=" + name
                    + " R=" + formatPositions(region)
                    + " B=" + formatPositions(boundary)
                    + " |R|=" + region.length
                    + " |B|=" + boundary.length
                    + " boundaryCells=" + spec.boundaryCells
                    + " localStates=" + spec.localStatesPerBoundary
                    + " totalJobs=" + spec.totalJobs);

            BigInteger maxJobs = BigInteger.valueOf(regionAtomTableMaxJobs);
            if (spec.totalJobs.compareTo(maxJobs) > 0) {
                table.localTableSkipped = true;
                System.out.println("[REGION_ATOM_TABLE] name=" + name
                        + " status=skipped reason=totalJobs>" + regionAtomTableMaxJobs);
                continue;
            }

            try {
                computeEmatOnlyTable(table, rcs);
            } catch (RuntimeException e) {
                table.localTableSkipped = true;
                System.err.println("[REGION_ATOM_TABLE] name=" + name
                        + " status=skipped reason=exception:" + e.getMessage());
            }

            // Phase 3: certified table for atoms named in branchmarkstar.regionAtom.certify.
            // Runs the deterministic local CCD oracle to tighten U_R(b). Requires
            // minimizingEcalc to be non-null for current branch-DP backends.
            if (regionAtomCertifyNames.contains(name)) {
                if (table.localTableSkipped) {
                    System.out.println("[REGION_ATOM_CERT] name=" + name
                            + " status=skipped reason=emat-table skipped (totalJobs too large)");
                } else {
                    try {
                        certifyRegionAtomTable(table, rcs);
                    } catch (RuntimeException e) {
                        System.err.println("[REGION_ATOM_CERT] name=" + name
                                + " status=failed reason=exception:" + e.getMessage());
                        e.printStackTrace(System.err);
                    }
                }
            }
        }

        // Phase 4 hook (regionAtomCertifyUseDP) is wired into enhanceWithDPBounds
        // callsites via tightenNodeWithCertifiedTables, so every search node that
        // covers an atom's R+B gets its upper bound tightened. No one-shot init
        // call is needed here — the hook runs lazily as nodes are created.
    }

    /** Parse "7;6,7,8;0,2,6,7,8"-style property into a list of sorted+deduped position arrays. */
    private List<int[]> parseRegionAtomRegions(String spec, RCs rcs) {
        List<int[]> out = new ArrayList<>();
        int numPos = rcs.getNumPos();
        String[] regionTokens = spec.split(";");
        for (String token : regionTokens) {
            String tok = token.trim();
            if (tok.isEmpty()) continue;
            String[] posTokens = tok.split(",");
            TreeSet<Integer> set = new TreeSet<>();
            boolean ok = true;
            for (String pt : posTokens) {
                String p = pt.trim();
                if (p.isEmpty()) continue;
                try {
                    int pos = Integer.parseInt(p);
                    if (pos < 0 || pos >= numPos) {
                        System.err.println("[REGION_ATOM] WARN region '" + tok
                                + "' contains out-of-range position " + pos
                                + " (numPos=" + numPos + "); skipping region");
                        ok = false;
                        break;
                    }
                    set.add(pos);
                } catch (NumberFormatException e) {
                    System.err.println("[REGION_ATOM] WARN region '" + tok
                            + "' contains non-integer token '" + p + "'; skipping region");
                    ok = false;
                    break;
                }
            }
            if (!ok || set.isEmpty()) continue;
            int[] arr = new int[set.size()];
            int i = 0;
            for (Integer v : set) arr[i++] = v;
            out.add(arr);
        }
        return out;
    }

    /** B = N_G(R) \ R against the active interaction graph (cut-filtered in SPARSE mode). */
    private int[] computeBoundary(int[] region, InteractionGraph graph) {
        Set<Integer> regionSet = new HashSet<>();
        for (int r : region) regionSet.add(r);
        TreeSet<Integer> boundary = new TreeSet<>();
        for (int r : region) {
            Set<Integer> neighbors = graph.getNeighbors(r);
            if (neighbors == null) continue;
            for (int n : neighbors) {
                if (!regionSet.contains(n)) boundary.add(n);
            }
        }
        int[] arr = new int[boundary.size()];
        int i = 0;
        for (Integer v : boundary) arr[i++] = v;
        return arr;
    }

    private static String formatPositions(int[] positions) {
        if (positions.length == 0) return "[]";
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < positions.length; i++) {
            if (i > 0) sb.append(",");
            sb.append(positions[i]);
        }
        sb.append("]");
        return sb.toString();
    }

    private static String makeAssignmentKey(int[] positions, int[] rcAtPos) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < positions.length; i++) {
            if (i > 0) sb.append(',');
            sb.append(rcAtPos[i]);
        }
        return sb.toString();
    }

    /**
     * Owned local energy for a region atom under boundary assignment b and
     * local assignment x_R. Owns ONLY:
     *   - one-body terms in R
     *   - R-R pair terms
     *   - R-B crossing pair terms
     * Excludes: one-body in B, B-B pairs, outside-only, B-to-outside terms,
     * and emat constant term.
     */
    private static double ownedLocalEnergy(EnergyMatrix emat,
                                           int[] region, int[] boundary,
                                           int[] xR, int[] xB) {
        double e = 0.0;
        for (int i = 0; i < region.length; i++) {
            int posI = region[i];
            int rcI = xR[i];
            e += emat.getOneBody(posI, rcI);
            // R-R
            for (int j = i + 1; j < region.length; j++) {
                int posJ = region[j];
                int rcJ = xR[j];
                e += emat.getPairwise(posI, rcI, posJ, rcJ);
            }
            // R-B crossings
            for (int j = 0; j < boundary.length; j++) {
                int posJ = boundary[j];
                int rcJ = xB[j];
                e += emat.getPairwise(posI, rcI, posJ, rcJ);
            }
        }
        return e;
    }

    /**
     * Phase 0 emat-only pre-min table: for each boundary assignment b enumerate
     * local x_R and accumulate
     *   L_R^pre(b) = sum_{x_R} exp(-E^+_R(x_R; b) / RT)   [rigid emat]
     *   U_R^pre(b) = sum_{x_R} exp(-E^-_R(x_R; b) / RT)   [minimizing emat]
     * Caches per-boundary aggregates in table.boundaryCells, totals in
     * table.localLower / localUpper / localGap. Records top per-(b,x_R) cells in
     * table.localCells to expose the highest pre-min gap slices.
     *
     * NOTE: This is a deterministic emat-level sandwich, NOT a certified table.
     */
    private void computeEmatOnlyTable(RegionAtomTable table, RCs rcs) {
        int[] region = table.spec.region;
        int[] boundary = table.spec.boundary;
        int[] xR = new int[region.length];
        int[] xB = new int[boundary.length];

        long t0 = System.currentTimeMillis();
        int[] bCellsContainer = new int[]{0};
        // Outer enumeration over boundary; for each boundary cell, inner enumeration over R.
        enumerateAssignments(boundary, 0, xB, rcs, () -> {
            String bKey = makeAssignmentKey(boundary, xB);
            // Track the largest pre-gap local cell within this boundary for [REGION_ATOM_LOCAL]
            String[] bestLocalKey = new String[1];
            double[] bestLocalPreGap = new double[]{Double.NEGATIVE_INFINITY};
            BigDecimal[] bestLocalGap = new BigDecimal[]{BigDecimal.ZERO};
            BigDecimal[] bestLocalLower = new BigDecimal[]{BigDecimal.ZERO};
            BigDecimal[] bestLocalUpper = new BigDecimal[]{BigDecimal.ZERO};

            enumerateAssignments(region, 0, xR, rcs, () -> {
                double ePlus = ownedLocalEnergy(branchRigidEmat, region, boundary, xR, xB);
                double eMinus = ownedLocalEnergy(branchMinimizingEmat, region, boundary, xR, xB);
                if (!Double.isFinite(ePlus) || !Double.isFinite(eMinus)) {
                    // +inf appears when sparse pruning marked this combination as pruned;
                    // skip — it is outside C_sparse.
                    return;
                }
                BigDecimal bL = bc.calc(ePlus);   // higher energy -> lower Z bound
                BigDecimal bU = bc.calc(eMinus);  // lower energy -> upper Z bound
                table.localLower = table.localLower.add(bL);
                table.localUpper = table.localUpper.add(bU);
                BigDecimal gap = bU.subtract(bL);
                if (gap.signum() > 0) {
                    table.localGap = table.localGap.add(gap);
                }
                // Per-cell (boundary) accumulator
                RegionAtomCell cell = table.boundaryCells.computeIfAbsent(bKey,
                        k -> new RegionAtomCell());
                cell.add(bL, bU, ePlus - eMinus);

                // Track largest-gap local within this boundary cell
                BigDecimal localGap = gap.signum() > 0 ? gap : BigDecimal.ZERO;
                if (bestLocalKey[0] == null || localGap.compareTo(bestLocalGap[0]) > 0) {
                    bestLocalKey[0] = makeAssignmentKey(region, xR);
                    bestLocalPreGap[0] = ePlus - eMinus;
                    bestLocalGap[0] = localGap;
                    bestLocalLower[0] = bL;
                    bestLocalUpper[0] = bU;
                }
            });

            // Save the best per-boundary local cell so we can later print [REGION_ATOM_LOCAL]
            if (bestLocalKey[0] != null) {
                RegionAtomCell lc = table.localCells.computeIfAbsent(bKey + "|" + bestLocalKey[0],
                        k -> new RegionAtomCell());
                lc.add(bestLocalLower[0], bestLocalUpper[0], bestLocalPreGap[0]);
            }
            bCellsContainer[0]++;
        });
        table.localBoundaryCellsComputed = bCellsContainer[0];
        table.localTableComputed = true;
        long t1 = System.currentTimeMillis();

        System.out.println("[REGION_ATOM_TABLE] name=" + table.spec.name
                + " status=computed"
                + " boundaryCells=" + table.localBoundaryCellsComputed
                + " sumL=" + formatBigExp(table.localLower)
                + " sumU=" + formatBigExp(table.localUpper)
                + " sumGap=" + formatBigExp(table.localGap)
                + " timeMs=" + (t1 - t0));

        // Top-K boundary cells by gap
        List<Map.Entry<String, RegionAtomCell>> sorted = new ArrayList<>(table.boundaryCells.entrySet());
        sorted.sort((a, b) -> b.getValue().gap.compareTo(a.getValue().gap));
        int topK = Math.min(regionAtomScoutTopCells, sorted.size());
        for (int i = 0; i < topK; i++) {
            Map.Entry<String, RegionAtomCell> e = sorted.get(i);
            RegionAtomCell cell = e.getValue();
            System.out.println("[REGION_ATOM_BOUNDARY] name=" + table.spec.name
                    + " phase=ematPre"
                    + " rank=" + (i + 1)
                    + " bKey=" + e.getKey()
                    + " localStates=" + cell.count
                    + " L=" + formatBigExp(cell.lower)
                    + " U=" + formatBigExp(cell.upper)
                    + " gap=" + formatBigExp(cell.gap)
                    + " meanPreGap=" + String.format("%.4f", cell.meanPreGap())
                    + " maxPreGap=" + (cell.preGapMax == Double.NEGATIVE_INFINITY
                            ? "-inf" : String.format("%.4f", cell.preGapMax)));
        }
        // Top per-boundary best local slice (already filtered to one per boundary)
        List<Map.Entry<String, RegionAtomCell>> sortedLocal = new ArrayList<>(table.localCells.entrySet());
        sortedLocal.sort((a, b) -> b.getValue().gap.compareTo(a.getValue().gap));
        int topL = Math.min(regionAtomScoutTopCells, sortedLocal.size());
        for (int i = 0; i < topL; i++) {
            Map.Entry<String, RegionAtomCell> e = sortedLocal.get(i);
            RegionAtomCell cell = e.getValue();
            System.out.println("[REGION_ATOM_LOCAL] name=" + table.spec.name
                    + " phase=ematPre"
                    + " rank=" + (i + 1)
                    + " key=" + e.getKey()
                    + " L=" + formatBigExp(cell.lower)
                    + " U=" + formatBigExp(cell.upper)
                    + " gap=" + formatBigExp(cell.gap)
                    + " preGap=" + String.format("%.4f", cell.preGapMax));
        }
    }

    private interface AssignmentSink { void visit(); }

    /** Recursive enumeration of RC assignments over {@code positions}; writes into {@code assign}. */
    private void enumerateAssignments(int[] positions, int idx, int[] assign, RCs rcs,
                                      AssignmentSink sink) {
        if (positions.length == 0) {
            sink.visit();
            return;
        }
        if (idx == positions.length) {
            sink.visit();
            return;
        }
        int pos = positions[idx];
        int[] rcList = rcs.get(pos);
        for (int r = 0; r < rcList.length; r++) {
            assign[idx] = rcList[r];
            enumerateAssignments(positions, idx + 1, assign, rcs, sink);
        }
    }

    // ========== Region-atom: Phase 1 queue scout (record leaf data) ==========

    /**
     * Record a leaf's current emat-level Z bounds and preGap into every active
     * region-atom table, keyed by the leaf's (boundary, local) RC assignment.
     * Idempotent per leaf conformation (uses {@link #regionAtomSeenLeafConfs}).
     *
     * Diagnostic only: does NOT modify the leaf, flatSumZ, queue, or DP tables.
     */
    private void recordScoutLeaf(DecompSearchNode leaf) {
        if (!regionAtomEnabled || regionAtomScoutComplete) return;
        if (leaf == null || leaf.isAggregate) return;
        if (!leaf.isLeaf() || leaf.minimized) return;
        int[] conf = leaf.partialConf;
        // Require fully assigned leaf (regular branch-DP leaves satisfy this).
        for (int i = 0; i < conf.length; i++) {
            if (conf[i] < 0) return;
        }
        String confKey = SimpleConfSpace.formatConfRCs(conf);
        if (!regionAtomSeenLeafConfs.add(confKey)) return;

        double eRigid = computeFullConfPairwiseEnergy(conf, branchRigidEmat);
        double eMin = computeFullConfPairwiseEnergy(conf, branchMinimizingEmat);
        if (!Double.isFinite(eRigid) || !Double.isFinite(eMin)) {
            // Conformations outside C_sparse should already be pruned at A*; guard anyway.
            return;
        }
        double preGap = eRigid - eMin;
        BigDecimal lowerZ = leaf.subtreeLowerBound;
        BigDecimal upperZ = leaf.subtreeUpperBound;

        for (RegionAtomTable table : regionAtomTables) {
            int[] region = table.spec.region;
            int[] boundary = table.spec.boundary;
            int[] xR = new int[region.length];
            int[] xB = new int[boundary.length];
            for (int i = 0; i < region.length; i++) xR[i] = conf[region[i]];
            for (int i = 0; i < boundary.length; i++) xB[i] = conf[boundary[i]];
            String bKey = makeAssignmentKey(boundary, xB);
            String lKey = makeAssignmentKey(region, xR);
            // Use suffix "@queue" on keys so they do not collide with Phase 0 emat-only cells.
            String bKeyQ = "queue|" + bKey;
            String lKeyQ = "queue|" + lKey;
            table.recordLeaf(bKeyQ, lKeyQ, lowerZ, upperZ, preGap);
        }

        regionAtomScoutLeavesSeen++;
        regionAtomScoutLower = regionAtomScoutLower.add(lowerZ);
        regionAtomScoutUpper = regionAtomScoutUpper.add(upperZ);
        BigDecimal gap = upperZ.subtract(lowerZ);
        if (gap.signum() > 0) regionAtomScoutGap = regionAtomScoutGap.add(gap);

        if (regionAtomScoutLeavesSeen >= regionAtomScoutMaxLeaves) {
            regionAtomScoutComplete = true;
        }
    }

    /**
     * Print Phase 1 queue-scout summary and Phase 2 what-if epsilon-drop projections.
     * Called once from {@link #compute(int)} after the search loop terminates
     * (normal completion or scoutOnly early-exit).
     */
    private void printRegionAtomSummary() {
        if (!regionAtomEnabled || regionAtomScoutSummaryPrinted) return;
        regionAtomScoutSummaryPrinted = true;

        BigDecimal zUpper = flatSumZUpper.signum() > 0 ? flatSumZUpper : BigDecimal.ONE;
        double epsBefore = computeEpsilon(flatSumZLower, flatSumZUpper);

        System.out.println("[REGION_ATOM_SCOUT] state=" + stateName
                + " leavesSeen=" + regionAtomScoutLeavesSeen
                + " maxLeaves=" + regionAtomScoutMaxLeaves
                + " sampledZLower=" + formatBigExp(regionAtomScoutLower)
                + " sampledZUpper=" + formatBigExp(regionAtomScoutUpper)
                + " sampledZGap=" + formatBigExp(regionAtomScoutGap)
                + " globalZLower=" + formatBigExp(flatSumZLower)
                + " globalZUpper=" + formatBigExp(flatSumZUpper)
                + " epsilon=" + String.format("%.6f", epsBefore));

        // Phase 4 DP-replacement stats (no-op when regionAtomCertifyUseDP=false).
        long certifiedTables = regionAtomTables.stream().filter(t -> t.certified).count();
        if (certifiedTables > 0) {
            System.out.println("[REGION_ATOM_DP] state=" + stateName
                    + " certifyUseDP=" + regionAtomCertifyUseDP
                    + " certifiedTables=" + certifiedTables
                    + " applied=" + regionAtomDPApplied
                    + " skippedPartial=" + regionAtomDPSkippedPartial
                    + " missingBoundary=" + regionAtomDPMissingBoundary
                    + " missingLocal=" + regionAtomDPMissingLocal);
        }

        for (RegionAtomTable table : regionAtomTables) {
            // Queue-scout boundary cells (suffix "queue|" key)
            List<Map.Entry<String, RegionAtomCell>> queueCells = new ArrayList<>();
            BigDecimal tableQueueGap = BigDecimal.ZERO;
            int tableQueueLeaves = 0;
            for (Map.Entry<String, RegionAtomCell> e : table.boundaryCells.entrySet()) {
                if (!e.getKey().startsWith("queue|")) continue;
                queueCells.add(e);
                tableQueueGap = tableQueueGap.add(e.getValue().gap);
                tableQueueLeaves += e.getValue().count;
            }
            queueCells.sort((a, b) -> b.getValue().gap.compareTo(a.getValue().gap));

            System.out.println("[REGION_ATOM_SCOUT] name=" + table.spec.name
                    + " phase=queue"
                    + " boundaryCells=" + queueCells.size()
                    + " leaves=" + tableQueueLeaves
                    + " sumQueueGap=" + formatBigExp(tableQueueGap)
                    + " ematPreGap=" + formatBigExp(table.localGap)
                    + " ematSkipped=" + table.localTableSkipped);

            int topK = Math.min(regionAtomScoutTopCells, queueCells.size());
            for (int i = 0; i < topK; i++) {
                Map.Entry<String, RegionAtomCell> e = queueCells.get(i);
                RegionAtomCell cell = e.getValue();
                System.out.println("[REGION_ATOM_BOUNDARY] name=" + table.spec.name
                        + " phase=queue"
                        + " rank=" + (i + 1)
                        + " bKey=" + e.getKey().substring("queue|".length())
                        + " leaves=" + cell.count
                        + " L=" + formatBigExp(cell.lower)
                        + " U=" + formatBigExp(cell.upper)
                        + " gap=" + formatBigExp(cell.gap)
                        + " meanPreGap=" + String.format("%.4f", cell.meanPreGap())
                        + " maxPreGap=" + (cell.preGapMax == Double.NEGATIVE_INFINITY
                                ? "-inf" : String.format("%.4f", cell.preGapMax)));
            }

            // Phase 2 what-if: simulate per-cell upper shrink by factor exp(-delta/RT).
            // We model each boundary cell as "if we certified this cell and shrank U by
            // factor f = exp(-delta/RT)" and report the resulting epsilon and jobs/drop ratio.
            // jobs = sum of (count of x_R observations) on this cell — used for jobs/drop ratio
            // approximation when comparing atoms.
            if (regionAtomWhatIfDeltas.length == 0 || queueCells.isEmpty()) continue;
            for (double delta : regionAtomWhatIfDeltas) {
                double rt = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
                double factor = Math.exp(-delta / rt);  // shrink U by this factor
                BigDecimal bigFactor = BigDecimal.valueOf(factor);

                // Apply shrink to the top topK boundary cells; cap the total shrink to
                // the queue-sampled upper mass attributable to those cells.
                BigDecimal upperDrop = BigDecimal.ZERO;
                BigDecimal lowerLift = BigDecimal.ZERO;
                int jobs = 0;
                int cellsUsed = Math.min(topK, queueCells.size());
                for (int i = 0; i < cellsUsed; i++) {
                    RegionAtomCell cell = queueCells.get(i).getValue();
                    BigDecimal newU = cell.upper.multiply(bigFactor);
                    BigDecimal newL = cell.lower; // lower stays in this model
                    if (newU.compareTo(cell.lower) < 0) newU = cell.lower;
                    upperDrop = upperDrop.add(cell.upper.subtract(newU));
                    lowerLift = lowerLift.add(newL.subtract(cell.lower));
                    jobs += cell.count;
                }
                BigDecimal newZUpper = flatSumZUpper.subtract(upperDrop);
                BigDecimal newZLower = flatSumZLower.add(lowerLift);
                if (newZUpper.signum() <= 0) newZUpper = BigDecimal.ONE;
                double epsAfter = computeEpsilon(newZLower, newZUpper);
                double epsDrop = epsBefore - epsAfter;
                double dropPerJob = jobs > 0 ? epsDrop / jobs : 0.0;
                BigDecimal sharedGap = upperDrop;
                double shareOfQueueGap = tableQueueGap.signum() > 0
                        ? MathTools.bigDivide(sharedGap, tableQueueGap, PartitionFunction.decimalPrecision).doubleValue()
                        : 0.0;
                System.out.println("[REGION_ATOM_WHATIF] name=" + table.spec.name
                        + " delta=" + String.format("%.4f", delta)
                        + " shrinkFactor=" + String.format("%.4e", factor)
                        + " topCells=" + cellsUsed
                        + " jobs=" + jobs
                        + " upperDrop=" + formatBigExp(upperDrop)
                        + " lowerLift=" + formatBigExp(lowerLift)
                        + " epsBefore=" + String.format("%.6f", epsBefore)
                        + " epsAfter=" + String.format("%.6f", epsAfter)
                        + " epsDrop=" + String.format("%.6f", epsDrop)
                        + " epsDropPerJob=" + String.format("%.3e", dropPerJob)
                        + " shareOfQueueGap=" + String.format("%.4f", shareOfQueueGap));
            }
        }
    }

    // ========== Region-atom: Phase 3 certified CPU oracle ==========

    /** Parse "7;6,7,8;0,2,6,7,8" into the matching {@code formatPositions} names. */
    private static Set<String> parseRegionAtomCertifyNames(String spec) {
        Set<String> out = new HashSet<>();
        if (spec == null || spec.trim().isEmpty()) return out;
        for (String token : spec.split(";")) {
            String tok = token.trim();
            if (tok.isEmpty()) continue;
            // Normalize by sorting the comma-separated positions so user can pass
            // any order; the canonical name uses ascending order.
            String[] parts = tok.split(",");
            TreeSet<Integer> sorted = new TreeSet<>();
            boolean ok = true;
            for (String p : parts) {
                try {
                    sorted.add(Integer.parseInt(p.trim()));
                } catch (NumberFormatException e) {
                    ok = false;
                    System.err.println("[REGION_ATOM] WARN certify name '" + tok
                            + "' contains non-integer '" + p + "'; ignoring");
                    break;
                }
            }
            if (!ok || sorted.isEmpty()) continue;
            int[] arr = sorted.stream().mapToInt(Integer::intValue).toArray();
            out.add(formatPositions(arr));
        }
        return out;
    }

    /**
     * Phase 3 certified table: for each (b, x_R) call the CCD energy oracle on the
     * tuple {@code (R ∪ B)} with R-owned residue interactions and accumulate
     *   L_R(b) = sum_{x_R} exp(-E^+_rigid(x_R; b) / RT)   [pre-min lower bound]
     *   U_R(b) = sum_{x_R} exp(-E_R^cert(x_R; b) / RT)    [CCD-min upper bound]
     *
     * Submits via {@code minimizingEcalc.calcEnergyAsync} so jobs run on the
     * task executor (CPU pool by default; CUDA-backed pool if the ecalc is
     * compiled GPU). Waits for completion before returning.
     *
     * NOTE: tighter U_R via joint R+B CCD with R-owned terms is sound (it
     * gives a smaller min energy than any constrained-DOF minimization,
     * hence a larger exp value, hence an UPPER bound on Z_R(b)).
     */
    private void certifyRegionAtomTable(RegionAtomTable table, RCs rcs) {
        int[] region = table.spec.region;
        int[] boundary = table.spec.boundary;
        int[] xR = new int[region.length];
        int[] xB = new int[boundary.length];

        long t0 = System.currentTimeMillis();
        final java.util.concurrent.atomic.AtomicLong jobsSubmitted = new java.util.concurrent.atomic.AtomicLong(0);
        final java.util.concurrent.atomic.AtomicLong jobsPruned = new java.util.concurrent.atomic.AtomicLong(0);
        final int totalLocalPerBoundary = table.spec.localStatesPerBoundary.intValueExact();

        // Optional top-K subsampling: certify only the K boundary cells with the
        // largest Phase 0 emat gap. The rest fall back to the pre-min sandwich
        // (still valid bounds, just not tightened). This makes Phase 3 feasible
        // for atoms with many boundary cells (e.g. {7} has 41040 cells).
        final Set<String> certifyAllowedBoundary;
        if (regionAtomCertifyTopK > 0
                && table.boundaryCells.size() > regionAtomCertifyTopK) {
            List<Map.Entry<String, RegionAtomCell>> sorted = new ArrayList<>(table.boundaryCells.entrySet());
            sorted.sort((a, b) -> b.getValue().gap.compareTo(a.getValue().gap));
            certifyAllowedBoundary = new HashSet<>();
            int k = Math.min(regionAtomCertifyTopK, sorted.size());
            for (int i = 0; i < k; i++) certifyAllowedBoundary.add(sorted.get(i).getKey());
            System.out.println("[REGION_ATOM_CERT] name=" + table.spec.name
                    + " topK=" + regionAtomCertifyTopK
                    + " totalBoundaryCells=" + table.boundaryCells.size()
                    + " certifying=" + certifyAllowedBoundary.size()
                    + " skipping=" + (table.boundaryCells.size() - certifyAllowedBoundary.size())
                    + " (others fall back to Phase 0 sandwich)");
        } else {
            certifyAllowedBoundary = null;  // certify all
        }

        enumerateAssignments(boundary, 0, xB, rcs, () -> {
            final String bKey = makeAssignmentKey(boundary, xB);
            if (certifyAllowedBoundary != null && !certifyAllowedBoundary.contains(bKey)) {
                return;  // skipped boundary cell — uses Phase 0 pre-min bound
            }
            final int[] xBSnapshot = xB.clone();
            final CertifiedBoundaryCell cell = table.certifiedCells.computeIfAbsent(bKey,
                    k -> new CertifiedBoundaryCell(bKey, xBSnapshot, totalLocalPerBoundary));

            enumerateAssignments(region, 0, xR, rcs, () -> {
                final int[] xRSnapshot = xR.clone();

                // Drop combinations pruned by sparse pfunc (rigid or minimizing emat returns +inf).
                double rigidLocal = ownedLocalEnergy(branchRigidEmat, region, boundary,
                        xRSnapshot, xBSnapshot);
                double emins = ownedLocalEnergy(branchMinimizingEmat, region, boundary,
                        xRSnapshot, xBSnapshot);
                if (!Double.isFinite(rigidLocal) || !Double.isFinite(emins)) {
                    cell.markPruned();
                    jobsPruned.incrementAndGet();
                    return;
                }

                final BigDecimal bL = bc.calc(rigidLocal);  // L stays at rigid (pre-min lower bound)

                RCTuple tuple = makeRegionAtomTuple(region, boundary, xRSnapshot, xBSnapshot);
                ResidueInteractions inters = makeRegionAtomInteractions(region, boundary,
                        xRSnapshot, xBSnapshot);

                long submitted = jobsSubmitted.incrementAndGet();
                if (submitted % 2000 == 0) {
                    System.out.println("[REGION_ATOM_CERT_PROGRESS] name=" + table.spec.name
                            + " submitted=" + submitted
                            + " elapsedMs=" + (System.currentTimeMillis() - t0));
                }
                final String localKey = makeAssignmentKey(region, xRSnapshot);
                minimizingEcalc.calcEnergyAsync(tuple, inters, epmol -> {
                    double eCert = epmol.energy;
                    BigDecimal bU = bc.calc(eCert);
                    cell.addLocal(localKey, bL, bU, eCert);
                });
            });
        });

        minimizingEcalc.tasks.waitForFinish();
        long t1 = System.currentTimeMillis();

        // Aggregate per-table sums
        BigDecimal sumL = BigDecimal.ZERO;
        BigDecimal sumU = BigDecimal.ZERO;
        BigDecimal sumGap = BigDecimal.ZERO;
        for (CertifiedBoundaryCell c : table.certifiedCells.values()) {
            sumL = sumL.add(c.lowerZ);
            sumU = sumU.add(c.upperZ);
            BigDecimal gap = c.upperZ.subtract(c.lowerZ);
            if (gap.signum() > 0) sumGap = sumGap.add(gap);
        }
        table.certifiedSumLower = sumL;
        table.certifiedSumUpper = sumU;
        table.certifiedSumGap = sumGap;
        table.certified = true;

        String gapShrinkStr = "n/a";
        if (table.localGap.signum() > 0) {
            double shrink = MathTools.bigDivide(sumGap, table.localGap,
                    PartitionFunction.decimalPrecision).doubleValue();
            gapShrinkStr = String.format("%.6f", shrink);
        }

        System.out.println("[REGION_ATOM_CERT] name=" + table.spec.name
                + " status=certified"
                + " boundaryCells=" + table.certifiedCells.size()
                + " jobs=" + jobsSubmitted.get()
                + " pruned=" + jobsPruned.get()
                + " sumL=" + formatBigExp(sumL)
                + " sumU=" + formatBigExp(sumU)
                + " sumGap=" + formatBigExp(sumGap)
                + " ematSumGap=" + formatBigExp(table.localGap)
                + " gapShrink=" + gapShrinkStr
                + " timeMs=" + (t1 - t0));

        // Top-K certified boundary cells by gap
        List<CertifiedBoundaryCell> cellList = new ArrayList<>(table.certifiedCells.values());
        cellList.sort((a, b) -> {
            BigDecimal ag = a.upperZ.subtract(a.lowerZ);
            BigDecimal bg = b.upperZ.subtract(b.lowerZ);
            return bg.compareTo(ag);
        });
        int topK = Math.min(regionAtomScoutTopCells, cellList.size());
        for (int i = 0; i < topK; i++) {
            CertifiedBoundaryCell c = cellList.get(i);
            BigDecimal gap = c.upperZ.subtract(c.lowerZ);
            if (gap.signum() < 0) gap = BigDecimal.ZERO;
            System.out.println("[REGION_ATOM_CERT_BOUNDARY] name=" + table.spec.name
                    + " rank=" + (i + 1)
                    + " bKey=" + c.bKey
                    + " localStatesDone=" + c.localStatesDone
                    + " pruned=" + c.localStatesPruned
                    + " L=" + formatBigExp(c.lowerZ)
                    + " U=" + formatBigExp(c.upperZ)
                    + " gap=" + formatBigExp(gap)
                    + " minE=" + String.format("%.4f", c.minCertEnergy)
                    + " maxE=" + String.format("%.4f", c.maxCertEnergy));
        }
    }

    /** Combine R ∪ B positions/assignments into a sorted RCTuple suitable for calcEnergy. */
    private static RCTuple makeRegionAtomTuple(int[] region, int[] boundary, int[] xR, int[] xB) {
        int n = region.length + boundary.length;
        int[][] entries = new int[n][2];
        int idx = 0;
        for (int i = 0; i < region.length; i++) {
            entries[idx][0] = region[i];
            entries[idx][1] = xR[i];
            idx++;
        }
        for (int i = 0; i < boundary.length; i++) {
            entries[idx][0] = boundary[i];
            entries[idx][1] = xB[i];
            idx++;
        }
        Arrays.sort(entries, (a, b) -> Integer.compare(a[0], b[0]));
        ArrayList<Integer> posList = new ArrayList<>(n);
        ArrayList<Integer> rcList = new ArrayList<>(n);
        for (int[] e : entries) {
            posList.add(e[0]);
            rcList.add(e[1]);
        }
        return new RCTuple(posList, rcList);
    }

    /**
     * R-owned ResidueInteractions:
     *   - one-body single for each R position (includes shell)
     *   - pair for every R-R combination
     *   - pair for every R-B crossing
     * Excludes B singles, B-B pairs, and all outside terms.
     */
    /**
     * R-owned ResidueInteractions matching the emat's {@link EnergyPartition} convention.
     *
     * The emat values include offsets (reference energies, residue entropy) and shell
     * interactions according to the chosen partition (Traditional puts shell+intra on
     * singles; AllOnPairs distributes intra+shell across pairs with weight 1/(n-1)).
     * To make calcEnergy(tuple, inters) return a value comparable to the emat sum, we
     * must use the SAME {@code epart.makeSingle} / {@code epart.makePair} that built
     * the emat — bare {@code addSingle}/{@code addPair} would silently drop offsets
     * and shell, producing energies that are systematically too negative (this was
     * the early-bug: certified U_R came out 10^3x larger than the pre-min sandwich).
     */
    private ResidueInteractions makeRegionAtomInteractions(int[] region, int[] boundary,
                                                           int[] xR, int[] xB) {
        ResidueInteractions inters = new ResidueInteractions();
        EnergyPartition epart = minimizingEcalc.epart;
        SimpleReferenceEnergies eref = minimizingEcalc.eref;
        boolean addEntropy = minimizingEcalc.addResEntropy;

        // R one-body
        for (int i = 0; i < region.length; i++) {
            inters.addAll(epart.makeSingle(confSpace, eref, addEntropy, region[i], xR[i]));
        }
        // R-R pairs
        for (int i = 0; i < region.length; i++) {
            for (int j = i + 1; j < region.length; j++) {
                inters.addAll(epart.makePair(confSpace, eref, addEntropy,
                        region[i], xR[i], region[j], xR[j]));
            }
        }
        // R-B crossings
        for (int i = 0; i < region.length; i++) {
            for (int j = 0; j < boundary.length; j++) {
                inters.addAll(epart.makePair(confSpace, eref, addEntropy,
                        region[i], xR[i], boundary[j], xB[j]));
            }
        }
        return inters;
    }

    // ========== Region-atom: Phase 4 DP factor replacement ==========

    /** Per-pfunc tally of Phase 4 hits across all node tighten attempts. */
    private long regionAtomDPApplied = 0;
    /** Accumulated flatSumZUpper delta from region-atom tightening (flushed once per round). */
    private BigDecimal regionAtomDPUpperDelta = BigDecimal.ZERO;
    private long regionAtomDPSkippedPartial = 0;
    private long regionAtomDPMissingBoundary = 0;
    private long regionAtomDPMissingLocal = 0;

    /**
     * DP-level integration: apply certified region-atom corrections directly into
     * DP table entries. For each lambda-edge whose lambda set covers R, and for
     * each (mIdx, lambdaIdx) entry, replace the R-owned energy from minimizingEmat
     * with the certified energy eCert. This tightens Z_upper at the DP source,
     * propagating to ALL nodes without per-node patching.
     *
     * Soundness: eCert >= eMin_R_owned (certified CCD min >= pairwise emat min),
     * so fullEnergyMin increases -> Boltzmann decreases -> Z_upper tighter.
     */
    private void applyRegionAtomToDP(double RT) {
        // Collect all lambda edges
        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);

        for (RegionAtomTable table : regionAtomTables) {
            if (!table.certified) continue;
            int[] region = table.spec.region;
            int[] boundary = table.spec.boundary;

            for (RootedTreeEdge edge : lambdaEdges) {
                int[] lambdaPos = edge.getLambdaPositionsSorted();
                int[] mPos = edge.getMPositionsSorted();

                System.out.println("[REGION_ATOM_DP_DEBUG] checking edge lambda="
                        + Arrays.toString(lambdaPos) + " m=" + Arrays.toString(mPos)
                        + " R=" + Arrays.toString(region) + " B=" + Arrays.toString(boundary));

                // Check: R must be a subset of lambda (so R is assigned at this edge)
                boolean rInLambda = true;
                for (int rp : region) {
                    boolean found = false;
                    for (int lp : lambdaPos) { if (lp == rp) { found = true; break; } }
                    if (!found) { rInLambda = false; break; }
                }
                if (!rInLambda) {
                    System.out.println("[REGION_ATOM_DP_DEBUG]   SKIP: R not in lambda");
                    continue;
                }

                // Check: B positions must all be in lambda ∪ M
                boolean bCovered = true;
                for (int bp : boundary) {
                    boolean found = false;
                    for (int lp : lambdaPos) { if (lp == bp) { found = true; break; } }
                    if (!found) {
                        for (int mp : mPos) { if (mp == bp) { found = true; break; } }
                    }
                    if (!found) {
                        System.out.println("[REGION_ATOM_DP_DEBUG]   SKIP: B pos " + bp + " not in lambda∪M");
                        bCovered = false; break;
                    }
                }
                if (!bCovered) continue;

                double[][] fullEnergyMin = edge.getFullEnergyMin();
                if (fullEnergyMin == null) continue;

                int corrected = 0;
                for (int mIdx = 0; mIdx < edge.getMArraySize(); mIdx++) {
                    int[] mRCs = edge.decodeMStatePublic(mIdx);
                    for (int lIdx = 0; lIdx < edge.getTotalLambdaStates(); lIdx++) {
                        if (Double.isNaN(fullEnergyMin[mIdx][lIdx])) continue;

                        int[] lambdaRCs = edge.decodeLambdaStatePublic(lIdx);

                        // Build xR and xB from lambda+M RCs
                        int[] xR = new int[region.length];
                        int[] xB = new int[boundary.length];
                        for (int i = 0; i < region.length; i++) {
                            xR[i] = resolveRC(region[i], lambdaPos, lambdaRCs, mPos, mRCs, edge.getRCs());
                        }
                        boolean valid = true;
                        for (int i = 0; i < boundary.length; i++) {
                            int rc = resolveRC(boundary[i], lambdaPos, lambdaRCs, mPos, mRCs, edge.getRCs());
                            if (rc < 0) { valid = false; break; }
                            xB[i] = rc;
                        }
                        if (!valid) continue;

                        // Look up certified energy
                        String bKey = makeAssignmentKey(boundary, xB);
                        CertifiedBoundaryCell cell = table.certifiedCells.get(bKey);
                        if (cell == null) continue;
                        String lKey = makeAssignmentKey(region, xR);
                        Double eCert = cell.certifiedEnergyPerLocal.get(lKey);
                        if (eCert == null) continue;

                        // Compute R-owned energy from minimizing emat
                        double eMinOwned = ownedLocalEnergy(branchMinimizingEmat, region, boundary, xR, xB);
                        if (!Double.isFinite(eMinOwned)) continue;

                        double correction = eCert - eMinOwned;
                        if (correction > 0) {
                            fullEnergyMin[mIdx][lIdx] += correction;
                            corrected++;
                        }
                    }
                }

                if (corrected > 0) {
                    // Re-sort lambda indices since energies changed
                    edge.initIncrementalEnumeration(
                            branchRigidEmat, branchMinimizingEmat, interactionGraph, RT, true);
                    // Re-apply corrections after re-init (which resets energies)
                    // Actually, let's apply corrections differently - modify after init
                    // For now, just apply corrections to the re-initialized table
                    fullEnergyMin = edge.getFullEnergyMin();
                    if (fullEnergyMin != null) {
                        int reapplied = 0;
                        for (int mIdx = 0; mIdx < edge.getMArraySize(); mIdx++) {
                            int[] mRCs = edge.decodeMStatePublic(mIdx);
                            for (int lIdx = 0; lIdx < edge.getTotalLambdaStates(); lIdx++) {
                                if (Double.isNaN(fullEnergyMin[mIdx][lIdx])) continue;
                                int[] lambdaRCs = edge.decodeLambdaStatePublic(lIdx);
                                int[] xR = new int[region.length];
                                int[] xB = new int[boundary.length];
                                for (int i = 0; i < region.length; i++)
                                    xR[i] = resolveRC(region[i], lambdaPos, lambdaRCs, mPos, mRCs, edge.getRCs());
                                boolean v2 = true;
                                for (int i = 0; i < boundary.length; i++) {
                                    int rc = resolveRC(boundary[i], lambdaPos, lambdaRCs, mPos, mRCs, edge.getRCs());
                                    if (rc < 0) { v2 = false; break; }
                                    xB[i] = rc;
                                }
                                if (!v2) continue;
                                String bKey = makeAssignmentKey(boundary, xB);
                                CertifiedBoundaryCell cell = table.certifiedCells.get(bKey);
                                if (cell == null) continue;
                                String lKey = makeAssignmentKey(region, xR);
                                Double eCert = cell.certifiedEnergyPerLocal.get(lKey);
                                if (eCert == null) continue;
                                double eMinOwned = ownedLocalEnergy(branchMinimizingEmat, region, boundary, xR, xB);
                                if (!Double.isFinite(eMinOwned)) continue;
                                double corr = eCert - eMinOwned;
                                if (corr > 0) {
                                    fullEnergyMin[mIdx][lIdx] += corr;
                                    reapplied++;
                                }
                            }
                        }
                        corrected = reapplied;
                    }

                    System.out.println("[REGION_ATOM_DP_INTEGRATE] edge lambda="
                            + Arrays.toString(lambdaPos) + " corrected=" + corrected
                            + "/" + (edge.getMArraySize() * edge.getTotalLambdaStates())
                            + " entries");
                }
            }
        }
    }

    /** Resolve a position's RC from lambda or M arrays. Returns the actual RC index. */
    private int resolveRC(int pos, int[] lambdaPos, int[] lambdaRCs,
                          int[] mPos, int[] mRCs, RCs edgeRCs) {
        for (int i = 0; i < lambdaPos.length; i++) {
            if (lambdaPos[i] == pos) return edgeRCs.get(pos, lambdaRCs[i]);
        }
        for (int i = 0; i < mPos.length; i++) {
            if (mPos[i] == pos) return edgeRCs.get(pos, mRCs[i]);
        }
        return -1;
    }

    /**
     * Region-atom Phase 4 hook: try to tighten a search node's bounds using any
     * applicable certified region-atom table. Sound replacement of the R-owned
     * slice of the leaf bound:
     *
     *   subtreeUpperBound_new = subtreeUpperBound_old * bc.calc(E_R_cert(x_R; b))
     *                                                / bc.calc(E_R_emat_min(x_R; b))
     *
     * Because E_R_cert ≥ E_R_emat_min (joint CCD min ≥ pair-wise emat min lower
     * bound), the ratio is ≤ 1, so the new upper bound is tighter while remaining
     * a valid upper bound on Z (E_R_cert ≤ E_R_true with B's DOFs free).
     *
     * Requires the node to have all positions in R ∪ B assigned. Partial nodes
     * are skipped (counted in regionAtomDPSkippedPartial).
     *
     * This is the per-leaf tightening flavor of DP factor replacement. A proper
     * DP factor replacement that aligns the atom's (R, B) with a branch-decomp
     * separator and replaces the edge's logZ table is a follow-up; this hook
     * already exercises the certified table end-to-end.
     */
    private void tightenNodeWithCertifiedTables(DecompSearchNode node) {
        // Per-node tightening disabled: DP-level integration (applyRegionAtomToDP)
        // bakes corrections into DP table entries before search starts.
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
                + "-Dbranchmarkstar.dp.cache=false for dense-only runs, or enable the shard-backed DP table once available.");
    }

    private void logDPMemoryPredictions(String context, double logTESS) {
        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);

        long totalFinalBytes = 0L;
        long totalCacheBytes = 0L;
        for (int edgeId = 0; edgeId < lambdaEdges.size(); edgeId++) {
            RootedTreeEdge edge = lambdaEdges.get(edgeId);
            long finalBytes = estimateDPTableBytes(edge.getMStateCount());
            long cacheBytes = (dpCacheEnabled && dpCacheMaxEntries > 0) ? finalBytes : 0L;
            totalFinalBytes = saturatingAdd(totalFinalBytes, finalBytes);
            totalCacheBytes = saturatingAdd(totalCacheBytes, cacheBytes);

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
                    + ", branchwidth=" + branchwidth
                    + ", logTESS=" + String.format(Locale.ROOT, "%.2f", logTESS));
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " DP memory prediction summary"
                + " context=" + context
                + ", state=" + stateName
                + ", lambdaEdges=" + lambdaEdges.size()
                + ", totalFinalTableBytes=" + formatBytes(totalFinalBytes)
                + ", totalCacheCopyBytes=" + formatBytes(totalCacheBytes)
                + ", dpCache=" + (dpCacheEnabled && dpCacheMaxEntries > 0)
                + ", dpCacheMaxTable=" + formatBytes(dpCacheMaxTableBytes)
                + ", dpCacheMaxTotal=" + formatBytes(dpCacheMaxTotalBytes)
                + ", branchwidth=" + branchwidth
                + ", logTESS=" + String.format(Locale.ROOT, "%.2f", logTESS));
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
                Iterator<String> it = DP_TABLE_CACHE.keySet().iterator();
                if (!it.hasNext()) break;
                String evictKey = it.next();
                CachedDPTable evicted = DP_TABLE_CACHE.get(evictKey);
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
                    regionAtomCertifyUseDP);
            DPCacheStats dpCacheStats = computeFullDPTables("initial", !regionAtomCertifyUseDP);

            // Region-atom DP-level integration: apply certified corrections
            // directly into the DP table entries BEFORE the full DP propagation
            // is used by the search. This tightens Z_upper at the source.
            if (regionAtomCertifyUseDP && !regionAtomTables.isEmpty()) {
                applyRegionAtomToDP(rtForDP);
                // Recompute full DP after modifying leaf-edge energies
                RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
            }

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

            // Phase 4 region-atom hook: lazily tighten bounds on every pulled node.
            // No-op when regionAtomCertifyUseDP=false or when no atom's R+B is fully
            // covered by this node's partialConf. Updates flatSumZ as needed.
            tightenNodeWithCertifiedTables(node);

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

        // Flush accumulated region-atom delta in one synchronized batch.
        if (regionAtomDPUpperDelta.signum() != 0) {
            synchronized (this) {
                flatSumZUpper = flatSumZUpper.add(regionAtomDPUpperDelta);
            }
            regionAtomDPUpperDelta = BigDecimal.ZERO;
        }

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

        // Region-atom queue scout: record every pulled leaf into the diagnostic
        // tables BEFORE the decision branch fires. This is diagnostic only —
        // it does NOT mutate the leaf, flatSumZ, or the queue.
        if (regionAtomEnabled && !leafNodes.isEmpty()) {
            for (DecompSearchNode leaf : leafNodes) {
                recordScoutLeaf(leaf);
            }
        }

        // Scout-only mode: skip leaf minimization entirely and force the internal-
        // expansion branch so the queue keeps producing new leaves. Pulled leaves
        // are dropped (their flatSumZ contribution becomes orphaned, which is
        // acceptable in diagnostic mode since scoutOnly does NOT certify a pfunc).
        if (regionAtomScoutOnly) {
            if (parallelInternal && parallelism.numThreads > 1 && internalNodes.size() > 1) {
                processInternalNodesParallel(internalNodes, drillDownThreshold);
            } else {
                List<DecompSearchNode> newNodes = new ArrayList<>();
                for (int i = 0; i < internalNodes.size(); i++) {
                    DecompSearchNode internal = internalNodes.get(i);
                    if (regionAtomScoutComplete) {
                        // Stop expanding; defer remaining internals back to queue
                        for (int j = i; j < internalNodes.size(); j++) {
                            decompQueue.add(internalNodes.get(j));
                        }
                        break;
                    }
                    if (!internal.isAggregate && shouldDrillDown(internal, drillDownThreshold)) {
                        drillDownFullEnum(internal, newNodes, drillDownThreshold);
                    } else {
                        expandDecompNode(internal);
                        updateBound();
                    }
                }
                for (DecompSearchNode n : newNodes) {
                    if (n.errorBound.signum() > 0) decompQueue.add(n);
                }
            }
            // Drop leafNodes: do not minimize, do not re-add.
            updateBound();
            return;
        }

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

        // Phase 4: try certified region-atom tables on this node (no-op when
        // regionAtomCertifyUseDP=false or when no atom's R+B covers this node).
        tightenNodeWithCertifiedTables(node);
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
            // Use the original (pre-tightening) upper bound for the aggregate
            // computation. If the node was tightened by region-atom, its
            // subtreeUpperBound is smaller than the DP-based children sum,
            // which would make the aggregate negative. Using the original
            // ensures children + aggregate = original (consistent with DP).
            BigDecimal parentUpperForAggregate = node.regionAtomOriginalUpper != null
                    ? node.regionAtomOriginalUpper : node.subtreeUpperBound;
            DecompSearchNode aggregate = DecompSearchNode.makeAggregate(
                    node, pendingIdx, lambdaAStar,
                    parentUpperForAggregate, childrenZUpper, bc);
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

            // Region-atom scoutOnly: terminate after enough leaves observed or
            // when the queue drains (diagnostic mode does NOT certify pfunc).
            if (regionAtomScoutOnly
                    && (regionAtomScoutComplete || decompQueue.isEmpty())) {
                System.out.println("[REGION_ATOM] scoutOnly early-exit: leavesSeen="
                        + regionAtomScoutLeavesSeen
                        + " maxLeaves=" + regionAtomScoutMaxLeaves
                        + " queueSize=" + decompQueue.size()
                        + " epsilon=" + String.format("%.6f", epsilonBound));
                break;
            }
        }

        System.out.println(BranchDpConfig.getBackendLogPrefix() + " Finished. epsilon=" + String.format("%.6f", epsilonBound)
                + " after " + totalEnumerationSteps + " enum steps, "
                + totalMinimizations + " minimizations, "
                + totalDrillDowns + " drill-downs, "
                + totalDrillDownEarlyStops + " early stops.");
        printLeafMinimizationProfile();
        printEdgeSelectionStats();
        printCorrectionAuditStats();
        printRegionAtomSummary();
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

        if (regionAtomScoutOnly) {
            // scoutOnly is diagnostic only — do NOT mark as Estimated; bounds are not certified.
            System.out.println(BranchDpConfig.getBackendLogPrefix() + " scoutOnly mode complete; pfunc NOT certified.");
        } else if (reportedEpsilon <= targetEpsilon) {
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
