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

package edu.duke.cs.osprey.markstar.framework;

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
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.branch.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.*;
import java.util.stream.IntStream;

/**
 * BranchMARK*: Partition function computation using branch decomposition.
 *
 * Top-down edge-level traversal of the decomposition tree.
 * Each "expand" assigns an entire edge's lambda-set.
 * MARK*-style priority queue decides between expand (internal) and minimize (leaf).
 * Z bounds maintained incrementally from h-scores (no DP needed).
 * O(1) bound updates per operation.
 */
public class BranchMARKStarBound extends MARKStarBound {

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
    private static final String ROOT_SPLIT_MAX_FSET_PROPERTY = "branchmarkstar.rootSplit.maxFset";
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

    private enum EdgeSelectionStrategy {
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
                    System.err.println("BranchMARK*: Unknown edge selection strategy '" + value
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
                    System.err.println("BranchMARK*: Unknown cutoff strategy '" + value
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
     *    "cut-but-large-interaction" full confs leak in. Triple corrections,
     *    GNN, and S8 subtree GNN are disabled because they are not certified
     *    against E_sparse / C_sparse.
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
                    System.err.println("BranchMARK*: Unknown energy mode '" + value
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
    /**
     * Energy matrices used by ALL pfunc-relevant operations under SPARSE mode:
     * scorers, A*, computeFullConfPairwiseEnergy, breakdownScoreByPosition,
     * tuple internal-energy queries, and leaf pairwise bounds. In SPARSE mode
     * these are copies of {@code rigidEmat}/{@code minimizingEmat} with cut
     * edges processed by {@link #makeSparseEnergyMatrix} (zeroed if accepted,
     * +inf if pruned). In FULL mode they alias the originals so behavior is
     * identical to the pre-sparse build.
     */
    private EnergyMatrix branchRigidEmat;
    private EnergyMatrix branchMinimizingEmat;
    private int sparsePrunedRcPairs = 0;
    private int sparseZeroedRcPairs = 0;
    private int sparseTotalCutPositionPairs = 0;

    /** True once {@link RootedTreeEdge#postOrderComputeFullDP} has populated logZ tables.
     * Only set in SPARSE mode; controls whether {@link #enhanceWithDPBounds} is applied. */
    private boolean dpTablesReady = false;
    /** RCs used to build the tree / DP tables; stashed for use by enhanceWithDPBounds. */
    private RCs searchRCs;

    // Whether branch decomposition is active
    private boolean useBranchDecomposition = true;

    // Conformation space
    private final SimpleConfSpace confSpace;

    // Branch decomposition state
    private InteractionGraph interactionGraph;
    private BranchDecomposition branchDecomposition;
    private int branchwidth;

    // Rooted tree
    private RootedTreeNode rootedRoot;
    private RootedTreeEdge rootedRootEdge;

    // Statistics
    private int totalEnumerationSteps = 0;
    private int totalMinimizations = 0;
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
    private BigDecimal flatSumZLower = BigDecimal.ZERO;
    private BigDecimal flatSumZUpper = BigDecimal.ZERO;

    // Scorers (created once, reused)
    private AStarScorer gScorerMin;
    private AStarScorer hScorerMin;
    private AStarScorer gScorerRigid;
    private AStarScorer hScorerNegRigid;

    // Set of already-minimized conformations to avoid duplicate work
    private final Set<String> minimizedConfs = new HashSet<>();

    // ========== GNN Pool (Strategy 7-style decoupled GNN for BranchMARK*) ==========

    private GNNConfEnergyCalculator gnnBatchCalc;

    // GNN pool parameters
    private int gpuBatchSize = 1000;    // fire ONNX when pool reaches this size
    private int gnnBudgetMax = 100;     // max simultaneous GNN-bounded nodes
    private double cpQ = 0.06;          // CP quantile (kcal/mol)
    private double cpAlpha = 0.001;
    private double cpDelta = 0.10;

    // GNN pool state (reset per pfunc via compute())
    private final List<DecompSearchNode> gnnPool = new ArrayList<>();
    private final Set<DecompSearchNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private final Map<DecompSearchNode, Double> gnnPredictions = new IdentityHashMap<>();
    private int gnnBudgetUsed = 0;

    // GNN statistics
    private int gnnBounded = 0;
    private int gnnCCDFromGNN = 0;
    private int gnnCCDFromOriginal = 0;
    private double gnnTimeMs = 0;
    private int gnnOnnxCalls = 0;
    private int gnnRoundCounter = 0;

    // ========== Subtree GNN (Strategy 8 for BranchMARK*) ==========

    private GNNSubtreeEnergyCalculator subtreeGNN;

    // Subtree GNN parameters
    private int subtreeBatchSize = 500;
    private double subtreeCpQ = 0.10;      // wider CP band for subtree (less accurate)
    private int subtreeBudgetMax = 50;

    // Subtree GNN pool state
    private final List<DecompSearchNode> subtreePool = new ArrayList<>();
    private final Set<DecompSearchNode> subtreeBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private int subtreeBudgetUsed = 0;

    // Subtree GNN statistics
    private int subtreeGNNBounded = 0;
    private double subtreeGNNTimeMs = 0;
    private int subtreeOnnxCalls = 0;

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

        RootingCandidate(RootedTreeNode root, int splitEdgeIndex, double logTESS,
                         int lambdaEdges, int maxFsetSize, int branchingEdges,
                         int totalFsetEdges, int rootFsetSize) {
            this.root = root;
            this.splitEdgeIndex = splitEdgeIndex;
            this.logTESS = logTESS;
            this.lambdaEdges = lambdaEdges;
            this.maxFsetSize = maxFsetSize;
            this.branchingEdges = branchingEdges;
            this.totalFsetEdges = totalFsetEdges;
            this.rootFsetSize = rootFsetSize;
        }
    }

    // ========== Constructor ==========

    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
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
            System.out.println("BranchMARK*: certifyFullBounds requested but disabled in SPARSE mode "
                    + "(sparse pfunc is a self-consistent Z over C_sparse, not an approximation of full Z).");
        }
        this.certifyFullBounds = requestedCertify && !useSparsePfunc;
        this.edgeSelectionStrategy = EdgeSelectionStrategy.fromProperty(
                getConfigProperty(EDGE_SELECTION_PROPERTY, null));
        this.edgeLookaheadMaxStates = Math.max(1,
                getConfigInteger(EDGE_LOOKAHEAD_MAX_STATES_PROPERTY, 2000));
        this.edgeLookaheadMaxPendingEdges = Math.max(2,
                getConfigInteger(EDGE_LOOKAHEAD_MAX_PENDING_EDGES_PROPERTY, 4));
        this.edgeLookaheadParallel = getConfigBoolean(
                EDGE_LOOKAHEAD_PARALLEL_PROPERTY, true);
        this.rootSplitStrategy = getConfigProperty(ROOT_SPLIT_PROPERTY, "legacy").trim();
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

        System.out.println("BranchMARK*: Cutoff strategy=" + cutoffStrategy
                + ", distCutoff=" + distCutoff
                + ", energyCutoff=" + energyCutoff
                + ", residualBudget=" + residualBudget
                + ", keepConnected=" + cutoffKeepConnected
                + ", certifyFullBounds=" + certifyFullBounds
                + ", energyMode=" + energyMode
                + ", sparsePruneThreshold=" + String.format("%.6f", sparsePruneThreshold));
        System.out.println("BranchMARK*: Edge selection=" + edgeSelectionStrategy
                + ", maxExactLookaheadStates=" + edgeLookaheadMaxStates
                + ", maxPendingLookaheadEdges=" + edgeLookaheadMaxPendingEdges
                + ", parallelLookahead=" + edgeLookaheadParallel
                + ", rootSplit=" + rootSplitStrategy
                + ", rootSplitMaxFset=" + rootSplitMaxFset
                + ", tripleCorrections=" + useHigherOrderCorrections
                + ", correctionAudit=" + correctionAudit
                + ", parallelInternal=" + parallelInternal + " [DEPRECATED]"
                + ", parallelEnumeration=" + parallelEnumeration + " [DEPRECATED]");
        if (trace) {
            System.out.println("BranchMARK*: trace enabled"
                    + " (rounds=" + traceRounds
                    + ", corrections=" + traceCorrectionLimit + ")");
        }
        if (parallelInternal) {
            System.out.println("BranchMARK*: parallelInternal is DEPRECATED; prefer pending-edge lookahead parallelism.");
        }
        if (parallelEnumeration) {
            System.out.println("BranchMARK*: parallelEnumeration is DEPRECATED; prefer pending-edge lookahead parallelism.");
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
            System.out.println("BranchMARK*: Strict SPARSE pfunc active. "
                    + "Cut-edge confs with |pair E| > " + String.format("%.6f", sparsePruneThreshold)
                    + " kcal/mol are pruned out of the conformation set. "
                    + "Accepted cut-edge pair energies are zeroed. "
                    + "Triple corrections / GNN / S8 are disabled.");
            System.out.println("BranchMARK*: Sparse emat stats: cutPositionPairs="
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

        System.out.println("BranchMARK*: Branch decomposition computed. Branchwidth=" + branchwidth
                + ", positions=" + interactionGraph.getNumPositions()
                + ", edges=" + branchDecomposition.getTree().getNumEdges());

        // Step 3: Root the tree. Legacy uses split edge 0; optional strategies
        // can choose a split that exposes multiple pending edges for lookahead.
        RootingCandidate rooting = selectRooting(rcs);
        rootedRoot = rooting != null ? rooting.root : null;
        if (rootedRoot == null) {
            System.out.println("BranchMARK*: Empty tree, falling back to standard MARK*.");
            useBranchDecomposition = false;
            return;
        }
        selectedRootSplitIndex = rooting.splitEdgeIndex;
        System.out.println("BranchMARK*: Root split edge=" + selectedRootSplitIndex
                + " (lambdaEdges=" + rooting.lambdaEdges
                + ", rootFset=" + rooting.rootFsetSize
                + ", maxFset=" + rooting.maxFsetSize
                + ", branchingEdges=" + rooting.branchingEdges
                + ", totalFsetEdges=" + rooting.totalFsetEdges + ")");

        rootedRootEdge = rootedRoot.getLeftChild().getChildOfEdge();

        // Step 4: TESS check
        double logTESS = rootedRootEdge.computeLogTESS();
        double logNaive = computeLogNaive(rcs);
        double tessRatio = Math.exp(logTESS - logNaive);

        System.out.println("BranchMARK*: Compact tree TESS: logTESS=" + String.format("%.2f", logTESS)
                + ", logNaive=" + String.format("%.2f", logNaive)
                + ", ratio=" + String.format("%.4f", tessRatio));

        if (tessRatio > TESS_FALLBACK_THRESHOLD) {
            if (useSparsePfunc) {
                // SPARSE mode uses branch* emat (pruned + zeroed); standard MARK* would
                // sum over the FULL conformation set with the original emat, which is a
                // different pfunc target. Stay on branch decomposition to keep Z = Z_sparse.
                System.out.println("BranchMARK*: TESS ratio " + String.format("%.4f", tessRatio)
                        + " > threshold " + TESS_FALLBACK_THRESHOLD
                        + ", continuing on branch decomposition to keep SPARSE pfunc self-consistent.");
            } else {
                System.out.println("BranchMARK*: TESS ratio " + String.format("%.4f", tessRatio)
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
            System.out.println("BranchMARK*: correction matrix attached"
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
            System.out.println("BranchMARK*: Sparse residual rho=0; branch graph matches the full pair graph.");
            return;
        }

        if (!Double.isFinite(rho)) {
            System.out.println("BranchMARK*: Sparse residual rho is infinite; no finite sparse-to-full Z certificate.");
            return;
        }

        double rt = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
        double log10Factor = rho / rt / Math.log(10.0);
        // Note: this residual bound assumes we sum over the FULL conformation set
        // with the sparse energy (FULL mode). Under SPARSE mode we additionally
        // remove cut-edge confs from C, so this is reported for context only.
        System.out.println("BranchMARK*: Sparse residual certificate: cutEdges="
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

        System.out.println("BranchMARK*: Sparse " + label + " emat: cutPositionPairs="
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

    private static String getConfigProperty(String key, String defaultValue) {
        String value = System.getProperty(key);
        if (value == null) {
            value = System.getProperty("osprey." + key);
        }
        return value != null ? value : defaultValue;
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

    private static double getConfigDouble(String key, double defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("BranchMARK*: Invalid double for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        String value = getConfigProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    private RootingCandidate selectRooting(RCs rcs) {
        int numEdges = branchDecomposition.getTree().getNumEdges();
        if (numEdges == 0) return null;

        String strategy = rootSplitStrategy.toLowerCase(Locale.ROOT);
        if (strategy.isEmpty() || strategy.equals("legacy") || strategy.equals("edge0")) {
            return evaluateRootSplit(rcs, 0, true);
        }

        try {
            int explicitSplit = Integer.parseInt(strategy);
            return evaluateRootSplit(rcs, explicitSplit, true);
        } catch (NumberFormatException ignored) {
            // fall through to named strategies
        }

        if (!strategy.equals("branching") && !strategy.equals("maxfset")
                && !strategy.equals("lookahead")) {
            System.err.println("BranchMARK*: Unknown root split strategy '" + rootSplitStrategy
                    + "', using legacy split edge 0.");
            return evaluateRootSplit(rcs, 0, true);
        }

        double logNaive = computeLogNaive(rcs);
        RootingCandidate best = null;
        for (int splitIdx = 0; splitIdx < numEdges; splitIdx++) {
            RootingCandidate candidate = evaluateRootSplit(rcs, splitIdx, false);
            if (candidate == null) continue;
            if (candidate.maxFsetSize > rootSplitMaxFset) continue;
            if (isBetterRooting(candidate, best, logNaive)) {
                best = candidate;
            }
        }
        return best != null ? evaluateRootSplit(rcs, best.splitEdgeIndex, true) : evaluateRootSplit(rcs, 0, true);
    }

    private RootingCandidate evaluateRootSplit(RCs rcs, int splitEdgeIndex,
                                               boolean initEnumerationArrays) {
        RootedTreeNode root = branchDecomposition.rootBranchTree(rcs, splitEdgeIndex);
        if (root == null) return null;

        RootedTreeEdge.postOrderCompLlambda(root, initEnumerationArrays);
        RootedTreeEdge rootEdge = root.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(root, lambdaEdges);

        int maxFsetSize = 0;
        int branchingEdges = 0;
        int totalFsetEdges = 0;
        for (RootedTreeEdge edge : lambdaEdges) {
            int fsetSize = edge.getFset() == null ? 0 : edge.getFset().size();
            maxFsetSize = Math.max(maxFsetSize, fsetSize);
            if (fsetSize > 1) {
                branchingEdges++;
            }
            totalFsetEdges += fsetSize;
        }

        int rootFsetSize = rootEdge.getFset() == null ? 0 : rootEdge.getFset().size();
        return new RootingCandidate(root, splitEdgeIndex, rootEdge.computeLogTESS(),
                lambdaEdges.size(), maxFsetSize, branchingEdges, totalFsetEdges, rootFsetSize);
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

    private double computeLogNaive(RCs rcs) {
        double logNaive = 0.0;
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            if (rcs.getNum(pos) > 0) {
                logNaive += Math.log(rcs.getNum(pos));
            }
        }
        return logNaive;
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
                    branchRigidEmat, branchMinimizingEmat, interactionGraph, rtForDP);
            RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
            dpTablesReady = true;
            long dpInitTime = System.currentTimeMillis() - dpInitStart;
            // Diagnostic: report root edge DP values so we can spot -inf / 0 issues.
            double[] rootLZL = rootedRootEdge.getLogZLower();
            double[] rootLZU = rootedRootEdge.getLogZUpper();
            String lzlStr = (rootLZL == null || rootLZL.length == 0) ? "null" : String.format("%.4f", rootLZL[0]);
            String lzuStr = (rootLZU == null || rootLZU.length == 0) ? "null" : String.format("%.4f", rootLZU[0]);
            System.out.println("BranchMARK*: DP tables computed in " + dpInitTime + " ms (SPARSE factored Z)"
                    + ", rootEdge logZ[0]: lower=" + lzlStr + ", upper=" + lzuStr
                    + ", rootEdge isLambdaEdge=" + rootedRootEdge.getIsLambdaEdge()
                    + ", mArraySize=" + (rootLZU == null ? -1 : rootLZU.length));
        } else {
            dpTablesReady = false;
            System.out.println("BranchMARK*: DP table precomputation skipped (FULL mode, DP would be invalid)");
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

        // Start with 1 minimization per round, increment like MARK*
        maxMinimizations = 1;

        // Initialize Z bounds from root node
        flatSumZLower = rootSearchNode.subtreeLowerBound;
        flatSumZUpper = rootSearchNode.subtreeUpperBound;

        // Compute initial epsilon
        updateBound();

        // Count lambda-edges and print decomposition stats
        java.util.List<edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge> lambdaEdges = new java.util.ArrayList<>();
        edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);
        int leafEdgeCount = 0;
        long maxLambdaStates = 0;
        for (edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge e : lambdaEdges) {
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
        System.out.println("BranchMARK*: " + lambdaEdges.size() + " lambda-edges ("
                + leafEdgeCount + " leaf), maxLambdaStates=" + maxLambdaStates
                + ", pendingEdges=" + rootSearchNode.pendingEdges.size());

        // [DEBUG] Identify isolated positions (not covered by any lambda-edge)
        java.util.Set<Integer> coveredPos = new java.util.HashSet<>();
        for (edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge e : lambdaEdges) {
            for (int p : e.getLambdaPositionsSorted()) coveredPos.add(p);
            for (int p : e.getMPositionsSorted()) coveredPos.add(p);
        }
        java.util.List<Integer> isolatedPos = new java.util.ArrayList<>();
        for (int p = 0; p < rcs.getNumPos(); p++) {
            if (!coveredPos.contains(p)) isolatedPos.add(p);
        }
        StringBuilder isoSb = new StringBuilder();
        isoSb.append("BranchMARK*: [DEBUG] isolated positions=").append(isolatedPos.size())
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
        degSb.append("BranchMARK*: [DEBUG] sparse-graph degree per position: [");
        for (int p = 0; p < rcs.getNumPos(); p++) {
            if (p > 0) degSb.append(",");
            degSb.append("p").append(p).append("=").append(degree[p]);
        }
        degSb.append("]");
        System.out.println(degSb.toString());

        System.out.println("BranchMARK*: Initialized. "
                + "Root node: upper=" + String.format("%.4e", flatSumZUpper.doubleValue())
                + ", lower=" + String.format("%.4e", flatSumZLower.doubleValue())
                + ", epsilon=" + String.format("%.6f", epsilonBound));
    }

    // ========== GNN configuration methods ==========

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
        if (useSparsePfunc) {
            // GNN predicts full-energy residuals against E_full, not E_sparse,
            // so its CP bounds are not valid corrections to the sparse pfunc.
            System.out.println("[BranchMARK*+GNN] disabled in SPARSE mode "
                    + "(model is not sparse-energy certified).");
            this.gnnBatchCalc = null;
            return;
        }
        this.gnnBatchCalc = calc;
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = size;
        System.out.println("[BranchMARK*+GNN] gpuBatchSize=" + size);
    }

    public void setCPParams(double alpha, double delta, double q) {
        this.cpAlpha = alpha;
        this.cpDelta = delta;
        this.cpQ = q;
        this.gnnBudgetMax = (int) Math.floor(delta / alpha);
        System.out.println("[BranchMARK*+GNN] CP params: alpha=" + alpha + ", delta=" + delta
                + ", q=" + String.format("%.6f", q) + " kcal/mol, budget=" + gnnBudgetMax);
    }

    // ========== Subtree GNN configuration methods ==========

    public void setSubtreeGNN(GNNSubtreeEnergyCalculator calc) {
        if (useSparsePfunc) {
            // Subtree GNN residuals are calibrated against E_full; in SPARSE
            // mode Z lives on E_sparse / C_sparse, so the CP correction is invalid.
            System.out.println("[BranchMARK*+S8] disabled in SPARSE mode "
                    + "(subtree GNN is not sparse-energy certified).");
            this.subtreeGNN = null;
            return;
        }
        this.subtreeGNN = calc;
        System.out.println("[BranchMARK*+S8] Subtree GNN loaded");
    }

    public void setSubtreeBatchSize(int size) {
        this.subtreeBatchSize = size;
    }

    public void setSubtreeCPParams(double q, int budget) {
        this.subtreeCpQ = q;
        this.subtreeBudgetMax = budget;
        System.out.println("[BranchMARK*+S8] Subtree CP: q=" + String.format("%.6f", q)
                + ", budget=" + budget);
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

        System.out.println("BranchMARK*: Numeric audit [" + context + "] "
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

    private static String formatDelta(BigDecimal before, BigDecimal after) {
        return formatBigExp(before.subtract(after));
    }

    private void trace(String message) {
        if (trace) {
            System.out.println("BranchMARK*: [TRACE] " + message);
        }
    }

    private boolean shouldPrintCorrectionDetail(long count) {
        return correctionAudit && (traceCorrectionLimit <= 0 || count <= traceCorrectionLimit);
    }

    private void correctionAuditLog(String message) {
        System.out.println("BranchMARK*: [CORRECTION_AUDIT] " + message);
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
        gnnRoundCounter++;
        if (gnnRoundCounter % 1000 == 0) {
            System.out.println("BranchMARK*: round=" + gnnRoundCounter
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

            // Higher-order corrections are disabled for BranchMARK* so the pfunc
            // certificate stays on the branch bound model.
            if (useHigherOrderCorrections && !node.isAggregate && !node.minimized) {
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
                                + " round=" + gnnRoundCounter
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
                                + " round=" + gnnRoundCounter
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
                                + " round=" + gnnRoundCounter
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

        if (trace && (gnnRoundCounter <= traceRounds || correctionsThisRound > 0)) {
            trace("round=" + gnnRoundCounter
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
            if (trace && (gnnRoundCounter <= traceRounds || correctionsThisRound > 0)) {
                trace("round=" + gnnRoundCounter + " decision=leaf"
                        + " leafZ=" + formatBigExp(leafZ)
                        + " internalZ=" + formatBigExp(internalZ));
            }
            // Leaf error dominates: minimize leaves, put internals back
            for (DecompSearchNode leaf : leafNodes) {
                minimizeDecompLeaf(leaf);
            }
            decompQueue.addAll(internalNodes);

            // Increment maxMinimizations toward full parallelism (like MARK*)
            if (maxMinimizations < parallelism.numThreads) {
                maxMinimizations++;
            }
        } else {
            if (trace && (gnnRoundCounter <= traceRounds || correctionsThisRound > 0)) {
                trace("round=" + gnnRoundCounter + " decision=internal"
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
            double[] logZL = pe.edge.getLogZLower();
            double[] logZU = pe.edge.getLogZUpper();
            if (logZL == null || logZU == null) return;
            if (pe.mIdx < 0 || pe.mIdx >= logZU.length || pe.mIdx >= logZL.length) return;
            double lzU = logZU[pe.mIdx];
            double lzL = logZL[pe.mIdx];
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
                throw new IllegalStateException("BranchMARK*: missing lambda enumeration chunk result");
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
        int[] fullConf = leaf.partialConf;

        // Check if already minimized
        String confKey = SimpleConfSpace.formatConfRCs(fullConf);
        if (minimizedConfs.contains(confKey)) {
            // Put back with current bounds
            decompQueue.add(leaf);
            return;
        }
        minimizedConfs.add(confKey);

        // Track GNN stats
        if (gnnBoundedNodes.remove(leaf)) {
            gnnBudgetUsed--;
            gnnCCDFromGNN++;
        } else {
            gnnCCDFromOriginal++;
        }

        // Use branch emats so pairwise scoring stays inside the sparse pfunc world.
        double eRigid = computeFullConfPairwiseEnergy(fullConf, branchRigidEmat);
        double eMin = computeFullConfPairwiseEnergy(fullConf, branchMinimizingEmat);
        double oldConfLower = leaf.confLowerBound;
        double oldConfUpper = leaf.confUpperBound;
        double epsilonBefore = computeEpsilon();

        // Minimize against the same pair set used by the bounds. In SPARSE mode
        // we limit ResidueInteractions to graph-kept pairs (cut pairs contribute
        // 0 to E_sparse, and A* has already pruned any conf with a +inf cut pair).
        // In FULL mode we use the standard ConfAnalyzer minimization.
        ConfSearch.ScoredConf scoredConf = new ConfSearch.ScoredConf(fullConf, eMin);
        double eTrue;
        if (useSparsePfunc) {
            eTrue = minimizingEcalc.calcEnergy(new RCTuple(fullConf), makeSparseFullConfInters(fullConf)).energy;
            if (useHigherOrderCorrections) {
                correctionSparseGenerationSkips++;
                if (shouldPrintCorrectionDetail(correctionSparseGenerationSkips)) {
                    correctionAuditLog("NO-GENERATE sparse leaf min#" + (totalMinimizations + 1)
                            + " state=" + stateName
                            + " reason=SPARSE uses sparse interactions and does not run ConfAnalyzer triple corrections"
                            + " conf=" + SimpleConfSpace.formatConfRCs(fullConf));
                }
            }
        } else {
            ConfAnalyzer.ConfAnalysis analysis = new ConfAnalyzer(minimizingEcalc).analyze(scoredConf);
            eTrue = analysis.epmol.energy;
            if (useHigherOrderCorrections) {
                computeTripleCorrections(analysis, scoredConf);
            }
        }
        totalMinimizations++;

        // Update Z bounds: remove old contribution, add new (exact) contribution
        BigDecimal oldLeafLower = leaf.subtreeLowerBound;
        BigDecimal oldLeafUpper = leaf.subtreeUpperBound;
        flatSumZLower = flatSumZLower.subtract(leaf.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(leaf.subtreeUpperBound);

        BigDecimal boltzTrue = bc.calc(eTrue);
        // After minimization, both lower and upper for this conformation become exact
        flatSumZLower = flatSumZLower.add(boltzTrue);
        flatSumZUpper = flatSumZUpper.add(boltzTrue);
        double epsilonAfter = computeEpsilon();

        // Mark as minimized (don't re-add to queue — it's exact now)
        leaf.minimized = true;
        leaf.minimizedEnergy = eTrue;

        recordLeafMinimizationProfile("BranchMARK*-" + edgeSelectionStrategy.name(),
                totalMinimizations, fullConf,
                oldConfLower, oldConfUpper, eTrue,
                oldLeafLower, oldLeafUpper, boltzTrue,
                epsilonBefore, epsilonAfter);

        if (totalMinimizations % 100 == 0 || totalMinimizations <= 5) {
            System.out.println("BranchMARK*: " + totalMinimizations
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

    // ========== GNN Pool: batch inference + CP bounding ==========

    /**
     * Flush the GNN pool: one big ONNX batch call, then apply CP bounds
     * to the lowest-energy predictions within the budget.
     */
    private void flushGNNPool() {
        if (gnnPool.isEmpty()) return;

        // Remove already-bounded nodes
        List<DecompSearchNode> candidates = new ArrayList<>();
        for (DecompSearchNode node : gnnPool) {
            if (!gnnBoundedNodes.contains(node) && !node.minimized) {
                candidates.add(node);
            }
        }
        gnnPool.clear();

        if (candidates.isEmpty()) return;

        // ONE big ONNX batch call
        int[][] confs = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            confs[i] = candidates.get(i).partialConf;
        }
        double[] gnnEnergies = gnnBatchCalc.calcEnergies(confs);
        gnnOnnxCalls++;

        // Sort by energy ascending (lowest energy = most important for Z)
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> Double.compare(gnnEnergies[a], gnnEnergies[b]));

        // Apply CP bounds to top budgetMax lowest-energy nodes
        int bounded = 0, skipped = 0;
        int slotsAvailable = gnnBudgetMax - gnnBudgetUsed;

        for (int rank = 0; rank < indices.length && bounded < slotsAvailable; rank++) {
            int idx = indices[rank];
            DecompSearchNode leaf = candidates.get(idx);
            double eGNN = gnnEnergies[idx];

            // Compute current pairwise bounds for this conf (branch emats keep us inside C_sparse)
            double eMin = computeFullConfPairwiseEnergy(leaf.partialConf, branchMinimizingEmat);
            double eRigid = computeFullConfPairwiseEnergy(leaf.partialConf, branchRigidEmat);

            // CP bounds: [eGNN - cpQ, eGNN + cpQ] intersected with [eMin, eRigid]
            double newLower = Math.max(eGNN - cpQ, eMin);
            double newUpper = Math.min(eGNN + cpQ, eRigid);

            if (newLower > newUpper) {
                skipped++;
                continue;
            }

            // Update Z bounds: remove old, add new tighter bounds
            flatSumZLower = flatSumZLower.subtract(leaf.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.subtract(leaf.subtreeUpperBound);

            BigDecimal newBoltzLower = bc.calc(newUpper);  // higher energy → lower Boltzmann
            BigDecimal newBoltzUpper = bc.calc(newLower);  // lower energy → higher Boltzmann
            leaf.subtreeLowerBound = newBoltzLower;
            leaf.subtreeUpperBound = newBoltzUpper;
            leaf.errorBound = newBoltzUpper.subtract(newBoltzLower);

            flatSumZLower = flatSumZLower.add(newBoltzLower);
            flatSumZUpper = flatSumZUpper.add(newBoltzUpper);

            // Re-add to queue with updated bounds
            if (leaf.errorBound.signum() > 0) {
                decompQueue.add(leaf);
            }

            gnnBoundedNodes.add(leaf);
            gnnPredictions.put(leaf, eGNN);
            gnnBudgetUsed++;
            gnnBounded++;
            bounded++;
        }

        System.out.println(String.format("[BranchMARK*+GNN r%d] GNN flush: candidates=%d bounded=%d skipped=%d budget=%d/%d",
                gnnRoundCounter, candidates.size(), bounded, skipped, gnnBudgetUsed, gnnBudgetMax));
    }

    /**
     * Flush the subtree GNN pool: predict mean residuals for internal nodes
     * and tighten their Z bounds using Conformal Prediction.
     *
     * Same logic as MARKStarBoundGNNS8.flushSubtreePool, adapted for DecompSearchNode.
     */
    private void flushSubtreePool() {
        if (subtreePool.isEmpty()) return;

        List<DecompSearchNode> candidates = new ArrayList<>();
        for (DecompSearchNode node : subtreePool) {
            if (!subtreeBoundedNodes.contains(node) && !node.isLeaf()) {
                candidates.add(node);
            }
        }
        subtreePool.clear();
        if (candidates.isEmpty()) return;

        // Build assignments array for subtree GNN
        int[][] assignments = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            assignments[i] = candidates.get(i).partialConf.clone();
        }

        // Predict mean residuals
        double[] meanResiduals = subtreeGNN.predictSubtreeResiduals(assignments);
        subtreeOnnxCalls++;

        // The GNN predicts correction c such that logZ_true ≈ logZ_emat - c/kT.
        // CP bound: |c_pred - c_true| ≤ subtreeCpQ
        final double kT = 0.5922;  // kcal/mol at 298K
        int bounded = 0;
        int slotsAvailable = subtreeBudgetMax - subtreeBudgetUsed;

        // Sort by error bound descending (tighten the loosest bounds first)
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> candidates.get(b).errorBound.compareTo(candidates.get(a).errorBound));

        for (int rank = 0; rank < indices.length && bounded < slotsAvailable; rank++) {
            int idx = indices[rank];
            DecompSearchNode node = candidates.get(idx);
            double c = meanResiduals[idx];

            BigDecimal oldZUpper = node.subtreeUpperBound;
            BigDecimal oldZLower = node.subtreeLowerBound;

            // Apply correction with CP bound
            double scaleUpper = Math.exp(-(c - subtreeCpQ) / kT);
            double scaleLower = Math.exp(-(c + subtreeCpQ) / kT);

            BigDecimal newZUpper = oldZUpper.multiply(new BigDecimal(scaleUpper));
            BigDecimal newZLower = oldZUpper.multiply(new BigDecimal(scaleLower));

            // Only tighten
            boolean tightenedUpper = newZUpper.compareTo(oldZUpper) < 0;
            boolean tightenedLower = newZLower.compareTo(oldZLower) > 0;
            if (!tightenedUpper && !tightenedLower) continue;

            // Sanity: new bounds must not cross
            BigDecimal finalLower = tightenedLower ? newZLower : oldZLower;
            BigDecimal finalUpper = tightenedUpper ? newZUpper : oldZUpper;
            if (finalLower.compareTo(finalUpper) > 0) continue;

            // Update Z totals
            flatSumZLower = flatSumZLower.subtract(oldZLower);
            flatSumZUpper = flatSumZUpper.subtract(oldZUpper);

            node.subtreeLowerBound = finalLower;
            node.subtreeUpperBound = finalUpper;
            node.errorBound = finalUpper.subtract(finalLower);

            flatSumZLower = flatSumZLower.add(finalLower);
            flatSumZUpper = flatSumZUpper.add(finalUpper);

            subtreeBoundedNodes.add(node);
            subtreeBudgetUsed++;
            subtreeGNNBounded++;
            bounded++;
        }

        if (gnnRoundCounter % 50 == 0 || bounded > 0) {
            System.out.println(String.format(
                "[BranchMARK*+S8 r%d] Subtree flush: candidates=%d bounded=%d budget=%d/%d",
                gnnRoundCounter, candidates.size(), bounded, subtreeBudgetUsed, subtreeBudgetMax));
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
                        synchronized (BranchMARKStarBound.this) {
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
                        synchronized (BranchMARKStarBound.this) {
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

        System.out.println("BranchMARK*: Edge selection stats: "
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

        System.out.println("BranchMARK*: Correction audit summary: "
                + "enabled=" + useHigherOrderCorrections
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
        if (useHigherOrderCorrections && useSparsePfunc && correctionSparseGenerationSkips > 0) {
            System.out.println("BranchMARK*: Correction audit note: SPARSE mode did not generate triple corrections; "
                    + "any applications would have to come from pre-existing correctionMatrix entries.");
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
        if (!useBranchDecomposition) {
            super.compute(maxNumConfs);
            return;
        }

        double lastEps = 1.0;

        while (epsilonBound > targetEpsilon) {
            tightenBoundInPhases();

            if (lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("BranchMARK*: Warning - bounds got looser. eps=" + epsilonBound);
            }
            lastEps = epsilonBound;
        }

        // Flush remaining GNN pool at end of pfunc
        if (gnnBatchCalc != null && !gnnPool.isEmpty()) {
            long gnnStart = System.nanoTime();
            flushGNNPool();
            gnnTimeMs += (System.nanoTime() - gnnStart) / 1e6;
            updateBound();
        }

        // Print GNN stats if used
        if (gnnBatchCalc != null) {
            System.out.println("[BranchMARK*+GNN] Done. eps=" + String.format("%.6f", epsilonBound)
                    + ", GNNbounded=" + gnnBounded
                    + ", CCD(fromGNN)=" + gnnCCDFromGNN
                    + ", CCD(fromOrig)=" + gnnCCDFromOriginal
                    + ", budget=" + gnnBudgetUsed + "/" + gnnBudgetMax
                    + ", gnnTime=" + String.format("%.1fms", gnnTimeMs)
                    + ", onnxCalls=" + gnnOnnxCalls);
        }

        // Print subtree GNN stats if used
        if (subtreeGNN != null) {
            System.out.println("[BranchMARK*+S8] Subtree GNN: bounded=" + subtreeGNNBounded
                    + ", time=" + String.format("%.1fms", subtreeGNNTimeMs)
                    + ", onnxCalls=" + subtreeOnnxCalls);
        }

        // Reset GNN per-pfunc state
        gnnPool.clear();
        gnnBoundedNodes.clear();
        gnnPredictions.clear();
        gnnBudgetUsed = 0;
        gnnBounded = 0;
        gnnCCDFromGNN = 0;
        gnnCCDFromOriginal = 0;
        gnnTimeMs = 0;
        gnnOnnxCalls = 0;
        gnnRoundCounter = 0;

        // Reset subtree GNN per-pfunc state
        subtreePool.clear();
        subtreeBoundedNodes.clear();
        subtreeBudgetUsed = 0;
        subtreeGNNBounded = 0;
        subtreeGNNTimeMs = 0;
        subtreeOnnxCalls = 0;

        System.out.println("BranchMARK*: Finished. epsilon=" + String.format("%.6f", epsilonBound)
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
            System.out.println("BranchMARK*: Sparse-to-full certificate epsilon="
                    + String.format("%.6f", reportedEpsilon)
                    + " > target=" + String.format("%.6f", targetEpsilon)
                    + "; leaving status as " + getStatus() + ".");
        }

        if (reportedEpsilon <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
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

        System.out.println("BranchMARK*: Applied sparse-to-full pfunc inflation: rho="
                + String.format("%.6f", rho)
                + ", lowerFactor=" + String.format("%.6e", lowerFactor.doubleValue())
                + ", upperFactor=" + String.format("%.6e", upperFactor.doubleValue())
                + ", sparseEpsilon=" + String.format("%.6f", epsilonBound)
                + ", certifiedEpsilon=" + String.format("%.6f", certifiedEpsilon));

        return certifiedEpsilon;
    }
}
