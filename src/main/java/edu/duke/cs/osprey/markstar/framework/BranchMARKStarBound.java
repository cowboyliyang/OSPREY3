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
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.branch.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.*;

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
    private static final double TESS_FALLBACK_THRESHOLD = 0.5;  // fallback if TESS/Naive > this

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

    // ========== Constructor ==========

    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        this.confSpace = confSpace;

        // Step 1: Build sparse interaction graph using BWM-style dual cutoff
        interactionGraph = InteractionGraph.buildWithDualCutoff(
                confSpace, rigidEmat, minimizingEmat, rcs,
                DEFAULT_DIST_CUTOFF, DEFAULT_ENERGY_CUTOFF);

        // Step 2: Compute branch decomposition
        branchDecomposition = new BranchDecomposition(interactionGraph);
        branchDecomposition.compute();
        branchwidth = branchDecomposition.getBranchwidth();

        System.out.println("BranchMARK*: Branch decomposition computed. Branchwidth=" + branchwidth
                + ", positions=" + interactionGraph.getNumPositions()
                + ", edges=" + branchDecomposition.getTree().getNumEdges());

        // Step 3: Root the tree
        rootedRoot = branchDecomposition.rootBranchTree(rcs);
        if (rootedRoot == null) {
            System.out.println("BranchMARK*: Empty tree, falling back to standard MARK*.");
            useBranchDecomposition = false;
            return;
        }

        // Step 4: Post-order traversal to compute L/lambda/F-sets
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);

        rootedRootEdge = rootedRoot.getLeftChild().getChildOfEdge();

        // Step 5: Compact the tree
        rootedRootEdge.compactTree();

        // Step 6: TESS check
        double logTESS = rootedRootEdge.computeLogTESS();
        double logNaive = 0.0;
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            if (rcs.getNum(pos) > 0) {
                logNaive += Math.log(rcs.getNum(pos));
            }
        }
        double tessRatio = Math.exp(logTESS - logNaive);

        System.out.println("BranchMARK*: Compact tree TESS: logTESS=" + String.format("%.2f", logTESS)
                + ", logNaive=" + String.format("%.2f", logNaive)
                + ", ratio=" + String.format("%.4f", tessRatio));

        if (tessRatio > TESS_FALLBACK_THRESHOLD) {
            System.out.println("BranchMARK*: TESS ratio " + String.format("%.4f", tessRatio)
                    + " > threshold " + TESS_FALLBACK_THRESHOLD + ", falling back to standard MARK*.");
            useBranchDecomposition = false;
            return;
        }

        // Initialize scorers and root search node
        initSearch(rcs);
    }

    /** Initialize search: create scorers and root search node */
    private void initSearch(RCs rcs) {
        int numPos = rcs.getNumPos();

        // Create scorers (same types as parent MARKStarBound uses)
        gScorerMin = new PairwiseGScorer(minimizingEmat);
        hScorerMin = new MPLPPairwiseHScorer(new EdgeUpdater(), minimizingEmat, 1, 0.0001);
        gScorerRigid = new PairwiseGScorer(rigidEmat);
        hScorerNegRigid = new TraditionalPairwiseHScorer(
                new NegatedEnergyMatrix(confSpace, rigidEmat), rcs);

        // DP table precomputation DISABLED — enhanceWithDPBounds is also disabled
        // (pruned interaction graph misses cut pairwise interactions, producing invalid Z bounds).
        // Skipping saves ~35s of startup time across all sequences.
        // long dpInitStart = System.currentTimeMillis();
        // RootedTreeEdge.postOrderInitIncremental(rootedRoot, rigidEmat, minimizingEmat,
        //         interactionGraph, RT);
        // RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
        // long dpInitTime = System.currentTimeMillis() - dpInitStart;
        // System.out.println("BranchMARK*: DP tables computed in " + dpInitTime + " ms");
        System.out.println("BranchMARK*: DP table precomputation skipped (disabled)");

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
            sb.append(", leaf=").append(!e.hasFsetChildren());
            System.out.println(sb.toString());
        }
        System.out.println("BranchMARK*: " + lambdaEdges.size() + " lambda-edges ("
                + leafEdgeCount + " leaf), maxLambdaStates=" + maxLambdaStates
                + ", pendingEdges=" + rootSearchNode.pendingEdges.size());

        System.out.println("BranchMARK*: Initialized. "
                + "Root node: upper=" + String.format("%.4e", flatSumZUpper.doubleValue())
                + ", lower=" + String.format("%.4e", flatSumZLower.doubleValue())
                + ", epsilon=" + String.format("%.6f", epsilonBound));
    }

    // ========== GNN configuration methods ==========

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
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
        if (flatSumZUpper.signum() <= 0) return Double.POSITIVE_INFINITY;
        if (flatSumZLower.signum() <= 0) return 1.0;

        double eps = MathTools.bigDivide(flatSumZUpper.subtract(flatSumZLower), flatSumZUpper,
                PartitionFunction.decimalPrecision).doubleValue();
        return Math.max(0.0, eps);
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
                    + ", drillDowns=" + totalDrillDowns);
        }

        // === Phase 1: Pull nodes (mirrors MARK*'s populateQueues) ===
        // Pull up to maxNodes internals and maxMinimizations leaves.
        // Corrected nodes are put back without classification.
        int maxNodes = 1000;
        List<DecompSearchNode> internalNodes = new ArrayList<>();
        List<DecompSearchNode> leafNodes = new ArrayList<>();
        List<DecompSearchNode> leftoverLeaves = new ArrayList<>();
        BigDecimal internalZ = BigDecimal.ZERO;
        BigDecimal leafZ = BigDecimal.ZERO;

        while (!decompQueue.isEmpty() && internalNodes.size() < maxNodes) {
            DecompSearchNode node = decompQueue.poll();

            if (node.errorBound.signum() <= 0) continue;

            // Apply correction matrix: triple corrections from previous minimizations
            // may tighten confLowerBound for any node, not just the minimized one.
            // Corrected nodes go back to queue (same as MARK*'s populateQueues).
            if (!node.isAggregate && !node.minimized) {
                double confCorrection = correctionMatrix.confE(node.partialConf);
                if (confCorrection > node.confLowerBound) {
                    flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
                    flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

                    node.confLowerBound = Math.min(confCorrection, node.confUpperBound);
                    node.recomputeZBounds(bc);

                    flatSumZLower = flatSumZLower.add(node.subtreeLowerBound);
                    flatSumZUpper = flatSumZUpper.add(node.subtreeUpperBound);

                    leftoverLeaves.add(node);
                    continue;
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

        // === Phase 2: Error-driven selection (mirrors MARK*'s decision logic) ===
        // Compare aggregate leaf error vs internal error; process only the dominant type.
        BigDecimal drillDownThreshold = new BigDecimal(1 - targetEpsilon);

        if (MathTools.isLessThan(internalZ, leafZ)) {
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
            // Internal error dominates: expand/drill-down internals, put leaves back
            List<DecompSearchNode> newNodes = new ArrayList<>();
            for (DecompSearchNode internal : internalNodes) {
                // High error check (same criterion as MARK*'s boundLowestBoundConfUnderNode):
                // lowerBound <= 1 AND upperBound/rootUpperBound > (1 - targetEpsilon)
                if (!MathTools.isGreaterThan(internal.subtreeLowerBound, BigDecimal.ONE) &&
                    MathTools.isGreaterThan(
                        MathTools.bigDivide(internal.subtreeUpperBound, flatSumZUpper,
                                PartitionFunction.decimalPrecision),
                        drillDownThreshold)) {
                    // High error: drill-down through all pending edges with full enumeration.
                    // All siblings become explicit children (no aggregates created).
                    drillDownFullEnum(internal, newNodes);
                } else {
                    // Low error: expand one pending edge with full enumeration (no aggregates)
                    expandDecompNode(internal);
                }
            }
            for (DecompSearchNode n : newNodes) {
                if (n.errorBound.signum() > 0) {
                    decompQueue.add(n);
                }
            }
            decompQueue.addAll(leafNodes);
        }

        updateBound();
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
        // DISABLED: DP table uses pruned interaction graph, missing cut pairwise
        // interactions. This can produce Z bounds that are NOT valid upper/lower bounds
        // for the full energy model, causing systematic underestimation of qstar.
        // TODO: fix by adding correction for cut interactions, or build decomposition
        // on the full interaction graph.
        return;
    }

    /** Number of lambda-states to produce per expand (top-k by energy). */
    private static final int LAMBDA_EXPAND_K = Integer.MAX_VALUE;

    /**
     * Expand a decomposition search node: assign one pending edge's lambda-set.
     *
     * Always uses LambdaAStar to lazily produce top-k lowest-energy lambda-states,
     * plus an aggregate node representing all remaining states. Never enumerates
     * the full lambda-set at once.
     */
    private void expandDecompNode(DecompSearchNode node) {
        if (node.isAggregate) {
            expandAggregateNode(node);
            return;
        }

        // Pick the pending edge to expand
        int pendingIdx = node.pickHighestErrorPendingEdge();
        DecompSearchNode.PendingEdge pe = node.pendingEdges.get(pendingIdx);
        RootedTreeEdge edge = pe.edge;

        // Remove parent's contribution from running Z totals
        flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

        totalEnumerationSteps++;

        // Lazy A* enumeration: produce top-k, aggregate the rest
        LambdaAStar lambdaAStar = new LambdaAStar(edge, pe.mIdx, minimizingEmat);

        int k = Math.min(LAMBDA_EXPAND_K, edge.getTotalLambdaStates());
        BigDecimal childrenZUpper = BigDecimal.ZERO;

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
     */
    private void drillDownFullEnum(DecompSearchNode node, List<DecompSearchNode> newNodes) {
        totalDrillDowns++;

        // Remove the original node's contribution (replaced by its children)
        flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

        DecompSearchNode current = node;

        while (!current.isLeaf()) {
            int pendingIdx = current.pickHighestErrorPendingEdge();
            DecompSearchNode.PendingEdge pe = current.pendingEdges.get(pendingIdx);
            RootedTreeEdge edge = pe.edge;

            totalEnumerationSteps++;

            // Enumerate ALL lambda-states for this edge (no aggregate)
            LambdaAStar lambdaAStar = new LambdaAStar(edge, pe.mIdx, minimizingEmat);
            DecompSearchNode bestChild = null;
            double bestChildLower = Double.POSITIVE_INFINITY;

            while (true) {
                int lambdaIdx = lambdaAStar.next();
                if (lambdaIdx < 0) break;

                DecompSearchNode child = DecompSearchNode.makeChild(
                        current, pendingIdx, lambdaIdx, RCs,
                        gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

                // Tighten Z bounds using DP table lookups for remaining pending edges
                enhanceWithDPBounds(child);

                // Add child to flat sum and newNodes
                flatSumZLower = flatSumZLower.add(child.subtreeLowerBound);
                flatSumZUpper = flatSumZUpper.add(child.subtreeUpperBound);
                newNodes.add(child);

                // Track best child (lowest confLowerBound) for continued drilling
                if (child.confLowerBound < bestChildLower) {
                    bestChild = child;
                    bestChildLower = child.confLowerBound;
                }
            }

            if (bestChild == null) return; // all pruned

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

        double eRigid = computeFullConfPairwiseEnergy(fullConf, rigidEmat);
        double eMin = computeFullConfPairwiseEnergy(fullConf, minimizingEmat);

        // Full minimization
        ConfSearch.ScoredConf scoredConf = new ConfSearch.ScoredConf(fullConf, eMin);
        ConfAnalyzer.ConfAnalysis analysis = new ConfAnalyzer(minimizingEcalc).analyze(scoredConf);
        double eTrue = analysis.epmol.energy;
        totalMinimizations++;

        // Compute triple corrections
        computeTripleCorrections(analysis, scoredConf);

        // Update Z bounds: remove old contribution, add new (exact) contribution
        BigDecimal oldLeafLower = leaf.subtreeLowerBound;
        BigDecimal oldLeafUpper = leaf.subtreeUpperBound;
        flatSumZLower = flatSumZLower.subtract(leaf.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(leaf.subtreeUpperBound);

        BigDecimal boltzTrue = bc.calc(eTrue);
        // After minimization, both lower and upper for this conformation become exact
        flatSumZLower = flatSumZLower.add(boltzTrue);
        flatSumZUpper = flatSumZUpper.add(boltzTrue);

        // Mark as minimized (don't re-add to queue — it's exact now)
        leaf.minimized = true;
        leaf.minimizedEnergy = eTrue;

        if (totalMinimizations % 100 == 0 || totalMinimizations <= 5) {
            System.out.println("BranchMARK*: " + totalMinimizations
                    + " minimizations, epsilon=" + String.format("%.6f", epsilonBound)
                    + ", flatSumZLower=" + String.format("%.6e", flatSumZLower.doubleValue())
                    + ", flatSumZUpper=" + String.format("%.6e", flatSumZUpper.doubleValue())
                    + ", eTrue=" + String.format("%.4f", eTrue)
                    + ", boltzTrue=" + String.format("%.6e", boltzTrue.doubleValue())
                    + ", leafOldLower=" + String.format("%.6e", oldLeafLower.doubleValue())
                    + ", leafOldUpper=" + String.format("%.6e", oldLeafUpper.doubleValue()));
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
        double eMin = computeFullConfPairwiseEnergy(leaf.partialConf, minimizingEmat);
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

            // Compute current pairwise bounds for this conf
            double eMin = computeFullConfPairwiseEnergy(leaf.partialConf, minimizingEmat);
            double eRigid = computeFullConfPairwiseEnergy(leaf.partialConf, rigidEmat);

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

        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(minimizingEmat);
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

        double tripleThreshold = 0.3;
        for (int[] pair : largePairs) {
            int pos1 = pair[0], pos2 = pair[1];
            for (int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos1 || pos3 == pos2) continue;
                RCTuple tuple = new RCTuple(
                        pos1, conf.getAssignments()[pos1],
                        pos2, conf.getAssignments()[pos2],
                        pos3, conf.getAssignments()[pos3]);
                double tupleBounds = rigidEmat.getInternalEnergy(tuple)
                        - minimizingEmat.getInternalEnergy(tuple);
                if (tupleBounds < tripleThreshold) continue;
                if (correctionMatrix.hasHigherOrderTermFor(tuple)) continue;

                minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
                    double tripleEnergy = minimizedTuple.energy;
                    double lowerbound = minimizingEmat.getInternalEnergy(tuple);
                    double correction = tripleEnergy - lowerbound;
                    if (correction != 0) {
                        correctionMatrix.setHigherOrder(tuple, correction);
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

    // ========== Override compute ==========

    @Override
    protected void updateBound() {
        if (useBranchDecomposition) {
            epsilonBound = computeEpsilon();
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
                + totalDrillDowns + " drill-downs.");

        // Set final partition function values
        PartitionFunction.Values vals = getValues();
        vals.qstar = flatSumZLower;
        vals.pstar = flatSumZUpper;
        vals.qprime = vals.pstar.subtract(vals.qstar);

        if (epsilonBound <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
    }
}
