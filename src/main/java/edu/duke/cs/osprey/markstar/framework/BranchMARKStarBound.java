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
 * Supports two compute modes:
 *
 * SUM_OF_PRODUCTS (original):
 *   Phase 1: One-shot bottom-up DP on branch decomposition tree.
 *   Phase 2: Full minimization with Z corrections.
 *
 * FLAT_SUM (new default):
 *   Top-down edge-level traversal of the decomposition tree.
 *   Each "expand" assigns an entire edge's lambda-set.
 *   MARK*-style priority queue decides between expand (internal) and minimize (leaf).
 *   Z bounds maintained incrementally from h-scores (no DP needed).
 *   O(1) bound updates per operation.
 */
public class BranchMARKStarBound extends MARKStarBound {

    // ========== Compute mode ==========

    public enum ComputeMode {
        /** Bottom-up DP on branch decomposition tree (original) */
        SUM_OF_PRODUCTS,
        /** Top-down decomposition-guided A* with flat sum updates (new default) */
        FLAT_SUM
    }

    // Configuration
    private static final double DEFAULT_DIST_CUTOFF = 8.0;      // Angstroms
    private static final double DEFAULT_ENERGY_CUTOFF = 0.1;    // kcal/mol
    private static final double TESS_FALLBACK_THRESHOLD = 0.5;  // fallback if TESS/Naive > this

    // Whether branch decomposition is active
    private boolean useBranchDecomposition = true;

    // Compute mode
    private final ComputeMode computeMode;

    // Conformation space
    private final SimpleConfSpace confSpace;

    // Branch decomposition state
    private InteractionGraph interactionGraph;
    private BranchDecomposition branchDecomposition;
    private int branchwidth;

    // Rooted tree
    private RootedTreeNode rootedRoot;
    private RootedTreeEdge rootedRootEdge;

    // All lambda-edges for scheduling (SUM_OF_PRODUCTS mode)
    private List<RootedTreeEdge> allLambdaEdges;

    // Log-space DP partition function bounds (SUM_OF_PRODUCTS mode)
    private double logDpLower = Double.NEGATIVE_INFINITY;
    private double logDpUpper = Double.NEGATIVE_INFINITY;

    // Phase 2: Z correction accumulators (SUM_OF_PRODUCTS mode)
    private BigDecimal corrLower = BigDecimal.ZERO;
    private BigDecimal corrUpper = BigDecimal.ZERO;

    // RT = R * T (gas constant x temperature)
    private final double RT;

    // Statistics
    private int totalEnumerationSteps = 0;
    private int totalMinimizations = 0;

    // Phase tracking (SUM_OF_PRODUCTS mode)
    private boolean phase1Complete = false;

    // ========== FLAT_SUM mode fields ==========

    // Priority queue for decomposition-guided search
    private PriorityQueue<DecompSearchNode> decompQueue;

    // Running Z bounds (FLAT_SUM mode)
    // Z_lower = sum of all nodes' subtreeLowerBound
    // Z_upper = sum of all nodes' subtreeUpperBound
    private BigDecimal flatSumZLower = BigDecimal.ZERO;
    private BigDecimal flatSumZUpper = BigDecimal.ZERO;

    // Scorers for FLAT_SUM mode (created once, reused)
    private AStarScorer gScorerMin;
    private AStarScorer hScorerMin;
    private AStarScorer gScorerRigid;
    private AStarScorer hScorerNegRigid;

    // Set of already-minimized conformations to avoid duplicate work
    private final Set<String> minimizedConfs = new HashSet<>();

    // ========== Constructors ==========

    /** Original constructor: SUM_OF_PRODUCTS mode */
    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism) {
        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism,
                ComputeMode.SUM_OF_PRODUCTS);
    }

    /** Constructor with explicit compute mode */
    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism, ComputeMode computeMode) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        this.confSpace = confSpace;
        this.computeMode = computeMode;
        this.RT = BoltzmannCalculator.RClassic * 298.0;

        // Step 1: Build sparse interaction graph using BWM-style dual cutoff
        interactionGraph = InteractionGraph.buildWithDualCutoff(
                confSpace, rigidEmat, minimizingEmat, rcs,
                DEFAULT_DIST_CUTOFF, DEFAULT_ENERGY_CUTOFF);

        // Step 2: Compute branch decomposition
        branchDecomposition = new BranchDecomposition(interactionGraph);
        branchDecomposition.compute();
        branchwidth = branchDecomposition.getBranchwidth();

        System.out.println("BranchMARK* [" + computeMode + "]: Branch decomposition computed. Branchwidth=" + branchwidth
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

        // Mode-specific initialization
        if (computeMode == ComputeMode.FLAT_SUM) {
            initFlatSumMode(rcs);
        } else {
            initSumOfProductsMode(rcs);
        }
    }

    /** Initialize SUM_OF_PRODUCTS mode: run full DP */
    private void initSumOfProductsMode(RCs rcs) {
        // Step 7: Initialize incremental enumeration (k=0, loose bounds)
        long initStart = System.currentTimeMillis();
        RootedTreeEdge.postOrderInitIncremental(rootedRoot, rigidEmat, minimizingEmat,
                interactionGraph, RT);
        long initTime = System.currentTimeMillis() - initStart;

        // Step 8: Collect all lambda-edges
        allLambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, allLambdaEdges);

        int leafEdgeCount = 0;
        for (RootedTreeEdge edge : allLambdaEdges) {
            if (!edge.hasFsetChildren()) leafEdgeCount++;
        }

        // Step 9: One-shot bottom-up DP
        long dpStart = System.currentTimeMillis();
        RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
        long dpTime = System.currentTimeMillis() - dpStart;

        // Get final bounds from root edge
        refreshGlobalBounds();
        updateBound();
        phase1Complete = true;

        System.out.println("BranchMARK*: One-shot DP complete in " + dpTime + " ms. "
                + allLambdaEdges.size() + " lambda-edges (" + leafEdgeCount + " leaf), "
                + "epsilon=" + String.format("%.6f", epsilonBound)
                + " (init=" + initTime + " ms)");
    }

    /** Initialize FLAT_SUM mode: create scorers and root search node */
    private void initFlatSumMode(RCs rcs) {
        int numPos = rcs.getNumPos();

        // Create scorers (same types as parent MARKStarBound uses)
        gScorerMin = new PairwiseGScorer(minimizingEmat);
        hScorerMin = new MPLPPairwiseHScorer(new EdgeUpdater(), minimizingEmat, 1, 0.0001);
        gScorerRigid = new PairwiseGScorer(rigidEmat);
        hScorerNegRigid = new TraditionalPairwiseHScorer(
                new NegatedEnergyMatrix(confSpace, rigidEmat), rcs);

        // Create root search node
        DecompSearchNode rootSearchNode = DecompSearchNode.makeRoot(
                rootedRootEdge, numPos, rcs,
                gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

        // Initialize priority queue
        decompQueue = new PriorityQueue<>();
        decompQueue.add(rootSearchNode);

        // Initialize Z bounds from root node
        flatSumZLower = rootSearchNode.subtreeLowerBound;
        flatSumZUpper = rootSearchNode.subtreeUpperBound;

        // Compute initial epsilon
        updateBound();

        System.out.println("BranchMARK* [FLAT_SUM]: Initialized. "
                + "Root node: upper=" + String.format("%.4e", flatSumZUpper.doubleValue())
                + ", lower=" + String.format("%.4e", flatSumZLower.doubleValue())
                + ", epsilon=" + String.format("%.6f", epsilonBound));
    }

    // ========== Log-space epsilon computation (SUM_OF_PRODUCTS) ==========

    private double computeEpsilon() {
        if (logDpUpper == Double.NEGATIVE_INFINITY) return Double.POSITIVE_INFINITY;

        if (corrLower.signum() == 0 && corrUpper.signum() == 0) {
            double logRatio = logDpLower - logDpUpper;
            if (logRatio > 0) logRatio = 0;
            return 1.0 - Math.exp(logRatio);
        }

        BigDecimal zLower = logToPositiveBigDecimal(logDpLower).add(corrLower);
        BigDecimal zUpper = logToPositiveBigDecimal(logDpUpper).add(corrUpper);

        if (zUpper.signum() <= 0) return Double.POSITIVE_INFINITY;
        if (zLower.signum() <= 0) return 1.0;

        double eps = MathTools.bigDivide(zUpper.subtract(zLower), zUpper,
                PartitionFunction.decimalPrecision).doubleValue();
        return Math.max(0.0, eps);
    }

    /** Compute epsilon for FLAT_SUM mode */
    private double computeFlatSumEpsilon() {
        if (flatSumZUpper.signum() <= 0) return Double.POSITIVE_INFINITY;
        if (flatSumZLower.signum() <= 0) return 1.0;

        double eps = MathTools.bigDivide(flatSumZUpper.subtract(flatSumZLower), flatSumZUpper,
                PartitionFunction.decimalPrecision).doubleValue();
        return Math.max(0.0, eps);
    }

    private void refreshGlobalBounds() {
        double[] rootLogZLower = rootedRootEdge.getLogZLower();
        double[] rootLogZUpper = rootedRootEdge.getLogZUpper();
        if (rootLogZLower != null && rootLogZLower.length > 0) {
            logDpLower = rootLogZLower[0];
            logDpUpper = rootLogZUpper[0];
        }
    }

    // ========== Bound tightening ==========

    @Override
    protected void tightenBoundInPhases() {
        if (computeMode == ComputeMode.FLAT_SUM && useBranchDecomposition) {
            tightenFlatSum();
            return;
        }

        if (!useBranchDecomposition) {
            super.tightenBoundInPhases();
            return;
        }

        // SUM_OF_PRODUCTS: Phase 2 minimization
        tightenPhase2();
    }

    // ========== FLAT_SUM: Top-down decomposition-guided expand/minimize ==========

    /**
     * One round of the FLAT_SUM tighten loop.
     * Mirrors MARK*'s tightenBoundInPhases: pull nodes, classify, decide, process.
     */
    private void tightenFlatSum() {
        // Pull a batch of nodes from the queue
        int batchSize = Math.max(1, Math.min(maxMinimizations, 100));
        List<DecompSearchNode> internalNodes = new ArrayList<>();
        List<DecompSearchNode> leafNodes = new ArrayList<>();

        BigDecimal internalZ = BigDecimal.ZERO;
        BigDecimal leafZ = BigDecimal.ZERO;

        int pulled = 0;
        while (pulled < batchSize && !decompQueue.isEmpty()) {
            DecompSearchNode node = decompQueue.poll();

            // Skip nodes with negligible error
            if (node.errorBound.signum() <= 0) continue;

            if (node.isLeaf()) {
                // Check correction matrix first (may skip minimization)
                if (!node.minimized) {
                    double confCorrection = correctionMatrix.confE(node.partialConf);
                    double eMin = computeFullConfPairwiseEnergy(node.partialConf, minimizingEmat);
                    if (confCorrection > eMin) {
                        // Correction tightens the lower bound without full minimization
                        applyCorrection(node, confCorrection);
                        continue;
                    }
                }
                leafNodes.add(node);
                leafZ = leafZ.add(node.errorBound);
            } else {
                internalNodes.add(node);
                internalZ = internalZ.add(node.errorBound);
            }
            pulled++;
        }

        // Decision: expand internals or minimize leaves?
        if (MathTools.isLessThan(internalZ, leafZ)) {
            // Minimize leaves (leaf error dominates)
            for (DecompSearchNode leaf : leafNodes) {
                minimizeDecompLeaf(leaf);
            }
            // Put internal nodes back
            decompQueue.addAll(internalNodes);
        } else {
            // Expand internals (internal error dominates)
            for (DecompSearchNode internal : internalNodes) {
                expandDecompNode(internal);
            }
            // Put leaf nodes back
            decompQueue.addAll(leafNodes);
        }

        updateBound();
    }

    /**
     * Expand a decomposition search node: assign one pending edge's lambda-set.
     * Creates children for each lambda-state and updates Z bounds.
     */
    private void expandDecompNode(DecompSearchNode node) {
        // Pick the pending edge to expand
        int pendingIdx = node.pickHighestErrorPendingEdge();
        DecompSearchNode.PendingEdge pe = node.pendingEdges.get(pendingIdx);
        RootedTreeEdge edge = pe.edge;

        // Remove parent's contribution from running Z totals
        flatSumZLower = flatSumZLower.subtract(node.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(node.subtreeUpperBound);

        // Enumerate all lambda-states for this edge at its M-state
        int numLambdaStates = edge.getTotalLambdaStates();
        totalEnumerationSteps++;

        for (int lambdaIdx = 0; lambdaIdx < numLambdaStates; lambdaIdx++) {
            DecompSearchNode child = DecompSearchNode.makeChild(
                    node, pendingIdx, lambdaIdx, RCs,
                    gScorerMin, hScorerMin, gScorerRigid, hScorerNegRigid, bc);

            // Add child's contribution to running Z totals
            flatSumZLower = flatSumZLower.add(child.subtreeLowerBound);
            flatSumZUpper = flatSumZUpper.add(child.subtreeUpperBound);

            // Add to queue if it has meaningful error
            if (child.errorBound.signum() > 0) {
                decompQueue.add(child);
            }
        }
    }

    /**
     * Minimize a leaf node (full conformation).
     * Same pipeline as SUM_OF_PRODUCTS Phase 2: calcEnergy + triple corrections.
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
        flatSumZLower = flatSumZLower.subtract(leaf.subtreeLowerBound);
        flatSumZUpper = flatSumZUpper.subtract(leaf.subtreeUpperBound);

        BigDecimal boltzTrue = bc.calc(eTrue);
        // After minimization, both lower and upper for this conformation become exact
        flatSumZLower = flatSumZLower.add(boltzTrue);
        flatSumZUpper = flatSumZUpper.add(boltzTrue);

        // Mark as minimized (don't re-add to queue — it's exact now)
        leaf.minimized = true;
        leaf.minimizedEnergy = eTrue;

        if (totalMinimizations % 100 == 0) {
            System.out.println("BranchMARK* [FLAT_SUM]: " + totalMinimizations
                    + " minimizations, epsilon=" + String.format("%.6f", epsilonBound));
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

    // ========== Phase 2: Full minimization with Z correction (SUM_OF_PRODUCTS) ==========

    private void tightenPhase2() {
        int batchSize = Math.max(1, Math.min(maxMinimizations, 10));

        for (int i = 0; i < batchSize; i++) {
            int[] fullConf = generateHighErrorConformation();
            if (fullConf == null) break;

            String confKey = SimpleConfSpace.formatConfRCs(fullConf);
            if (minimizedConfs.contains(confKey)) continue;
            minimizedConfs.add(confKey);

            double eRigid = computeFullConfPairwiseEnergy(fullConf, rigidEmat);
            double eMin = computeFullConfPairwiseEnergy(fullConf, minimizingEmat);

            double confCorrection = correctionMatrix.confE(fullConf);
            if (confCorrection > eMin) {
                BigDecimal boltzCorr = bc.calc(confCorrection);
                BigDecimal boltzRigid = bc.calc(eRigid);
                corrLower = corrLower.add(boltzCorr.subtract(boltzRigid));
                continue;
            }

            ConfSearch.ScoredConf scoredConf = new ConfSearch.ScoredConf(fullConf, eMin);
            ConfAnalyzer.ConfAnalysis analysis = new ConfAnalyzer(minimizingEcalc).analyze(scoredConf);
            double eTrue = analysis.epmol.energy;
            totalMinimizations++;

            computeTripleCorrections(analysis, scoredConf);

            BigDecimal boltzTrue = bc.calc(eTrue);
            BigDecimal boltzRigid = bc.calc(eRigid);
            BigDecimal boltzMin = bc.calc(eMin);

            corrLower = corrLower.add(boltzTrue.subtract(boltzRigid));
            corrUpper = corrUpper.add(boltzTrue.subtract(boltzMin));
        }

        updateBound();

        if (totalMinimizations > 0 && totalMinimizations % 100 == 0) {
            System.out.println("BranchMARK*: Phase 2 progress: " + totalMinimizations
                    + " minimizations, epsilon=" + String.format("%.6f", epsilonBound));
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

    // ========== SUM_OF_PRODUCTS: conformation generation ==========

    private int[] generateHighErrorConformation() {
        int numPos = RCs.getNumPos();
        int[] conf = new int[numPos];
        Arrays.fill(conf, -1);
        generateConfRecursive(rootedRootEdge, 0, conf);
        for (int i = 0; i < numPos; i++) {
            if (conf[i] == -1) {
                System.err.println("BranchMARK*: Warning - position " + i + " not assigned in generated conf");
                return null;
            }
        }
        return conf;
    }

    private void generateConfRecursive(RootedTreeEdge edge, int mIdx, int[] conf) {
        int lambdaIdx = edge.getBestErrorLambdaState(mIdx);
        if (lambdaIdx < 0) lambdaIdx = 0;

        int[] mPositions = edge.getMPositionsSorted();
        int[] mGlobalRCs = edge.getMGlobalRCs(mIdx);
        for (int i = 0; i < mPositions.length; i++) {
            conf[mPositions[i]] = mGlobalRCs[i];
        }

        int[] lambdaPositions = edge.getLambdaPositionsSorted();
        int[] lambdaGlobalRCs = edge.getLambdaGlobalRCs(lambdaIdx);
        for (int i = 0; i < lambdaPositions.length; i++) {
            conf[lambdaPositions[i]] = lambdaGlobalRCs[i];
        }

        if (edge.getFset() != null) {
            int[] mRCsLocal = edge.decodeMStatePublic(mIdx);
            int[] lambdaRCsLocal = edge.decodeLambdaStatePublic(lambdaIdx);
            for (RootedTreeEdge fEdge : edge.getFset()) {
                int[] fMRCs = edge.getMstateForFullState(mRCsLocal, lambdaRCsLocal, fEdge);
                int fIndex = fEdge.computeIndexInA(fMRCs);
                generateConfRecursive(fEdge, fIndex, conf);
            }
        }
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
        if (computeMode == ComputeMode.FLAT_SUM && useBranchDecomposition) {
            epsilonBound = computeFlatSumEpsilon();
            return;
        }

        if (!useBranchDecomposition || logDpUpper == Double.NEGATIVE_INFINITY) {
            super.updateBound();
            return;
        }

        epsilonBound = computeEpsilon();
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

        System.out.println("BranchMARK* [" + computeMode + "]: Finished. epsilon=" + String.format("%.6f", epsilonBound)
                + " after " + totalEnumerationSteps + " enum steps, "
                + totalMinimizations + " minimizations.");

        // Set final partition function values
        PartitionFunction.Values vals = getValues();
        if (computeMode == ComputeMode.FLAT_SUM) {
            vals.qstar = flatSumZLower;
            vals.pstar = flatSumZUpper;
        } else {
            if (corrLower.signum() == 0 && corrUpper.signum() == 0) {
                vals.qstar = logToPositiveBigDecimal(logDpLower);
                vals.pstar = logToPositiveBigDecimal(logDpUpper);
            } else {
                vals.qstar = logToPositiveBigDecimal(logDpLower).add(corrLower);
                vals.pstar = logToPositiveBigDecimal(logDpUpper).add(corrUpper);
            }
        }
        vals.qprime = vals.pstar.subtract(vals.qstar);

        if (epsilonBound <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
    }

    /** Convert log-space value to BigDecimal, handling overflow/underflow. */
    private static BigDecimal logToPositiveBigDecimal(double logVal) {
        if (logVal == Double.NEGATIVE_INFINITY) return BigDecimal.ZERO;
        if (logVal > 709) {
            double logBase = logVal / Math.log(10);
            int exponent = (int) Math.floor(logBase);
            double mantissa = Math.pow(10, logBase - exponent);
            return BigDecimal.valueOf(mantissa).scaleByPowerOfTen(exponent);
        }
        return BigDecimal.valueOf(Math.exp(logVal));
    }
}
