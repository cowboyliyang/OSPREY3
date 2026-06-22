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

import edu.duke.cs.osprey.astar.conf.RCs;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.*;

/**
 * Computes a branch decomposition of a residue interaction graph.
 * Uses the heuristic algorithm from Hicks (2005): star graph → push → 2-sep → 3-sep → eigenvector.
 *
 * Adapted from BWM's BranchDecomposition.BranchDecompositionH.
 * Key changes: int position indices instead of String vertex names,
 * InteractionGraph input instead of file I/O,
 * Apache Commons Math 3 instead of JAMA for eigenvalue decomposition.
 */
public class BranchDecomposition {

    public enum Strategy {
        HICKS,
        GREEDY_MERGE,
        EXACT,
        EXACT_IMPROVE,
        AUTO;

        public static Strategy fromProperty(String value) {
            if (value == null || value.trim().isEmpty()) {
                return HICKS;
            }
            switch (value.trim().toLowerCase(Locale.ROOT)) {
                case "hicks":
                case "hicks2005":
                case "legacy":
                    return HICKS;
                case "greedy":
                case "greedymerge":
                case "greedy_merge":
                case "edgegreedy":
                case "edge_greedy":
                    return GREEDY_MERGE;
                case "exact":
                case "exactbranchwidth":
                case "exact_branchwidth":
                case "optimal":
                    return EXACT;
                case "exactimprove":
                case "exact_improve":
                case "improve":
                case "improve_exact":
                case "targeted":
                case "exactfast":
                case "exact_fast":
                    return EXACT_IMPROVE;
                case "auto":
                case "best":
                    return AUTO;
                default:
                    System.err.println("BranchDecomposition: unknown strategy '" + value
                            + "', using hicks.");
                    return HICKS;
            }
        }
    }

    private BranchTree bt;
    private Set<Integer> graphVertices;
    private int numPositions;
    private InteractionGraph interactionGraph;
    private final Strategy strategy;
    private boolean computed = false;
    private static final int DEFAULT_EXACT_MAX_POSITIONS = 20;
    private static final int DEFAULT_EXACT_IMPROVE_MAX_DROP = 1;
    private static final long DEFAULT_EXACT_IMPROVE_MAX_MILLIS = 120_000L;

    /**
     * Build a branch decomposition from an interaction graph.
     * Call compute() to run the algorithm.
     */
    public BranchDecomposition(InteractionGraph graph) {
        this(graph, Strategy.HICKS);
    }

    public BranchDecomposition(InteractionGraph graph, Strategy strategy) {
        this.interactionGraph = graph;
        this.strategy = strategy == null ? Strategy.HICKS : strategy;
        this.numPositions = graph.getNumPositions();
        resetTree();
        if (this.strategy == Strategy.HICKS) {
            initializeHicksTree();
        }
    }

    private void resetTree() {
        this.bt = new BranchTree();
        this.graphVertices = new LinkedHashSet<>();
    }

    private void initializeHicksTree() {
        resetTree();
        // Build initial leaf nodes from interaction graph edges (analogous to readGraphFile)
        Set<Integer> gvDup = new LinkedHashSet<>();
        Set<Integer> seen = new LinkedHashSet<>();

        for (int[] edge : interactionGraph.getEdgeList()) {
            int pos1 = edge[0];
            int pos2 = edge[1];
            bt.addNode(new BranchNode(true, pos1, pos2));

            if (!seen.add(pos1)) gvDup.add(pos1);
            if (!seen.add(pos2)) gvDup.add(pos2);
        }
        graphVertices.addAll(seen);

        // Build star graph
        constructStarGraph(gvDup);
    }

    /** Run the branch decomposition algorithm. */
    public void compute() {
        if (computed) {
            return;
        }

        switch (strategy) {
            case HICKS:
                computeHicks();
                break;
            case GREEDY_MERGE:
                computeGreedyMerge();
                break;
            case EXACT:
                computeExact();
                break;
            case EXACT_IMPROVE:
                computeExactImprove();
                break;
            case AUTO:
                computeAuto();
                break;
            default:
                throw new IllegalStateException("Unhandled branch decomposition strategy " + strategy);
        }
        computed = true;
    }

    private void computeHicks() {
        boolean done = false;
        boolean done2sep = false;
        boolean done3sep = false;

        while (!done) {
            LinkedHashSet<BranchNode> bigNodes = getNodesToSplit();

            if (!bigNodes.isEmpty()) {
                LinkedHashSet<BranchEdge> sx = new LinkedHashSet<>();
                LinkedHashSet<BranchEdge> sy = new LinkedHashSet<>();

                BranchNode pn = findPushNode(bigNodes);

                if (pn != null) {
                    pushNode(pn, sx, sy);
                } else if (!done2sep) {
                    boolean found = false;
                    for (BranchNode bn : bigNodes) {
                        if (find2sep(bn, sx, sy)) {
                            pn = bn;
                            found = true;
                            break;
                        }
                    }
                    if (!found) done2sep = true;
                } else if (!done3sep) {
                    boolean found = false;
                    for (BranchNode bn : bigNodes) {
                        if (find3sep(bn, sx, sy)) {
                            pn = bn;
                            found = true;
                            break;
                        }
                    }
                    if (!found) done3sep = true;
                } else {
                    pn = bigNodes.iterator().next();
                    doEigen(pn, sx, sy);
                }

                if (!sx.isEmpty()) {
                    performSplit(pn, sx, sy);
                }
            } else {
                done = true;
            }
        }
    }

    public BranchTree getTree() { return bt; }
    public int getBranchwidth() { return bt.getBranchwidth(); }
    public Strategy getStrategy() { return strategy; }

    private void computeAuto() {
        BranchDecomposition hicks = new BranchDecomposition(interactionGraph, Strategy.HICKS);
        hicks.compute();

        BranchDecomposition greedy = new BranchDecomposition(interactionGraph, Strategy.GREEDY_MERGE);
        greedy.compute();

        BranchDecomposition selected = greedy.getBranchwidth() < hicks.getBranchwidth() ? greedy : hicks;
        this.bt = selected.bt;
        this.graphVertices = selected.graphVertices;
        System.out.println("BranchDecomposition auto: hicks bw=" + hicks.getBranchwidth()
                + ", greedyMerge bw=" + greedy.getBranchwidth()
                + ", selected=" + (selected == greedy ? "greedyMerge" : "hicks"));
    }

    private void computeExact() {
        int activePositions = countActiveGraphPositions();
        int maxPositions = getExactMaxPositions();
        if (activePositions > maxPositions) {
            System.out.println("BranchDecomposition exact: activePositions=" + activePositions
                    + " > maxPositions=" + maxPositions + ", falling back to hicks.");
            computeHicksFallback();
            return;
        }

        BranchDecomposition hicks = new BranchDecomposition(interactionGraph, Strategy.HICKS);
        hicks.compute();
        int hicksWidth = hicks.getBranchwidth();
        if (interactionGraph.getEdgeList().isEmpty()) {
            this.bt = hicks.bt;
            this.graphVertices = hicks.graphVertices;
            return;
        }

        int lowerBound = 0;
        try {
            lowerBound = ExactTreewidth.compute(interactionGraph).branchwidthLowerBound;
        } catch (RuntimeException e) {
            System.err.println("BranchDecomposition exact: exact treewidth lower bound failed: "
                    + e.getMessage() + ". Starting search at bw=0.");
        }

        ExactBranchwidth.Result exact = ExactBranchwidth.compute(
                interactionGraph, Math.min(lowerBound, hicksWidth), hicksWidth);
        if (exact.branchwidth <= hicksWidth) {
            this.bt = exact.tree;
            resetGraphVerticesFromEdges();
            System.out.println("BranchDecomposition exact: bw=" + exact.branchwidth
                    + ", hicksBw=" + hicksWidth
                    + ", lowerBound=" + exact.lowerBound
                    + ", activePositions=" + exact.activePositions
                    + ", graphEdges=" + exact.graphEdges
                    + ", elapsedMs=" + String.format(Locale.ROOT, "%.1f", exact.elapsedMillis()));
        } else {
            this.bt = hicks.bt;
            this.graphVertices = hicks.graphVertices;
            System.out.println("BranchDecomposition exact: reconstructed bw=" + exact.branchwidth
                    + " > hicksBw=" + hicksWidth + ", using hicks.");
        }
    }

    private void computeExactImprove() {
        int activePositions = countActiveGraphPositions();
        int maxPositions = getExactMaxPositions();
        BranchDecomposition hicks = new BranchDecomposition(interactionGraph, Strategy.HICKS);
        hicks.compute();
        int hicksWidth = hicks.getBranchwidth();

        if (activePositions > maxPositions) {
            this.bt = hicks.bt;
            this.graphVertices = hicks.graphVertices;
            System.out.println("BranchDecomposition exact_improve: activePositions=" + activePositions
                    + " > maxPositions=" + maxPositions + ", using hicks bw=" + hicksWidth + ".");
            return;
        }
        if (interactionGraph.getEdgeList().isEmpty() || hicksWidth <= 0) {
            this.bt = hicks.bt;
            this.graphVertices = hicks.graphVertices;
            return;
        }

        int lowerBound = 0;
        try {
            lowerBound = ExactTreewidth.compute(interactionGraph).branchwidthLowerBound;
        } catch (RuntimeException e) {
            System.err.println("BranchDecomposition exact_improve: exact treewidth lower bound failed: "
                    + e.getMessage() + ". Using lowerBound=0.");
        }

        int maxDrop = getExactImproveMaxDrop();
        int minTarget = Math.max(lowerBound, hicksWidth - maxDrop);
        long targetTimeoutNanos = getExactImproveMaxMillis() * 1_000_000L;

        for (int targetWidth = hicksWidth - 1; targetWidth >= minTarget; targetWidth--) {
            long started = System.nanoTime();
            try {
                ExactBranchwidth.Result exact = ExactBranchwidth.computeAtMost(
                        interactionGraph, targetWidth, targetTimeoutNanos);
                if (exact == null) {
                    System.out.println("BranchDecomposition exact_improve: no decomposition with bw<="
                            + targetWidth + " (hicksBw=" + hicksWidth
                            + ", elapsedMs=" + String.format(Locale.ROOT, "%.1f",
                            (System.nanoTime() - started) / 1_000_000.0) + "), using hicks.");
                    break;
                }

                this.bt = exact.tree;
                resetGraphVerticesFromEdges();
                System.out.println("BranchDecomposition exact_improve: bw=" + exact.branchwidth
                        + ", hicksBw=" + hicksWidth
                        + ", targetBw=" + targetWidth
                        + ", lowerBound=" + lowerBound
                        + ", activePositions=" + exact.activePositions
                        + ", graphEdges=" + exact.graphEdges
                        + ", elapsedMs=" + String.format(Locale.ROOT, "%.1f", exact.elapsedMillis()));
                return;
            } catch (ExactBranchwidth.SearchTimeoutException e) {
                System.out.println("BranchDecomposition exact_improve: targetBw=" + targetWidth
                        + " timed out after " + getExactImproveMaxMillis()
                        + " ms (hicksBw=" + hicksWidth + "), using hicks.");
                break;
            }
        }

        this.bt = hicks.bt;
        this.graphVertices = hicks.graphVertices;
    }

    private void computeHicksFallback() {
        BranchDecomposition hicks = new BranchDecomposition(interactionGraph, Strategy.HICKS);
        hicks.compute();
        this.bt = hicks.bt;
        this.graphVertices = hicks.graphVertices;
    }

    private int countActiveGraphPositions() {
        boolean[] active = new boolean[numPositions];
        int count = 0;
        for (int[] edge : interactionGraph.getEdgeList()) {
            if (!active[edge[0]]) {
                active[edge[0]] = true;
                count++;
            }
            if (!active[edge[1]]) {
                active[edge[1]] = true;
                count++;
            }
        }
        return count;
    }

    private void resetGraphVerticesFromEdges() {
        graphVertices = new LinkedHashSet<>();
        for (int[] edge : interactionGraph.getEdgeList()) {
            graphVertices.add(edge[0]);
            graphVertices.add(edge[1]);
        }
    }

    private static int getExactMaxPositions() {
        String value = System.getProperty("branchdp.decomp.exact.maxPositions");
        if (value == null || value.trim().isEmpty()) {
            value = System.getProperty("packstar.decomp.exact.maxPositions");
        }
        if (value == null || value.trim().isEmpty()) {
            return DEFAULT_EXACT_MAX_POSITIONS;
        }
        try {
            return Math.max(0, Integer.parseInt(value.trim()));
        } catch (NumberFormatException e) {
            System.err.println("BranchDecomposition exact: invalid maxPositions '" + value
                    + "', using " + DEFAULT_EXACT_MAX_POSITIONS + ".");
            return DEFAULT_EXACT_MAX_POSITIONS;
        }
    }

    private static int getExactImproveMaxDrop() {
        String value = System.getProperty("branchdp.decomp.exactImprove.maxDrop");
        if (value == null || value.trim().isEmpty()) {
            value = System.getProperty("packstar.decomp.exactImprove.maxDrop");
        }
        if (value == null || value.trim().isEmpty()) {
            return DEFAULT_EXACT_IMPROVE_MAX_DROP;
        }
        try {
            return Math.max(1, Integer.parseInt(value.trim()));
        } catch (NumberFormatException e) {
            System.err.println("BranchDecomposition exact_improve: invalid maxDrop '" + value
                    + "', using " + DEFAULT_EXACT_IMPROVE_MAX_DROP + ".");
            return DEFAULT_EXACT_IMPROVE_MAX_DROP;
        }
    }

    private static long getExactImproveMaxMillis() {
        String value = System.getProperty("branchdp.decomp.exactImprove.maxMillis");
        if (value == null || value.trim().isEmpty()) {
            value = System.getProperty("packstar.decomp.exactImprove.maxMillis");
        }
        if (value == null || value.trim().isEmpty()) {
            return DEFAULT_EXACT_IMPROVE_MAX_MILLIS;
        }
        try {
            return Math.max(1L, Long.parseLong(value.trim()));
        } catch (NumberFormatException e) {
            System.err.println("BranchDecomposition exact_improve: invalid maxMillis '" + value
                    + "', using " + DEFAULT_EXACT_IMPROVE_MAX_MILLIS + ".");
            return DEFAULT_EXACT_IMPROVE_MAX_MILLIS;
        }
    }

    // ========== Star graph construction ==========

    private void constructStarGraph(Set<Integer> gvDup) {
        BranchNode starNode = new BranchNode(false, -1, -1);
        bt.addNode(starNode);
        int starNodeIndex = starNode.getIndex();

        for (int i = 0; i < bt.getNumNodes(); i++) {
            if (i != starNodeIndex) {
                int eInd = bt.addEdge(i, starNodeIndex);
                BranchEdge be = bt.getEdge(eInd);
                BranchNode bn = bt.getNode(i);
                int pos1 = bn.getPos1();
                int pos2 = bn.getPos2();

                if (gvDup.contains(pos1)) be.addM(pos1);
                if (gvDup.contains(pos2)) be.addM(pos2);
            }
        }
    }

    // ========== Greedy edge-cluster branch decomposition ==========

    private static class GreedyCluster {
        final BranchNode rootNode;
        final BitSet edgeSet;
        final int leafCount;
        final int boundarySize;
        final long serial;

        GreedyCluster(BranchNode rootNode, BitSet edgeSet, int leafCount,
                      int boundarySize, long serial) {
            this.rootNode = rootNode;
            this.edgeSet = edgeSet;
            this.leafCount = leafCount;
            this.boundarySize = boundarySize;
            this.serial = serial;
        }
    }

    private static class GreedyPair {
        final int i;
        final int j;
        final int mergeWidth;
        final int unionBoundarySize;
        final int balance;
        final int unionLeafCount;
        final long serialSum;

        GreedyPair(int i, int j, int mergeWidth, int unionBoundarySize,
                   int balance, int unionLeafCount, long serialSum) {
            this.i = i;
            this.j = j;
            this.mergeWidth = mergeWidth;
            this.unionBoundarySize = unionBoundarySize;
            this.balance = balance;
            this.unionLeafCount = unionLeafCount;
            this.serialSum = serialSum;
        }

        boolean isBetterThan(GreedyPair other) {
            if (other == null) return true;
            if (mergeWidth != other.mergeWidth) return mergeWidth < other.mergeWidth;
            if (unionBoundarySize != other.unionBoundarySize) return unionBoundarySize < other.unionBoundarySize;
            if (balance != other.balance) return balance < other.balance;
            if (unionLeafCount != other.unionLeafCount) return unionLeafCount < other.unionLeafCount;
            return serialSum < other.serialSum;
        }
    }

    private void computeGreedyMerge() {
        resetTree();

        List<int[]> graphEdges = interactionGraph.getEdgeList();
        int numGraphEdges = graphEdges.size();
        if (numGraphEdges == 0) {
            return;
        }

        List<GreedyCluster> clusters = new ArrayList<>();
        long nextSerial = 0;
        for (int edgeIndex = 0; edgeIndex < numGraphEdges; edgeIndex++) {
            int[] edge = graphEdges.get(edgeIndex);
            BranchNode leaf = new BranchNode(true, edge[0], edge[1]);
            bt.addNode(leaf);
            graphVertices.add(edge[0]);
            graphVertices.add(edge[1]);

            BitSet edgeSet = new BitSet(numGraphEdges);
            edgeSet.set(edgeIndex);
            clusters.add(new GreedyCluster(leaf, edgeSet, 1,
                    middleSetSize(edgeSet, graphEdges), nextSerial++));
        }

        while (clusters.size() > 2) {
            GreedyPair pair = findGreedyPair(clusters, graphEdges, numGraphEdges);
            GreedyCluster merged = mergeGreedyClusters(
                    clusters.get(pair.i), clusters.get(pair.j), graphEdges, nextSerial++);

            clusters.remove(pair.j);
            clusters.remove(pair.i);
            clusters.add(merged);
        }

        if (clusters.size() == 1) {
            BranchNode root = new BranchNode(false, -1, -1);
            bt.addNode(root);
            int edgeIndex = bt.addEdge(root.getIndex(), clusters.get(0).rootNode.getIndex());
            if (edgeIndex >= 0) {
                bt.getEdge(edgeIndex).setM(middleSet(clusters.get(0).edgeSet, graphEdges));
            }
        } else {
            GreedyCluster left = clusters.get(0);
            GreedyCluster right = clusters.get(1);
            int edgeIndex = bt.addEdge(left.rootNode.getIndex(), right.rootNode.getIndex());
            if (edgeIndex >= 0) {
                bt.getEdge(edgeIndex).setM(middleSet(left.edgeSet, graphEdges));
            }
        }
    }

    private GreedyPair findGreedyPair(List<GreedyCluster> clusters, List<int[]> graphEdges,
                                      int numGraphEdges) {
        GreedyPair best = null;
        for (int i = 0; i < clusters.size(); i++) {
            GreedyCluster a = clusters.get(i);
            for (int j = i + 1; j < clusters.size(); j++) {
                GreedyCluster b = clusters.get(j);
                BitSet union = unionEdgeSets(a.edgeSet, b.edgeSet);
                int unionBoundarySize = middleSetSize(union, graphEdges);
                int mergeWidth = Math.max(Math.max(a.boundarySize, b.boundarySize), unionBoundarySize);
                int unionLeafCount = a.leafCount + b.leafCount;
                int balance = Math.abs(numGraphEdges - 2 * unionLeafCount);
                GreedyPair candidate = new GreedyPair(i, j, mergeWidth, unionBoundarySize,
                        balance, unionLeafCount, a.serial + b.serial);
                if (candidate.isBetterThan(best)) {
                    best = candidate;
                }
            }
        }
        return best;
    }

    private GreedyCluster mergeGreedyClusters(GreedyCluster left, GreedyCluster right,
                                              List<int[]> graphEdges, long serial) {
        BranchNode parent = new BranchNode(false, -1, -1);
        bt.addNode(parent);

        int leftEdgeIndex = bt.addEdge(parent.getIndex(), left.rootNode.getIndex());
        if (leftEdgeIndex >= 0) {
            bt.getEdge(leftEdgeIndex).setM(middleSet(left.edgeSet, graphEdges));
        }

        int rightEdgeIndex = bt.addEdge(parent.getIndex(), right.rootNode.getIndex());
        if (rightEdgeIndex >= 0) {
            bt.getEdge(rightEdgeIndex).setM(middleSet(right.edgeSet, graphEdges));
        }

        BitSet union = unionEdgeSets(left.edgeSet, right.edgeSet);
        return new GreedyCluster(parent, union, left.leafCount + right.leafCount,
                middleSetSize(union, graphEdges), serial);
    }

    private BitSet unionEdgeSets(BitSet a, BitSet b) {
        BitSet union = (BitSet) a.clone();
        union.or(b);
        return union;
    }

    private int middleSetSize(BitSet edgeSet, List<int[]> graphEdges) {
        return middleSet(edgeSet, graphEdges).size();
    }

    private LinkedHashSet<Integer> middleSet(BitSet edgeSet, List<int[]> graphEdges) {
        boolean[] inSubset = new boolean[numPositions];
        boolean[] outSubset = new boolean[numPositions];

        for (int edgeIndex = 0; edgeIndex < graphEdges.size(); edgeIndex++) {
            int[] edge = graphEdges.get(edgeIndex);
            boolean[] side = edgeSet.get(edgeIndex) ? inSubset : outSubset;
            side[edge[0]] = true;
            side[edge[1]] = true;
        }

        LinkedHashSet<Integer> middle = new LinkedHashSet<>();
        for (int pos = 0; pos < numPositions; pos++) {
            if (inSubset[pos] && outSubset[pos]) {
                middle.add(pos);
            }
        }
        return middle;
    }

    // ========== Node splitting logic ==========

    private LinkedHashSet<BranchNode> getNodesToSplit() {
        LinkedHashSet<BranchNode> pn = new LinkedHashSet<>();
        for (int i = 0; i < bt.getNumNodes(); i++) {
            BranchNode bn = bt.getNode(i);
            if (!bn.getIsLeaf() && bn.getNumEdges() > 3)
                pn.add(bn);
        }
        return pn;
    }

    private BranchNode findPushNode(LinkedHashSet<BranchNode> bigNodes) {
        for (BranchNode cn : bigNodes) {
            BranchEdge[] ce = bt.getEdgesForNode(cn.getIndex());
            for (int j = 0; j < cn.getNumEdges(); j++) {
                if (!ce[j].isChecked(cn))
                    return cn;
            }
        }
        return null;
    }

    // ========== Push operation ==========

    private void pushNode(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        BranchEdge[] pe = bt.getEdgesForNode(pn.getIndex());

        // Track which positions appear in which edge M-sets
        boolean[][] vInMe = new boolean[numPositions][pe.length];
        for (int k = 0; k < pe.length; k++) {
            for (int pos : pe[k].getM()) {
                if (pos >= 0 && pos < numPositions)
                    vInMe[pos][k] = true;
            }
        }

        for (int i = 0; i < pe.length; i++) {
            if (!pe[i].isChecked(pn)) {
                pe[i].setChecked(pn, true);

                for (int j = 0; j < pe.length; j++) {
                    if (i != j) {
                        LinkedHashSet<Integer> uMij = new LinkedHashSet<>(pe[i].getM());
                        uMij.addAll(pe[j].getM());

                        // Remove vertices only in M(i) ∪ M(j) but no other edges
                        LinkedHashSet<Integer> toRemove = new LinkedHashSet<>();
                        for (int pos : uMij) {
                            if (pos < 0 || pos >= numPositions) continue;
                            boolean found = false;
                            for (int l = 0; l < pe.length; l++) {
                                if (l != i && l != j && vInMe[pos][l]) {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) toRemove.add(pos);
                        }
                        uMij.removeAll(toRemove);

                        // Check pushing inequality
                        if (uMij.size() <= Math.max(pe[i].getM().size(), pe[j].getM().size())) {
                            sx.add(pe[i]);
                            sx.add(pe[j]);
                            for (int k = 0; k < pe.length; k++) {
                                if (k != i && k != j)
                                    sy.add(pe[k]);
                            }
                            return;
                        }
                    }
                }
            }
        }
    }

    // ========== 2-separation ==========

    private boolean find2sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        BranchEdge[] D = bt.getEdgesForNode(pn.getIndex());
        LinkedHashSet<Integer> S = new LinkedHashSet<>();
        for (BranchEdge e : D) S.addAll(e.getM());
        Object[] Sa = S.toArray();
        BranchTree H = constructGraphH(Sa, D);
        return find23sepHelper(H, sx, sy, D, true);
    }

    // ========== 3-separation ==========

    private boolean find3sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        BranchEdge[] D = bt.getEdgesForNode(pn.getIndex());
        LinkedHashSet<Integer> S = new LinkedHashSet<>();
        for (BranchEdge e : D) S.addAll(e.getM());
        Object[] Sa = S.toArray();
        BranchTree H = constructGraphH(Sa, D);

        // Check for 3-separation via triangle
        for (int i = 0; i < H.getNumNodes(); i++) {
            BranchNode x1 = H.getNode(i);
            if (x1.getNumEdges() == 3) {
                for (int j = 0; j < H.getNumNodes(); j++) {
                    if (j != i) {
                        BranchNode x2 = H.getNode(j);
                        if (H.edgeExists(x1.getIndex(), x2.getIndex())) {
                            for (int k = j + 1; k < H.getNumNodes(); k++) {
                                if (k != i) {
                                    BranchNode x3 = H.getNode(k);
                                    if (H.edgeExists(x1.getIndex(), x3.getIndex())
                                            && H.edgeExists(x2.getIndex(), x3.getIndex())) {
                                        // Found a triangle
                                        LinkedHashSet<BranchEdge> sxtmp = new LinkedHashSet<>();
                                        LinkedHashSet<BranchEdge> sytmp = new LinkedHashSet<>();
                                        for (BranchEdge d : D) {
                                            sytmp.add(d);
                                            if (d.getM().contains(x1.getPos1()))
                                                sxtmp.add(d);
                                        }
                                        sytmp.removeAll(sxtmp);
                                        if (sytmp.size() >= 2) {
                                            sx.addAll(sxtmp);
                                            sy.addAll(sytmp);
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Check for 3-separations with L and R sets of size >= 2
        return find23sepHelper(H, sx, sy, D, false);
    }

    // ========== 2/3-separation helper ==========

    private boolean find23sepHelper(BranchTree H, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy,
                                     BranchEdge[] D, boolean do2sep) {
        // Find node with minimum degree
        BranchNode nn = null;
        int minDegree = Integer.MAX_VALUE;
        for (int i = 0; i < H.getNumNodes(); i++) {
            if (H.getNode(i).getNumEdges() < minDegree) {
                nn = H.getNode(i);
                minDegree = H.getNode(i).getNumEdges();
            }
        }
        if (nn == null) return false;

        // Find neighbors of nn
        BranchEdge[] cDedges = H.getEdgesForNode(nn.getIndex());
        BranchNode[] nnNeighborsN = new BranchNode[cDedges.length];
        LinkedHashSet<Integer> nnNeighbors = new LinkedHashSet<>();
        for (int i = 0; i < cDedges.length; i++) {
            if (cDedges[i].getn1().getIndex() == nn.getIndex()) {
                nnNeighbors.add(cDedges[i].getn2().getIndex());
                nnNeighborsN[i] = cDedges[i].getn2();
            } else {
                nnNeighbors.add(cDedges[i].getn1().getIndex());
                nnNeighborsN[i] = cDedges[i].getn1();
            }
        }

        // Check all non-neighbors of nn
        for (int i = 0; i < H.getNumNodes(); i++) {
            if (!nnNeighbors.contains(H.getNode(i).getIndex()) && H.getNode(i).getIndex() != nn.getIndex()) {
                BranchEdge e = new BranchEdge(new BranchNode(false, -1, -1), new BranchNode(false, -1, -1));

                if (do2sep && find2sepHelper(H, nn, H.getNode(i), sx, sy, D, e))
                    return true;
                if (!do2sep && findKsepPairV(H, nn, H.getNode(i), 3, sx, sy, D, e))
                    return true;
            }
        }

        // Check pairs of neighbors of nn
        for (int i = 0; i < nnNeighborsN.length; i++) {
            for (int j = i + 1; j < nnNeighborsN.length; j++) {
                if (!H.edgeExists(nnNeighborsN[i].getIndex(), nnNeighborsN[j].getIndex())) {
                    BranchEdge e = new BranchEdge(new BranchNode(false, -1, -1), new BranchNode(false, -1, -1));

                    if (do2sep && find2sepHelper(H, nnNeighborsN[i], nnNeighborsN[j], sx, sy, D, e))
                        return true;
                    if (!do2sep && findKsepPairV(H, nnNeighborsN[i], nnNeighborsN[j], 3, sx, sy, D, e))
                        return true;
                }
            }
        }

        return false;
    }

    private boolean find2sepHelper(BranchTree H, BranchNode v, BranchNode w,
                                    LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy,
                                    BranchEdge[] D, BranchEdge e) {
        if (!H.edgeExists(v.getIndex(), w.getIndex())) {
            MinVertexCut mvc = new MinVertexCut();
            mvc.findMinVertexCut(H, v.getIndex(), w.getIndex());
            e.setM(mvc.getCutSet());

            if (mvc.getCutSet().size() == 2 && mvc.getLset().size() > 0 && mvc.getRset().size() > 0) {
                completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
                return true;
            }
        }
        return false;
    }

    // ========== k-separation with vertex contraction ==========

    private boolean findKsepPairV(BranchTree H, BranchNode v, BranchNode w, int k,
                                   LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy,
                                   BranchEdge[] D, BranchEdge e) {
        if (!H.edgeExists(v.getIndex(), w.getIndex())) {
            MinVertexCut mvc = new MinVertexCut();
            mvc.findMinVertexCut(H, v.getIndex(), w.getIndex());
            e.setM(mvc.getCutSet());

            if (mvc.getCutSet().size() == k) {
                // Check if both L and R have size >= 2, or M is not independent in original graph
                if ((mvc.getLset().size() >= 2 && mvc.getRset().size() >= 2)
                        || !isIndependentSetGraph(mvc.getCutSet())) {
                    completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
                    return true;
                }

                if (v.getNumEdges() <= k && w.getNumEdges() <= k) {
                    // Contract v and w to their neighbors
                    LinkedHashSet<BranchNode> vN = H.getNeighborsForNode(v, false);
                    LinkedHashSet<BranchNode> wN = H.getNeighborsForNode(w, false);
                    Object[] vNa = vN.toArray();
                    Object[] wNa = wN.toArray();

                    for (Object vNeighObj : vNa) {
                        BranchNode vNeigh = (BranchNode) vNeighObj;
                        if (!wN.contains(vNeigh)) {
                            for (Object wNeighObj : wNa) {
                                BranchNode wNeigh = (BranchNode) wNeighObj;
                                if (!vN.contains(wNeigh)) {
                                    if (!H.edgeExists(vNeigh.getIndex(), wNeigh.getIndex())) {
                                        int[] xyInd = new int[2];
                                        BranchTree modH = contractV(H, v, vNeigh, vNa, w, wNeigh, wNa, xyInd);
                                        mvc = new MinVertexCut();
                                        mvc.findMinVertexCut(modH, xyInd[0], xyInd[1]);
                                        if (mvc.getCutSet().size() <= k) {
                                            mvc.addToLset(v.getPos1());
                                            mvc.addToRset(w.getPos1());
                                            completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
                                            e.setM(mvc.getCutSet());
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (v.getNumEdges() <= k || w.getNumEdges() <= k) {
                    // One has degree <= k, the other > k
                    BranchNode sn = v.getNumEdges() <= k ? v : w;
                    BranchNode ln = v.getNumEdges() <= k ? w : v;
                    LinkedHashSet<BranchNode> sN = H.getNeighborsForNode(sn, false);
                    Object[] sNa = sN.toArray();

                    for (Object sNeighObj : sNa) {
                        BranchNode sNeigh = (BranchNode) sNeighObj;
                        if (!H.edgeExists(sNeigh.getIndex(), ln.getIndex())) {
                            int[] xyInd = new int[2];
                            BranchTree modH = contractV(H, sn, sNeigh, sNa, xyInd, ln);
                            mvc = new MinVertexCut();
                            mvc.findMinVertexCut(modH, xyInd[0], xyInd[1]);
                            if (mvc.getCutSet().size() <= k) {
                                mvc.addToLset(sn.getPos1());
                                completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
                                e.setM(mvc.getCutSet());
                                return true;
                            }
                        }
                    }
                } else {
                    // Both have degree > k
                    completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
                    e.setM(mvc.getCutSet());
                    return true;
                }
            }
        }
        return false;
    }

    /** Check if positions in set s are independent in the original interaction graph. */
    private boolean isIndependentSetGraph(LinkedHashSet<Integer> s) {
        Integer[] sA = s.toArray(new Integer[0]);
        for (int i = 0; i < sA.length; i++) {
            for (int j = i + 1; j < sA.length; j++) {
                if (sA[i] >= 0 && sA[j] >= 0 && interactionGraph.hasEdge(sA[i], sA[j]))
                    return false;
            }
        }
        return true;
    }

    // ========== Vertex contraction ==========

    /** Contract nodes v→x and w→y in graph H. */
    private BranchTree contractV(BranchTree h, BranchNode v, BranchNode x, Object[] vNa,
                                  BranchNode w, BranchNode y, Object[] wNa, int[] xyInd) {
        BranchTree modH = h.deepCopy();
        BranchNode vMod = modH.getNode(v.getIndex());
        BranchNode xMod = modH.getNode(x.getIndex());
        BranchNode wMod = modH.getNode(w.getIndex());
        BranchNode yMod = modH.getNode(y.getIndex());
        unifyEdges(modH, vNa, xMod);
        unifyEdges(modH, wNa, yMod);
        modH.deleteNode(vMod.getIndex());
        modH.deleteNode(wMod.getIndex());
        xyInd[0] = xMod.getIndex();
        xyInd[1] = yMod.getIndex();
        return modH;
    }

    /** Contract node v→x in graph H (single contraction). */
    private BranchTree contractV(BranchTree h, BranchNode v, BranchNode x, Object[] vNa,
                                  int[] xyInd, BranchNode ln) {
        BranchTree modH = h.deepCopy();
        BranchNode vMod = modH.getNode(v.getIndex());
        BranchNode xMod = modH.getNode(x.getIndex());
        BranchNode lnMod = modH.getNode(ln.getIndex());
        unifyEdges(modH, vNa, xMod);
        modH.deleteNode(vMod.getIndex());
        xyInd[0] = xMod.getIndex();
        xyInd[1] = lnMod.getIndex();
        return modH;
    }

    private void unifyEdges(BranchTree h, Object[] vNa, BranchNode x) {
        for (Object nodeObj : vNa) {
            BranchNode node = (BranchNode) nodeObj;
            if (x.getIndex() != node.getIndex())
                h.addEdge(x.getIndex(), node.getIndex());
        }
    }

    // ========== Eigenvector heuristic ==========

    private void doEigen(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        BranchEdge[] D = bt.getEdgesForNode(pn.getIndex());

        // Union of M-sets
        LinkedHashSet<Integer> S = new LinkedHashSet<>();
        for (BranchEdge e : D) S.addAll(e.getM());
        Object[] Sa = S.toArray();

        // Nv[i] = edges in D that contain vertex Sa[i] in their M-set
        BranchEdge[][] Nv = new BranchEdge[Sa.length][];
        for (int i = 0; i < Sa.length; i++) {
            List<BranchEdge> edges = new ArrayList<>();
            for (BranchEdge d : D) {
                if (d.getM().contains(Sa[i]))
                    edges.add(d);
            }
            Nv[i] = edges.toArray(new BranchEdge[0]);
        }

        // Build F matrix (Laplacian-like)
        double[][] F = new double[D.length][D.length];
        for (int i = 0; i < D.length; i++) {
            for (int j = 0; j < D.length; j++) {
                F[i][j] = 0.0;
                if (i != j) {
                    for (int k = 0; k < Nv.length; k++) {
                        if (isElement(Nv[k], D[i]) && isElement(Nv[k], D[j]))
                            F[i][j] += -1.0 / (Nv[k].length - 1);
                    }
                } else {
                    for (int k = 0; k < Nv.length; k++) {
                        if (isElement(Nv[k], D[i]))
                            F[i][j]++;
                    }
                }
            }
        }

        // Eigenvalue decomposition using Apache Commons Math 3
        RealMatrix Fm = new Array2DRowRealMatrix(F);
        EigenDecomposition ed = new EigenDecomposition(Fm);
        double[] eigenvalues = ed.getRealEigenvalues();

        // Find second smallest eigenvalue
        double[] sorted = Arrays.copyOf(eigenvalues, eigenvalues.length);
        Arrays.sort(sorted);
        double secondSmallest = sorted.length > 1 ? sorted[1] : sorted[0];
        int evlind = -1;
        for (int i = 0; i < eigenvalues.length; i++) {
            if (eigenvalues[i] == secondSmallest) {
                evlind = i;
                break;
            }
        }

        // Get eigenvector for second smallest eigenvalue
        double[] secondEvc = ed.getEigenvector(evlind).toArray();

        // Sort D by eigenvector component
        double[] secondEvcSorted = Arrays.copyOf(secondEvc, secondEvc.length);
        Arrays.sort(secondEvcSorted);
        BranchEdge[] Dsorted = new BranchEdge[D.length];
        boolean[] used = new boolean[D.length];
        for (int i = 0; i < secondEvcSorted.length; i++) {
            for (int j = 0; j < secondEvc.length; j++) {
                if (!used[j] && secondEvc[j] == secondEvcSorted[i]) {
                    Dsorted[i] = D[j];
                    used[j] = true;
                    break;
                }
            }
        }

        // Partition: first |D|/3 to sx, last |D|/3 to sy
        int numE = (int) Math.ceil(D.length / 3.0);
        for (int i = 0; i < numE; i++) {
            sx.add(Dsorted[i]);
            sy.add(Dsorted[Dsorted.length - 1 - i]);
        }

        // Build helper graph H and find min vertex cut
        BranchTree H = constructGraphH(Sa, D);
        int[] vi = new int[2];
        H = modifyGraphH(H, sx, sy, vi);

        MinVertexCut mvc = new MinVertexCut();
        mvc.findMinVertexCut(H.deepCopy(), vi[0], vi[1]);
        mvc.getCutSet().remove(-1);
        mvc.getLset().remove(-1);
        mvc.getRset().remove(-1);

        // Complete the partition
        completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
    }

    // ========== Helper graph construction ==========

    /**
     * Construct graph H with nodes = vertices in Sa (M-set positions),
     * edges between positions that co-occur in some M-set.
     */
    private BranchTree constructGraphH(Object[] Sa, BranchEdge[] D) {
        BranchTree H = new BranchTree();
        for (Object posObj : Sa) {
            int pos = (Integer) posObj;
            H.addNode(new BranchNode(true, pos, -1)); // vertex identity in pos1
        }
        for (int i = 0; i < H.getNumNodes(); i++) {
            int v1 = H.getNode(i).getPos1();
            for (int j = i + 1; j < H.getNumNodes(); j++) {
                int v2 = H.getNode(j).getPos1();
                for (BranchEdge d : D) {
                    if (d.getM().contains(v1) && d.getM().contains(v2)) {
                        H.addEdge(i, j);
                        break;
                    }
                }
            }
        }
        return H;
    }

    /**
     * Modify graph H by adding source and sink nodes connected to vertices in sx/sy M-sets.
     */
    private BranchTree modifyGraphH(BranchTree H, LinkedHashSet<BranchEdge> sx,
                                     LinkedHashSet<BranchEdge> sy, int[] vi) {
        BranchNode va = new BranchNode(true, -1, -1); // source sentinel
        H.addNode(va);
        addEdgesToArtificialNode(H, va, sx);

        BranchNode vb = new BranchNode(true, -1, -1); // sink sentinel
        H.addNode(vb);
        addEdgesToArtificialNode(H, vb, sy);

        vi[0] = va.getIndex();
        vi[1] = vb.getIndex();
        return H;
    }

    /** Add edges from artificial node vn to all positions in M-sets of the given edges. */
    private void addEdgesToArtificialNode(BranchTree H, BranchNode vn, LinkedHashSet<BranchEdge> edgeSet) {
        LinkedHashSet<Integer> positions = new LinkedHashSet<>();
        for (BranchEdge be : edgeSet) {
            positions.addAll(be.getM());
        }
        for (int pos : positions) {
            for (int j = 0; j < H.getNumNodes(); j++) {
                if (H.getNode(j).getPos1() == pos) {
                    H.addEdge(vn.getIndex(), H.getNode(j).getIndex());
                    break;
                }
            }
        }
    }

    // ========== Edge partition completion ==========

    /**
     * Complete the partition of edges D into sx and sy based on L-set and R-set.
     * For each unassigned edge in D:
     *   - if M intersects only L: add to sx
     *   - if M intersects only R: add to sy
     *   - if M intersects neither: add to sy (default)
     */
    private void completeXY(LinkedHashSet<Integer> Lset, LinkedHashSet<Integer> Rset,
                             BranchEdge[] D, LinkedHashSet<BranchEdge> sx,
                             LinkedHashSet<BranchEdge> sy, boolean inX) {
        for (BranchEdge d : D) {
            if (!sx.contains(d) && !sy.contains(d)) {
                LinkedHashSet<Integer> eM = d.getM();
                LinkedHashSet<Integer> intersectL = new LinkedHashSet<>(Lset);
                intersectL.retainAll(eM);
                LinkedHashSet<Integer> intersectR = new LinkedHashSet<>(Rset);
                intersectR.retainAll(eM);

                if (intersectR.isEmpty() && intersectL.isEmpty()) {
                    if (inX) sx.add(d);
                    else sy.add(d);
                } else if (intersectR.isEmpty()) {
                    sx.add(d);
                } else if (intersectL.isEmpty()) {
                    sy.add(d);
                } else {
                    // Edge connects both sides -- add to the larger set
                    if (sx.size() <= sy.size()) sx.add(d);
                    else sy.add(d);
                }
            }
        }
    }

    // ========== Split execution ==========

    private void performSplit(BranchNode bn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        BranchNode bx = new BranchNode(false, -1, -1);
        bt.addNode(bx);
        BranchNode by = new BranchNode(false, -1, -1);
        bt.addNode(by);

        BranchEdge[] bnEdges = bt.getEdgesForNode(bn.getIndex());

        for (BranchEdge bnEdge : bnEdges) {
            BranchNode otherNode;
            if (bnEdge.getn1().getIndex() == bn.getIndex())
                otherNode = bnEdge.getn2();
            else
                otherNode = bnEdge.getn1();

            int eInd;
            if (sx.contains(bnEdge))
                eInd = bt.addEdge(bx.getIndex(), otherNode.getIndex());
            else
                eInd = bt.addEdge(by.getIndex(), otherNode.getIndex());

            BranchEdge ce = bt.getEdge(eInd);
            ce.setM(bnEdge.getM());
            ce.setChecked(otherNode, bnEdge.isChecked(otherNode));
            if (sx.contains(bnEdge))
                ce.setChecked(bx, bnEdge.isChecked(bn));
            else
                ce.setChecked(by, bnEdge.isChecked(bn));
        }

        bt.deleteNode(bn.getIndex());

        int bxbyInd = bt.addEdge(bx.getIndex(), by.getIndex());
        BranchEdge bxby = bt.getEdge(bxbyInd);
        compMLRsets(bxby, sx, sy);
    }

    /** Compute M-set for new edge = intersection of union(M-sets in sx) and union(M-sets in sy). */
    private void compMLRsets(BranchEdge e, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy) {
        LinkedHashSet<Integer> uMsx = new LinkedHashSet<>();
        for (BranchEdge edge : sx) uMsx.addAll(edge.getM());

        LinkedHashSet<Integer> uMsy = new LinkedHashSet<>();
        for (BranchEdge edge : sy) uMsy.addAll(edge.getM());

        LinkedHashSet<Integer> xyM = new LinkedHashSet<>(uMsx);
        xyM.retainAll(uMsy);
        e.setM(xyM);
    }

    // ========== Utility ==========

    private boolean isElement(Object[] a, Object e) {
        for (Object obj : a) {
            if (obj == e) return true;
        }
        return false;
    }

    /** Get the index of a node in H that stores the given position. */
    private int getIndOfVertex(BranchTree H, int pos) {
        for (int i = 0; i < H.getNumNodes(); i++) {
            if (H.getNode(i).getPos1() == pos)
                return H.getNode(i).getIndex();
        }
        return -1;
    }

    // ========== Tree Rooting (BWM*-style) ==========

    /**
     * Convert the unrooted BranchTree into a rooted binary tree.
     * Following BWM*'s transformRootedTree():
     * 1. Create virtual root (name=-3) and s-node (name=-2)
     * 2. Root edge (root→s) has empty M-set, isRootEdge=true
     * 3. Use edge[0] as the split edge: its two endpoints become s's children
     * 4. Recursively assign children for internal nodes
     *
     * @return the virtual root RootedTreeNode
     */
    public RootedTreeNode rootBranchTree(RCs rcs) {
        return rootBranchTree(rcs, 0);
    }

    public RootedTreeNode rootBranchTree(RCs rcs, int splitEdgeIndex) {
        if (bt.getNumEdges() == 0) {
            return null;
        }
        if (bt.getNumEdges() == 1) {
            return rootSingleEdgeTree(rcs);
        }
        if (splitEdgeIndex < 0 || splitEdgeIndex >= bt.getNumEdges()) {
            throw new IllegalArgumentException("Invalid branch decomposition split edge index "
                    + splitEdgeIndex + " for " + bt.getNumEdges() + " edges");
        }

        // Map from BranchNode index → RootedTreeNode
        Map<Integer, RootedTreeNode> nodeMap = new HashMap<>();
        for (int i = 0; i < bt.getNumNodes(); i++) {
            BranchNode bn = bt.getNode(i);
            nodeMap.put(i, new RootedTreeNode(i, bn.getIsLeaf(), bn.getPos1(), bn.getPos2()));
        }

        // Create virtual root and s-node
        RootedTreeNode root = new RootedTreeNode(-3, false, -1, -1);
        RootedTreeNode s = new RootedTreeNode(-2, false, -1, -1);

        // Root edge: root → s, M = empty, isRootEdge = true
        RootedTreeEdge rootEdge = new RootedTreeEdge(root, s, new LinkedHashSet<>(), true, rcs);
        s.setChildOfEdge(rootEdge);
        s.setParent(root);
        root.setLeftChild(s);

        // Use the requested split edge. The legacy behavior is splitEdgeIndex=0.
        BranchEdge splitEdge = bt.getEdge(splitEdgeIndex);
        LinkedHashSet<Integer> splitM = splitEdge.getM();
        int n1Idx = splitEdge.getn1().getIndex();
        int n2Idx = splitEdge.getn2().getIndex();

        RootedTreeNode leftNode = nodeMap.get(n1Idx);
        RootedTreeNode rightNode = nodeMap.get(n2Idx);

        // s → leftNode (with split M-set)
        RootedTreeEdge leftEdge = new RootedTreeEdge(s, leftNode, splitM, false, rcs);
        leftNode.setChildOfEdge(leftEdge);
        leftNode.setParent(s);
        s.setLeftChild(leftNode);

        // s → rightNode (with split M-set)
        RootedTreeEdge rightEdge = new RootedTreeEdge(s, rightNode, splitM, false, rcs);
        rightNode.setChildOfEdge(rightEdge);
        rightNode.setParent(s);
        s.setRightChild(rightNode);

        // Track which edges have been used (the split edge is used)
        Set<Integer> usedEdges = new HashSet<>();
        usedEdges.add(splitEdgeIndex);

        // Recursively build subtrees
        if (!leftNode.getIsLeaf()) {
            buildRootedHelper(leftNode, nodeMap, usedEdges, rcs);
        }
        if (!rightNode.getIsLeaf()) {
            buildRootedHelper(rightNode, nodeMap, usedEdges, rcs);
        }

        return root;
    }

    private RootedTreeNode rootSingleEdgeTree(RCs rcs) {
        BranchEdge onlyEdge = bt.getEdge(0);
        BranchNode n1 = onlyEdge.getn1();
        BranchNode n2 = onlyEdge.getn2();
        BranchNode leaf = n1.getIsLeaf() ? n1 : (n2.getIsLeaf() ? n2 : null);
        if (leaf == null) {
            throw new RuntimeException("Single-edge branch tree has no leaf endpoint");
        }

        RootedTreeNode root = new RootedTreeNode(-3, false, -1, -1);
        RootedTreeNode leafNode = new RootedTreeNode(
                leaf.getIndex(), true, leaf.getPos1(), leaf.getPos2());
        RootedTreeEdge rootEdge = new RootedTreeEdge(
                root, leafNode, new LinkedHashSet<>(), true, rcs);

        leafNode.setChildOfEdge(rootEdge);
        leafNode.setParent(root);
        root.setLeftChild(leafNode);
        return root;
    }

    /**
     * For an internal node that already has its parent set,
     * find the 2 remaining incident edges in the unrooted tree
     * and assign them as left/right children.
     */
    private void buildRootedHelper(RootedTreeNode node, Map<Integer, RootedTreeNode> nodeMap,
                                    Set<Integer> usedEdges, RCs rcs) {
        int nodeIdx = node.getName();

        // Find the 2 unused edges incident to this node
        List<Integer> childEdgeIndices = new ArrayList<>();
        for (int i = 0; i < bt.getNumEdges(); i++) {
            if (usedEdges.contains(i)) continue;
            BranchEdge be = bt.getEdge(i);
            if (be.getn1().getIndex() == nodeIdx || be.getn2().getIndex() == nodeIdx) {
                childEdgeIndices.add(i);
            }
        }

        // Degenerate decompositions for tiny sparse graphs can produce unary
        // internal nodes. Keep those rooted instead of failing before DP starts.
        if (childEdgeIndices.size() > 2) {
            throw new RuntimeException("Internal node " + nodeIdx + " has " + childEdgeIndices.size()
                    + " remaining edges (expected <= 2)");
        }

        for (int side = 0; side < childEdgeIndices.size(); side++) {
            int edgeIdx = childEdgeIndices.get(side);
            usedEdges.add(edgeIdx);
            BranchEdge be = bt.getEdge(edgeIdx);

            // The child is the node on the other end
            int childIdx = (be.getn1().getIndex() == nodeIdx) ?
                    be.getn2().getIndex() : be.getn1().getIndex();
            RootedTreeNode childNode = nodeMap.get(childIdx);

            RootedTreeEdge edge = new RootedTreeEdge(node, childNode, be.getM(), false, rcs);
            childNode.setChildOfEdge(edge);
            childNode.setParent(node);

            if (side == 0) node.setLeftChild(childNode);
            else node.setRightChild(childNode);

            if (!childNode.getIsLeaf()) {
                buildRootedHelper(childNode, nodeMap, usedEdges, rcs);
            }
        }
    }

    // ========== TESS (Tree Estimation of State Space) ==========

    /**
     * Result of TESS computation, stored in log-space to avoid overflow.
     * logTess = log(TESS), logNaive = log(Naive).
     * ratio = TESS / Naive; in log-space: logTess - logNaive.
     * If logTess >= logNaive (ratio >= 1), decomposition is not beneficial.
     */
    public static class TESSResult {
        public final double logTess;
        public final double logNaive;

        public TESSResult(double logTess, double logNaive) {
            this.logTess = logTess;
            this.logNaive = logNaive;
        }

        /** TESS / Naive ratio. Values >= 1 mean decomposition is not beneficial. */
        public double getRatio() {
            return Math.exp(logTess - logNaive);
        }

        @Override
        public String toString() {
            return String.format("TESS(logTess=%.2f, logNaive=%.2f, ratio=%.4f)", logTess, logNaive, getRatio());
        }
    }

    /**
     * Compute TESS: the estimated total state space after branch decomposition.
     *
     * TESS = Σ (for each edge side with non-empty lambda) Π(numRCs[pos] for pos in sidePositions)
     * Naive = Π(numRCs[pos] for all positions)
     *
     * Uses log-space arithmetic: logTESS = log(Σ exp(logSideSpace_i))
     *
     * An edge side is only counted if its lambda (= sidePositions \ M) is non-empty,
     * matching BWM*'s compactTree logic which skips edges with empty lambda.
     */
    public TESSResult computeTESS(RCs rcs) {
        // Compute logNaive = Σ log(numRCs[pos]) for all positions
        double logNaive = 0.0;
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            int numRC = rcs.getNum(pos);
            if (numRC > 0) {
                logNaive += Math.log(numRC);
            }
        }

        // Accumulate TESS using log-sum-exp for numerical stability
        List<Double> logTerms = new ArrayList<>();

        for (int eIdx = 0; eIdx < bt.getNumEdges(); eIdx++) {
            BranchEdge edge = bt.getEdge(eIdx);
            Set<Integer> mSet = edge.getM();
            BranchNode n1 = edge.getn1();
            BranchNode n2 = edge.getn2();

            // Process each side of the edge
            for (int side = 0; side < 2; side++) {
                BranchNode startNode = (side == 0) ? n1 : n2;
                BranchNode excludeNode = (side == 0) ? n2 : n1;

                Set<Integer> sidePositions = collectPositionsForTESS(startNode, excludeNode);

                // Lambda = sidePositions \ M
                Set<Integer> lambda = new LinkedHashSet<>(sidePositions);
                lambda.removeAll(mSet);

                // Only count this side if lambda is non-empty (matches BWM* compactTree)
                if (!lambda.isEmpty()) {
                    double logSideSpace = 0.0;
                    for (int pos : sidePositions) {
                        int numRC = rcs.getNum(pos);
                        if (numRC > 0) {
                            logSideSpace += Math.log(numRC);
                        }
                    }
                    logTerms.add(logSideSpace);
                }
            }
        }

        // logTESS = log-sum-exp of all terms
        double logTess;
        if (logTerms.isEmpty()) {
            logTess = 0.0;
        } else {
            // log-sum-exp: log(Σ exp(x_i)) = max + log(Σ exp(x_i - max))
            double maxLogTerm = Collections.max(logTerms);
            double sumExp = 0.0;
            for (double logTerm : logTerms) {
                sumExp += Math.exp(logTerm - maxLogTerm);
            }
            logTess = maxLogTerm + Math.log(sumExp);
        }

        return new TESSResult(logTess, logNaive);
    }

    /** Collect all positions reachable from startNode without going through excludeNode via BFS. */
    private Set<Integer> collectPositionsForTESS(BranchNode startNode, BranchNode excludeNode) {
        Set<Integer> positions = new LinkedHashSet<>();
        Set<Integer> visited = new HashSet<>();
        Queue<BranchNode> bfsQueue = new LinkedList<>();

        visited.add(startNode.getIndex());
        visited.add(excludeNode.getIndex());
        bfsQueue.add(startNode);

        while (!bfsQueue.isEmpty()) {
            BranchNode current = bfsQueue.poll();

            if (current.getIsLeaf()) {
                if (current.getPos1() >= 0) positions.add(current.getPos1());
                if (current.getPos2() >= 0) positions.add(current.getPos2());
            }

            BranchEdge[] edges = bt.getEdgesForNode(current.getIndex());
            for (BranchEdge e : edges) {
                BranchNode neighbor = (e.getn1().getIndex() == current.getIndex()) ? e.getn2() : e.getn1();
                if (visited.add(neighbor.getIndex())) {
                    bfsQueue.add(neighbor);
                }
            }
        }

        return positions;
    }
}
