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
** You should have received a copy of OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.branchdp;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;

/**
 * Exact branchwidth for small simple interaction graphs.
 *
 * <p>The decision DP follows the vertex-subset recurrence of Kaneda,
 * Kobayashi, and Tamaki (2026). It is exponential in the number of active
 * graph vertices, so callers should keep this behind a small-position guard.</p>
 */
final class ExactBranchwidth {

    static final class SearchTimeoutException extends RuntimeException {
        SearchTimeoutException(int targetWidth) {
            super("Exact branchwidth target bw<=" + targetWidth + " timed out");
        }
    }

    static final class Result {
        final int branchwidth;
        final BranchTree tree;
        final int lowerBound;
        final int upperBound;
        final int activePositions;
        final int graphEdges;
        final long elapsedNanos;

        private Result(int branchwidth, BranchTree tree, int lowerBound, int upperBound,
                       int activePositions, int graphEdges, long elapsedNanos) {
            this.branchwidth = branchwidth;
            this.tree = tree;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
            this.activePositions = activePositions;
            this.graphEdges = graphEdges;
            this.elapsedNanos = elapsedNanos;
        }

        double elapsedMillis() {
            return elapsedNanos / 1_000_000.0;
        }
    }

    private static final class Shape {
        final int edgeIndex;
        final Shape left;
        final Shape right;
        final BitSet edges;

        private Shape(int edgeIndex, Shape left, Shape right, BitSet edges) {
            this.edgeIndex = edgeIndex;
            this.left = left;
            this.right = right;
            this.edges = edges;
        }

        static Shape leaf(int edgeIndex) {
            BitSet edges = new BitSet();
            edges.set(edgeIndex);
            return new Shape(edgeIndex, null, null, edges);
        }

        static Shape internal(Shape left, Shape right) {
            BitSet edges = new BitSet();
            edges.or(left.edges);
            edges.or(right.edges);
            return new Shape(-1, left, right, edges);
        }

        boolean isLeaf() {
            return edgeIndex >= 0;
        }
    }

    private static final class SplitChoice {
        final int x1;
        final int x2;
        final long score;

        SplitChoice(int x1, int x2, long score) {
            this.x1 = x1;
            this.x2 = x2;
            this.score = score;
        }
    }

    private final InteractionGraph graph;
    private final List<int[]> originalEdges;
    private final int originalNumPositions;
    private final int[] compactToOriginal;
    private final int[] compactU;
    private final int[] compactV;
    private final int[] compactEdgeMask;
    private final int[] adjacency;
    private final int activePositions;
    private final int graphEdges;
    private final int fullMask;
    private final int[] boundaryMask;
    private final int[] edgeCount;

    private ExactBranchwidth(InteractionGraph graph) {
        this.graph = graph;
        this.originalEdges = graph.getEdgeList();
        this.originalNumPositions = graph.getNumPositions();
        this.graphEdges = originalEdges.size();

        boolean[] active = new boolean[originalNumPositions];
        for (int[] edge : originalEdges) {
            active[edge[0]] = true;
            active[edge[1]] = true;
        }

        int[] originalToCompact = new int[originalNumPositions];
        java.util.Arrays.fill(originalToCompact, -1);
        List<Integer> compact = new ArrayList<>();
        for (int pos = 0; pos < originalNumPositions; pos++) {
            if (active[pos]) {
                originalToCompact[pos] = compact.size();
                compact.add(pos);
            }
        }

        this.activePositions = compact.size();
        this.compactToOriginal = new int[activePositions];
        for (int i = 0; i < activePositions; i++) {
            compactToOriginal[i] = compact.get(i);
        }

        this.compactU = new int[graphEdges];
        this.compactV = new int[graphEdges];
        this.compactEdgeMask = new int[graphEdges];
        this.adjacency = new int[activePositions];
        for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
            int[] edge = originalEdges.get(edgeIndex);
            int u = originalToCompact[edge[0]];
            int v = originalToCompact[edge[1]];
            compactU[edgeIndex] = u;
            compactV[edgeIndex] = v;
            compactEdgeMask[edgeIndex] = (1 << u) | (1 << v);
            adjacency[u] |= 1 << v;
            adjacency[v] |= 1 << u;
        }

        this.fullMask = activePositions == 0 ? 0 : (1 << activePositions) - 1;
        int numMasks = 1 << activePositions;
        this.boundaryMask = new int[numMasks];
        this.edgeCount = new int[numMasks];
        precomputeMasks();
    }

    static Result compute(InteractionGraph graph, int lowerBound, int upperBound) {
        long started = System.nanoTime();
        ExactBranchwidth exact = new ExactBranchwidth(graph);

        if (exact.graphEdges == 0) {
            return new Result(0, new BranchTree(), 0, 0,
                    exact.activePositions, exact.graphEdges, System.nanoTime() - started);
        }

        int lower = Math.max(0, lowerBound);
        int upper = Math.max(lower, upperBound);
        BranchTree bestTree = null;
        int bestWidth = Integer.MAX_VALUE;

        for (int k = lower; k <= upper; k++) {
            FixedKSearch search = exact.new FixedKSearch(k);
            if (!search.feasible(exact.fullMask)) {
                continue;
            }

            Shape shape = exact.buildShape(exact.fullMask, search);
            BranchTree tree = exact.toBranchTree(shape);
            exact.recomputeMiddleSets(tree);
            int width = tree.getBranchwidth();
            if (width < bestWidth) {
                bestWidth = width;
                bestTree = tree;
            }
            if (width <= k) {
                return new Result(width, tree, lower, upper,
                        exact.activePositions, exact.graphEdges, System.nanoTime() - started);
            }
        }

        if (bestTree == null) {
            throw new IllegalStateException("Exact branchwidth DP found no decomposition up to width " + upper);
        }
        return new Result(bestWidth, bestTree, lower, upper,
                exact.activePositions, exact.graphEdges, System.nanoTime() - started);
    }

    static Result computeAtMost(InteractionGraph graph, int targetWidth, long timeoutNanos) {
        long started = System.nanoTime();
        ExactBranchwidth exact = new ExactBranchwidth(graph);

        if (exact.graphEdges == 0) {
            return new Result(0, new BranchTree(), 0, targetWidth,
                    exact.activePositions, exact.graphEdges, System.nanoTime() - started);
        }

        int target = Math.max(0, targetWidth);
        long deadlineNanos = timeoutNanos > 0 ? started + timeoutNanos : 0L;
        FixedKSearch search = exact.new FixedKSearch(target, deadlineNanos);
        if (!search.feasible(exact.fullMask)) {
            return null;
        }

        Shape shape = exact.buildShape(exact.fullMask, search);
        BranchTree tree = exact.toBranchTree(shape);
        exact.recomputeMiddleSets(tree);
        int width = tree.getBranchwidth();
        if (width > target) {
            return null;
        }

        return new Result(width, tree, target, target,
                exact.activePositions, exact.graphEdges, System.nanoTime() - started);
    }

    private void precomputeMasks() {
        for (int mask = 0; mask <= fullMask; mask++) {
            int outside = fullMask & ~mask;
            int boundary = 0;
            for (int v = 0; v < activePositions; v++) {
                int bit = 1 << v;
                if ((mask & bit) != 0 && (adjacency[v] & outside) != 0) {
                    boundary |= bit;
                }
            }
            boundaryMask[mask] = boundary;

            int count = 0;
            for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
                int edgeMask = compactEdgeMask[edgeIndex];
                if ((edgeMask & ~mask) == 0) {
                    count++;
                }
            }
            edgeCount[mask] = count;
        }
    }

    private final class FixedKSearch {
        private final int k;
        private final long deadlineNanos;
        private long steps = 0L;
        private final byte[] state;
        private final int[] split1;
        private final int[] split2;

        FixedKSearch(int k) {
            this(k, 0L);
        }

        FixedKSearch(int k, long deadlineNanos) {
            this.k = k;
            this.deadlineNanos = deadlineNanos;
            int numMasks = 1 << activePositions;
            this.state = new byte[numMasks];
            this.split1 = new int[numMasks];
            this.split2 = new int[numMasks];
            java.util.Arrays.fill(split1, -1);
            java.util.Arrays.fill(split2, -1);
        }

        boolean feasible(int mask) {
            checkTimeout();
            byte cached = state[mask];
            if (cached != 0) {
                return cached == 2;
            }
            if (edgeCount[mask] <= 1) {
                state[mask] = 2;
                return true;
            }
            if (Integer.bitCount(boundaryMask[mask]) > k) {
                state[mask] = 1;
                return false;
            }

            SplitChoice choice = findSplit(mask);
            if (choice == null) {
                state[mask] = 1;
                return false;
            }
            split1[mask] = choice.x1;
            split2[mask] = choice.x2;
            state[mask] = 2;
            return true;
        }

        private SplitChoice findSplit(int mask) {
            SplitChoice best = null;
            for (int x1 = (mask - 1) & mask; x1 != 0; x1 = (x1 - 1) & mask) {
                int missingFromX1 = mask & ~x1;
                if (missingFromX1 == 0 || Integer.bitCount(boundaryMask[x1]) > k) {
                    continue;
                }

                for (int overlap = x1; ; overlap = (overlap - 1) & x1) {
                    checkTimeout();
                    int x2 = missingFromX1 | overlap;
                    if (x2 != mask && x1 <= x2 && Integer.bitCount(boundaryMask[x2]) <= k) {
                        SplitChoice candidate = evaluateSplit(mask, x1, x2);
                        if (candidate != null && feasible(x1) && feasible(x2)
                                && (best == null || candidate.score < best.score)) {
                            best = candidate;
                        }
                    }
                    if (overlap == 0) {
                        break;
                    }
                }
            }
            return best;
        }

        private void checkTimeout() {
            if (deadlineNanos == 0L) {
                return;
            }
            if ((steps++ & 8191L) == 0L && System.nanoTime() > deadlineNanos) {
                throw new SearchTimeoutException(k);
            }
        }

        private SplitChoice evaluateSplit(int mask, int x1, int x2) {
            int intersection = x1 & x2;
            if (intersection != (boundaryMask[x1] & boundaryMask[x2])) {
                return null;
            }

            int fPlusVertices = fPlusVertexMask(mask, x1, x2);
            int separator = boundaryMask[mask] | fPlusVertices;
            int separatorSize = Integer.bitCount(separator);
            if (separatorSize > k) {
                return null;
            }

            int overlapEdges = commonEdgeCount(x1, x2);
            int maxChildEdges = Math.max(edgeCount[x1], edgeCount[x2]);
            int balance = Math.abs(edgeCount[x1] - edgeCount[x2]);

            long score = (long) separatorSize * 1_000_000_000L
                    + (long) overlapEdges * 10_000_000L
                    + (long) maxChildEdges * 10_000L
                    + balance;
            return new SplitChoice(x1, x2, score);
        }
    }

    private int fPlusVertexMask(int mask, int x1, int x2) {
        int fPlus = 0;
        for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
            int edgeMask = compactEdgeMask[edgeIndex];
            if ((edgeMask & ~mask) == 0
                    && (edgeMask & ~x1) != 0
                    && (edgeMask & ~x2) != 0) {
                fPlus |= edgeMask;
            }
        }
        return fPlus;
    }

    private int commonEdgeCount(int x1, int x2) {
        int count = 0;
        for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
            int edgeMask = compactEdgeMask[edgeIndex];
            if ((edgeMask & ~x1) == 0 && (edgeMask & ~x2) == 0) {
                count++;
            }
        }
        return count;
    }

    private Shape buildShape(int mask, FixedKSearch search) {
        int edges = edgeCount[mask];
        if (edges == 0) {
            return null;
        }
        if (edges == 1) {
            return Shape.leaf(findSingleEdge(mask));
        }

        int x1 = search.split1[mask];
        int x2 = search.split2[mask];
        if (x1 < 0 || x2 < 0) {
            throw new IllegalStateException("No exact branchwidth split recorded for mask " + mask);
        }

        Shape left = buildShape(x1, search);
        Shape right = buildShape(x2, search);
        if (left != null && right != null) {
            right = pruneEdges(right, left.edges);
        }

        Shape shape = compose(left, right);
        BitSet present = shape == null ? new BitSet(graphEdges) : (BitSet) shape.edges.clone();
        BitSet missing = edgesInMask(mask);
        missing.andNot(present);
        shape = compose(shape, balancedLeaves(missing));
        return shape;
    }

    private int findSingleEdge(int mask) {
        for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
            int edgeMask = compactEdgeMask[edgeIndex];
            if ((edgeMask & ~mask) == 0) {
                return edgeIndex;
            }
        }
        throw new IllegalStateException("Mask has no edge");
    }

    private BitSet edgesInMask(int mask) {
        BitSet edges = new BitSet(graphEdges);
        for (int edgeIndex = 0; edgeIndex < graphEdges; edgeIndex++) {
            int edgeMask = compactEdgeMask[edgeIndex];
            if ((edgeMask & ~mask) == 0) {
                edges.set(edgeIndex);
            }
        }
        return edges;
    }

    private Shape pruneEdges(Shape shape, BitSet remove) {
        if (shape == null) {
            return null;
        }
        if (shape.isLeaf()) {
            return remove.get(shape.edgeIndex) ? null : shape;
        }
        return compose(pruneEdges(shape.left, remove), pruneEdges(shape.right, remove));
    }

    private Shape compose(Shape left, Shape right) {
        if (left == null) return right;
        if (right == null) return left;
        return Shape.internal(left, right);
    }

    private Shape balancedLeaves(BitSet edges) {
        List<Shape> leaves = new ArrayList<>();
        for (int edgeIndex = edges.nextSetBit(0); edgeIndex >= 0; edgeIndex = edges.nextSetBit(edgeIndex + 1)) {
            leaves.add(Shape.leaf(edgeIndex));
        }
        return combineBalanced(leaves);
    }

    private Shape combineBalanced(List<Shape> shapes) {
        if (shapes.isEmpty()) {
            return null;
        }
        List<Shape> work = new ArrayList<>(shapes);
        work.sort(Comparator.comparingInt(a -> a.edges.cardinality()));
        while (work.size() > 1) {
            Shape a = work.remove(0);
            Shape b = work.remove(0);
            Shape merged = Shape.internal(a, b);
            int insert = Collections.binarySearch(work, merged,
                    Comparator.comparingInt(s -> s.edges.cardinality()));
            if (insert < 0) {
                insert = -insert - 1;
            }
            work.add(insert, merged);
        }
        return work.get(0);
    }

    private BranchTree toBranchTree(Shape shape) {
        BranchTree tree = new BranchTree();
        if (shape == null) {
            return tree;
        }

        int rootIndex = addShapeNode(tree, shape);
        if (shape.isLeaf()) {
            BranchNode root = new BranchNode(false, -1, -1);
            tree.addNode(root);
            tree.addEdge(root.getIndex(), rootIndex);
        }
        return tree;
    }

    private int addShapeNode(BranchTree tree, Shape shape) {
        int nodeIndex;
        if (shape.isLeaf()) {
            int[] edge = originalEdges.get(shape.edgeIndex);
            BranchNode leaf = new BranchNode(true, edge[0], edge[1]);
            tree.addNode(leaf);
            nodeIndex = leaf.getIndex();
        } else {
            BranchNode internal = new BranchNode(false, -1, -1);
            tree.addNode(internal);
            nodeIndex = internal.getIndex();
            int leftIndex = addShapeNode(tree, shape.left);
            int rightIndex = addShapeNode(tree, shape.right);
            tree.addEdge(nodeIndex, leftIndex);
            tree.addEdge(nodeIndex, rightIndex);
        }
        return nodeIndex;
    }

    private void recomputeMiddleSets(BranchTree tree) {
        Map<Long, Integer> edgeIndexByKey = new HashMap<>();
        for (int edgeIndex = 0; edgeIndex < originalEdges.size(); edgeIndex++) {
            int[] edge = originalEdges.get(edgeIndex);
            edgeIndexByKey.put(edgeKey(edge[0], edge[1]), edgeIndex);
        }

        for (int branchEdgeIndex = 0; branchEdgeIndex < tree.getNumEdges(); branchEdgeIndex++) {
            BranchEdge edge = tree.getEdge(branchEdgeIndex);
            BitSet sideEdges = collectLeafEdgesOnSide(tree,
                    edge.getn1().getIndex(), edge.getn2().getIndex(), edgeIndexByKey);
            edge.setM(middleSet(sideEdges));
        }
    }

    private BitSet collectLeafEdgesOnSide(BranchTree tree, int startNode, int blockedNode,
                                          Map<Long, Integer> edgeIndexByKey) {
        BitSet sideEdges = new BitSet(graphEdges);
        boolean[] visited = new boolean[tree.getNumNodes()];
        Queue<Integer> queue = new ArrayDeque<>();
        visited[blockedNode] = true;
        visited[startNode] = true;
        queue.add(startNode);

        while (!queue.isEmpty()) {
            int nodeIndex = queue.poll();
            BranchNode node = tree.getNode(nodeIndex);
            if (node.getIsLeaf()) {
                Integer edgeIndex = edgeIndexByKey.get(edgeKey(node.getPos1(), node.getPos2()));
                if (edgeIndex == null) {
                    throw new IllegalStateException("Exact branch tree leaf is not an interaction edge: "
                            + node.getPos1() + "," + node.getPos2());
                }
                sideEdges.set(edgeIndex);
            }

            for (BranchEdge incident : tree.getEdgesForNode(nodeIndex)) {
                int neighbor = incident.getn1().getIndex() == nodeIndex
                        ? incident.getn2().getIndex()
                        : incident.getn1().getIndex();
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.add(neighbor);
                }
            }
        }
        return sideEdges;
    }

    private LinkedHashSet<Integer> middleSet(BitSet edgeSet) {
        boolean[] inSubset = new boolean[originalNumPositions];
        boolean[] outSubset = new boolean[originalNumPositions];

        for (int edgeIndex = 0; edgeIndex < originalEdges.size(); edgeIndex++) {
            int[] edge = originalEdges.get(edgeIndex);
            boolean[] side = edgeSet.get(edgeIndex) ? inSubset : outSubset;
            side[edge[0]] = true;
            side[edge[1]] = true;
        }

        LinkedHashSet<Integer> middle = new LinkedHashSet<>();
        for (int pos = 0; pos < originalNumPositions; pos++) {
            if (inSubset[pos] && outSubset[pos]) {
                middle.add(pos);
            }
        }
        return middle;
    }

    private long edgeKey(int a, int b) {
        int lo = Math.min(a, b);
        int hi = Math.max(a, b);
        return (((long) lo) << 32) ^ (hi & 0xffffffffL);
    }
}
