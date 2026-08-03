package edu.duke.cs.osprey.branchdp;

import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestBranchDecompositionStrategies {

    @Test
    public void greedyMergeBuildsRootableBranchTree() {
        InteractionGraph graph = InteractionGraph.buildComplete(5);
        BranchDecomposition bd = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.GREEDY_MERGE);

        bd.compute();

        assertEquals(graph.getNumEdges(), countLeaves(bd.getTree()));
        assertTrue(bd.getBranchwidth() <= graph.getNumPositions());
        assertNotNull(bd.rootBranchTree(null));
    }

    @Test
    public void autoNeverChoosesHigherBranchwidthThanHicks() {
        InteractionGraph graph = InteractionGraph.buildComplete(6);

        BranchDecomposition hicks = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.HICKS);
        hicks.compute();

        BranchDecomposition auto = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.AUTO);
        auto.compute();

        assertTrue(auto.getBranchwidth() <= hicks.getBranchwidth());
        assertEquals(graph.getNumEdges(), countLeaves(auto.getTree()));
    }

    @Test
    public void weightedHicksBuildsRootableTreeWithoutIncreasingWeightedMaximum() {
        InteractionGraph graph = InteractionGraph.buildComplete(6);
        int[] rcCounts = {97, 2, 3, 5, 7, 11};

        BranchDecomposition hicks = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.HICKS, rcCounts);
        hicks.compute();
        BranchDecomposition weighted = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.WEIGHTED_HICKS, rcCounts);
        weighted.compute();

        assertEquals(graph.getNumEdges(), countLeaves(weighted.getTree()));
        assertNotNull(weighted.rootBranchTree(null));
        assertTrue(maxBoundaryStates(weighted.getTree(), rcCounts)
                        .compareTo(maxBoundaryStates(hicks.getTree(), rcCounts)) <= 0,
                "weighted Hicks/NNI must not increase the maximum separator state product");
    }

    @Test
    public void adaptiveStandaloneUsesRootableWeightedFirstStage() {
        InteractionGraph graph = InteractionGraph.buildComplete(6);
        int[] rcCounts = {97, 2, 3, 5, 7, 11};

        BranchDecomposition weighted = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.WEIGHTED_HICKS, rcCounts);
        weighted.compute();
        BranchDecomposition adaptive = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.ADAPTIVE, rcCounts);
        adaptive.compute();

        assertEquals(BranchDecomposition.Strategy.ADAPTIVE,
                BranchDecomposition.Strategy.fromProperty("adaptive_exact"));
        assertEquals(graph.getNumEdges(), countLeaves(adaptive.getTree()));
        assertNotNull(adaptive.rootBranchTree(null));
        assertEquals(maxBoundaryStates(weighted.getTree(), rcCounts),
                maxBoundaryStates(adaptive.getTree(), rcCounts));
    }

    @Test
    public void exactTreewidthMatchesSmallGraphs() {
        ExactTreewidth.Result path = ExactTreewidth.compute(InteractionGraph.buildFromEdges(5, Arrays.asList(
                new int[] { 0, 1 },
                new int[] { 1, 2 },
                new int[] { 2, 3 },
                new int[] { 3, 4 })));
        assertEquals(1, path.treewidth);
        assertEquals(1, path.branchwidthLowerBound);

        ExactTreewidth.Result cycle = ExactTreewidth.compute(InteractionGraph.buildFromEdges(4, Arrays.asList(
                new int[] { 0, 1 },
                new int[] { 1, 2 },
                new int[] { 2, 3 },
                new int[] { 3, 0 })));
        assertEquals(2, cycle.treewidth);
        assertEquals(2, cycle.branchwidthLowerBound);

        ExactTreewidth.Result clique = ExactTreewidth.compute(InteractionGraph.buildComplete(6));
        assertEquals(5, clique.treewidth);
        assertEquals(4, clique.branchwidthLowerBound);
    }

    @Test
    public void exactBranchwidthMatchesSmallGraphs() {
        assertExactBranchwidth(InteractionGraph.buildFromEdges(4, Arrays.asList(
                new int[] { 0, 1 },
                new int[] { 0, 2 },
                new int[] { 0, 3 })), 1);

        assertExactBranchwidth(InteractionGraph.buildFromEdges(4, Arrays.asList(
                new int[] { 0, 1 },
                new int[] { 1, 2 },
                new int[] { 2, 3 },
                new int[] { 3, 0 })), 2);

        assertExactBranchwidth(InteractionGraph.buildComplete(4), 3);
        assertExactBranchwidth(InteractionGraph.buildComplete(6), 4);
    }

    private static int countLeaves(BranchTree tree) {
        int count = 0;
        for (int i = 0; i < tree.getNumNodes(); i++) {
            if (tree.getNode(i).getIsLeaf()) {
                count++;
            }
        }
        return count;
    }

    private static BigInteger maxBoundaryStates(BranchTree tree, int[] rcCounts) {
        BigInteger max = BigInteger.ZERO;
        for (int edgeIndex = 0; edgeIndex < tree.getNumEdges(); edgeIndex++) {
            BigInteger states = BigInteger.ONE;
            for (int pos : tree.getEdge(edgeIndex).getM()) {
                states = states.multiply(BigInteger.valueOf(rcCounts[pos]));
            }
            max = max.max(states);
        }
        return max;
    }

    private static void assertExactBranchwidth(InteractionGraph graph, int expectedBranchwidth) {
        BranchDecomposition exact = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.EXACT);
        exact.compute();

        assertEquals(expectedBranchwidth, exact.getBranchwidth());
        assertEquals(graph.getNumEdges(), countLeaves(exact.getTree()));
        assertNotNull(exact.rootBranchTree(null));
    }
}
