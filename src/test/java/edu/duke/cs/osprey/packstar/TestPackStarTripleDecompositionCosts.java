package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.BranchDpBackend;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.branchdp.RootedTreeEdge;
import edu.duke.cs.osprey.branchdp.RootedTreeNode;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class TestPackStarTripleDecompositionCosts {
    private static long edge(int a, int b) {
        return ((long) Math.min(a, b) << 32) | Math.max(a, b);
    }

    private static PackStarTripleDecompositionCosts.Cost cost(long work, long table) {
        return new PackStarTripleDecompositionCosts.Cost(2, table / 16, table / 16,
                table, table, table, 0, BigInteger.valueOf(work));
    }

    private static PackStarTripleDecompositionCosts.Limits unlimited() {
        return new PackStarTripleDecompositionCosts.Limits(Long.MAX_VALUE, Long.MAX_VALUE,
                Long.MAX_VALUE, Long.MAX_VALUE, Long.MAX_VALUE);
    }

    private static final class Fixture {
        final RCs rcs;
        final InteractionGraph graph;
        final PackStarProposalLearning.Data data;

        Fixture(List<int[]> edges, boolean sharedFill) {
            int positions = 9, n = 1 << positions;
            int[][] allowed = new int[positions][];
            for (int p = 0; p < positions; p++) allowed[p] = new int[]{0, 1};
            rcs = new RCs(allowed);
            graph = InteractionGraph.buildFromEdges(positions, edges);
            int[][] conf = new int[n][positions];
            double[] residual = new double[n];
            for (int i = 0; i < n; i++) {
                for (int p = 0; p < positions; p++) conf[i][p] = (i >>> p) & 1;
                residual[i] = 1.2 * parity(conf[i], 0, 1, 2)
                        + 0.8 * (sharedFill ? parity(conf[i], 0, 1, 3) : parity(conf[i], 3, 4, 5));
            }
            data = new PackStarProposalLearning.Data(conf, residual, new double[n]);
        }

        PackStarTripleEtaCorrections fit(int triples,
                PackStarTripleDecompositionCosts.Previewer preview,
                PackStarTripleDecompositionCosts.Limits limits) {
            return PackStarTripleEtaCorrections.emptyFitted(9).fitSelectedSecondMoment(
                    rcs, graph, data, new PackStarProposalLearning.Data[]{data, data},
                    1, triples, 3, 1, 0, 10, 1000, preview, limits);
        }
    }

    private static double parity(int[] conf, int a, int b, int c) {
        return (conf[a] ^ conf[b] ^ conf[c]) == 0 ? -1 : 1;
    }

    @Test void equalFillCountsCanHaveDifferentStructuralCostsAndChangeWinner() {
        List<int[]> path = new ArrayList<>();
        for (int p = 0; p < 8; p++) path.add(new int[]{p, p + 1});
        Fixture f = new Fixture(path, false);
        assertEquals("0,1,2", f.fit(1, null, null).positionTriples().get(0).toString());
        PackStarTripleDecompositionCosts.Previewer preview = fill ->
                cost(fill.contains(edge(0, 2)) ? 100000 : 100, 160);
        var selected = f.fit(1, preview, unlimited());
        assertEquals("3,4,5", selected.positionTriples().get(0).toString());
        assertEquals(1, selected.requiredFillEdges(f.graph).size());
    }

    @Test void hardLimitsRejectEvenAHighGainCandidateAndAllowPairOnly() {
        List<int[]> path = new ArrayList<>();
        for (int p = 0; p < 8; p++) path.add(new int[]{p, p + 1});
        Fixture f = new Fixture(path, false);
        var limits = new PackStarTripleDecompositionCosts.Limits(100, 160, 160, 160, 0);
        assertTrue(f.fit(3, fill -> cost(100, 320), limits).positionTriples().isEmpty());
        assertNull(limits.rejection(cost(100, 160))); // equality is admitted
        assertEquals("max-table-bytes", limits.rejection(cost(100, 176)));
        assertEquals("host-bytes", limits.rejection(new PackStarTripleDecompositionCosts.Cost(
                2, 10, 10, 160, 160, 161, 0, BigInteger.TEN)));
        assertEquals("file-bytes", limits.rejection(new PackStarTripleDecompositionCosts.Cost(
                2, 10, 10, 160, 160, 100, 1, BigInteger.TEN)));
    }

    @Test void selectedTriplesShareFillUnionAndDoNotPayForTheSameEdgeTwice() {
        List<int[]> edges = new ArrayList<>();
        for (int a = 0; a < 9; a++) for (int b = a + 1; b < 9; b++)
            if (!(a == 0 && b == 1)) edges.add(new int[]{a, b});
        Fixture f = new Fixture(edges, true);
        List<Set<Long>> requested = new ArrayList<>();
        var selected = f.fit(2, fill -> {
            requested.add(Set.copyOf(fill));
            return cost(100, 160);
        }, unlimited());
        Set<String> triples = new HashSet<>();
        for (var triple : selected.positionTriples()) triples.add(triple.toString());
        assertEquals(Set.of("0,1,2", "0,1,3"), triples);
        assertEquals(1, selected.requiredFillEdges(f.graph).size());
        assertTrue(requested.stream().allMatch(fill -> fill.isEmpty() || fill.equals(Set.of(edge(0, 1)))));
    }

    @Test void noStatisticalSignalDoesNotBuildAnyPreview() {
        var rcs = new RCs(new int[][]{{0, 1}, {0, 1}, {0, 1}});
        var graph = InteractionGraph.buildFromEdges(3, List.of(new int[]{0, 1}, new int[]{1, 2}));
        var data = new PackStarProposalLearning.Data(new int[][]{{0, 0, 0}, {1, 1, 1}},
                new double[2], new double[2]);
        var result = PackStarTripleEtaCorrections.emptyFitted(3).fitSelectedSecondMoment(
                rcs, graph, data, new PackStarProposalLearning.Data[]{data, data},
                1, 3, 3, 1, 0, 10, 1000,
                fill -> { fail("no signal must not trigger a decomposition"); return null; }, unlimited());
        assertTrue(result.positionTriples().isEmpty());
    }

    @Test void zeroGrowthHasFinitePriorityAndMemoryGrowthIsNotFreeWhenWorkDrops() {
        var base = cost(100, 160);
        assertEquals(0, base.growthFrom(base));
        assertEquals(9, cost(10, 1600).growthFrom(base));
        assertEquals(0.1, 1.0 / (1 + cost(10, 1600).growthFrom(base)), 1e-12);
        var mapped = new PackStarTripleDecompositionCosts.Cost(2, 10, 10, 160, 160,
                40, 160, BigInteger.valueOf(100));
        assertEquals(1, mapped.growthFrom(base));
    }

    @Test void realPreviewMatchesProductionRootingCachesUnionAndBoundsCacheSize() throws Exception {
        RCs rcs = new RCs(new int[][]{{0, 1}, {0, 1, 2}, {0, 1}, {0, 1}});
        var graph = InteractionGraph.buildFromEdges(4, List.of(new int[]{0, 1}, new int[]{1, 2}, new int[]{2, 3}));
        var cache = new PackStarTripleDecompositionCosts.Cache(rcs, graph, null, 1);
        var base = cache.preview(Set.of());
        var fill = new LinkedHashSet<>(List.of(edge(0, 2), edge(1, 3)));
        var preview = cache.preview(fill);
        assertSame(preview, cache.preview(new LinkedHashSet<>(List.of(edge(1, 3), edge(0, 2)))));
        assertEquals(2, cache.builds);
        List<int[]> edges = new ArrayList<>(graph.getEdgeList());
        edges.add(new int[]{0, 2}); edges.add(new int[]{1, 3});
        var production = PackStarTripleDecompositionCosts.rootProposalGraph(
                InteractionGraph.buildFromEdges(4, edges), rcs, null, true);
        RootedTreeNode root = production.selected.root;
        try {
            assertEquals(production.cost().toString(), preview.toString());
        } finally {
            RootedTreeEdge.postOrderReleaseLargeMemory(root);
        }
        var structural = PackStarTripleDecompositionCosts.rootProposalGraph(
                InteractionGraph.buildFromEdges(4, edges), rcs, null, false);
        root = structural.selected.root;
        try {
            List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
            RootedTreeEdge.collectLambdaEdges(root, lambdaEdges);
            for (RootedTreeEdge e : lambdaEdges) {
                assertFalse(e.hasDPTable());
                var field = RootedTreeEdge.class.getDeclaredField("enumeratedCount");
                field.setAccessible(true);
                assertNull(field.get(e));
                assertNull(e.getFullEnergyMin());
                assertNull(e.getFullEnergyRigid());
            }
        } finally {
            RootedTreeEdge.postOrderReleaseLargeMemory(root);
        }
        cache.preview(Set.of(edge(0, 2)));
        assertSame(base, cache.preview(Set.of())); // base is pinned
        cache.preview(fill);
        assertEquals(4, cache.builds); // the one-entry LRU evicted the first union
    }

    private static final class Properties implements AutoCloseable {
        private final java.util.Map<String, String> old = new java.util.HashMap<>();
        Properties(String... pairs) {
            for (int i = 0; i < pairs.length; i += 2) {
                old.put(pairs[i], System.getProperty(pairs[i]));
                System.setProperty(pairs[i], pairs[i + 1]);
            }
        }
        @Override public void close() {
            old.forEach((key, value) -> {
                if (value == null) System.clearProperty(key); else System.setProperty(key, value);
            });
        }
    }

    @Test void gpuRootOptimizationImprovesOverEdgeZeroAndMatchesMaterializedProposal() {
        try (var properties = new Properties("packstar.rootSplit", "gpubytes",
                "packstar.rootSplit.hostBudgetBytes", "1GiB",
                "packstar.rootSplit.gpuBudgetBytes", "1GiB")) {
            int[] counts = {2, 7, 3, 11, 2, 5};
            int[][] allowed = new int[counts.length][];
            for (int p = 0; p < counts.length; p++)
                allowed[p] = java.util.stream.IntStream.range(0, counts[p]).toArray();
            RCs rcs = new RCs(allowed);
            var edges = new ArrayList<int[]>();
            for (int p = 0; p < counts.length - 1; p++) edges.add(new int[]{p, p + 1});
            var base = InteractionGraph.buildFromEdges(counts.length, edges);
            var cache = new PackStarTripleDecompositionCosts.Cache(rcs, base, null, 8);
            var fill = Set.of(edge(0, 2), edge(2, 4), edge(1, 4));
            var preview = cache.preview(fill);
            edges.add(new int[]{0, 2}); edges.add(new int[]{1, 4}); edges.add(new int[]{2, 4});
            var graph = InteractionGraph.buildFromEdges(counts.length, edges);
            var production = PackStarTripleDecompositionCosts.rootProposalGraph(graph, rcs, null, true);
            try {
                assertEquals(production.cost().toString(), preview.toString());
                assertNotEquals(0, preview.rootSplitIndex, "regression: proposal must not silently use edge 0");
            } finally { RootedTreeEdge.postOrderReleaseLargeMemory(production.selected.root); }
            System.setProperty("packstar.rootSplit", "0");
            var legacy = PackStarTripleDecompositionCosts.rootProposalGraph(graph, rcs, null, false);
            try {
                assertTrue(preview.work.compareTo(legacy.cost().work) < 0,
                        "selected=" + preview + ", edge0=" + legacy.cost());
            } finally { RootedTreeEdge.postOrderReleaseLargeMemory(legacy.selected.root); }
        }
    }

    @Test void configuredExplicitRootAndHeapRejectionApplyToPreviewAndExecution() {
        try (var properties = new Properties("packstar.rootSplit", "1",
                "packstar.rootSplit.hostBudgetBytes", "1GiB")) {
            var rcs = new RCs(new int[][]{{0, 1}, {0, 1, 2}, {0, 1}, {0, 1}});
            var base = InteractionGraph.buildFromEdges(4,
                    List.of(new int[]{0, 1}, new int[]{1, 2}, new int[]{2, 3}));
            var graph = InteractionGraph.buildFromEdges(4,
                    List.of(new int[]{0, 1}, new int[]{1, 2}, new int[]{2, 3}, new int[]{0, 2}));
            var preview = new PackStarTripleDecompositionCosts.Cache(rcs, base, null, 2)
                    .preview(Set.of(edge(0, 2)));
            var production = PackStarTripleDecompositionCosts.rootProposalGraph(graph, rcs, null, true);
            try {
                assertEquals(1, preview.rootSplitIndex);
                assertEquals(production.cost().toString(), preview.toString());
            } finally { RootedTreeEdge.postOrderReleaseLargeMemory(production.selected.root); }
            System.setProperty("packstar.rootSplit.hostBudgetBytes", "1");
            var cache = new PackStarTripleDecompositionCosts.Cache(rcs, base, null, 2);
            var rejected = cache.preview(Set.of(edge(0, 2)));
            assertEquals("size-overflow", unlimited().rejection(rejected));
            assertSame(rejected, cache.preview(Set.of(edge(0, 2))));
            assertThrows(BranchDpBackend.RootSelectionException.class,
                    () -> PackStarTripleDecompositionCosts.rootProposalGraph(graph, rcs, null, true));
        }
    }

    @Test void mmapUsesAuxiliaryHeapAndHonorsPackstarAliasWithoutAllocatingFiles() {
        String key = "packstar.dp.tableMode";
        String old = System.getProperty(key);
        try {
            System.setProperty(key, "mmap");
            var rcs = new RCs(new int[][]{{0, 1}, {0, 1}, {0, 1}});
            var graph = InteractionGraph.buildFromEdges(3,
                    List.of(new int[]{0, 1}, new int[]{1, 2}, new int[]{0, 2}));
            var mapped = new PackStarTripleDecompositionCosts.Cache(rcs, graph, null, 1).preview(Set.of());
            assertTrue(mapped.fileBytes > 0);
            assertEquals(mapped.totalTableBytes, mapped.fileBytes);
            System.setProperty(key, "auto");
            var heap = new PackStarTripleDecompositionCosts.Cache(rcs, graph, null, 1).preview(Set.of());
            assertEquals(0, heap.fileBytes);
            assertEquals(heap.hostBytes, mapped.hostBytes + mapped.fileBytes);
        } finally {
            if (old == null) System.clearProperty(key); else System.setProperty(key, old);
        }
    }

    @Test void backendIndexingLimitIsCachedAsAnInfeasibleCandidate() {
        String key = "packstar.dp.maxMStates";
        String old = System.getProperty(key);
        try {
            System.setProperty(key, "1");
            var rcs = new RCs(new int[][]{{0, 1}, {0, 1}, {0, 1}});
            var graph = InteractionGraph.buildFromEdges(3,
                    List.of(new int[]{0, 1}, new int[]{1, 2}));
            var cache = new PackStarTripleDecompositionCosts.Cache(rcs, graph, null, 2);
            var rejected = cache.preview(Set.of(edge(0, 2)));
            assertEquals("size-overflow", unlimited().rejection(rejected));
            assertSame(rejected, cache.preview(Set.of(edge(0, 2))));
            assertEquals(1, cache.builds);
        } finally {
            if (old == null) System.clearProperty(key); else System.setProperty(key, old);
        }
    }
}
