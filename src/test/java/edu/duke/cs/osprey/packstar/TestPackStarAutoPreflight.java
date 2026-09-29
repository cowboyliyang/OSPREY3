package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.*;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.math.BigInteger;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.jupiter.api.Assertions.*;

class TestPackStarAutoPreflight {
    @AfterEach void clear() { BranchDpAdmission.clearExactPolicies(); }

    private static BranchDpAdmission.Prediction prediction(String key, BigInteger work, int drop) {
        return new BranchDpAdmission.Prediction("Complex", key, 8 - drop, 2,
                work, BigInteger.ZERO, true, 0,
                new BranchDpAdmission.Hardware(0, 0, 0), 1,
                8, work, Double.POSITIVE_INFINITY, drop > 0, drop > 0);
    }

    private static PackStarAutoPreflight.Settings settings(int rounds, int threads) {
        return new PackStarAutoPreflight.Settings(3, rounds, threads, .01, 5, 60);
    }

    @Test void noHoursOrHardwarePlateauStillTriesLargerBudgets() {
        AtomicInteger calls = new AtomicInteger(), checkpoints = new AtomicInteger();
        List<Long> budgets = new ArrayList<>();
        var state = new PackStarAutoPreflight.State(() -> {
            calls.incrementAndGet();
            var policy = BranchDpAdmission.getExactPolicy("same");
            if (policy != null) budgets.add(policy.maxMillis);
            return prediction("same", BigInteger.valueOf(1000), 0);
        });
        var result = PackStarAutoPreflight.optimize(List.of(state), settings(12, 1),
                false, checkpoints::incrementAndGet);
        assertEquals("plateau", result.stopReason);
        assertEquals(3, result.rounds);
        assertEquals(4, calls.get());
        assertEquals(List.of(5L, 20L, 60L), budgets);
        assertEquals(4, checkpoints.get());
        assertFalse(result.summary.withinSla());
        assertFalse(state.best.hasFinitePrediction());
        assertNull(BranchDpAdmission.getExactPolicy("same"));
    }

    @Test void retainsBestDespiteRegressionAndReplaysWorkCeiling(@TempDir Path dir) {
        AtomicInteger calls = new AtomicInteger();
        var state = new PackStarAutoPreflight.State(() -> {
            int call = calls.incrementAndGet();
            return prediction("winner", BigInteger.valueOf(call == 2 ? 500 : 1000), call == 2 ? 1 : 0);
        });
        var result = PackStarAutoPreflight.optimize(List.of(state), settings(12, 1), false, () -> {});
        assertEquals(BigInteger.valueOf(500), result.finalWork);
        assertEquals(4, result.rounds);
        var file = dir.resolve("policy.tsv").toFile();
        BranchDpAdmission.writeExactPolicies(file);
        BranchDpAdmission.clearExactPolicies();
        BranchDpAdmission.loadExactPolicies(file);
        assertEquals(BigInteger.valueOf(500), BranchDpAdmission.getExactPolicy("winner").maxTotalGpuWork);
        assertDoesNotThrow(() -> BranchDpAdmission.enforceRetainedPredictionCeiling(state.best, "test"));
        assertThrows(IllegalStateException.class, () -> BranchDpAdmission.enforceRetainedPredictionCeiling(
                prediction("winner", BigInteger.valueOf(501), 0), "test"));
    }

    @Test void smallImprovementsAreSavedButDoNotKeepSearchAliveForever() {
        AtomicInteger calls = new AtomicInteger();
        var state = new PackStarAutoPreflight.State(() -> prediction("tiny",
                BigInteger.valueOf(1000 - calls.getAndIncrement()), 1));
        var result = PackStarAutoPreflight.optimize(List.of(state), settings(12, 1), false, () -> {});
        assertEquals("plateau", result.stopReason);
        assertEquals(BigInteger.valueOf(997), result.finalWork);
        assertEquals(3, result.rounds);
    }

    @Test void keepsOptimizingAllStatesAndReportsSearchLimit() {
        List<PackStarAutoPreflight.State> states = new ArrayList<>();
        for (String key : List.of("a", "b", "c")) {
            AtomicInteger calls = new AtomicInteger();
            states.add(new PackStarAutoPreflight.State(() -> prediction(key,
                    BigInteger.valueOf(1000 >> calls.getAndIncrement()), 1)));
        }
        var result = PackStarAutoPreflight.optimize(states, settings(2, 3), false, () -> {});
        assertEquals("search-limit", result.stopReason);
        assertEquals(BigInteger.valueOf(750), result.finalWork);
        assertEquals(3, BranchDpAdmission.exactPolicyCount());
    }

    @Test void failedAndWrongIdentityCandidatesRollbackWithoutFalseCleanPlateau() {
        AtomicInteger calls = new AtomicInteger();
        var state = new PackStarAutoPreflight.State(() -> {
            int call = calls.incrementAndGet();
            if (call == 1) return prediction("original", BigInteger.TEN, 0);
            if (call == 2) return prediction("other", BigInteger.ONE, 1);
            throw new IllegalStateException("failed trial");
        });
        var result = PackStarAutoPreflight.optimize(List.of(state), settings(12, 1), false, () -> {});
        assertEquals("plateau-with-preview-errors", result.stopReason);
        assertEquals(BigInteger.TEN, result.finalWork);
        assertEquals(0, BranchDpAdmission.exactPolicyCount());
    }

    @Test void initialFailureAndCheckpointFailureAreNotSwallowed() {
        var broken = new PackStarAutoPreflight.State(() -> null);
        assertThrows(IllegalStateException.class, () -> PackStarAutoPreflight.optimize(
                List.of(broken), settings(12, 1), false, () -> {}));
        var valid = new PackStarAutoPreflight.State(() -> prediction("ok", BigInteger.ONE, 0));
        assertThrows(IllegalStateException.class, () -> PackStarAutoPreflight.optimize(
                List.of(valid), settings(12, 1), false,
                () -> { throw new IllegalStateException("disk full"); }));
    }

    @Test void lockedPolicyIsOnlyPreviewedOnce() {
        AtomicInteger calls = new AtomicInteger();
        var policy = new BranchDpAdmission.ExactPolicy(1, 1, 10).withWorkCeiling(BigInteger.TEN);
        BranchDpAdmission.putExactPolicy("locked", policy);
        var state = new PackStarAutoPreflight.State(() -> {
            calls.incrementAndGet();
            return prediction("locked", BigInteger.TEN, 1);
        });
        var result = PackStarAutoPreflight.optimize(List.of(state), settings(12, 1), true, () -> {});
        assertEquals("loaded-policy", result.stopReason);
        assertEquals(1, calls.get());
        assertSame(policy, BranchDpAdmission.getExactPolicy("locked"));
    }

    @Test void hugeWorkDoesNotOverflowConvergenceMetric() {
        BigInteger before = BigInteger.TEN.pow(500);
        assertEquals(.5, PackStarAutoPreflight.relativeGain(before, before.divide(BigInteger.TWO)), 0);
    }

    @Test void automaticRootIsExhaustiveAndIndependentOfThroughput() {
        var saved = new java.util.HashMap<String, String>();
        String[] keys = {"packstar.admission.mode", "packstar.rootSplit",
                "packstar.admission.gpuWorkPerSecondPerGpu", "packstar.admission.gpuCount"};
        for (String key : keys) saved.put(key, System.getProperty(key));
        try (var scope = BranchDpConfig.enterPackStarAliasScope()) {
            int[] counts = {2, 7, 3, 11, 2, 5};
            int[][] allowed = new int[counts.length][];
            for (int p = 0; p < counts.length; p++)
                allowed[p] = java.util.stream.IntStream.range(0, counts[p]).toArray();
            RCs rcs = new RCs(allowed);
            var edges = new ArrayList<int[]>();
            for (int p = 0; p < counts.length - 1; p++) edges.add(new int[]{p, p + 1});
            var graph = InteractionGraph.buildFromEdges(counts.length, edges);
            var tree = new BranchDecomposition(graph, BranchDecomposition.Strategy.WEIGHTED_HICKS, counts);
            tree.compute();
            System.setProperty("packstar.admission.mode", "sla");
            BigInteger expected = null;
            for (int root = 0; root < tree.getTree().getNumEdges(); root++) {
                System.setProperty("packstar.rootSplit", Integer.toString(root));
                var candidate = BranchDpBackend.selectConfiguredRoot(tree, graph, rcs, null, false);
                if (expected == null || candidate.gpuWork.compareTo(expected) < 0) expected = candidate.gpuWork;
                RootedTreeEdge.postOrderReleaseLargeMemory(candidate.root);
            }
            System.setProperty("packstar.admission.mode", "auto");
            for (String rate : List.of("0", "9999999999999")) {
                System.setProperty("packstar.admission.gpuWorkPerSecondPerGpu", rate);
                System.setProperty("packstar.admission.gpuCount", "4");
                var candidate = BranchDpBackend.selectConfiguredRoot(tree, graph, rcs, null, false);
                assertEquals(expected, candidate.gpuWork);
                List<RootedTreeEdge> lambda = new ArrayList<>();
                RootedTreeEdge.collectLambdaEdges(candidate.root, lambda);
                for (var edge : lambda) assertFalse(edge.hasDPTable());
                RootedTreeEdge.postOrderReleaseLargeMemory(candidate.root);
            }
        } finally {
            saved.forEach((key, value) -> {
                if (value == null) System.clearProperty(key); else System.setProperty(key, value);
            });
        }
    }
}
