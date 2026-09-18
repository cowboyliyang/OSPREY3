package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import org.junit.jupiter.api.Test;
import java.math.BigInteger;
import java.util.*;
import static org.junit.jupiter.api.Assertions.*;

class TestPackStarForwardTripleSelection {
    private static PackStarTripleDecompositionCosts.Cost cost(long work) {
        return new PackStarTripleDecompositionCosts.Cost(2, 8, 24, 128, 384, 384, 0, BigInteger.valueOf(work));
    }
    private static PackStarTripleDecompositionCosts.Limits limits(long work) {
        return new PackStarTripleDecompositionCosts.Limits(1000, 10000, 10000, 10000, 0, BigInteger.valueOf(work));
    }
    private static class Fixture {
        final RCs rcs;
        final InteractionGraph graph;
        final PackStarProposalLearning.Data data;
        Fixture(double third, boolean missingEdges) {
            int positions = 9, n = 1 << positions;
            int[][] allowed = new int[positions][];
            for (int p = 0; p < positions; p++) allowed[p] = new int[]{0, 1};
            rcs = new RCs(allowed);
            List<int[]> edges = new ArrayList<>();
            for (int b = 0; b < positions; b += 3) {
                edges.add(new int[]{b, b+1}); edges.add(new int[]{b+1, b+2});
                if (!missingEdges) edges.add(new int[]{b, b+2});
            }
            graph = InteractionGraph.buildFromEdges(positions, edges);
            int[][] conf = new int[n][positions];
            double[] r = new double[n];
            for (int i = 0; i < n; i++) {
                for (int p = 0; p < positions; p++) conf[i][p] = (i >>> p) & 1;
                r[i] = parity(conf[i], 0) + 0.7 * parity(conf[i], 3) + third * parity(conf[i], 6);
            }
            data = new PackStarProposalLearning.Data(conf, r, new double[n]);
        }
        PackStarTripleEtaCorrections.MomentPath fit(int k, int fill,
                PackStarTripleDecompositionCosts.Previewer preview, long work) {
            return PackStarTripleEtaCorrections.emptyFitted(9).fitSecondMomentPath(rcs, graph,
                    data, new PackStarProposalLearning.Data[]{data, data}, 1, k, fill, 1, 0, 3, 1000,
                    preview, limits(work));
        }
    }
    private static double parity(int[] c, int p) { return (c[p] ^ c[p+1] ^ c[p+2]) == 0 ? -1 : 1; }
    private static long edge(int a, int b) { return ((long) a << 32) | b; }

    @Test void retainsEveryPrefixAndStopsWhenRemainingCliqueHasNoSignal() {
        Fixture f = new Fixture(0, false);
        var path = f.fit(3, 0, fill -> cost(100), 1000);
        assertEquals(3, path.models.size(), path.stopReason); // K=0,1,2
        assertEquals("no-feasible-reliable-improvement", path.stopReason);
        for (int k = 0; k < path.models.size(); k++) assertEquals(k, path.models.get(k).positionTripleCount);
        assertSame(path.models.get(2), path.atMost(3));
        double[] h = new double[f.data.conf.length];
        for (int i = 0; i < h.length; i++) h[i] = path.models.get(2).scoreResidual(f.data.conf[i], (a, ra, b, rb) -> 0);
        assertEquals(0, PackStarProposalLearning.logRho(f.data, h, -1, 1), 1e-7);
        assertNotSame(path.models.get(1).entries.get(0), path.models.get(2).entries.get(0));
    }

    @Test void canRetainAllFourModels() {
        var path = new Fixture(0.45, false).fit(3, 0, fill -> cost(100), 1000);
        assertEquals(4, path.models.size(), path.stopReason);
        assertEquals("maximum-K", path.stopReason);
    }

    @Test void gainRanksFeasibleCandidatesAndWorkBudgetIsAHardGate() {
        Fixture f = new Fixture(0, true);
        PackStarTripleDecompositionCosts.Previewer preview = fill -> cost(fill.contains(edge(0, 2)) ? 100000 : 100);
        assertEquals("0,1,2", f.fit(1, 1, preview, 100000).models.get(1).positionTriples().get(0).toString());
        assertEquals("3,4,5", f.fit(1, 1, preview, 1000).models.get(1).positionTriples().get(0).toString());
        assertEquals(1, f.fit(3, 1, fill -> cost(1001), 1000).models.size());
        assertNull(limits(100).rejection(cost(100)));
        assertEquals("work", limits(100).rejection(cost(101)));
    }

    @Test void jointFitCoordinatesOverlappingFactorsAndHonorsCap() {
        // Two overlapping factors have correlated features. Fitting the first
        // alone absorbs part of the second; joint fitting must revise it.
        int[][] features = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};
        double[] r = {-2, 0, 0, 2}, a = {Math.log(4), 0, 0, Math.log(4)};
        var first = PackStarProposalLearning.jointMomentFit(new int[][]{{0},{0},{1},{1}},
                r, a, new double[2], new double[2], 1, 3, 1000);
        double[] sequential = {first.parameters[0], first.parameters[1], -1, 1};
        var fit = PackStarProposalLearning.jointMomentFit(features, r, a, new double[4], sequential, 1, 3, 1000);
        assertTrue(fit.converged);
        assertEquals(0, fit.objective, 1e-9);
        assertTrue(Math.abs(first.parameters[0] - fit.parameters[0]) > 0.05);
        var capped = PackStarProposalLearning.jointMomentFit(features, r, a, new double[4], new double[4], 1, 0.2, 1000);
        assertTrue(capped.converged);
        for (double h : capped.parameters) assertTrue(Math.abs(h) <= 0.2);
        assertTrue(capped.objective > fit.objective + 0.01);
    }

    @Test void jointObjectiveGradientAndSourceMultiplicityAgree() {
        int[][] f = {{0,2}, {1,2}, {1,3}};
        double[] r = {-1, 0.5, 2}, a = {Math.log(2), 0, 0};
        double[] penalty = {0.1, 0.2, 0.3, 0.4}, h = {-0.4, 0.1, 0.2, 0.5}, g = new double[4];
        PackStarProposalLearning.momentObjective(f, r, a, penalty, h, 0.6, g);
        for (int j = 0; j < h.length; j++) {
            double[] plus = h.clone(), minus = h.clone(); plus[j] += 1e-5; minus[j] -= 1e-5;
            double derivative = (PackStarProposalLearning.momentObjective(f, r, a, penalty, plus, 0.6, null)
                    - PackStarProposalLearning.momentObjective(f, r, a, penalty, minus, 0.6, null)) / 2e-5;
            assertEquals(derivative, g[j], 1e-8);
        }
        var weighted = PackStarProposalLearning.jointMomentFit(f, r, a, penalty, h, 0.6, 3, 1000);
        var repeated = PackStarProposalLearning.jointMomentFit(new int[][]{{0,2},{0,2},{1,2},{1,3}},
                new double[]{-1,-1,0.5,2}, new double[]{10000,10000,10000,10000}, penalty, h, 0.6, 3, 1000);
        assertTrue(weighted.converged); assertTrue(repeated.converged);
        assertArrayEquals(weighted.parameters, repeated.parameters, 1e-5);
    }

    @Test void contradictoryInnerEvidenceAndInsufficientSupportStopBeforeAdding() {
        Fixture f = new Fixture(0, false);
        PackStarProposalLearning.Data[] contradictory = new PackStarProposalLearning.Data[2];
        for (int fold = 0; fold < 2; fold++) {
            double[] residual = f.data.residual.clone();
            for (int i = 0; i < residual.length; i++)
                if (f.data.folds[i] == fold) residual[i] = -residual[i];
            contradictory[fold] = new PackStarProposalLearning.Data(f.data.conf, residual, f.data.logWeight);
        }
        var path = PackStarTripleEtaCorrections.emptyFitted(9).fitSecondMomentPath(f.rcs, f.graph,
                f.data, contradictory, 1, 3, 0, 1, 0, 3, 1000, fill -> cost(100), limits(1000));
        assertEquals(1, path.models.size());
        var unsupported = PackStarTripleEtaCorrections.emptyFitted(9).fitSecondMomentPath(f.rcs, f.graph,
                f.data, new PackStarProposalLearning.Data[]{f.data, f.data}, 1, 3, 0,
                1000, 4, 3, 1000, fill -> {
                    fail("unsupported candidates must not build a decomposition"); return null;
                }, limits(1000));
        assertEquals(1, unsupported.models.size());
        var one = f.fit(1, 0, fill -> cost(100), 1000).models.get(1);
        assertEquals(0, one.scoreResidual(new int[]{2,0,0,0,0,0,0,0,0}, (a,ra,b,rb) -> 0));
    }

    @Test void fivePercentUsesGlobalMinimumAndPrefersSmallerK() {
        record Model(int k, long n, int structure) {}
        Comparator<Model> simpler = Comparator.comparingInt(Model::k).thenComparingInt(Model::structure);
        Model pair = new Model(0, 600, 0), triple = new Model(3, 590, 10);
        assertSame(pair, PackStarProposalLearning.chooseWithinFivePercent(List.of(triple, pair), Model::n, simpler));
        Model good = new Model(3, 700, 10);
        assertSame(good, PackStarProposalLearning.chooseWithinFivePercent(List.of(new Model(0, 2000, 0), good), Model::n, simpler));
        // 105 and 110 are close, but 110 is outside 5% of the global best 100.
        Model best = new Model(3, 100, 3), middle = new Model(2, 105, 2), bad = new Model(0, 110, 0);
        List<Model> order = new ArrayList<>(List.of(best, middle, bad));
        for (int i = 0; i < 6; i++) {
            Collections.shuffle(order, new Random(i));
            assertSame(middle, PackStarProposalLearning.chooseWithinFivePercent(order, Model::n, simpler));
        }
        assertNull(PackStarProposalLearning.chooseWithinFivePercent(List.<Model>of(), Model::n, simpler));
    }

    @Test void parallelFitsPreserveEveryPrefixAndPreviewOrder() {
        String key = "packstar.pac.frequencySeverity.tripleFitThreads";
        String previous = System.getProperty(key);
        try {
            Fixture fixture = new Fixture(0.45, true);
            List<Set<Long>> serialPreviews = new ArrayList<>(), parallelPreviews = new ArrayList<>();
            System.setProperty(key, "1");
            var serial = fixture.fit(3, 3, fill -> {
                serialPreviews.add(new LinkedHashSet<>(fill));
                return cost(fill.contains(edge(0, 2)) ? 200 : 100);
            }, 1000);
            System.setProperty(key, "4");
            var parallel = fixture.fit(3, 3, fill -> {
                parallelPreviews.add(new LinkedHashSet<>(fill));
                return cost(fill.contains(edge(0, 2)) ? 200 : 100);
            }, 1000);
            assertEquals(serialPreviews, parallelPreviews);
            assertEquals(serial.stopReason, parallel.stopReason);
            assertEquals(serial.gains, parallel.gains);
            assertEquals(serial.models.size(), parallel.models.size());
            for (int k = 0; k < serial.models.size(); k++) {
                assertEquals(serial.models.get(k).signatureSha256, parallel.models.get(k).signatureSha256);
                for (int[] conf : fixture.data.conf)
                    assertEquals(serial.models.get(k).scoreResidual(conf, (a,ra,b,rb) -> 0),
                            parallel.models.get(k).scoreResidual(conf, (a,ra,b,rb) -> 0), 0.0);
            }
        } finally {
            if (previous == null) System.clearProperty(key); else System.setProperty(key, previous);
        }
    }

    @Test void limitingFoldWorkPreservesEveryRequiredPrefix() {
        Fixture fixture = new Fixture(0.45, true);
        var complete = fixture.fit(3, 3, fill -> cost(100), 1000);
        assertEquals(4, complete.models.size());
        for (int required = 0; required < 3; required++) {
            final boolean mayPreview = required > 0;
            var limited = fixture.fit(required, 3, fill -> {
                assertTrue(mayPreview, "K=0 must not fit or preview any triple");
                return cost(100);
            }, 1000);
            assertEquals(required + 1, limited.models.size());
            for (int k = 0; k <= required; k++) {
                assertEquals(complete.models.get(k).signatureSha256, limited.models.get(k).signatureSha256);
                assertEquals(complete.gains.get(k), limited.gains.get(k));
            }
        }
    }

    @Test void shiftedMomentSumsAgreeWithSequentialReference() {
        Random random = new Random(179);
        for (int trial = 0; trial < 20; trial++) {
            int n = 300 + trial * 17, p = 37;
            int[][] features = new int[n][];
            double[] residual = new double[n], weights = new double[n];
            double[] penalty = new double[p], point = new double[p], actualGradient = new double[p];
            double rt = 0.6;
            for (int j = 0; j < p; j++) { penalty[j] = random.nextDouble(); point[j] = 6 * random.nextDouble() - 3; }
            for (int i = 0; i < n; i++) {
                features[i] = new int[i % 4];
                for (int j = 0; j < features[i].length; j++) features[i][j] = random.nextInt(p);
                residual[i] = 60 * random.nextDouble() - 30;
                weights[i] = 10000 - 1000 * random.nextDouble();
            }
            // Independent oracle: the original sequential logAdd reduction.
            double shift = Arrays.stream(weights).max().orElseThrow();
            double residualShift = residual[0] / rt;
            double logA = Double.NEGATIVE_INFINITY, logB = logA, first = logA;
            double[] a = new double[n], b = new double[n], expectedGradient = new double[p];
            for (int i = 0; i < n; i++) {
                double h = 0;
                for (int j : features[i]) h += point[j] / rt;
                double source = weights[i] - shift, r = residual[i] / rt - residualShift;
                a[i] = source - h; b[i] = source - 2 * r + h;
                logA = PackStarProposalLearning.logAdd(logA, a[i]);
                logB = PackStarProposalLearning.logAdd(logB, b[i]);
                first = PackStarProposalLearning.logAdd(first, source - r);
            }
            double expected = logA + logB - 2 * first;
            for (int i = 0; i < n; i++) {
                double g = (Math.exp(b[i] - logB) - Math.exp(a[i] - logA)) / rt;
                for (int j : features[i]) expectedGradient[j] += g;
            }
            for (int j = 0; j < p; j++) {
                expected += 0.5 * penalty[j] * Math.pow(point[j] / rt, 2);
                expectedGradient[j] += penalty[j] * point[j] / (rt * rt);
            }
            double actual = new PackStarProposalLearning.MomentObjective(features, residual, weights, penalty, rt, true)
                    .value(point, actualGradient);
            assertEquals(expected, actual, 1e-10);
            assertArrayEquals(expectedGradient, actualGradient, 1e-10);
        }
    }

    @Test void objectiveScratchCanBeReusedAcrossBacktrackingPoints() {
        int[][] features = {{0, 2}, {1, 2}, {1, 3}};
        double[] residual = {-1, 0.5, 2}, weights = {10000, 9999, 10001};
        double[] penalty = {0.1, 0.2, 0.3, 0.4};
        double[][] points = {{-.4, .1, .2, .5}, {2, -2, 0, 1}, {-.4, .1, .2, .5}};
        for (boolean shifted : new boolean[]{false, true}) {
            var cached = new PackStarProposalLearning.MomentObjective(features, residual, weights, penalty, 0.6, shifted);
            double tolerance = shifted ? 1e-12 : 0.0;
            for (double[] point : points) {
                double[] expectedGradient = new double[4], actualGradient = new double[4];
                Arrays.fill(actualGradient, Double.NaN);
                double expected = PackStarProposalLearning.momentObjective(
                        features, residual, weights, penalty, point, 0.6, expectedGradient);
                assertEquals(expected, cached.value(point, null), tolerance);
                cached.gradientAtLastValue(point, actualGradient);
                assertArrayEquals(expectedGradient, actualGradient, tolerance);
                // Repeated gradient reads must not consume or mutate cached masses.
                cached.gradientAtLastValue(point, actualGradient);
                assertArrayEquals(expectedGradient, actualGradient, tolerance);
                assertEquals(expected, cached.value(point, actualGradient), tolerance);
                assertArrayEquals(expectedGradient, actualGradient, tolerance);
            }
        }
    }
}
