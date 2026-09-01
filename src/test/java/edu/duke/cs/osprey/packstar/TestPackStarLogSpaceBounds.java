package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.tools.MathTools;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;

/**
 * Log-space PAC epsilon (design 1.3), the doc's remaining [TODO]:
 * "PAC log-space epsilon tests with logZ=+/-1e6".
 *
 * PackStarEstimator.epsilonFromLogBounds computes
 *   epsilon = 1 - exp(logZLower - logZUpper)
 * from the LOG ratio so it stays finite where the old direct-double path
 * (exp(logZ) then a BigDecimal ratio) collapsed to Z=0 and produced a fake
 * epsilon = NaN/1 (e.g. 5dc4). These tests pin that: epsilon depends only on the
 * gap between the bounds and is finite/in-range even at |logZ| ~ 1e6.
 */
public class TestPackStarLogSpaceBounds {

    private static double eps(double lo, double up) {
        return PackStarEstimator.epsilonFromLogBounds(lo, up);
    }

    @Test
    public void finiteAndCorrectAtExtremeMagnitudes() {
        // |logZ| ~ 1e6: exp(logZ) underflows/overflows a double, but the gap is small.
        // The bound INPUTS only resolve to ~ulp(1e6) ~ 1e-10, so the achievable
        // gap (and thus epsilon) carries ~1e-10 quantization noise. The point is
        // that epsilon stays finite and ~correct -- not the 0/1 collapse the old
        // exp(logZ) path produced -- so tolerate the input quantization (~1e-8).
        double e1 = eps(-1_000_000.0, -999_999.9);   // gap ~0.1
        assertTrue(Double.isFinite(e1));
        assertEquals(1.0 - Math.exp(-0.1), e1, 1e-8, "logZ ~ -1e6");

        double e2 = eps(999_999.5, 1_000_000.0);      // gap ~0.5
        assertTrue(Double.isFinite(e2));
        assertEquals(1.0 - Math.exp(-0.5), e2, 1e-8, "logZ ~ +1e6");
    }

    @Test
    public void dependsOnlyOnGap_underflowProof() {
        // Same gap at tiny vs huge magnitude => same epsilon. The OLD path would
        // give 0/0 = NaN at the huge magnitude; the log-space path does not.
        double gap = 0.25;
        double small = eps(-gap, 0.0);
        double huge = eps(-1_000_000.0 - gap, -1_000_000.0);
        assertEquals(small, huge, 1e-9, "epsilon must depend only on the bound gap");
        assertEquals(1.0 - Math.exp(-gap), huge, 1e-9);
    }

    @Test
    public void equalBoundsGiveZero() {
        assertEquals(0.0, eps(-1_000_000.0, -1_000_000.0), 0.0, "equal bounds -> eps 0");
        assertEquals(0.0, eps(0.0, 0.0), 0.0);
    }

    @Test
    public void tinyFpExcursionsAreSnapped() {
        // lower a hair above upper (FP noise) -> eps slightly negative -> snapped to 0
        double e = eps(Math.nextUp(5.0), 5.0);
        assertEquals(0.0, e, 0.0, "tiny negative epsilon must snap to exactly 0.0");
        assertTrue(e >= 0.0);
    }

    @Test
    public void hugeGapApproachesOneWithoutExceeding() {
        double e = eps(-1_000_000.0, 1_000_000.0); // exp(-2e6) underflows to 0
        assertEquals(1.0, e, 0.0, "huge gap -> exactly 1.0, never > 1");
        assertTrue(e <= 1.0);
    }

    @Test
    public void nonFiniteBoundsGiveEpsilonOne() {
        assertEquals(1.0, eps(Double.NEGATIVE_INFINITY, 5.0), 0.0);
        assertEquals(1.0, eps(5.0, Double.POSITIVE_INFINITY), 0.0);
        assertEquals(1.0, eps(Double.NaN, 5.0), 0.0);
        assertEquals(1.0, eps(5.0, Double.NaN), 0.0);
    }

    @Test
    public void alwaysInUnitInterval_forFiniteOrderedBounds() {
        double[] mags = {-1e6, -1e3, -1.0, 0.0, 1.0, 1e3, 1e6};
        double[] gaps = {0.0, 1e-6, 1e-3, 0.5, 5.0, 50.0, 1e6};
        for (double up : mags) {
            for (double g : gaps) {
                double lo = up - g; // lower <= upper
                double e = eps(lo, up);
                assertTrue(Double.isFinite(e), "finite at up=" + up + " gap=" + g);
                assertTrue(e >= 0.0 && e <= 1.0, "in [0,1] at up=" + up + " gap=" + g + " got " + e);
            }
        }
    }

    @Test
    public void materiallyReversedBoundsFailClosed() {
        assertEquals(1.0, eps(6.0, 5.0), 0.0);
    }

    @Test
    public void exactBulkTailPointwiseDecomposition() {
        double rt = 0.593050165;
        double[] gaps = {0.0, 0.1, 0.5, 2.0, 10.0};
        double[] etas = {-3.0, -0.2, 0.0, 0.4, 2.0, 8.0};
        double[] logCaps = {-2.0, -0.2, 0.0, 0.8, 3.0};
        for (double gap : gaps) {
            for (double eta : etas) {
                for (double logCap : logCaps) {
                    double x = PackStarEstimator.clippedBulkNormalized(
                            gap, eta, rt, logCap);
                    double h = PackStarEstimator.tailRemainderH(
                            gap, eta, rt, logCap);
                    double bulkRelativeToEm = Math.exp(logCap - eta / rt) * x;
                    double trueRelativeToEm = Math.exp(-gap / rt);
                    assertEquals(trueRelativeToEm, bulkRelativeToEm + h,
                            2e-13 * Math.max(1.0, trueRelativeToEm),
                            "pointwise identity gap=" + gap + " eta=" + eta
                                    + " logCap=" + logCap);
                    assertTrue(x >= 0.0 && x <= 1.0);
                    assertTrue(h >= 0.0 && h <= 1.0);
                    if (h > 0.0) {
                        assertTrue(eta > rt * logCap,
                                "positive tail must imply eta > kappa");
                    }
                }
            }
        }
    }

    @Test
    public void tiltedTailNormalizationIsAlwaysBounded() {
        double rt = 0.593050165;
        double[] lambdas = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0};
        for (double gap : new double[]{0.0, 0.01, 0.2, 1.0, 5.0}) {
            for (double eta : new double[]{-2.0, 0.0, 0.2, 0.8, 2.0, 10.0}) {
                for (double logCap : new double[]{-1.0, 0.0, 0.5, 2.0}) {
                    double h = PackStarEstimator.tailRemainderH(
                            gap, eta, rt, logCap);
                    for (double lambda : lambdas) {
                        double y = PackStarEstimator.tiltedTailNormalized(
                                gap, eta, rt, logCap, lambda);
                        assertTrue(Double.isFinite(y));
                        assertTrue(y >= 0.0 && y <= 1.0,
                                "normalized tail in [0,1]");
                        if (lambda == 0.0) {
                            assertEquals(h, y, 1e-15,
                                    "lambda=0 must be the p_m fallback");
                        }
                    }
                }
            }
        }
    }

    @Test
    public void logAddExpHandlesZeroComponentsAndExtremeLogs() {
        assertEquals(-3.0, PackStarEstimator.logAddExp(
                Double.NEGATIVE_INFINITY, -3.0), 0.0);
        assertEquals(1000.0 + Math.log1p(Math.exp(-2000.0)),
                PackStarEstimator.logAddExp(1000.0, -1000.0), 0.0);
        assertEquals(Math.log(2.0),
                PackStarEstimator.logAddExp(0.0, 0.0), 1e-15);
    }

    @Test
    public void unitIntervalSizingUsesAValidWorstCaseSampleVariance() {
        assertEquals(0.5,
                PackStarEstimator.maximumUnitIntervalSampleVariance(2), 0.0);
        assertEquals(1.0 / 3.0,
                PackStarEstimator.maximumUnitIntervalSampleVariance(4), 1e-15);
        assertTrue(PackStarEstimator.maximumUnitIntervalSampleVariance(256)
                > 0.25);
        assertTrue(PackStarEstimator.maximumUnitIntervalSampleVariance(256)
                < PackStarEstimator.maximumUnitIntervalSampleVariance(16));
    }

    @Test
    public void dynamicBHasClosedFormRadiusBalance() {
        double logZEta = 2.0;
        double logZLambda = 8.0;
        double lambda = 2.0;
        double bulkDelta = 0.1;
        double tailDelta = 0.2;
        PackStarEstimator.RadiusBalancedThreshold design =
                PackStarEstimator.balanceLogCap(
                        logZEta, logZLambda, lambda,
                        bulkDelta, tailDelta);
        assertNotNull(design);

        double expected = (Math.log(lambda)
                + logZLambda + Math.log(tailDelta)
                - logZEta - Math.log(bulkDelta)) / (lambda + 1.0);
        assertEquals(expected, design.logCap, 1e-15);

        double logBulkRadius = logZEta + Math.log(bulkDelta)
                + design.logCap;
        double logTailRadius = logZLambda + Math.log(tailDelta)
                - lambda * design.logCap;
        assertEquals(Math.log(lambda),
                logBulkRadius - logTailRadius, 1e-14,
                "the interior optimum satisfies bulkRadius=lambda*tailRadius");
    }

    @Test
    public void dynamicBUsesNonnegativeBoundaryForLambdaZeroOrSmallTail() {
        PackStarEstimator.RadiusBalancedThreshold lambdaZero =
                PackStarEstimator.balanceLogCap(
                        4.0, 10.0, 0.0, 0.2, 0.2);
        assertNotNull(lambdaZero);
        assertEquals(0.0, lambdaZero.logCap, 0.0);

        PackStarEstimator.RadiusBalancedThreshold boundary =
                PackStarEstimator.balanceLogCap(
                        10.0, 0.0, 1.0, 0.2, 0.2);
        assertNotNull(boundary);
        assertEquals(0.0, boundary.logCap, 0.0,
                "B is constrained to be nonnegative");
    }

    @Test
    public void dynamicBRejectsInvalidInputs() {
        assertNull(PackStarEstimator.balanceLogCap(
                Double.NaN, 1.0, 1.0, 0.1, 0.1));
        assertNull(PackStarEstimator.balanceLogCap(
                1.0, 1.0, -1.0, 0.1, 0.1));
        assertNull(PackStarEstimator.balanceLogCap(
                1.0, 1.0, 1.0, 0.0, 0.1));
    }

    @Test
    public void tighterTailIntervalImprovesJointSizingPrediction() {
        double loose = PackStarEstimator.combinedBulkTailEpsilon(
                0.0, 0.8, 1.0,
                0.0, 0.3, 0.2,
                Math.log(10.0));
        double tight = PackStarEstimator.combinedBulkTailEpsilon(
                0.0, 0.8, 1.0,
                0.0, 0.3, 0.05,
                Math.log(10.0));
        assertTrue(tight < loose,
                "more tail samples must improve the predicted combined epsilon");

        double capped = PackStarEstimator.combinedBulkTailEpsilon(
                0.0, 0.8, 1.0,
                0.0, 0.3, 0.2,
                Math.log(1.2));
        assertTrue(capped <= loose,
                "the deterministic q_m cap can only tighten the upper endpoint");
    }

    @Test
    public void certificateValidationFailsClosed() {
        assertTrue(PackStarEstimator.isValidCertificate(
                BigDecimal.ONE, BigDecimal.TEN, 0.9));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.ZERO, BigDecimal.ZERO, 0.0),
                "[0,0] cannot be an Estimated K* partition-function certificate");
        assertTrue(PackStarEstimator.isValidCertificate(
                BigDecimal.ZERO, BigDecimal.ONE, 1.0));
        assertFalse(PackStarEstimator.isValidCertificate(
                MathTools.BigNaN, BigDecimal.ONE, 0.5));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.ONE, MathTools.BigPositiveInfinity, 0.5));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.TEN, BigDecimal.ONE, -9.0));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.valueOf(-1), BigDecimal.ONE, 0.5));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.ONE, BigDecimal.TEN, Double.NaN));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.ONE, BigDecimal.TEN, -1e-3));
        assertFalse(PackStarEstimator.isValidCertificate(
                BigDecimal.ONE, BigDecimal.TEN, 1.001));
    }

    @Test
    public void finiteLogBoundsRetainArbitraryDecimalExponents() {
        for (double expectedLog10 : new double[]{-169.590919, -701.0}) {
            BigDecimal value = PackStarEstimator.bigExpFromLog(
                    expectedLog10 * Math.log(10.0));
            assertTrue(value.signum() > 0,
                    "finite logZ must never be represented as zero");
            assertEquals(expectedLog10, MathTools.log10(value), 1.0e-12);
        }
        assertEquals(0.0,
                PackStarEstimator.bigExpFromLog(
                        -701.0 * Math.log(10.0)).doubleValue(),
                0.0,
                "the BigDecimal must stay positive even when double underflows");
        assertEquals(BigDecimal.ZERO,
                PackStarEstimator.bigExpFromLog(
                        Double.NEGATIVE_INFINITY));
    }

    @Test
    public void conditionalResidualBoundGridIsCanonicalAndCappedBySelection() {
        assertArrayEquals(new double[]{1.0, 1.5, 2.0},
                PackStarEstimator.parsePositiveGrid(
                        "2, 1, 1.5, 1", "test B grid"), 0.0);

        double[] grid = {1.0, 1.5, 2.0};
        assertEquals(1.0,
                PackStarEstimator.selectSmallestResidualBound(-0.1, grid), 0.0);
        assertEquals(1.0,
                PackStarEstimator.selectSmallestResidualBound(1.0, grid), 0.0);
        assertEquals(1.5,
                PackStarEstimator.selectSmallestResidualBound(1.000001, grid), 0.0);
        assertEquals(1.5,
                PackStarEstimator.selectSmallestResidualBound(1.5, grid), 0.0);
        assertEquals(2.0,
                PackStarEstimator.selectSmallestResidualBound(1.8, grid), 0.0);
        assertTrue(Double.isNaN(
                PackStarEstimator.selectSmallestResidualBound(2.000001, grid)));
        assertTrue(Double.isNaN(
                PackStarEstimator.selectSmallestResidualBound(Double.NaN, grid)));
    }

    @Test
    public void conditionalResidualBoundGridRejectsInvalidValues() {
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parsePositiveGrid("", "test B grid"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parsePositiveGrid("0,1", "test B grid"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parsePositiveGrid("-1,1", "test B grid"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parsePositiveGrid("NaN,1", "test B grid"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parsePositiveGrid("bad,1", "test B grid"));
    }

    @Test
    public void conditionalClipAndExcessGiveExactBoundedWeightDecomposition() {
        double rt = 0.593050165;
        double boundKcal = 2.0;
        double logBound = boundKcal / rt;
        for (double logCap : new double[]{0.0, 0.5, 1.5, logBound}) {
            double excessRange = PackStarEstimator.conditionalExcessRange(
                    boundKcal, rt, logCap);
            assertTrue(Double.isFinite(excessRange));
            assertTrue(excessRange >= 0.0);
            for (double logWeight : new double[]{-10.0, -0.2, 0.0,
                    0.4, 1.2, logBound}) {
                double clipped = PackStarEstimator.conditionalClippedNormalized(
                        logWeight, logCap);
                double excess = PackStarEstimator.conditionalExcessNormalized(
                        logWeight, logCap);
                assertTrue(clipped >= 0.0 && clipped <= 1.0);
                assertTrue(excess >= 0.0 && excess <= excessRange + 1e-12);
                assertEquals(Math.exp(logWeight - logCap), clipped + excess,
                        2e-14 * Math.max(1.0,
                                Math.exp(logWeight - logCap)));
            }
        }
        assertEquals(Math.expm1(boundKcal / rt),
                PackStarEstimator.conditionalExcessRange(
                        boundKcal, rt, 0.0), 1e-14);
        assertEquals(0.0,
                PackStarEstimator.conditionalExcessRange(
                        boundKcal, rt, logBound), 1e-14);
    }

    @Test
    public void conditionalDpRichIdentityProposalPreservesOrdinaryStatistics() {
        double[] logRho = {0.0, 0.0, 0.0, 0.0};
        PackStarEstimator.ConditionalImportanceWeights normalized =
                PackStarEstimator.normalizeConditionalImportanceWeights(
                        logRho);
        assertNotNull(normalized);
        assertArrayEquals(new double[]{0.25, 0.25, 0.25, 0.25},
                normalized.normalized, 0.0);
        assertEquals(4.0, normalized.effectiveSampleSize, 0.0);

        double[] values = {1.0, 2.0, 4.0, 8.0};
        double mean = PackStarEstimator.conditionalWeightedMean(
                values, normalized.normalized);
        assertEquals(3.75, mean, 0.0);
        assertEquals(28.75 / 3.0,
                PackStarEstimator.conditionalWeightedSampleVariance(
                        values, normalized.normalized, mean), 1e-14);
        assertEquals(6.2,
                PackStarEstimator.conditionalWeightedQuantile(
                        values, normalized.normalized, 0.85), 1e-14,
                "equal weights must reproduce legacy i/(n-1) interpolation");
    }

    @Test
    public void conditionalDpRichTwoStateReweightingIsExact() {
        double rt = 1.0;
        double sourceLogZ = Math.log(2.0);
        double targetLogZ = Math.log(4.0 / 3.0);
        double[] logRho = {
                PackStarEstimator.conditionalImportanceLogRatio(
                        0.0, 0.0, sourceLogZ, targetLogZ, rt),
                PackStarEstimator.conditionalImportanceLogRatio(
                        0.0, Math.log(3.0),
                        sourceLogZ, targetLogZ, rt)
        };
        PackStarEstimator.ConditionalImportanceWeights normalized =
                PackStarEstimator.normalizeConditionalImportanceWeights(
                        logRho);
        assertNotNull(normalized);
        assertArrayEquals(new double[]{0.75, 0.25},
                normalized.normalized, 1e-15,
                "uniform source must reweight to the exact target Boltzmann law");
        assertEquals(0.25,
                PackStarEstimator.conditionalWeightedMean(
                        new double[]{0.0, 1.0}, normalized.normalized),
                1e-15);
        assertEquals(1.6, normalized.effectiveSampleSize, 1e-14);
    }

    @Test
    public void conditionalDpRichLogNormalizationStaysFinite() {
        PackStarEstimator.ConditionalImportanceWeights normalized =
                PackStarEstimator.normalizeConditionalImportanceWeights(
                        new double[]{1000.0, 999.0, -1000.0});
        assertNotNull(normalized);
        double sum = 0.0;
        for (double weight : normalized.normalized) {
            assertTrue(Double.isFinite(weight));
            assertTrue(weight >= 0.0 && weight <= 1.0);
            sum += weight;
        }
        assertEquals(1.0, sum, 1e-15);
        assertTrue(Double.isFinite(normalized.effectiveSampleSize));
        assertTrue(normalized.effectiveSampleSize > 1.0);
        assertTrue(normalized.effectiveSampleSize < 2.0);
    }

    @Test
    public void conditionalDpRichLowEssFailsClosed() {
        PackStarEstimator.ConditionalImportanceWeights collapsed =
                PackStarEstimator.normalizeConditionalImportanceWeights(
                        new double[]{1000.0, 0.0, 0.0, 0.0});
        PackStarEstimator.ConditionalImportanceWeights identity =
                PackStarEstimator.normalizeConditionalImportanceWeights(
                        new double[]{0.0, 0.0, 0.0, 0.0});
        assertNotNull(collapsed);
        assertNotNull(identity);
        assertFalse(PackStarEstimator.conditionalDpRichEssPasses(
                collapsed.effectiveSampleSize, 4, 0.25));
        assertTrue(PackStarEstimator.conditionalDpRichEssPasses(
                identity.effectiveSampleSize, 4, 0.25));
        assertFalse(PackStarEstimator.conditionalDpRichEssPasses(
                Double.NaN, 100, 0.25));
        assertFalse(PackStarEstimator.conditionalDpRichEssPasses(
                100.0, 100, 0.0));
    }

    @Test
    public void conditionalDpRichCandidateTieBreakIsDeterministic() {
        assertEquals(12.0,
                PackStarEstimator.conditionalLogRangeProxy(
                        10.0, 1.0, 0.5), 0.0);
        assertTrue(PackStarEstimator.isBetterConditionalDpRichCandidate(
                12.0, 1.0, 0.5,
                Double.NaN, Double.NaN, Double.NaN));
        assertTrue(PackStarEstimator.isBetterConditionalDpRichCandidate(
                11.0, 2.0, 0.1,
                12.0, 1.0, 1.0));
        assertTrue(PackStarEstimator.isBetterConditionalDpRichCandidate(
                12.0, 1.0, 0.1,
                12.0, 1.5, 1.0));
        assertTrue(PackStarEstimator.isBetterConditionalDpRichCandidate(
                12.0, 1.0, 0.9,
                12.0, 1.0, 0.8));
        assertFalse(PackStarEstimator.isBetterConditionalDpRichCandidate(
                12.0, 1.0, 0.7,
                12.0, 1.0, 0.8));
    }

    @Test
    public void conditionalDpRichZeroEtaIsIdentityAgainstMinProposal() {
        double sourceEnergy = -17.25;
        double sourceLogZ = 31.5;
        assertEquals(0.0,
                PackStarEstimator.conditionalImportanceLogRatio(
                        sourceEnergy, sourceEnergy,
                        sourceLogZ, sourceLogZ, 0.593050165),
                0.0,
                "eta=0 makes source and target both equal the E_m proposal");
        assertTrue(Double.isNaN(
                PackStarEstimator.conditionalImportanceLogRatio(
                        sourceEnergy, sourceEnergy,
                        sourceLogZ, sourceLogZ, 0.0)));
    }

    @Test
    public void conditionalEtaV3ChiSquareProxyIsGaugeInvariant() {
        double[] uniform = {0.5, 0.5};
        double expected = 1.0 / 9.0;
        assertEquals(expected,
                PackStarEstimator.conditionalWeightChiSquareProxy(
                        new double[]{0.0, Math.log(2.0)}, uniform),
                1e-15);
        assertEquals(expected,
                PackStarEstimator.conditionalWeightChiSquareProxy(
                        new double[]{-1000.0, -1000.0 + Math.log(2.0)},
                        uniform),
                1e-13,
                "a conformation-independent energy shift cannot change q or chi-square");
        assertEquals(0.0,
                PackStarEstimator.conditionalWeightChiSquareProxy(
                        new double[]{17.0, 17.0, 17.0},
                        new double[]{0.2, 0.3, 0.5}),
                1e-15);
    }

    @Test
    public void conditionalEtaV3CandidateOrderingRequiresReachability() {
        assertFalse(PackStarEstimator.isBetterConditionalEtaV3Candidate(
                false, 100, 0.0, 1.0, 0.1,
                -1, Double.NaN, Double.NaN, Double.NaN));
        assertTrue(PackStarEstimator.isBetterConditionalEtaV3Candidate(
                true, 800, 0.4, 10.0, 0.6,
                -1, Double.NaN, Double.NaN, Double.NaN));
        assertTrue(PackStarEstimator.isBetterConditionalEtaV3Candidate(
                true, 700, 100.0, 20.0, 0.65,
                800, 0.4, 10.0, 0.6),
                "predicted certificate cost is the primary objective");
        assertTrue(PackStarEstimator.isBetterConditionalEtaV3Candidate(
                true, 700, 0.2, 30.0, 0.65,
                700, 0.3, 10.0, 0.5),
                "chi-square breaks equal-cost ties");
        assertFalse(PackStarEstimator.isBetterConditionalEtaV3Candidate(
                true, 700, Double.NaN, 1.0, 0.5,
                -1, Double.NaN, Double.NaN, Double.NaN));
    }

    @Test
    public void conditionalEtaV3ProofBoundMustComeFromFrozenGrid() {
        double[] grid = {1.0, 1.5, 2.0};
        assertTrue(PackStarEstimator.gridContains(grid, 2.0));
        assertTrue(PackStarEstimator.gridContains(grid,
                Math.nextAfter(1.5, 2.0)));
        assertFalse(PackStarEstimator.gridContains(grid, 2.1));
        assertFalse(PackStarEstimator.gridContains(grid, Double.NaN));
    }
}
