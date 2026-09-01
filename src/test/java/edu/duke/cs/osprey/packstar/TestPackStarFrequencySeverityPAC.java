package edu.duke.cs.osprey.packstar;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

public class TestPackStarFrequencySeverityPAC {

    @Test
    public void clopperPearsonZeroCountUsesExactClosedForm() {
        int n = 100;
        double delta = 0.025;
        double expected = -Math.expm1(Math.log(delta) / n);
        assertEquals(expected,
                PackStarFrequencySeverityPAC.clopperPearsonUpper(
                        0, n, delta), 1.0e-15);
    }

    @Test
    public void frequencySeverityIntervalMatchesItsDecomposition() {
        double[] logRelative = {
                Math.log(1.0), Math.log(2.0),
                Math.log(4.0), Math.log(0.5)
        };
        double logClip = Math.log(2.0);
        double severity = 3.0;
        double bulkDelta = 0.025;
        double frequencyDelta = 0.025;

        PackStarFrequencySeverityPAC.Interval interval =
                PackStarFrequencySeverityPAC.evaluate(
                        logRelative, logClip, severity,
                        bulkDelta, frequencyDelta);

        // Y=min(R/C,1) = [0.5,1,1,0.25].  Only R=4 is a strict tail.
        assertEquals(4, interval.sampleCount);
        assertEquals(1, interval.tailCount);
        assertEquals(0.6875, interval.bulkMean, 0.0);
        assertEquals(0.25, interval.tailProbabilityEmpirical, 0.0);
        assertEquals(0.25, interval.empiricalTailMean, 1.0e-15);
        assertEquals(1.0, interval.empiricalConditionalSeverity, 1.0e-15);
        assertEquals(1.0, interval.observedMaxConditionalSeverity, 1.0e-15);
        assertEquals(interval.tailProbabilityUpper * severity,
                interval.inducedTailMeanUpper, 0.0);
        assertEquals(interval.bulkUpper + interval.inducedTailMeanUpper,
                interval.normalizedMeanUpper, 0.0);
        assertEquals(1.0 - interval.normalizedMeanLower
                        / interval.normalizedMeanUpper,
                interval.epsilon, 1.0e-15);
    }

    @Test
    public void addingARawEnergyGaugeDoesNotChangeRelativeInterval() {
        double[] rawLogWeights = {-1.5, -0.2, 0.1, 1.7, 2.2};
        double logMu = 0.3;
        double gauge = 19.0;
        double[] relativeBefore = new double[rawLogWeights.length];
        double[] relativeAfter = new double[rawLogWeights.length];
        for (int i = 0; i < rawLogWeights.length; i++) {
            relativeBefore[i] = rawLogWeights[i] - logMu;
            relativeAfter[i] = rawLogWeights[i] + gauge
                    - (logMu + gauge);
        }
        PackStarFrequencySeverityPAC.Interval before =
                PackStarFrequencySeverityPAC.evaluate(
                        relativeBefore, 1.0, 20.0, 0.025, 0.025);
        PackStarFrequencySeverityPAC.Interval after =
                PackStarFrequencySeverityPAC.evaluate(
                        relativeAfter, 1.0, 20.0, 0.025, 0.025);
        assertEquals(before.tailCount, after.tailCount);
        assertEquals(before.bulkMean, after.bulkMean, 1.0e-15);
        assertEquals(before.bulkLower, after.bulkLower, 1.0e-15);
        assertEquals(before.normalizedMeanUpper,
                after.normalizedMeanUpper, 1.0e-15);
        assertEquals(before.epsilon, after.epsilon, 1.0e-15);
    }

    @Test
    public void pooledCrossfitMomentsRecoverUniformCombinedSample() {
        PackStarFrequencySeverityPAC.PooledCrossfitMoments pooled =
                PackStarFrequencySeverityPAC.poolCrossfitMoments(
                        new int[]{2, 2},
                        new double[]{2.0, 2.0},
                        new double[]{0.5, 0.5},
                        new double[]{0.5, 0.5},
                        new double[]{0.0, 1.0});
        assertEquals(4, pooled.sampleCount);
        assertEquals(4.0, pooled.effectiveSampleSize, 1.0e-15);
        assertEquals(1.0, pooled.effectiveSampleFraction, 1.0e-15);
        assertEquals(0.5, pooled.moments.bulkMean, 1.0e-15);
        assertEquals(1.0 / 3.0,
                pooled.moments.bulkVariance, 1.0e-15);
        assertEquals(0.5, pooled.moments.tailProbability, 1.0e-15);
    }

    @Test
    public void pooledCrossfitEssAccountsForFoldWeightConcentration() {
        PackStarFrequencySeverityPAC.PooledCrossfitMoments pooled =
                PackStarFrequencySeverityPAC.poolCrossfitMoments(
                        new int[]{100, 100},
                        new double[]{25.0, 100.0},
                        new double[]{0.4, 0.6},
                        new double[]{0.02, 0.02},
                        new double[]{0.1, 0.3});
        assertEquals(80.0, pooled.effectiveSampleSize, 1.0e-12);
        assertEquals(0.4, pooled.effectiveSampleFraction, 1.0e-15);
        assertEquals(0.5, pooled.moments.bulkMean, 1.0e-15);
        assertEquals(0.2, pooled.moments.tailProbability, 1.0e-15);
        assertThrows(IllegalArgumentException.class,
                () -> PackStarFrequencySeverityPAC.poolCrossfitMoments(
                        new int[]{10}, new double[]{11.0},
                        new double[]{0.5}, new double[]{0.1},
                        new double[]{0.0}));
    }

    @Test
    public void sizingUsesTheUnreachableCapWhenTargetCannotBeReached() {
        PackStarFrequencySeverityPAC.Moments collapsed =
                new PackStarFrequencySeverityPAC.Moments(
                        0.01, 0.2, 0.5);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        collapsed, 4000, 400,
                        0.683, 0.9, 20.0,
                        0.025, 0.025);
        assertFalse(sizing.reachableAtMax);
        assertEquals(400, sizing.finalSamples);
        assertTrue(sizing.epsilonAtMaxSamples > 0.683);
    }

    @Test
    public void sizingFindsAReachableFrozenCount() {
        PackStarFrequencySeverityPAC.Moments stable =
                new PackStarFrequencySeverityPAC.Moments(
                        0.8, 0.01, 0.0);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        stable, 4000, 400,
                        0.683, 0.9, 1.0,
                        0.025, 0.025);
        assertTrue(sizing.reachableAtMax);
        assertTrue(sizing.finalSamples >= 2);
        assertTrue(sizing.finalSamples <= 4000);
        assertTrue(sizing.epsilonAtFinalSamples
                <= 0.683 * 0.9 + 1.0e-12);
    }

    @Test
    public void independentSeverityMonitorCanRejectButNotValidate() {
        double[] hugeTail = new double[10];
        for (int i = 0; i < hugeTail.length; i++) {
            hugeTail[i] = Math.log1p(100.0);
        }
        PackStarFrequencySeverityPAC.SeverityTest rejected =
                PackStarFrequencySeverityPAC.testConditionalSeverity(
                        hugeTail, 0.0, 1.0, 0.05);
        assertTrue(rejected.sufficientTailSamples);
        assertTrue(rejected.rejected);
        assertTrue(rejected.pValueUpper <= 0.05);

        double[] noTail = {-2.0, -1.0, -0.5};
        PackStarFrequencySeverityPAC.SeverityTest inconclusive =
                PackStarFrequencySeverityPAC.testConditionalSeverity(
                        noTail, 0.0, 1.0, 0.05);
        assertFalse(inconclusive.sufficientTailSamples);
        assertFalse(inconclusive.rejected);
    }

    @Test
    public void productionCandidateScalesAreFrozenAtOne() {
        PackStarEstimator.FrequencySeverityShrinkPair[] shrink =
                PackStarEstimator.parseFrequencySeverityShrinkGrid(
                        "0:0,2:5,5:10");
        assertEquals(3, shrink.length);
        assertEquals(0.0, shrink[0].unary, 0.0);
        assertEquals(0.0, shrink[0].pair, 0.0);
        assertArrayEquals(new double[]{1.0},
                PackStarEstimator.parseFixedOneGrid(
                        "1", "test grid"), 0.0);
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parseFrequencySeverityShrinkGrid(
                        "2:5"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parseFixedOneGrid(
                        "0,1", "test grid"));
        assertThrows(IllegalArgumentException.class,
                () -> PackStarEstimator.parseFixedOneGrid(
                        "0.999", "test grid"));
    }

    @Test
    public void candidateIdentityRecordsTheFixedOneScales() {
        PackStarEstimator.FrequencySeverityShrinkPair raw =
                new PackStarEstimator.FrequencySeverityShrinkPair(0.0, 0.0);
        assertEquals("ku-0.00000-kp-0.00000-alpha-1.00000-pair-only",
                PackStarEstimator.frequencySeverityCandidateId(
                        raw, 1.0, 0.0));
        assertEquals(
                "ku-0.00000-kp-0.00000-alpha-1.00000-plus-triple-eta-gamma-1.00000",
                PackStarEstimator.frequencySeverityCandidateId(
                        raw, 1.0, 1.0));
    }

    @Test
    public void sourceAwareWeightRetainsTheSourceNormalizer() {
        double rt = 0.6;
        double targetEnergy = 13.0;
        double sourceEnergy = 11.5;
        double sourceLogZ = 27.0;
        assertEquals(-(targetEnergy - sourceEnergy) / rt + sourceLogZ,
                PackStarEstimator.frequencySeveritySourceLogWeight(
                        targetEnergy, sourceEnergy, sourceLogZ, rt),
                0.0);
    }

    @Test
    public void sourceAwareWeightIsInvariantToASourceEnergyGauge() {
        double rt = 0.6;
        double targetEnergy = 13.0;
        double sourceEnergy = 11.5;
        double sourceLogZ = 27.0;
        double gauge = 4.25;
        double before = PackStarEstimator.frequencySeveritySourceLogWeight(
                targetEnergy, sourceEnergy, sourceLogZ, rt);
        // Adding g to every source energy changes logZ_source by -g/RT.
        double after = PackStarEstimator.frequencySeveritySourceLogWeight(
                targetEnergy, sourceEnergy + gauge,
                sourceLogZ - gauge / rt, rt);
        assertEquals(before, after, 1.0e-14);
    }

    @Test
    public void zeroBulkLowerWithReachablePointSizingExtendsFrozenDiscovery() {
        double logClip = 1.0;
        double[] logRelative = new double[100];
        for (int i = 0; i < logRelative.length; i++) {
            // Twenty observations at the clipping boundary and eighty nearly
            // zero bulk values mimic a high-variance, tail-free pilot.
            logRelative[i] = i < 20
                    ? logClip : logClip + Math.log(1.0e-12);
        }
        PackStarFrequencySeverityPAC.Interval interval =
                PackStarFrequencySeverityPAC.evaluate(
                        logRelative, logClip, 20.0, 0.025, 0.025);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        new PackStarFrequencySeverityPAC.Moments(
                                interval.bulkMean,
                                interval.bulkVariance * 1.3,
                                interval.tailProbabilityEmpirical),
                        4000, 400, 0.683, 0.9, 20.0,
                        0.025, 0.025);

        assertFalse(interval.hasPositiveBulkLower());
        assertTrue(sizing.reachableAtMax);
        assertTrue(PackStarEstimator.shouldExtendFrequencySeverityDiscovery(
                100, 400, interval, sizing, 0.683));
        assertFalse(PackStarEstimator.shouldExtendFrequencySeverityDiscovery(
                400, 400, interval, sizing, 0.683));
    }

    @Test
    public void discoveryExtensionDoesNotMaskAnUnreachableProposal() {
        double[] logRelative = new double[100];
        Arrays.fill(logRelative, -20.0);
        PackStarFrequencySeverityPAC.Interval interval =
                PackStarFrequencySeverityPAC.evaluate(
                        logRelative, 1.0, 20.0, 0.025, 0.025);
        PackStarFrequencySeverityPAC.Sizing unreachable =
                PackStarFrequencySeverityPAC.size(
                        new PackStarFrequencySeverityPAC.Moments(
                                0.01, 0.2, 0.5),
                        4000, 400, 0.683, 0.9, 20.0,
                        0.025, 0.025);

        assertFalse(unreachable.reachableAtMax);
        assertFalse(PackStarEstimator.shouldExtendFrequencySeverityDiscovery(
                100, 400, interval, unreachable, 0.683));
    }
}
