package edu.duke.cs.osprey.branchdp;

import edu.duke.cs.osprey.astar.conf.RCs;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class TestBranchDpAdmission {

    @AfterEach
    public void clearPolicies() {
        BranchDpAdmission.clearExactPolicies();
    }

    @Test
    public void predictionAddsComputeAndOutOfCoreTime() {
        BranchDpAdmission.Hardware hardware =
                new BranchDpAdmission.Hardware(100.0, 2, 50.0);
        BranchDpAdmission.Prediction prediction = prediction(
                "Complex", "state", BigInteger.valueOf(20_000),
                BigInteger.valueOf(5_000), true, 0, hardware);

        assertEquals(100.0, prediction.computeSeconds, 1e-12);
        assertEquals(100.0, prediction.oocSeconds, 1e-12);
        assertEquals(200.0, prediction.predictedSeconds, 1e-12);
    }

    @Test
    public void predictionMultipliesBothWorkAndTrafficByDpSweeps() {
        BranchDpAdmission.Hardware hardware =
                new BranchDpAdmission.Hardware(100.0, 2, 50.0);
        BranchDpAdmission.Prediction prediction =
                new BranchDpAdmission.Prediction(
                        "Complex", "state", 12, 0,
                        BigInteger.valueOf(20_000), BigInteger.valueOf(5_000),
                        true, 0, hardware, 6,
                        12, BigInteger.valueOf(20_000), 200.0,
                        false, false);

        assertEquals(6, prediction.dpSweeps);
        assertEquals(BigInteger.valueOf(120_000), prediction.totalGpuWork());
        assertEquals(BigInteger.valueOf(30_000), prediction.totalOocTrafficBytes());
        assertEquals(1200.0, prediction.predictedSeconds, 1e-12);
    }

    @Test
    public void safetyFactorConservativelyInflatesBothLegs() {
        BranchDpAdmission.Hardware hardware =
                new BranchDpAdmission.Hardware(100.0, 2, 50.0, 1.25);
        BranchDpAdmission.Prediction prediction = prediction(
                "Complex", "state", BigInteger.valueOf(20_000),
                BigInteger.valueOf(5_000), true, 0, hardware);

        assertEquals(125.0, prediction.computeSeconds, 1e-12);
        assertEquals(125.0, prediction.oocSeconds, 1e-12);
        assertEquals(250.0, prediction.predictedSeconds, 1e-12);
        assertEquals(1.25, prediction.safetyFactor, 0.0);
    }

    @Test
    public void missingCalibrationOrTrafficFailsClosed() {
        BranchDpAdmission.Hardware noComputeRate =
                new BranchDpAdmission.Hardware(0.0, 8, 1_000.0);
        assertFalse(prediction("Complex", "a", BigInteger.ONE,
                BigInteger.ZERO, true, 0, noComputeRate)
                .hasFinitePrediction());

        BranchDpAdmission.Hardware noTrafficRate =
                new BranchDpAdmission.Hardware(100.0, 8, 0.0);
        assertFalse(prediction("Complex", "b", BigInteger.ONE,
                BigInteger.ONE, true, 0, noTrafficRate)
                .hasFinitePrediction());
        assertFalse(prediction("Complex", "c", BigInteger.ONE,
                BigInteger.ZERO, false, 0, noTrafficRate)
                .hasFinitePrediction());
        assertFalse(prediction("Complex", "d", BigInteger.ONE,
                BigInteger.ZERO, true, 1, noTrafficRate)
                .hasFinitePrediction());
    }

    @Test
    public void wholeCaseSlaIsSumNotPerStateLimit() {
        BranchDpAdmission.Hardware hardware =
                new BranchDpAdmission.Hardware(1.0, 1, 1.0);
        BranchDpAdmission.Prediction first = prediction(
                "Protein", "p", BigInteger.valueOf(2000),
                BigInteger.ZERO, true, 0, hardware);
        BranchDpAdmission.Prediction second = prediction(
                "Complex", "c", BigInteger.valueOf(2000),
                BigInteger.ZERO, true, 0, hardware);

        // Each state is below one hour, but the complete case is above one hour.
        assertTrue(first.predictedHours() < 1.0);
        assertTrue(second.predictedHours() < 1.0);
        BranchDpAdmission.CaseSummary rejected =
                new BranchDpAdmission.CaseSummary(
                        Arrays.asList(first, second), 1.0);
        assertEquals(4000.0 / 3600.0, rejected.totalHours(), 1e-12);
        assertFalse(rejected.withinSla());

        BranchDpAdmission.CaseSummary admitted =
                new BranchDpAdmission.CaseSummary(
                        Arrays.asList(first, second), 2.0);
        assertTrue(admitted.withinSla());
    }

    @Test
    public void contributorsAreSortedByWholeStatePrediction() {
        BranchDpAdmission.Hardware hardware =
                new BranchDpAdmission.Hardware(1.0, 1, 1.0);
        BranchDpAdmission.Prediction small = prediction(
                "small", "s", BigInteger.ONE, BigInteger.ZERO,
                true, 0, hardware);
        BranchDpAdmission.Prediction large = prediction(
                "large", "l", BigInteger.TEN, BigInteger.ZERO,
                true, 0, hardware);

        List<BranchDpAdmission.Prediction> sorted =
                new BranchDpAdmission.CaseSummary(
                        Arrays.asList(small, large), 1.0)
                        .contributorsDescending();
        assertSame(large, sorted.get(0));
        assertSame(small, sorted.get(1));
    }

    @Test
    public void exactPolicyUsesExactRcIdentityAndRetainsCeiling() {
        RCs first = new RCs(new int[][] {{0, 2}, {1}});
        RCs second = new RCs(new int[][] {{0, 3}, {1}});
        String firstKey = BranchDpAdmission.stateKey("Complex", first);
        String secondKey = BranchDpAdmission.stateKey("Complex", second);
        assertNotEquals(firstKey, secondKey);

        BranchDpAdmission.ExactPolicy policy =
                new BranchDpAdmission.ExactPolicy(2, 2, 1234)
                        .withPredictionCeiling(42.0);
        BranchDpAdmission.putExactPolicy(firstKey, policy);
        assertSame(policy, BranchDpAdmission.getExactPolicy(firstKey));
        assertNull(BranchDpAdmission.getExactPolicy(secondKey));
        assertEquals(1, BranchDpAdmission.exactPolicyCount());
        assertEquals(42.0, policy.maxPredictedSeconds, 0.0);
    }

    @Test
    public void exactImproveOptionsClampToBoundedValidRange() {
        BranchDecomposition.ExactImproveOptions options =
                new BranchDecomposition.ExactImproveOptions(3, 1, 0);
        assertEquals(3, options.minDrop);
        assertEquals(3, options.maxDrop);
        assertEquals(1L, options.maxMillis);
    }

    @Test
    public void exactPolicyDumpIsKeySortedAndRoundTrips(
            @TempDir Path tempDir) throws IOException {
        BranchDpAdmission.putExactPolicy("z-state",
                new BranchDpAdmission.ExactPolicy(2, 3, 200L, 20.0));
        BranchDpAdmission.putExactPolicy("a-state",
                new BranchDpAdmission.ExactPolicy(1, 1, 100L, 10.0));
        File policyFile = tempDir.resolve("policy.tsv").toFile();

        assertEquals(2,
                BranchDpAdmission.writeExactPolicies(policyFile));
        List<String> rows = new ArrayList<>();
        for (String line : Files.readAllLines(policyFile.toPath())) {
            if (!line.isEmpty() && !line.startsWith("#")) rows.add(line);
        }
        assertEquals(2, rows.size());
        assertTrue(rows.get(0).endsWith("\ta-state"));
        assertTrue(rows.get(1).endsWith("\tz-state"));

        BranchDpAdmission.clearExactPolicies();
        assertEquals(2,
                BranchDpAdmission.loadExactPolicies(policyFile));
        assertEquals(1,
                BranchDpAdmission.getExactPolicy("a-state").minDrop);
        assertEquals(3,
                BranchDpAdmission.getExactPolicy("z-state").maxDrop);
    }

    private static BranchDpAdmission.Prediction prediction(
            String stateName, String key, BigInteger work,
            BigInteger traffic, boolean trafficAvailable,
            int unsupportedEdges,
            BranchDpAdmission.Hardware hardware) {
        return new BranchDpAdmission.Prediction(
                stateName, key, 10, 0, work, traffic,
                trafficAvailable, unsupportedEdges, hardware, 1,
                10, work, hardware.totalSeconds(work, traffic,
                trafficAvailable, unsupportedEdges), false, false);
    }
}
