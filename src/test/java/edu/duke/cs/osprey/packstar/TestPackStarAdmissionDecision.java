package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestPackStarAdmissionDecision {

    private static final BranchDpAdmission.Hardware HARDWARE =
            new BranchDpAdmission.Hardware(1.0, 1, 1.0);

    @BeforeEach
    public void clearPoliciesBeforeTest() {
        BranchDpAdmission.clearExactPolicies();
    }

    @AfterEach
    public void clearPoliciesAfterTest() {
        BranchDpAdmission.clearExactPolicies();
    }

    @Test
    public void initialWholeCasePassAdmitsWithoutFinalRetries() {
        AtomicInteger firstPreviews = new AtomicInteger();
        AtomicInteger secondPreviews = new AtomicInteger();
        List<PackStarAdmissionDecision.State> states = List.of(
                new PackStarAdmissionDecision.State("first-sequence", () -> {
                    firstPreviews.incrementAndGet();
                    return prediction("Protein", "protein|first",
                            0.30, 10, 10, false, false);
                }),
                new PackStarAdmissionDecision.State("second-sequence", () -> {
                    secondPreviews.incrementAndGet();
                    return prediction("Complex", "complex|second",
                            0.45, 11, 11, false, false);
                }));

        PackStarAdmissionDecision.Outcome outcome =
                PackStarAdmissionDecision.decide(states, settings(1.0, 3));

        assertTrue(outcome.finalSummary.withinSla());
        assertEquals(0.75, outcome.finalSummary.totalHours(), 1.0e-12);
        assertSame(outcome.initialSummary, outcome.finalSummary);
        assertFalse(outcome.usedFinalPass());
        assertEquals(0, outcome.attemptedFinalStates);
        assertEquals(0, outcome.acceptedFinalStates);
        assertEquals(1, firstPreviews.get());
        assertEquals(1, secondPreviews.get());
        assertEquals(0, BranchDpAdmission.exactPolicyCount());
    }

    @Test
    public void finalPassImprovesLargestContributorAndAdmitsCase() {
        String heavyKey = "complex|heavy";
        AtomicInteger heavyPreviews = new AtomicInteger();
        AtomicInteger smallPreviews = new AtomicInteger();
        PackStarAdmissionDecision.State heavy =
                new PackStarAdmissionDecision.State("heavy-sequence", () -> {
                    int call = heavyPreviews.incrementAndGet();
                    if (call == 1) {
                        assertNull(BranchDpAdmission.getExactPolicy(heavyKey));
                        return prediction("Complex", heavyKey,
                                2.0, 12, 12, false, false);
                    }
                    BranchDpAdmission.ExactPolicy trial =
                            BranchDpAdmission.getExactPolicy(heavyKey);
                    assertNotNull(trial);
                    assertEquals(1, trial.minDrop);
                    assertEquals(1, trial.maxDrop);
                    assertEquals(54_321L, trial.maxMillis);
                    assertTrue(Double.isInfinite(trial.maxPredictedSeconds));
                    return prediction("Complex", heavyKey,
                            0.40, 11, 12, true, true);
                });
        PackStarAdmissionDecision.State small =
                new PackStarAdmissionDecision.State("small-sequence", () -> {
                    smallPreviews.incrementAndGet();
                    return prediction("Protein", "protein|small",
                            0.40, 9, 9, false, false);
                });

        PackStarAdmissionDecision.Outcome outcome =
                PackStarAdmissionDecision.decide(
                        List.of(heavy, small), settings(1.0, 2));

        assertFalse(outcome.initialSummary.withinSla());
        assertEquals(2.40, outcome.initialSummary.totalHours(), 1.0e-12);
        assertTrue(outcome.finalSummary.withinSla());
        assertEquals(0.80, outcome.finalSummary.totalHours(), 1.0e-12);
        assertTrue(outcome.usedFinalPass());
        assertEquals(1, outcome.attemptedFinalStates);
        assertEquals(1, outcome.acceptedFinalStates);
        assertEquals(2, heavyPreviews.get());
        assertEquals(1, smallPreviews.get());
        assertEquals(0.40, heavy.prediction().predictedHours(), 1.0e-12);

        BranchDpAdmission.ExactPolicy retained =
                BranchDpAdmission.getExactPolicy(heavyKey);
        assertNotNull(retained);
        assertEquals(1, retained.minDrop);
        assertEquals(1, retained.maxDrop);
        assertEquals(54_321L, retained.maxMillis);
        assertEquals(0.40 * 3600.0,
                retained.maxPredictedSeconds, 1.0e-12);
    }

    @Test
    public void distributedShardForcesFinalPassEvenWhenLocalShardFits() {
        String stateKey = "complex|distributed-shard";
        AtomicInteger previews = new AtomicInteger();
        PackStarAdmissionDecision.State state =
                new PackStarAdmissionDecision.State("distributed-sequence", () -> {
                    int call = previews.incrementAndGet();
                    if (call == 1) {
                        return prediction("Complex", stateKey,
                                0.25, 12, 12, false, false);
                    }
                    assertNotNull(BranchDpAdmission.getExactPolicy(stateKey));
                    return prediction("Complex", stateKey,
                            0.10, 11, 12, true, true);
                });

        PackStarAdmissionDecision.Settings distributedSettings =
                new PackStarAdmissionDecision.Settings(
                        1.0, 1, 54_321L, 1, 1,
                        false, true,
                        new BranchDecomposition.ExactImproveOptions(
                                1, 1, 12_345L));
        PackStarAdmissionDecision.Outcome outcome =
                PackStarAdmissionDecision.decide(
                        List.of(state), distributedSettings);

        assertTrue(outcome.initialSummary.withinSla());
        assertTrue(outcome.finalSummary.withinSla());
        assertTrue(outcome.usedFinalPass());
        assertEquals(1, outcome.attemptedFinalStates);
        assertEquals(1, outcome.acceptedFinalStates);
        assertEquals(2, previews.get());
    }

    @Test
    public void loadedPolicyIsPreviewedOnceWithoutPolicyResearch() {
        String stateKey = "complex|frozen";
        BranchDpAdmission.ExactPolicy frozen =
                new BranchDpAdmission.ExactPolicy(2, 2, 54_321L, 900.0);
        BranchDpAdmission.putExactPolicy(stateKey, frozen);
        AtomicInteger previews = new AtomicInteger();
        PackStarAdmissionDecision.State state =
                new PackStarAdmissionDecision.State("frozen-sequence", () -> {
                    previews.incrementAndGet();
                    assertSame(frozen,
                            BranchDpAdmission.getExactPolicy(stateKey));
                    return prediction("Complex", stateKey,
                            0.20, 10, 12, true, true);
                });

        PackStarAdmissionDecision.Settings lockedSettings =
                new PackStarAdmissionDecision.Settings(
                        1.0, 3, 54_321L, 2, 1,
                        true, false,
                        new BranchDecomposition.ExactImproveOptions(
                                1, 1, 12_345L));
        PackStarAdmissionDecision.Outcome outcome =
                PackStarAdmissionDecision.decide(
                        List.of(state), lockedSettings);

        assertTrue(outcome.finalSummary.withinSla());
        assertFalse(outcome.usedFinalPass());
        assertEquals(1, previews.get());
        assertSame(frozen, BranchDpAdmission.getExactPolicy(stateKey));
    }

    @Test
    public void loadedPolicyOverSlaRejectsWithoutPolicyResearch() {
        String stateKey = "complex|frozen-over-sla";
        BranchDpAdmission.ExactPolicy frozen =
                new BranchDpAdmission.ExactPolicy(1, 1, 54_321L, 7200.0);
        BranchDpAdmission.putExactPolicy(stateKey, frozen);
        AtomicInteger previews = new AtomicInteger();
        PackStarAdmissionDecision.State state =
                new PackStarAdmissionDecision.State("frozen-heavy", () -> {
                    previews.incrementAndGet();
                    return prediction("Complex", stateKey,
                            2.0, 12, 12, true, true);
                });
        PackStarAdmissionDecision.Settings lockedSettings =
                new PackStarAdmissionDecision.Settings(
                        1.0, 3, 54_321L, 2, 1,
                        true, false,
                        new BranchDecomposition.ExactImproveOptions(
                                1, 1, 12_345L));

        IllegalStateException rejected = assertThrows(
                IllegalStateException.class,
                () -> PackStarAdmissionDecision.decide(
                        List.of(state), lockedSettings));

        assertTrue(rejected.getMessage().contains(
                "loaded exact-policy plan still exceeds SLA"));
        assertEquals(1, previews.get());
        assertSame(frozen, BranchDpAdmission.getExactPolicy(stateKey));
    }

    @Test
    public void stillOverSlaAbortsBeforeFormalMaterialization() {
        String heavyKey = "complex|unimproved";
        AtomicInteger heavyPreviews = new AtomicInteger();
        AtomicInteger formalMaterializations = new AtomicInteger();
        PackStarAdmissionDecision.State heavy =
                new PackStarAdmissionDecision.State("heavy-sequence", () -> {
                    heavyPreviews.incrementAndGet();
                    return prediction("Complex", heavyKey,
                            2.0, 12, 12, false, false);
                });
        PackStarAdmissionDecision.State small =
                new PackStarAdmissionDecision.State("small-sequence", () ->
                        prediction("Protein", "protein|small",
                                0.20, 9, 9, false, false));

        IllegalStateException rejected = assertThrows(
                IllegalStateException.class, () -> {
                    PackStarAdmissionDecision.decide(
                            List.of(heavy, small), settings(1.0, 1));
                    formalMaterializations.incrementAndGet();
                });

        assertTrue(rejected.getMessage().contains(
                "rejected before DP-table materialization"));
        assertTrue(rejected.getMessage().contains(
                "whole case still exceeds SLA"));
        assertTrue(rejected.getMessage().contains("complex|unimproved")
                || rejected.getMessage().contains("heavy-sequence"));
        assertEquals(2, heavyPreviews.get());
        assertEquals(0, formalMaterializations.get());
        assertEquals(0, BranchDpAdmission.exactPolicyCount(),
                "a rejected exact candidate must restore the previous policy");
    }

    @Test
    public void acceptedExactPlanIsLockedAndSlowerReproductionFailsBeforeAllocation() {
        String stateKey = "complex|soft-exact";
        BranchDecomposition.ExactImproveOptions initialOptions =
                new BranchDecomposition.ExactImproveOptions(1, 2, 116_000L);
        PackStarAdmissionDecision.State state =
                new PackStarAdmissionDecision.State("soft-exact-sequence", () ->
                        prediction("Complex", stateKey,
                                0.25, 11, 12, true, true));

        PackStarAdmissionDecision.Outcome outcome =
                PackStarAdmissionDecision.decide(List.of(state),
                        new PackStarAdmissionDecision.Settings(
                                1.0, 3, 54_321L, initialOptions));

        assertTrue(outcome.finalSummary.withinSla());
        BranchDpAdmission.ExactPolicy retained =
                BranchDpAdmission.getExactPolicy(stateKey);
        assertNotNull(retained);
        assertEquals(1, retained.minDrop);
        assertEquals(2, retained.maxDrop);
        assertEquals(116_000L, retained.maxMillis);
        assertEquals(0.25 * 3600.0,
                retained.maxPredictedSeconds, 1.0e-12);

        BranchDpAdmission.Prediction reproduced = prediction(
                "Complex", stateKey, 0.25, 11, 12, true, true);
        assertDoesNotThrow(() ->
                BranchDpAdmission.enforceRetainedPredictionCeiling(
                        reproduced, "PACK*:"));

        AtomicInteger formalMaterializations = new AtomicInteger();
        BranchDpAdmission.Prediction slower = prediction(
                "Complex", stateKey, 0.50, 11, 12, true, true);
        IllegalStateException mismatch = assertThrows(
                IllegalStateException.class, () -> {
                    BranchDpAdmission.enforceRetainedPredictionCeiling(
                            slower, "PACK*:");
                    formalMaterializations.incrementAndGet();
                });
        assertTrue(mismatch.getMessage().contains(
                "retained whole-case decomposition could not be reproduced"));
        assertTrue(mismatch.getMessage().contains(
                "before DP-table materialization"));
        assertEquals(0, formalMaterializations.get());
    }

    private static PackStarAdmissionDecision.Settings settings(
            double slaHours, int finalMaxStates) {
        return new PackStarAdmissionDecision.Settings(
                slaHours, finalMaxStates, 54_321L,
                new BranchDecomposition.ExactImproveOptions(
                        1, 1, 12_345L));
    }

    private static BranchDpAdmission.Prediction prediction(
            String stateName, String stateKey, double hours,
            int branchwidth, int firstStageBranchwidth,
            boolean adaptiveAttempted, boolean adaptiveAccepted) {
        long seconds = Math.round(hours * 3600.0);
        BigInteger work = BigInteger.valueOf(seconds);
        return new BranchDpAdmission.Prediction(
                stateName, stateKey, branchwidth, 0,
                work, BigInteger.ZERO, true, 0,
                HARDWARE, 1,
                firstStageBranchwidth, work, seconds,
                adaptiveAttempted, adaptiveAccepted);
    }
}
