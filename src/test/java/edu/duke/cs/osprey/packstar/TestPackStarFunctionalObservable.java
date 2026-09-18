package edu.duke.cs.osprey.packstar;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class TestPackStarFunctionalObservable {

    @Test
    public void selfNormalizedProbabilityUsesFinalImportanceWeights() {
        PackStarFunctionalObservableAccumulator accumulator =
                new PackStarFunctionalObservableAccumulator("contact", 1.0);
        accumulator.add(true, Math.log(3.0));
        accumulator.add(false, 0.0);

        PackStarFunctionalObservableResult result = accumulator.finish();

        assertEquals(PackStarFunctionalObservableResult.Status.RESOLVED,
                result.getStatus());
        assertEquals(0.75, result.getProbability(), 1.0e-15);
        assertEquals(-Math.log(0.75),
                result.getRestrictionFreeEnergy(), 1.0e-15);
        assertEquals(1.6, result.getEffectiveSampleSize(), 1.0e-15);
        assertEquals(0.75, result.getMaxNormalizedWeight(), 1.0e-15);
        assertEquals(1, result.getHitCount());
        assertEquals(0.5, result.getRawHitFraction(), 1.0e-15);
        assertTrue(result.isResolved());
    }

    @Test
    public void zeroHitIsReportedAsUnresolved() {
        PackStarFunctionalObservableAccumulator accumulator =
                new PackStarFunctionalObservableAccumulator("rare", 0.6);
        accumulator.add(false, 0.0);
        accumulator.add(false, Math.log(4.0));

        PackStarFunctionalObservableResult result = accumulator.finish();

        assertEquals(PackStarFunctionalObservableResult.Status.UNRESOLVED_NO_HITS,
                result.getStatus());
        assertEquals(0.0, result.getProbability(), 0.0);
        assertTrue(Double.isNaN(result.getRestrictionFreeEnergy()));
        assertEquals(25.0 / 17.0, result.getEffectiveSampleSize(), 1.0e-15);
        assertFalse(result.isResolved());
    }

    @Test
    public void invalidWeightsDoNotProduceAResolvedFreeEnergy() {
        PackStarFunctionalObservableAccumulator accumulator =
                new PackStarFunctionalObservableAccumulator("bad", 0.6);
        accumulator.add(true, 0.0);
        accumulator.add(false, Double.NaN);

        PackStarFunctionalObservableResult result = accumulator.finish();

        assertEquals(PackStarFunctionalObservableResult.Status.UNRESOLVED_INVALID_SAMPLES,
                result.getStatus());
        assertEquals(1.0, result.getProbability(), 1.0e-15);
        assertTrue(Double.isNaN(result.getRestrictionFreeEnergy()));
        assertEquals(1, result.getInvalidSampleCount());
        assertFalse(result.isResolved());
    }

    @Test
    public void discreteAndComposedEventsInspectTheAssignment() {
        PackStarFunctionalEvent.Sample sample =
                new PackStarFunctionalEvent.Sample(
                        0, new int[]{2, 7, 4}, 0.0, 0.0, null);
        PackStarFunctionalEvent event = PackStarFunctionalEvents.allOf(
                PackStarFunctionalEvents.conformationAt(0, 2),
                PackStarFunctionalEvents.conformationAt(1, 7));

        assertTrue(event.test(sample));
        assertTrue(PackStarFunctionalEvents.conformationEquals(
                new int[]{2, 7, 4}).test(sample));
        assertFalse(PackStarFunctionalEvents.conformationAt(2, 9).test(sample));
    }
}
