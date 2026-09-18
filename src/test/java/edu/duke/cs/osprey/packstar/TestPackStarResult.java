package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class TestPackStarResult {
    @Test
    public void kstarRetainsFunctionalEstimateAndBackendCounters() {
        var accumulator = new PackStarFunctionalObservableAccumulator("contact", 0.6);
        accumulator.add(true, Math.log(3));
        accumulator.add(false, 0);
        var original = PartitionFunction.Result.makeAborted().setStat("fullCCD", 7);
        var result = new PackStarResult(original, accumulator.finish());
        var score = new KStarScore(result, result, result);
        var retained = (PackStarResult) score.complex;
        assertEquals(7, retained.getStat("fullCCD"));
        assertEquals(original.status, retained.status);
        assertEquals(0.75, retained.getFunctionalObservableResult().getProbability(), 1e-12);
        assertTrue(retained.getFunctionalObservableResult().isResolved());
    }

    @Test
    public void formalProtocolRejectsRuntimeDrift() {
        PackStar4wyuFormalCase.validateRuntimeSettings(0.1, 32);
        assertThrows(IllegalArgumentException.class,
                () -> PackStar4wyuFormalCase.validateRuntimeSettings(0.2, 32));
    }
}
