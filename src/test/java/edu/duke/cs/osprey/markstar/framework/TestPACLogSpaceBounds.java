package edu.duke.cs.osprey.markstar.framework;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/**
 * Log-space PAC epsilon (design 1.3), the doc's remaining [TODO]:
 * "PAC log-space epsilon tests with logZ=+/-1e6".
 *
 * PACPartitionFunction.epsilonFromLogBounds computes
 *   epsilon = 1 - exp(logZLower - logZUpper)
 * from the LOG ratio so it stays finite where the old direct-double path
 * (exp(logZ) then a BigDecimal ratio) collapsed to Z=0 and produced a fake
 * epsilon = NaN/1 (e.g. 5dc4). These tests pin that: epsilon depends only on the
 * gap between the bounds and is finite/in-range even at |logZ| ~ 1e6.
 */
public class TestPACLogSpaceBounds {

    private static double eps(double lo, double up) {
        return PACPartitionFunction.epsilonFromLogBounds(lo, up);
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
}
