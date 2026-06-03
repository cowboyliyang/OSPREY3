package edu.duke.cs.osprey.markstar.framework;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/**
 * Byte-limited DP cache accounting (design 1.1), the doc's remaining [TODO]:
 * "Cache skip/eviction byte-accounting test".
 *
 * Covers the pure decisions that drive the cache's memory safety:
 *   - estimateDPTableBytes: 2 * mStates * 8, with a saturating clamp instead of
 *     a silent long overflow (an overflow here would let a huge table look tiny
 *     and slip past the budget -> OOM, exactly the 3k3q failure mode);
 *   - dpCacheShouldSkip: skip if M-state count reaches skipIfMStates, or the
 *     table exceeds the per-table or total byte budget (boundary-exact).
 */
public class TestDPCacheBudget {

    private static final long MiB = 1024L * 1024L;
    private static final long GiB = 1024L * MiB;

    // recommended defaults from the design doc
    private static final long MAX_TABLE = 256 * MiB;
    private static final long MAX_TOTAL = 4 * GiB;
    private static final long SKIP_M    = 8_000_000L;

    // ---- estimateDPTableBytes --------------------------------------------

    @Test
    public void byteEstimateIsTwoArraysOfDoubles() {
        assertEquals(0L, BranchMARKStarBound.estimateDPTableBytes(0));
        assertEquals(16L, BranchMARKStarBound.estimateDPTableBytes(1));     // 2 * 1 * 8
        assertEquals(160L, BranchMARKStarBound.estimateDPTableBytes(10));
        assertEquals(2L * Integer.MAX_VALUE * 8L,
                BranchMARKStarBound.estimateDPTableBytes(Integer.MAX_VALUE));
    }

    @Test
    public void byteEstimateSaturatesInsteadOfOverflowing() {
        long threshold = Long.MAX_VALUE / 16L; // largest mStates that does not overflow 2*8*m
        assertEquals(2L * threshold * 8L, BranchMARKStarBound.estimateDPTableBytes(threshold),
                "at the threshold the exact product is still returned");
        assertEquals(Long.MAX_VALUE, BranchMARKStarBound.estimateDPTableBytes(threshold + 1L),
                "just past the threshold must clamp, not wrap negative");
        assertEquals(Long.MAX_VALUE, BranchMARKStarBound.estimateDPTableBytes(Long.MAX_VALUE));
        assertTrue(BranchMARKStarBound.estimateDPTableBytes(Long.MAX_VALUE) > 0,
                "clamped value must stay positive (overflow would go negative)");
    }

    // ---- dpCacheShouldSkip ------------------------------------------------

    private static boolean skip(long len) {
        return BranchMARKStarBound.dpCacheShouldSkip(
                len, BranchMARKStarBound.estimateDPTableBytes(len), SKIP_M, MAX_TABLE, MAX_TOTAL);
    }

    @Test
    public void keepsSmallTables() {
        assertFalse(skip(1));
        assertFalse(skip(1000));
        assertFalse(skip(1_000_000)); // 16 MB, well under 256 MiB
    }

    @Test
    public void skipIfMStatesBoundaryIsInclusive() {
        // skip uses >=, so exactly SKIP_M skips and SKIP_M-1 does not (bytes are tiny here:
        // SKIP_M states == 128 MB < 256 MiB, so only the M-state rule can trigger).
        assertFalse(BranchMARKStarBound.dpCacheShouldSkip(
                SKIP_M - 1, estBytes(SKIP_M - 1), SKIP_M, MAX_TABLE, MAX_TOTAL));
        assertTrue(BranchMARKStarBound.dpCacheShouldSkip(
                SKIP_M, estBytes(SKIP_M), SKIP_M, MAX_TABLE, MAX_TOTAL));
    }

    @Test
    public void perTableByteBudgetBoundaryIsStrict() {
        // tableBytes == maxTableBytes is allowed; +1 over is skipped. Disable the other
        // two rules so only the per-table budget can fire.
        long big = Long.MAX_VALUE; // generous skipIfMStates / maxTotal
        assertFalse(BranchMARKStarBound.dpCacheShouldSkip(0, MAX_TABLE, big, MAX_TABLE, big),
                "exactly at the per-table budget must be cacheable");
        assertTrue(BranchMARKStarBound.dpCacheShouldSkip(0, MAX_TABLE + 1, big, MAX_TABLE, big),
                "one byte over the per-table budget must skip");
    }

    @Test
    public void totalByteBudgetIsEnforced() {
        long big = Long.MAX_VALUE;
        // a single table larger than the whole-cache budget is skipped even if maxTable is huge
        assertTrue(BranchMARKStarBound.dpCacheShouldSkip(0, MAX_TOTAL + 1, big, big, MAX_TOTAL));
        assertFalse(BranchMARKStarBound.dpCacheShouldSkip(0, MAX_TOTAL, big, big, MAX_TOTAL));
    }

    @Test
    public void hugeTableAlwaysSkips() {
        // a table that overflowed naive byte math must still be caught (saturated to MAX_VALUE)
        assertTrue(skip(Long.MAX_VALUE / 4L));
    }

    private static long estBytes(long m) {
        return BranchMARKStarBound.estimateDPTableBytes(m);
    }
}
