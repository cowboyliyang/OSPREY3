package edu.duke.cs.osprey.branchdp;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

/**
 * Correctness tests for the sharded final-DP table abstraction.
 *
 * This covers the lynchpin of the Phase B sharding work in
 * BranchMARKStar_PAC_Reliability_Design.md: dense-vs-sharded value equivalence
 * and the long shard-routing math, which must hold exactly so the sharded path
 * cannot silently return wrong DP values on large cases.
 *
 * These are pure unit tests over the DPTable implementations: no conf space /
 * protein setup required. (The edge-level mixed-radix index math in
 * RootedTreeEdge.computeIndexInA is an instance method tied to a constructed
 * edge; it is left to an integration-level test.)
 */
public class TestDPTable {

	// deterministic, distinct lower/upper values per index
	private static double lowerVal(long idx) { return -1.5 * (idx + 1); }
	private static double upperVal(long idx) { return  3.25 * (idx + 1) + 0.5; }

	// ---- dense vs sharded equivalence -------------------------------------

	private void checkDenseShardedEquivalent(int size, int shardSize) {
		DPTable dense = new DenseDPTable(size);
		DPTable sharded = new ShardedDPTable(size, shardSize);

		assertEquals(size, dense.size());
		assertEquals(size, sharded.size());
		assertTrue(dense.isDenseArrayBacked());
		assertFalse(sharded.isDenseArrayBacked());

		// write through both, in several strides to stress the routing math
		for (int step : new int[] { 1, 3, 7 }) {
			for (long idx = 0; idx < size; idx += step) {
				dense.set(idx, lowerVal(idx), upperVal(idx));
				sharded.set(idx, lowerVal(idx), upperVal(idx));
			}
		}

		for (long idx = 0; idx < size; idx++) {
			String where = " (size=" + size + ", shardSize=" + shardSize + ", idx=" + idx + ")";
			assertEquals(lowerVal(idx), dense.lower(idx), 0.0, "dense lower" + where);
			assertEquals(upperVal(idx), dense.upper(idx), 0.0, "dense upper" + where);
			assertEquals(dense.lower(idx), sharded.lower(idx), 0.0, "sharded vs dense lower" + where);
			assertEquals(dense.upper(idx), sharded.upper(idx), 0.0, "sharded vs dense upper" + where);
		}
	}

	@Test
	public void denseShardedEquivalence_variousShapes() {
		// exact multiples, non-multiples, size==1, and shardSize >= size (single shard)
		int[][] cases = {
			{1, 1}, {1, 8}, {7, 1}, {8, 4}, {8, 3}, {100, 16},
			{1000, 7}, {1003, 7}, {1024, 1024}, {1025, 1024}, {5000, 97},
		};
		for (int[] c : cases) {
			checkDenseShardedEquivalent(c[0], c[1]);
		}
	}

	@Test
	public void shardBoundaryRouting() {
		// 13 states with shardSize 4 -> shards of (4,4,4,1); check straddling indices
		int shardSize = 4;
		int size = 13;
		ShardedDPTable t = new ShardedDPTable(size, shardSize);
		for (long idx = 0; idx < size; idx++) {
			t.set(idx, lowerVal(idx), upperVal(idx));
		}
		for (long idx : new long[] {0, 3, 4, 7, 8, 11, 12}) {
			assertEquals(lowerVal(idx), t.lower(idx), 0.0, "lower at boundary " + idx);
			assertEquals(upperVal(idx), t.upper(idx), 0.0, "upper at boundary " + idx);
		}
	}

	@Test
	public void manyShardsRoutingIsContiguous() {
		// shardSize 1 forces one shard per index: every routing decision exercised
		int size = 257;
		ShardedDPTable t = new ShardedDPTable(size, 1);
		for (long idx = 0; idx < size; idx++) {
			t.set(idx, lowerVal(idx), upperVal(idx));
		}
		for (long idx = 0; idx < size; idx++) {
			assertEquals(lowerVal(idx), t.lower(idx), 0.0, "lower at idx " + idx);
			assertEquals(upperVal(idx), t.upper(idx), 0.0, "upper at idx " + idx);
		}
	}

	@Test
	public void fillSetsAllEntries_bothBackings() {
		DPTable dense = new DenseDPTable(50);
		DPTable sharded = new ShardedDPTable(50, 7);
		dense.fill(-9.0, 9.0);
		sharded.fill(-9.0, 9.0);
		for (long idx = 0; idx < 50; idx++) {
			assertEquals(-9.0, dense.lower(idx), 0.0);
			assertEquals(9.0, dense.upper(idx), 0.0);
			assertEquals(-9.0, sharded.lower(idx), 0.0);
			assertEquals(9.0, sharded.upper(idx), 0.0);
		}
	}

	@Test
	public void outOfBoundsRejected_sharded() {
		ShardedDPTable t = new ShardedDPTable(10, 4);
		assertThrows(IndexOutOfBoundsException.class, () -> t.lower(10));
		assertThrows(IndexOutOfBoundsException.class, () -> t.lower(-1));
		assertThrows(IndexOutOfBoundsException.class, () -> t.set(10, 0, 0));
	}

	@Test
	public void estimatedBytesAndFlags() {
		assertEquals(2L * 10 * 8L, new DenseDPTable(10).estimatedBytes());
		assertEquals(2L * 10 * 8L, new ShardedDPTable(10, 4).estimatedBytes());
	}

	@Test
	public void denseRejectsOverIntLimit() {
		assertThrows(IllegalArgumentException.class,
			() -> new DenseDPTable((long) Integer.MAX_VALUE + 1L));
	}

	@Test
	public void shardedUnsafeArrayAccessThrows() {
		ShardedDPTable t = new ShardedDPTable(10, 4);
		assertThrows(IllegalStateException.class, t::lowerArrayUnsafe);
		assertThrows(IllegalStateException.class, t::upperArrayUnsafe);
	}

	@Test
	public void denseUnsafeArrayMatchesAccessors() {
		DenseDPTable t = new DenseDPTable(8);
		for (long idx = 0; idx < 8; idx++) {
			t.set(idx, lowerVal(idx), upperVal(idx));
		}
		double[] lo = t.lowerArrayUnsafe();
		double[] hi = t.upperArrayUnsafe();
		for (int i = 0; i < 8; i++) {
			assertEquals(t.lower(i), lo[i], 0.0);
			assertEquals(t.upper(i), hi[i], 0.0);
		}
	}
}
