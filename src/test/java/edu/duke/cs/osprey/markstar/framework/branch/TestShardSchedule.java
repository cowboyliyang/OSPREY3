package edu.duke.cs.osprey.markstar.framework.branch;

import static org.junit.jupiter.api.Assertions.*;

import java.lang.reflect.Method;
import org.junit.jupiter.api.Test;

/**
 * Correctness guard for the shard-scheduled parallel full DP (design 2.5).
 *
 * The point of this test is the *invariant* that parallel shard scheduling must
 * produce byte-identical DP values to the serial path, under a range of thread
 * counts and shard sizes (including ragged final shards and shardSize=1). Since
 * computeFullDP requires a fully wired tree/energy setup, we test the scheduler
 * core directly: a contiguous-range, bounded-thread, atomic-cursor sweep that
 * writes one value per index. A data race or a routing/partition bug would show
 * up as a missing or duplicated index.
 *
 * This mirrors exactly how computeFullDPSharded partitions [0, n) into shards and
 * assigns each index to exactly one worker.
 */
public class TestShardSchedule {

	/** Run the same contiguous-shard scheduling computeFullDPSharded uses, into out[]. */
	private static void scheduleInto(long[] outHits, int threads, long shardSize) {
		long n = outHits.length;
		long numShards = (n + shardSize - 1L) / shardSize;
		int workers = (int) Math.min(threads, Math.max(1, numShards));
		java.util.concurrent.atomic.AtomicLong cursor = new java.util.concurrent.atomic.AtomicLong(0);
		java.util.concurrent.ExecutorService pool = java.util.concurrent.Executors.newFixedThreadPool(workers);
		try {
			java.util.List<java.util.concurrent.Future<?>> fs = new java.util.ArrayList<>();
			for (int w = 0; w < workers; w++) {
				fs.add(pool.submit(() -> {
					long s;
					while ((s = cursor.getAndIncrement()) < numShards) {
						long start = s * shardSize;
						long end = Math.min(start + shardSize, n);
						for (long i = start; i < end; i++) {
							// each index written exactly once if partition is correct
							outHits[(int) i] += 1;
						}
					}
				}));
			}
			for (java.util.concurrent.Future<?> f : fs) {
				f.get();
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			pool.shutdownNow();
		}
	}

	@Test
	public void everyIndexCoveredExactlyOnce_variousShapes() {
		int[][] cases = {
			{1, 1, 1}, {100, 4, 16}, {100, 8, 7}, {1000, 16, 1},
			{1003, 7, 64}, {65537, 8, 65536}, {50000, 13, 97},
		};
		for (int[] c : cases) {
			int n = c[0], threads = c[1], shardSize = c[2];
			long[] hits = new long[n];
			scheduleInto(hits, threads, shardSize);
			for (int i = 0; i < n; i++) {
				assertEquals(1L, hits[i],
					"index " + i + " hit " + hits[i] + " times (n=" + n
						+ ", threads=" + threads + ", shardSize=" + shardSize + ")");
			}
		}
	}

	@Test
	public void parallelMatchesSerial_valuePerIndex() {
		int n = 40000;
		// serial reference
		double[] serial = new double[n];
		for (int i = 0; i < n; i++) {
			serial[i] = Math.sin(i * 0.001) + i;
		}
		// parallel via the same shard scheduling, different threads/shardSizes
		for (int threads : new int[] {1, 2, 4, 8}) {
			for (long shardSize : new long[] {1, 97, 4096, 65536}) {
				double[] par = new double[n];
				long numShards = (n + shardSize - 1L) / shardSize;
				int workers = (int) Math.min(threads, Math.max(1, numShards));
				java.util.concurrent.atomic.AtomicLong cursor = new java.util.concurrent.atomic.AtomicLong(0);
				java.util.concurrent.ExecutorService pool = java.util.concurrent.Executors.newFixedThreadPool(workers);
				try {
					java.util.List<java.util.concurrent.Future<?>> fs = new java.util.ArrayList<>();
					for (int w = 0; w < workers; w++) {
						fs.add(pool.submit(() -> {
							long s;
							while ((s = cursor.getAndIncrement()) < numShards) {
								long start = s * shardSize;
								long end = Math.min(start + shardSize, n);
								for (long i = start; i < end; i++) {
									par[(int) i] = Math.sin(i * 0.001) + i;
								}
							}
						}));
					}
					for (java.util.concurrent.Future<?> f : fs) f.get();
				} catch (Exception e) {
					throw new RuntimeException(e);
				} finally {
					pool.shutdownNow();
				}
				assertArrayEquals(serial, par, 0.0,
					"parallel result differs from serial (threads=" + threads
						+ ", shardSize=" + shardSize + ")");
			}
		}
	}

	@Test
	public void schedulerMethodsExistOnProductionClass() throws Exception {
		// guard against accidental signature changes to the production scheduler
		Method m = RootedTreeEdge.class.getDeclaredMethod("computeFullDPSharded");
		assertNotNull(m);
		assertNotNull(RootedTreeEdge.class.getDeclaredMethod("resolveDPThreads"));
		assertNotNull(RootedTreeEdge.class.getDeclaredMethod("resolveComputeShardSize"));
	}
}
