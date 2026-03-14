package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Phase 5: Partial Start Cache Integration
 *
 * Provides warm start DOF values for CCD minimization by caching
 * minimized conformations and reusing DOF values from similar tuples.
 *
 * Key Features:
 * - Stores 7-tuple results as 5-subset and 6-subset entries
 * - Query: O(28) hash lookups instead of O(N) linear scan
 * - 6-subset match = 6 positions overlap (better warm start)
 * - 5-subset match = 5 positions overlap (good warm start)
 * - Non-overlapping positions use default DOF values
 */
public class PartialStartIntegration {

    // Global feature flag
    public static boolean ENABLE_PARTIALSTART_CACHE = false;  // Disabled by default

    // Global cache instances (one per confSpace)
    private static final Map<SimpleConfSpace, PartialStartCache> globalCaches =
            new ConcurrentHashMap<>();

    /**
     * Get or create PartialStartCache for a conformation space
     */
    public static PartialStartCache getOrCreateCache(SimpleConfSpace confSpace, ParametricMolecule pmol) {
        if (!ENABLE_PARTIALSTART_CACHE) {
            return null;
        }

        return globalCaches.computeIfAbsent(confSpace, cs -> {
            PartialStartCache cache = new PartialStartCache();
            cache.initialize(pmol);
            return cache;
        });
    }

    /**
     * Get warm start DOF values for a tuple
     * Returns null if cache is disabled or no match found
     */
    public static double[] getWarmStart(SimpleConfSpace confSpace, ParametricMolecule pmol,
                                        RCTuple tuple, DoubleMatrix1D defaultDOFs) {
        if (!ENABLE_PARTIALSTART_CACHE) {
            return null;
        }

        PartialStartCache cache = getOrCreateCache(confSpace, pmol);
        if (cache == null) {
            return null;
        }

        return cache.getWarmStart(tuple, defaultDOFs);
    }

    /**
     * Get warm start information including unmatched DOF indices
     * Returns null if cache is disabled or no match found
     */
    public static PartialStartCache.WarmStartResult getWarmStartWithInfo(
            SimpleConfSpace confSpace,
            ParametricMolecule pmol,
            RCTuple tuple,
            DoubleMatrix1D defaultDOFs) {

        if (!ENABLE_PARTIALSTART_CACHE) {
            return null;
        }

        PartialStartCache cache = getOrCreateCache(confSpace, pmol);
        if (cache == null) {
            return null;
        }

        return cache.getWarmStartWithInfo(tuple, defaultDOFs);
    }

    /**
     * Cache minimization result
     */
    public static void cacheResult(SimpleConfSpace confSpace, ParametricMolecule pmol,
                                   RCTuple tuple, DoubleMatrix1D dofValues) {
        if (!ENABLE_PARTIALSTART_CACHE) {
            return;
        }

        PartialStartCache cache = getOrCreateCache(confSpace, pmol);
        if (cache != null) {
            cache.put(tuple, dofValues);
        }
    }

    /**
     * Print statistics for all cache instances
     */
    public static void printGlobalStats() {
        if (globalCaches.isEmpty()) {
            System.out.println("No PartialStartCache instances active");
            return;
        }

        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIAL START CACHE (PHASE 5) GLOBAL STATISTICS");
        System.out.println("=".repeat(100));
        System.out.println("Total cache instances: " + globalCaches.size());
        System.out.println();

        long total6Hits = 0;
        long total5Hits = 0;
        long totalMisses = 0;
        int totalUniqueEntries = 0;
        int totalIndexEntries = 0;

        int instanceNum = 1;
        for (Map.Entry<SimpleConfSpace, PartialStartCache> entry : globalCaches.entrySet()) {
            PartialStartCache cache = entry.getValue();
            PartialStartCache.CacheStats stats = cache.getStats();

            System.out.println("Cache Instance #" + instanceNum + " (ConfSpace: " +
                entry.getKey().positions.size() + " positions)");
            System.out.println("  " + stats);
            System.out.println("  Unique entries: " + cache.size());
            System.out.println("  Total index entries: " + cache.totalIndexEntries());
            System.out.println();

            total6Hits += stats.sixSubsetHits;
            total5Hits += stats.fiveSubsetHits;
            totalMisses += stats.misses;
            totalUniqueEntries += cache.size();
            totalIndexEntries += cache.totalIndexEntries();

            instanceNum++;
        }

        // Summary
        long totalQueries = total6Hits + total5Hits + totalMisses;
        long totalHits = total6Hits + total5Hits;
        double hitRate = totalQueries > 0 ? 100.0 * totalHits / totalQueries : 0.0;

        System.out.println("=".repeat(100));
        System.out.println("SUMMARY:");
        System.out.println(String.format("  Total queries: %d", totalQueries));
        System.out.println(String.format("  6-subset hits: %d (%.1f%%)", total6Hits,
            totalQueries > 0 ? 100.0 * total6Hits / totalQueries : 0.0));
        System.out.println(String.format("  5-subset hits: %d (%.1f%%)", total5Hits,
            totalQueries > 0 ? 100.0 * total5Hits / totalQueries : 0.0));
        System.out.println(String.format("  Misses: %d (%.1f%%)", totalMisses,
            totalQueries > 0 ? 100.0 * totalMisses / totalQueries : 0.0));
        System.out.println(String.format("  Overall hit rate: %.1f%%", hitRate));
        System.out.println(String.format("  Unique entries: %d", totalUniqueEntries));
        System.out.println(String.format("  Total index entries: %d", totalIndexEntries));
        System.out.println("=".repeat(100) + "\n");
    }

    /**
     * Clear all cache instances
     */
    public static void clearAllCaches() {
        for (PartialStartCache cache : globalCaches.values()) {
            cache.clear();
        }
        globalCaches.clear();
    }

    /**
     * Get aggregated statistics
     */
    public static AggregatedStats getAggregatedStats() {
        long total6Hits = 0;
        long total5Hits = 0;
        long totalMisses = 0;
        int totalSize = 0;

        for (PartialStartCache cache : globalCaches.values()) {
            PartialStartCache.CacheStats stats = cache.getStats();
            total6Hits += stats.sixSubsetHits;
            total5Hits += stats.fiveSubsetHits;
            totalMisses += stats.misses;
            totalSize += cache.size();
        }

        return new AggregatedStats(total6Hits, total5Hits, totalMisses, totalSize);
    }

    /**
     * Aggregated statistics across all cache instances
     */
    public static class AggregatedStats {
        public final long sixSubsetHits;
        public final long fiveSubsetHits;
        public final long misses;
        public final int cacheSize;

        public AggregatedStats(long sixHits, long fiveHits, long misses, int size) {
            this.sixSubsetHits = sixHits;
            this.fiveSubsetHits = fiveHits;
            this.misses = misses;
            this.cacheSize = size;
        }

        public long totalQueries() {
            return sixSubsetHits + fiveSubsetHits + misses;
        }

        public long totalHits() {
            return sixSubsetHits + fiveSubsetHits;
        }

        public double hitRate() {
            long total = totalQueries();
            return total > 0 ? (double) totalHits() / total : 0.0;
        }

        @Override
        public String toString() {
            return String.format("Queries: %d, 6-hits: %d, 5-hits: %d, Misses: %d, HitRate: %.1f%%, Size: %d",
                totalQueries(), sixSubsetHits, fiveSubsetHits, misses, hitRate() * 100, cacheSize);
        }
    }

    // ========== Legacy instance-based API (for backward compatibility) ==========

    private final PartialStartCache cache;
    private final ObjectiveFunction objFunc;
    private final ParametricMolecule pmol;
    private boolean enableCache;

    public PartialStartIntegration(ObjectiveFunction objFunc, ParametricMolecule pmol) {
        this.cache = new PartialStartCache();
        this.objFunc = objFunc;
        this.pmol = pmol;
        this.enableCache = true;
        cache.initialize(pmol);
    }

    /**
     * Minimize with automatic warm start from cache
     */
    public Minimizer.Result minimize(RCTuple tuple) {
        Minimizer minimizer = new SimpleCCDMinimizer(objFunc);
        Minimizer.Result result;

        if (enableCache) {
            DoubleMatrix1D centerDOFs = objFunc.getDOFsCenter();
            double[] warmStart = cache.getWarmStart(tuple, centerDOFs);

            if (warmStart != null) {
                result = minimizer.minimizeFrom(convertToMatrix(warmStart));
            } else {
                result = minimizer.minimizeFromCenter();
            }

            cache.put(tuple, result.dofValues);
        } else {
            result = minimizer.minimizeFromCenter();
        }

        return result;
    }

    public void setEnableCache(boolean enable) {
        this.enableCache = enable;
    }

    public PartialStartCache.CacheStats getStats() {
        return cache.getStats();
    }

    public void clearCache() {
        cache.clear();
    }

    public int getCacheSize() {
        return cache.size();
    }

    private static DoubleMatrix1D convertToMatrix(double[] arr) {
        DoubleMatrix1D matrix = cern.colt.matrix.DoubleFactory1D.dense.make(arr.length);
        for (int i = 0; i < arr.length; i++) {
            matrix.set(i, arr[i]);
        }
        return matrix;
    }
}
