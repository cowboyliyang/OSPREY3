package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

/**
 * PartialFixCache Integration Helper
 *
 * Provides easy integration of PartialFixCache (Phase 4 internally) into existing MARK* workflow.
 * Uses BWM*-inspired branch decomposition with L-set/M-set separation for fast caching.
 */
public class PartialFixIntegration {

    // Feature flag for PartialFixCache
    public static boolean ENABLE_PARTIALFIX_CACHE = true;

    // Global PartialFixCache instance (one per confSpace)
    private static java.util.Map<SimpleConfSpace, PartialFixCache> globalCaches =
            new java.util.concurrent.ConcurrentHashMap<>();

    /**
     * Get or create PartialFixCache for a conformation space
     */
    public static PartialFixCache getOrCreateCache(SimpleConfSpace confSpace) {
        if (!ENABLE_PARTIALFIX_CACHE) {
            return null;
        }

        return globalCaches.computeIfAbsent(confSpace, cs -> {
            // Create branch decomposition
            BranchDecomposition branchDecomp = new BranchDecomposition(cs);

            // Optimization: If tree is single-leaf, PartialFixCache provides no benefit
            // Return null to avoid overhead of cache creation and queries
            if (branchDecomp.root.isLeaf) {
                System.out.println("PartialFixCache SKIPPED for confSpace with " +
                                 cs.positions.size() + " positions (tree is single-leaf, no caching benefit)");
                return null;
            }

            // Create PartialFixCache
            PartialFixCache cache = new PartialFixCache(branchDecomp, cs);

            System.out.println("PartialFixCache initialized for confSpace with " +
                             cs.positions.size() + " positions");
            System.out.println("  Branch width: " + branchDecomp.branchWidth);

            return cache;
        });
    }

    /**
     * Print statistics for all PartialFixCache instances
     */
    public static void printAllStatistics() {
        if (globalCaches.isEmpty()) {
            System.out.println("No PartialFixCache instances active");
            return;
        }

        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIALFIX CACHE GLOBAL STATISTICS");
        System.out.println("=".repeat(100));
        System.out.println("Total cache instances: " + globalCaches.size());
        System.out.println();

        int instanceNum = 1;
        for (java.util.Map.Entry<SimpleConfSpace, PartialFixCache> entry : globalCaches.entrySet()) {
            System.out.println("Cache Instance #" + instanceNum + " (ConfSpace: " +
                entry.getKey().positions.size() + " positions)");
            entry.getValue().printStatistics();
            instanceNum++;
        }
    }

    /**
     * Print compact one-line statistics for all caches
     */
    public static void printCompactStatistics() {
        if (globalCaches.isEmpty()) {
            return;
        }

        System.out.println("\n[PartialFixCache Global Summary]");
        int totalQueries = 0;
        long totalHits = 0;
        long totalMisses = 0;

        for (PartialFixCache cache : globalCaches.values()) {
            totalQueries += cache.getTotalQueries();
            totalHits += cache.getCacheHits();
            totalMisses += cache.getCacheMisses();
            cache.printCompactStatistics();
        }

        double overallHitRate = (totalHits + totalMisses) > 0 ?
            100.0 * totalHits / (totalHits + totalMisses) : 0.0;

        System.out.println(String.format("[Overall] Queries: %d | Hits: %d (%.1f%%) | Misses: %d",
            totalQueries, totalHits, overallHitRate, totalMisses));
    }

    /**
     * Get aggregated statistics across all caches
     */
    public static AggregatedStatistics getAggregatedStatistics() {
        long totalQueries = 0;
        long totalHits = 0;
        long totalMisses = 0;
        int totalCacheSize = 0;

        for (PartialFixCache cache : globalCaches.values()) {
            totalQueries += cache.getTotalQueries();
            totalHits += cache.getCacheHits();
            totalMisses += cache.getCacheMisses();
            totalCacheSize += cache.getCacheSize();
        }

        return new AggregatedStatistics(totalQueries, totalHits, totalMisses, totalCacheSize);
    }

    /**
     * Clear all PartialFixCache instances
     */
    public static void clearAllCaches() {
        for (PartialFixCache cache : globalCaches.values()) {
            cache.clearCache();
        }
        globalCaches.clear();
    }

    /**
     * Reset statistics for all caches without clearing cache content
     */
    public static void resetAllStatistics() {
        for (PartialFixCache cache : globalCaches.values()) {
            cache.resetStatistics();
        }
    }

    /**
     * Aggregated statistics across all cache instances
     */
    public static class AggregatedStatistics {
        public final long totalQueries;
        public final long totalHits;
        public final long totalMisses;
        public final int totalCacheSize;

        public AggregatedStatistics(long queries, long hits, long misses, int cacheSize) {
            this.totalQueries = queries;
            this.totalHits = hits;
            this.totalMisses = misses;
            this.totalCacheSize = cacheSize;
        }

        public double getHitRate() {
            long totalAccess = totalHits + totalMisses;
            return totalAccess > 0 ? (double) totalHits / totalAccess : 0.0;
        }

        @Override
        public String toString() {
            return String.format("Queries: %d, Hits: %d (%.1f%%), Misses: %d, Cache size: %d",
                totalQueries, totalHits, getHitRate() * 100, totalMisses, totalCacheSize);
        }
    }
}
