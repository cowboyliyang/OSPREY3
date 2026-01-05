package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.minimization.Minimizer;

import java.util.*;

/**
 * Simplified Phase 2: Minimization Result Cache
 *
 * This is a lightweight wrapper that caches minimization results for identical conformations.
 * Unlike the full SubtreeDOFCache with branch decomposition, this simply caches complete
 * conformation minimizations.
 *
 * Benefits:
 * - Simple to integrate (just wraps existing minimizer calls)
 * - No modifications to minimization code needed
 * - Still provides speedup when same conformations are minimized multiple times
 *
 * Limitations:
 * - Doesn't exploit subtree structure (requires identical full conformations)
 * - Less powerful than full Phase 2 with branch decomposition
 *
 * This serves as a stepping stone toward full Phase 2 integration.
 */
public class MinimizationCache {

    // Cache key: RCTuple -> Minimization result
    private final Map<RCTuple, CachedResult> cache;

    // Statistics
    private long cacheHits = 0;
    private long cacheMisses = 0;
    private long totalQueries = 0;

    // Configuration
    private static final int MAX_CACHE_SIZE = 50000;
    public static boolean ENABLE_CACHE = false; // Global enable/disable

    public MinimizationCache() {
        this.cache = new LinkedHashMap<RCTuple, CachedResult>(16, 0.75f, true) {
            @Override
            protected boolean removeEldestEntry(Map.Entry<RCTuple, CachedResult> eldest) {
                return size() > MAX_CACHE_SIZE;
            }
        };
    }

    /**
     * Minimization wrapper with caching
     *
     * @param conf The conformation being minimized
     * @param minimizer The minimizer to use (if not cached)
     * @param initialDOFs Initial DOF values
     * @return Minimization result (from cache or freshly computed)
     */
    public Minimizer.Result minimizeWithCache(
            RCTuple conf,
            Minimizer minimizer,
            DoubleMatrix1D initialDOFs) {

        if (!ENABLE_CACHE) {
            // Cache disabled, just minimize
            return minimizer.minimizeFrom(initialDOFs);
        }

        totalQueries++;

        // Check cache
        CachedResult cached = cache.get(conf);
        if (cached != null) {
            cacheHits++;
            return new Minimizer.Result(cached.dofValues.copy(), cached.energy);
        }

        // Cache miss - perform minimization
        cacheMisses++;
        Minimizer.Result result = minimizer.minimizeFrom(initialDOFs);

        // Store in cache (use conf as-is, assuming it's normalized)
        cache.put(conf, new CachedResult(result.dofValues, result.energy));

        return result;
    }

    /**
     * Alternative wrapper: just check cache, don't minimize if not found
     * Useful when you want to control minimization timing
     */
    public CachedResult getCached(RCTuple conf) {
        if (!ENABLE_CACHE) {
            return null;
        }

        totalQueries++;
        CachedResult cached = cache.get(conf);

        if (cached != null) {
            cacheHits++;
        } else {
            cacheMisses++;
        }

        return cached;
    }

    /**
     * Manually add a result to the cache
     */
    public void put(RCTuple conf, DoubleMatrix1D dofValues, double energy) {
        if (ENABLE_CACHE) {
            cache.put(conf, new CachedResult(dofValues, energy));
        }
    }

    // Statistics

    public void printStats() {
        if (!ENABLE_CACHE) {
            System.out.println("\n=== Minimization Cache: DISABLED ===\n");
            return;
        }

        System.out.println("\n=== Minimization Cache Statistics ===");
        System.out.println("Total queries: " + totalQueries);
        System.out.println("Cache hits:    " + cacheHits);
        System.out.println("Cache misses:  " + cacheMisses);

        if (totalQueries > 0) {
            double hitRate = 100.0 * cacheHits / totalQueries;
            System.out.println("Hit rate:      " + String.format("%.1f%%", hitRate));
        }

        System.out.println("Cache size:    " + cache.size() + " / " + MAX_CACHE_SIZE);

        // Estimate memory usage
        long estimatedBytes = cache.size() * 500; // Rough estimate: 500 bytes per entry
        System.out.println("Est. memory:   " + String.format("%.1f MB", estimatedBytes / 1024.0 / 1024.0));

        System.out.println("======================================\n");
    }

    public void clearCache() {
        cache.clear();
        cacheHits = 0;
        cacheMisses = 0;
        totalQueries = 0;
    }

    public long getCacheHits() {
        return cacheHits;
    }

    public long getCacheMisses() {
        return cacheMisses;
    }

    public double getHitRate() {
        if (totalQueries == 0) return 0.0;
        return 100.0 * cacheHits / totalQueries;
    }

    /**
     * Cached minimization result
     */
    public static class CachedResult {
        public final DoubleMatrix1D dofValues;
        public final double energy;

        public CachedResult(DoubleMatrix1D dofValues, double energy) {
            this.dofValues = dofValues.copy(); // Deep copy to avoid external modification
            this.energy = energy;
        }
    }
}
