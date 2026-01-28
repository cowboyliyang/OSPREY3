package edu.duke.cs.osprey.markstar;

/**
 * Global timer to track minimization time across all phases of MARK*.
 *
 * This class provides centralized timing statistics for:
 * - Phase 1: Triple correction computation
 * - Phase 2: A* search minimization (through SubtreeDOFCache)
 *
 * All times are tracked in nanoseconds for precision.
 */
public class MinimizationTimer {

    // Phase 1: Triple Correction Minimization
    private static long phase1MinimizationTimeNs = 0;
    private static long phase1MinimizationCount = 0;

    // Phase 2: A* Search Minimization (through cache)
    private static long phase2MinimizationTimeNs = 0;
    private static long phase2MinimizationCount = 0;

    // Cache lookup overhead
    private static long cacheLookupTimeNs = 0;

    // Total program runtime for percentage calculation
    private static long programStartTimeNs = 0;
    private static long programEndTimeNs = 0;

    /**
     * Reset all timing statistics. Call this at the beginning of each test run.
     */
    public static synchronized void reset() {
        phase1MinimizationTimeNs = 0;
        phase1MinimizationCount = 0;
        phase2MinimizationTimeNs = 0;
        phase2MinimizationCount = 0;
        cacheLookupTimeNs = 0;
        programStartTimeNs = 0;
        programEndTimeNs = 0;
    }

    /**
     * Start the global timer (call at beginning of MARK* run)
     */
    public static synchronized void startProgram() {
        programStartTimeNs = System.nanoTime();
    }

    /**
     * Stop the global timer (call at end of MARK* run)
     */
    public static synchronized void stopProgram() {
        programEndTimeNs = System.nanoTime();
    }

    // Phase 1 recording methods

    public static synchronized void recordPhase1Minimization(long durationNs) {
        phase1MinimizationTimeNs += durationNs;
        phase1MinimizationCount++;
    }

    // Phase 2 recording methods

    public static synchronized void recordPhase2Minimization(long durationNs) {
        phase2MinimizationTimeNs += durationNs;
        phase2MinimizationCount++;
    }

    public static synchronized void recordCacheLookup(long durationNs) {
        cacheLookupTimeNs += durationNs;
    }

    // Getter methods

    public static synchronized long getPhase1MinimizationTimeNs() {
        return phase1MinimizationTimeNs;
    }

    public static synchronized long getPhase1MinimizationCount() {
        return phase1MinimizationCount;
    }

    public static synchronized long getPhase2MinimizationTimeNs() {
        return phase2MinimizationTimeNs;
    }

    public static synchronized long getPhase2MinimizationCount() {
        return phase2MinimizationCount;
    }

    public static synchronized long getCacheLookupTimeNs() {
        return cacheLookupTimeNs;
    }

    public static synchronized long getTotalProgramTimeNs() {
        if (programEndTimeNs > 0) {
            return programEndTimeNs - programStartTimeNs;
        } else {
            return System.nanoTime() - programStartTimeNs;
        }
    }

    /**
     * Print comprehensive timing statistics
     */
    public static synchronized void printStatistics() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("MINIMIZATION TIMING STATISTICS (All Phases)");
        System.out.println("=".repeat(100));

        long totalMinTimeNs = phase1MinimizationTimeNs + phase2MinimizationTimeNs;
        long totalMinCount = phase1MinimizationCount + phase2MinimizationCount;
        long totalProgramTimeNs = getTotalProgramTimeNs();

        // Convert to seconds
        double phase1MinSec = phase1MinimizationTimeNs / 1e9;
        double phase2MinSec = phase2MinimizationTimeNs / 1e9;
        double totalMinSec = totalMinTimeNs / 1e9;
        double cacheLookupSec = cacheLookupTimeNs / 1e9;
        double totalProgramSec = totalProgramTimeNs / 1e9;

        System.out.println("\n--- Phase 1: Triple Correction Minimization ---");
        System.out.println("Minimization calls:     " + phase1MinimizationCount);
        System.out.println("Total time:             " + String.format("%.2f s", phase1MinSec));
        if (phase1MinimizationCount > 0) {
            double avgMs = phase1MinimizationTimeNs / (phase1MinimizationCount * 1e6);
            System.out.println("Average per call:       " + String.format("%.2f ms", avgMs));
        }
        if (totalProgramTimeNs > 0) {
            double pct = 100.0 * phase1MinimizationTimeNs / totalProgramTimeNs;
            System.out.println("% of total runtime:     " + String.format("%.1f%%", pct));
        }

        System.out.println("\n--- Phase 2: A* Search Minimization (via Cache) ---");
        System.out.println("Minimization calls:     " + phase2MinimizationCount);
        System.out.println("Total time:             " + String.format("%.2f s", phase2MinSec));
        if (phase2MinimizationCount > 0) {
            double avgMs = phase2MinimizationTimeNs / (phase2MinimizationCount * 1e6);
            System.out.println("Average per call:       " + String.format("%.2f ms", avgMs));
        }
        if (totalProgramTimeNs > 0) {
            double pct = 100.0 * phase2MinimizationTimeNs / totalProgramTimeNs;
            System.out.println("% of total runtime:     " + String.format("%.1f%%", pct));
        }

        System.out.println("\n--- Cache Overhead ---");
        System.out.println("Cache lookup time:      " + String.format("%.2f s", cacheLookupSec));
        if (totalProgramTimeNs > 0) {
            double pct = 100.0 * cacheLookupTimeNs / totalProgramTimeNs;
            System.out.println("% of total runtime:     " + String.format("%.1f%%", pct));
        }

        System.out.println("\n--- TOTAL MINIMIZATION ---");
        System.out.println("Total minimizations:    " + totalMinCount);
        System.out.println("Total min time:         " + String.format("%.2f s", totalMinSec));
        if (totalMinCount > 0) {
            double avgMs = totalMinTimeNs / (totalMinCount * 1e6);
            System.out.println("Average per call:       " + String.format("%.2f ms", avgMs));
        }

        if (totalProgramTimeNs > 0) {
            System.out.println("\n--- OVERALL BREAKDOWN ---");
            System.out.println("Total program time:     " + String.format("%.2f s", totalProgramSec));

            double minPct = 100.0 * totalMinTimeNs / totalProgramTimeNs;
            double cachePct = 100.0 * cacheLookupTimeNs / totalProgramTimeNs;
            double otherPct = 100.0 - minPct - cachePct;

            System.out.println("Time in minimization:   " + String.format("%.2f s (%.1f%%)", totalMinSec, minPct));
            System.out.println("Time in cache lookup:   " + String.format("%.2f s (%.1f%%)", cacheLookupSec, cachePct));
            System.out.println("Time in other ops:      " + String.format("%.2f s (%.1f%%)",
                (totalProgramTimeNs - totalMinTimeNs - cacheLookupTimeNs) / 1e9, otherPct));

            System.out.println("\nKEY INSIGHT:");
            if (minPct > 50) {
                System.out.println("- Minimization dominates runtime (" + String.format("%.1f%%", minPct) +
                    ") - optimization focus should be on minimization efficiency");
            } else if (minPct > 25) {
                System.out.println("- Minimization is significant (" + String.format("%.1f%%", minPct) +
                    ") but not dominant - balanced optimization needed");
            } else {
                System.out.println("- Minimization is only " + String.format("%.1f%%", minPct) +
                    " of runtime - other operations dominate (energy calculation, search, etc.)");
            }
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }
}
