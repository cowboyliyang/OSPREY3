package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.structure.Residue;

import java.util.*;

/**
 * Phase 5: Partial Start Cache
 *
 * Provides warm start DOF values for CCD minimization by caching
 * minimized conformations and reusing DOF values from similar tuples.
 *
 * Key Idea:
 * - Store 7-tuple minimization results as 5-subset and 6-subset entries
 * - When querying, find matching subsets and use cached DOF values
 * - For overlapping positions: use cached DOF values (warm start)
 * - For non-overlapping positions: use default DOF values
 *
 * Index Structure:
 * - 6-subset indices: C(7,6) = 7 combinations (higher priority, 6 positions overlap)
 * - 5-subset indices: C(7,5) = 21 combinations (lower priority, 5 positions overlap)
 *
 * Performance:
 * - Lookup: O(7 + 21) = O(28) hash lookups instead of O(N) linear scan
 * - Storage: Decreases over time as subsets fill up
 */
public class PartialStartCache {

    /**
     * Result of warm start query with unmatched DOF information
     */
    public static class WarmStartResult {
        public final double[] dofValues;
        public final int[] unmatchedDOFIndices;  // Indices of DOFs for unmatched positions
        public final int subsetSize;  // 5 or 6

        public WarmStartResult(double[] dofValues, int[] unmatchedDOFIndices, int subsetSize) {
            this.dofValues = dofValues;
            this.unmatchedDOFIndices = unmatchedDOFIndices;
            this.subsetSize = subsetSize;
        }
    }

    // 6-subset indices: C(7,6) = 7 combinations (query first, more overlap)
    private final List<Map<String, CachedEntry>> sixSubsetIndices;
    private final List<int[]> sixSubsetCombinations;

    // 5-subset indices: C(7,5) = 21 combinations (query second)
    private final List<Map<String, CachedEntry>> fiveSubsetIndices;
    private final List<int[]> fiveSubsetCombinations;

    // Configuration
    private static final int TUPLE_SIZE = 7;
    private static final int NUM_6_SUBSETS = 7;   // C(7,6)
    private static final int NUM_5_SUBSETS = 21;  // C(7,5)

    // Position → DOF indices mapping (pre-computed for performance)
    private Map<Integer, List<Integer>> positionToDOFIndices;
    private ParametricMolecule pmol;

    // Statistics
    private final CacheStats stats;

    public PartialStartCache() {
        // Generate combinations
        this.sixSubsetCombinations = generateCombinations(TUPLE_SIZE, 6);
        this.fiveSubsetCombinations = generateCombinations(TUPLE_SIZE, 5);

        // Initialize indices
        this.sixSubsetIndices = new ArrayList<>(NUM_6_SUBSETS);
        for (int i = 0; i < NUM_6_SUBSETS; i++) {
            this.sixSubsetIndices.add(new HashMap<>());
        }

        this.fiveSubsetIndices = new ArrayList<>(NUM_5_SUBSETS);
        for (int i = 0; i < NUM_5_SUBSETS; i++) {
            this.fiveSubsetIndices.add(new HashMap<>());
        }

        this.positionToDOFIndices = null;
        this.stats = new CacheStats();
    }

    /**
     * Generate all C(n, k) combinations
     */
    private List<int[]> generateCombinations(int n, int k) {
        List<int[]> result = new ArrayList<>();
        generateCombinationsHelper(n, k, 0, new int[k], 0, result);
        return result;
    }

    private void generateCombinationsHelper(int n, int k, int start, int[] current, int idx, List<int[]> result) {
        if (idx == k) {
            result.add(current.clone());
            return;
        }
        for (int i = start; i <= n - (k - idx); i++) {
            current[idx] = i;
            generateCombinationsHelper(n, k, i + 1, current, idx + 1, result);
        }
    }

    /**
     * Initialize cache with ParametricMolecule
     * MUST be called before use to build position→DOF mapping
     */
    public void initialize(ParametricMolecule pmol) {
        this.pmol = pmol;
        buildPositionToDOFMapping();
    }

    /**
     * Get warm start DOF values for a 7-tuple
     *
     * Strategy:
     * 1. Try 6-subset indices first (7 lookups, 6 positions overlap)
     * 2. Try 5-subset indices second (21 lookups, 5 positions overlap)
     *
     * Complexity: O(28) hash lookups
     */
    public double[] getWarmStart(RCTuple tuple, DoubleMatrix1D defaultDOFs) {
        if (pmol == null || positionToDOFIndices == null) {
            throw new IllegalStateException("PartialStartCache not initialized. Call initialize() first.");
        }

        if (tuple.size() != TUPLE_SIZE) {
            stats.misses++;
            return null;  // Only support 7-tuples
        }

        // Priority 1: Try 6-subset indices (more overlap = better warm start)
        for (int i = 0; i < NUM_6_SUBSETS; i++) {
            int[] subsetPositions = sixSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            CachedEntry cached = sixSubsetIndices.get(i).get(subsetKey);
            if (cached != null) {
                stats.sixSubsetHits++;
                double[] result = buildWarmStartDOFs(tuple, cached, subsetPositions, defaultDOFs);
                System.out.println("[WARM_START] 6-subset hit for tuple: " + tuple.stringListing() + " (6/7 overlap)");
                return result;
            }
        }

        // Priority 2: Try 5-subset indices
        for (int i = 0; i < NUM_5_SUBSETS; i++) {
            int[] subsetPositions = fiveSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            CachedEntry cached = fiveSubsetIndices.get(i).get(subsetKey);
            if (cached != null) {
                stats.fiveSubsetHits++;
                double[] result = buildWarmStartDOFs(tuple, cached, subsetPositions, defaultDOFs);
                System.out.println("[WARM_START] 5-subset hit for tuple: " + tuple.stringListing() + " (5/7 overlap) **LOW OVERLAP**");
                return result;
            }
        }

        // No match found
        stats.misses++;
        return null;
    }

    /**
     * Get warm start with unmatched DOF information (for priority optimization)
     *
     * Returns WarmStartResult containing:
     * - DOF values (same as getWarmStart)
     * - Indices of DOFs for unmatched positions
     * - Subset size (5 or 6)
     */
    public WarmStartResult getWarmStartWithInfo(RCTuple tuple, DoubleMatrix1D defaultDOFs) {
        if (pmol == null || positionToDOFIndices == null) {
            throw new IllegalStateException("PartialStartCache not initialized. Call initialize() first.");
        }

        if (tuple.size() != TUPLE_SIZE) {
            stats.misses++;
            return null;  // Only support 7-tuples
        }

        // Priority 1: Try 6-subset indices (more overlap = better warm start)
        for (int i = 0; i < NUM_6_SUBSETS; i++) {
            int[] subsetPositions = sixSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            CachedEntry cached = sixSubsetIndices.get(i).get(subsetKey);
            if (cached != null) {
                stats.sixSubsetHits++;
                double[] warmStart = buildWarmStartDOFs(tuple, cached, subsetPositions, defaultDOFs);
                int[] unmatchedDOFs = findUnmatchedDOFs(tuple, subsetPositions);
                System.out.println("[WARM_START] 6-subset hit for tuple: " + tuple.stringListing()
                    + " (6/7 overlap) UnmatchedDOFs=" + unmatchedDOFs.length);
                return new WarmStartResult(warmStart, unmatchedDOFs, 6);
            }
        }

        // Priority 2: Try 5-subset indices
        for (int i = 0; i < NUM_5_SUBSETS; i++) {
            int[] subsetPositions = fiveSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            CachedEntry cached = fiveSubsetIndices.get(i).get(subsetKey);
            if (cached != null) {
                stats.fiveSubsetHits++;
                double[] warmStart = buildWarmStartDOFs(tuple, cached, subsetPositions, defaultDOFs);
                int[] unmatchedDOFs = findUnmatchedDOFs(tuple, subsetPositions);
                System.out.println("[WARM_START] 5-subset hit for tuple: " + tuple.stringListing()
                    + " (5/7 overlap) **LOW OVERLAP** UnmatchedDOFs=" + unmatchedDOFs.length);
                return new WarmStartResult(warmStart, unmatchedDOFs, 5);
            }
        }

        // No match found
        stats.misses++;
        return null;
    }

    /**
     * Find DOF indices for unmatched positions
     *
     * Returns array of DOF indices corresponding to positions NOT in the subset
     */
    private int[] findUnmatchedDOFs(RCTuple tuple, int[] subsetPositions) {
        // Create set for fast lookup
        Set<Integer> subsetPosSet = new HashSet<>();
        for (int pos : subsetPositions) {
            subsetPosSet.add(pos);
        }

        List<Integer> unmatchedDOFs = new ArrayList<>();

        // For each position in the tuple
        for (int pos = 0; pos < tuple.size(); pos++) {
            if (!subsetPosSet.contains(pos)) {
                // This position is NOT in the subset - its DOFs are unmatched
                List<Integer> dofIndices = positionToDOFIndices.get(pos);
                if (dofIndices != null) {
                    unmatchedDOFs.addAll(dofIndices);
                }
            }
        }

        return unmatchedDOFs.stream().mapToInt(Integer::intValue).toArray();
    }

    /**
     * Cache a minimized 7-tuple
     *
     * Strategy: Split into 6-subsets and 5-subsets, store if not already present
     */
    public void put(RCTuple tuple, DoubleMatrix1D dofValues) {
        if (pmol == null || positionToDOFIndices == null) {
            return;  // Not initialized, skip caching
        }

        if (tuple.size() != TUPLE_SIZE) {
            return;  // Only support 7-tuples
        }

        // Extract DOF values by position (do this once)
        Map<Integer, double[]> dofsByPosition = extractDOFsByPosition(tuple, dofValues);

        // Create cache entry
        CachedEntry entry = new CachedEntry(tuple, dofsByPosition);

        int stored6 = 0;
        int stored5 = 0;

        // Store in 6-subset indices (skip if already has entry)
        for (int i = 0; i < NUM_6_SUBSETS; i++) {
            int[] subsetPositions = sixSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            Map<String, CachedEntry> index = sixSubsetIndices.get(i);
            if (!index.containsKey(subsetKey)) {
                index.put(subsetKey, entry);
                stored6++;
            }
        }

        // Store in 5-subset indices (skip if already has entry)
        for (int i = 0; i < NUM_5_SUBSETS; i++) {
            int[] subsetPositions = fiveSubsetCombinations.get(i);
            String subsetKey = buildSubsetKey(tuple, subsetPositions);

            Map<String, CachedEntry> index = fiveSubsetIndices.get(i);
            if (!index.containsKey(subsetKey)) {
                index.put(subsetKey, entry);
                stored5++;
            }
        }

        stats.totalStored++;
        stats.sixSubsetEntriesAdded += stored6;
        stats.fiveSubsetEntriesAdded += stored5;
    }

    /**
     * Build subset key from tuple using specified positions
     */
    private String buildSubsetKey(RCTuple tuple, int[] positions) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < positions.length; i++) {
            if (i > 0) sb.append(':');
            int pos = positions[i];
            sb.append(tuple.RCs.get(pos));
        }
        return sb.toString();
    }

    /**
     * Build warm start DOF array from cached entry
     *
     * For positions in subset: use cached DOF values (guaranteed to match)
     * For positions outside subset: use default DOF values
     */
    private double[] buildWarmStartDOFs(
        RCTuple queryTuple,
        CachedEntry cachedEntry,
        int[] subsetPositions,
        DoubleMatrix1D defaultDOFs
    ) {
        double[] warmStart = toArray(defaultDOFs);  // Start from default

        // Create set for fast lookup
        Set<Integer> subsetPosSet = new HashSet<>();
        for (int pos : subsetPositions) {
            subsetPosSet.add(pos);
        }

        // For each position in the tuple
        for (int pos = 0; pos < queryTuple.size(); pos++) {
            if (subsetPosSet.contains(pos)) {
                // This position is in the subset - RCs are guaranteed to match
                // Use cached DOF values
                double[] cachedDOFsForPos = cachedEntry.dofsByPosition.get(pos);
                List<Integer> dofIndices = positionToDOFIndices.get(pos);

                if (cachedDOFsForPos != null && dofIndices != null) {
                    for (int i = 0; i < Math.min(dofIndices.size(), cachedDOFsForPos.length); i++) {
                        int dofIdx = dofIndices.get(i);
                        if (dofIdx < warmStart.length) {
                            warmStart[dofIdx] = cachedDOFsForPos[i];
                        }
                    }
                }
            }
            // else: Position not in subset - keep default value
        }

        return warmStart;
    }

    /**
     * Extract DOF values organized by position
     */
    private Map<Integer, double[]> extractDOFsByPosition(
        RCTuple tuple,
        DoubleMatrix1D dofValues
    ) {
        Map<Integer, double[]> result = new HashMap<>();

        for (int pos = 0; pos < tuple.size(); pos++) {
            List<Integer> dofIndices = positionToDOFIndices.get(pos);

            if (dofIndices != null && !dofIndices.isEmpty()) {
                double[] posDOFs = new double[dofIndices.size()];

                for (int i = 0; i < dofIndices.size(); i++) {
                    int dofIdx = dofIndices.get(i);
                    if (dofIdx < dofValues.size()) {
                        posDOFs[i] = dofValues.get(dofIdx);
                    }
                }

                result.put(pos, posDOFs);
            }
        }

        return result;
    }

    /**
     * Build mapping from position index to DOF indices
     */
    private void buildPositionToDOFMapping() {
        positionToDOFIndices = new HashMap<>();

        if (pmol == null || pmol.dofs == null || pmol.dofs.isEmpty()) {
            return;
        }

        for (int dofIdx = 0; dofIdx < pmol.dofs.size(); dofIdx++) {
            DegreeOfFreedom dof = pmol.dofs.get(dofIdx);
            Integer pos = getDOFPosition(dof);

            if (pos != null) {
                positionToDOFIndices
                    .computeIfAbsent(pos, k -> new ArrayList<>())
                    .add(dofIdx);
            }
        }
    }

    /**
     * Extract position index from a DOF
     */
    private Integer getDOFPosition(DegreeOfFreedom dof) {
        if (dof instanceof FreeDihedral) {
            FreeDihedral dihedralDOF = (FreeDihedral) dof;
            Residue residue = dihedralDOF.getResidue();

            if (residue != null) {
                return residue.indexInMolecule;
            }
        }
        return null;
    }

    /**
     * Get cache statistics
     */
    public CacheStats getStats() {
        return stats;
    }

    /**
     * Clear all cached entries
     */
    public void clear() {
        for (Map<String, CachedEntry> index : sixSubsetIndices) {
            index.clear();
        }
        for (Map<String, CachedEntry> index : fiveSubsetIndices) {
            index.clear();
        }
        stats.reset();
    }

    /**
     * Get total number of unique cached entries
     */
    public int size() {
        Set<CachedEntry> uniqueEntries = new HashSet<>();
        for (Map<String, CachedEntry> index : sixSubsetIndices) {
            uniqueEntries.addAll(index.values());
        }
        for (Map<String, CachedEntry> index : fiveSubsetIndices) {
            uniqueEntries.addAll(index.values());
        }
        return uniqueEntries.size();
    }

    /**
     * Get total index entries count
     */
    public int totalIndexEntries() {
        int total = 0;
        for (Map<String, CachedEntry> index : sixSubsetIndices) {
            total += index.size();
        }
        for (Map<String, CachedEntry> index : fiveSubsetIndices) {
            total += index.size();
        }
        return total;
    }

    // Helper: Convert DoubleMatrix1D to array
    private static double[] toArray(DoubleMatrix1D matrix) {
        double[] arr = new double[matrix.size()];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = matrix.get(i);
        }
        return arr;
    }

    /**
     * Cached entry containing minimized DOF values by position
     */
    private static class CachedEntry {
        final RCTuple tuple;
        final Map<Integer, double[]> dofsByPosition;

        CachedEntry(RCTuple tuple, Map<Integer, double[]> dofsByPosition) {
            // Defensive copy
            this.tuple = new RCTuple(
                new ArrayList<>(tuple.pos),
                new ArrayList<>(tuple.RCs)
            );
            this.dofsByPosition = dofsByPosition;
        }
    }

    /**
     * Cache statistics
     */
    public static class CacheStats {
        public long sixSubsetHits = 0;
        public long fiveSubsetHits = 0;
        public long misses = 0;
        public long totalStored = 0;
        public long sixSubsetEntriesAdded = 0;
        public long fiveSubsetEntriesAdded = 0;

        public void reset() {
            sixSubsetHits = 0;
            fiveSubsetHits = 0;
            misses = 0;
            totalStored = 0;
            sixSubsetEntriesAdded = 0;
            fiveSubsetEntriesAdded = 0;
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
            return String.format(
                "PartialStartCache Stats: queries=%d, 6-hits=%d, 5-hits=%d, miss=%d, hitRate=%.1f%%, " +
                "stored=%d, 6-entries=%d, 5-entries=%d",
                totalQueries(), sixSubsetHits, fiveSubsetHits, misses, hitRate() * 100,
                totalStored, sixSubsetEntriesAdded, fiveSubsetEntriesAdded
            );
        }
    }
}
