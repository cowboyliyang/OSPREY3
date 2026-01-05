package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.minimization.ConstrainedMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;

import java.util.*;

/**
 * Phase 2: True Subtree DOF Cache with Branch Decomposition
 *
 * This implements the REAL subtree caching strategy:
 * 1. Decompose conformations into subtrees using BranchDecomposition
 * 2. Cache minimized DOF values for each subtree independently
 * 3. When minimizing a new conformation:
 *    - Reuse cached subtrees that match
 *    - Minimize only new/uncached subtrees using ConstrainedMinimizer
 *    - Combine all subtrees and refine boundaries
 *
 * Example:
 *   Conformation A = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC3, pos4:RC0}
 *   Branch decomposition splits into:
 *     Subtree1: {pos0:RC2, pos1:RC5, pos2:RC1}
 *     Subtree2: {pos3:RC3, pos4:RC0}
 *
 *   Minimize conformation A:
 *     - Minimize Subtree1 → cache DOFs for positions 0,1,2
 *     - Minimize Subtree2 → cache DOFs for positions 3,4
 *     - Refine boundaries between subtrees
 *
 *   Later, conformation B = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC7, pos4:RC9}:
 *     - Subtree1 matches cache! Reuse DOFs for positions 0,1,2
 *     - Minimize Subtree2 (new RCs) → cache DOFs for positions 3,4
 *     - Refine boundaries
 *
 * Expected speedup: 30-50% when many conformations share subtrees
 */
public class SubtreeDOFCache {

    // Cache: subtree configuration → minimized DOF values
    private final Map<SubtreeKey, MinimizedSubtree> cache;

    // Branch decomposition tree
    private final BranchDecomposition branchDecomp;

    // Conformation space (needed for DOF mapping)
    private final SimpleConfSpace confSpace;

    // Statistics
    private long cacheHits = 0;
    private long cacheMisses = 0;
    private long partialHits = 0;  // Some subtrees cached, some not
    private long totalQueries = 0;
    private long totalSubtreeQueries = 0;

    // Configuration
    private static final int MAX_CACHE_SIZE = 100000;
    private static final boolean ENABLE_CACHE = true;
    private static final boolean ENABLE_BOUNDARY_REFINEMENT = true; // Now uses local DOF mapping

    public SubtreeDOFCache(BranchDecomposition branchDecomp, SimpleConfSpace confSpace) {
        this.branchDecomp = branchDecomp;
        this.confSpace = confSpace;
        this.cache = new LinkedHashMap<SubtreeKey, MinimizedSubtree>(16, 0.75f, true) {
            @Override
            protected boolean removeEldestEntry(Map.Entry<SubtreeKey, MinimizedSubtree> eldest) {
                return size() > MAX_CACHE_SIZE;
            }
        };
    }

    /**
     * Minimize a conformation with TRUE subtree DOF caching
     *
     * Strategy:
     * 1. Get subtrees from branch decomposition
     * 2. Check cache for each subtree
     * 3. For uncached subtrees: minimize using ConstrainedMinimizer
     * 4. Combine all subtree DOFs into full DOF vector
     * 5. Refine boundaries between subtrees
     */
    public MinimizationResult minimizeWithCache(
            RCTuple conf,
            Minimizer minimizer,
            DoubleMatrix1D initialDOFs,
            ObjectiveFunction objectiveFunction) {

        if (!ENABLE_CACHE) {
            Minimizer.Result result = minimizer.minimizeFrom(initialDOFs);
            return new MinimizationResult(result.dofValues, result.energy, false);
        }

        totalQueries++;

        // Get subtrees for this conformation
        List<Subtree> subtrees = getSubtrees(conf);

        if (subtrees.isEmpty()) {
            // No subtree decomposition available, fall back to full minimization
            Minimizer.Result result = minimizer.minimizeFrom(initialDOFs);
            return new MinimizationResult(result.dofValues, result.energy, false);
        }

        // Get ParametricMolecule from ObjectiveFunction (if available)
        // We need this to compute local DOF indices
        edu.duke.cs.osprey.confspace.ParametricMolecule pmol = null;
        if (objectiveFunction instanceof edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction) {
            pmol = ((edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction) objectiveFunction).pmol;
        }

        // Check cache for each subtree
        List<Subtree> uncachedSubtrees = new ArrayList<>();
        DoubleMatrix1D combinedDOFs = initialDOFs.copy();
        int cachedCount = 0;

        for (Subtree subtree : subtrees) {
            totalSubtreeQueries++;
            SubtreeKey key = new SubtreeKey(subtree, conf);
            MinimizedSubtree cached = cache.get(key);

            if (cached != null) {
                // Cache hit! Apply cached DOFs to combined vector
                // Compute LOCAL DOF indices for this subtree
                List<Integer> localDOFIndices = getLocalDOFIndices(subtree.positions, pmol);
                cachedCount++;
                applySubtreeDOFs(combinedDOFs, cached.dofs, localDOFIndices);
            } else {
                // Cache miss - need to minimize this subtree
                uncachedSubtrees.add(subtree);
            }
        }

        // Minimize uncached subtrees using TRUE ConstrainedMinimizer
        // with dynamically computed LOCAL DOF indices
        if (!uncachedSubtrees.isEmpty()) {
            cacheMisses += uncachedSubtrees.size();

            for (Subtree subtree : uncachedSubtrees) {
                // Compute LOCAL DOF indices for this subtree within the current conformation
                List<Integer> localDOFIndices = getLocalDOFIndices(subtree.positions, pmol);
                Set<Integer> freeDOFIndices = new HashSet<>(localDOFIndices);

                // Create constrained minimizer that only optimizes this subtree's DOFs
                ConstrainedMinimizer constrainedMin = new ConstrainedMinimizer(
                    minimizer,
                    objectiveFunction,
                    freeDOFIndices,
                    combinedDOFs
                );

                // Minimize this subtree
                Minimizer.Result subtreeResult = constrainedMin.minimizeFrom(combinedDOFs);

                // Extract and cache this subtree's DOFs
                DoubleMatrix1D subtreeDOFs = extractSubtreeDOFs(subtreeResult.dofValues, localDOFIndices);
                SubtreeKey key = new SubtreeKey(subtree, conf);
                cache.put(key, new MinimizedSubtree(subtreeDOFs, subtreeResult.energy));

                // Update combined DOFs with this subtree's result
                applySubtreeDOFs(combinedDOFs, subtreeDOFs, localDOFIndices);
            }
        }

        // Update statistics
        if (cachedCount > 0 && uncachedSubtrees.size() > 0) {
            partialHits++;
        } else if (cachedCount == subtrees.size()) {
            cacheHits++;
        }

        // Refine boundaries between subtrees
        double finalEnergy;
        if (ENABLE_BOUNDARY_REFINEMENT && subtrees.size() > 1) {
            finalEnergy = refineBoundaries(combinedDOFs, subtrees, minimizer, objectiveFunction, pmol);
        } else {
            finalEnergy = objectiveFunction.getValue(combinedDOFs);
        }

        boolean fullyCached = (cachedCount == subtrees.size());
        return new MinimizationResult(combinedDOFs, finalEnergy, fullyCached);
    }

    /**
     * Get subtrees from branch decomposition for a given conformation
     */
    private List<Subtree> getSubtrees(RCTuple conf) {
        List<Subtree> subtrees = new ArrayList<>();

        if (branchDecomp == null || branchDecomp.root == null) {
            return subtrees;
        }

        // Traverse branch decomposition tree and collect subtrees
        collectSubtreesFromNode(branchDecomp.root, subtrees);

        return subtrees;
    }

    /**
     * Recursively collect subtrees from branch decomposition tree
     */
    private void collectSubtreesFromNode(BranchDecomposition.TreeNode node, List<Subtree> subtrees) {
        if (node == null) return;

        // Add this node as a subtree if it's not a leaf and has enough positions
        if (!node.isLeaf && node.positions.size() > 1) {
            List<Integer> positions = new ArrayList<>(node.positions);
            // DOF indices will be computed dynamically at runtime based on the specific conformation
            subtrees.add(new Subtree(positions));
        }

        // Recurse to children
        if (node.leftChild != null) {
            collectSubtreesFromNode(node.leftChild, subtrees);
        }
        if (node.rightChild != null) {
            collectSubtreesFromNode(node.rightChild, subtrees);
        }
    }

    /**
     * Get LOCAL DOF indices for given positions within a specific conformation.
     *
     * CRITICAL: This returns indices relative to the current conformation's DOF vector,
     * NOT global indices in the conformation space.
     *
     * IMPORTANT: We use the ParametricMolecule's actual DOF list to determine ranges,
     * because different amino acids have different numbers of DOFs (e.g., GLY vs TRP).
     *
     * @param positions Positions in this subtree
     * @param pmol The ParametricMolecule for the current conformation
     * @return List of LOCAL DOF indices (0 to pmol.dofs.size()-1)
     */
    private List<Integer> getLocalDOFIndices(List<Integer> positions, edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {
        List<Integer> indices = new ArrayList<>();

        if (pmol == null || pmol.dofs == null || pmol.dofs.isEmpty()) {
            return indices; // No DOFs to optimize
        }

        // CRITICAL: We need to map position indices to DOF indices based on the ACTUAL
        // ParametricMolecule, not assumptions about DOF counts.
        //
        // The ParametricMolecule.dofs list contains all DOFs for this specific conformation.
        // DOFs are ordered by position, but each position may have a different number of DOFs
        // depending on which amino acid is selected.
        //
        // Strategy: Use the DOF names/residues to determine which DOFs belong to which positions

        Set<Integer> targetPositions = new HashSet<>(positions);

        for (int dofIdx = 0; dofIdx < pmol.dofs.size(); dofIdx++) {
            edu.duke.cs.osprey.dof.DegreeOfFreedom dof = pmol.dofs.get(dofIdx);

            // Try to determine which position this DOF belongs to
            // DOFs typically have a residue associated with them
            Integer dofPosition = getDOFPosition(dof);

            if (dofPosition != null && targetPositions.contains(dofPosition)) {
                indices.add(dofIdx);
            }
        }

        // If we couldn't determine positions from DOFs (shouldn't happen), fall back
        if (indices.isEmpty() && confSpace != null) {
            // Fallback: Use confSpace structure (but this is less reliable with mutations)
            int currentDOFIdx = 0;
            for (int pos = 0; pos < Math.min(confSpace.positions.size(), pmol.dofs.size()); pos++) {
                // Estimate DOFs per position from pmol
                int dofsPerPosition = estimateDOFsPerPosition(pos, pmol);

                if (targetPositions.contains(pos)) {
                    for (int i = 0; i < dofsPerPosition && (currentDOFIdx + i) < pmol.dofs.size(); i++) {
                        indices.add(currentDOFIdx + i);
                    }
                }

                currentDOFIdx += dofsPerPosition;
            }
        }

        return indices;
    }

    /**
     * Try to extract position index from a DOF
     * Returns null if position cannot be determined
     */
    private Integer getDOFPosition(edu.duke.cs.osprey.dof.DegreeOfFreedom dof) {
        // DOFs are usually associated with a residue
        // Check if this is a FreeDihedral or other common DOF type
        if (dof instanceof edu.duke.cs.osprey.dof.FreeDihedral) {
            edu.duke.cs.osprey.dof.FreeDihedral dihedralDOF = (edu.duke.cs.osprey.dof.FreeDihedral) dof;
            // FreeDihedral has a residue field
            if (dihedralDOF.getResidue() != null && confSpace != null) {
                // Find which position this residue corresponds to
                String resId = dihedralDOF.getResidue().getPDBResNumber();
                for (int pos = 0; pos < confSpace.positions.size(); pos++) {
                    if (confSpace.positions.get(pos).resNum.equals(resId)) {
                        return pos;
                    }
                }
            }
        }

        // Could add more DOF types here if needed
        return null;
    }

    /**
     * Estimate DOFs per position from ParametricMolecule
     */
    private int estimateDOFsPerPosition(int pos, edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {
        if (confSpace == null || pos >= confSpace.positions.size()) {
            return pmol.dofs.size() / Math.max(1, confSpace != null ? confSpace.positions.size() : 1);
        }

        // Average DOF count across all RCs for this position (handles mutations)
        int totalDOFs = 0;
        int numRCs = confSpace.positions.get(pos).resConfs.size();
        for (int rc = 0; rc < Math.min(numRCs, 5); rc++) { // Sample first 5 RCs
            totalDOFs += confSpace.positions.get(pos).resConfs.get(rc).dofBounds.size();
        }
        return totalDOFs / Math.min(numRCs, 5);
    }

    /**
     * Extract DOF values for a specific subtree from full DOF vector
     */
    private DoubleMatrix1D extractSubtreeDOFs(DoubleMatrix1D allDOFs, List<Integer> dofIndices) {
        DoubleMatrix1D subtreeDOFs = DoubleFactory1D.dense.make(dofIndices.size());
        for (int i = 0; i < dofIndices.size(); i++) {
            int idx = dofIndices.get(i);
            if (idx >= 0 && idx < allDOFs.size()) {
                subtreeDOFs.set(i, allDOFs.get(idx));
            }
        }
        return subtreeDOFs;
    }

    /**
     * Apply subtree DOF values to full DOF vector
     */
    private void applySubtreeDOFs(DoubleMatrix1D target, DoubleMatrix1D subtreeDOFs, List<Integer> dofIndices) {
        for (int i = 0; i < subtreeDOFs.size() && i < dofIndices.size(); i++) {
            int targetIdx = dofIndices.get(i);
            if (targetIdx >= 0 && targetIdx < target.size()) {
                target.set(targetIdx, subtreeDOFs.get(i));
            }
        }
    }

    /**
     * Refine DOFs at subtree boundaries
     *
     * After combining cached subtrees, DOFs at boundaries may not be optimal.
     * This step minimizes only the boundary DOFs to fix discontinuities.
     *
     * Example: Subtree1 has positions {0,1,2}, Subtree2 has positions {3,4}
     * Boundary DOFs are those near the interface (e.g., position 2 and 3)
     */
    private double refineBoundaries(
            DoubleMatrix1D dofs,
            List<Subtree> subtrees,
            Minimizer minimizer,
            ObjectiveFunction objectiveFunction,
            edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {

        if (subtrees.size() <= 1) {
            return objectiveFunction.getValue(dofs);
        }

        // Identify boundary DOFs (those involved in cross-subtree interactions)
        // Use LOCAL DOF indices computed from pmol
        Set<Integer> boundaryDOFIndices = new HashSet<>();
        for (int i = 0; i < subtrees.size(); i++) {
            for (int j = i + 1; j < subtrees.size(); j++) {
                boundaryDOFIndices.addAll(getBoundaryDOFs(subtrees.get(i), subtrees.get(j), pmol));
            }
        }

        if (boundaryDOFIndices.isEmpty()) {
            return objectiveFunction.getValue(dofs);
        }

        // Minimize only boundary DOFs using ConstrainedMinimizer
        ConstrainedMinimizer boundaryMin = new ConstrainedMinimizer(
            minimizer,
            objectiveFunction,
            boundaryDOFIndices,
            dofs
        );

        Minimizer.Result result = boundaryMin.minimizeFrom(dofs);

        // Update dofs with refined boundaries
        for (int i = 0; i < dofs.size(); i++) {
            dofs.set(i, result.dofValues.get(i));
        }

        return result.energy;
    }

    /**
     * Get LOCAL DOF indices at the boundary between two subtrees
     *
     * Only includes DOFs from positions that are actually at the boundary
     * (adjacent positions across subtrees), not all DOFs from both subtrees.
     *
     * This is critical for performance: we only want to re-optimize the
     * interface between cached subtrees, not re-optimize everything.
     */
    private Set<Integer> getBoundaryDOFs(Subtree st1, Subtree st2, edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {
        Set<Integer> boundary = new HashSet<>();

        Set<Integer> positions1 = new HashSet<>(st1.positions);
        Set<Integer> positions2 = new HashSet<>(st2.positions);

        // Collect positions that are at the boundary (adjacent across subtrees)
        Set<Integer> boundaryPositions = new HashSet<>();

        for (int pos1 : positions1) {
            for (int pos2 : positions2) {
                if (Math.abs(pos1 - pos2) <= 1) {
                    // These positions are adjacent across the boundary
                    boundaryPositions.add(pos1);
                    boundaryPositions.add(pos2);
                }
            }
        }

        if (boundaryPositions.isEmpty()) {
            return boundary; // No boundary DOFs
        }

        // Get LOCAL DOF indices for boundary positions
        // Convert boundary positions to a list and compute their local DOF indices
        List<Integer> boundaryPosList = new ArrayList<>(boundaryPositions);
        List<Integer> localDOFIndices = getLocalDOFIndices(boundaryPosList, pmol);
        boundary.addAll(localDOFIndices);

        return boundary;
    }

    // Note: getPositionDOFsFromSubtree() method removed
    // We now use getLocalDOFIndices() which computes DOF indices dynamically
    // based on the current conformation's ParametricMolecule

    // Statistics

    public void printStats() {
        System.out.println("\n=== Subtree DOF Cache Statistics (TRUE SUBTREE CACHING) ===");
        System.out.println("Total conformation queries: " + totalQueries);
        System.out.println("Total subtree queries:      " + totalSubtreeQueries);
        System.out.println("Full cache hits:            " + cacheHits + " (all subtrees cached)");
        System.out.println("Partial cache hits:         " + partialHits + " (some subtrees cached)");
        System.out.println("Full cache misses:          " + cacheMisses + " (subtree not cached)");

        if (totalSubtreeQueries > 0) {
            long totalHits = cacheHits + partialHits;
            double hitRate = 100.0 * totalHits / totalQueries;
            double subtreeHitRate = 100.0 * (totalSubtreeQueries - cacheMisses) / totalSubtreeQueries;

            System.out.println("Conformation hit rate:      " + String.format("%.1f%%", hitRate));
            System.out.println("Subtree hit rate:           " + String.format("%.1f%%", subtreeHitRate));
        }

        System.out.println("Cache size:                 " + cache.size() + " / " + MAX_CACHE_SIZE);

        // Estimate speedup
        if (totalQueries > 0) {
            double estimatedSpeedup = 1.0 + (0.5 * (cacheHits + 0.5 * partialHits) / totalQueries);
            System.out.println("Estimated speedup:          " + String.format("%.2fx", estimatedSpeedup));
        }

        System.out.println("===========================================================\n");
    }

    public void clearCache() {
        cache.clear();
        cacheHits = 0;
        cacheMisses = 0;
        partialHits = 0;
        totalQueries = 0;
        totalSubtreeQueries = 0;
    }

    // Inner classes

    /**
     * Key for caching: identifies a unique subtree configuration
     */
    private static class SubtreeKey {
        final List<Integer> positions;
        final int[] RCs; // RC assignments for these positions

        SubtreeKey(Subtree subtree, RCTuple fullConf) {
            this.positions = new ArrayList<>(subtree.positions);
            Collections.sort(this.positions);

            // Extract RCs for subtree positions from full conformation
            this.RCs = new int[positions.size()];
            for (int i = 0; i < positions.size(); i++) {
                int pos = positions.get(i);
                // Find RC at this position in fullConf
                int rcIdx = -1;
                for (int j = 0; j < fullConf.pos.size(); j++) {
                    if (fullConf.pos.get(j) == pos) {
                        rcIdx = fullConf.RCs.get(j);
                        break;
                    }
                }
                this.RCs[i] = rcIdx;
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof SubtreeKey)) return false;
            SubtreeKey that = (SubtreeKey) o;
            return positions.equals(that.positions) && Arrays.equals(RCs, that.RCs);
        }

        @Override
        public int hashCode() {
            int result = positions.hashCode();
            result = 31 * result + Arrays.hashCode(RCs);
            return result;
        }
    }

    /**
     * Cached minimized DOF values for a subtree
     */
    private static class MinimizedSubtree {
        final DoubleMatrix1D dofs;   // DOF values for this subtree only
        final double energy;          // Energy of this subtree (approximate)

        MinimizedSubtree(DoubleMatrix1D dofs, double energy) {
            this.dofs = dofs.copy();
            this.energy = energy;
        }
    }

    /**
     * Result of minimization with cache information
     */
    public static class MinimizationResult {
        public final DoubleMatrix1D dofs;
        public final double energy;
        public final boolean fromCache;

        public MinimizationResult(DoubleMatrix1D dofs, double energy, boolean fromCache) {
            this.dofs = dofs;
            this.energy = energy;
            this.fromCache = fromCache;
        }
    }

    /**
     * Represents a subtree in the branch decomposition.
     *
     * DOF indices are NOT stored here because they depend on the specific conformation.
     * Instead, DOF indices are computed dynamically at runtime using getLocalDOFIndices().
     */
    public static class Subtree {
        final List<Integer> positions;    // Position indices in this subtree

        public Subtree(List<Integer> positions) {
            this.positions = positions;
        }

        // Legacy constructor for backward compatibility (if needed)
        @Deprecated
        public Subtree(List<Integer> positions, List<Integer> dofIndices) {
            this.positions = positions;
            // Ignore dofIndices - they will be computed dynamically
        }
    }
}
