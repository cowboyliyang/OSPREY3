package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Constrained DOF Minimizer for Phase 4
 *
 * This minimizer optimizes only a subset of DOFs while keeping others fixed.
 * Used for:
 * 1. L-set minimization: optimize L-set DOFs, keep M-set + rest fixed
 * 2. M-set quick optimization: optimize M-set DOFs, keep L-set fixed
 *
 * Implementation wraps existing minimizer (CCD) and constrains DOF bounds
 * to fix certain positions.
 */
public class ConstrainedDOFMinimizer {

    private final Minimizer baseMinimizer;
    private final Set<Integer> positionsToOptimize;
    private final ObjectiveFunction objectiveFunction;
    private final edu.duke.cs.osprey.confspace.ParametricMolecule pmol;
    private int maxIterations = 30; // Default CCD iterations

    /**
     * Create constrained minimizer
     *
     * @param baseMinimizer Base minimizer (e.g., SimpleCCDMinimizer)
     * @param positionsToOptimize Set of position indices to optimize (others fixed)
     * @param objectiveFunction Energy function
     * @param pmol ParametricMolecule for DOF mapping
     */
    public ConstrainedDOFMinimizer(Minimizer baseMinimizer,
                                   Set<Integer> positionsToOptimize,
                                   ObjectiveFunction objectiveFunction,
                                   edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {
        this.baseMinimizer = baseMinimizer;
        this.positionsToOptimize = positionsToOptimize;
        this.objectiveFunction = objectiveFunction;
        this.pmol = pmol;
    }

    /**
     * Set maximum iterations (for quick M-set optimization)
     */
    public void setMaxIterations(int maxIterations) {
        this.maxIterations = maxIterations;
    }

    /**
     * Minimize with fixed DOFs
     *
     * Strategy:
     * 1. Create wrapper ObjectiveFunction that only exposes free DOFs
     * 2. Set fixed DOFs to have min=max (zero range)
     * 3. Run base minimizer
     * 4. Restore full DOF vector
     */
    public DoubleMatrix1D minimizeFixedDOFs(DoubleMatrix1D initialDOFs) {
        // Get DOF indices to optimize
        List<Integer> freeDOFIndices = getFreeDOFIndices();

        if (freeDOFIndices.isEmpty()) {
            // Nothing to optimize
            return initialDOFs.copy();
        }

        // Create constrained objective function
        ConstrainedObjectiveFunction constrainedObjFunc =
                new ConstrainedObjectiveFunction(objectiveFunction, freeDOFIndices, initialDOFs);

        // Create minimizer for constrained function
        Minimizer minimizer;
        if (baseMinimizer instanceof SimpleCCDMinimizer) {
            minimizer = createConstrainedCCDMinimizer(constrainedObjFunc);
        } else {
            // Fall back to SimpleCCDMinimizer
            minimizer = new SimpleCCDMinimizer(constrainedObjFunc);
        }

        // Minimize constrained DOFs (starting from center of bounds)
        Minimizer.Result result = minimizer.minimizeFromCenter();

        // Restore full DOF vector
        DoubleMatrix1D fullDOFs = initialDOFs.copy();
        for (int i = 0; i < freeDOFIndices.size(); i++) {
            int dofIdx = freeDOFIndices.get(i);
            fullDOFs.set(dofIdx, result.dofValues.get(i));
        }

        return fullDOFs;
    }

    /**
     * Create CCD minimizer with iteration limit
     */
    private Minimizer createConstrainedCCDMinimizer(ObjectiveFunction objFunc) {
        SimpleCCDMinimizer ccd = new SimpleCCDMinimizer(objFunc);

        // Set iteration limit if needed (by modifying SimpleCCDMinimizer.MaxIterations)
        // For now, use default
        // TODO: Add API to SimpleCCDMinimizer for setting max iterations

        return ccd;
    }

    /**
     * Get DOF indices that should be optimized (free DOFs)
     */
    private List<Integer> getFreeDOFIndices() {
        List<Integer> freeDOFs = new ArrayList<>();

        // For each position to optimize, add its DOF indices
        for (int pos : positionsToOptimize) {
            List<Integer> posDOFs = getDOFIndicesForPosition(pos);
            freeDOFs.addAll(posDOFs);
        }

        return freeDOFs;
    }

    /**
     * Get DOF indices for a position using ParametricMolecule
     * This is the CORRECT way to map positions to DOF indices.
     */
    private List<Integer> getDOFIndicesForPosition(int pos) {
        List<Integer> indices = new ArrayList<>();

        if (pmol == null || pmol.dofs == null || pmol.dofs.isEmpty()) {
            return indices; // No DOFs to optimize
        }

        // CRITICAL: Map position indices to DOF indices based on the ACTUAL ParametricMolecule
        // DOFs are ordered by position, but each position may have different number of DOFs
        for (int dofIdx = 0; dofIdx < pmol.dofs.size(); dofIdx++) {
            edu.duke.cs.osprey.dof.DegreeOfFreedom dof = pmol.dofs.get(dofIdx);

            // Determine which position this DOF belongs to
            Integer dofPosition = getDOFPosition(dof);

            if (dofPosition != null && dofPosition == pos) {
                indices.add(dofIdx);
            }
        }

        return indices;
    }

    /**
     * Extract position index from a DOF
     * Returns null if position cannot be determined
     */
    private Integer getDOFPosition(edu.duke.cs.osprey.dof.DegreeOfFreedom dof) {
        // DOFs are usually associated with a residue
        if (dof instanceof edu.duke.cs.osprey.dof.FreeDihedral) {
            edu.duke.cs.osprey.dof.FreeDihedral dihedralDOF = (edu.duke.cs.osprey.dof.FreeDihedral) dof;
            if (dihedralDOF.getResidue() != null) {
                // The residue number directly corresponds to position index
                // This assumes residues are numbered sequentially starting from 0
                // TODO: May need to map PDB residue numbers to position indices via confSpace
                return dihedralDOF.getResidue().indexInMolecule;
            }
        }
        return null;
    }

    /**
     * Wrapper objective function that only exposes free DOFs
     */
    private static class ConstrainedObjectiveFunction implements ObjectiveFunction {

        private final ObjectiveFunction baseFunc;
        private final List<Integer> freeDOFIndices;
        private final DoubleMatrix1D fixedDOFs;
        private final int numFreeDOFs;

        ConstrainedObjectiveFunction(ObjectiveFunction baseFunc,
                                    List<Integer> freeDOFIndices,
                                    DoubleMatrix1D fixedDOFs) {
            this.baseFunc = baseFunc;
            this.freeDOFIndices = freeDOFIndices;
            this.fixedDOFs = fixedDOFs.copy();
            this.numFreeDOFs = freeDOFIndices.size();
        }

        @Override
        public int getNumDOFs() {
            return numFreeDOFs;
        }

        @Override
        public DoubleMatrix1D[] getConstraints() {
            // Get base constraints
            DoubleMatrix1D[] baseConstraints = baseFunc.getConstraints();

            // Extract constraints for free DOFs only
            DoubleMatrix1D minConstraints = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);
            DoubleMatrix1D maxConstraints = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);

            for (int i = 0; i < numFreeDOFs; i++) {
                int dofIdx = freeDOFIndices.get(i);
                minConstraints.set(i, baseConstraints[0].get(dofIdx));
                maxConstraints.set(i, baseConstraints[1].get(dofIdx));
            }

            return new DoubleMatrix1D[]{minConstraints, maxConstraints};
        }

        @Override
        public void setDOFs(DoubleMatrix1D dofs) {
            // Expand to full DOF vector
            DoubleMatrix1D fullDOFs = createFullDOFs(dofs);
            baseFunc.setDOFs(fullDOFs);
        }

        @Override
        public void setDOF(int dof, double val) {
            // Map to full DOF index
            int fullDOFIdx = freeDOFIndices.get(dof);
            baseFunc.setDOF(fullDOFIdx, val);
        }

        @Override
        public double getValue(DoubleMatrix1D dofs) {
            DoubleMatrix1D fullDOFs = createFullDOFs(dofs);
            return baseFunc.getValue(fullDOFs);
        }

        @Override
        public double getValForDOF(int dof, double val) {
            // Map to full DOF index
            int fullDOFIdx = freeDOFIndices.get(dof);
            return baseFunc.getValForDOF(fullDOFIdx, val);
        }

        @Override
        public double getInitStepSize(int dof) {
            int fullDOFIdx = freeDOFIndices.get(dof);
            return baseFunc.getInitStepSize(fullDOFIdx);
        }

        @Override
        public boolean isDOFAngle(int dof) {
            int fullDOFIdx = freeDOFIndices.get(dof);
            return baseFunc.isDOFAngle(fullDOFIdx);
        }

        @Override
        public DoubleMatrix1D getDOFsCenter() {
            // Get center for free DOFs only
            DoubleMatrix1D fullCenter = baseFunc.getDOFsCenter();
            DoubleMatrix1D freeCenter = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);

            for (int i = 0; i < numFreeDOFs; i++) {
                int fullDOFIdx = freeDOFIndices.get(i);
                freeCenter.set(i, fullCenter.get(fullDOFIdx));
            }

            return freeCenter;
        }

        /**
         * Create full DOF vector from free DOF values
         */
        private DoubleMatrix1D createFullDOFs(DoubleMatrix1D freeDOFs) {
            DoubleMatrix1D fullDOFs = fixedDOFs.copy();

            for (int i = 0; i < numFreeDOFs; i++) {
                int fullDOFIdx = freeDOFIndices.get(i);
                fullDOFs.set(fullDOFIdx, freeDOFs.get(i));
            }

            return fullDOFs;
        }
    }
}
