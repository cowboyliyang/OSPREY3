package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

import java.util.Arrays;
import java.util.Set;

/**
 * Constrained Minimizer: Minimize only a subset of DOFs while keeping others fixed
 *
 * This is the key component for true subtree DOF caching in Phase 2.
 * It wraps an existing Minimizer and adds the ability to optimize only
 * specified DOF indices while keeping the rest fixed.
 *
 * Implementation Strategy:
 * 1. Extract free DOFs (those to optimize) from full DOF vector
 * 2. Create a constrained objective function that only varies free DOFs
 * 3. Minimize using the underlying minimizer
 * 4. Reconstruct full DOF vector with optimized free DOFs
 *
 * This enables subtree minimization:
 * - Minimize subtree A: optimize only DOFs for positions {0,1,2}
 * - Minimize subtree B: optimize only DOFs for positions {3,4}
 * - Combine: merge both results, then refine boundaries
 */
public class ConstrainedMinimizer implements Minimizer {

    private final Minimizer delegate;
    private final ObjectiveFunction fullObjective;
    private final Set<Integer> freeDOFIndices; // Indices of DOFs to optimize
    private final DoubleMatrix1D fixedDOFs;    // Full DOF vector with fixed values

    /**
     * Create a constrained minimizer
     *
     * @param delegate The underlying minimizer to use
     * @param fullObjective The objective function over all DOFs
     * @param freeDOFIndices Indices of DOFs to optimize (others are fixed)
     * @param fixedDOFs Initial DOF values (free DOFs will be extracted, fixed DOFs kept)
     */
    public ConstrainedMinimizer(
            Minimizer delegate,
            ObjectiveFunction fullObjective,
            Set<Integer> freeDOFIndices,
            DoubleMatrix1D fixedDOFs) {

        this.delegate = delegate;
        this.fullObjective = fullObjective;
        this.freeDOFIndices = freeDOFIndices;
        this.fixedDOFs = fixedDOFs.copy(); // Deep copy to avoid external modification
    }

    @Override
    public Result minimize() {
        return minimizeFromCenter();
    }

    @Override
    public Result minimizeFromCenter() {
        // Start from current fixed DOF values
        return minimizeFrom(fixedDOFs);
    }

    @Override
    public Result minimizeFrom(DoubleMatrix1D x) {
        if (freeDOFIndices == null || freeDOFIndices.isEmpty()) {
            // No free DOFs - return current state without optimization
            double energy = fullObjective.getValue(x);
            return new Result(x.copy(), energy);
        }

        // Special case: if no free DOFs, just evaluate energy and return
        if (freeDOFIndices.isEmpty()) {
            double energy = fullObjective.getValue(x);
            return new Result(x.copy(), energy);
        }

        if (freeDOFIndices.size() == x.size()) {
            // All DOFs are free - use standard minimization
            return delegate.minimizeFrom(x);
        }

        // Extract free DOFs from full vector
        int numFreeDOFs = freeDOFIndices.size();
        DoubleMatrix1D freeDOFs = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);
        Integer[] sortedIndices = freeDOFIndices.toArray(new Integer[0]);
        Arrays.sort(sortedIndices); // Sort for consistent ordering

        for (int i = 0; i < numFreeDOFs; i++) {
            int fullIdx = sortedIndices[i];
            if (fullIdx >= 0 && fullIdx < x.size()) {
                freeDOFs.set(i, x.get(fullIdx));
            }
        }

        // Create constrained objective function
        ConstrainedObjectiveFunction constrainedObj =
            new ConstrainedObjectiveFunction(fullObjective, sortedIndices, x.copy());

        // Create temporary minimizer for free DOFs only
        Minimizer freeMinimizer = createFreeDOFMinimizer(constrainedObj);

        // Minimize free DOFs
        Result freeResult = freeMinimizer.minimizeFrom(freeDOFs);

        // Check if minimization failed
        if (freeResult == null) {
            // Minimization failed, return current DOFs with evaluated energy
            double energy = fullObjective.getValue(x);
            return new Result(x.copy(), energy);
        }

        // Reconstruct full DOF vector
        DoubleMatrix1D fullDOFs = x.copy();
        for (int i = 0; i < numFreeDOFs; i++) {
            int fullIdx = sortedIndices[i];
            if (fullIdx >= 0 && fullIdx < fullDOFs.size()) {
                fullDOFs.set(fullIdx, freeResult.dofValues.get(i));
            }
        }

        // Return result with full DOF vector
        return new Result(fullDOFs, freeResult.energy);
    }

    @Override
    public void clean() {
        if (delegate instanceof NeedsCleanup) {
            ((NeedsCleanup) delegate).clean();
        }
    }

    /**
     * Create a minimizer for free DOFs only
     * Uses the same type as the delegate minimizer
     */
    private Minimizer createFreeDOFMinimizer(ObjectiveFunction constrainedObj) {
        // Try to create same type of minimizer as delegate
        if (delegate instanceof CCDMinimizer) {
            return new CCDMinimizer(constrainedObj, false);
        } else {
            // Default: use CCD minimizer
            return new CCDMinimizer(constrainedObj, false);
        }
    }

    /**
     * Constrained objective function: varies only free DOFs, keeps others fixed
     */
    private static class ConstrainedObjectiveFunction implements ObjectiveFunction {

        private final ObjectiveFunction fullObjective;
        private final Integer[] freeDOFIndices; // Sorted indices of free DOFs
        private final DoubleMatrix1D fixedDOFs;  // Full DOF vector with fixed values

        public ConstrainedObjectiveFunction(
                ObjectiveFunction fullObjective,
                Integer[] freeDOFIndices,
                DoubleMatrix1D fixedDOFs) {

            this.fullObjective = fullObjective;
            this.freeDOFIndices = freeDOFIndices;
            this.fixedDOFs = fixedDOFs;
        }

        @Override
        public int getNumDOFs() {
            return freeDOFIndices.length;
        }

        @Override
        public DoubleMatrix1D[] getConstraints() {
            // Extract constraints for free DOFs only
            DoubleMatrix1D[] fullConstraints = fullObjective.getConstraints();
            if (fullConstraints == null || fullConstraints.length != 2) {
                return null;
            }

            DoubleMatrix1D[] freeConstraints = new DoubleMatrix1D[2];
            int numFreeDOFs = freeDOFIndices.length;

            freeConstraints[0] = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);
            freeConstraints[1] = cern.colt.matrix.DoubleFactory1D.dense.make(numFreeDOFs);

            for (int i = 0; i < numFreeDOFs; i++) {
                int fullIdx = freeDOFIndices[i];
                if (fullIdx >= 0 && fullIdx < fullConstraints[0].size()) {
                    freeConstraints[0].set(i, fullConstraints[0].get(fullIdx));
                    freeConstraints[1].set(i, fullConstraints[1].get(fullIdx));
                }
            }

            return freeConstraints;
        }

        @Override
        public void setDOFs(DoubleMatrix1D freeDOFValues) {
            // Reconstruct full DOF vector: copy fixed DOFs, update free DOFs
            DoubleMatrix1D fullDOFs = fixedDOFs.copy();

            for (int i = 0; i < freeDOFIndices.length && i < freeDOFValues.size(); i++) {
                int fullIdx = freeDOFIndices[i];
                if (fullIdx >= 0 && fullIdx < fullDOFs.size()) {
                    fullDOFs.set(fullIdx, freeDOFValues.get(i));
                }
            }

            // Update full objective function
            fullObjective.setDOFs(fullDOFs);

            // Update our fixed DOF cache
            for (int i = 0; i < fixedDOFs.size(); i++) {
                fixedDOFs.set(i, fullDOFs.get(i));
            }
        }

        public void setDOFsNoCopy(DoubleMatrix1D freeDOFValues) {
            setDOFs(freeDOFValues);
        }

        @Override
        public double getValue(DoubleMatrix1D freeDOFValues) {
            // Reconstruct full DOF vector
            DoubleMatrix1D fullDOFs = fixedDOFs.copy();

            for (int i = 0; i < freeDOFIndices.length && i < freeDOFValues.size(); i++) {
                int fullIdx = freeDOFIndices[i];
                if (fullIdx >= 0 && fullIdx < fullDOFs.size()) {
                    fullDOFs.set(fullIdx, freeDOFValues.get(i));
                }
            }

            // Evaluate full objective
            return fullObjective.getValue(fullDOFs);
        }

        @Override
        public double getValForDOF(int freeDOFIndex, double val) {
            // Map free DOF index to full DOF index
            if (freeDOFIndex >= 0 && freeDOFIndex < freeDOFIndices.length) {
                int fullIdx = freeDOFIndices[freeDOFIndex];
                return fullObjective.getValForDOF(fullIdx, val);
            }
            return Double.NaN;
        }

        @Override
        public void setDOF(int freeDOFIndex, double val) {
            // Map free DOF index to full DOF index and set it
            if (freeDOFIndex >= 0 && freeDOFIndex < freeDOFIndices.length) {
                int fullIdx = freeDOFIndices[freeDOFIndex];
                fullObjective.setDOF(fullIdx, val);
                // Update our fixed DOF cache
                if (fullIdx >= 0 && fullIdx < fixedDOFs.size()) {
                    fixedDOFs.set(fullIdx, val);
                }
            }
        }

        @Override
        public double getInitStepSize(int freeDOFIndex) {
            if (freeDOFIndex >= 0 && freeDOFIndex < freeDOFIndices.length) {
                int fullIdx = freeDOFIndices[freeDOFIndex];
                return fullObjective.getInitStepSize(fullIdx);
            }
            return 0.0;
        }

        @Override
        public boolean isDOFAngle(int freeDOFIndex) {
            if (freeDOFIndex >= 0 && freeDOFIndex < freeDOFIndices.length) {
                int fullIdx = freeDOFIndices[freeDOFIndex];
                return fullObjective.isDOFAngle(fullIdx);
            }
            return false;
        }
    }
}
