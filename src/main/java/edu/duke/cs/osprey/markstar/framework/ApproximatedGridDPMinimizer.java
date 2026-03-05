package edu.duke.cs.osprey.markstar.framework;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;

import java.util.HashMap;
import java.util.Map;

/**
 * GridDP minimizer variant that evaluates one-body/pairwise grid tables
 * with a confspace-specific approximator matrix.
 *
 * Falls back to forcefield energy if an approximator term is missing.
 */
public class ApproximatedGridDPMinimizer extends GridDPMinimizer {

    private final ApproximatorMatrix amat;
    private final boolean fallbackToForcefield;

    private final Map<Long, ApproximatedObjectiveFunction.Approximator.Addable> oneBodyApproximators = new HashMap<>();
    private final Map<Long, ApproximatedObjectiveFunction.Approximator.Addable> pairApproximators = new HashMap<>();

    public ApproximatedGridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                                       RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                                       SimpleReferenceEnergies eref, boolean useCache, int numThreads,
                                       ApproximatorMatrix amat, boolean fallbackToForcefield) {
        super(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, useCache, numThreads);
        this.amat = amat;
        this.fallbackToForcefield = fallbackToForcefield;
    }

    public ApproximatedGridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                                       RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                                       SimpleReferenceEnergies eref, boolean useCache,
                                       ApproximatorMatrix amat, boolean fallbackToForcefield) {
        this(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, useCache, 1, amat, fallbackToForcefield);
    }

    @Override
    protected double evalOneBodyEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                       int pos, int rc, int gridState, double[] posDofValues,
                                       double erefOffset) {

        ApproximatedObjectiveFunction.Approximator.Addable approximator = getOneBodyApproximator(pos, rc);
        if (approximator == null) {
            return fallbackOneBody(mol, efunc, pos, rc, gridState, posDofValues, erefOffset);
        }

        if (approximator.numDofs() != posDofValues.length) {
            if (fallbackToForcefield) {
                return fallbackOneBody(mol, efunc, pos, rc, gridState, posDofValues, erefOffset);
            }
            throw new IllegalStateException(String.format(
                    "One-body approximator DOF mismatch at pos=%d rc=%d: expected %d, got %d",
                    pos, rc, approximator.numDofs(), posDofValues.length
            ));
        }

        DoubleMatrix1D x = DoubleFactory1D.dense.make(posDofValues);
        return approximator.getValue(x) + erefOffset;
    }

    @Override
    protected double evalPairEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                    int pos1, int rc1, int gridState1, double[] pos1DofValues,
                                    int pos2, int rc2, int gridState2, double[] pos2DofValues) {

        ApproximatedObjectiveFunction.Approximator.Addable approximator = getPairApproximator(pos1, rc1, pos2, rc2);
        if (approximator == null) {
            return fallbackPair(mol, efunc, pos1, rc1, gridState1, pos1DofValues, pos2, rc2, gridState2, pos2DofValues);
        }

        double[] pairDofs = concat(pos1DofValues, pos2DofValues);
        if (approximator.numDofs() != pairDofs.length) {
            if (fallbackToForcefield) {
                return fallbackPair(mol, efunc, pos1, rc1, gridState1, pos1DofValues, pos2, rc2, gridState2, pos2DofValues);
            }
            throw new IllegalStateException(String.format(
                    "Pair approximator DOF mismatch at (%d,%d)-(%d,%d): expected %d, got %d",
                    pos1, rc1, pos2, rc2, approximator.numDofs(), pairDofs.length
            ));
        }

        DoubleMatrix1D x = DoubleFactory1D.dense.make(pairDofs);
        return approximator.getValue(x);
    }

    private double fallbackOneBody(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                   int pos, int rc, int gridState, double[] posDofValues,
                                   double erefOffset) {
        return super.evalOneBodyEnergy(mol, efunc, pos, rc, gridState, posDofValues, erefOffset);
    }

    private double fallbackPair(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                int pos1, int rc1, int gridState1, double[] pos1DofValues,
                                int pos2, int rc2, int gridState2, double[] pos2DofValues) {
        return super.evalPairEnergy(mol, efunc, pos1, rc1, gridState1, pos1DofValues, pos2, rc2, gridState2, pos2DofValues);
    }

    private ApproximatedObjectiveFunction.Approximator.Addable getOneBodyApproximator(int pos, int rc) {
        if (amat == null) return null;
        long key = oneBodyKey(pos, rc);
        if (!oneBodyApproximators.containsKey(key)) {
            oneBodyApproximators.put(key, amat.get(pos, rc));
        }
        return oneBodyApproximators.get(key);
    }

    private ApproximatedObjectiveFunction.Approximator.Addable getPairApproximator(int pos1, int rc1, int pos2, int rc2) {
        if (amat == null) return null;
        long key = pairKey(pos1, rc1, pos2, rc2);
        if (!pairApproximators.containsKey(key)) {
            pairApproximators.put(key, amat.get(pos1, rc1, pos2, rc2));
        }
        return pairApproximators.get(key);
    }

    private static long oneBodyKey(int pos, int rc) {
        return ((long) pos << 32) | (rc & 0xFFFFFFFFL);
    }

    private static long pairKey(int pos1, int rc1, int pos2, int rc2) {
        int aPos = Math.min(pos1, pos2);
        int bPos = Math.max(pos1, pos2);
        int aRc = (aPos == pos1) ? rc1 : rc2;
        int bRc = (aPos == pos1) ? rc2 : rc1;
        return ((long) aPos << 48) | ((long) aRc << 32) | ((long) bPos << 16) | (bRc & 0xFFFFL);
    }

    private static double[] concat(double[] a, double[] b) {
        double[] out = new double[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }
}

