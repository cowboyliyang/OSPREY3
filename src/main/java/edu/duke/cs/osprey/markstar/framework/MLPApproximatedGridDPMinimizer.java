package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.approximation.MLPEnergyModel;
import edu.duke.cs.osprey.energy.approximation.MLPSurrogateMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;

/**
 * GridDP variant that uses per-task MLP models for one-body and pairwise energies.
 */
public class MLPApproximatedGridDPMinimizer extends GridDPMinimizer {

    private final MLPSurrogateMatrix surrogate;
    private final boolean fallbackToForcefield;

    public MLPApproximatedGridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                                          RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                                          SimpleReferenceEnergies eref, boolean useCache, int numThreads,
                                          MLPSurrogateMatrix surrogate, boolean fallbackToForcefield) {
        super(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, useCache, numThreads);
        this.surrogate = surrogate;
        this.fallbackToForcefield = fallbackToForcefield;
    }

    public MLPApproximatedGridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                                          RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                                          SimpleReferenceEnergies eref, boolean useCache,
                                          MLPSurrogateMatrix surrogate, boolean fallbackToForcefield) {
        this(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, useCache, 1, surrogate, fallbackToForcefield);
    }

    @Override
    protected double evalOneBodyEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                       int pos, int rc, int gridState, double[] posDofValues,
                                       double erefOffset) {

        MLPEnergyModel model = surrogate == null ? null : surrogate.getOneBody(pos, rc);
        if (model == null) {
            return fallbackOneBody(mol, efunc, pos, rc, gridState, posDofValues, erefOffset);
        }
        if (model.numInputs() != posDofValues.length) {
            if (fallbackToForcefield) {
                return fallbackOneBody(mol, efunc, pos, rc, gridState, posDofValues, erefOffset);
            }
            throw new IllegalStateException(String.format(
                    "MLP one-body input mismatch at pos=%d rc=%d: model=%d dofs=%d",
                    pos, rc, model.numInputs(), posDofValues.length
            ));
        }
        return model.predict(posDofValues) + erefOffset;
    }

    @Override
    protected double evalPairEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                    int pos1, int rc1, int gridState1, double[] pos1DofValues,
                                    int pos2, int rc2, int gridState2, double[] pos2DofValues) {

        boolean swap = pos1 > pos2;
        int aPos = swap ? pos2 : pos1;
        int aRc = swap ? rc2 : rc1;
        int bPos = swap ? pos1 : pos2;
        int bRc = swap ? rc1 : rc2;

        MLPEnergyModel model = surrogate == null ? null : surrogate.getPair(aPos, aRc, bPos, bRc);
        if (model == null) {
            return fallbackPair(mol, efunc, pos1, rc1, gridState1, pos1DofValues, pos2, rc2, gridState2, pos2DofValues);
        }

        double[] in = swap ? concat(pos2DofValues, pos1DofValues) : concat(pos1DofValues, pos2DofValues);
        if (model.numInputs() != in.length) {
            if (fallbackToForcefield) {
                return fallbackPair(mol, efunc, pos1, rc1, gridState1, pos1DofValues, pos2, rc2, gridState2, pos2DofValues);
            }
            throw new IllegalStateException(String.format(
                    "MLP pair input mismatch at (%d,%d)-(%d,%d): model=%d dofs=%d",
                    pos1, rc1, pos2, rc2, model.numInputs(), in.length
            ));
        }
        return model.predict(in);
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

    private static double[] concat(double[] a, double[] b) {
        double[] out = new double[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }
}

