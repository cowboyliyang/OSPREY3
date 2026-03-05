package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DofInfo;
import edu.duke.cs.osprey.energy.approximation.MLPEnergyModel;
import edu.duke.cs.osprey.energy.approximation.MLPSurrogateMatrix;

import java.util.ArrayList;

/**
 * Objective function that evaluates a tuple energy using one-body and pairwise MLP surrogates only.
 */
public class MLPSurrogateObjectiveFunction implements ObjectiveFunction {

    private static final long serialVersionUID = 839588757735124424L;

    public final ParametricMolecule pmol;
    private final DofInfo dofInfo;
    private final MLPSurrogateMatrix surrogate;
    private final boolean requireFullApproximation;
    private final double erefOffset;
    private final DoubleMatrix1D curDOFVals;

    public MLPSurrogateObjectiveFunction(SimpleConfSpace confSpace, RCTuple tuple, MLPSurrogateMatrix surrogate,
                                         boolean requireFullApproximation, double erefOffset) {
        this.pmol = confSpace.makeMolecule(tuple);
        this.dofInfo = confSpace.makeDofInfo(tuple);
        this.surrogate = surrogate;
        this.requireFullApproximation = requireFullApproximation;
        this.erefOffset = erefOffset;
        this.curDOFVals = DoubleFactory1D.dense.make(pmol.dofs.size());
    }

    @Override
    public int getNumDOFs() {
        return pmol.dofs.size();
    }

    @Override
    public DoubleMatrix1D[] getConstraints() {
        return pmol.dofBounds.getBounds();
    }

    @Override
    public void setDOF(int d, double val) {
        curDOFVals.set(d, val);
        pmol.dofs.get(d).apply(val);
    }

    @Override
    public double getValForDOF(int d, double val) {
        setDOF(d, val);
        return getValue(curDOFVals);
    }

    @Override
    public void setDOFs(DoubleMatrix1D x) {
        curDOFVals.assign(x);
        for (int d = 0; d < x.size(); d++) {
            pmol.dofs.get(d).apply(x.get(d));
        }
    }

    @Override
    public double getValue(DoubleMatrix1D x) {
        setDOFs(x);
        return evalCurrent();
    }

    @Override
    public double getInitStepSize(int d) {
        return MoleculeModifierAndScorer.getInitStepSize(pmol.dofs.get(d));
    }

    @Override
    public boolean isDOFAngle(int d) {
        return MoleculeModifierAndScorer.isDOFAngle(pmol.dofs.get(d));
    }

    @Override
    public ArrayList<Integer> getInitFixableDOFs() {
        throw new UnsupportedOperationException("implement me!");
    }

    private double evalCurrent() {

        double e = erefOffset;

        int nBlocks = dofInfo.size();
        for (int bi = 0; bi < nBlocks; bi++) {
            int pos = dofInfo.positions.get(bi).index;
            int rc = dofInfo.resConfs.get(bi).index;
            MLPEnergyModel one = surrogate.getOneBody(pos, rc);
            if (one == null) {
                if (requireFullApproximation) {
                    throw new IllegalStateException(String.format("missing one-body MLP model for (%d,%d)", pos, rc));
                }
                continue;
            }
            double[] x = blockValues(bi);
            if (one.numInputs() != x.length) {
                throw new IllegalStateException(String.format(
                        "one-body input mismatch for (%d,%d): model=%d dofs=%d",
                        pos, rc, one.numInputs(), x.length
                ));
            }
            e += one.predict(x);
        }

        for (int bi = 0; bi < nBlocks; bi++) {
            int posI = dofInfo.positions.get(bi).index;
            int rcI = dofInfo.resConfs.get(bi).index;
            double[] xI = blockValues(bi);
            for (int bj = bi + 1; bj < nBlocks; bj++) {
                int posJ = dofInfo.positions.get(bj).index;
                int rcJ = dofInfo.resConfs.get(bj).index;
                double[] xJ = blockValues(bj);

                boolean swap = posI > posJ;
                int aPos = swap ? posJ : posI;
                int aRc = swap ? rcJ : rcI;
                int bPos = swap ? posI : posJ;
                int bRc = swap ? rcI : rcJ;

                MLPEnergyModel pair = surrogate.getPair(aPos, aRc, bPos, bRc);
                if (pair == null) {
                    if (requireFullApproximation) {
                        throw new IllegalStateException(String.format(
                                "missing pair MLP model for (%d,%d)-(%d,%d)", aPos, aRc, bPos, bRc
                        ));
                    }
                    continue;
                }

                double[] in = swap ? concat(xJ, xI) : concat(xI, xJ);
                if (pair.numInputs() != in.length) {
                    throw new IllegalStateException(String.format(
                            "pair input mismatch for (%d,%d)-(%d,%d): model=%d dofs=%d",
                            aPos, aRc, bPos, bRc, pair.numInputs(), in.length
                    ));
                }
                e += pair.predict(in);
            }
        }

        return e;
    }

    private double[] blockValues(int block) {
        int offset = dofInfo.offsets.get(block);
        int count = dofInfo.counts.get(block);
        double[] vals = new double[count];
        for (int d = 0; d < count; d++) {
            vals[d] = curDOFVals.get(offset + d);
        }
        return vals;
    }

    private static double[] concat(double[] a, double[] b) {
        double[] out = new double[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }
}

