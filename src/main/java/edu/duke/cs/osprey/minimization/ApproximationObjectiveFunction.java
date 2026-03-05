package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction;

import java.util.ArrayList;

/**
 * Objective function backed entirely by an approximator model.
 *
 * This keeps CCD unchanged while replacing forcefield energy calls with
 * confspace-specific surrogate evaluation.
 */
public class ApproximationObjectiveFunction implements ObjectiveFunction {

    private static final long serialVersionUID = -1105233457882286755L;

    public final ParametricMolecule pmol;
    public final ApproximatedObjectiveFunction.Approximator approximator;
    public final DoubleMatrix1D curDOFVals;

    public ApproximationObjectiveFunction(ParametricMolecule pmol,
                                          ApproximatedObjectiveFunction.Approximator approximator) {

        this.pmol = pmol;
        this.approximator = approximator;
        this.curDOFVals = DoubleFactory1D.dense.make(pmol.dofBounds.size());

        if (pmol.dofs.size() != approximator.numDofs()) {
            throw new IllegalArgumentException(String.format(
                    "Approximator DOF count mismatch: pmol=%d approximator=%d",
                    pmol.dofs.size(), approximator.numDofs()
            ));
        }
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
        return approximator.getValForDOF(d, val, curDOFVals);
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
        return approximator.getValue(x);
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
}

