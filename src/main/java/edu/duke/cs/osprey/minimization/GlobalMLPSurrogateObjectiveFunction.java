package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DofInfo;
import edu.duke.cs.osprey.energy.approximation.TaskGlobalMLPSurrogate;

import java.util.ArrayList;

/**
 * Objective function backed by a single confspace-global MLP surrogate.
 */
public class GlobalMLPSurrogateObjectiveFunction implements ObjectiveFunction {

    private static final long serialVersionUID = 8305561746088099840L;

    public final ParametricMolecule pmol;
    private final DofInfo dofInfo;
    private final TaskGlobalMLPSurrogate surrogate;
    private final DoubleMatrix1D curDOFVals;
    private final double[] featureScratch;

    public GlobalMLPSurrogateObjectiveFunction(SimpleConfSpace confSpace, RCTuple tuple, TaskGlobalMLPSurrogate surrogate) {
        this.pmol = confSpace.makeMolecule(tuple);
        this.dofInfo = confSpace.makeDofInfo(tuple);
        this.surrogate = surrogate;
        this.curDOFVals = DoubleFactory1D.dense.make(pmol.dofs.size());
        this.featureScratch = new double[surrogate.getInputDim()];
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
        return surrogate.predict(dofInfo.tuple, dofInfo, curDOFVals, featureScratch);
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
