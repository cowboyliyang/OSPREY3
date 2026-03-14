package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DofInfo;

import java.util.Arrays;

/**
 * A single MLP surrogate for an entire conf space task.
 *
 * Input features are concatenated as:
 * - one-hot RC identity for each position
 * - padded per-position DOF values (up to max DOF count seen at that position)
 */
public class TaskGlobalMLPSurrogate {

    public final SimpleConfSpace confSpace;
    public final MLPEnergyModel model;

    private final int[] rcOffsetsByPos;
    private final int[] dofOffsetsByPos;
    private final int[] maxDofsByPos;
    private final int oneHotDim;
    private final int inputDim;

    public TaskGlobalMLPSurrogate(SimpleConfSpace confSpace, MLPEnergyModel model) {
        this.confSpace = confSpace;
        this.model = model;

        int nPos = confSpace.positions.size();
        this.rcOffsetsByPos = new int[nPos];
        this.dofOffsetsByPos = new int[nPos];
        this.maxDofsByPos = new int[nPos];

        int cursor = 0;
        for (SimpleConfSpace.Position pos : confSpace.positions) {
            rcOffsetsByPos[pos.index] = cursor;
            cursor += pos.resConfs.size();
        }
        this.oneHotDim = cursor;

        for (SimpleConfSpace.Position pos : confSpace.positions) {
            int maxDofs = 0;
            for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
                maxDofs = Math.max(maxDofs, rc.dofBounds.size());
            }
            maxDofsByPos[pos.index] = maxDofs;
            dofOffsetsByPos[pos.index] = cursor;
            cursor += maxDofs;
        }
        this.inputDim = cursor;

        if (model != null && model.numInputs() != inputDim) {
            throw new IllegalArgumentException(String.format(
                    "global MLP input mismatch: model=%d expected=%d",
                    model.numInputs(), inputDim
            ));
        }
    }

    public int getInputDim() {
        return inputDim;
    }

    public int getOneHotDim() {
        return oneHotDim;
    }

    public void encodeInto(RCTuple tuple, DofInfo dofInfo, DoubleMatrix1D x, double[] out) {
        if (out.length != inputDim) {
            throw new IllegalArgumentException(String.format("feature buffer has %d dims, expected %d", out.length, inputDim));
        }

        Arrays.fill(out, 0.0);

        for (int i = 0; i < tuple.size(); i++) {
            int pos = tuple.pos.get(i);
            int rc = tuple.RCs.get(i);
            out[rcOffsetsByPos[pos] + rc] = 1.0;
        }

        for (int b = 0; b < dofInfo.size(); b++) {
            SimpleConfSpace.Position pos = dofInfo.positions.get(b);
            if (pos == null) {
                continue;
            }

            int posIdx = pos.index;
            int srcOffset = dofInfo.offsets.get(b);
            int count = Math.min(dofInfo.counts.get(b), maxDofsByPos[posIdx]);
            int dstOffset = dofOffsetsByPos[posIdx];

            for (int d = 0; d < count; d++) {
                out[dstOffset + d] = x.get(srcOffset + d);
            }
        }
    }

    public double predict(RCTuple tuple, DofInfo dofInfo, DoubleMatrix1D x, double[] scratchFeatures) {
        if (model == null) {
            throw new IllegalStateException("global MLP model is not set");
        }
        encodeInto(tuple, dofInfo, x, scratchFeatures);
        return model.predict(scratchFeatures);
    }
}
