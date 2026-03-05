package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.tools.IOable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

/**
 * Per-RC one-body and pairwise MLP energy models for a conf space.
 */
public class MLPSurrogateMatrix implements IOable {

    private static final int FORMAT_VERSION = 1;

    public final SimpleConfSpace confSpace;
    private final TupleMatrixGeneric<MLPEnergyModel> models;

    public MLPSurrogateMatrix(SimpleConfSpace confSpace) {
        this.confSpace = confSpace;
        this.models = new TupleMatrixGeneric<>(confSpace);
    }

    public MLPEnergyModel getOneBody(int pos, int rc) {
        return models.getOneBody(pos, rc);
    }

    public void setOneBody(int pos, int rc, MLPEnergyModel model) {
        models.setOneBody(pos, rc, model);
    }

    public MLPEnergyModel getPair(int pos1, int rc1, int pos2, int rc2) {
        return models.getPairwise(pos1, rc1, pos2, rc2);
    }

    public void setPair(int pos1, int rc1, int pos2, int rc2, MLPEnergyModel model) {
        models.setPairwise(pos1, rc1, pos2, rc2, model);
    }

    @Override
    public void writeTo(DataOutput out) throws IOException {
        out.writeInt(FORMAT_VERSION);

        // singles
        for (SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                writeModel(getOneBody(pos1.index, rc1.index), out);
            }
        }

        // pairs
        for (SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
                    for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {
                        writeModel(getPair(pos1.index, rc1.index, pos2.index, rc2.index), out);
                    }
                }
            }
        }
    }

    @Override
    public void readFrom(DataInput in) throws IOException {
        int version = in.readInt();
        if (version != FORMAT_VERSION) {
            throw new IOException("unsupported MLPSurrogateMatrix format version: " + version);
        }

        // singles
        for (SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                setOneBody(pos1.index, rc1.index, readModel(in));
            }
        }

        // pairs
        for (SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
                    for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {
                        setPair(
                                pos1.index, rc1.index, pos2.index, rc2.index,
                                readModel(in)
                        );
                    }
                }
            }
        }
    }

    private static void writeModel(MLPEnergyModel model, DataOutput out) throws IOException {
        out.writeBoolean(model != null);
        if (model != null) {
            model.writeTo(out);
        }
    }

    private static MLPEnergyModel readModel(DataInput in) throws IOException {
        boolean present = in.readBoolean();
        if (!present) {
            return null;
        }
        return MLPEnergyModel.readFrom(in);
    }
}
