package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DofInfo;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.tools.Progress;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

import static edu.duke.cs.osprey.tools.Log.log;

/**
 * Train or load a single confspace-specific MLP surrogate for full-tuple energies.
 */
public class ConfSpaceSpecificGlobalMLPSurrogateCalculator {

    private static final int FORMAT_VERSION = 1;

    public final ConfEnergyCalculator confEcalc;
    public final MLPEnergyModel.TrainingConfig trainingConfig = new MLPEnergyModel.TrainingConfig();

    private int numTupleSamples = 256;
    private int samplesPerTuple = 32;
    private File cacheFile = null;

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator(ConfEnergyCalculator confEcalc) {
        this.confEcalc = confEcalc;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setNumTupleSamples(int val) {
        numTupleSamples = val;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setSamplesPerTuple(int val) {
        samplesPerTuple = val;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setCacheFile(File val) {
        cacheFile = val;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setHiddenSizes(int h1, int h2) {
        trainingConfig.hidden1 = h1;
        trainingConfig.hidden2 = h2;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setEpochs(int epochs) {
        trainingConfig.epochs = epochs;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setBatchSize(int batchSize) {
        trainingConfig.batchSize = batchSize;
        return this;
    }

    public ConfSpaceSpecificGlobalMLPSurrogateCalculator setLearningRate(double lr) {
        trainingConfig.learningRate = lr;
        return this;
    }

    public TaskGlobalMLPSurrogate calc() {

        TaskGlobalMLPSurrogate layout = new TaskGlobalMLPSurrogate(confEcalc.confSpace, null);

        if (cacheFile != null && cacheFile.exists()) {
            MLPEnergyModel cached = readModel(cacheFile);
            if (cached.numInputs() == layout.getInputDim()) {
                log("read global MLP surrogate model from file: %s", cacheFile.getAbsolutePath());
                return new TaskGlobalMLPSurrogate(confEcalc.confSpace, cached);
            }
            log(
                    "global MLP cache input mismatch (cache=%d expected=%d), retraining: %s",
                    cached.numInputs(),
                    layout.getInputDim(),
                    cacheFile.getAbsolutePath()
            );
        }

        int tuples = Math.max(1, numTupleSamples);
        int perTuple = Math.max(1, samplesPerTuple);
        int totalSamples = tuples * perTuple;

        int inputDim = layout.getInputDim();
        double[][] x = new double[totalSamples][inputDim];
        double[] y = new double[totalSamples];

        long seed = ((long) ConfSpaceSpecificApproximatorCache.fingerprint(confEcalc.confSpace).hashCode() << 32)
                ^ 0x9e3779b97f4a7c15L;
        Random rand = new Random(seed);

        Progress progress = new Progress(tuples);
        log("training 1 global MLP model (%d tuples x %d samples = %d points) ...", tuples, perTuple, totalSamples);

        int idx = 0;
        for (int t = 0; t < tuples; t++) {
            RCTuple tuple = randomFullTuple(confEcalc.confSpace, rand);
            DofInfo dofInfo = confEcalc.confSpace.makeDofInfo(tuple);
            ParametricMolecule pmol = confEcalc.confSpace.makeMolecule(tuple);
            ResidueInteractions inters = confEcalc.makeFragInters(tuple);

            try (EnergyFunction ff = confEcalc.ecalc.makeEnergyFunction(pmol, inters)) {

                MoleculeObjectiveFunction f = new MoleculeObjectiveFunction(pmol, ff);
                DoubleMatrix1D[] bounds = f.getConstraints();
                int numDofs = f.getNumDOFs();
                DoubleMatrix1D xi = DoubleFactory1D.dense.make(numDofs);
                double[] feat = new double[inputDim];

                for (int s = 0; s < perTuple; s++) {
                    if (s == 0) {
                        xi.assign(f.getDOFsCenter());
                    } else {
                        for (int d = 0; d < numDofs; d++) {
                            double lo = bounds[0].get(d);
                            double hi = bounds[1].get(d);
                            xi.set(d, lo + rand.nextDouble() * (hi - lo));
                        }
                    }

                    layout.encodeInto(tuple, dofInfo, xi, feat);
                    x[idx] = Arrays.copyOf(feat, feat.length);
                    y[idx] = f.getValue(xi);
                    idx++;
                }
            }

            progress.incrementProgress();
        }

        if (idx < x.length) {
            x = Arrays.copyOf(x, idx);
            y = Arrays.copyOf(y, idx);
        }

        double[] inputMin = new double[inputDim];
        double[] inputMax = new double[inputDim];
        Arrays.fill(inputMin, Double.POSITIVE_INFINITY);
        Arrays.fill(inputMax, Double.NEGATIVE_INFINITY);
        for (double[] sample : x) {
            for (int d = 0; d < inputDim; d++) {
                inputMin[d] = Math.min(inputMin[d], sample[d]);
                inputMax[d] = Math.max(inputMax[d], sample[d]);
            }
        }

        for (int d = 0; d < layout.getOneHotDim(); d++) {
            inputMin[d] = 0.0;
            inputMax[d] = 1.0;
        }
        for (int d = 0; d < inputDim; d++) {
            if (!Double.isFinite(inputMin[d]) || !Double.isFinite(inputMax[d])) {
                inputMin[d] = 0.0;
                inputMax[d] = 0.0;
            }
        }

        MLPEnergyModel model = MLPEnergyModel.train(
                x,
                y,
                inputMin,
                inputMax,
                trainingConfig,
                seed ^ 0x5deece66dL
        );

        TaskGlobalMLPSurrogate trained = new TaskGlobalMLPSurrogate(confEcalc.confSpace, model);

        if (cacheFile != null) {
            writeModel(cacheFile, model);
            log("wrote global MLP surrogate model to file: %s", cacheFile.getAbsolutePath());
        }

        return trained;
    }

    private static RCTuple randomFullTuple(SimpleConfSpace confSpace, Random rand) {
        int[] conf = new int[confSpace.positions.size()];
        for (SimpleConfSpace.Position pos : confSpace.positions) {
            conf[pos.index] = rand.nextInt(pos.resConfs.size());
        }
        return new RCTuple(conf);
    }

    private static void writeModel(File file, MLPEnergyModel model) {
        File parent = file.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new RuntimeException("Failed to create global MLP cache directory: " + parent.getAbsolutePath());
        }

        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)))) {
            out.writeInt(FORMAT_VERSION);
            model.writeTo(out);
        } catch (IOException ex) {
            throw new RuntimeException("can't write global MLP model to file: " + file.getAbsolutePath(), ex);
        }
    }

    private static MLPEnergyModel readModel(File file) {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(file)))) {
            int version = in.readInt();
            if (version != FORMAT_VERSION) {
                throw new RuntimeException("unsupported global MLP surrogate format version: " + version);
            }
            return MLPEnergyModel.readFrom(in);
        } catch (IOException ex) {
            throw new RuntimeException("can't read global MLP model from file: " + file.getAbsolutePath(), ex);
        }
    }
}
