package edu.duke.cs.osprey.energy.approximation.branch;

import ai.onnxruntime.*;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.io.File;
import java.nio.LongBuffer;
import java.util.Collections;

/**
 * Conformation energy calculator using a trained GNN model via ONNX Runtime.
 *
 * Prediction: E_pred(conf) = E_emat(conf) + GNN_residual(conf)
 *
 * The GNN model is loaded from an ONNX file exported by gnn/train.py.
 */
public class GNNConfEnergyCalculator implements AutoCloseable {

    private final OrtEnvironment env;
    private final OrtSession session;
    private final EnergyMatrix ematMinimized;
    private final int numPositions;

    public GNNConfEnergyCalculator(File onnxModel, EnergyMatrix ematMinimized, int numPositions) {
        this.ematMinimized = ematMinimized;
        this.numPositions = numPositions;
        try {
            this.env = OrtEnvironment.getEnvironment();
            OrtSession.SessionOptions opts = new OrtSession.SessionOptions();
            opts.setOptimizationLevel(OrtSession.SessionOptions.OptLevel.ALL_OPT);
            this.session = env.createSession(onnxModel.getAbsolutePath(), opts);
        } catch (OrtException e) {
            throw new RuntimeException("Failed to load ONNX model: " + onnxModel, e);
        }
    }

    /**
     * Predict CCD-quality energy for a single conformation.
     */
    public double calcEnergy(int[] conf) {
        double ematEnergy = ematMinimized.confE(conf);
        double residual = predictResidual(conf);
        return ematEnergy + residual;
    }

    /**
     * Predict CCD-quality energy for a single conformation given as RCTuple.
     */
    public double calcEnergy(RCTuple tuple) {
        int[] conf = new int[numPositions];
        for (int i = 0; i < tuple.pos.size(); i++) {
            conf[tuple.pos.get(i)] = tuple.RCs.get(i);
        }
        return calcEnergy(conf);
    }

    /**
     * Batch predict energies for multiple conformations.
     */
    public double[] calcEnergies(int[][] confs) {
        double[] ematEnergies = new double[confs.length];
        for (int i = 0; i < confs.length; i++) {
            ematEnergies[i] = ematMinimized.confE(confs[i]);
        }
        double[] residuals = predictResidualBatch(confs);
        double[] results = new double[confs.length];
        for (int i = 0; i < confs.length; i++) {
            results[i] = ematEnergies[i] + residuals[i];
        }
        return results;
    }

    private double predictResidual(int[] conf) {
        return predictResidualBatch(new int[][]{conf})[0];
    }

    private double[] predictResidualBatch(int[][] confs) {
        int batchSize = confs.length;
        try {
            // Prepare input tensor: (batch, numPositions) as int64
            long[] flatConfs = new long[batchSize * numPositions];
            for (int i = 0; i < batchSize; i++) {
                for (int j = 0; j < numPositions; j++) {
                    flatConfs[i * numPositions + j] = confs[i][j];
                }
            }
            OnnxTensor inputTensor = OnnxTensor.createTensor(
                    env,
                    LongBuffer.wrap(flatConfs),
                    new long[]{batchSize, numPositions}
            );

            // Run inference
            OrtSession.Result result = session.run(
                    Collections.singletonMap("confs", inputTensor));
            Object rawOutput = result.get(0).getValue();

            double[] residuals = new double[batchSize];
            if (rawOutput instanceof float[][]) {
                float[][] output2d = (float[][]) rawOutput;
                for (int i = 0; i < batchSize; i++) {
                    residuals[i] = output2d[i][0];
                }
            } else if (rawOutput instanceof float[]) {
                float[] output1d = (float[]) rawOutput;
                for (int i = 0; i < batchSize; i++) {
                    residuals[i] = output1d[i];
                }
            } else {
                throw new RuntimeException("Unexpected ONNX output type: " + rawOutput.getClass().getName());
            }

            inputTensor.close();
            result.close();

            return residuals;
        } catch (OrtException e) {
            throw new RuntimeException("ONNX inference failed", e);
        }
    }

    @Override
    public void close() {
        try {
            session.close();
        } catch (OrtException e) {
            // ignore
        }
    }
}
