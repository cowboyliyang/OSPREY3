package edu.duke.cs.osprey.energy.approximation.branch;

import ai.onnxruntime.*;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.io.File;
import java.nio.LongBuffer;
import java.util.*;

/**
 * Conformation energy calculator using a trained GNN model via ONNX Runtime.
 *
 * Prediction: E_pred(conf) = E_emat(conf) + GNN_residual(conf)
 *
 * The GNN model is loaded from an ONNX file exported by gnn/train.py.
 * Supports RC index mapping when inference confspace differs from training confspace.
 */
public class GNNConfEnergyCalculator implements AutoCloseable {

    private final OrtEnvironment env;
    private final OrtSession session;
    private final EnergyMatrix ematMinimized;
    private final int numPositions;

    // Optional RC mapping: rcMap[pos][inferenceRC] = trainingRC (-1 = unmapped)
    private int[][] rcMap = null;

    public GNNConfEnergyCalculator(File onnxModel, EnergyMatrix ematMinimized, int numPositions) {
        this.ematMinimized = ematMinimized;
        this.numPositions = numPositions;
        try {
            this.env = OrtEnvironment.getEnvironment();
            OrtSession.SessionOptions opts = new OrtSession.SessionOptions();
            opts.setOptimizationLevel(OrtSession.SessionOptions.OptLevel.ALL_OPT);
            // Try GPU first, fall back to CPU
            try {
                opts.addCUDA(0);
                System.out.println("[GNN] ONNX using CUDA provider (GPU 0)");
            } catch (Exception e) {
                System.out.println("[GNN] CUDA provider not available, using CPU: " + e.getMessage());
            }
            this.session = env.createSession(onnxModel.getAbsolutePath(), opts);
        } catch (OrtException e) {
            throw new RuntimeException("Failed to load ONNX model: " + onnxModel, e);
        }
    }

    /**
     * Set RC index mapping for cross-confspace inference.
     * rcMap[pos][inferenceRC] = trainingRC. -1 means unmapped (will be clamped to 0).
     */
    public void setRCMapping(int[][] rcMap) {
        this.rcMap = rcMap;
        int mapped = 0, unmapped = 0;
        for (int p = 0; p < rcMap.length; p++) {
            for (int r = 0; r < rcMap[p].length; r++) {
                if (rcMap[p][r] >= 0) mapped++; else unmapped++;
            }
        }
        System.out.println("[GNN] RC mapping set: " + mapped + " mapped, " + unmapped + " unmapped across " + rcMap.length + " positions");
    }

    /**
     * Build RC mapping from inference confspace to training confspace.
     * Matches by (position resNum, AA template name, rotamer index within that AA).
     * Both confspaces must have the same positions (same resNums in same order).
     */
    public static int[][] buildRCMapping(SimpleConfSpace inferenceCS, SimpleConfSpace trainingCS) {
        int numPos = inferenceCS.positions.size();
        if (numPos != trainingCS.positions.size()) {
            throw new IllegalArgumentException("Position count mismatch: inference=" + numPos
                    + " training=" + trainingCS.positions.size());
        }

        int[][] rcMap = new int[numPos][];

        for (int p = 0; p < numPos; p++) {
            SimpleConfSpace.Position infPos = inferenceCS.positions.get(p);
            SimpleConfSpace.Position trainPos = trainingCS.positions.get(p);

            if (!infPos.resNum.equals(trainPos.resNum)) {
                throw new IllegalArgumentException("Position " + p + " resNum mismatch: "
                        + infPos.resNum + " vs " + trainPos.resNum);
            }

            // Build training lookup: (AA name, rotamer index within AA) → training RC index
            Map<String, Integer> trainAACounter = new LinkedHashMap<>();
            Map<String, List<Integer>> trainAARCs = new LinkedHashMap<>();
            for (int r = 0; r < trainPos.resConfs.size(); r++) {
                SimpleConfSpace.ResidueConf rc = trainPos.resConfs.get(r);
                String aaName = rc.template.name;
                trainAARCs.computeIfAbsent(aaName, k -> new ArrayList<>()).add(r);
            }

            // Map each inference RC
            Map<String, Integer> infAACounter = new LinkedHashMap<>();
            rcMap[p] = new int[infPos.resConfs.size()];
            int mappedCount = 0;
            for (int r = 0; r < infPos.resConfs.size(); r++) {
                SimpleConfSpace.ResidueConf rc = infPos.resConfs.get(r);
                String aaName = rc.template.name;
                int rotIdx = infAACounter.merge(aaName, 0, (a, b) -> a + 1);

                List<Integer> trainRCs = trainAARCs.get(aaName);
                if (trainRCs != null && rotIdx < trainRCs.size()) {
                    rcMap[p][r] = trainRCs.get(rotIdx);
                    mappedCount++;
                } else {
                    rcMap[p][r] = -1; // unmapped
                }
            }
            System.out.println("[GNN] pos " + p + " (" + infPos.resNum + "): "
                    + infPos.resConfs.size() + " inference RCs, "
                    + mappedCount + " mapped to training RCs");
        }

        return rcMap;
    }

    /** Apply RC mapping to a single conf. Returns mapped conf (or original if no mapping). */
    private int[] mapConf(int[] conf) {
        if (rcMap == null) return conf;
        int[] mapped = new int[conf.length];
        for (int p = 0; p < conf.length; p++) {
            int rc = conf[p];
            if (rc >= 0 && rc < rcMap[p].length && rcMap[p][rc] >= 0) {
                mapped[p] = rcMap[p][rc];
            } else {
                mapped[p] = 0; // fallback for unmapped
            }
        }
        return mapped;
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
        // DEBUG: print first few
        if (debugCount <= DEBUG_LIMIT) {
            for (int i = 0; i < Math.min(confs.length, 3); i++) {
                System.out.println("[GNN_DEBUG] calcEnergies: conf=" + java.util.Arrays.toString(confs[i])
                    + " ematE=" + ematEnergies[i] + " residual=" + residuals[i]
                    + " total=" + results[i]);
            }
        }
        return results;
    }

    private double predictResidual(int[] conf) {
        return predictResidualBatch(new int[][]{conf})[0];
    }

    private static int debugCount = 0;
    private static final int DEBUG_LIMIT = 10;

    private double[] predictResidualBatch(int[][] confs) {
        int batchSize = confs.length;
        try {
            // Prepare input tensor: (batch, numPositions) as int64, applying RC mapping if set
            long[] flatConfs = new long[batchSize * numPositions];
            for (int i = 0; i < batchSize; i++) {
                int[] mapped = mapConf(confs[i]);
                for (int j = 0; j < numPositions; j++) {
                    flatConfs[i * numPositions + j] = mapped[j];
                }
            }

            // DEBUG: print first few batches
            if (debugCount < DEBUG_LIMIT) {
                System.out.println("[GNN_DEBUG] batch " + debugCount + ": batchSize=" + batchSize
                    + " numPositions=" + numPositions);
                for (int i = 0; i < Math.min(batchSize, 3); i++) {
                    System.out.println("[GNN_DEBUG]   conf[" + i + "] orig=" + java.util.Arrays.toString(confs[i])
                        + " (len=" + confs[i].length + ")");
                    long[] row = new long[numPositions];
                    System.arraycopy(flatConfs, i * numPositions, row, 0, numPositions);
                    System.out.println("[GNN_DEBUG]   conf[" + i + "] flat=" + java.util.Arrays.toString(row));
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

            // DEBUG: print raw residuals
            if (debugCount < DEBUG_LIMIT) {
                for (int i = 0; i < Math.min(batchSize, 3); i++) {
                    System.out.println("[GNN_DEBUG]   conf[" + i + "] raw_residual=" + residuals[i]);
                }
                debugCount++;
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
