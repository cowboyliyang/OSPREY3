package edu.duke.cs.osprey.energy.approximation.branch;

import ai.onnxruntime.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.io.File;
import java.nio.LongBuffer;
import java.util.HashMap;
import java.util.Map;

/**
 * Subtree-level GNN energy calculator via ONNX Runtime.
 *
 * Unlike GNNConfEnergyCalculator (which predicts residuals for fully-assigned conformations),
 * this predicts the free energy correction ΔF = F_CCD - F_emat for partial assignments (internal nodes).
 *
 * Input:  (fixed_rcs, free_mask) where free positions have assignments = -1
 * Output: ΔF such that logZ_CCD ≈ logZ_emat - ΔF/kT
 *
 * The ONNX model has two inputs:
 *   - "fixed_rcs": (batch, num_pos) int64 — RC indices, 0 for free positions
 *   - "free_mask":  (batch, num_pos) bool  — True for free positions
 */
public class GNNSubtreeEnergyCalculator implements AutoCloseable {

    private final OrtEnvironment env;
    private final OrtSession session;
    private final EnergyMatrix ematMinimized;
    private final int numPositions;

    private int[][] rcMap = null;

    private static int debugCount = 0;
    private static final int DEBUG_LIMIT = 5;

    public GNNSubtreeEnergyCalculator(File onnxModel, EnergyMatrix ematMinimized, int numPositions) {
        this(onnxModel, ematMinimized, numPositions, 0);
    }

    public GNNSubtreeEnergyCalculator(File onnxModel, EnergyMatrix ematMinimized, int numPositions, int gpuDeviceId) {
        this.ematMinimized = ematMinimized;
        this.numPositions = numPositions;
        try {
            this.env = OrtEnvironment.getEnvironment();
            OrtSession.SessionOptions opts = new OrtSession.SessionOptions();
            opts.setOptimizationLevel(OrtSession.SessionOptions.OptLevel.ALL_OPT);
            try {
                opts.addCUDA(gpuDeviceId);
                System.out.println("[GNN-Subtree] ONNX using CUDA provider (GPU " + gpuDeviceId + ")");
            } catch (Exception e) {
                System.out.println("[GNN-Subtree] CUDA not available, using CPU: " + e.getMessage());
            }
            this.session = env.createSession(onnxModel.getAbsolutePath(), opts);
            System.out.println("[GNN-Subtree] Loaded model: " + onnxModel.getAbsolutePath());
        } catch (OrtException e) {
            throw new RuntimeException("Failed to load subtree ONNX model: " + onnxModel, e);
        }
    }

    public void setRCMapping(int[][] rcMap) {
        this.rcMap = rcMap;
    }

    /**
     * Predict free energy correction ΔF for a batch of partial assignments.
     *
     * @param assignments  array of (numPositions,) int arrays; -1 = free position
     * @return ΔF = F_CCD - F_emat (kcal/mol) for each partial assignment
     */
    public double[] predictSubtreeResiduals(int[][] assignments) {
        int batchSize = assignments.length;
        try {
            long[] flatRCs = new long[batchSize * numPositions];
            boolean[][] boolMask = new boolean[batchSize][numPositions];  // true = free

            for (int i = 0; i < batchSize; i++) {
                for (int j = 0; j < numPositions; j++) {
                    int rc = assignments[i][j];
                    if (rc < 0) {
                        // Free position
                        flatRCs[i * numPositions + j] = 0;  // placeholder
                        boolMask[i][j] = true;  // free = True
                    } else {
                        // Fixed position — apply RC mapping if set
                        int mapped = rc;
                        if (rcMap != null && j < rcMap.length && rc < rcMap[j].length && rcMap[j][rc] >= 0) {
                            mapped = rcMap[j][rc];
                        }
                        flatRCs[i * numPositions + j] = mapped;
                        boolMask[i][j] = false;  // fixed = False
                    }
                }
            }

            if (debugCount < DEBUG_LIMIT) {
                System.out.println("[GNN-Subtree] batch: size=" + batchSize);
                for (int i = 0; i < Math.min(batchSize, 2); i++) {
                    StringBuilder sb = new StringBuilder("[GNN-Subtree]   assignments[" + i + "]=[");
                    for (int j = 0; j < numPositions; j++) {
                        if (j > 0) sb.append(",");
                        sb.append(assignments[i][j] < 0 ? "*" : String.valueOf(assignments[i][j]));
                    }
                    sb.append("]");
                    System.out.println(sb.toString());
                }
                debugCount++;
            }

            OnnxTensor rcTensor = OnnxTensor.createTensor(
                    env, LongBuffer.wrap(flatRCs), new long[]{batchSize, numPositions});
            OnnxTensor maskTensor = OnnxTensor.createTensor(env, boolMask);

            Map<String, OnnxTensor> inputs = new HashMap<>();
            inputs.put("fixed_rcs", rcTensor);
            inputs.put("free_mask", maskTensor);

            OrtSession.Result result = session.run(inputs);
            Object rawOutput = result.get(0).getValue();

            double[] residuals = new double[batchSize];
            if (rawOutput instanceof float[][]) {
                float[][] out2d = (float[][]) rawOutput;
                for (int i = 0; i < batchSize; i++) residuals[i] = out2d[i][0];
            } else if (rawOutput instanceof float[]) {
                float[] out1d = (float[]) rawOutput;
                for (int i = 0; i < batchSize; i++) residuals[i] = out1d[i];
            } else {
                throw new RuntimeException("Unexpected output type: " + rawOutput.getClass().getName());
            }

            rcTensor.close();
            maskTensor.close();
            result.close();

            return residuals;
        } catch (OrtException e) {
            throw new RuntimeException("Subtree ONNX inference failed", e);
        }
    }

    /**
     * Compute emat energy for a partial assignment (sum of one-body + fixed-fixed pairwise).
     */
    public double partialEmatEnergy(int[] assignments) {
        double e = 0;
        for (int i = 0; i < assignments.length; i++) {
            if (assignments[i] >= 0) {
                e += ematMinimized.getOneBody(i, assignments[i]);
                for (int j = i + 1; j < assignments.length; j++) {
                    if (assignments[j] >= 0) {
                        e += ematMinimized.getPairwise(i, assignments[i], j, assignments[j]);
                    }
                }
            }
        }
        return e;
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
