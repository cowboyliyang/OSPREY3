package edu.duke.cs.osprey.energy.approximation.branch;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * One shared linear model over hashed branch features:
 *   residual(conf) = w dot psi(conf)
 */
public class BranchSharedRidgeModel {

    private static final int FORMAT_VERSION = 1;

    public static class TrainingData {
        public final int hashDim;
        public final int[][] indices;
        public final double[][] values;
        public final double[] targets;

        public TrainingData(int hashDim, int[][] indices, double[][] values, double[] targets) {
            this.hashDim = hashDim;
            this.indices = indices;
            this.values = values;
            this.targets = targets;

            if (indices.length != values.length || values.length != targets.length) {
                throw new IllegalArgumentException("indices/values/targets length mismatch");
            }
        }

        public int size() {
            return targets.length;
        }

        public static TrainingData fromLists(
                int hashDim,
                List<BranchResidualFeatureEncoder.SparseFeatures> features,
                List<Double> targets
        ) {
            if (features.size() != targets.size()) {
                throw new IllegalArgumentException("features/targets size mismatch");
            }
            int n = features.size();
            int[][] idx = new int[n][];
            double[][] val = new double[n][];
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                BranchResidualFeatureEncoder.SparseFeatures sf = features.get(i);
                idx[i] = sf.indices;
                val[i] = sf.values;
                y[i] = targets.get(i);
            }
            return new TrainingData(hashDim, idx, val, y);
        }
    }

    public final int hashDim;
    public final double[] weights;

    public BranchSharedRidgeModel(int hashDim, double[] weights) {
        this.hashDim = hashDim;
        this.weights = weights;
        if (weights.length != hashDim) {
            throw new IllegalArgumentException(String.format(
                    "weight dimension mismatch: weights=%d hashDim=%d", weights.length, hashDim
            ));
        }
    }

    public double predict(BranchResidualFeatureEncoder.SparseFeatures features) {
        return predict(features.indices, features.values);
    }

    public double predict(int[] indices, double[] values) {
        double out = 0.0;
        for (int i = 0; i < indices.length; i++) {
            int idx = indices[i];
            if (idx >= 0 && idx < weights.length) {
                out += values[i] * weights[idx];
            }
        }
        return out;
    }

    public static BranchSharedRidgeModel fit(
            TrainingData data,
            double lambda,
            int maxIterations,
            double tolerance
    ) {
        return fitWithPrior(data, lambda, null, maxIterations, tolerance);
    }

    /**
     * Solve:
     *   argmin ||Xw - y||^2 + lambda ||w - prior||^2
     */
    public static BranchSharedRidgeModel fitWithPrior(
            TrainingData data,
            double lambda,
            BranchSharedRidgeModel prior,
            int maxIterations,
            double tolerance
    ) {

        if (data == null) {
            throw new IllegalArgumentException("data must not be null");
        }
        if (data.hashDim <= 0) {
            throw new IllegalArgumentException("data.hashDim must be > 0");
        }

        double[] priorWeights = null;
        if (prior != null) {
            if (prior.hashDim != data.hashDim) {
                throw new IllegalArgumentException("prior hashDim does not match training data");
            }
            priorWeights = prior.weights;
        }

        double[] w = solveByConjugateGradient(
                data,
                Math.max(0.0, lambda),
                priorWeights,
                Math.max(1, maxIterations),
                Math.max(1e-14, tolerance)
        );

        return new BranchSharedRidgeModel(data.hashDim, w);
    }

    private static double[] solveByConjugateGradient(
            TrainingData data,
            double lambda,
            double[] prior,
            int maxIterations,
            double tolerance
    ) {
        int dim = data.hashDim;

        double[] b = new double[dim];
        if (prior != null && lambda > 0.0) {
            for (int d = 0; d < dim; d++) {
                b[d] = lambda * prior[d];
            }
        }

        for (int i = 0; i < data.size(); i++) {
            int[] idx = data.indices[i];
            double[] val = data.values[i];
            double y = data.targets[i];
            for (int k = 0; k < idx.length; k++) {
                int d = idx[k];
                if (d >= 0 && d < dim) {
                    b[d] += val[k] * y;
                }
            }
        }

        double[] x = (prior != null) ? Arrays.copyOf(prior, prior.length) : new double[dim];
        double[] Ax = new double[dim];
        applyNormalEq(data, lambda, x, Ax);

        double[] r = new double[dim];
        for (int d = 0; d < dim; d++) {
            r[d] = b[d] - Ax[d];
        }

        double rNorm = l2Norm(r);
        if (rNorm <= tolerance) {
            return x;
        }

        double[] p = Arrays.copyOf(r, r.length);
        double[] Ap = new double[dim];
        double rsOld = dot(r, r);

        for (int it = 0; it < maxIterations; it++) {

            Arrays.fill(Ap, 0.0);
            applyNormalEq(data, lambda, p, Ap);

            double denom = dot(p, Ap);
            if (Math.abs(denom) < 1e-20) {
                break;
            }

            double alpha = rsOld / denom;
            for (int d = 0; d < dim; d++) {
                x[d] += alpha * p[d];
                r[d] -= alpha * Ap[d];
            }

            double rsNew = dot(r, r);
            if (Math.sqrt(rsNew) <= tolerance) {
                break;
            }

            double beta = rsNew / rsOld;
            for (int d = 0; d < dim; d++) {
                p[d] = r[d] + beta * p[d];
            }
            rsOld = rsNew;
        }

        return x;
    }

    /**
     * out = (X^T X + lambda I) * in
     */
    private static void applyNormalEq(TrainingData data, double lambda, double[] in, double[] out) {
        int dim = data.hashDim;

        if (lambda > 0.0) {
            for (int d = 0; d < dim; d++) {
                out[d] += lambda * in[d];
            }
        }

        for (int i = 0; i < data.size(); i++) {
            int[] idx = data.indices[i];
            double[] val = data.values[i];

            double xdotv = 0.0;
            for (int k = 0; k < idx.length; k++) {
                int d = idx[k];
                if (d >= 0 && d < dim) {
                    xdotv += val[k] * in[d];
                }
            }

            for (int k = 0; k < idx.length; k++) {
                int d = idx[k];
                if (d >= 0 && d < dim) {
                    out[d] += val[k] * xdotv;
                }
            }
        }
    }

    private static double dot(double[] a, double[] b) {
        double out = 0.0;
        for (int i = 0; i < a.length; i++) {
            out += a[i] * b[i];
        }
        return out;
    }

    private static double l2Norm(double[] x) {
        return Math.sqrt(dot(x, x));
    }

    public void writeTo(DataOutputStream out) throws IOException {
        out.writeInt(FORMAT_VERSION);
        out.writeInt(hashDim);
        out.writeInt(weights.length);
        for (double weight : weights) {
            out.writeDouble(weight);
        }
    }

    public static BranchSharedRidgeModel readFrom(DataInputStream in) throws IOException {
        int version = in.readInt();
        if (version != FORMAT_VERSION) {
            throw new IOException("Unsupported BranchSharedRidgeModel format version: " + version);
        }

        int hashDim = in.readInt();
        int len = in.readInt();
        if (len != hashDim) {
            throw new IOException(String.format(
                    "Invalid serialized ridge model: hashDim=%d weightLen=%d", hashDim, len
            ));
        }

        double[] w = new double[len];
        for (int i = 0; i < len; i++) {
            w[i] = in.readDouble();
        }
        return new BranchSharedRidgeModel(hashDim, w);
    }

    public void writeToFile(File file) {
        File parent = file.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new RuntimeException("Failed to create ridge model directory: " + parent.getAbsolutePath());
        }

        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)))) {
            writeTo(out);
        } catch (IOException ex) {
            throw new RuntimeException("Failed to write ridge model: " + file.getAbsolutePath(), ex);
        }
    }

    public static BranchSharedRidgeModel readFromFile(File file) {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(file)))) {
            return readFrom(in);
        } catch (IOException ex) {
            throw new RuntimeException("Failed to read ridge model: " + file.getAbsolutePath(), ex);
        }
    }
}
