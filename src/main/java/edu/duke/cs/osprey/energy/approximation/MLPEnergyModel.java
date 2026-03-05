package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleMatrix1D;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

/**
 * Lightweight fully-connected MLP regressor for forcefield energy approximation.
 *
 * Architecture: input -> tanh(hidden1) -> tanh(hidden2) -> linear output
 */
public class MLPEnergyModel {

    public static class TrainingConfig {
        public int hidden1 = 32;
        public int hidden2 = 32;
        public int epochs = 300;
        public int batchSize = 64;
        public double learningRate = 1e-3;
        public double l2 = 1e-8;
        public double validationFraction = 0.2;
    }

    public final int inputDim;
    public final int hidden1;
    public final int hidden2;

    public final double[] inputMin;
    public final double[] inputMax;

    public final double targetMean;
    public final double targetStd;

    public final double[] w1; // [hidden1 * inputDim]
    public final double[] b1; // [hidden1]
    public final double[] w2; // [hidden2 * hidden1]
    public final double[] b2; // [hidden2]
    public final double[] w3; // [hidden2]
    public final double b3;

    public final double maxAbsError;

    private MLPEnergyModel(
            int inputDim,
            int hidden1,
            int hidden2,
            double[] inputMin,
            double[] inputMax,
            double targetMean,
            double targetStd,
            double[] w1,
            double[] b1,
            double[] w2,
            double[] b2,
            double[] w3,
            double b3,
            double maxAbsError
    ) {
        this.inputDim = inputDim;
        this.hidden1 = hidden1;
        this.hidden2 = hidden2;
        this.inputMin = inputMin;
        this.inputMax = inputMax;
        this.targetMean = targetMean;
        this.targetStd = targetStd;
        this.w1 = w1;
        this.b1 = b1;
        this.w2 = w2;
        this.b2 = b2;
        this.w3 = w3;
        this.b3 = b3;
        this.maxAbsError = maxAbsError;
    }

    public int numInputs() {
        return inputDim;
    }

    public double error() {
        return maxAbsError;
    }

    public double predict(DoubleMatrix1D x) {
        double[] in = new double[inputDim];
        for (int d = 0; d < inputDim; d++) {
            in[d] = x.get(d);
        }
        return predict(in);
    }

    public double predict(double[] xRaw) {
        if (xRaw.length != inputDim) {
            throw new IllegalArgumentException(String.format("x has %d dims, expected %d", xRaw.length, inputDim));
        }

        double[] x = new double[inputDim];
        normalizeInput(xRaw, x);

        double[] a1 = new double[hidden1];
        for (int i = 0; i < hidden1; i++) {
            double z = b1[i];
            int row = i * inputDim;
            for (int d = 0; d < inputDim; d++) {
                z += w1[row + d] * x[d];
            }
            a1[i] = Math.tanh(z);
        }

        double[] a2 = new double[hidden2];
        for (int j = 0; j < hidden2; j++) {
            double z = b2[j];
            int row = j * hidden1;
            for (int i = 0; i < hidden1; i++) {
                z += w2[row + i] * a1[i];
            }
            a2[j] = Math.tanh(z);
        }

        double yNorm = b3;
        for (int j = 0; j < hidden2; j++) {
            yNorm += w3[j] * a2[j];
        }
        return yNorm * targetStd + targetMean;
    }

    public void writeTo(DataOutput out) throws IOException {
        out.writeInt(inputDim);
        out.writeInt(hidden1);
        out.writeInt(hidden2);
        writeArray(out, inputMin);
        writeArray(out, inputMax);
        out.writeDouble(targetMean);
        out.writeDouble(targetStd);
        writeArray(out, w1);
        writeArray(out, b1);
        writeArray(out, w2);
        writeArray(out, b2);
        writeArray(out, w3);
        out.writeDouble(b3);
        out.writeDouble(maxAbsError);
    }

    public static MLPEnergyModel readFrom(DataInput in) throws IOException {
        int inputDim = in.readInt();
        int hidden1 = in.readInt();
        int hidden2 = in.readInt();
        double[] inputMin = readArray(in);
        double[] inputMax = readArray(in);
        double targetMean = in.readDouble();
        double targetStd = in.readDouble();
        double[] w1 = readArray(in);
        double[] b1 = readArray(in);
        double[] w2 = readArray(in);
        double[] b2 = readArray(in);
        double[] w3 = readArray(in);
        double b3 = in.readDouble();
        double maxAbsError = in.readDouble();
        return new MLPEnergyModel(
                inputDim, hidden1, hidden2, inputMin, inputMax, targetMean, targetStd,
                w1, b1, w2, b2, w3, b3, maxAbsError
        );
    }

    public static int numParams(int inputDim, int hidden1, int hidden2) {
        return hidden1 * inputDim + hidden1
                + hidden2 * hidden1 + hidden2
                + hidden2 + 1;
    }

    public static MLPEnergyModel train(
            double[][] xRaw,
            double[] yRaw,
            double[] inputMin,
            double[] inputMax,
            TrainingConfig cfg,
            long seed
    ) {
        if (xRaw.length != yRaw.length) {
            throw new IllegalArgumentException("x and y sample counts do not match");
        }
        if (xRaw.length == 0) {
            throw new IllegalArgumentException("no samples");
        }
        int n = xRaw.length;
        int inputDim = xRaw[0].length;

        double targetMean = mean(yRaw);
        double targetStd = std(yRaw, targetMean);
        if (targetStd < 1e-12) {
            targetStd = 1.0;
        }

        double[] yNorm = new double[n];
        for (int i = 0; i < n; i++) {
            yNorm[i] = (yRaw[i] - targetMean) / targetStd;
        }

        double[][] xNorm = new double[n][inputDim];
        for (int i = 0; i < n; i++) {
            normalizeInput(xRaw[i], xNorm[i], inputMin, inputMax);
        }

        if (n == 1 || inputDim == 0) {
            return constant(inputDim, inputMin, inputMax, targetMean, maxAbsError(yRaw, targetMean));
        }

        int[] order = new int[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        Random splitRand = new Random(seed ^ 0x5deece66dL);
        shuffle(order, splitRand);
        int nVal = Math.max(1, (int) Math.round(n * cfg.validationFraction));
        if (nVal >= n) {
            nVal = n - 1;
        }
        int nTrain = n - nVal;
        int[] trainIdx = Arrays.copyOfRange(order, 0, nTrain);
        int[] valIdx = Arrays.copyOfRange(order, nTrain, n);

        int h1 = Math.max(1, cfg.hidden1);
        int h2 = Math.max(1, cfg.hidden2);
        Random rand = new Random(seed);

        double[] w1 = new double[h1 * inputDim];
        double[] b1 = new double[h1];
        double[] w2 = new double[h2 * h1];
        double[] b2 = new double[h2];
        double[] w3 = new double[h2];
        double b3 = 0.0;

        initWeights(w1, inputDim, rand);
        initWeights(w2, h1, rand);
        initWeights(w3, h2, rand);

        double[] mw1 = new double[w1.length], vw1 = new double[w1.length];
        double[] mb1 = new double[b1.length], vb1 = new double[b1.length];
        double[] mw2 = new double[w2.length], vw2 = new double[w2.length];
        double[] mb2 = new double[b2.length], vb2 = new double[b2.length];
        double[] mw3 = new double[w3.length], vw3 = new double[w3.length];
        double mb3 = 0.0, vb3 = 0.0;

        double[] gw1 = new double[w1.length];
        double[] gb1 = new double[b1.length];
        double[] gw2 = new double[w2.length];
        double[] gb2 = new double[b2.length];
        double[] gw3 = new double[w3.length];
        double gb3;

        double[] a1 = new double[h1];
        double[] a2 = new double[h2];
        double[] dz2 = new double[h2];
        double[] da1 = new double[h1];
        double[] dz1 = new double[h1];

        int batchSize = Math.max(1, cfg.batchSize);
        double bestValLoss = Double.POSITIVE_INFINITY;
        double[] bestW1 = Arrays.copyOf(w1, w1.length);
        double[] bestB1 = Arrays.copyOf(b1, b1.length);
        double[] bestW2 = Arrays.copyOf(w2, w2.length);
        double[] bestB2 = Arrays.copyOf(b2, b2.length);
        double[] bestW3 = Arrays.copyOf(w3, w3.length);
        double bestB3 = b3;

        long t = 0;
        Random epochRand = new Random(seed ^ 0x9e3779b97f4a7c15L);
        for (int epoch = 0; epoch < cfg.epochs; epoch++) {
            shuffle(trainIdx, epochRand);

            for (int start = 0; start < nTrain; start += batchSize) {
                int end = Math.min(nTrain, start + batchSize);
                int batchN = end - start;

                Arrays.fill(gw1, 0.0);
                Arrays.fill(gb1, 0.0);
                Arrays.fill(gw2, 0.0);
                Arrays.fill(gb2, 0.0);
                Arrays.fill(gw3, 0.0);
                gb3 = 0.0;

                for (int ii = start; ii < end; ii++) {
                    int idx = trainIdx[ii];
                    double[] x = xNorm[idx];
                    double y = yNorm[idx];

                    for (int i = 0; i < h1; i++) {
                        double z = b1[i];
                        int row = i * inputDim;
                        for (int d = 0; d < inputDim; d++) {
                            z += w1[row + d] * x[d];
                        }
                        a1[i] = Math.tanh(z);
                    }

                    for (int j = 0; j < h2; j++) {
                        double z = b2[j];
                        int row = j * h1;
                        for (int i = 0; i < h1; i++) {
                            z += w2[row + i] * a1[i];
                        }
                        a2[j] = Math.tanh(z);
                    }

                    double yhat = b3;
                    for (int j = 0; j < h2; j++) {
                        yhat += w3[j] * a2[j];
                    }
                    double dy = 2.0 * (yhat - y);

                    for (int j = 0; j < h2; j++) {
                        gw3[j] += dy * a2[j];
                    }
                    gb3 += dy;

                    Arrays.fill(da1, 0.0);
                    for (int j = 0; j < h2; j++) {
                        dz2[j] = dy * w3[j] * (1.0 - a2[j] * a2[j]);
                        gb2[j] += dz2[j];
                        int row = j * h1;
                        for (int i = 0; i < h1; i++) {
                            gw2[row + i] += dz2[j] * a1[i];
                            da1[i] += dz2[j] * w2[row + i];
                        }
                    }

                    for (int i = 0; i < h1; i++) {
                        dz1[i] = da1[i] * (1.0 - a1[i] * a1[i]);
                        gb1[i] += dz1[i];
                        int row = i * inputDim;
                        for (int d = 0; d < inputDim; d++) {
                            gw1[row + d] += dz1[i] * x[d];
                        }
                    }
                }

                double invBatch = 1.0 / batchN;
                scaleAndRegularize(gw1, w1, invBatch, cfg.l2);
                scaleAndRegularize(gb1, b1, invBatch, 0.0);
                scaleAndRegularize(gw2, w2, invBatch, cfg.l2);
                scaleAndRegularize(gb2, b2, invBatch, 0.0);
                scaleAndRegularize(gw3, w3, invBatch, cfg.l2);
                gb3 *= invBatch;

                t++;
                adamUpdate(w1, gw1, mw1, vw1, cfg.learningRate, t);
                adamUpdate(b1, gb1, mb1, vb1, cfg.learningRate, t);
                adamUpdate(w2, gw2, mw2, vw2, cfg.learningRate, t);
                adamUpdate(b2, gb2, mb2, vb2, cfg.learningRate, t);
                adamUpdate(w3, gw3, mw3, vw3, cfg.learningRate, t);
                double[] mb3Arr = new double[]{mb3};
                double[] vb3Arr = new double[]{vb3};
                double[] gb3Arr = new double[]{gb3};
                double[] b3Arr = new double[]{b3};
                adamUpdate(b3Arr, gb3Arr, mb3Arr, vb3Arr, cfg.learningRate, t);
                b3 = b3Arr[0];
                mb3 = mb3Arr[0];
                vb3 = vb3Arr[0];
            }

            double valLoss = computeLoss(xNorm, yNorm, valIdx, w1, b1, w2, b2, w3, b3, inputDim, h1, h2);
            if (valLoss < bestValLoss) {
                bestValLoss = valLoss;
                System.arraycopy(w1, 0, bestW1, 0, w1.length);
                System.arraycopy(b1, 0, bestB1, 0, b1.length);
                System.arraycopy(w2, 0, bestW2, 0, w2.length);
                System.arraycopy(b2, 0, bestB2, 0, b2.length);
                System.arraycopy(w3, 0, bestW3, 0, w3.length);
                bestB3 = b3;
            }
        }

        MLPEnergyModel model = new MLPEnergyModel(
                inputDim, h1, h2,
                Arrays.copyOf(inputMin, inputMin.length),
                Arrays.copyOf(inputMax, inputMax.length),
                targetMean, targetStd,
                bestW1, bestB1, bestW2, bestB2, bestW3, bestB3,
                0.0
        );
        double err = model.maxAbsErrorOnDataset(xRaw, yRaw);
        return new MLPEnergyModel(
                model.inputDim, model.hidden1, model.hidden2,
                model.inputMin, model.inputMax, model.targetMean, model.targetStd,
                model.w1, model.b1, model.w2, model.b2, model.w3, model.b3, err
        );
    }

    private static MLPEnergyModel constant(int inputDim, double[] inputMin, double[] inputMax, double value, double err) {
        int h1 = 1;
        int h2 = 1;
        double[] w1 = new double[h1 * inputDim];
        double[] b1 = new double[h1];
        double[] w2 = new double[h2 * h1];
        double[] b2 = new double[h2];
        double[] w3 = new double[h2];
        return new MLPEnergyModel(
                inputDim, h1, h2,
                Arrays.copyOf(inputMin, inputMin.length),
                Arrays.copyOf(inputMax, inputMax.length),
                value, 1.0, w1, b1, w2, b2, w3, 0.0, err
        );
    }

    private double maxAbsErrorOnDataset(double[][] xRaw, double[] yRaw) {
        double err = 0.0;
        for (int i = 0; i < xRaw.length; i++) {
            double pred = predict(xRaw[i]);
            err = Math.max(err, Math.abs(pred - yRaw[i]));
        }
        return err;
    }

    private static void normalizeInput(double[] raw, double[] out, double[] min, double[] max) {
        for (int d = 0; d < raw.length; d++) {
            double lo = min[d];
            double hi = max[d];
            double range = hi - lo;
            if (Math.abs(range) < 1e-12) {
                out[d] = 0.0;
            } else {
                out[d] = (2.0 * (raw[d] - lo) / range) - 1.0;
            }
        }
    }

    private void normalizeInput(double[] raw, double[] out) {
        normalizeInput(raw, out, inputMin, inputMax);
    }

    private static double mean(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s / x.length;
    }

    private static double std(double[] x, double mean) {
        if (x.length < 2) {
            return 0.0;
        }
        double s = 0.0;
        for (double v : x) {
            double d = v - mean;
            s += d * d;
        }
        return Math.sqrt(s / x.length);
    }

    private static double maxAbsError(double[] y, double value) {
        double e = 0.0;
        for (double v : y) {
            e = Math.max(e, Math.abs(v - value));
        }
        return e;
    }

    private static void initWeights(double[] w, int fanIn, Random rand) {
        double scale = Math.sqrt(2.0 / Math.max(1, fanIn));
        for (int i = 0; i < w.length; i++) {
            w[i] = rand.nextGaussian() * scale;
        }
    }

    private static void scaleAndRegularize(double[] grad, double[] param, double scale, double l2) {
        for (int i = 0; i < grad.length; i++) {
            grad[i] = grad[i] * scale + l2 * param[i];
        }
    }

    private static void adamUpdate(double[] p, double[] g, double[] m, double[] v, double lr, long t) {
        final double beta1 = 0.9;
        final double beta2 = 0.999;
        final double eps = 1e-8;
        double b1t = 1.0 - Math.pow(beta1, t);
        double b2t = 1.0 - Math.pow(beta2, t);
        for (int i = 0; i < p.length; i++) {
            m[i] = beta1 * m[i] + (1.0 - beta1) * g[i];
            v[i] = beta2 * v[i] + (1.0 - beta2) * g[i] * g[i];
            double mhat = m[i] / b1t;
            double vhat = v[i] / b2t;
            p[i] -= lr * mhat / (Math.sqrt(vhat) + eps);
        }
    }

    private static double computeLoss(
            double[][] xNorm,
            double[] yNorm,
            int[] idx,
            double[] w1, double[] b1,
            double[] w2, double[] b2,
            double[] w3, double b3,
            int inputDim, int h1, int h2
    ) {
        double loss = 0.0;
        double[] a1 = new double[h1];
        double[] a2 = new double[h2];
        for (int id : idx) {
            double[] x = xNorm[id];
            for (int i = 0; i < h1; i++) {
                double z = b1[i];
                int row = i * inputDim;
                for (int d = 0; d < inputDim; d++) {
                    z += w1[row + d] * x[d];
                }
                a1[i] = Math.tanh(z);
            }
            for (int j = 0; j < h2; j++) {
                double z = b2[j];
                int row = j * h1;
                for (int i = 0; i < h1; i++) {
                    z += w2[row + i] * a1[i];
                }
                a2[j] = Math.tanh(z);
            }
            double yhat = b3;
            for (int j = 0; j < h2; j++) {
                yhat += w3[j] * a2[j];
            }
            double d = yhat - yNorm[id];
            loss += d * d;
        }
        return loss / idx.length;
    }

    private static void shuffle(int[] a, Random rand) {
        for (int i = a.length - 1; i > 0; i--) {
            int j = rand.nextInt(i + 1);
            int tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }

    private static void writeArray(DataOutput out, double[] x) throws IOException {
        out.writeInt(x.length);
        for (double v : x) {
            out.writeDouble(v);
        }
    }

    private static double[] readArray(DataInput in) throws IOException {
        int n = in.readInt();
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = in.readDouble();
        }
        return out;
    }
}
