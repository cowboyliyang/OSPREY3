package edu.duke.cs.osprey.packstar;

import java.util.Random;

/** Isolate fitter cost; synthetic timings are not full protein speedups. */
public final class MomentFitSpeedBench {
    public static void main(String[] args) {
        System.out.println("samples\tscopes\tsequential_s\tgrouped_s\tspeedup\tmax_parameter_delta");
        for (int scopes = 1; scopes <= 3; scopes++) {
            int n = 1000, cells = 16, p = cells * scopes, repeats = 30;
            int[][] features = new int[n][scopes];
            double[] residual = new double[n], weights = new double[n];
            double[] penalty = new double[p], initial = new double[p];
            Random random = new Random(9871 + scopes);
            for (int j = 0; j < p; j++) penalty[j] = 0.002;
            for (int i = 0; i < n; i++) {
                for (int t = 0; t < scopes; t++) {
                    int cell = random.nextInt(cells);
                    features[i][t] = t * cells + cell;
                    residual[i] += (cell - 7.5) * (0.2 / (t + 1));
                }
                residual[i] += random.nextGaussian() * 0.03;
                weights[i] = random.nextDouble() * 4 - 2;
            }
            for (int warm = 0; warm < 5; warm++) for (boolean grouped : new boolean[]{false, true})
                PackStarProposalLearning.jointMomentFit(features, residual, weights, penalty, initial, 0.6, 3, 400, grouped);
            double[] seconds = new double[2];
            double maximumDelta = 0;
            for (int repeat = 0; repeat < repeats; repeat++) {
                PackStarProposalLearning.MomentFit[] fits = new PackStarProposalLearning.MomentFit[2];
                for (int order = 0; order < 2; order++) {
                    int mode = (order + repeat) % 2;
                    long start = System.nanoTime();
                    fits[mode] = PackStarProposalLearning.jointMomentFit(features, residual, weights, penalty, initial, 0.6, 3, 400, mode == 1);
                    seconds[mode] += (System.nanoTime() - start) / 1e9;
                }
                if (fits[0].converged != fits[1].converged) throw new AssertionError("convergence changed");
                for (int j = 0; j < p; j++) maximumDelta = Math.max(maximumDelta, Math.abs(fits[0].parameters[j] - fits[1].parameters[j]));
                if (maximumDelta > 1e-8) throw new AssertionError("fit changed by " + maximumDelta);
            }
            System.out.printf(java.util.Locale.ROOT, "%d\t%d\t%.6f\t%.6f\t%.3f\t%.4g%n",
                    n, scopes, seconds[0], seconds[1], seconds[0]/seconds[1], maximumDelta);
        }
    }
}
