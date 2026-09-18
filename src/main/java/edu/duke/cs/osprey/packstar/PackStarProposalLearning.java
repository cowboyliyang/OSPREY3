package edu.duke.cs.osprey.packstar;

import java.util.*;

/** Source-weighted fitting utilities. All inputs are adaptation data only. */
final class PackStarProposalLearning {
    private PackStarProposalLearning() {}

    static double logAdd(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        double max = Math.max(a, b);
        return max + Math.log1p(Math.exp(Math.min(a, b) - max));
    }

    // Group repeated conformations so they never cross the inner fit/score split.
    static int fold(int[] conf) {
        long h = 0xcbf29ce484222325L;
        for (int rc : conf) { h ^= rc; h *= 0x100000001b3L; }
        h ^= h >>> 33; h *= 0xff51afd7ed558ccdL; h ^= h >>> 33;
        return (int) (h & 1L);
    }

    static final class Data {
        final int[][] conf;
        final double[] residual, logWeight;
        final String[] contexts;
        final int[] folds;
        Data(int[][] conf, double[] residual, double[] logWeight) {
            if (conf.length != residual.length || conf.length != logWeight.length)
                throw new IllegalArgumentException("proposal sample shape mismatch");
            contexts = new String[conf.length];
            folds = new int[conf.length];
            for (int i = 0; i < conf.length; i++) {
                if (conf[i] == null || !Double.isFinite(residual[i])
                        || !Double.isFinite(logWeight[i]))
                    throw new IllegalArgumentException("invalid proposal sample/source weight: " + i);
                contexts[i] = Arrays.toString(conf[i]);
                folds[i] = fold(conf[i]);
            }
            this.conf = conf; this.residual = residual; this.logWeight = logWeight;
        }
    }

    static final class JointFit {
        final double[] delta;
        final double offset;
        JointFit(double[] delta, double offset) { this.delta = delta; this.offset = offset; }
    }

    /** Weighted joint residual regression around the termwise prior.
     * Penalty per cell = its probability mass * shrink / distinct contexts.
     * Sparse cells and unseen cells therefore stay close to the existing model.
     */
    static JointFit jointFit(int[][] features, double[] residual, double[] logWeights,
                             double[] shrinkRatios, double cap, int sweeps) {
        int n = residual.length, p = shrinkRatios.length;
        if (n == 0 || features.length != n || logWeights.length != n
                || !(cap > 0) || !Double.isFinite(cap) || sweeps < 1)
            throw new IllegalArgumentException("invalid joint fit inputs");
        double logSum = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(logWeights[i]) || !Double.isFinite(residual[i]))
                throw new IllegalArgumentException("nonfinite joint fit sample");
            logSum = logAdd(logSum, logWeights[i]);
        }
        List<List<Integer>> rows = new ArrayList<>(p);
        for (int j = 0; j < p; j++) {
            if (!(shrinkRatios[j] >= 0) || !Double.isFinite(shrinkRatios[j]))
                throw new IllegalArgumentException("invalid joint regularizer");
            rows.add(new ArrayList<>());
        }
        double[] weights = new double[n], errors = residual.clone(), delta = new double[p];
        double[] mass = new double[p];
        for (int i = 0; i < n; i++) {
            weights[i] = Math.exp(logWeights[i] - logSum);
            for (int j : features[i]) { rows.get(j).add(i); mass[j] += weights[i]; }
        }
        double offset = 0;
        for (int sweep = 0; sweep < sweeps; sweep++) {
            double shift = 0, maxChange = 0;
            for (int i = 0; i < n; i++) shift += weights[i] * errors[i];
            offset += shift;
            for (int i = 0; i < n; i++) errors[i] -= shift;
            for (int j = 0; j < p; j++) {
                if (!(mass[j] > 0)) continue;
                double numerator = 0;
                for (int i : rows.get(j)) numerator += weights[i] * (errors[i] + delta[j]);
                double value = Math.max(-cap, Math.min(cap,
                        numerator / (mass[j] * (1 + shrinkRatios[j]))));
                double change = value - delta[j];
                delta[j] = value;
                for (int i : rows.get(j)) errors[i] -= change;
                maxChange = Math.max(maxChange, Math.abs(change));
            }
            if (maxChange < 1e-8 && Math.abs(shift) < 1e-8) break;
        }
        return new JointFit(delta, offset);
    }

    static final class MomentCell {
        double logMass = Double.NEGATIVE_INFINITY;
        double logMassSquared = Double.NEGATIVE_INFINITY;
        double logMoment = Double.NEGATIVE_INFINITY;
        final Set<String> contexts = new HashSet<>();
        void add(int[] conf, double a, double r, double rt) {
            add(Arrays.toString(conf), a, r, rt);
        }
        void add(String context, double a, double r, double rt) {
            logMass = logAdd(logMass, a);
            logMassSquared = logAdd(logMassSquared, 2 * a);
            logMoment = logAdd(logMoment, a - 2 * r / rt);
            contexts.add(context);
        }
        double correction(double globalLogMoment, double rt, int minContexts,
                          double strength, double cap) {
            if (contexts.size() < minContexts) return 0;
            double ess = Math.exp(Math.min(700, 2 * logMass - logMassSquared));
            double support = Math.min(contexts.size(), ess);
            double fraction = support / (support + strength);
            double logRatio = logMoment - logMass - globalLogMoment;
            // Shrink the conditional moment toward the global moment, so an
            // unobserved cell has exactly zero correction in the same gauge.
            double shrunk = fraction == 1 ? logRatio : logAdd(
                    Math.log(fraction) + logRatio, Math.log1p(-fraction));
            return Math.max(-cap, Math.min(cap, -0.5 * rt * shrunk));
        }
    }

    /** Anchor the tolerance to the global minimum, avoiding order-dependent
     * chains of pairwise 5% comparisons. The comparator orders complexity first.
     */
    static <T> T chooseWithinFivePercent(List<T> candidates,
            java.util.function.ToLongFunction<T> samples, Comparator<T> complexity) {
        long minimum = Long.MAX_VALUE;
        for (T candidate : candidates) minimum = Math.min(minimum, samples.applyAsLong(candidate));
        T best = null;
        for (T candidate : candidates) {
            if (samples.applyAsLong(candidate) > 1.05 * minimum + 1e-10) continue;
            if (best == null || complexity.compare(candidate, best) < 0) best = candidate;
        }
        return best;
    }

    static final class MomentFit {
        final double[] parameters;
        final double objective;
        final int iterations;
        final boolean converged;

        MomentFit(double[] parameters, double objective, int iterations, boolean converged) {
            this.parameters = parameters;
            this.objective = objective;
            this.iterations = iterations;
            this.converged = converged;
        }
    }

    /** Convex, box-constrained joint log-rho fit. Parameters are in kcal/mol.
     * The fixed pair proposal defines a and r. The penalty is
     * 0.5 * sum_j penalty[j] * (h_j/RT)^2; unsupported cells have no parameter.
     * Projected accelerated gradient with backtracking uses sparse cell indices.
     */
    static MomentFit jointMomentFit(int[][] features, double[] residual, double[] logWeights,
            double[] penalty, double[] initial, double rt, double cap, int maxIterations) {
        return jointMomentFit(features, residual, logWeights, penalty, initial, rt, cap, maxIterations, false);
    }

    static MomentFit jointMomentFit(int[][] features, double[] residual, double[] logWeights,
            double[] penalty, double[] initial, double rt, double cap, int maxIterations,
            boolean shiftedSums) {
        int n = residual.length, p = penalty.length;
        if (n < 2 || features.length != n || logWeights.length != n || initial.length != p
                || !(rt > 0) || !Double.isFinite(rt) || !(cap > 0)
                || !Double.isFinite(cap) || maxIterations < 1)
            throw new IllegalArgumentException("invalid joint moment fit inputs");
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(residual[i]) || !Double.isFinite(logWeights[i]))
                throw new IllegalArgumentException("nonfinite joint moment sample");
            for (int j : features[i]) if (j < 0 || j >= p)
                throw new IllegalArgumentException("invalid moment feature");
        }
        for (int j = 0; j < p; j++)
            if (!(penalty[j] >= 0) || !Double.isFinite(penalty[j]) || !Double.isFinite(initial[j]))
                throw new IllegalArgumentException("invalid moment regularizer/initial value");
        double[] x = initial.clone();
        for (int j = 0; j < p; j++) x[j] = Math.max(-cap, Math.min(cap, x[j]));
        double[] y = x.clone(), gradient = new double[p];
        MomentObjective objective = new MomentObjective(features, residual, logWeights, penalty, rt, shiftedSums);
        double value = objective.value(x, null);
        double momentum = 1, lipschitz = 1 / (rt * rt);
        for (int iteration = 1; iteration <= maxIterations; iteration++) {
            double atY = objective.value(y, gradient);
            double[] next = new double[p];
            double nextValue;
            while (true) {
                double linear = 0, norm = 0;
                for (int j = 0; j < p; j++) {
                    next[j] = Math.max(-cap, Math.min(cap, y[j] - gradient[j] / lipschitz));
                    double d = next[j] - y[j];
                    linear += gradient[j] * d;
                    norm += d * d;
                }
                nextValue = objective.value(next, null);
                if (Double.isFinite(nextValue) && nextValue <= atY + linear + 0.5 * lipschitz * norm + 1e-12) break;
                lipschitz *= 2;
                if (!Double.isFinite(lipschitz)) throw new IllegalStateException("moment line search failed");
            }
            if (nextValue > value + 1e-12) {
                y = x.clone();
                momentum = 1;
                continue;
            }
            // Check the projected gradient at the accepted point, not merely
            // objective stagnation (which is unreliable for sparse cells).
            // The accepted line-search evaluation already owns the identical
            // per-sample terms and normalizers. Only its gradient is missing.
            objective.gradientAtLastValue(next, gradient);
            double stationarity = 0;
            for (int j = 0; j < p; j++) stationarity = Math.max(stationarity,
                    Math.abs(next[j] - Math.max(-cap, Math.min(cap, next[j] - gradient[j]))));
            if (stationarity <= 1e-7)
                return new MomentFit(next, nextValue, iteration, true);
            double nextMomentum = 0.5 * (1 + Math.sqrt(1 + 4 * momentum * momentum));
            for (int j = 0; j < p; j++)
                y[j] = next[j] + (momentum - 1) / nextMomentum * (next[j] - x[j]);
            x = next;
            value = nextValue;
            momentum = nextMomentum;
        }
        return new MomentFit(x, value, maxIterations, false);
    }

    static double momentObjective(int[][] features, double[] residual, double[] logWeights,
            double[] penalty, double[] parameters, double rt, double[] gradient) {
        return new MomentObjective(features, residual, logWeights, penalty, rt).value(parameters, gradient);
    }

    /** One optimizer invocation owns its scratch arrays. Centered two-pass
     * log sums use one logarithm per reduction rather than one per sample.
     * All fitting, line-search and convergence thresholds remain unchanged.
     */
    static final class MomentObjective {
        private final int[][] features;
        private final double[] penalty, source, secondBase, a, b, scaledParameters;
        private final double rt, logFirst;
        private double logA, logB;
        private double sumA, sumB;
        private final boolean shiftedSums;

        MomentObjective(int[][] features, double[] residual, double[] logWeights,
                        double[] penalty, double rt) {
            this(features, residual, logWeights, penalty, rt, false);
        }

        MomentObjective(int[][] features, double[] residual, double[] logWeights,
                        double[] penalty, double rt, boolean shiftedSums) {
            this.shiftedSums = shiftedSums;
            this.penalty = penalty;
            this.rt = rt;
            scaledParameters = new double[penalty.length];
            int n = residual.length;
            double[] rawSource = new double[n], rawSecond = new double[n];
            double sourceShift = Double.NEGATIVE_INFINITY;
            for (double weight : logWeights) sourceShift = Math.max(sourceShift, weight);
            double residualShift = residual[0] / rt;
            double first = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < n; i++) {
                rawSource[i] = logWeights[i] - sourceShift;
                double r = residual[i] / rt - residualShift;
                rawSecond[i] = rawSource[i] - 2 * r;
                first = logAdd(first, rawSource[i] - r);
            }
            logFirst = first;
            if (shiftedSums) {
                // Rows with identical active cells have the same h for EVERY
                // parameter vector. Their two weighted masses are sufficient
                // for the objective and gradient. Preserve both masses rather
                // than averaging residuals or dropping draw multiplicity.
                Map<java.nio.IntBuffer, Integer> groups = new LinkedHashMap<>();
                List<int[]> groupedFeatures = new ArrayList<>();
                double[] groupedSource = new double[n], groupedSecond = new double[n];
                for (int i = 0; i < n; i++) {
                    java.nio.IntBuffer key = java.nio.IntBuffer.wrap(features[i]);
                    Integer group = groups.get(key);
                    if (group == null) {
                        group = groups.size();
                        groups.put(key, group);
                        groupedFeatures.add(features[i]);
                        groupedSource[group] = rawSource[i];
                        groupedSecond[group] = rawSecond[i];
                    } else {
                        groupedSource[group] = logAdd(groupedSource[group], rawSource[i]);
                        groupedSecond[group] = logAdd(groupedSecond[group], rawSecond[i]);
                    }
                }
                this.features = groupedFeatures.toArray(new int[0][]);
                source = Arrays.copyOf(groupedSource, groups.size());
                secondBase = Arrays.copyOf(groupedSecond, groups.size());
            } else {
                this.features = features;
                source = rawSource;
                secondBase = rawSecond;
            }
            a = new double[source.length];
            b = new double[source.length];
        }

        double value(double[] parameters, double[] gradient) {
            logA = Double.NEGATIVE_INFINITY;
            logB = logA;
            // Each parameter may occur in many rows. Preserve division and
            // feature summation order, but divide each parameter only once.
            for (int j = 0; j < parameters.length; j++) scaledParameters[j] = parameters[j] / rt;
            for (int i = 0; i < source.length; i++) {
                double h = 0;
                for (int j : features[i]) h += scaledParameters[j];
                a[i] = source[i] - h;
                b[i] = secondBase[i] + h;
                logA = shiftedSums ? Math.max(logA, a[i]) : logAdd(logA, a[i]);
                logB = shiftedSums ? Math.max(logB, b[i]) : logAdd(logB, b[i]);
            }
            if (shiftedSums) {
                sumA = 0; sumB = 0;
                double compensationA = 0, compensationB = 0;
                for (int i = 0; i < source.length; i++) {
                    // Kahan accumulation retains small source-weight terms in a
                    // long reduction. Max shifting prevents exponential overflow.
                    // Retain the shifted masses for the gradient, including
                    // the gradient requested after an accepted line search.
                    a[i] = Math.exp(a[i] - logA);
                    b[i] = Math.exp(b[i] - logB);
                    double termA = a[i] - compensationA;
                    double nextA = sumA + termA;
                    compensationA = (nextA - sumA) - termA;
                    sumA = nextA;
                    double termB = b[i] - compensationB;
                    double nextB = sumB + termB;
                    compensationB = (nextB - sumB) - termB;
                    sumB = nextB;
                }
                logA += Math.log(sumA);
                logB += Math.log(sumB);
            }
            if (gradient != null) gradientAtLastValue(parameters, gradient);
            double objective = logA + logB - 2 * logFirst;
            for (int j = 0; j < parameters.length; j++) {
                objective += 0.5 * penalty[j] * Math.pow(parameters[j] / rt, 2);
            }
            return objective;
        }

        /** parameters must be the unchanged point of the most recent value().
         * Avoid a second pair of log-sum-exp reductions at the same point.
         */
        void gradientAtLastValue(double[] parameters, double[] gradient) {
            Arrays.fill(gradient, 0);
            for (int i = 0; i < source.length; i++) {
                double g = shiftedSums
                        ? (b[i] / sumB - a[i] / sumA) / rt
                        : (Math.exp(b[i] - logB) - Math.exp(a[i] - logA)) / rt;
                for (int j : features[i]) gradient[j] += g;
            }
            for (int j = 0; j < parameters.length; j++)
                gradient[j] += penalty[j] * parameters[j] / (rt * rt);
        }
    }

    /** Empirical log rho under the candidate, including its changed normalizer.
     * log rho = log sum exp(a-h/RT) + log sum exp(a-2r/RT+h/RT)
     *           - 2 log sum exp(a-r/RT).
     * The source constants and draw multiplicity must remain in a.
     */
    static double logRho(Data data, double[] correction, int heldOutFold, double rt) {
        double mass = Double.NEGATIVE_INFINITY, second = mass, first = mass;
        int count = 0;
        for (int i = 0; i < data.conf.length; i++) {
            if (heldOutFold >= 0 && data.folds[i] != heldOutFold) continue;
            double a = data.logWeight[i], r = data.residual[i] / rt, h = correction[i] / rt;
            mass = logAdd(mass, a - h);
            second = logAdd(second, a - 2 * r + h);
            first = logAdd(first, a - r);
            count++;
        }
        return count < 2 ? Double.NaN : mass + second - 2 * first;
    }
}
