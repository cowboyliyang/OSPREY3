package edu.duke.cs.osprey.packstar;

import java.util.Arrays;

/**
 * Exact finite-state algebra for the variance-targeted eta-v4 proposal.
 *
 * <p>For an unnormalised positive target {@code f} and proposal {@code q},
 * the importance-sampling second moment is {@code sum(f^2/q)}.  Its value
 * divided by {@code Z^2}, minus one, is both the squared coefficient of
 * variation of the weights and {@code chi^2(p || q)} for
 * {@code p=f/Z}.  Unlike an energy-regression loss, this is the objective that
 * directly controls the variance of the PACK* partition-function estimator.
 *
 * <p>The routines here deliberately contain no sampling or mutable DP state.
 * They provide auditable log-space identities for exact toy systems and for
 * the later adaptation/control-flow integration.</p>
 */
final class PackStarEtaV4Importance {

    private static final double NORMALIZATION_TOLERANCE = 1.0e-10;
    private static final double BOUND_TOLERANCE = 1.0e-12;

    private PackStarEtaV4Importance() {
    }

    static final class Objective {
        final double logSecondMoment;
        final double chiSquare;
        /** Gradient of log(second moment) in the supplied coordinates. */
        final double[] gradient;

        Objective(double logSecondMoment, double chiSquare,
                  double[] gradient) {
            this.logSecondMoment = logSecondMoment;
            this.chiSquare = chiSquare;
            this.gradient = gradient;
        }
    }

    static final class SparseGradient {
        final double[] gradient;
        final double tailEffectiveSampleSize;
        final double minLogImportanceWeight;
        final double maxLogImportanceWeight;

        SparseGradient(double[] gradient,
                       double tailEffectiveSampleSize,
                       double minLogImportanceWeight,
                       double maxLogImportanceWeight) {
            this.gradient = gradient;
            this.tailEffectiveSampleSize = tailEffectiveSampleSize;
            this.minLogImportanceWeight = minLogImportanceWeight;
            this.maxLogImportanceWeight = maxLogImportanceWeight;
        }
    }

    /**
     * Log density of {@code qMix = sum_k alpha[k] q[k]} at one state.
     * Mixture weights must already be frozen and normalized; this method does
     * not silently renormalize malformed inputs.
     */
    static double logMixtureProbability(
            double[] logComponentProbabilities, double[] mixtureWeights) {
        if (!validMixtureInputs(logComponentProbabilities, mixtureWeights)) {
            return Double.NaN;
        }
        double max = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < mixtureWeights.length; k++) {
            if (mixtureWeights[k] == 0.0) continue;
            double term = Math.log(mixtureWeights[k])
                    + logComponentProbabilities[k];
            max = Math.max(max, term);
        }
        if (!Double.isFinite(max)) return Double.NEGATIVE_INFINITY;
        double sum = 0.0;
        for (int k = 0; k < mixtureWeights.length; k++) {
            if (mixtureWeights[k] == 0.0) continue;
            sum += Math.exp(Math.log(mixtureWeights[k])
                    + logComponentProbabilities[k] - max);
        }
        if (!Double.isFinite(sum) || !(sum > 0.0)) return Double.NaN;
        return max + Math.log(sum);
    }

    static double mixtureLogImportanceWeight(
            double logUnnormalisedTarget,
            double[] logComponentProbabilities,
            double[] mixtureWeights) {
        if (!Double.isFinite(logUnnormalisedTarget)) return Double.NaN;
        double logQMix = logMixtureProbability(
                logComponentProbabilities, mixtureWeights);
        if (!Double.isFinite(logQMix)) return Double.NaN;
        return logUnnormalisedTarget - logQMix;
    }

    /**
     * Conditional pointwise range inherited from one anchor component.
     *
     * <p>If {@code E_anchor(c)-E_true(c) <= B} globally, then
     * {@code f/q_anchor <= Z_anchor exp(B/RT)}.  Since
     * {@code q_mix >= alpha_anchor q_anchor}, the mixture weight is bounded by
     * {@code Z_anchor exp(B/RT)/alpha_anchor}.</p>
     */
    static double anchorLogRange(double anchorLogZ, double boundKcal,
                                 double rt, double anchorWeight) {
        if (!Double.isFinite(anchorLogZ)
                || !Double.isFinite(boundKcal) || boundKcal < 0.0
                || !Double.isFinite(rt) || !(rt > 0.0)
                || !Double.isFinite(anchorWeight)
                || !(anchorWeight > 0.0) || anchorWeight > 1.0) {
            return Double.NaN;
        }
        double result = anchorLogZ + boundKcal / rt
                - Math.log(anchorWeight);
        return Double.isFinite(result) ? result : Double.NaN;
    }

    /**
     * Pointwise range inherited from a complementary cover of components.
     *
     * <p>If, for every state, at least one positive-weight component
     * {@code k} satisfies {@code E_k-E_true <= B_k}, then
     *
     * <pre>
     * f/q_mix <= max_k Z_k exp(B_k/RT) / alpha_k.
     * </pre>
     *
     * Zero-weight components cannot contribute to the cover and are ignored.
     */
    static double componentCoverLogRange(
            double[] componentLogZ, double[] boundsKcal,
            double rt, double[] mixtureWeights) {
        if (componentLogZ == null || boundsKcal == null
                || mixtureWeights == null
                || componentLogZ.length == 0
                || componentLogZ.length != boundsKcal.length
                || componentLogZ.length != mixtureWeights.length
                || !Double.isFinite(rt) || !(rt > 0.0)
                || !validNormalizedWeights(mixtureWeights)) {
            return Double.NaN;
        }
        double logRange = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < componentLogZ.length; k++) {
            if (mixtureWeights[k] == 0.0) continue;
            double candidate = anchorLogRange(
                    componentLogZ[k], boundsKcal[k], rt,
                    mixtureWeights[k]);
            if (!Double.isFinite(candidate)) return Double.NaN;
            logRange = Math.max(logRange, candidate);
        }
        return Double.isFinite(logRange) ? logRange : Double.NaN;
    }

    /** Whether at least one positive-weight component covers this state. */
    static boolean componentCoverContains(
            double[] componentOverCorrectionsKcal,
            double[] boundsKcal, double[] mixtureWeights) {
        if (componentOverCorrectionsKcal == null || boundsKcal == null
                || mixtureWeights == null
                || componentOverCorrectionsKcal.length == 0
                || componentOverCorrectionsKcal.length
                != boundsKcal.length
                || componentOverCorrectionsKcal.length
                != mixtureWeights.length
                || !validNormalizedWeights(mixtureWeights)) {
            return false;
        }
        for (int k = 0; k < mixtureWeights.length; k++) {
            if (mixtureWeights[k] == 0.0) continue;
            if (!Double.isFinite(componentOverCorrectionsKcal[k])
                    || !Double.isFinite(boundsKcal[k])
                    || boundsKcal[k] < 0.0) return false;
            if (componentOverCorrectionsKcal[k]
                    <= boundsKcal[k] + BOUND_TOLERANCE) {
                return true;
            }
        }
        return false;
    }

    /**
     * Mixture weights that minimize the cover range for fixed components:
     * {@code alpha_k proportional to Z_k exp(B_k/RT)}.  At these weights every
     * positive component contributes the same cap and the common range is the
     * sum of the component caps.
     */
    static double[] rangeBalancedCoverWeights(
            double[] componentLogZ, double[] boundsKcal, double rt) {
        if (componentLogZ == null || boundsKcal == null
                || componentLogZ.length == 0
                || componentLogZ.length != boundsKcal.length
                || !Double.isFinite(rt) || !(rt > 0.0)) {
            return null;
        }
        double max = Double.NEGATIVE_INFINITY;
        double[] logCaps = new double[componentLogZ.length];
        for (int k = 0; k < componentLogZ.length; k++) {
            if (!Double.isFinite(componentLogZ[k])
                    || !Double.isFinite(boundsKcal[k])
                    || boundsKcal[k] < 0.0) return null;
            logCaps[k] = componentLogZ[k] + boundsKcal[k] / rt;
            max = Math.max(max, logCaps[k]);
        }
        if (!Double.isFinite(max)) return null;
        double sum = 0.0;
        double[] weights = new double[logCaps.length];
        for (int k = 0; k < logCaps.length; k++) {
            weights[k] = Math.exp(logCaps[k] - max);
            sum += weights[k];
        }
        if (!Double.isFinite(sum) || !(sum > 0.0)) return null;
        for (int k = 0; k < weights.length; k++) weights[k] /= sum;
        return validNormalizedWeights(weights) ? weights : null;
    }

    /**
     * Floating-point-safe gate for the intended absolute tail ESS boundary.
     * The tolerance is deliberately tiny relative to the threshold; proposal
     * candidates must still pass the independent held-out objective.
     */
    static boolean tailEssPasses(double observedEss, double minimumEss) {
        if (!Double.isFinite(observedEss)
                || !Double.isFinite(minimumEss)
                || !(minimumEss > 0.0)) return false;
        double tolerance = Math.max(1.0e-6, 1.0e-4 * minimumEss);
        return observedEss >= minimumEss - tolerance;
    }

    /**
     * Scale one mixture importance weight to its certified anchor range.
     * Values above one by more than floating-point tolerance fail closed.
     */
    static double boundedNormalizedMixtureWeight(
            double logUnnormalisedTarget,
            double[] logComponentProbabilities,
            double[] mixtureWeights,
            double anchorLogRange) {
        if (!Double.isFinite(anchorLogRange)) return Double.NaN;
        double logWeight = mixtureLogImportanceWeight(
                logUnnormalisedTarget,
                logComponentProbabilities, mixtureWeights);
        if (!Double.isFinite(logWeight)) return Double.NaN;
        double logNormalized = logWeight - anchorLogRange;
        if (logNormalized > Math.log1p(BOUND_TOLERANCE)) {
            return Double.NaN;
        }
        if (logNormalized >= 0.0) return 1.0;
        double normalized = Math.exp(logNormalized);
        return Double.isFinite(normalized) ? normalized : Double.NaN;
    }

    /**
     * Exact chi-square objective and gradient for one exponential-family
     * proposal over an enumerated finite state space.
     *
     * <p>For {@code q_theta(x)=exp(theta*T(x)-A(theta))},
     *
     * <pre>
     * grad_theta log sum_x f(x)^2/q_theta(x)
     *     = E_q[T] - E_r[T],  r(x) proportional to f(x)^2/q_theta(x).
     * </pre>
     *
     * A gradient-descent step therefore moves proposal mass toward states
     * dominating the squared importance weights.</p>
     */
    static Objective exactExponentialFamilyObjective(
            double[] logUnnormalisedTarget,
            double[] logProposalProbabilities,
            double[][] sufficientStatistics) {
        int states = validateFiniteStateInputs(
                logUnnormalisedTarget, logProposalProbabilities,
                sufficientStatistics);
        if (states < 1) return null;
        int dimension = sufficientStatistics[0].length;

        double logZ = logSumExp(logUnnormalisedTarget);
        double[] secondMomentTerms = new double[states];
        for (int state = 0; state < states; state++) {
            secondMomentTerms[state] = 2.0 * logUnnormalisedTarget[state]
                    - logProposalProbabilities[state];
        }
        double logSecondMoment = logSumExp(secondMomentTerms);
        if (!Double.isFinite(logZ) || !Double.isFinite(logSecondMoment)) {
            return null;
        }

        double[] proposalMean = new double[dimension];
        double[] squaredWeightMean = new double[dimension];
        for (int state = 0; state < states; state++) {
            double q = Math.exp(logProposalProbabilities[state]);
            double r = Math.exp(secondMomentTerms[state]
                    - logSecondMoment);
            for (int feature = 0; feature < dimension; feature++) {
                double statistic = sufficientStatistics[state][feature];
                proposalMean[feature] += q * statistic;
                squaredWeightMean[feature] += r * statistic;
            }
        }
        double[] gradient = new double[dimension];
        for (int feature = 0; feature < dimension; feature++) {
            gradient[feature] = proposalMean[feature]
                    - squaredWeightMean[feature];
            if (!Double.isFinite(gradient[feature])) return null;
        }
        return new Objective(logSecondMoment,
                chiSquareFromLogs(logSecondMoment, logZ), gradient);
    }

    /**
     * Sparse sample analogue of the exponential-family objective gradient.
     *
     * <p>{@code baseWeights} describe the empirical proposal measure.  They
     * are normalized internally, so an IID on-policy batch uses all ones.
     * For {@code tailPower=2}, tail weights are proportional to
     * {@code baseWeight * importanceWeight^2} and the returned vector is
     * {@code E_q[T]-E_tail[T]}.  Smaller powers are permitted only as an
     * adaptation homotopy; candidate acceptance must still use the untampered
     * chi-square objective on fresh data.</p>
     */
    static SparseGradient sparseChiSquareGradient(
            int dimension, int[][] featureIndices,
            double[] baseWeights, double[] logImportanceWeights,
            double tailPower) {
        if (dimension <= 0 || featureIndices == null
                || baseWeights == null || logImportanceWeights == null
                || featureIndices.length == 0
                || featureIndices.length != baseWeights.length
                || featureIndices.length != logImportanceWeights.length
                || !Double.isFinite(tailPower) || !(tailPower > 0.0)
                || tailPower > 2.0) return null;

        double sumBase = 0.0;
        double maxTailLog = Double.NEGATIVE_INFINITY;
        double minLogWeight = Double.POSITIVE_INFINITY;
        double maxLogWeight = Double.NEGATIVE_INFINITY;
        for (int sample = 0; sample < featureIndices.length; sample++) {
            if (featureIndices[sample] == null
                    || !Double.isFinite(baseWeights[sample])
                    || baseWeights[sample] < 0.0
                    || !Double.isFinite(logImportanceWeights[sample])) {
                return null;
            }
            for (int index : featureIndices[sample]) {
                if (index < 0 || index >= dimension) return null;
            }
            sumBase += baseWeights[sample];
            minLogWeight = Math.min(minLogWeight,
                    logImportanceWeights[sample]);
            maxLogWeight = Math.max(maxLogWeight,
                    logImportanceWeights[sample]);
            if (baseWeights[sample] > 0.0) {
                maxTailLog = Math.max(maxTailLog,
                        Math.log(baseWeights[sample])
                                + tailPower
                                * logImportanceWeights[sample]);
            }
        }
        if (!Double.isFinite(sumBase) || !(sumBase > 0.0)
                || !Double.isFinite(maxTailLog)) return null;

        double sumTailScaled = 0.0;
        double[] tailScaled = new double[featureIndices.length];
        for (int sample = 0; sample < featureIndices.length; sample++) {
            if (baseWeights[sample] == 0.0) continue;
            tailScaled[sample] = Math.exp(
                    Math.log(baseWeights[sample])
                            + tailPower * logImportanceWeights[sample]
                            - maxTailLog);
            sumTailScaled += tailScaled[sample];
        }
        if (!Double.isFinite(sumTailScaled) || !(sumTailScaled > 0.0)) {
            return null;
        }

        double[] gradient = new double[dimension];
        double tailWeightSquares = 0.0;
        for (int sample = 0; sample < featureIndices.length; sample++) {
            double qWeight = baseWeights[sample] / sumBase;
            double tailWeight = tailScaled[sample] / sumTailScaled;
            tailWeightSquares += tailWeight * tailWeight;
            double difference = qWeight - tailWeight;
            for (int index : featureIndices[sample]) {
                gradient[index] += difference;
            }
        }
        if (!Double.isFinite(tailWeightSquares)
                || !(tailWeightSquares > 0.0)) return null;
        for (double value : gradient) {
            if (!Double.isFinite(value)) return null;
        }
        return new SparseGradient(gradient,
                1.0 / tailWeightSquares,
                minLogWeight, maxLogWeight);
    }

    /**
     * Exact objective and gradient with respect to frozen mixture weights.
     * For fixed components, {@code sum f^2/(sum alpha*q)} is convex on the
     * simplex.  The returned gradient is for its logarithm; optimization must
     * still project onto the simplex and preserve any anchor-weight floor.
     */
    static Objective exactMixtureObjective(
            double[] logUnnormalisedTarget,
            double[][] logComponentProbabilities,
            double[] mixtureWeights) {
        if (logUnnormalisedTarget == null
                || logComponentProbabilities == null
                || mixtureWeights == null
                || logUnnormalisedTarget.length == 0
                || logComponentProbabilities.length == 0
                || logComponentProbabilities.length
                != mixtureWeights.length
                || !validNormalizedWeights(mixtureWeights)) {
            return null;
        }
        int states = logUnnormalisedTarget.length;
        for (double value : logUnnormalisedTarget) {
            if (!Double.isFinite(value)) return null;
        }
        for (double[] component : logComponentProbabilities) {
            if (component == null || component.length != states
                    || !isNormalizedLogProbability(component)) return null;
        }

        double[] logMixture = new double[states];
        double[] secondMomentTerms = new double[states];
        double[] atState = new double[mixtureWeights.length];
        for (int state = 0; state < states; state++) {
            for (int component = 0;
                 component < mixtureWeights.length; component++) {
                atState[component] =
                        logComponentProbabilities[component][state];
            }
            logMixture[state] = logMixtureProbability(
                    atState, mixtureWeights);
            if (!Double.isFinite(logMixture[state])) return null;
            secondMomentTerms[state] =
                    2.0 * logUnnormalisedTarget[state]
                            - logMixture[state];
        }
        double logZ = logSumExp(logUnnormalisedTarget);
        double logSecondMoment = logSumExp(secondMomentTerms);
        if (!Double.isFinite(logZ) || !Double.isFinite(logSecondMoment)) {
            return null;
        }

        double[] gradient = new double[mixtureWeights.length];
        for (int state = 0; state < states; state++) {
            double r = Math.exp(secondMomentTerms[state]
                    - logSecondMoment);
            for (int component = 0;
                 component < mixtureWeights.length; component++) {
                double componentOverMixture = Math.exp(
                        logComponentProbabilities[component][state]
                                - logMixture[state]);
                gradient[component] -= r * componentOverMixture;
            }
        }
        for (double value : gradient) {
            if (!Double.isFinite(value)) return null;
        }
        return new Objective(logSecondMoment,
                chiSquareFromLogs(logSecondMoment, logZ), gradient);
    }

    private static int validateFiniteStateInputs(
            double[] logUnnormalisedTarget,
            double[] logProposalProbabilities,
            double[][] sufficientStatistics) {
        if (logUnnormalisedTarget == null
                || logProposalProbabilities == null
                || sufficientStatistics == null
                || logUnnormalisedTarget.length == 0
                || logUnnormalisedTarget.length
                != logProposalProbabilities.length
                || logUnnormalisedTarget.length
                != sufficientStatistics.length
                || !isNormalizedLogProbability(
                logProposalProbabilities)) return -1;
        int dimension = sufficientStatistics[0] == null
                ? -1 : sufficientStatistics[0].length;
        if (dimension < 1) return -1;
        for (int state = 0;
             state < logUnnormalisedTarget.length; state++) {
            if (!Double.isFinite(logUnnormalisedTarget[state])
                    || sufficientStatistics[state] == null
                    || sufficientStatistics[state].length != dimension) {
                return -1;
            }
            for (double statistic : sufficientStatistics[state]) {
                if (!Double.isFinite(statistic)) return -1;
            }
        }
        return logUnnormalisedTarget.length;
    }

    private static double chiSquareFromLogs(
            double logSecondMoment, double logZ) {
        double logRatio = logSecondMoment - 2.0 * logZ;
        if (logRatio > Math.log(Double.MAX_VALUE)) {
            return Double.POSITIVE_INFINITY;
        }
        double result = Math.expm1(logRatio);
        if (!Double.isFinite(result)) return Double.POSITIVE_INFINITY;
        // Exact arithmetic is nonnegative; tolerate cancellation near zero.
        return Math.max(0.0, result);
    }

    private static boolean validMixtureInputs(
            double[] logComponentProbabilities,
            double[] mixtureWeights) {
        if (logComponentProbabilities == null || mixtureWeights == null
                || logComponentProbabilities.length == 0
                || logComponentProbabilities.length
                != mixtureWeights.length
                || !validNormalizedWeights(mixtureWeights)) return false;
        for (double logProbability : logComponentProbabilities) {
            if (!(Double.isFinite(logProbability)
                    || logProbability == Double.NEGATIVE_INFINITY)) {
                return false;
            }
        }
        return true;
    }

    private static boolean validNormalizedWeights(double[] weights) {
        if (weights == null || weights.length == 0) return false;
        double sum = 0.0;
        for (double weight : weights) {
            if (!Double.isFinite(weight) || weight < 0.0
                    || weight > 1.0) return false;
            sum += weight;
        }
        return Double.isFinite(sum)
                && Math.abs(sum - 1.0) <= NORMALIZATION_TOLERANCE;
    }

    private static boolean isNormalizedLogProbability(
            double[] logProbabilities) {
        if (logProbabilities == null
                || logProbabilities.length == 0) return false;
        for (double value : logProbabilities) {
            if (!(Double.isFinite(value)
                    || value == Double.NEGATIVE_INFINITY)) return false;
        }
        double logSum = logSumExp(logProbabilities);
        return Double.isFinite(logSum)
                && Math.abs(logSum) <= NORMALIZATION_TOLERANCE;
    }

    private static double logSumExp(double[] values) {
        if (values == null || values.length == 0) return Double.NaN;
        double max = Arrays.stream(values).max()
                .orElse(Double.NEGATIVE_INFINITY);
        if (!Double.isFinite(max)) return max;
        double sum = 0.0;
        for (double value : values) {
            if (!(Double.isFinite(value)
                    || value == Double.NEGATIVE_INFINITY)) {
                return Double.NaN;
            }
            sum += Math.exp(value - max);
        }
        return max + Math.log(sum);
    }
}
