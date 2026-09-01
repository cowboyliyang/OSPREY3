/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
*/

package edu.duke.cs.osprey.packstar;

import org.apache.commons.math3.distribution.BetaDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * Frequency/severity PAC confidence interval for importance weights.
 *
 * <p>For a positive relative importance weight {@code R}, a fixed relative
 * clipping threshold {@code C}, and {@code U=max(R/C-1,0)}, the interval uses
 * an empirical-Bernstein interval for {@code Y=min(R/C,1)} and a one-sided
 * exact binomial limit for {@code P(U>0)}.  The tail factorization uses the
 * externally supplied conditional severity premise
 * {@code E[U | U>0] <= S0}.</p>
 *
 * <p>This class deliberately works on the relative scale.  A caller may use
 * any positive, adaptation-only normalizer {@code mu} and set {@code R=W/mu};
 * the corresponding partition-function scale is then
 * {@code Zproposal * mu * C}.</p>
 */
public final class PackStarFrequencySeverityPAC {

    private static final double[] DEFAULT_BETTING_LAMBDAS = {
            1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3,
            1.0e-2, 3.0e-2, 0.1, 0.3, 0.5, 0.8
    };

    private PackStarFrequencySeverityPAC() {}

    static double[] severityBettingLambdas() {
        return DEFAULT_BETTING_LAMBDAS.clone();
    }

    /** Point estimates used only for adaptation-time sizing and ranking. */
    public static final class Moments {
        public final double bulkMean;
        public final double bulkVariance;
        public final double tailProbability;

        public Moments(double bulkMean, double bulkVariance,
                       double tailProbability) {
            requireProbability(bulkMean, "bulk mean");
            if (!Double.isFinite(bulkVariance) || bulkVariance < 0.0) {
                throw new IllegalArgumentException(
                        "bulk variance must be finite and nonnegative");
            }
            requireProbability(tailProbability, "tail probability");
            this.bulkMean = bulkMean;
            this.bulkVariance = bulkVariance;
            this.tailProbability = tailProbability;
        }
    }

    /**
     * Aggregated cross-fit diagnostics for fold-specific, self-normalized
     * proposal scores.  The folds may use different fitted eta values, so this
     * object describes their sample-size-weighted mixture and is adaptation
     * evidence only; it is never an issued confidence interval.
     */
    static final class PooledCrossfitMoments {
        final int sampleCount;
        final double effectiveSampleSize;
        final double effectiveSampleFraction;
        final Moments moments;

        private PooledCrossfitMoments(int sampleCount,
                                      double effectiveSampleSize,
                                      Moments moments) {
            this.sampleCount = sampleCount;
            this.effectiveSampleSize = effectiveSampleSize;
            this.effectiveSampleFraction =
                    effectiveSampleSize / sampleCount;
            this.moments = moments;
        }
    }

    /**
     * Pool normalized weighted moments from disjoint cross-fit folds.
     *
     * <p>Within fold {@code f}, weights sum to one and have
     * {@code sum(w^2)=1/ESS_f}.  Giving the fold mass {@code n_f/N} produces
     * exact pooled weighted mean, variance, and ESS identities without
     * pretending that either fold-specific eta was fitted on its validation
     * observations.</p>
     */
    static PooledCrossfitMoments poolCrossfitMoments(
            int[] sampleCounts,
            double[] effectiveSampleSizes,
            double[] bulkMeans,
            double[] bulkVariances,
            double[] tailProbabilities) {
        if (sampleCounts == null || sampleCounts.length == 0
                || effectiveSampleSizes == null
                || bulkMeans == null || bulkVariances == null
                || tailProbabilities == null
                || effectiveSampleSizes.length != sampleCounts.length
                || bulkMeans.length != sampleCounts.length
                || bulkVariances.length != sampleCounts.length
                || tailProbabilities.length != sampleCounts.length) {
            throw new IllegalArgumentException(
                    "cross-fit moment arrays must have one common positive length");
        }

        int total = 0;
        for (int count : sampleCounts) {
            if (count <= 1 || total > Integer.MAX_VALUE - count) {
                throw new IllegalArgumentException(
                        "cross-fit fold sample counts must exceed one and have a finite sum");
            }
            total += count;
        }

        double pooledMean = 0.0;
        double pooledSecond = 0.0;
        double pooledTail = 0.0;
        double pooledWeightSquares = 0.0;
        for (int fold = 0; fold < sampleCounts.length; fold++) {
            int count = sampleCounts[fold];
            double ess = effectiveSampleSizes[fold];
            double mean = bulkMeans[fold];
            double variance = bulkVariances[fold];
            double tail = tailProbabilities[fold];
            if (!Double.isFinite(ess) || !(ess > 0.0)
                    || ess > count * (1.0 + 1.0e-9)
                    || !Double.isFinite(mean) || mean < 0.0 || mean > 1.0
                    || !Double.isFinite(variance) || variance < 0.0
                    || !Double.isFinite(tail) || tail < 0.0 || tail > 1.0) {
                throw new IllegalArgumentException(
                        "invalid cross-fit weighted fold moments");
            }
            ess = Math.min(ess, count);
            double foldMass = (double) count / total;
            double foldWeightSquares = 1.0 / ess;
            double centeredNumerator =
                    variance * Math.max(0.0, 1.0 - foldWeightSquares);
            pooledMean += foldMass * mean;
            pooledSecond += foldMass
                    * (centeredNumerator + mean * mean);
            pooledTail += foldMass * tail;
            pooledWeightSquares += foldMass * foldMass
                    * foldWeightSquares;
        }
        if (!(pooledWeightSquares > 0.0)
                || !(pooledWeightSquares < 1.0)
                || !Double.isFinite(pooledWeightSquares)) {
            throw new IllegalArgumentException(
                    "cross-fit pooled weights have invalid concentration");
        }
        double varianceNumerator = Math.max(0.0,
                pooledSecond - pooledMean * pooledMean);
        double pooledVariance = varianceNumerator
                / (1.0 - pooledWeightSquares);
        Moments moments = new Moments(
                Math.max(0.0, Math.min(1.0, pooledMean)),
                pooledVariance,
                Math.max(0.0, Math.min(1.0, pooledTail)));
        return new PooledCrossfitMoments(
                total, 1.0 / pooledWeightSquares, moments);
    }

    /** A finite-sample interval in units of the frozen clipping scale. */
    public static final class Interval {
        public final int sampleCount;
        public final int tailCount;
        public final double bulkMean;
        public final double bulkVariance;
        public final double bulkRadius;
        public final double bulkLower;
        public final double bulkUpper;
        public final double tailProbabilityEmpirical;
        public final double tailProbabilityUpper;
        public final double conditionalSeverityCap;
        public final double inducedTailMeanUpper;
        public final double normalizedMeanLower;
        public final double normalizedMeanUpper;
        public final double epsilon;
        public final double empiricalTailMean;
        public final double empiricalConditionalSeverity;
        public final double observedMaxConditionalSeverity;

        private Interval(int sampleCount, int tailCount,
                         double bulkMean, double bulkVariance,
                         double bulkRadius, double bulkLower,
                         double bulkUpper,
                         double tailProbabilityEmpirical,
                         double tailProbabilityUpper,
                         double conditionalSeverityCap,
                         double inducedTailMeanUpper,
                         double normalizedMeanLower,
                         double normalizedMeanUpper,
                         double epsilon,
                         double empiricalTailMean,
                         double empiricalConditionalSeverity,
                         double observedMaxConditionalSeverity) {
            this.sampleCount = sampleCount;
            this.tailCount = tailCount;
            this.bulkMean = bulkMean;
            this.bulkVariance = bulkVariance;
            this.bulkRadius = bulkRadius;
            this.bulkLower = bulkLower;
            this.bulkUpper = bulkUpper;
            this.tailProbabilityEmpirical = tailProbabilityEmpirical;
            this.tailProbabilityUpper = tailProbabilityUpper;
            this.conditionalSeverityCap = conditionalSeverityCap;
            this.inducedTailMeanUpper = inducedTailMeanUpper;
            this.normalizedMeanLower = normalizedMeanLower;
            this.normalizedMeanUpper = normalizedMeanUpper;
            this.epsilon = epsilon;
            this.empiricalTailMean = empiricalTailMean;
            this.empiricalConditionalSeverity =
                    empiricalConditionalSeverity;
            this.observedMaxConditionalSeverity =
                    observedMaxConditionalSeverity;
        }

        public boolean hasPositiveBulkLower() {
            return bulkLower > 0.0;
        }

        public boolean isFinite() {
            return Double.isFinite(normalizedMeanUpper)
                    && normalizedMeanUpper > 0.0
                    && normalizedMeanLower >= 0.0
                    && normalizedMeanLower <= normalizedMeanUpper
                    && Double.isFinite(epsilon);
        }

        public Moments moments() {
            return new Moments(bulkMean, bulkVariance,
                    tailProbabilityEmpirical);
        }
    }

    /** Adaptation-only projected sample size.  It is not a certificate. */
    public static final class Sizing {
        public final int finalSamples;
        public final boolean reachableAtMax;
        public final double epsilonAtFinalSamples;
        public final double epsilonAtMaxSamples;
        public final double sizingTarget;

        private Sizing(int finalSamples, boolean reachableAtMax,
                       double epsilonAtFinalSamples,
                       double epsilonAtMaxSamples,
                       double sizingTarget) {
            this.finalSamples = finalSamples;
            this.reachableAtMax = reachableAtMax;
            this.epsilonAtFinalSamples = epsilonAtFinalSamples;
            this.epsilonAtMaxSamples = epsilonAtMaxSamples;
            this.sizingTarget = sizingTarget;
        }
    }

    /** Independent monitor result that can reject, but never validate, S0. */
    public static final class SeverityTest {
        public final int tailCount;
        public final boolean sufficientTailSamples;
        public final boolean logicalViolation;
        public final boolean rejected;
        public final double logEValue;
        public final double pValueUpper;

        private SeverityTest(int tailCount,
                             boolean sufficientTailSamples,
                             boolean logicalViolation,
                             boolean rejected,
                             double logEValue,
                             double pValueUpper) {
            this.tailCount = tailCount;
            this.sufficientTailSamples = sufficientTailSamples;
            this.logicalViolation = logicalViolation;
            this.rejected = rejected;
            this.logEValue = logEValue;
            this.pValueUpper = pValueUpper;
        }
    }

    /** Evaluate a fresh IID sample of log relative weights. */
    public static Interval evaluate(double[] logRelativeWeights,
                                    double logRelativeClip,
                                    double conditionalSeverityCap,
                                    double bulkDelta,
                                    double frequencyDelta) {
        requireLogWeights(logRelativeWeights);
        if (!Double.isFinite(logRelativeClip)) {
            throw new IllegalArgumentException(
                    "relative log clip must be finite");
        }
        requireSeverity(conditionalSeverityCap);
        requireDelta(bulkDelta, "bulk delta");
        requireDelta(frequencyDelta, "frequency delta");

        int n = logRelativeWeights.length;
        double[] bulk = new double[n];
        int tailCount = 0;
        double sumTail = 0.0;
        double maxTail = 0.0;
        boolean infiniteTail = false;

        for (int i = 0; i < n; i++) {
            double relative = logRelativeWeights[i] - logRelativeClip;
            if (relative > 0.0) {
                bulk[i] = 1.0;
                tailCount++;
                double excess = safeExpm1(relative);
                if (Double.isFinite(excess) && !infiniteTail) {
                    sumTail += excess;
                    if (!Double.isFinite(sumTail)) infiniteTail = true;
                } else {
                    infiniteTail = true;
                }
                maxTail = Math.max(maxTail, excess);
            } else {
                bulk[i] = Math.exp(relative);
            }
        }

        double mean = mean(bulk);
        double variance = sampleVariance(bulk, mean);
        double radius = empiricalBernsteinRadius(
                n, variance, 1.0, bulkDelta);
        double lower = Math.max(0.0, mean - radius);
        double upper = Math.min(1.0, mean + radius);
        double pHat = (double) tailCount / n;
        double pUpper = clopperPearsonUpper(
                tailCount, n, frequencyDelta);
        double tailUpper = pUpper * conditionalSeverityCap;
        double meanUpper = upper + tailUpper;
        double epsilon = intervalEpsilon(lower, meanUpper);
        double empiricalTailMean = infiniteTail
                ? Double.POSITIVE_INFINITY : sumTail / n;
        double empiricalSeverity = tailCount > 0
                ? empiricalTailMean / pHat : 0.0;

        return new Interval(
                n, tailCount, mean, variance, radius, lower, upper,
                pHat, pUpper, conditionalSeverityCap, tailUpper,
                lower, meanUpper, epsilon, empiricalTailMean,
                empiricalSeverity, maxTail);
    }

    /**
     * Project an interval from pilot/cross-fit point estimates.  This is a
     * design diagnostic only; only {@link #evaluate} on fresh IID samples is
     * the issued confidence interval.
     */
    public static Interval project(int sampleCount, Moments moments,
                                   double conditionalSeverityCap,
                                   double bulkDelta,
                                   double frequencyDelta) {
        if (sampleCount <= 1) {
            throw new IllegalArgumentException(
                    "projected sample count must exceed one");
        }
        if (moments == null) {
            throw new IllegalArgumentException("moments are required");
        }
        requireSeverity(conditionalSeverityCap);
        requireDelta(bulkDelta, "bulk delta");
        requireDelta(frequencyDelta, "frequency delta");

        double radius = empiricalBernsteinRadius(
                sampleCount, moments.bulkVariance, 1.0, bulkDelta);
        double lower = Math.max(0.0, moments.bulkMean - radius);
        double upper = Math.min(1.0, moments.bulkMean + radius);
        int pseudoTailCount = (int) Math.min(sampleCount,
                Math.ceil(moments.tailProbability * sampleCount));
        double pUpper = clopperPearsonUpper(
                pseudoTailCount, sampleCount, frequencyDelta);
        double tailUpper = pUpper * conditionalSeverityCap;
        double meanUpper = upper + tailUpper;
        double epsilon = intervalEpsilon(lower, meanUpper);

        return new Interval(
                sampleCount, pseudoTailCount, moments.bulkMean,
                moments.bulkVariance, radius, lower, upper,
                moments.tailProbability, pUpper,
                conditionalSeverityCap, tailUpper,
                lower, meanUpper, epsilon,
                Double.NaN, Double.NaN, Double.NaN);
    }

    /** Choose a pilot-frozen final N, including the unreachable diagnostic cap. */
    public static Sizing size(Moments moments,
                              int maxSamples,
                              int unreachableSamples,
                              double targetEpsilon,
                              double safetyFraction,
                              double conditionalSeverityCap,
                              double bulkDelta,
                              double frequencyDelta) {
        if (maxSamples <= 1 || unreachableSamples <= 1
                || unreachableSamples > maxSamples) {
            throw new IllegalArgumentException(
                    "invalid frequency/severity maximum/unreachable sample counts");
        }
        if (!Double.isFinite(targetEpsilon)
                || !(targetEpsilon > 0.0) || targetEpsilon >= 1.0) {
            throw new IllegalArgumentException(
                    "target epsilon must be finite and in (0,1)");
        }
        if (!Double.isFinite(safetyFraction)
                || !(safetyFraction > 0.0) || safetyFraction > 1.0) {
            throw new IllegalArgumentException(
                    "sizing safety fraction must be in (0,1]");
        }

        double sizingTarget = targetEpsilon * safetyFraction;
        Interval atMax = project(maxSamples, moments,
                conditionalSeverityCap, bulkDelta, frequencyDelta);
        if (atMax.epsilon > targetEpsilon) {
            Interval unreachable = project(unreachableSamples, moments,
                    conditionalSeverityCap, bulkDelta, frequencyDelta);
            return new Sizing(unreachableSamples, false,
                    unreachable.epsilon, atMax.epsilon, sizingTarget);
        }
        if (atMax.epsilon > sizingTarget) {
            return new Sizing(maxSamples, true, atMax.epsilon,
                    atMax.epsilon, sizingTarget);
        }

        int low = 2;
        int high = maxSamples;
        while (low < high) {
            int middle = (low + high) >>> 1;
            double epsilon = project(middle, moments,
                    conditionalSeverityCap, bulkDelta,
                    frequencyDelta).epsilon;
            if (epsilon <= sizingTarget) high = middle;
            else low = middle + 1;
        }
        Interval selected = project(low, moments,
                conditionalSeverityCap, bulkDelta, frequencyDelta);
        return new Sizing(low, true, selected.epsilon,
                atMax.epsilon, sizingTarget);
    }

    /** One-sided exact Clopper-Pearson upper confidence limit. */
    public static double clopperPearsonUpper(int successes, int trials,
                                             double delta) {
        if (trials <= 0 || successes < 0 || successes > trials) {
            throw new IllegalArgumentException(
                    "invalid binomial successes/trials");
        }
        requireDelta(delta, "binomial delta");
        if (successes == trials) return 1.0;
        if (successes == 0) {
            return -Math.expm1(Math.log(delta) / trials);
        }
        BetaDistribution beta = new BetaDistribution(
                successes + 1.0, trials - successes);
        double value = beta.inverseCumulativeProbability(1.0 - delta);
        return Math.max(0.0, Math.min(1.0, value));
    }

    /** Two-sided Maurer-Pontil radius for values in an interval of this range. */
    public static double empiricalBernsteinRadius(int sampleCount,
                                                  double sampleVariance,
                                                  double range,
                                                  double delta) {
        if (sampleCount <= 1) return Double.POSITIVE_INFINITY;
        if (!Double.isFinite(sampleVariance) || sampleVariance < 0.0
                || !Double.isFinite(range) || range < 0.0) {
            throw new IllegalArgumentException(
                    "invalid variance/range for empirical Bernstein");
        }
        requireDelta(delta, "empirical-Bernstein delta");
        double logTerm = Math.log(4.0 / delta);
        return Math.sqrt(2.0 * sampleVariance * logTerm / sampleCount)
                + 7.0 * range * logTerm
                / (3.0 * (sampleCount - 1));
    }

    /**
     * Test the conditional severity premise on an independent monitor batch.
     * Failure to reject is never interpreted as validation.
     */
    public static SeverityTest testConditionalSeverity(
            double[] logRelativeWeights,
            double logRelativeClip,
            double conditionalSeverityCap,
            double testAlpha) {
        requireLogWeights(logRelativeWeights);
        if (!Double.isFinite(logRelativeClip)) {
            throw new IllegalArgumentException(
                    "relative log clip must be finite");
        }
        requireSeverity(conditionalSeverityCap);
        requireDelta(testAlpha, "severity-test alpha");

        List<Double> logExcesses = new ArrayList<>();
        for (double logWeight : logRelativeWeights) {
            double relative = logWeight - logRelativeClip;
            if (relative <= 0.0) continue;
            if (relative > 40.0) {
                logExcesses.add(relative + Math.log1p(-Math.exp(-relative)));
            } else {
                logExcesses.add(Math.log(Math.expm1(relative)));
            }
        }
        int tailCount = logExcesses.size();
        if (conditionalSeverityCap == 0.0) {
            boolean violation = tailCount > 0;
            return new SeverityTest(tailCount, tailCount > 1,
                    violation, violation,
                    violation ? Double.POSITIVE_INFINITY : 0.0,
                    violation ? 0.0 : 1.0);
        }
        if (tailCount <= 1) {
            return new SeverityTest(tailCount, false, false,
                    false, Double.NaN, Double.NaN);
        }

        double logCap = Math.log(conditionalSeverityCap);
        List<Double> components = new ArrayList<>();
        components.add(logSumExp(logExcesses)
                - Math.log(tailCount) - logCap);
        for (double lambda : DEFAULT_BETTING_LAMBDAS) {
            double sum = 0.0;
            double logOneMinus = Math.log1p(-lambda);
            double logLambda = Math.log(lambda);
            for (double logExcess : logExcesses) {
                sum += logAddExp(logOneMinus,
                        logLambda + logExcess - logCap);
            }
            components.add(sum);
        }
        double logE = logSumExp(components) - Math.log(components.size());
        double p = logE > 0.0 ? Math.min(1.0, Math.exp(-logE)) : 1.0;
        boolean rejected = logE >= -Math.log(testAlpha);
        return new SeverityTest(tailCount, true, false,
                rejected, logE, p);
    }

    private static double intervalEpsilon(double lower, double upper) {
        if (!(upper > 0.0) || !Double.isFinite(upper)
                || lower < 0.0 || lower > upper) return 1.0;
        return Math.max(0.0, Math.min(1.0, 1.0 - lower / upper));
    }

    private static double mean(double[] values) {
        double sum = 0.0;
        for (double value : values) sum += value;
        return sum / values.length;
    }

    private static double sampleVariance(double[] values, double mean) {
        double sum = 0.0;
        for (double value : values) {
            double delta = value - mean;
            sum += delta * delta;
        }
        return Math.max(0.0, sum / (values.length - 1));
    }

    private static double safeExpm1(double value) {
        if (value > Math.log(Double.MAX_VALUE)) {
            return Double.POSITIVE_INFINITY;
        }
        return Math.expm1(value);
    }

    private static double logSumExp(List<Double> values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) max = Math.max(max, value);
        if (max == Double.POSITIVE_INFINITY) return max;
        if (max == Double.NEGATIVE_INFINITY) return max;
        double sum = 0.0;
        for (double value : values) sum += Math.exp(value - max);
        return max + Math.log(sum);
    }

    private static double logAddExp(double first, double second) {
        if (first == Double.NEGATIVE_INFINITY) return second;
        if (second == Double.NEGATIVE_INFINITY) return first;
        double max = Math.max(first, second);
        return max + Math.log1p(Math.exp(Math.min(first, second) - max));
    }

    private static void requireLogWeights(double[] logWeights) {
        if (logWeights == null || logWeights.length <= 1) {
            throw new IllegalArgumentException(
                    "at least two log relative weights are required");
        }
        for (double value : logWeights) {
            if (!Double.isFinite(value)) {
                throw new IllegalArgumentException(
                        "log relative weights must be finite");
            }
        }
    }

    private static void requireSeverity(double severity) {
        if (!Double.isFinite(severity) || severity < 0.0) {
            throw new IllegalArgumentException(
                    "conditional severity cap must be finite and nonnegative");
        }
    }

    private static void requireDelta(double delta, String label) {
        if (!Double.isFinite(delta) || !(delta > 0.0) || delta >= 1.0) {
            throw new IllegalArgumentException(label + " must be in (0,1)");
        }
    }

    private static void requireProbability(double value, String label) {
        if (!Double.isFinite(value) || value < 0.0 || value > 1.0) {
            throw new IllegalArgumentException(label + " must be in [0,1]");
        }
    }
}
