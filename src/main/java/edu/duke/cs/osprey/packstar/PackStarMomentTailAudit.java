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

import java.io.File;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Diagnostic-only audit of importance-weight second moments and tails.
 *
 * <p>The input samples must be fresh IID draws from the frozen proposal.  The
 * reported normalized second moment and ESS are invariant to a constant shift
 * in every log weight.  They are empirical diagnostics, not certified bounds
 * on an unseen population moment or tail.</p>
 */
final class PackStarMomentTailAudit {

    static final String ASSIGNMENT_SCHEMA =
            "packstar-moment-tail-samples-v1";
    static final String SUMMARY_SCHEMA =
            "packstar-moment-tail-summary-v1";

    static final class Sample {
        final int[] assignment;
        final double eMinKcal;
        final double eTrueKcal;
        final double logWeight;
        final double pointwiseOverCorrectionKcal;

        Sample(int[] assignment, double eMinKcal, double eTrueKcal,
               double logWeight, double pointwiseOverCorrectionKcal) {
            this.assignment = assignment.clone();
            this.eMinKcal = eMinKcal;
            this.eTrueKcal = eTrueKcal;
            this.logWeight = logWeight;
            this.pointwiseOverCorrectionKcal =
                    pointwiseOverCorrectionKcal;
        }
    }

    static final class Result {
        final List<Sample> samples;
        final double rt;
        final double[] thresholdsKcal;
        final double[] logWeightKcal;
        final int uniqueAssignments;
        final int lowerBoundViolations;
        final double logMeanWeight;
        final double logMeanSquaredWeight;
        final double normalizedSecondMoment;
        final double cvSquared;
        final double effectiveSampleSize;
        final double effectiveSampleSizeFraction;
        final double maxFirstMomentShare;
        final double maxSecondMomentShare;
        final double prefixQuarterNormalizedSecondMoment;
        final double prefixHalfNormalizedSecondMoment;
        final double batchMinNormalizedSecondMoment;
        final double batchMedianNormalizedSecondMoment;
        final double batchMaxNormalizedSecondMoment;
        final int[] exceedCounts;
        final double[] firstMomentTailFractions;
        final double[] secondMomentTailFractions;

        Result(List<Sample> samples, double rt,
               double[] thresholdsKcal, double[] logWeightKcal,
               int uniqueAssignments, int lowerBoundViolations,
               double logMeanWeight, double logMeanSquaredWeight,
               double normalizedSecondMoment, double cvSquared,
               double effectiveSampleSize,
               double effectiveSampleSizeFraction,
               double maxFirstMomentShare,
               double maxSecondMomentShare,
               double prefixQuarterNormalizedSecondMoment,
               double prefixHalfNormalizedSecondMoment,
               double batchMinNormalizedSecondMoment,
               double batchMedianNormalizedSecondMoment,
               double batchMaxNormalizedSecondMoment,
               int[] exceedCounts,
               double[] firstMomentTailFractions,
               double[] secondMomentTailFractions) {
            this.samples = samples;
            this.rt = rt;
            this.thresholdsKcal = thresholdsKcal;
            this.logWeightKcal = logWeightKcal;
            this.uniqueAssignments = uniqueAssignments;
            this.lowerBoundViolations = lowerBoundViolations;
            this.logMeanWeight = logMeanWeight;
            this.logMeanSquaredWeight = logMeanSquaredWeight;
            this.normalizedSecondMoment = normalizedSecondMoment;
            this.cvSquared = cvSquared;
            this.effectiveSampleSize = effectiveSampleSize;
            this.effectiveSampleSizeFraction =
                    effectiveSampleSizeFraction;
            this.maxFirstMomentShare = maxFirstMomentShare;
            this.maxSecondMomentShare = maxSecondMomentShare;
            this.prefixQuarterNormalizedSecondMoment =
                    prefixQuarterNormalizedSecondMoment;
            this.prefixHalfNormalizedSecondMoment =
                    prefixHalfNormalizedSecondMoment;
            this.batchMinNormalizedSecondMoment =
                    batchMinNormalizedSecondMoment;
            this.batchMedianNormalizedSecondMoment =
                    batchMedianNormalizedSecondMoment;
            this.batchMaxNormalizedSecondMoment =
                    batchMaxNormalizedSecondMoment;
            this.exceedCounts = exceedCounts;
            this.firstMomentTailFractions = firstMomentTailFractions;
            this.secondMomentTailFractions = secondMomentTailFractions;
        }
    }

    private PackStarMomentTailAudit() {}

    static Result analyze(List<Sample> inputSamples,
                          double[] inputThresholdsKcal,
                          double rt) {
        if (inputSamples == null || inputSamples.size() < 2) {
            throw new IllegalArgumentException(
                    "moment/tail audit requires at least two samples");
        }
        if (!Double.isFinite(rt) || !(rt > 0.0)) {
            throw new IllegalArgumentException(
                    "moment/tail audit RT must be finite and positive");
        }
        double[] thresholds = inputThresholdsKcal.clone();
        if (thresholds.length == 0) {
            throw new IllegalArgumentException(
                    "moment/tail threshold grid is empty");
        }
        for (int i = 0; i < thresholds.length; i++) {
            if (!Double.isFinite(thresholds[i])
                    || (i > 0 && !(thresholds[i] > thresholds[i - 1]))) {
                throw new IllegalArgumentException(
                        "moment/tail thresholds must be finite and strictly increasing");
            }
        }

        int n = inputSamples.size();
        List<Sample> samples = new ArrayList<>(n);
        double[] logWeights = new double[n];
        double[] logWeightKcal = new double[n];
        Set<String> assignments = new HashSet<>();
        int lowerBoundViolations = 0;
        for (int i = 0; i < n; i++) {
            Sample sample = inputSamples.get(i);
            if (sample == null || sample.assignment.length == 0
                    || !Double.isFinite(sample.eMinKcal)
                    || !Double.isFinite(sample.eTrueKcal)
                    || !Double.isFinite(sample.logWeight)
                    || !Double.isFinite(
                    sample.pointwiseOverCorrectionKcal)) {
                throw new IllegalArgumentException(
                        "moment/tail sample " + i
                                + " is incomplete or non-finite");
            }
            samples.add(new Sample(sample.assignment, sample.eMinKcal,
                    sample.eTrueKcal, sample.logWeight,
                    sample.pointwiseOverCorrectionKcal));
            logWeights[i] = sample.logWeight;
            logWeightKcal[i] = rt * sample.logWeight;
            assignments.add(assignmentText(sample.assignment));
            if (sample.eMinKcal > sample.eTrueKcal) {
                lowerBoundViolations++;
            }
        }

        double logSumW = logSumExp(logWeights, 1.0, 0, n);
        double logSumW2 = logSumExp(logWeights, 2.0, 0, n);
        double logN = Math.log(n);
        double logMeanWeight = logSumW - logN;
        double logMeanSquaredWeight = logSumW2 - logN;
        double normalizedSecondMoment = Math.exp(
                logN + logSumW2 - 2.0 * logSumW);
        double effectiveSampleSize = Math.exp(
                2.0 * logSumW - logSumW2);
        double maxLogWeight = Arrays.stream(logWeights)
                .max().orElseThrow(IllegalStateException::new);
        double maxFirstMomentShare = Math.exp(
                maxLogWeight - logSumW);
        double maxSecondMomentShare = Math.exp(
                2.0 * maxLogWeight - logSumW2);

        int quarter = Math.max(2, n / 4);
        int half = Math.max(2, n / 2);
        double prefixQuarter = normalizedSecondMoment(
                logWeights, 0, quarter);
        double prefixHalf = normalizedSecondMoment(
                logWeights, 0, half);

        int batchCount = Math.min(8, n / 2);
        double[] batchMoments = new double[batchCount];
        for (int batch = 0; batch < batchCount; batch++) {
            int from = batch * n / batchCount;
            int to = (batch + 1) * n / batchCount;
            batchMoments[batch] = normalizedSecondMoment(
                    logWeights, from, to);
        }
        Arrays.sort(batchMoments);

        int[] exceedCounts = new int[thresholds.length];
        double[] logTailW = new double[thresholds.length];
        double[] logTailW2 = new double[thresholds.length];
        Arrays.fill(logTailW, Double.NEGATIVE_INFINITY);
        Arrays.fill(logTailW2, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < thresholds.length; j++) {
                if (logWeightKcal[i] > thresholds[j]) {
                    exceedCounts[j]++;
                    logTailW[j] = logAddExp(
                            logTailW[j], logWeights[i]);
                    logTailW2[j] = logAddExp(
                            logTailW2[j], 2.0 * logWeights[i]);
                }
            }
        }
        double[] firstTail = new double[thresholds.length];
        double[] secondTail = new double[thresholds.length];
        for (int j = 0; j < thresholds.length; j++) {
            firstTail[j] = tailFraction(logTailW[j], logSumW);
            secondTail[j] = tailFraction(logTailW2[j], logSumW2);
        }

        return new Result(samples, rt, thresholds, logWeightKcal,
                assignments.size(), lowerBoundViolations,
                logMeanWeight, logMeanSquaredWeight,
                normalizedSecondMoment,
                Math.max(0.0, normalizedSecondMoment - 1.0),
                effectiveSampleSize, effectiveSampleSize / n,
                maxFirstMomentShare, maxSecondMomentShare,
                prefixQuarter, prefixHalf,
                batchMoments[0], quantileSorted(batchMoments, 0.5),
                batchMoments[batchMoments.length - 1],
                exceedCounts, firstTail, secondTail);
    }

    static void writeAssignmentTsv(File output, Result result,
                                   String scheme, String stage)
            throws Exception {
        refuseInvalidOutput(output);
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tscheme\tstage\tsampleIndex"
                    + "\tassignment\tE_m_kcal\tE_true_kcal"
                    + "\tlogWeight\tlogWeightEquivalentKcal"
                    + "\tpointwiseOverCorrectionKcal"
                    + "\tlowerBoundViolation");
            for (int i = 0; i < result.samples.size(); i++) {
                Sample sample = result.samples.get(i);
                writer.printf(Locale.ROOT,
                        "%s\t%s\t%s\t%d\t%s"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%s%n",
                        ASSIGNMENT_SCHEMA, sanitize(scheme),
                        sanitize(stage), i,
                        assignmentText(sample.assignment),
                        sample.eMinKcal, sample.eTrueKcal,
                        sample.logWeight, result.logWeightKcal[i],
                        sample.pointwiseOverCorrectionKcal,
                        Boolean.toString(
                                sample.eMinKcal > sample.eTrueKcal));
            }
        }
    }

    static void writeSummaryTsv(File output, Result result,
                                String scheme, String stage)
            throws Exception {
        refuseInvalidOutput(output);
        double qMin = quantile(result.logWeightKcal, 0.0);
        double q50 = quantile(result.logWeightKcal, 0.5);
        double q90 = quantile(result.logWeightKcal, 0.9);
        double q95 = quantile(result.logWeightKcal, 0.95);
        double q99 = quantile(result.logWeightKcal, 0.99);
        double q999 = quantile(result.logWeightKcal, 0.999);
        double qMax = quantile(result.logWeightKcal, 1.0);
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tscheme\tstage\tguarantee"
                    + "\tsamples\tuniqueAssignments\tRT"
                    + "\tlowerBoundViolations"
                    + "\tlogMeanWeight\tlogMeanSquaredWeight"
                    + "\tnormalizedSecondMoment\tcvSquared"
                    + "\teffectiveSampleSize\teffectiveSampleSizeFraction"
                    + "\tmaxFirstMomentShare\tmaxSecondMomentShare"
                    + "\tprefixQuarterNormalizedSecondMoment"
                    + "\tprefixHalfNormalizedSecondMoment"
                    + "\tbatchMinNormalizedSecondMoment"
                    + "\tbatchMedianNormalizedSecondMoment"
                    + "\tbatchMaxNormalizedSecondMoment"
                    + "\tlogWeightKcalMin\tlogWeightKcalP50"
                    + "\tlogWeightKcalP90\tlogWeightKcalP95"
                    + "\tlogWeightKcalP99\tlogWeightKcalP999"
                    + "\tlogWeightKcalMax\tthresholdKcal"
                    + "\texceedCount\texceedFraction"
                    + "\tfirstMomentTailFraction"
                    + "\tsecondMomentTailFraction");
            for (int j = 0; j < result.thresholdsKcal.length; j++) {
                writer.printf(Locale.ROOT,
                        "%s\t%s\t%s\tempirical-fresh-IID-diagnostic-only"
                                + "\t%d\t%d\t%.17g\t%d"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%d\t%.17g\t%.17g\t%.17g%n",
                        SUMMARY_SCHEMA, sanitize(scheme), sanitize(stage),
                        result.samples.size(), result.uniqueAssignments,
                        result.rt, result.lowerBoundViolations,
                        result.logMeanWeight,
                        result.logMeanSquaredWeight,
                        result.normalizedSecondMoment,
                        result.cvSquared,
                        result.effectiveSampleSize,
                        result.effectiveSampleSizeFraction,
                        result.maxFirstMomentShare,
                        result.maxSecondMomentShare,
                        result.prefixQuarterNormalizedSecondMoment,
                        result.prefixHalfNormalizedSecondMoment,
                        result.batchMinNormalizedSecondMoment,
                        result.batchMedianNormalizedSecondMoment,
                        result.batchMaxNormalizedSecondMoment,
                        qMin, q50, q90, q95, q99, q999, qMax,
                        result.thresholdsKcal[j],
                        result.exceedCounts[j],
                        result.exceedCounts[j]
                                / (double) result.samples.size(),
                        result.firstMomentTailFractions[j],
                        result.secondMomentTailFractions[j]);
            }
        }
    }

    private static double normalizedSecondMoment(
            double[] logWeights, int from, int to) {
        int n = to - from;
        if (n < 2) return Double.NaN;
        double logSum = logSumExp(logWeights, 1.0, from, to);
        double logSum2 = logSumExp(logWeights, 2.0, from, to);
        return Math.exp(Math.log(n) + logSum2 - 2.0 * logSum);
    }

    private static double logSumExp(double[] values, double scale,
                                    int from, int to) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = from; i < to; i++) {
            max = Math.max(max, scale * values[i]);
        }
        double sum = 0.0;
        for (int i = from; i < to; i++) {
            sum += Math.exp(scale * values[i] - max);
        }
        return max + Math.log(sum);
    }

    private static double logAddExp(double first, double second) {
        if (first == Double.NEGATIVE_INFINITY) return second;
        if (second == Double.NEGATIVE_INFINITY) return first;
        double max = Math.max(first, second);
        return max + Math.log1p(Math.exp(Math.min(first, second) - max));
    }

    private static double tailFraction(double logTail, double logTotal) {
        return logTail == Double.NEGATIVE_INFINITY
                ? 0.0 : Math.exp(logTail - logTotal);
    }

    private static double quantile(double[] values, double q) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        return quantileSorted(sorted, q);
    }

    private static double quantileSorted(double[] sorted, double q) {
        double position = Math.max(0.0, Math.min(1.0, q))
                * (sorted.length - 1);
        int low = (int) Math.floor(position);
        int high = (int) Math.ceil(position);
        if (low == high) return sorted[low];
        double fraction = position - low;
        return sorted[low] * (1.0 - fraction)
                + sorted[high] * fraction;
    }

    private static String assignmentText(int[] assignment) {
        StringBuilder out = new StringBuilder();
        for (int i = 0; i < assignment.length; i++) {
            if (i > 0) out.append(',');
            out.append(assignment[i]);
        }
        return out.toString();
    }

    private static String sanitize(String value) {
        return value == null ? ""
                : value.replace('\t', ' ')
                .replace('\n', ' ').replace('\r', ' ');
    }

    private static void refuseInvalidOutput(File output) {
        File parent = output.getParentFile();
        if (parent == null || !parent.isDirectory()) {
            throw new IllegalArgumentException(
                    "moment/tail output parent does not exist: "
                            + (parent == null ? "null" : parent));
        }
        if (output.exists()) {
            throw new IllegalArgumentException(
                    "moment/tail audit refuses to overwrite output: "
                            + output);
        }
    }
}
