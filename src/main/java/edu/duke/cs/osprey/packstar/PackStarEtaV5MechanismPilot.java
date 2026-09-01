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
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Diagnostic-only analysis for the eta-v5 protein/10 mechanism pilot.
 *
 * Samples are IID from p_m.  The analysis asks whether the canonical
 * importance residual eta-g actually exceeds each frozen threshold, or
 * whether the exact Phase-A set {eta > tau} was mostly a conservative false
 * positive caused by the weak lower bound g >= 0.
 */
final class PackStarEtaV5MechanismPilot {

    static final String ASSIGNMENT_SCHEMA =
            "packstar-eta-v5-pm-mechanism-samples-v1";
    static final String SUMMARY_SCHEMA =
            "packstar-eta-v5-pm-mechanism-summary-v1";

    static final class Sample {
        final int[] assignment;
        final double eMinKcal;
        final double eProposalRawKcal;
        final double eTrueKcal;

        Sample(int[] assignment, double eMinKcal,
               double eProposalRawKcal, double eTrueKcal) {
            this.assignment = assignment.clone();
            this.eMinKcal = eMinKcal;
            this.eProposalRawKcal = eProposalRawKcal;
            this.eTrueKcal = eTrueKcal;
        }
    }

    static final class Result {
        final List<Sample> samples;
        final double canonicalShiftKcal;
        final double[] thresholdsKcal;
        final double[] gapKcal;
        final double[] rawEtaKcal;
        final double[] canonicalEtaKcal;
        final double[] rawProposalErrorKcal;
        final double[] canonicalResidualKcal;
        final boolean[][] exceeds;
        final double[][] h;
        final int[] exceedCounts;
        final double[] meanH;
        final double[] varianceH;
        final double[] maxH;
        final double pearsonRawEtaGap;
        final boolean pearsonRawEtaGapDefined;
        final double spearmanRawEtaGap;
        final boolean spearmanRawEtaGapDefined;

        Result(List<Sample> samples, double canonicalShiftKcal,
               double[] thresholdsKcal, double[] gapKcal,
               double[] rawEtaKcal, double[] canonicalEtaKcal,
               double[] rawProposalErrorKcal,
               double[] canonicalResidualKcal, boolean[][] exceeds,
               double[][] h, int[] exceedCounts, double[] meanH,
               double[] varianceH, double[] maxH,
               double pearsonRawEtaGap,
               boolean pearsonRawEtaGapDefined,
               double spearmanRawEtaGap,
               boolean spearmanRawEtaGapDefined) {
            this.samples = samples;
            this.canonicalShiftKcal = canonicalShiftKcal;
            this.thresholdsKcal = thresholdsKcal;
            this.gapKcal = gapKcal;
            this.rawEtaKcal = rawEtaKcal;
            this.canonicalEtaKcal = canonicalEtaKcal;
            this.rawProposalErrorKcal = rawProposalErrorKcal;
            this.canonicalResidualKcal = canonicalResidualKcal;
            this.exceeds = exceeds;
            this.h = h;
            this.exceedCounts = exceedCounts;
            this.meanH = meanH;
            this.varianceH = varianceH;
            this.maxH = maxH;
            this.pearsonRawEtaGap = pearsonRawEtaGap;
            this.pearsonRawEtaGapDefined = pearsonRawEtaGapDefined;
            this.spearmanRawEtaGap = spearmanRawEtaGap;
            this.spearmanRawEtaGapDefined = spearmanRawEtaGapDefined;
        }
    }

    private PackStarEtaV5MechanismPilot() {}

    static Result analyze(List<Sample> inputSamples,
                          double canonicalShiftKcal,
                          double[] inputThresholdsKcal,
                          double rt) {
        if (inputSamples == null || inputSamples.size() < 2) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism pilot requires at least two samples");
        }
        if (!Double.isFinite(canonicalShiftKcal)) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism pilot canonical shift must be finite");
        }
        if (!Double.isFinite(rt) || !(rt > 0.0)) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism pilot RT must be finite and positive");
        }
        double[] thresholdsKcal = inputThresholdsKcal.clone();
        if (thresholdsKcal.length == 0) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism pilot threshold grid is empty");
        }
        for (int i = 0; i < thresholdsKcal.length; i++) {
            if (!Double.isFinite(thresholdsKcal[i])
                    || (i > 0 && !(thresholdsKcal[i]
                    > thresholdsKcal[i - 1]))) {
                throw new IllegalArgumentException(
                        "eta-v5 mechanism pilot thresholds must be finite"
                                + " and strictly increasing");
            }
        }

        List<Sample> samples = new ArrayList<>(inputSamples.size());
        int n = inputSamples.size();
        int k = thresholdsKcal.length;
        double[] gapKcal = new double[n];
        double[] rawEtaKcal = new double[n];
        double[] canonicalEtaKcal = new double[n];
        double[] rawProposalErrorKcal = new double[n];
        double[] canonicalResidualKcal = new double[n];
        boolean[][] exceeds = new boolean[n][k];
        double[][] h = new double[n][k];
        int[] exceedCounts = new int[k];
        double[] sumH = new double[k];
        double[] sumH2 = new double[k];
        double[] maxH = new double[k];

        for (int i = 0; i < n; i++) {
            Sample sample = inputSamples.get(i);
            if (sample == null || sample.assignment.length == 0
                    || !Double.isFinite(sample.eMinKcal)
                    || !Double.isFinite(sample.eProposalRawKcal)
                    || !Double.isFinite(sample.eTrueKcal)) {
                throw new IllegalArgumentException(
                        "eta-v5 mechanism pilot sample " + i
                                + " is incomplete or non-finite");
            }
            if (sample.eMinKcal > sample.eTrueKcal) {
                throw new IllegalArgumentException(
                        "eta-v5 mechanism pilot lower-bound violation at"
                                + " sample " + i + ": E_m="
                                + sample.eMinKcal + " E_true="
                                + sample.eTrueKcal);
            }
            samples.add(new Sample(sample.assignment, sample.eMinKcal,
                    sample.eProposalRawKcal, sample.eTrueKcal));
            gapKcal[i] = sample.eTrueKcal - sample.eMinKcal;
            rawEtaKcal[i] = sample.eProposalRawKcal - sample.eMinKcal;
            canonicalEtaKcal[i] = rawEtaKcal[i] - canonicalShiftKcal;
            rawProposalErrorKcal[i] =
                    sample.eProposalRawKcal - sample.eTrueKcal;
            canonicalResidualKcal[i] =
                    canonicalEtaKcal[i] - gapKcal[i];

            for (int j = 0; j < k; j++) {
                exceeds[i][j] = canonicalResidualKcal[i]
                        > thresholdsKcal[j];
                if (exceeds[i][j]) exceedCounts[j]++;
                double value = PackStarEstimator.tailRemainderH(
                        gapKcal[i], canonicalEtaKcal[i], rt,
                        thresholdsKcal[j] / rt);
                if (!Double.isFinite(value) || value < 0.0
                        || value > 1.0) {
                    throw new IllegalArgumentException(
                            "eta-v5 mechanism pilot invalid h at sample "
                                    + i + " threshold "
                                    + thresholdsKcal[j] + ": " + value);
                }
                h[i][j] = value;
                sumH[j] += value;
                sumH2[j] += value * value;
                maxH[j] = Math.max(maxH[j], value);
            }
        }

        double[] meanH = new double[k];
        double[] varianceH = new double[k];
        for (int j = 0; j < k; j++) {
            meanH[j] = sumH[j] / n;
            varianceH[j] = Math.max(0.0,
                    (sumH2[j] - n * meanH[j] * meanH[j]) / (n - 1));
        }

        double rawPearson = pearson(rawEtaKcal, gapKcal);
        double rawSpearman = spearman(rawEtaKcal, gapKcal);
        boolean pearsonRawEtaGapDefined = Double.isFinite(rawPearson);
        boolean spearmanRawEtaGapDefined = Double.isFinite(rawSpearman);
        // A highly concentrated p_m can return the same assignment in all
        // 200 IID draws.  That makes correlation mathematically undefined,
        // but the residual/exceedance diagnostic remains fully informative.
        // Preserve that distinction explicitly instead of aborting the run.
        double pearsonRawEtaGap = pearsonRawEtaGapDefined
                ? rawPearson : 0.0;
        double spearmanRawEtaGap = spearmanRawEtaGapDefined
                ? rawSpearman : 0.0;

        return new Result(samples, canonicalShiftKcal, thresholdsKcal,
                gapKcal, rawEtaKcal, canonicalEtaKcal,
                rawProposalErrorKcal, canonicalResidualKcal,
                exceeds, h, exceedCounts, meanH, varianceH, maxH,
                pearsonRawEtaGap, pearsonRawEtaGapDefined,
                spearmanRawEtaGap, spearmanRawEtaGapDefined);
    }

    static void writeAssignmentTsv(File output, Result result,
                                   String proposalSource,
                                   String randomStreamIdentity)
            throws Exception {
        refuseInvalidOutput(output);
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.print("schema\tsampleIndex\tproposalSource"
                    + "\trandomStreamIdentity\tassignment"
                    + "\tcanonicalShiftKcal\tE_m_kcal"
                    + "\tE_proposal_raw_kcal"
                    + "\tE_proposal_canonical_kcal\tE_true_kcal"
                    + "\tg_kcal\trawEtaKcal\tcanonicalEtaKcal"
                    + "\trawProposalErrorKcal\tcanonicalResidualKcal");
            for (double threshold : result.thresholdsKcal) {
                String token = thresholdToken(threshold);
                writer.print("\texceed_" + token + "\th_" + token);
            }
            writer.println();

            for (int i = 0; i < result.samples.size(); i++) {
                Sample sample = result.samples.get(i);
                writer.printf(Locale.ROOT,
                        "%s\t%d\t%s\t%s\t%s"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g",
                        ASSIGNMENT_SCHEMA, i,
                        sanitize(proposalSource),
                        sanitize(randomStreamIdentity),
                        assignmentText(sample.assignment),
                        result.canonicalShiftKcal,
                        sample.eMinKcal,
                        sample.eProposalRawKcal,
                        sample.eProposalRawKcal
                                - result.canonicalShiftKcal,
                        sample.eTrueKcal,
                        result.gapKcal[i],
                        result.rawEtaKcal[i],
                        result.canonicalEtaKcal[i],
                        result.rawProposalErrorKcal[i],
                        result.canonicalResidualKcal[i]);
                for (int j = 0; j < result.thresholdsKcal.length; j++) {
                    writer.printf(Locale.ROOT, "\t%s\t%.17g",
                            Boolean.toString(result.exceeds[i][j]),
                            result.h[i][j]);
                }
                writer.println();
            }
        }
    }

    static void writeSummaryTsv(File output, Result result,
                                String proposalSource,
                                String randomStreamIdentity)
            throws Exception {
        refuseInvalidOutput(output);
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tproposalSource\trandomStreamIdentity"
                    + "\tthresholdGauge\tsamples\tcanonicalShiftKcal"
                    + "\tcanonicalThresholdKcal\tactualExceedCount"
                    + "\tactualExceedFraction\tmeanH\tvarianceH\tmaxH"
                    + "\tmeanRawEtaKcal\tmeanGapKcal"
                    + "\tpearsonRawEtaGap\tpearsonRawEtaGapDefined"
                    + "\tspearmanRawEtaGap\tspearmanRawEtaGapDefined"
                    + "\trawProposalErrorMinKcal"
                    + "\trawProposalErrorMedianKcal"
                    + "\trawProposalErrorP90Kcal"
                    + "\trawProposalErrorP95Kcal"
                    + "\trawProposalErrorMaxKcal"
                    + "\tcanonicalResidualMinKcal"
                    + "\tcanonicalResidualMedianKcal"
                    + "\tcanonicalResidualP90Kcal"
                    + "\tcanonicalResidualP95Kcal"
                    + "\tcanonicalResidualMaxKcal"
                    + "\tlowerBoundViolations");

            double meanRawEta = mean(result.rawEtaKcal);
            double meanGap = mean(result.gapKcal);
            double rawErrorMin = quantile(result.rawProposalErrorKcal, 0.0);
            double rawErrorMedian = quantile(
                    result.rawProposalErrorKcal, 0.5);
            double rawErrorP90 = quantile(
                    result.rawProposalErrorKcal, 0.9);
            double rawErrorP95 = quantile(
                    result.rawProposalErrorKcal, 0.95);
            double rawErrorMax = quantile(result.rawProposalErrorKcal, 1.0);
            double residualMin = quantile(
                    result.canonicalResidualKcal, 0.0);
            double residualMedian = quantile(
                    result.canonicalResidualKcal, 0.5);
            double residualP90 = quantile(
                    result.canonicalResidualKcal, 0.9);
            double residualP95 = quantile(
                    result.canonicalResidualKcal, 0.95);
            double residualMax = quantile(
                    result.canonicalResidualKcal, 1.0);

            for (int j = 0; j < result.thresholdsKcal.length; j++) {
                writer.printf(Locale.ROOT,
                        "%s\t%s\t%s\tcanonical-normalizer"
                                + "\t%d\t%.17g\t%.17g\t%d\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%s\t%.17g\t%s"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t0%n",
                        SUMMARY_SCHEMA,
                        sanitize(proposalSource),
                        sanitize(randomStreamIdentity),
                        result.samples.size(),
                        result.canonicalShiftKcal,
                        result.thresholdsKcal[j],
                        result.exceedCounts[j],
                        result.exceedCounts[j]
                                / (double) result.samples.size(),
                        result.meanH[j],
                        result.varianceH[j],
                        result.maxH[j],
                        meanRawEta,
                        meanGap,
                        result.pearsonRawEtaGap,
                        Boolean.toString(
                                result.pearsonRawEtaGapDefined),
                        result.spearmanRawEtaGap,
                        Boolean.toString(
                                result.spearmanRawEtaGapDefined),
                        rawErrorMin,
                        rawErrorMedian,
                        rawErrorP90,
                        rawErrorP95,
                        rawErrorMax,
                        residualMin,
                        residualMedian,
                        residualP90,
                        residualP95,
                        residualMax);
            }
        }
    }

    private static void refuseInvalidOutput(File output) {
        File parent = output.getParentFile();
        if (parent == null || !parent.isDirectory()) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism output parent does not exist: "
                            + (parent == null ? "null" : parent));
        }
        if (output.exists()) {
            throw new IllegalArgumentException(
                    "eta-v5 mechanism pilot refuses to overwrite output: "
                            + output);
        }
    }

    static double pearson(double[] x, double[] y) {
        if (x.length != y.length || x.length < 2) return Double.NaN;
        double meanX = mean(x);
        double meanY = mean(y);
        double covariance = 0.0;
        double varianceX = 0.0;
        double varianceY = 0.0;
        for (int i = 0; i < x.length; i++) {
            double dx = x[i] - meanX;
            double dy = y[i] - meanY;
            covariance += dx * dy;
            varianceX += dx * dx;
            varianceY += dy * dy;
        }
        if (!(varianceX > 0.0) || !(varianceY > 0.0)) {
            return Double.NaN;
        }
        return clampCorrelation(covariance
                / Math.sqrt(varianceX * varianceY));
    }

    static double spearman(double[] x, double[] y) {
        if (x.length != y.length || x.length < 2) return Double.NaN;
        return pearson(ranks(x), ranks(y));
    }

    private static double[] ranks(double[] values) {
        Integer[] order = new Integer[values.length];
        for (int i = 0; i < values.length; i++) order[i] = i;
        Arrays.sort(order, Comparator.comparingDouble(i -> values[i]));
        double[] ranks = new double[values.length];
        int start = 0;
        while (start < order.length) {
            int end = start + 1;
            while (end < order.length
                    && Double.compare(values[order[start]],
                    values[order[end]]) == 0) {
                end++;
            }
            double averageRank = 0.5 * (start + end - 1) + 1.0;
            for (int i = start; i < end; i++) {
                ranks[order[i]] = averageRank;
            }
            start = end;
        }
        return ranks;
    }

    private static double mean(double[] values) {
        double sum = 0.0;
        for (double value : values) sum += value;
        return sum / values.length;
    }

    private static double quantile(double[] values, double q) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        double position = Math.max(0.0, Math.min(1.0, q))
                * (sorted.length - 1);
        int low = (int) Math.floor(position);
        int high = (int) Math.ceil(position);
        if (low == high) return sorted[low];
        double fraction = position - low;
        return sorted[low] * (1.0 - fraction)
                + sorted[high] * fraction;
    }

    private static double clampCorrelation(double value) {
        if (value > 1.0 && value < 1.0 + 1.0e-12) return 1.0;
        if (value < -1.0 && value > -1.0 - 1.0e-12) return -1.0;
        return value;
    }

    private static String assignmentText(int[] assignment) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < assignment.length; i++) {
            if (i > 0) builder.append(',');
            builder.append(assignment[i]);
        }
        return builder.toString();
    }

    private static String thresholdToken(double threshold) {
        String text = BigDecimalText.format(threshold);
        return text.replace("-", "m")
                .replace("+", "p")
                .replace(".", "p");
    }

    private static String sanitize(String value) {
        if (value == null) return "";
        return value.replace('\t', ' ')
                .replace('\n', ' ')
                .replace('\r', ' ');
    }

    private static final class BigDecimalText {
        static String format(double value) {
            return java.math.BigDecimal.valueOf(value)
                    .stripTrailingZeros().toPlainString();
        }
    }
}
