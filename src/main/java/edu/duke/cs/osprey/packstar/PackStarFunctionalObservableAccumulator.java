package edu.duke.cs.osprey.packstar;

import java.util.ArrayList;
import java.util.List;

/** Numerically stable accumulator for a self-normalized event probability. */
final class PackStarFunctionalObservableAccumulator {

    private final String eventName;
    private final double rt;
    private final List<Double> finiteLogWeights = new ArrayList<>();
    private int sampleCount;
    private int invalidSampleCount;
    private int hitCount;
    private final List<Boolean> finiteHits = new ArrayList<>();

    PackStarFunctionalObservableAccumulator(String eventName, double rt) {
        if (!Double.isFinite(rt) || !(rt > 0.0)) {
            throw new IllegalArgumentException("RT must be finite and positive");
        }
        this.eventName = eventName;
        this.rt = rt;
    }

    void add(boolean hit, double logWeight) {
        sampleCount++;
        if (hit) {
            hitCount++;
        }
        if (Double.isNaN(logWeight) || logWeight == Double.POSITIVE_INFINITY) {
            invalidSampleCount++;
            return;
        }
        if (Double.isFinite(logWeight)) {
            finiteLogWeights.add(logWeight);
            finiteHits.add(hit);
        }
    }

    PackStarFunctionalObservableResult finish() {
        int weightedSampleCount = finiteLogWeights.size();
        double rawHitFraction = sampleCount > 0
                ? ((double) hitCount) / sampleCount : Double.NaN;
        if (weightedSampleCount == 0) {
            return PackStarFunctionalObservableResult.fromSamples(
                    eventName,
                    PackStarFunctionalObservableResult.Status.UNRESOLVED_NO_FINITE_WEIGHT,
                    sampleCount, 0, invalidSampleCount, hitCount,
                    Double.NaN, Double.NaN, 0.0, Double.NaN,
                    rawHitFraction, rt,
                    "no sample had a finite nonzero importance weight");
        }

        double maxLogWeight = Double.NEGATIVE_INFINITY;
        for (double logWeight : finiteLogWeights) {
            maxLogWeight = Math.max(maxLogWeight, logWeight);
        }
        double totalScaledWeight = 0.0;
        double hitScaledWeight = 0.0;
        double squaredScaledWeight = 0.0;
        for (int i = 0; i < finiteLogWeights.size(); i++) {
            double scaled = Math.exp(finiteLogWeights.get(i) - maxLogWeight);
            totalScaledWeight += scaled;
            squaredScaledWeight += scaled * scaled;
            if (finiteHits.get(i)) {
                hitScaledWeight += scaled;
            }
        }

        double probability = hitScaledWeight / totalScaledWeight;
        probability = Math.max(0.0, Math.min(1.0, probability));
        double effectiveSampleSize = squaredScaledWeight > 0.0
                ? totalScaledWeight * totalScaledWeight / squaredScaledWeight
                : 0.0;
        double maxNormalizedWeight = 1.0 / totalScaledWeight;
        PackStarFunctionalObservableResult.Status status;
        double restrictionFreeEnergy;
        String diagnostic;
        if (invalidSampleCount > 0) {
            status = PackStarFunctionalObservableResult.Status.UNRESOLVED_INVALID_SAMPLES;
            restrictionFreeEnergy = Double.NaN;
            diagnostic = "one or more final sample weights were non-finite";
        } else if (!(hitScaledWeight > 0.0)) {
            status = PackStarFunctionalObservableResult.Status.UNRESOLVED_NO_HITS;
            restrictionFreeEnergy = Double.NaN;
            diagnostic = "the final weighted sample contained no event hit";
        } else {
            status = PackStarFunctionalObservableResult.Status.RESOLVED;
            restrictionFreeEnergy = -rt * Math.log(probability);
            diagnostic = "resolved from the frozen final proposal sample";
        }

        return PackStarFunctionalObservableResult.fromSamples(
                eventName, status, sampleCount, weightedSampleCount,
                invalidSampleCount, hitCount, probability,
                restrictionFreeEnergy, effectiveSampleSize,
                maxNormalizedWeight, rawHitFraction, rt, diagnostic);
    }
}
