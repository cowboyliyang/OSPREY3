/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
*/

package edu.duke.cs.osprey.packstar;

import java.util.Locale;

/**
 * Weighted final-sample estimate of a PACK* functional event.
 *
 * <p>The probability is a self-normalized importance estimate.  The
 * restriction free energy is reported only when the event has a positive
 * weighted hit mass and every sample weight was valid.  A zero-hit estimate is
 * kept as probability zero but its free energy is NaN and its status is
 * unresolved, because a finite sample cannot establish that the event is
 * physically impossible.</p>
 */
public final class PackStarFunctionalObservableResult {

    public enum Status {
        NOT_CONFIGURED,
        NOT_COMPUTED,
        RESOLVED,
        UNRESOLVED_NO_FINITE_WEIGHT,
        UNRESOLVED_NO_HITS,
        UNRESOLVED_INVALID_SAMPLES,
        UNRESOLVED_EVENT_ERROR
    }

    private final String eventName;
    private final Status status;
    private final int sampleCount;
    private final int weightedSampleCount;
    private final int invalidSampleCount;
    private final int hitCount;
    private final double probability;
    private final double restrictionFreeEnergy;
    private final double effectiveSampleSize;
    private final double maxNormalizedWeight;
    private final double rawHitFraction;
    private final double rt;
    private final String diagnostic;

    private PackStarFunctionalObservableResult(
            String eventName, Status status,
            int sampleCount, int weightedSampleCount,
            int invalidSampleCount, int hitCount,
            double probability, double restrictionFreeEnergy,
            double effectiveSampleSize, double maxNormalizedWeight,
            double rawHitFraction, double rt, String diagnostic) {
        this.eventName = eventName;
        this.status = status;
        this.sampleCount = sampleCount;
        this.weightedSampleCount = weightedSampleCount;
        this.invalidSampleCount = invalidSampleCount;
        this.hitCount = hitCount;
        this.probability = probability;
        this.restrictionFreeEnergy = restrictionFreeEnergy;
        this.effectiveSampleSize = effectiveSampleSize;
        this.maxNormalizedWeight = maxNormalizedWeight;
        this.rawHitFraction = rawHitFraction;
        this.rt = rt;
        this.diagnostic = diagnostic == null ? "" : diagnostic;
    }

    static PackStarFunctionalObservableResult notConfigured() {
        return new PackStarFunctionalObservableResult(
                null, Status.NOT_CONFIGURED,
                0, 0, 0, 0,
                Double.NaN, Double.NaN,
                Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, "no functional event configured");
    }

    static PackStarFunctionalObservableResult notComputed(String eventName) {
        return notComputed(eventName, "final sample has not been evaluated");
    }

    static PackStarFunctionalObservableResult notComputed(String eventName, String reason) {
        return new PackStarFunctionalObservableResult(
                eventName, Status.NOT_COMPUTED,
                0, 0, 0, 0,
                Double.NaN, Double.NaN,
                Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, reason);
    }

    static PackStarFunctionalObservableResult fromSamples(
            String eventName, Status status,
            int sampleCount, int weightedSampleCount,
            int invalidSampleCount, int hitCount,
            double probability, double restrictionFreeEnergy,
            double effectiveSampleSize, double maxNormalizedWeight,
            double rawHitFraction, double rt, String diagnostic) {
        return new PackStarFunctionalObservableResult(
                eventName, status, sampleCount, weightedSampleCount,
                invalidSampleCount, hitCount, probability,
                restrictionFreeEnergy, effectiveSampleSize,
                maxNormalizedWeight, rawHitFraction, rt, diagnostic);
    }

    static PackStarFunctionalObservableResult eventError(
            String eventName, int sampleCount, int hitCount, double rt,
            String diagnostic) {
        return new PackStarFunctionalObservableResult(
                eventName, Status.UNRESOLVED_EVENT_ERROR,
                sampleCount, 0, 0, hitCount,
                Double.NaN, Double.NaN,
                Double.NaN, Double.NaN,
                Double.NaN,
                rt, diagnostic);
    }

    public String getEventName() {
        return eventName;
    }

    public Status getStatus() {
        return status;
    }

    public boolean isResolved() {
        return status == Status.RESOLVED;
    }

    public int getSampleCount() {
        return sampleCount;
    }

    /** Number of samples with a finite, nonzero log weight. */
    public int getWeightedSampleCount() {
        return weightedSampleCount;
    }

    public int getInvalidSampleCount() {
        return invalidSampleCount;
    }

    public int getHitCount() {
        return hitCount;
    }

    /** Estimated functional-conformation probability P_C. */
    public double getProbability() {
        return probability;
    }

    /** Estimated restriction free energy in the PACK* energy units. */
    public double getRestrictionFreeEnergy() {
        return restrictionFreeEnergy;
    }

    /** Effective sample size from the final importance weights. */
    public double getEffectiveSampleSize() {
        return effectiveSampleSize;
    }

    /** Largest normalized final importance weight. */
    public double getMaxNormalizedWeight() {
        return maxNormalizedWeight;
    }

    /** Unweighted hit fraction, useful for diagnosing rare-event coverage. */
    public double getRawHitFraction() {
        return rawHitFraction;
    }

    /** The thermal energy RT used for the restriction free energy. */
    public double getRT() {
        return rt;
    }

    public String getDiagnostic() {
        return diagnostic;
    }

    @Override
    public String toString() {
        return String.format(Locale.ROOT,
                "PackStarFunctionalObservableResult[event=%s,status=%s,P_C=%.9g,G_restrict=%.9g,"
                        + "samples=%d,weightedSamples=%d,hits=%d,ESS=%.6g,maxWeight=%.6g]",
                eventName, status, probability, restrictionFreeEnergy,
                sampleCount, weightedSampleCount, hitCount,
                effectiveSampleSize, maxNormalizedWeight);
    }
}
