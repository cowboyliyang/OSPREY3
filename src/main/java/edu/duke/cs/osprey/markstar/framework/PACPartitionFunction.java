/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
**
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.branch.*;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.IntConsumer;

/**
 * PAC (Probably Approximately Correct) Partition Function Estimation
 * via Rao-Blackwellized Importance Sampling.
 *
 * Algorithm:
 *   Z = Z_min * E_p[phi(c)]
 *   where phi(c) = exp(-g(c)/kT), g(c) = E_true(c) - E_min(c)
 *
 *   Rao-Blackwellization: decompose g(c) = f_pair(c) + residual(c)
 *   where f_pair is pair-decomposable (eta corrections learned from CCD samples).
 *
 *   Z = Z_corrected * E_p'[psi(c)]
 *   where psi(c) = exp(-residual(c)/kT), with much lower variance than phi.
 *
 *   Bernstein inequality gives PAC bound on E[psi].
 */
public class PACPartitionFunction {

    // Configuration
    private static final String PAC_SAMPLES_PROPERTY = "branchmarkstar.pac.samples";
    private static final String PAC_CONFIDENCE_PROPERTY = "branchmarkstar.pac.confidence";
    private static final String PAC_ADAPTIVE_PROPERTY = "branchmarkstar.pac.adaptive";
    private static final String PAC_MIN_SAMPLES_PROPERTY = "branchmarkstar.pac.minSamples";
    private static final String PAC_MAX_SAMPLES_PROPERTY = "branchmarkstar.pac.maxSamples";
    private static final String PAC_BATCH_SIZE_PROPERTY = "branchmarkstar.pac.batchSize";
    private static final String PAC_TARGET_EPSILON_PROPERTY = "branchmarkstar.pac.targetEpsilon";
    private static final String PAC_SAMPLING_BATCHED_PROPERTY = "branchmarkstar.pac.sampling.batched";
    private static final String PAC_SAMPLING_PARALLEL_PROPERTY = "branchmarkstar.pac.sampling.parallel";
    private static final String PAC_SAMPLING_THREADS_PROPERTY = "branchmarkstar.pac.sampling.threads";
    private static final String PAC_SAMPLING_LARGE_LAMBDA_PROPERTY = "branchmarkstar.pac.sampling.largeLambdaThreshold";
    private static final String PAC_SAMPLING_PROGRESS_PROPERTY = "branchmarkstar.pac.sampling.progress";
    private static final String DP_PARALLEL_THREADS_PROPERTY = "branchmarkstar.dp.parallel.threads";
    private static final String DP_CACHE_ENABLED_PROPERTY = "branchmarkstar.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "branchmarkstar.dp.cache.maxEntries";
    private static final String DP_CACHE_MAX_TABLE_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTableBytes";
    private static final String DP_CACHE_MAX_TOTAL_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTotalBytes";
    private static final String DP_CACHE_SKIP_IF_M_STATES_PROPERTY = "branchmarkstar.dp.cache.skipIfMStates";
    private static final int DEFAULT_SAMPLES = 500;
    private static final double DEFAULT_CONFIDENCE = 0.05; // delta = 0.05 => 95% confidence
    private static final int DEFAULT_ADAPTIVE_MIN_SAMPLES = 100;
    private static final int DEFAULT_ADAPTIVE_BATCH_SIZE = 200;
    private static final double DEFAULT_TARGET_EPSILON = 0.683;
    private static final int DEFAULT_SAMPLING_LARGE_LAMBDA = 65_536;
    private static final int DEFAULT_DP_CACHE_MAX_ENTRIES = 20000;
    private static final long DEFAULT_DP_CACHE_MAX_TABLE_BYTES = 256L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_MAX_TOTAL_BYTES = 4L * 1024L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_SKIP_IF_M_STATES = 8_000_000L;

    private static final Object CORRECTED_DP_CACHE_LOCK = new Object();
    private static final LinkedHashMap<String, CachedDPTable> CORRECTED_DP_CACHE =
            new LinkedHashMap<>(1024, 0.75f, true);
    private static long correctedDPCacheBytes = 0L;

    // Inputs
    private final RootedTreeNode rootedRoot;
    private final RootedTreeEdge rootedRootEdge;
    private final EnergyMatrix branchMinimizingEmat;
    private final EnergyMatrix branchRigidEmat;
    private final InteractionGraph interactionGraph;
    private final ConfEnergyCalculator minimizingEcalc;
    private final RCs rcs;
    private final SimpleConfSpace confSpace;
    private final double RT;

    // Configuration
    private final int numSamples;
    private final double delta; // confidence parameter
    private final boolean adaptiveSampling;
    private final int minSamples;
    private final int maxSamples;
    private final int batchSize;
    private final double targetEpsilon;
    private final boolean batchedSampling;
    private final int samplingThreads;
    private final int samplingLargeLambdaThreshold;
    private final boolean samplingProgress;
    private final boolean dpCacheEnabled;
    private final int dpCacheMaxEntries;
    private final long dpCacheMaxTableBytes;
    private final long dpCacheMaxTotalBytes;
    private final long dpCacheSkipIfMStates;

    // Results
    private BigDecimal zLower;
    private BigDecimal zUpper;
    private double epsilon;
    private double logZLowerPAC;
    private double logZUpperPAC;
    private int totalCCDCalls;

    // Original logZ_min before Phase 4 modifies DP tables
    private double logZMinOriginal;

    // Statistics
    private double meanPsi;
    private double varPsi;
    private double cvPsi; // coefficient of variation
    private double meanResidual;
    private double stdResidual;

    private final BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private final LocalRCMap[] localRCByGlobalRC;

    public PACPartitionFunction(RootedTreeNode rootedRoot,
                                 RootedTreeEdge rootedRootEdge,
                                 EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace) {
        this(rootedRoot, rootedRootEdge,
                branchMinimizingEmat, branchRigidEmat,
                interactionGraph, minimizingEcalc,
                rcs, confSpace,
                Double.NaN);
    }

    public PACPartitionFunction(RootedTreeNode rootedRoot,
                                 RootedTreeEdge rootedRootEdge,
                                 EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace,
                                 double requestedTargetEpsilon) {
        this.rootedRoot = rootedRoot;
        this.rootedRootEdge = rootedRootEdge;
        this.branchMinimizingEmat = branchMinimizingEmat;
        this.branchRigidEmat = branchRigidEmat;
        this.interactionGraph = interactionGraph;
        this.minimizingEcalc = minimizingEcalc;
        this.rcs = rcs;
        this.confSpace = confSpace;
        this.RT = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;

        this.numSamples = Math.max(1, getConfigInteger(PAC_SAMPLES_PROPERTY, DEFAULT_SAMPLES));
        String deltaStr = getConfigProperty(PAC_CONFIDENCE_PROPERTY, null);
        this.delta = (deltaStr != null) ? Double.parseDouble(deltaStr) : DEFAULT_CONFIDENCE;
        this.adaptiveSampling = getConfigBoolean(PAC_ADAPTIVE_PROPERTY, false);
        this.maxSamples = Math.max(1,
                getConfigInteger(PAC_MAX_SAMPLES_PROPERTY, this.numSamples));
        this.minSamples = Math.max(1, Math.min(this.maxSamples,
                getConfigInteger(PAC_MIN_SAMPLES_PROPERTY,
                        Math.min(DEFAULT_ADAPTIVE_MIN_SAMPLES, this.maxSamples))));
        this.batchSize = Math.max(1,
                getConfigInteger(PAC_BATCH_SIZE_PROPERTY,
                        Math.min(DEFAULT_ADAPTIVE_BATCH_SIZE, this.maxSamples)));
        double defaultTarget = Double.isFinite(requestedTargetEpsilon) && requestedTargetEpsilon > 0.0
                ? requestedTargetEpsilon
                : DEFAULT_TARGET_EPSILON;
        this.targetEpsilon = Math.max(0.0,
                getConfigDouble(PAC_TARGET_EPSILON_PROPERTY, defaultTarget));
        this.batchedSampling = getConfigBoolean(PAC_SAMPLING_BATCHED_PROPERTY, true);
        this.samplingThreads = resolveConfiguredSamplingThreads();
        this.samplingLargeLambdaThreshold = Math.max(1,
                getConfigInteger(PAC_SAMPLING_LARGE_LAMBDA_PROPERTY,
                        DEFAULT_SAMPLING_LARGE_LAMBDA));
        this.samplingProgress = getConfigBoolean(PAC_SAMPLING_PROGRESS_PROPERTY, true);
        this.dpCacheEnabled = getConfigBoolean(DP_CACHE_ENABLED_PROPERTY, true);
        this.dpCacheMaxEntries = Math.max(0,
                getConfigInteger(DP_CACHE_MAX_ENTRIES_PROPERTY, DEFAULT_DP_CACHE_MAX_ENTRIES));
        this.dpCacheMaxTableBytes = Math.max(0L,
                getConfigBytes(DP_CACHE_MAX_TABLE_BYTES_PROPERTY, DEFAULT_DP_CACHE_MAX_TABLE_BYTES));
        this.dpCacheMaxTotalBytes = Math.max(0L,
                getConfigBytes(DP_CACHE_MAX_TOTAL_BYTES_PROPERTY, DEFAULT_DP_CACHE_MAX_TOTAL_BYTES));
        this.dpCacheSkipIfMStates = Math.max(0L,
                getConfigLong(DP_CACHE_SKIP_IF_M_STATES_PROPERTY, DEFAULT_DP_CACHE_SKIP_IF_M_STATES));
        this.localRCByGlobalRC = buildLocalRCMaps(rcs);
    }

    private static class CachedDPTable {
        final double[] lower;
        final double[] upper;
        final long bytes;

        CachedDPTable(double[] lower, double[] upper, long bytes) {
            this.lower = lower;
            this.upper = upper;
            this.bytes = bytes;
        }
    }

    private static class DPCacheStats {
        int hits = 0;
        int misses = 0;
        int stores = 0;
        int skippedLarge = 0;
        int evictions = 0;
    }

    private static class CorrectedDPResult {
        final double logZCorrected;
        final DPCacheStats cacheStats;

        CorrectedDPResult(double logZCorrected, DPCacheStats cacheStats) {
            this.logZCorrected = logZCorrected;
            this.cacheStats = cacheStats;
        }
    }

    private static class LocalRCMap {
        final int[] dense;
        final Map<Integer, Integer> sparse;

        LocalRCMap(int[] dense, Map<Integer, Integer> sparse) {
            this.dense = dense;
            this.sparse = sparse;
        }

        int get(int globalRC) {
            if (dense != null) {
                if (globalRC >= 0 && globalRC < dense.length) {
                    return dense[globalRC];
                }
                return -1;
            }
            Integer local = sparse.get(globalRC);
            return local == null ? -1 : local;
        }
    }

    private static class SampleGroup {
        final long mIdx;
        int[] sampleIndices = new int[4];
        int count = 0;

        SampleGroup(long mIdx) {
            this.mIdx = mIdx;
        }

        void add(int sampleIndex) {
            if (count >= sampleIndices.length) {
                sampleIndices = Arrays.copyOf(sampleIndices, sampleIndices.length * 2);
            }
            sampleIndices[count++] = sampleIndex;
        }
    }

    private static String getConfigProperty(String key, String defaultValue) {
        String value = System.getProperty(key);
        if (value == null) {
            value = System.getProperty("osprey." + key);
        }
        return value != null ? value : defaultValue;
    }

    private static int getConfigInteger(String key, int defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("[PAC] Invalid integer for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static long getConfigLong(String key, long defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Long.parseLong(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("[PAC] Invalid long for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static double getConfigDouble(String key, double defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("[PAC] Invalid double for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static long getConfigBytes(String key, long defaultValue) {
        String value = getConfigProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return parseByteCount(value.trim());
        } catch (NumberFormatException e) {
            System.err.println("[PAC] Invalid byte count for '" + key
                    + "': '" + value + "', using " + defaultValue + ".");
            return defaultValue;
        }
    }

    private static long parseByteCount(String value) {
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace("_", "");
        long multiplier = 1L;
        if (normalized.endsWith("kib")) {
            multiplier = 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("kb") || normalized.endsWith("k")) {
            multiplier = 1024L;
            normalized = normalized.replaceAll("kb?$", "");
        } else if (normalized.endsWith("mib")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("mb") || normalized.endsWith("m")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.replaceAll("mb?$", "");
        } else if (normalized.endsWith("gib")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("gb") || normalized.endsWith("g")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("gb?$", "");
        } else if (normalized.endsWith("tib")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("tb") || normalized.endsWith("t")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("tb?$", "");
        } else if (normalized.endsWith("b")) {
            normalized = normalized.substring(0, normalized.length() - 1);
        }

        double amount = Double.parseDouble(normalized.trim());
        if (amount < 0 || amount > Long.MAX_VALUE / (double) multiplier) {
            throw new NumberFormatException(value);
        }
        return (long) (amount * multiplier);
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        String value = getConfigProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    private static int resolveConfiguredSamplingThreads() {
        if (!getConfigBoolean(PAC_SAMPLING_PARALLEL_PROPERTY, true)) {
            return 1;
        }

        int configured = getConfigInteger(PAC_SAMPLING_THREADS_PROPERTY, 0);
        if (configured <= 0) {
            configured = getConfigInteger(DP_PARALLEL_THREADS_PROPERTY, 0);
        }
        if (configured <= 0) {
            configured = Runtime.getRuntime().availableProcessors();
        }
        return Math.max(1, configured);
    }

    private static LocalRCMap[] buildLocalRCMaps(RCs rcs) {
        LocalRCMap[] maps = new LocalRCMap[rcs.getNumPos()];
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            int numRCs = rcs.getNum(pos);
            int maxRC = -1;
            for (int local = 0; local < numRCs; local++) {
                maxRC = Math.max(maxRC, rcs.get(pos, local));
            }

            if (maxRC >= 0 && maxRC <= Math.max(4096, numRCs * 8)) {
                int[] dense = new int[maxRC + 1];
                Arrays.fill(dense, -1);
                for (int local = 0; local < numRCs; local++) {
                    dense[rcs.get(pos, local)] = local;
                }
                maps[pos] = new LocalRCMap(dense, null);
            } else {
                Map<Integer, Integer> sparse = new HashMap<>(Math.max(16, numRCs * 2));
                for (int local = 0; local < numRCs; local++) {
                    sparse.put(rcs.get(pos, local), local);
                }
                maps[pos] = new LocalRCMap(null, sparse);
            }
        }
        return maps;
    }

    private static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }

    private static ThreadFactory daemonThreadFactory(String namePrefix) {
        AtomicInteger counter = new AtomicInteger(0);
        return runnable -> {
            Thread thread = new Thread(runnable, namePrefix + "-" + counter.getAndIncrement());
            thread.setDaemon(true);
            return thread;
        };
    }

    /**
     * Run the full PAC estimation pipeline.
     * Returns the epsilon achieved.
     */
    public double compute() {
        long startTime = System.currentTimeMillis();

        // Phase 0: Z_min from existing DP (already computed)
        double logZMin = rootedRootEdge.getLogZUpper(0);
        double logZRigid = rootedRootEdge.getLogZLower(0);
        logZMinOriginal = logZMin; // save before Phase 4 modifies DP tables
        System.out.println("[PAC] Phase 0: logZ_min=" + String.format("%.4f", logZMin)
                + ", logZ_rigid=" + String.format("%.4f", logZRigid)
                + ", gap=" + String.format("%.4f", logZMin - logZRigid));

        if (logZMin == Double.NEGATIVE_INFINITY) {
            setZeroBounds("initial DP upper bound is zero");
            System.out.println("[PAC] Degenerate zero-mass state; skipping sampling and returning epsilon=1.0");
            return epsilon;
        }
        if (!Double.isFinite(logZMin)) {
            setZeroBounds("initial DP upper bound is non-finite: " + logZMin);
            System.out.println("[PAC] Degenerate non-finite state; skipping sampling and returning epsilon=1.0");
            return epsilon;
        }

        if (adaptiveSampling) {
            runAdaptiveSamplingPAC(startTime, logZRigid);
        } else {
            runFixedSamplingPAC(startTime, logZRigid);
        }

        return epsilon;
    }

    private void runFixedSamplingPAC(long startTime, double logZRigid) {
        Random rng = new Random(42); // reproducible

        // Phase 1: Sample conformations from p(c) proportional to exp(-E_min/kT)
        long phase1Start = System.currentTimeMillis();
        List<int[]> sampledConfs = sampleConformationsFromDP(numSamples, rng);
        long phase1Time = System.currentTimeMillis() - phase1Start;
        System.out.println("[PAC] Phase 1: Sampled " + sampledConfs.size()
                + " conformations in " + phase1Time + " ms");

        // Phase 2: Parallel CCD minimization
        long phase2Start = System.currentTimeMillis();
        List<CCDResult> ccdResults = runParallelCCD(sampledConfs);
        long phase2Time = System.currentTimeMillis() - phase2Start;
        totalCCDCalls = ccdResults.size();
        System.out.println("[PAC] Phase 2: " + totalCCDCalls + " CCD minimizations in "
                + phase2Time + " ms");

        runCorrectionAndBoundPhases(ccdResults, logZRigid, "Phase");
        printFinalSummary(startTime);
    }

    private void runAdaptiveSamplingPAC(long startTime, double logZRigid) {
        int adaptiveMax = Math.max(maxSamples, minSamples);
        int adaptiveBatch = Math.min(batchSize, adaptiveMax);
        System.out.println("[PAC] Adaptive sampling enabled: minSamples=" + minSamples
                + ", maxSamples=" + adaptiveMax
                + ", batchSize=" + adaptiveBatch
                + ", targetEpsilon=" + String.format("%.6f", targetEpsilon));

        Random rng = new Random(42);
        List<CCDResult> allCCDResults = new ArrayList<>(adaptiveMax);
        int requestedSamples = 0;
        int round = 0;
        long totalSamplingTime = 0L;
        long totalCCDTime = 0L;
        boolean proposalDPIsOriginal = true;

        while (requestedSamples < adaptiveMax) {
            round++;
            int targetValidSamples = allCCDResults.isEmpty()
                    ? Math.min(minSamples, adaptiveMax)
                    : Math.min(allCCDResults.size() + adaptiveBatch, adaptiveMax);
            int batchRequest = Math.max(1, targetValidSamples - allCCDResults.size());
            batchRequest = Math.min(batchRequest, adaptiveMax - requestedSamples);

            if (!proposalDPIsOriginal) {
                restoreOriginalProposalDPForAdaptiveSampling();
                proposalDPIsOriginal = true;
            }

            long phase1Start = System.currentTimeMillis();
            List<int[]> batch = sampleConformationsFromDP(batchRequest, rng);
            long phase1Time = System.currentTimeMillis() - phase1Start;
            totalSamplingTime += phase1Time;
            requestedSamples += batchRequest;
            System.out.println("[PAC] Adaptive round " + round + ": sampled "
                    + batch.size() + "/" + batchRequest
                    + " candidate conformations in " + phase1Time + " ms"
                    + ", requested=" + requestedSamples + "/" + adaptiveMax);

            if (batch.isEmpty()) {
                System.out.println("[PAC] Adaptive stop: no valid sampled conformations");
                break;
            }

            long phase2Start = System.currentTimeMillis();
            List<CCDResult> batchResults = runParallelCCD(batch);
            long phase2Time = System.currentTimeMillis() - phase2Start;
            totalCCDTime += phase2Time;
            allCCDResults.addAll(batchResults);
            totalCCDCalls = allCCDResults.size();
            System.out.println("[PAC] Adaptive round " + round + ": "
                    + batchResults.size() + " new CCD minimizations, "
                    + totalCCDCalls + " total, phase2=" + phase2Time + " ms");

            runCorrectionAndBoundPhases(allCCDResults, logZRigid,
                    "Adaptive round " + round);
            proposalDPIsOriginal = false;

            if (totalCCDCalls >= minSamples && epsilon <= targetEpsilon) {
                System.out.println("[PAC] Adaptive stop: epsilon="
                        + String.format("%.6f", epsilon)
                        + " <= target=" + String.format("%.6f", targetEpsilon)
                        + " with " + totalCCDCalls + " CCD calls");
                break;
            }

            if (batchResults.isEmpty()) {
                System.out.println("[PAC] Adaptive stop: no valid new CCD samples");
                break;
            }
        }

        System.out.println("[PAC] Adaptive sampling summary: requestedSamples="
                + requestedSamples + ", validCCD=" + totalCCDCalls
                + ", samplingMs=" + totalSamplingTime
                + ", ccdMs=" + totalCCDTime);
        printFinalSummary(startTime);
    }

    private void restoreOriginalProposalDPForAdaptiveSampling() {
        long restoreStart = System.currentTimeMillis();
        RootedTreeEdge.postOrderInitIncremental(rootedRoot,
                branchRigidEmat, branchMinimizingEmat, interactionGraph, RT);
        RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
        System.out.println("[PAC] Adaptive restored original proposal DP in "
                + (System.currentTimeMillis() - restoreStart) + " ms");
    }

    private void runCorrectionAndBoundPhases(List<CCDResult> ccdResults,
                                             double logZRigid,
                                             String labelPrefix) {
        long phase3Start = System.currentTimeMillis();
        EtaCorrections eta = extractEtaCorrections(ccdResults);
        long phase3Time = System.currentTimeMillis() - phase3Start;
        System.out.println("[PAC] " + labelPrefix + " 3: Eta extraction in " + phase3Time + " ms"
                + ", oneBodyTerms=" + eta.oneBodyCount
                + ", pairTerms=" + eta.pairCount);

        long phase4Start = System.currentTimeMillis();
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta);
        double logZCorrected = correctedDP.logZCorrected;
        long phase4Time = System.currentTimeMillis() - phase4Start;
        System.out.println("[PAC] " + labelPrefix + " 4: logZ_corrected="
                + String.format("%.4f", logZCorrected)
                + " (improvement over logZ_rigid: "
                + String.format("%.4f", logZCorrected - logZRigid) + ")"
                + " in " + phase4Time + " ms"
                + formatDPCacheStats(correctedDP.cacheStats));

        long phase5Start = System.currentTimeMillis();
        computePACBound(ccdResults, correctedEmat, logZCorrected);
        long phase5Time = System.currentTimeMillis() - phase5Start;
        System.out.println("[PAC] " + labelPrefix + " 5: epsilon="
                + String.format("%.6f", epsilon)
                + ", meanPsi=" + String.format("%.6f", meanPsi)
                + ", cvPsi=" + String.format("%.4f", cvPsi)
                + ", meanResidual=" + String.format("%.4f", meanResidual) + " kcal/mol"
                + ", stdResidual=" + String.format("%.4f", stdResidual) + " kcal/mol"
                + " in " + phase5Time + " ms");
    }

    private void printFinalSummary(long startTime) {
        long totalTime = System.currentTimeMillis() - startTime;
        System.out.println("[PAC] Total: " + totalTime + " ms, " + totalCCDCalls
                + " CCD calls, epsilon=" + String.format("%.6f", epsilon)
                + ", confidence=" + String.format("%.2f%%", (1.0 - delta) * 100));
        System.out.println("[PAC] Z bounds: lower=" + String.format("%.6e", zLower.doubleValue())
                + ", upper=" + String.format("%.6e", zUpper.doubleValue())
                + ", log10Lower=" + formatLog10(logZLowerPAC)
                + ", log10Upper=" + formatLog10(logZUpperPAC));
    }

    // ========== Phase 1: Ancestral sampling from DP ==========

    /**
     * Sample conformations from p(c) ∝ exp(-E_min(c)/kT) using ancestral sampling
     * on the tree decomposition.
     *
     * Top-down traversal: for each edge, sample lambda-state proportional to
     * exp(-fullEnergyMin[mIdx][lIdx]/RT), conditioned on the M-state from parent.
     */
    private List<int[]> sampleConformationsFromDP(int n, Random rng) {
        if (!batchedSampling) {
            return sampleConformationsFromDPLegacy(n, rng);
        }

        int numPos = rcs.getNumPos();
        int[][] confs = new int[n][numPos];
        for (int[] conf : confs) {
            Arrays.fill(conf, -1);
        }

        SplittableRandom[] sampleRngs = new SplittableRandom[n];
        long masterSeed = rng.nextLong();
        for (int s = 0; s < n; s++) {
            sampleRngs[s] = new SplittableRandom(mix64(masterSeed + 0x9E3779B97F4A7C15L * (long) s));
        }

        List<RootedTreeEdge> topDownOrder = buildTopDownOrder(rootedRootEdge);
        ExecutorService samplingPool = samplingThreads > 1
                ? Executors.newFixedThreadPool(samplingThreads,
                        daemonThreadFactory("branchmarkstar-pac-sample"))
                : null;

        try {
            for (RootedTreeEdge edge : topDownOrder) {
                if (!edge.getIsLambdaEdge()) continue;

                long edgeStart = System.currentTimeMillis();
                Map<Long, SampleGroup> groups = groupSamplesByMIdx(edge, confs);
                sampleEdgeGroups(edge, confs, sampleRngs, groups, samplingPool);
                logSamplingEdgeProgress(edge, n, groups, System.currentTimeMillis() - edgeStart);
            }
        } finally {
            if (samplingPool != null) {
                samplingPool.shutdownNow();
            }
        }

        List<int[]> conformations = new ArrayList<>(n);
        for (int[] conf : confs) {
            boolean valid = true;
            for (int p = 0; p < numPos; p++) {
                if (conf[p] < 0) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                conformations.add(conf);
            }
        }

        return conformations;
    }

    private List<int[]> sampleConformationsFromDPLegacy(int n, Random rng) {
        int numPos = rcs.getNumPos();
        List<int[]> conformations = new ArrayList<>(n);
        List<RootedTreeEdge> topDownOrder = buildTopDownOrder(rootedRootEdge);

        for (int s = 0; s < n; s++) {
            int[] conf = new int[numPos];
            Arrays.fill(conf, -1);

            for (RootedTreeEdge edge : topDownOrder) {
                if (!edge.getIsLambdaEdge()) continue;
                int[] mPositions = edge.getMPositionsSorted();
                int[] lambdaPositions = edge.getLambdaPositionsSorted();

                int[] mLocalRCs = new int[mPositions.length];
                for (int i = 0; i < mPositions.length; i++) {
                    int globalRC = conf[mPositions[i]];
                    mLocalRCs[i] = findLocalRCIndex(mPositions[i], globalRC);
                }
                long mIdx = edge.computeIndexInA(mLocalRCs);

                int lIdx = sampleLambdaState(edge, mIdx, rng);

                int[] lambdaLocalRCs = edge.decodeLambdaStatePublic(lIdx);
                for (int i = 0; i < lambdaPositions.length; i++) {
                    conf[lambdaPositions[i]] = rcs.get(lambdaPositions[i], lambdaLocalRCs[i]);
                }
            }

            // Verify all positions assigned
            boolean valid = true;
            for (int p = 0; p < numPos; p++) {
                if (conf[p] < 0) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                conformations.add(conf);
            }
        }

        return conformations;
    }

    private Map<Long, SampleGroup> groupSamplesByMIdx(RootedTreeEdge edge, int[][] confs) {
        Map<Long, SampleGroup> groups = new HashMap<>();
        int[] mPositions = edge.getMPositionsSorted();
        int[] mLocalRCs = new int[mPositions.length];

        for (int s = 0; s < confs.length; s++) {
            for (int i = 0; i < mPositions.length; i++) {
                int pos = mPositions[i];
                int globalRC = confs[s][pos];
                if (globalRC < 0) {
                    throw new IllegalStateException("PAC sampling reached edge before M position "
                            + pos + " was assigned");
                }
                mLocalRCs[i] = findLocalRCIndex(pos, globalRC);
            }
            long mIdx = edge.computeIndexInA(mLocalRCs);
            groups.computeIfAbsent(mIdx, SampleGroup::new).add(s);
        }

        return groups;
    }

    private void sampleEdgeGroups(RootedTreeEdge edge, int[][] confs,
                                  SplittableRandom[] sampleRngs,
                                  Map<Long, SampleGroup> groups,
                                  ExecutorService samplingPool) {
        int totalLambda = edge.getTotalLambdaStates();
        if (totalLambda <= 0) {
            throw new IllegalStateException("Lambda edge has no lambda states");
        }
        if (totalLambda == 1) {
            for (SampleGroup group : groups.values()) {
                for (int i = 0; i < group.count; i++) {
                    writeLambdaState(edge, 0, confs[group.sampleIndices[i]]);
                }
            }
            return;
        }

        if (samplingPool == null || samplingThreads <= 1) {
            for (SampleGroup group : groups.values()) {
                sampleGroupSerial(edge, group, confs, sampleRngs);
            }
            return;
        }

        if (totalLambda >= samplingLargeLambdaThreshold) {
            processSingletonGroupsParallel(edge, groups.values(), confs, sampleRngs, samplingPool);
            for (SampleGroup group : groups.values()) {
                if (group.count >= 2) {
                    sampleGroupWithParallelCDF(edge, group, confs, sampleRngs, samplingPool);
                }
            }
        } else {
            processSmallGroupsParallel(edge, groups.values(), confs, sampleRngs, samplingPool);
        }
    }

    private void processSingletonGroupsParallel(RootedTreeEdge edge,
                                                Collection<SampleGroup> groups,
                                                int[][] confs,
                                                SplittableRandom[] sampleRngs,
                                                ExecutorService samplingPool) {
        List<Future<?>> futures = new ArrayList<>();
        for (SampleGroup group : groups) {
            if (group.count != 1) continue;
            futures.add(samplingPool.submit(() -> {
                int sampleIndex = group.sampleIndices[0];
                int lIdx = sampleLambdaStateStreaming(edge, group.mIdx, sampleRngs[sampleIndex]);
                writeLambdaState(edge, lIdx, confs[sampleIndex]);
            }));
        }
        waitForSamplingFutures(futures);
    }

    private void processSmallGroupsParallel(RootedTreeEdge edge,
                                            Collection<SampleGroup> groups,
                                            int[][] confs,
                                            SplittableRandom[] sampleRngs,
                                            ExecutorService samplingPool) {
        List<Future<?>> futures = new ArrayList<>(groups.size());
        for (SampleGroup group : groups) {
            futures.add(samplingPool.submit(() ->
                    sampleGroupSerial(edge, group, confs, sampleRngs)));
        }
        waitForSamplingFutures(futures);
    }

    private void sampleGroupSerial(RootedTreeEdge edge, SampleGroup group,
                                   int[][] confs, SplittableRandom[] sampleRngs) {
        if (group.count == 1) {
            int sampleIndex = group.sampleIndices[0];
            int lIdx = sampleLambdaStateStreaming(edge, group.mIdx, sampleRngs[sampleIndex]);
            writeLambdaState(edge, lIdx, confs[sampleIndex]);
            return;
        }

        double[] cdf = buildConditionalCDFSerial(edge, group.mIdx);
        sampleGroupFromCDF(edge, group, confs, sampleRngs, cdf);
    }

    private void sampleGroupWithParallelCDF(RootedTreeEdge edge, SampleGroup group,
                                            int[][] confs, SplittableRandom[] sampleRngs,
                                            ExecutorService samplingPool) {
        double[] cdf = buildConditionalCDFParallel(edge, group.mIdx, samplingPool);
        sampleGroupFromCDF(edge, group, confs, sampleRngs, cdf);
    }

    private void sampleGroupFromCDF(RootedTreeEdge edge, SampleGroup group,
                                    int[][] confs, SplittableRandom[] sampleRngs,
                                    double[] cdf) {
        double total = cdf[cdf.length - 1];
        if (!(total > 0.0) || !Double.isFinite(total)) {
            for (int i = 0; i < group.count; i++) {
                int sampleIndex = group.sampleIndices[i];
                int lIdx = sampleLambdaStateStreaming(edge, group.mIdx, sampleRngs[sampleIndex]);
                writeLambdaState(edge, lIdx, confs[sampleIndex]);
            }
            return;
        }

        for (int i = 0; i < group.count; i++) {
            int sampleIndex = group.sampleIndices[i];
            int lIdx = sampleFromCDF(cdf, total, sampleRngs[sampleIndex]);
            writeLambdaState(edge, lIdx, confs[sampleIndex]);
        }
    }

    private double[] buildConditionalCDFSerial(RootedTreeEdge edge, long mIdx) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] cdf = new double[totalLambda];
        int[] mRCs = edge.hasFsetChildren() ? edge.decodeMStatePublic(mIdx) : null;

        double maxLog = Double.NEGATIVE_INFINITY;
        for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
            double logWeight = computeLambdaLogWeight(edge, mIdx, mRCs, lIdx);
            cdf[lIdx] = logWeight;
            if (logWeight > maxLog) {
                maxLog = logWeight;
            }
        }

        if (!Double.isFinite(maxLog)) {
            return cdf;
        }

        double running = 0.0;
        for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
            running += Math.exp(cdf[lIdx] - maxLog);
            cdf[lIdx] = running;
        }
        return cdf;
    }

    private double[] buildConditionalCDFParallel(RootedTreeEdge edge, long mIdx,
                                                 ExecutorService samplingPool) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] cdf = new double[totalLambda];
        int chunks = Math.min(samplingThreads, totalLambda);
        int[] mRCs = edge.hasFsetChildren() ? edge.decodeMStatePublic(mIdx) : null;
        double[] chunkMax = new double[chunks];
        Arrays.fill(chunkMax, Double.NEGATIVE_INFINITY);

        runParallelChunks(samplingPool, chunks, chunk -> {
            int start = chunkStart(totalLambda, chunks, chunk);
            int end = chunkStart(totalLambda, chunks, chunk + 1);
            double max = Double.NEGATIVE_INFINITY;
            for (int lIdx = start; lIdx < end; lIdx++) {
                double logWeight = computeLambdaLogWeight(edge, mIdx, mRCs, lIdx);
                cdf[lIdx] = logWeight;
                if (logWeight > max) {
                    max = logWeight;
                }
            }
            chunkMax[chunk] = max;
        });

        double maxLog = Double.NEGATIVE_INFINITY;
        for (double max : chunkMax) {
            if (max > maxLog) {
                maxLog = max;
            }
        }
        if (!Double.isFinite(maxLog)) {
            return cdf;
        }

        final double globalMaxLog = maxLog;
        double[] chunkSums = new double[chunks];
        runParallelChunks(samplingPool, chunks, chunk -> {
            int start = chunkStart(totalLambda, chunks, chunk);
            int end = chunkStart(totalLambda, chunks, chunk + 1);
            double sum = 0.0;
            for (int lIdx = start; lIdx < end; lIdx++) {
                double weight = Math.exp(cdf[lIdx] - globalMaxLog);
                cdf[lIdx] = weight;
                sum += weight;
            }
            chunkSums[chunk] = sum;
        });

        double[] chunkOffsets = new double[chunks];
        double offset = 0.0;
        for (int chunk = 0; chunk < chunks; chunk++) {
            chunkOffsets[chunk] = offset;
            offset += chunkSums[chunk];
        }

        runParallelChunks(samplingPool, chunks, chunk -> {
            int start = chunkStart(totalLambda, chunks, chunk);
            int end = chunkStart(totalLambda, chunks, chunk + 1);
            double running = chunkOffsets[chunk];
            for (int lIdx = start; lIdx < end; lIdx++) {
                running += cdf[lIdx];
                cdf[lIdx] = running;
            }
        });

        return cdf;
    }

    private int sampleLambdaStateStreaming(RootedTreeEdge edge, long mIdx,
                                           SplittableRandom rng) {
        int totalLambda = edge.getTotalLambdaStates();
        if (totalLambda <= 1) {
            return 0;
        }

        int[] mRCs = edge.hasFsetChildren() ? edge.decodeMStatePublic(mIdx) : null;
        int bestIdx = totalLambda - 1;
        double bestScore = Double.NEGATIVE_INFINITY;

        for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
            double logWeight = computeLambdaLogWeight(edge, mIdx, mRCs, lIdx);
            if (Double.isNaN(logWeight)) {
                continue;
            }
            double score = logWeight + sampleGumbel(rng);
            if (score > bestScore) {
                bestScore = score;
                bestIdx = lIdx;
            }
        }
        return bestIdx;
    }

    private double computeLambdaLogWeight(RootedTreeEdge edge, long mIdx,
                                          int[] mRCs, int lIdx) {
        if (!edge.hasFsetChildren()) {
            return -edge.getFullEnergyMin(mIdx, lIdx) / RT;
        }

        int[] lambdaRCs = edge.decodeLambdaStatePublic(lIdx);
        double localEnergy = computeLocalEnergy(edge, mRCs, lambdaRCs, branchMinimizingEmat);
        double fSumUpper = 0.0;
        if (edge.getFset() != null) {
            for (RootedTreeEdge fEdge : edge.getFset()) {
                int[] fM = edge.getMstateForFullState(mRCs, lambdaRCs, fEdge);
                long fIdx = fEdge.computeIndexInA(fM);
                fSumUpper += fEdge.getLogZUpper(fIdx);
            }
        }
        return -localEnergy / RT + fSumUpper;
    }

    private void writeLambdaState(RootedTreeEdge edge, int lIdx, int[] conf) {
        int[] lambdaPositions = edge.getLambdaPositionsSorted();
        int[] lambdaLocalRCs = edge.decodeLambdaStatePublic(lIdx);
        for (int i = 0; i < lambdaPositions.length; i++) {
            conf[lambdaPositions[i]] = rcs.get(lambdaPositions[i], lambdaLocalRCs[i]);
        }
    }

    private int sampleFromCDF(double[] cdf, double total, SplittableRandom rng) {
        double target = rng.nextDouble() * total;
        int lo = 0;
        int hi = cdf.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (target < cdf[mid]) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        return lo;
    }

    private double sampleGumbel(SplittableRandom rng) {
        double u = rng.nextDouble();
        if (u <= 0.0) {
            u = Double.MIN_VALUE;
        }
        return -Math.log(-Math.log(u));
    }

    private int chunkStart(int total, int chunks, int chunk) {
        return (int) (((long) total * (long) chunk) / (long) chunks);
    }

    private void runParallelChunks(ExecutorService pool, int chunks, IntConsumer task) {
        List<Future<?>> futures = new ArrayList<>(chunks);
        for (int chunk = 0; chunk < chunks; chunk++) {
            final int c = chunk;
            futures.add(pool.submit(() -> task.accept(c)));
        }
        waitForSamplingFutures(futures);
    }

    private void waitForSamplingFutures(List<Future<?>> futures) {
        try {
            for (Future<?> future : futures) {
                future.get();
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("PAC sampling interrupted", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new RuntimeException("PAC sampling failed", cause);
        }
    }

    private void logSamplingEdgeProgress(RootedTreeEdge edge, int requestedSamples,
                                         Map<Long, SampleGroup> groups,
                                         long elapsedMs) {
        if (!samplingProgress) {
            return;
        }

        int totalLambda = edge.getTotalLambdaStates();
        if (totalLambda < samplingLargeLambdaThreshold && elapsedMs < 1000L) {
            return;
        }

        int singletonGroups = 0;
        int maxGroup = 0;
        for (SampleGroup group : groups.values()) {
            if (group.count == 1) {
                singletonGroups++;
            }
            maxGroup = Math.max(maxGroup, group.count);
        }

        System.out.println("[PAC] Phase 1 edge: lambdaStates=" + totalLambda
                + ", samples=" + requestedSamples
                + ", distinctMIdx=" + groups.size()
                + ", singletonGroups=" + singletonGroups
                + ", maxGroup=" + maxGroup
                + ", threads=" + samplingThreads
                + ", elapsed=" + elapsedMs + " ms");
    }

    /**
     * Find the local RC index for a position given its global RC.
     */
    private int findLocalRCIndex(int pos, int globalRC) {
        int local = localRCByGlobalRC[pos].get(globalRC);
        if (local >= 0) {
            return local;
        }
        throw new IllegalStateException("Global RC " + globalRC + " not found at position " + pos);
    }

    /**
     * Sample a lambda-state proportional to its Boltzmann weight exp(-fullEnergyMin/RT).
     * For leaf edges, uses pre-computed fullEnergyMin.
     * For non-leaf edges, uses local energy + child logZ.
     */
    private int sampleLambdaState(RootedTreeEdge edge, long mIdx, Random rng) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] logWeights = new double[totalLambda];

        if (!edge.hasFsetChildren()) {
            // Leaf edge: weight = exp(-fullEnergyMin/RT). The energy table may be
            // streamed/on-demand instead of materialized as [M x lambda].
            for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
                logWeights[lIdx] = -edge.getFullEnergyMin(mIdx, lIdx) / RT;
            }
        } else {
            // Non-leaf edge: weight = exp(-localEnergy/RT) * product of child Z_upper
            int[] mRCs = edge.decodeMStatePublic(mIdx);
            for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
                int[] lambdaRCs = edge.decodeLambdaStatePublic(lIdx);
                // Compute local energy (lambda-only + lambda-M interactions)
                double localEnergy = computeLocalEnergy(edge, mRCs, lambdaRCs, branchMinimizingEmat);
                double fSumUpper = 0.0;
                if (edge.getFset() != null) {
                    for (RootedTreeEdge fEdge : edge.getFset()) {
                        int[] fM = edge.getMstateForFullState(mRCs, lambdaRCs, fEdge);
                        long fIdx = fEdge.computeIndexInA(fM);
                        fSumUpper += fEdge.getLogZUpper(fIdx);
                    }
                }
                logWeights[lIdx] = -localEnergy / RT + fSumUpper;
            }
        }

        return sampleFromLogWeights(logWeights, rng);
    }

    /**
     * Sample an index from log-space weights using Gumbel-max trick.
     */
    private int sampleFromLogWeights(double[] logWeights, Random rng) {
        // Find max for numerical stability
        double maxLog = Double.NEGATIVE_INFINITY;
        for (double lw : logWeights) {
            if (lw > maxLog) maxLog = lw;
        }

        // Compute normalized probabilities
        double[] probs = new double[logWeights.length];
        double sum = 0;
        for (int i = 0; i < logWeights.length; i++) {
            probs[i] = Math.exp(logWeights[i] - maxLog);
            sum += probs[i];
        }

        // Sample
        double r = rng.nextDouble() * sum;
        double cumulative = 0;
        for (int i = 0; i < probs.length; i++) {
            cumulative += probs[i];
            if (r <= cumulative) return i;
        }
        return probs.length - 1; // fallback
    }

    /**
     * Build top-down traversal order for the compact tree.
     */
    private List<RootedTreeEdge> buildTopDownOrder(RootedTreeEdge rootEdge) {
        List<RootedTreeEdge> order = new ArrayList<>();
        Queue<RootedTreeEdge> queue = new LinkedList<>();
        queue.add(rootEdge);

        while (!queue.isEmpty()) {
            RootedTreeEdge edge = queue.poll();
            if (edge.getIsLambdaEdge()) {
                order.add(edge);
            }
            // Add compact children
            RootedTreeNode leftChild = edge.getCompactLeftChild();
            RootedTreeNode rightChild = edge.getCompactRightChild();
            if (leftChild != null && leftChild.getChildOfEdge() != null) {
                queue.add(leftChild.getChildOfEdge());
            }
            if (rightChild != null && rightChild.getChildOfEdge() != null) {
                queue.add(rightChild.getChildOfEdge());
            }
        }

        return order;
    }

    /**
     * Compute local energy for a (M, lambda) pair on an edge.
     * Local = lambda-only + lambda-M interactions.
     */
    private double computeLocalEnergy(RootedTreeEdge edge, int[] mRCs, int[] lambdaRCs,
                                       EnergyMatrix emat) {
        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        double energy = 0.0;

        // Lambda-only: one-body + lambda-lambda pairs
        for (int i = 0; i < lambdaPositions.length; i++) {
            int posI = lambdaPositions[i];
            int rcI = rcs.get(posI, lambdaRCs[i]);
            energy += emat.getOneBody(posI, rcI);
            for (int j = i + 1; j < lambdaPositions.length; j++) {
                int posJ = lambdaPositions[j];
                if (interactionGraph.hasEdge(posI, posJ)) {
                    int rcJ = rcs.get(posJ, lambdaRCs[j]);
                    energy += emat.getPairwise(posI, rcI, posJ, rcJ);
                }
            }
        }

        // Lambda-M interactions
        for (int i = 0; i < lambdaPositions.length; i++) {
            int posI = lambdaPositions[i];
            int rcI = rcs.get(posI, lambdaRCs[i]);
            for (int j = 0; j < mPositions.length; j++) {
                int posJ = mPositions[j];
                if (interactionGraph.hasEdge(posI, posJ)) {
                    int rcJ = rcs.get(posJ, mRCs[j]);
                    energy += emat.getPairwise(posI, rcI, posJ, rcJ);
                }
            }
        }

        return energy;
    }

    // ========== Phase 2: Parallel CCD minimization ==========

    private static class CCDResult {
        final int[] conf;
        final double eTrue;
        final double eMin;
        final EnergyCalculator.EnergiedParametricMolecule epmol;

        CCDResult(int[] conf, double eTrue, double eMin,
                  EnergyCalculator.EnergiedParametricMolecule epmol) {
            this.conf = conf;
            this.eTrue = eTrue;
            this.eMin = eMin;
            this.epmol = epmol;
        }
    }

    /**
     * Run CCD minimizations in parallel via minimizingEcalc.tasks.
     */
    private List<CCDResult> runParallelCCD(List<int[]> conformations) {
        ConcurrentLinkedQueue<CCDResult> results = new ConcurrentLinkedQueue<>();
        AtomicInteger completed = new AtomicInteger(0);
        int total = conformations.size();

        for (int[] conf : conformations) {
            ResidueInteractions inters = makeSparseFullConfInters(conf);
            RCTuple tuple = new RCTuple(conf);

            minimizingEcalc.calcEnergyAsync(tuple, inters, epmol -> {
                double eTrue = epmol.energy;
                double eMin = computeFullConfPairwiseEnergy(conf, branchMinimizingEmat);
                results.add(new CCDResult(conf, eTrue, eMin, epmol));

                int done = completed.incrementAndGet();
                if (done % 100 == 0 || done == total) {
                    System.out.println("[PAC] CCD progress: " + done + "/" + total);
                }
            });
        }

        minimizingEcalc.tasks.waitForFinish();
        return new ArrayList<>(results);
    }

    // ========== Phase 3: Extract per-term eta corrections ==========

    private static class EtaCorrections {
        // eta_oneBody[pos][rc] = mean correction for one-body term
        final double[][] oneBody;
        final int[][] oneBodyCounts; // per-(pos,rc) counts, NOT per-pos
        // eta_pair[pos1][rc1][pos2][rc2] = mean correction for pair term
        // Stored as map to avoid sparse array issues
        final Map<Long, double[]> pairSums; // key = pack(pos1,rc1,pos2,rc2), [0]=sum, [1]=count
        int oneBodyCount = 0;
        int pairCount = 0;

        EtaCorrections(int numPos, int[] numRCsPerPos) {
            oneBody = new double[numPos][];
            oneBodyCounts = new int[numPos][];
            for (int p = 0; p < numPos; p++) {
                oneBody[p] = new double[numRCsPerPos[p]];
                oneBodyCounts[p] = new int[numRCsPerPos[p]];
            }
            pairSums = new HashMap<>();
        }

        static long packPairKey(int pos1, int rc1, int pos2, int rc2) {
            // Ensure pos1 < pos2
            if (pos1 > pos2) {
                int tmp = pos1; pos1 = pos2; pos2 = tmp;
                tmp = rc1; rc1 = rc2; rc2 = tmp;
            }
            return ((long) pos1 << 48) | ((long) rc1 << 32) | ((long) pos2 << 16) | rc2;
        }

        void addOneBodySample(int pos, int rc, double correction) {
            oneBody[pos][rc] += correction;
            oneBodyCounts[pos][rc]++;
        }

        void addPairSample(int pos1, int rc1, int pos2, int rc2, double correction) {
            long key = packPairKey(pos1, rc1, pos2, rc2);
            double[] val = pairSums.computeIfAbsent(key, k -> new double[2]);
            val[0] += correction;
            val[1] += 1;
        }

        double getOneBodyEta(int pos, int rc) {
            if (oneBodyCounts[pos][rc] == 0) return 0.0;
            return oneBody[pos][rc] / oneBodyCounts[pos][rc];
        }

        double getPairEta(int pos1, int rc1, int pos2, int rc2) {
            long key = packPairKey(pos1, rc1, pos2, rc2);
            double[] val = pairSums.get(key);
            if (val == null || val[1] == 0) return 0.0;
            return val[0] / val[1];
        }
    }

    /**
     * Extract per-term corrections from CCD results.
     *
     * For each minimized conformation, evaluate individual energy terms at the
     * CCD-optimized coordinates and compare with emat_min entries.
     *
     * Uses the forcefield breakdown approach from ConfAnalyzer.
     */
    private EtaCorrections extractEtaCorrections(List<CCDResult> ccdResults) {
        int numPos = rcs.getNumPos();
        int[] numRCs = new int[numPos];
        for (int p = 0; p < numPos; p++) {
            numRCs[p] = branchMinimizingEmat.getNumConfAtPos(p);
        }
        EtaCorrections eta = new EtaCorrections(numPos, numRCs);

        for (CCDResult result : ccdResults) {
            int[] conf = result.conf;

            // Create the forcefield energy function at CCD-optimized coordinates
            ResidueForcefieldEnergy efunc =
                    (ResidueForcefieldEnergy) minimizingEcalc.ecalc.makeEnergyFunction(result.epmol);

            // Evaluate per-term energies at optimized coordinates
            for (int pos = 0; pos < numPos; pos++) {
                if (conf[pos] < 0) continue;
                // Single-position energy at optimized DOFs
                ResidueInteractions singleInters = minimizingEcalc.makeSingleInters(pos, conf[pos]);
                double actualSingle = efunc.makeSubset(singleInters).getEnergy();
                double ematSingle = branchMinimizingEmat.getOneBody(pos, conf[pos]);
                eta.addOneBodySample(pos, conf[pos], actualSingle - ematSingle);
            }

            for (int pos1 = 0; pos1 < numPos; pos1++) {
                if (conf[pos1] < 0) continue;
                for (int pos2 = pos1 + 1; pos2 < numPos; pos2++) {
                    if (conf[pos2] < 0) continue;
                    if (!interactionGraph.hasEdge(pos1, pos2)) continue;

                    ResidueInteractions pairInters =
                            minimizingEcalc.makePairInters(pos1, conf[pos1], pos2, conf[pos2]);
                    double actualPair = efunc.makeSubset(pairInters).getEnergy();
                    double ematPair = branchMinimizingEmat.getPairwise(pos1, conf[pos1], pos2, conf[pos2]);
                    eta.addPairSample(pos1, conf[pos1], pos2, conf[pos2], actualPair - ematPair);
                }
            }

            // Clean up
            efunc = null;
        }

        // Count unique (pos,rc) terms with data
        for (int p = 0; p < numPos; p++) {
            for (int rc = 0; rc < numRCs[p]; rc++) {
                if (eta.oneBodyCounts[p][rc] > 0) eta.oneBodyCount++;
            }
        }
        eta.pairCount = eta.pairSums.size();

        return eta;
    }

    // ========== Phase 4: Build corrected emat and recompute DP ==========

    private EnergyMatrix buildCorrectedEmat(EtaCorrections eta) {
        EnergyMatrix corrected = new EnergyMatrix(branchMinimizingEmat);
        int numPos = rcs.getNumPos();

        // Add eta corrections to one-body terms
        for (int pos = 0; pos < numPos; pos++) {
            int numRC = corrected.getNumConfAtPos(pos);
            for (int rc = 0; rc < numRC; rc++) {
                double original = corrected.getOneBody(pos, rc);
                double correction = eta.getOneBodyEta(pos, rc);
                corrected.setOneBody(pos, rc, original + correction);
            }
        }

        // Add eta corrections to pair terms
        for (int pos1 = 0; pos1 < numPos; pos1++) {
            for (int pos2 = pos1 + 1; pos2 < numPos; pos2++) {
                if (!interactionGraph.hasEdge(pos1, pos2)) continue;
                int numRC1 = corrected.getNumConfAtPos(pos1);
                int numRC2 = corrected.getNumConfAtPos(pos2);
                for (int rc1 = 0; rc1 < numRC1; rc1++) {
                    for (int rc2 = 0; rc2 < numRC2; rc2++) {
                        double original = corrected.getPairwise(pos1, rc1, pos2, rc2);
                        double correction = eta.getPairEta(pos1, rc1, pos2, rc2);
                        corrected.setPairwise(pos1, rc1, pos2, rc2, original + correction);
                    }
                }
            }
        }

        return corrected;
    }

    /**
     * Recompute DP with corrected emat.
     * The corrected emat replaces the minimizing emat (upper bound side).
     * Returns logZ_corrected.
     */
    private CorrectedDPResult recomputeDP(EnergyMatrix correctedEmat, EtaCorrections eta) {
        // Re-initialize DP with corrected emat as the "min" side
        // Keep rigid emat the same (lower bound)
        RootedTreeEdge.postOrderInitIncremental(rootedRoot,
                branchRigidEmat, correctedEmat, interactionGraph, RT);
        DPCacheStats stats = computeCorrectedDPTables(eta);

        return new CorrectedDPResult(rootedRootEdge.getLogZUpper(0), stats);
    }

    private DPCacheStats computeCorrectedDPTables(EtaCorrections eta) {
        if (!dpCacheEnabled || dpCacheMaxEntries == 0) {
            RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
            return null;
        }

        DPCacheStats stats = new DPCacheStats();
        IdentityHashMap<RootedTreeEdge, String> edgeKeys = new IdentityHashMap<>();
        String namespace = makeCorrectedDPCacheNamespace();
        postOrderComputeCorrectedDPWithCache(rootedRoot, namespace, eta, edgeKeys, stats);
        return stats;
    }

    private void postOrderComputeCorrectedDPWithCache(RootedTreeNode node, String namespace,
                                                      EtaCorrections eta,
                                                      IdentityHashMap<RootedTreeEdge, String> edgeKeys,
                                                      DPCacheStats stats) {
        if (node == null) return;
        postOrderComputeCorrectedDPWithCache(node.getLeftChild(), namespace, eta, edgeKeys, stats);
        postOrderComputeCorrectedDPWithCache(node.getRightChild(), namespace, eta, edgeKeys, stats);

        RootedTreeEdge edge = node.getChildOfEdge();
        if (edge == null || !edge.getIsLambdaEdge()) return;

        if (!edge.hasDenseDPTable()) {
            edge.computeFullDP();
            edgeKeys.put(edge, buildCorrectedDPCacheKey(namespace, edge, eta, edgeKeys));
            stats.skippedLarge++;
            return;
        }

        String key = buildCorrectedDPCacheKey(namespace, edge, eta, edgeKeys);
        CachedDPTable cached = getCorrectedCachedDPTable(key, edge.getMArraySize());
        if (cached != null) {
            System.arraycopy(cached.lower, 0, edge.getLogZLower(), 0, cached.lower.length);
            System.arraycopy(cached.upper, 0, edge.getLogZUpper(), 0, cached.upper.length);
            stats.hits++;
        } else {
            edge.computeFullDP();
            if (putCorrectedCachedDPTable(key, edge.getLogZLower(), edge.getLogZUpper(), stats)) {
                stats.stores++;
            }
            stats.misses++;
        }
        edgeKeys.put(edge, key);
    }

    private CachedDPTable getCorrectedCachedDPTable(String key, int expectedLength) {
        synchronized (CORRECTED_DP_CACHE_LOCK) {
            CachedDPTable cached = CORRECTED_DP_CACHE.get(key);
            if (cached == null) return null;
            if (cached.lower.length != expectedLength || cached.upper.length != expectedLength) {
                removeCorrectedCachedDPTable(key, cached);
                return null;
            }
            return cached;
        }
    }

    private boolean putCorrectedCachedDPTable(String key, double[] lower, double[] upper,
                                              DPCacheStats stats) {
        long tableBytes = estimateDPTableBytes(lower.length);
        if (lower.length >= dpCacheSkipIfMStates
                || tableBytes > dpCacheMaxTableBytes
                || tableBytes > dpCacheMaxTotalBytes) {
            stats.skippedLarge++;
            return false;
        }

        synchronized (CORRECTED_DP_CACHE_LOCK) {
            CachedDPTable previous = CORRECTED_DP_CACHE.remove(key);
            if (previous != null) {
                correctedDPCacheBytes -= previous.bytes;
            }

            CORRECTED_DP_CACHE.put(key, new CachedDPTable(
                    Arrays.copyOf(lower, lower.length),
                    Arrays.copyOf(upper, upper.length),
                    tableBytes));
            correctedDPCacheBytes += tableBytes;

            while (CORRECTED_DP_CACHE.size() > dpCacheMaxEntries
                    || correctedDPCacheBytes > dpCacheMaxTotalBytes) {
                Iterator<String> it = CORRECTED_DP_CACHE.keySet().iterator();
                if (!it.hasNext()) break;
                String evictKey = it.next();
                CachedDPTable evicted = CORRECTED_DP_CACHE.get(evictKey);
                it.remove();
                if (evicted != null) {
                    correctedDPCacheBytes -= evicted.bytes;
                }
                stats.evictions++;
            }
        }
        return true;
    }

    private static void removeCorrectedCachedDPTable(String key, CachedDPTable cached) {
        CORRECTED_DP_CACHE.remove(key);
        correctedDPCacheBytes -= cached.bytes;
    }

    private static long estimateDPTableBytes(int mStates) {
        return 2L * (long) mStates * (long) Double.BYTES;
    }

    private static String formatBytes(long bytes) {
        if (bytes < 1024L) {
            return bytes + " B";
        }
        double value = bytes / 1024.0;
        String[] units = {"KiB", "MiB", "GiB", "TiB"};
        int unit = 0;
        while (value >= 1024.0 && unit < units.length - 1) {
            value /= 1024.0;
            unit++;
        }
        return String.format(Locale.ROOT, "%.1f %s", value, units[unit]);
    }

    private String formatDPCacheStats(DPCacheStats stats) {
        if (stats == null) return "";
        int size;
        long bytes;
        synchronized (CORRECTED_DP_CACHE_LOCK) {
            size = CORRECTED_DP_CACHE.size();
            bytes = correctedDPCacheBytes;
        }
        return ", correctedDPCacheHits=" + stats.hits
                + ", correctedDPCacheMisses=" + stats.misses
                + ", correctedDPCacheStores=" + stats.stores
                + ", correctedDPCacheSkippedLarge=" + stats.skippedLarge
                + ", correctedDPCacheEvictions=" + stats.evictions
                + ", correctedDPCacheSize=" + size
                + ", correctedDPCacheBytes=" + formatBytes(bytes);
    }

    private String makeCorrectedDPCacheNamespace() {
        StringBuilder sb = new StringBuilder();
        sb.append("corrected")
                .append("|conf=").append(System.identityHashCode(confSpace))
                .append("|rigid=").append(System.identityHashCode(branchRigidEmat))
                .append("|minBase=").append(System.identityHashCode(branchMinimizingEmat))
                .append("|rt=").append(Double.toHexString(RT))
                .append("|rootM=");
        appendPositionsOnly(sb, rootedRootEdge.getMPositionsSorted());
        sb.append("|rootL=");
        appendPositionsOnly(sb, rootedRootEdge.getLambdaPositionsSorted());
        sb.append("|graph=");
        for (int[] edge : interactionGraph.getEdgeList()) {
            sb.append(edge[0]).append('-').append(edge[1]).append(';');
        }
        return sb.toString();
    }

    private String buildCorrectedDPCacheKey(String namespace, RootedTreeEdge edge,
                                            EtaCorrections eta,
                                            IdentityHashMap<RootedTreeEdge, String> edgeKeys) {
        StringBuilder sb = new StringBuilder(namespace);
        sb.append("|M=");
        appendPositionsWithRCs(sb, edge.getMPositionsSorted());
        sb.append("|L=");
        appendPositionsWithRCs(sb, edge.getLambdaPositionsSorted());
        sb.append("|eta=");
        appendLocalEtaSignature(sb, edge, eta);
        sb.append("|F=");
        if (edge.getFset() != null) {
            List<String> childKeys = new ArrayList<>();
            for (RootedTreeEdge child : edge.getFset()) {
                String childKey = edgeKeys.get(child);
                childKeys.add(childKey == null ? "missing" : childKey);
            }
            Collections.sort(childKeys);
            for (String childKey : childKeys) {
                sb.append('{').append(childKey).append("},");
            }
        }
        return sb.toString();
    }

    private void appendPositionsOnly(StringBuilder sb, int[] positions) {
        for (int pos : positions) {
            sb.append(pos).append(',');
        }
    }

    private void appendPositionsWithRCs(StringBuilder sb, int[] positions) {
        for (int pos : positions) {
            sb.append(pos).append('[');
            int n = rcs.getNum(pos);
            for (int i = 0; i < n; i++) {
                if (i > 0) sb.append(',');
                sb.append(rcs.get(pos, i));
            }
            sb.append("];");
        }
    }

    private void appendLocalEtaSignature(StringBuilder sb, RootedTreeEdge edge, EtaCorrections eta) {
        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        sb.append("ob=");
        for (int pos : lambdaPositions) {
            sb.append(pos).append('[');
            int n = rcs.getNum(pos);
            for (int i = 0; i < n; i++) {
                if (i > 0) sb.append(',');
                int rc = rcs.get(pos, i);
                appendDoubleBits(sb, eta.getOneBodyEta(pos, rc));
            }
            sb.append("];");
        }

        sb.append("ll=");
        for (int i = 0; i < lambdaPositions.length; i++) {
            int pos1 = lambdaPositions[i];
            for (int j = i + 1; j < lambdaPositions.length; j++) {
                int pos2 = lambdaPositions[j];
                if (!interactionGraph.hasEdge(pos1, pos2)) continue;
                appendPairEtaSignature(sb, eta, pos1, pos2);
            }
        }

        sb.append("lm=");
        for (int pos1 : lambdaPositions) {
            for (int pos2 : mPositions) {
                if (!interactionGraph.hasEdge(pos1, pos2)) continue;
                appendPairEtaSignature(sb, eta, pos1, pos2);
            }
        }
    }

    private void appendPairEtaSignature(StringBuilder sb, EtaCorrections eta, int pos1, int pos2) {
        int a = Math.min(pos1, pos2);
        int b = Math.max(pos1, pos2);
        sb.append(a).append('-').append(b).append('[');
        int n1 = rcs.getNum(a);
        int n2 = rcs.getNum(b);
        for (int i = 0; i < n1; i++) {
            int rc1 = rcs.get(a, i);
            if (i > 0) sb.append(';');
            for (int j = 0; j < n2; j++) {
                if (j > 0) sb.append(',');
                int rc2 = rcs.get(b, j);
                appendDoubleBits(sb, eta.getPairEta(a, rc1, b, rc2));
            }
        }
        sb.append("];");
    }

    private void appendDoubleBits(StringBuilder sb, double value) {
        sb.append(Long.toHexString(Double.doubleToLongBits(value))).append(':');
    }

    // ========== Phase 5: Compute residuals and PAC bound ==========

    /**
     * Unbiased estimator with control variate for variance reduction.
     *
     * Basic estimator: Z = Z_min * E_p[phi(c)], phi(c) = exp(-g(c)/kT), samples c ~ p ∝ exp(-E_min/kT)
     *
     * Control variate: alpha(c) = exp(-eta(c)/kT) has known expectation
     *   E_p[alpha] = Z_corrected / Z_min  (computable from DP)
     *
     * Variance-reduced estimator:
     *   phi_cv(c) = phi(c) - beta * (alpha(c) - E_p[alpha])
     *   E[phi_cv] = E[phi]  (unbiased)
     *   beta* = Cov(phi, alpha) / Var(alpha)
     *
     * Then Z = Z_min * mean(phi_cv), with Bernstein bound on phi_cv.
     */
    private void computePACBound(List<CCDResult> ccdResults,
                                  EnergyMatrix correctedEmat,
                                  double logZCorrected) {
        int N = ccdResults.size();
        if (N == 0) {
            setZeroBounds("PAC has no valid CCD samples");
            return;
        }

        double[] logPhiValues = new double[N];    // log phi(c) = -g/kT
        double[] logAlphaValues = new double[N];  // log alpha(c) = -eta/kT
        double[] residuals = new double[N];    // residual = g - eta (for diagnostics)

        for (int i = 0; i < N; i++) {
            CCDResult result = ccdResults.get(i);
            double g = result.eTrue - result.eMin; // g(c) = E_true - E_min
            double eta = computeFullConfPairwiseEnergy(result.conf, correctedEmat)
                    - result.eMin;  // eta(c) = E_corrected - E_min
            residuals[i] = g - eta;
            logPhiValues[i] = -g / RT;
            logAlphaValues[i] = -eta / RT;
        }

        // E_p[alpha] = Z_corrected / Z_min (exact from DP)
        // logZMinOriginal was saved before Phase 4 modified the DP tables.
        double logExpectedAlpha = logZCorrected - logZMinOriginal;
        if (!Double.isFinite(logExpectedAlpha)) {
            setZeroBounds("corrected/min DP ratio is NaN");
            return;
        }

        double weightShift = logExpectedAlpha;
        for (int i = 0; i < N; i++) {
            weightShift = Math.max(weightShift, logPhiValues[i]);
            weightShift = Math.max(weightShift, logAlphaValues[i]);
        }
        if (!Double.isFinite(weightShift)) {
            setZeroBounds("PAC sample weights are non-finite");
            return;
        }
        double expectedAlpha = Math.exp(logExpectedAlpha - weightShift);

        double[] phiValues = new double[N];    // shifted by exp(-weightShift)
        double[] alphaValues = new double[N];  // shifted by exp(-weightShift)
        for (int i = 0; i < N; i++) {
            phiValues[i] = Math.exp(logPhiValues[i] - weightShift);
            alphaValues[i] = Math.exp(logAlphaValues[i] - weightShift);
        }

        // Compute optimal beta for control variate: beta* = Cov(phi, alpha) / Var(alpha)
        double sumPhi = 0, sumAlpha = 0;
        for (int i = 0; i < N; i++) {
            sumPhi += phiValues[i];
            sumAlpha += alphaValues[i];
        }
        double meanPhi = sumPhi / N;
        double meanAlpha = sumAlpha / N;

        double covPhiAlpha = 0, varAlpha = 0;
        for (int i = 0; i < N; i++) {
            double dPhi = phiValues[i] - meanPhi;
            double dAlpha = alphaValues[i] - meanAlpha;
            covPhiAlpha += dPhi * dAlpha;
            varAlpha += dAlpha * dAlpha;
        }
        if (N > 1) {
            covPhiAlpha /= (N - 1);
            varAlpha /= (N - 1);
        } else {
            covPhiAlpha = 0.0;
            varAlpha = 0.0;
        }

        double beta = (varAlpha > 1e-30) ? covPhiAlpha / varAlpha : 0.0;

        // Compute control-variate-adjusted values: phi_cv = phi - beta*(alpha - E[alpha])
        double[] phiCV = new double[N];
        double sumCV = 0, sumCV2 = 0;
        double minCV = Double.MAX_VALUE, maxCV = -Double.MAX_VALUE;
        for (int i = 0; i < N; i++) {
            phiCV[i] = phiValues[i] - beta * (alphaValues[i] - expectedAlpha);
            sumCV += phiCV[i];
            sumCV2 += phiCV[i] * phiCV[i];
            minCV = Math.min(minCV, phiCV[i]);
            maxCV = Math.max(maxCV, phiCV[i]);
        }

        double meanCV = sumCV / N;
        double varCV = N > 1 ? Math.max(0, (sumCV2 - N * meanCV * meanCV) / (N - 1)) : 0.0;
        double rangeCV = maxCV - minCV;

        // Residual diagnostics (for reporting)
        double sumRes = 0, sumRes2 = 0;
        for (int i = 0; i < N; i++) {
            sumRes += residuals[i];
            sumRes2 += residuals[i] * residuals[i];
        }
        meanResidual = sumRes / N;
        stdResidual = N > 1
                ? Math.sqrt(Math.max(0, (sumRes2 - N * meanResidual * meanResidual) / (N - 1)))
                : 0.0;

        // Also compute basic phi stats for reporting. meanPsi is unshifted when
        // representable as a double; log10 bounds above are the authoritative
        // diagnostics for extreme cases.
        meanPsi = shiftedMeanToDouble(meanPhi, weightShift);
        double varPhi = Math.max(0, (phiValues.length > 1 ?
                (Arrays.stream(phiValues).map(x -> (x - meanPhi) * (x - meanPhi)).sum()) / (N - 1) : 0));
        cvPsi = (meanPhi > 0) ? Math.sqrt(varPhi) / meanPhi : Double.MAX_VALUE;

        // Variance reduction ratio
        double vrRatio = (varPhi > 1e-30) ? varCV / varPhi : 1.0;

        // Bernstein bound on phi_cv (one-sided: we only need lower bound for Z_lower)
        double deltaPer = delta / 3.0; // union bound over 3 pfuncs
        double boundDelta = solveBernsteinDelta(N, varCV, rangeCV, deltaPer);
        if (!Double.isFinite(meanCV) || !Double.isFinite(varCV)
                || !Double.isFinite(rangeCV) || !Double.isFinite(boundDelta)) {
            setZeroBounds("PAC statistics are non-finite");
            return;
        }

        // Z = Z_min * E_p[phi] = Z_min * E_p[phi_cv]  (unbiased)
        // The confidence interval is still in the shifted scale. Convert by
        // adding weightShift in log-space instead of multiplying tiny doubles.
        double cvBarLower = Math.max(0, meanCV - boundDelta);
        double cvBarUpper = meanCV + boundDelta;
        if (!Double.isFinite(cvBarLower) || !Double.isFinite(cvBarUpper) || cvBarUpper <= 0) {
            setZeroBounds("PAC confidence interval is non-finite");
            return;
        }

        logZLowerPAC = cvBarLower > 0
                ? logZMinOriginal + weightShift + Math.log(cvBarLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = logZMinOriginal + weightShift + Math.log(cvBarUpper);

        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);

        if (Double.isFinite(logZLowerPAC) && Double.isFinite(logZUpperPAC)) {
            epsilon = 1.0 - Math.exp(logZLowerPAC - logZUpperPAC);
            if (epsilon < 0.0 && epsilon > -1e-12) epsilon = 0.0;
            if (epsilon > 1.0 && epsilon < 1.0 + 1e-12) epsilon = 1.0;
        } else {
            epsilon = 1.0;
        }

        System.out.println("[PAC] Bernstein: Delta=" + String.format("%.6f", boundDelta)
                + ", meanCV=" + String.format("%.6f", meanCV)
                + ", rangeCV=" + String.format("%.6f", rangeCV)
                + ", varCV=" + String.format("%.6f", varCV)
                + ", varPhi=" + String.format("%.6f", varPhi)
                + ", vrRatio=" + String.format("%.4f", vrRatio)
                + ", beta=" + String.format("%.4f", beta)
                + ", logEAlpha=" + String.format("%.6f", logExpectedAlpha)
                + ", weightShift=" + String.format("%.6f", weightShift)
                + ", log10ZLower=" + formatLog10(logZLowerPAC)
                + ", log10ZUpper=" + formatLog10(logZUpperPAC)
                + ", N=" + N);
    }

    /**
     * Solve one-sided Bernstein inequality for Delta:
     * exp(-N*Delta^2 / (2*s^2 + 2*Delta*range/3)) = delta_target
     *
     * One-sided (not two-sided) since we only need the lower bound on E[phi_cv].
     */
    private double solveBernsteinDelta(int N, double s2, double range, double deltaTarget) {
        double logRHS = Math.log(1.0 / deltaTarget); // one-sided: ln(1/delta), not ln(2/delta)

        // Bisection on Delta
        double lo = 0, hi = range + 10 * Math.sqrt(s2); // generous upper bound
        for (int iter = 0; iter < 100; iter++) {
            double mid = (lo + hi) / 2;
            double lhs = (double) N * mid * mid / (2 * s2 + 2 * mid * range / 3);
            if (lhs < logRHS) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return hi;
    }

    // ========== Utility methods ==========

    private double shiftedMeanToDouble(double shiftedMean, double shift) {
        if (shiftedMean <= 0.0) return 0.0;
        return Math.exp(shift + Math.log(shiftedMean));
    }

    private static String formatLog10(double logZ) {
        if (logZ == Double.NEGATIVE_INFINITY) return "-inf";
        if (logZ == Double.POSITIVE_INFINITY) return "+inf";
        if (Double.isNaN(logZ)) return "NaN";
        return String.format(Locale.ROOT, "%.6f", logZ / Math.log(10.0));
    }

    /**
     * Compute full conformation pairwise energy from an energy matrix.
     * Only sums over pairs in the interaction graph (sparse mode).
     */
    private double computeFullConfPairwiseEnergy(int[] conf, EnergyMatrix emat) {
        double energy = emat.getConstTerm();
        int numPos = conf.length;
        for (int i = 0; i < numPos; i++) {
            if (conf[i] < 0) continue;
            energy += emat.getOneBody(i, conf[i]);
            for (int j = i + 1; j < numPos; j++) {
                if (conf[j] < 0) continue;
                if (!interactionGraph.hasEdge(i, j)) continue;
                energy += emat.getPairwise(i, conf[i], j, conf[j]);
            }
        }
        return energy;
    }

    /**
     * Build sparse residue interactions for a full conformation.
     */
    private ResidueInteractions makeSparseFullConfInters(int[] conf) {
        ResidueInteractions inters = new ResidueInteractions();
        for (int pos = 0; pos < conf.length; pos++) {
            if (conf[pos] < 0) continue;
            inters.addAll(minimizingEcalc.makeSingleInters(pos, conf[pos]));
        }
        for (int pos1 = 0; pos1 < conf.length; pos1++) {
            if (conf[pos1] < 0) continue;
            for (int pos2 = pos1 + 1; pos2 < conf.length; pos2++) {
                if (conf[pos2] < 0 || !interactionGraph.hasEdge(pos1, pos2)) continue;
                inters.addAll(minimizingEcalc.makePairInters(pos1, conf[pos1], pos2, conf[pos2]));
            }
        }
        if (minimizingEcalc.addShellInters) {
            inters.addAll(minimizingEcalc.makeShellInters());
        }
        return inters;
    }

    /**
     * Convert a log-space value to BigDecimal.
     */
    private BigDecimal bigExpFromLog(double logVal) {
        if (Double.isInfinite(logVal) && logVal < 0) return BigDecimal.ZERO;
        // exp(logVal) = exp(logVal/ln10 * ln10) = 10^(logVal/ln10)
        // Use BoltzmannCalculator's precision
        return bc.calc(-logVal * RT); // bc.calc(energy) = exp(-energy/RT)
        // Actually, we want exp(logVal) directly. bc.calc(e) = exp(-e/RT)
        // So bc.calc(-logVal * RT) = exp(-(-logVal*RT)/RT) = exp(logVal). Correct.
    }

    private void setZeroBounds(String reason) {
        System.out.println("[PAC] WARNING: " + reason + "; returning zero Z bounds");
        zLower = BigDecimal.ZERO;
        zUpper = BigDecimal.ZERO;
        logZLowerPAC = Double.NEGATIVE_INFINITY;
        logZUpperPAC = Double.NEGATIVE_INFINITY;
        epsilon = 1.0;
        meanPsi = 0.0;
        varPsi = 0.0;
        cvPsi = Double.MAX_VALUE;
        meanResidual = 0.0;
        stdResidual = 0.0;
    }

    // ========== Result getters ==========

    public BigDecimal getZLower() { return zLower; }
    public BigDecimal getZUpper() { return zUpper; }
    public double getEpsilon() { return epsilon; }
    public int getTotalCCDCalls() { return totalCCDCalls; }
    public double getMeanPsi() { return meanPsi; }
    public double getCvPsi() { return cvPsi; }
    public double getMeanResidual() { return meanResidual; }
    public double getStdResidual() { return stdResidual; }
}
