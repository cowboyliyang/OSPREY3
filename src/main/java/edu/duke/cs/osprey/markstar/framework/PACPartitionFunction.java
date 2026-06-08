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
import java.io.File;
import java.io.PrintWriter;

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
    private static final String PAC_TARGET_EPSILON_PROPERTY = "branchmarkstar.pac.targetEpsilon";
    private static final String PAC_SAMPLING_BATCHED_PROPERTY = "branchmarkstar.pac.sampling.batched";
    private static final String PAC_SAMPLING_PARALLEL_PROPERTY = "branchmarkstar.pac.sampling.parallel";
    private static final String PAC_SAMPLING_THREADS_PROPERTY = "branchmarkstar.pac.sampling.threads";
    private static final String PAC_SAMPLING_LARGE_LAMBDA_PROPERTY = "branchmarkstar.pac.sampling.largeLambdaThreshold";
    private static final String PAC_SAMPLING_PROGRESS_PROPERTY = "branchmarkstar.pac.sampling.progress";
    private static final String PAC_SAMPLING_GPU_PROPERTY = "branchmarkstar.pac.sampling.gpu";
    // Tier-1 two-stage estimator: learn eta on a train split (from p_m), solve the
    // corrected DP ONCE, then sample the estimation set from p_eta and reweight only
    // the residual xi = E_true - E_eta. Sample size N* is computed up-front from a
    // pilot by inverting empirical Bernstein (no per-round DP re-solve / restore).
    private static final String PAC_TRAIN_SAMPLES_PROPERTY = "branchmarkstar.pac.trainSamples";
    private static final String PAC_PILOT_SAMPLES_PROPERTY = "branchmarkstar.pac.pilotSamples";
    private static final String PAC_MAX_EST_SAMPLES_PROPERTY = "branchmarkstar.pac.maxEstSamples";
    private static final String PAC_NSTAR_INFLATE_PROPERTY = "branchmarkstar.pac.nstarInflate";
    private static final String PAC_RESIDUAL_BOUND_PROPERTY = "branchmarkstar.pac.residualBound";
    private static final String DP_PARALLEL_THREADS_PROPERTY = "branchmarkstar.dp.parallel.threads";
    private static final String DP_CACHE_ENABLED_PROPERTY = "branchmarkstar.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "branchmarkstar.dp.cache.maxEntries";
    private static final String DP_CACHE_MAX_TABLE_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTableBytes";
    private static final String DP_CACHE_MAX_TOTAL_BYTES_PROPERTY = "branchmarkstar.dp.cache.maxTotalBytes";
    private static final String DP_CACHE_SKIP_IF_M_STATES_PROPERTY = "branchmarkstar.dp.cache.skipIfMStates";
    private static final int DEFAULT_SAMPLES = 500;
    private static final double DEFAULT_CONFIDENCE = 0.05; // delta = 0.05 => 95% confidence
    private static final double DEFAULT_TARGET_EPSILON = 0.683;
    private static final int DEFAULT_SAMPLING_LARGE_LAMBDA = 65_536;
    private static final double DEFAULT_TRAIN_FRACTION = 0.5;
    private static final double DEFAULT_PILOT_FRACTION = 0.1;
    private static final double DEFAULT_NSTAR_INFLATE = 1.3;
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
    private final double targetEpsilon;
    private final boolean batchedSampling;
    private final boolean gpuSampling;
    private final int samplingThreads;
    private final int samplingLargeLambdaThreshold;
    private final boolean samplingProgress;
    private final int trainSamples;
    private final int pilotSamples;
    private final int maxEstSamples;
    private double clipLogCap = Double.NaN; // clip threshold fixed on the pilot (PAC: chosen before S2)
    private double logZMinDet = Double.POSITIVE_INFINITY; // log q_m: deterministic, assumption-free upper bound on log q
    private final double nstarInflate;
    private final double residualBoundKcal;
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
        double defaultTarget = Double.isFinite(requestedTargetEpsilon) && requestedTargetEpsilon > 0.0
                ? requestedTargetEpsilon
                : DEFAULT_TARGET_EPSILON;
        this.targetEpsilon = Math.max(0.0,
                getConfigDouble(PAC_TARGET_EPSILON_PROPERTY, defaultTarget));
        this.batchedSampling = getConfigBoolean(PAC_SAMPLING_BATCHED_PROPERTY, true);
        this.gpuSampling = getConfigBoolean(PAC_SAMPLING_GPU_PROPERTY, false);
        this.samplingThreads = resolveConfiguredSamplingThreads();
        this.samplingLargeLambdaThreshold = Math.max(1,
                getConfigInteger(PAC_SAMPLING_LARGE_LAMBDA_PROPERTY,
                        DEFAULT_SAMPLING_LARGE_LAMBDA));
        this.samplingProgress = getConfigBoolean(PAC_SAMPLING_PROGRESS_PROPERTY, true);
        int twoStageBudget = Math.max(6, this.numSamples);
        int defaultTrainSamples = Math.max(2, (int) Math.floor(twoStageBudget * DEFAULT_TRAIN_FRACTION));
        int defaultPilotSamples = Math.max(2, (int) Math.floor(twoStageBudget * DEFAULT_PILOT_FRACTION));
        int defaultMaxEstSamples = Math.max(2, twoStageBudget - defaultTrainSamples - defaultPilotSamples);
        this.trainSamples = Math.max(2, getConfigInteger(PAC_TRAIN_SAMPLES_PROPERTY, defaultTrainSamples));
        this.pilotSamples = Math.max(2, getConfigInteger(PAC_PILOT_SAMPLES_PROPERTY, defaultPilotSamples));
        this.maxEstSamples = Math.max(2,
                getConfigInteger(PAC_MAX_EST_SAMPLES_PROPERTY, defaultMaxEstSamples));
        this.nstarInflate = Math.max(1.0, getConfigDouble(PAC_NSTAR_INFLATE_PROPERTY, DEFAULT_NSTAR_INFLATE));
        double configuredResidualBound = getConfigDouble(PAC_RESIDUAL_BOUND_PROPERTY, Double.NaN);
        this.residualBoundKcal = configuredResidualBound >= 0.0 ? configuredResidualBound : Double.NaN;
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

        System.out.println("[PAC-2stage] Strict pilot holdout is enabled"
                + ", samples(train/pilot/maxEst)="
                + trainSamples + "/" + pilotSamples + "/" + maxEstSamples
                + ", residualBound="
                + (Double.isFinite(residualBoundKcal)
                ? String.format("%.4f kcal/mol", residualBoundKcal)
                : "not configured"));
        this.logZMinDet = logZMin; // q_m >= q always (E_m <= E_true): assumption-free upper-bound fallback
        runTwoStagePAC(startTime, logZRigid);

        return epsilon;
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

    // ========== Tier-1 two-stage estimator (p_eta proposal + residual reweight) ==========

    /**
     * Two-stage PAC estimator.
     *
     *   Z = q_eta * E_{p_eta}[ exp(-xi/RT) ],   xi = E_true - E_eta   (exact identity)
     *
     * Stage A (train, from p_m): learn eta, build the corrected emat and solve the
     * corrected DP ONCE -> p_eta (proposal) and q_eta = exp(logZCorrected).
     * Stage B (est, from p_eta): draw a held-out pilot, estimate the residual-leg
     * variance, invert empirical Bernstein to get the required N* up-front, then
     * draw a fresh estimation set of that fixed size in a single shot. No per-round
     * DP re-solve and no proposal restore.
     *
     * Sample-splitting: eta/p_eta are functions of the train split only, so the
     * estimation samples (and their residual weights) are independent of eta and of
     * the pilot given the train split, so conditioned on the pilot the final N is
     * fixed, exactly as the empirical-Bernstein PAC guarantee requires.
     */
    private void runTwoStagePAC(long startTime, double logZRigid) {
        Random rng = new Random(42); // reproducible

        // ---- Stage A: train from p_m -> eta, corrected DP (p_eta, q_eta) ----
        long tA = System.currentTimeMillis();
        List<int[]> trainConfs = sampleConformationsFromDP(trainSamples, rng);
        List<CCDResult> trainCCD = runParallelCCD(trainConfs);
        System.out.println("[PAC-2stage] Stage A train: " + trainCCD.size()
                + " CCD from p_m in " + (System.currentTimeMillis() - tA) + " ms");
        if (trainCCD.isEmpty()) {
            setZeroBounds("two-stage: no valid train samples");
            printFinalSummary(startTime);
            return;
        }

        EtaCorrections eta = extractEtaCorrections(trainCCD);
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        maybeDumpTrainingSamples(trainCCD, eta, correctedEmat);
        long tDP = System.currentTimeMillis();
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta); // loads p_eta into the tree
        double logZCorrected = correctedDP.logZCorrected;
        System.out.println("[PAC-2stage] corrected DP solved ONCE: logZ_corrected="
                + String.format("%.4f", logZCorrected)
                + " (improvement over logZ_rigid: " + String.format("%.4f", logZCorrected - logZRigid) + ")"
                + " in " + (System.currentTimeMillis() - tDP) + " ms"
                + ", oneBodyTerms=" + eta.oneBodyCount + ", pairTerms=" + eta.pairCount);
        if (!Double.isFinite(logZCorrected)) {
            setZeroBounds("two-stage: corrected logZ is non-finite");
            printFinalSummary(startTime);
            return;
        }

        // ---- Stage B pilot from p_eta -> estimate residual-leg stats -> N* ----
        // The tree now holds the corrected DP, so sampleConformationsFromDP draws p_eta.
        long tPilot = System.currentTimeMillis();
        int nPilot = Math.min(pilotSamples, maxEstSamples);
        List<CCDResult> pilotCCD = new ArrayList<>(runParallelCCD(sampleConformationsFromDP(nPilot, rng)));
        System.out.println("[PAC-2stage] Stage B pilot: " + pilotCCD.size()
                + " CCD from p_eta in " + (System.currentTimeMillis() - tPilot) + " ms");
        if (pilotCCD.isEmpty()) {
            setZeroBounds("two-stage: no valid pilot samples");
            printFinalSummary(startTime);
            return;
        }

        // ---- distribution-shift refinement (smart iteration) ----
        // eta is learned on p_m, but the corrected p_eta can concentrate mass on
        // rotamers p_m never sampled (eta untrained there -> systematic residual,
        // meanW collapse). When the p_eta pilot reveals this, re-learn eta on
        // train(p_m) + pilot(p_eta) -- which now covers the shifted region -- and
        // re-solve the DP once. The old pilot is reused (no new train batch), so the
        // only added cost is one DP solve plus a fresh pilot. The final estimation
        // set is still drawn fresh from the refined p_eta, so sample-splitting (and
        // thus the PAC guarantee) is preserved.
        int extraRefineCCD = 0;
        boolean iterate = getConfigBoolean("branchmarkstar.pac.iterate", true);
        double collapseThresh = getConfigDouble("branchmarkstar.pac.iterate.meanWThreshold", 0.3);
        double driftFracThresh = getConfigDouble("branchmarkstar.pac.iterate.driftFraction", 0.2);
        int minTrainCount = getConfigInteger("branchmarkstar.pac.iterate.minTrainCount", 5);
        int maxRounds = getConfigInteger("branchmarkstar.pac.iterate.maxRounds", 4);
        if (iterate) {
            // Multi-round (EM-style) distribution-shift refinement. eta is learned on p_m,
            // but the corrected p_eta can concentrate mass on rotamers p_m never sampled
            // (eta untrained there -> systematic residual, meanW collapse). Each round folds
            // the current pilot (drawn from the current p_eta) into a growing TRAINING pool,
            // re-learns eta on the pool -- which now covers the shifted region -- re-solves the
            // DP, and draws a FRESH pilot from the refined p_eta. Because re-learning eta moves
            // p_eta (a moving target), one round rarely suffices; we iterate until drift/collapse
            // clears or maxRounds. Two OR-combined detectors:
            //  (1) drift: fraction of pilot mass on (pos,rc) cells eta trained < minTrainCount
            //      times -- the root cause, robust to rare huge-weight (negative-xi) outliers.
            //  (2) collapse: absolute meanW near 0 (backstop for pure positive-xi cases).
            // The final estimation set is still drawn fresh from the final p_eta, so
            // sample-splitting (and the PAC guarantee) is preserved.
            List<CCDResult> trainPool = new ArrayList<>(trainCCD);
            for (int round = 1; round <= maxRounds; round++) {
                double meanWPilot = estimateMeanW(pilotCCD, correctedEmat);
                double driftFrac = pilotUndertrainedFraction(pilotCCD, eta, minTrainCount);
                boolean collapse = Double.isFinite(meanWPilot) && meanWPilot < collapseThresh;
                boolean drift = driftFrac >= driftFracThresh;
                System.out.println("[PAC-2stage] refinement check (round " + round + "/" + maxRounds
                        + "): meanWPilot=" + String.format("%.4f", meanWPilot)
                        + ", driftFrac=" + String.format("%.3f", driftFrac)
                        + " (collapse=" + collapse + ", drift=" + drift + ")");
                if (!(collapse || drift)) break;

                // fold the current pilot (from the current p_eta) into the training pool
                trainPool.addAll(pilotCCD);
                extraRefineCCD += pilotCCD.size();
                EtaCorrections eta2 = extractEtaCorrections(trainPool);
                EnergyMatrix correctedEmat2 = buildCorrectedEmat(eta2);
                CorrectedDPResult correctedDP2 = recomputeDP(correctedEmat2, eta2); // loads refined p_eta into the tree
                if (!Double.isFinite(correctedDP2.logZCorrected)) {
                    System.out.println("[PAC-2stage] refinement round " + round
                            + ": DP non-finite, keeping previous eta");
                    break;
                }
                eta = eta2;
                correctedEmat = correctedEmat2;
                logZCorrected = correctedDP2.logZCorrected;
                long tRef = System.currentTimeMillis();
                List<CCDResult> pilotNew =
                        new ArrayList<>(runParallelCCD(sampleConformationsFromDP(nPilot, rng)));
                double meanWNew = estimateMeanW(pilotNew, correctedEmat);
                System.out.println("[PAC-2stage] distribution-shift refinement round " + round
                        + ": meanW " + String.format("%.4f", meanWPilot) + " -> " + String.format("%.4f", meanWNew)
                        + ", logZcorr -> " + String.format("%.4f", logZCorrected)
                        + ", trainPool=" + trainPool.size()
                        + " (1 DP solve + fresh pilot in " + (System.currentTimeMillis() - tRef) + " ms)");
                if (pilotNew.isEmpty()) break;
                pilotCCD = pilotNew;
            }
        }

        int nStar = computeRequiredEstSamples(pilotCCD, correctedEmat);
        nStar = Math.max(2, Math.min(nStar, maxEstSamples));

        long tEst = System.currentTimeMillis();
        List<CCDResult> estCCD = runParallelCCD(sampleConformationsFromDP(nStar, rng));
        System.out.println("[PAC-2stage] Stage B fresh estimation: " + estCCD.size() + "/" + nStar
                + " CCD from p_eta in " + (System.currentTimeMillis() - tEst) + " ms"
                + " (pilot held out)");
        totalCCDCalls = trainCCD.size() + pilotCCD.size() + estCCD.size() + extraRefineCCD;

        // ---- clip threshold (if enabled) is fixed on the held-out pilot, BEFORE the
        // estimation weights S2 are observed, so the clipped weights stay i.i.d. given c ----
        if (getConfigBoolean("branchmarkstar.pac.clip", false)) {
            double clipQuantile = getConfigDouble("branchmarkstar.pac.clipQuantile", 0.9);
            double[] plogW = new double[pilotCCD.size()];
            for (int i = 0; i < plogW.length; i++) {
                CCDResult r = pilotCCD.get(i);
                plogW[i] = -(r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat)) / RT;
            }
            clipLogCap = quantile(plogW, clipQuantile);
            System.out.println("[PAC-2stage] clip threshold from pilot: logCap="
                    + String.format("%.4f", clipLogCap) + " (clipQ=" + clipQuantile
                    + ", nPilot=" + plogW.length + ")");
        }

        // ---- residual-leg empirical-Bernstein bound: Z = q_eta * E_{p_eta}[exp(-xi/RT)] ----
        computePACBoundResidual(estCCD, correctedEmat, logZCorrected);
        printFinalSummary(startTime);
    }

    /**
     * Compute the estimation sample size N* needed to reach targetEpsilon, by
     * inverting the empirical Bernstein bound on the residual-leg weight
     * w = exp(-xi/RT) using pilot estimates of its mean/variance and a
     * deterministic residual-weight range.
     *
     * epsilon = 2*Delta/(meanW + Delta)  =>  target Delta = eps*meanW/(2-eps).
     * Returns the smallest N in [2, maxEstSamples] whose Bernstein Delta <= target.
     * The pilot variance is inflated by nstarInflate for a conservative choice.
     */
    private int computeRequiredEstSamples(List<CCDResult> pilot, EnergyMatrix correctedEmat) {
        int n = pilot.size();
        if (n <= 1) return maxEstSamples;
        boolean clip = getConfigBoolean("branchmarkstar.pac.clip", false);
        double shift = Double.NEGATIVE_INFINITY;
        double[] logW = new double[n];
        for (int i = 0; i < n; i++) {
            CCDResult r = pilot.get(i);
            double xi = r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat);
            logW[i] = -xi / RT;
            shift = Math.max(shift, logW[i]);
        }
        if (!Double.isFinite(shift)) return maxEstSamples;
        // Under clipping, mirror the final estimator: shift by the clip threshold
        // logCap (a logW quantile) and use the deterministic range = 1.
        if (clip) {
            double clipQuantile = getConfigDouble("branchmarkstar.pac.clipQuantile", 0.9);
            shift = quantile(logW, clipQuantile);
            if (!Double.isFinite(shift)) return maxEstSamples;
        }
        double sum = 0, sum2 = 0, mn = Double.MAX_VALUE, mx = -Double.MAX_VALUE;
        double sumU = 0, sumU2 = 0; int nClip = 0;
        for (int i = 0; i < n; i++) {
            double w = clip ? Math.exp(Math.min(logW[i], shift) - shift) : Math.exp(logW[i] - shift);
            sum += w; sum2 += w * w;
            mn = Math.min(mn, w); mx = Math.max(mx, w);
            if (clip) {
                double u = Math.max(0.0, Math.exp(logW[i] - shift) - 1.0); // clip-excess in units of c
                sumU += u; sumU2 += u * u;
                if (logW[i] > shift) nClip++;
            }
        }
        double mean = sum / n;
        if (!(mean > 0)) return maxEstSamples;
        double var = Math.max(0.0, (sum2 - n * mean * mean) / (n - 1)) * nstarInflate;
        double sampleRange = mx - mn;
        double deltaPer = delta / 6.0; // 2 tails x 3 pfuncs

        int nStar;
        if (clip) {
            // Predict the FULL clipped epsilon at candidate N, INCLUDING the clip-bias
            // upper term; otherwise N is undersized and biasUpper dominates the real eps.
            double meanU = sumU / n;
            double varU = (n > 1 ? Math.max(0.0, (sumU2 - n * meanU * meanU) / (n - 1)) : 0.0) * nstarInflate;
            double pHat = (double) nClip / n;
            double uRange = (Double.isFinite(residualBoundKcal) && residualBoundKcal / RT - shift > 0)
                    ? Math.expm1(residualBoundKcal / RT - shift) : 0.0;
            nStar = solveClipSizing(mean, var, meanU, varU, uRange, pHat, deltaPer);
            System.out.println("[PAC-2stage] N* calc(clip): pilot=" + n
                    + ", meanWc=" + String.format("%.6f", mean)
                    + ", varWc=" + String.format("%.6f", var)
                    + ", meanExcess=" + String.format("%.6f", meanU)
                    + ", uRange=" + String.format("%.4f", uRange)
                    + ", pClip=" + String.format("%.3f", pHat)
                    + ", targetEps=" + String.format("%.4f", targetEpsilon)
                    + ", predEps@N*=" + String.format("%.4f",
                        predictClipEps(nStar, mean, var, meanU, varU, uRange, pHat, deltaPer))
                    + " -> N*_est=" + nStar + " (cap " + maxEstSamples + ")");
            return nStar;
        }

        double range = resolveResidualRangeForBernstein(shift, sampleRange, "N* sizing");
        double targetDelta = targetEpsilon * mean / (2.0 - targetEpsilon);
        if (solveBernsteinDelta(maxEstSamples, var, range, deltaPer) > targetDelta) {
            nStar = maxEstSamples; // even the cap cannot reach target
        } else {
            int lo = 2, hi = maxEstSamples;
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (solveBernsteinDelta(mid, var, range, deltaPer) <= targetDelta) hi = mid;
                else lo = mid + 1;
            }
            nStar = lo;
        }
        System.out.println("[PAC-2stage] N* calc: pilot=" + n
                + ", meanW=" + String.format("%.6f", mean)
                + ", varW(inflated)=" + String.format("%.6f", var)
                + ", rangeW=" + String.format("%.6f", range)
                + ", sampleRangeW=" + String.format("%.6f", sampleRange)
                + ", targetEps=" + String.format("%.4f", targetEpsilon)
                + ", targetDelta=" + String.format("%.6f", targetDelta)
                + " -> N*_est=" + nStar + " (cap " + maxEstSamples + ")");
        return nStar;
    }

    // Smallest N whose predicted clipped epsilon <= a SAFETY-margin'd target. Two reasons
    // the raw point prediction undersizes: the pilot (meanU/varU) is noisy, and aiming for
    // predEps == target leaves zero headroom, so the realized eps lands just over target.
    // We therefore size to sizeTarget = targetEpsilon * sizeSafety (< target). If even maxEst
    // cannot reach the REAL target, sampling will never converge this pfunc (heavy tail), so we
    // draw only a modest unreachableCap and let the run report honestly loose / fall back to q_m,
    // instead of wastefully drawing the full maxEst.
    // predictClipEps is monotone decreasing in N (Delta and bias both shrink), so binary search is valid.
    private int solveClipSizing(double mean, double var, double meanU, double varU,
                                double uRange, double pHat, double deltaPer) {
        double sizeSafety = getConfigDouble("branchmarkstar.pac.sizeSafety", 0.9);
        double sizeTarget = targetEpsilon * sizeSafety;
        int unreachableCap = getConfigInteger("branchmarkstar.pac.unreachableCap",
                Math.min(maxEstSamples, Math.max(2 * pilotSamples, 400)));
        double epsAtCap = predictClipEps(maxEstSamples, mean, var, meanU, varU, uRange, pHat, deltaPer);
        if (epsAtCap > targetEpsilon) {
            return Math.min(maxEstSamples, unreachableCap); // unreachable by sampling -> don't waste the full budget
        }
        if (epsAtCap > sizeTarget) {
            return maxEstSamples; // reachable at target but the safety margin needs the whole budget
        }
        int lo = 2, hi = maxEstSamples;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (predictClipEps(mid, mean, var, meanU, varU, uRange, pHat, deltaPer) <= sizeTarget) hi = mid;
            else lo = mid + 1;
        }
        return lo;
    }

    // Sampling-only clipped epsilon at sample size N (ignores the q_m cap, which only tightens):
    //   eps = 1 - max(0,mean-D) / (mean + D + bias),
    //   D = Bernstein(N,var,range=1), bias = min(uRange*(pHat+Hoeffding), meanU + Bernstein(N,varU,uRange)).
    private double predictClipEps(int N, double mean, double var, double meanU, double varU,
                                  double uRange, double pHat, double deltaPer) {
        double d = solveBernsteinDelta(N, var, 1.0, deltaPer);
        double bias = 0.0;
        if (uRange > 0) {
            double biasWorst = uRange * Math.min(1.0, pHat + hoeffdingHalfWidth(N, deltaPer));
            double biasEmp = meanU + solveBernsteinDelta(N, varU, uRange, deltaPer);
            bias = Math.min(biasWorst, biasEmp);
        }
        double lower = Math.max(0.0, mean - d);
        double upper = mean + d + bias;
        if (!(upper > 0)) return 1.0;
        return 1.0 - lower / upper;
    }

    private double resolveResidualRangeForBernstein(double shift, double sampleRange, String context) {
        double rangeBound = shiftedResidualWeightRangeBound(shift);
        if (Double.isFinite(rangeBound) && rangeBound >= sampleRange) {
            return rangeBound;
        } else if (Double.isFinite(rangeBound)) {
            System.out.println("[PAC-2stage] WARNING: branchmarkstar.pac.residualBound is smaller than "
                    + "the observed residual-weight range; " + context
                    + " uses the sample range and is empirical, not a strict PAC certificate.");
            return sampleRange;
        }
        System.out.println("[PAC-2stage] WARNING: no finite branchmarkstar.pac.residualBound is configured; "
                + context + " uses the sample range and is empirical, not a strict PAC certificate.");
        return sampleRange;
    }

    // meanW = E_{p_eta}[exp(-xi/RT)] estimated on a sample, computed via the
    // log-sum-exp shift to avoid overflow. Used as a collapse detector: meanW near 0
    // means q_eta badly overestimates q (eta trained off the p_eta distribution).
    private double estimateMeanW(List<CCDResult> samples, EnergyMatrix correctedEmat) {
        int n = samples.size();
        if (n == 0) return 0.0;
        double shift = Double.NEGATIVE_INFINITY;
        double[] logW = new double[n];
        for (int i = 0; i < n; i++) {
            CCDResult r = samples.get(i);
            double xi = r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat);
            logW[i] = -xi / RT;
            shift = Math.max(shift, logW[i]);
        }
        if (!Double.isFinite(shift)) return 0.0;
        double sum = 0;
        for (int i = 0; i < n; i++) sum += Math.exp(logW[i] - shift);
        double logMeanW = shift + Math.log(sum / n);
        return logMeanW < 50.0 ? Math.exp(logMeanW) : Double.MAX_VALUE;
    }

    // Fraction of pilot samples that touch an under-trained one-body cell: a
    // (pos,rc) the training distribution p_m sampled fewer than minCount times.
    // This is the root-cause signal for distribution-shift collapse -- p_eta has
    // concentrated mass on rotamers eta never learned -- and unlike absolute
    // meanW it is robust to rare huge-weight (negative-xi) outliers. Re-learning
    // eta on train+pilot covers exactly these cells, so a high fraction also
    // predicts the refinement will help.
    private double pilotUndertrainedFraction(List<CCDResult> pilot, EtaCorrections eta, int minCount) {
        if (pilot.isEmpty()) return 0.0;
        int hit = 0;
        for (CCDResult r : pilot) {
            for (int pos = 0; pos < r.conf.length; pos++) {
                int rc = r.conf[pos];
                if (rc < 0) continue;
                if (eta.oneBodyCounts[pos][rc] < minCount) { hit++; break; }
            }
        }
        return (double) hit / pilot.size();
    }

    private double shiftedResidualWeightRangeBound(double shift) {
        if (!Double.isFinite(residualBoundKcal)) return Double.NaN;
        double x = residualBoundKcal / RT;
        double hi = Math.exp(x - shift);
        double lo = Math.exp(-x - shift);
        double range = hi - lo;
        return range >= 0.0 ? range : Double.NaN;
    }

    // ===== Diagnostic dump for multi-body residual attribution =====
    // Enabled by -Dbranchmarkstar.pac.dumpDir=<dir>. Writes, per pfunc, the
    // per-sample residual xi = E_true - E_eta (pairwise part already removed) plus
    // the branch-decomposition bags, so offline analysis can locate the 3-body
    // interactions that the pairwise eta cannot absorb.
    private static final AtomicInteger PAC_DUMP_COUNTER = new AtomicInteger(0);

    private void maybeDumpResidualSamples(List<CCDResult> est, double[] xi,
                                          double logZCorrected,
                                          double meanXi, double stdXi) {
        String dir = getConfigProperty("branchmarkstar.pac.dumpDir", null);
        if (dir == null || dir.trim().isEmpty()) return;
        try {
            File d = new File(dir.trim());
            d.mkdirs();
            int n = PAC_DUMP_COUNTER.incrementAndGet();
            int numPos = rcs.getNumPos();
            File rf = new File(d, String.format("residual_%04d.tsv", n));
            try (PrintWriter w = new PrintWriter(rf)) {
                w.printf("# logZcorr=%.6f N=%d numPos=%d meanXi=%.6f stdXi=%.6f%n",
                        logZCorrected, est.size(), numPos, meanXi, stdXi);
                w.println("# idx\tconf(comma)\teTrue\teEta\txi");
                for (int i = 0; i < est.size(); i++) {
                    CCDResult r = est.get(i);
                    StringBuilder c = new StringBuilder();
                    for (int p = 0; p < r.conf.length; p++) {
                        if (p > 0) c.append(',');
                        c.append(r.conf[p]);
                    }
                    double eEta = r.eTrue - xi[i];
                    w.printf("%d\t%s\t%.6f\t%.6f\t%.6f%n", i, c, r.eTrue, eEta, xi[i]);
                }
            }
            File bf = new File(d, String.format("bags_%04d.tsv", n));
            try (PrintWriter w = new PrintWriter(bf)) {
                w.println("# edgeIdx\tlambda(comma)\tM(comma)");
                List<RootedTreeEdge> order = buildTopDownOrder(rootedRootEdge);
                int e = 0;
                for (RootedTreeEdge edge : order) {
                    w.printf("%d\t%s\t%s%n", e++,
                            joinInts(edge.getLambda()), joinInts(edge.getM()));
                }
            }
            System.out.println("[PAC-2stage] residual dump -> " + rf.getName()
                    + " (N=" + est.size() + ", logZcorr=" + String.format("%.4f", logZCorrected)
                    + ", meanXi=" + String.format("%.4f", meanXi) + ")");
        } catch (Exception ex) {
            System.err.println("[PAC-2stage] residual dump failed: " + ex.getMessage());
        }
    }

    private static final AtomicInteger PAC_TRAIN_COUNTER = new AtomicInteger(0);

    private void maybeDumpTrainingSamples(List<CCDResult> train, EtaCorrections eta,
                                          EnergyMatrix correctedEmat) {
        String dir = getConfigProperty("branchmarkstar.pac.dumpDir", null);
        if (dir == null || dir.trim().isEmpty()) return;
        try {
            File d = new File(dir.trim());
            d.mkdirs();
            int n = PAC_TRAIN_COUNTER.incrementAndGet();
            int numPos = rcs.getNumPos();
            File tf = new File(d, String.format("train_%04d.tsv", n));
            try (PrintWriter w = new PrintWriter(tf)) {
                w.printf("# trainN=%d numPos=%d (p_m samples)%n", train.size(), numPos);
                w.println("# idx\tconf(comma)\teTrue\teMin\teEta\txiTrain");
                for (int i = 0; i < train.size(); i++) {
                    CCDResult r = train.get(i);
                    StringBuilder c = new StringBuilder();
                    for (int p = 0; p < r.conf.length; p++) {
                        if (p > 0) c.append(',');
                        c.append(r.conf[p]);
                    }
                    double eEta = computeFullConfPairwiseEnergy(r.conf, correctedEmat);
                    w.printf("%d\t%s\t%.6f\t%.6f\t%.6f\t%.6f%n",
                            i, c, r.eTrue, r.eMin, eEta, r.eTrue - eEta);
                }
            }
            File ef = new File(d, String.format("eta_%04d.tsv", n));
            try (PrintWriter w = new PrintWriter(ef)) {
                w.println("# pos\trc\tetaOneBody\ttrainCount");
                for (int p = 0; p < eta.oneBody.length; p++) {
                    for (int rc = 0; rc < eta.oneBody[p].length; rc++) {
                        int cnt = eta.oneBodyCounts[p][rc];
                        if (cnt > 0) {
                            w.printf("%d\t%d\t%.6f\t%d%n", p, rc, eta.oneBody[p][rc] / cnt, cnt);
                        }
                    }
                }
            }
            System.out.println("[PAC-2stage] training dump -> " + tf.getName()
                    + " (trainN=" + train.size() + ")");
        } catch (Exception ex) {
            System.err.println("[PAC-2stage] training dump failed: " + ex.getMessage());
        }
    }

    private static String joinInts(java.util.Collection<Integer> s) {
        StringBuilder b = new StringBuilder();
        boolean first = true;
        for (Integer v : s) {
            if (!first) b.append(',');
            b.append(v);
            first = false;
        }
        return b.toString();
    }

    /**
     * Final residual-leg bound: estimand w = exp(-xi/RT), xi = E_true - E_eta, with
     * samples drawn from p_eta. Z = q_eta * E_{p_eta}[w], so
     *   logZ = logZCorrected + log(meanW +/- Delta),  Delta = empirical Bernstein.
     * The deterministic factor q_eta needs no bound; only the residual leg does.
     */
    private void computePACBoundResidual(List<CCDResult> estResults,
                                         EnergyMatrix correctedEmat,
                                         double logZCorrected) {
        int N = estResults.size();
        if (N == 0) {
            setZeroBounds("two-stage: no estimation samples");
            return;
        }
        if (getConfigBoolean("branchmarkstar.pac.clip", false)) {
            double clipQuantile = getConfigDouble("branchmarkstar.pac.clipQuantile", 0.9);
            computePACBoundResidualClipped(estResults, correctedEmat, logZCorrected, clipQuantile);
            return;
        }
        double shift = Double.NEGATIVE_INFINITY;
        double[] logW = new double[N];
        double[] xi = new double[N];
        for (int i = 0; i < N; i++) {
            CCDResult r = estResults.get(i);
            xi[i] = r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat);
            logW[i] = -xi[i] / RT;
            shift = Math.max(shift, logW[i]);
        }
        if (!Double.isFinite(shift)) {
            setZeroBounds("two-stage: residual weights are non-finite");
            return;
        }
        double sum = 0, sum2 = 0, mn = Double.MAX_VALUE, mx = -Double.MAX_VALUE;
        for (int i = 0; i < N; i++) {
            double w = Math.exp(logW[i] - shift);
            sum += w; sum2 += w * w;
            mn = Math.min(mn, w); mx = Math.max(mx, w);
        }
        double mean = sum / N;
        double var = N > 1 ? Math.max(0.0, (sum2 - N * mean * mean) / (N - 1)) : 0.0;
        double sampleRange = mx - mn;
        double range = resolveResidualRangeForBernstein(shift, sampleRange, "final bound");

        double sumRes = 0, sumRes2 = 0;
        for (int i = 0; i < N; i++) { sumRes += xi[i]; sumRes2 += xi[i] * xi[i]; }
        meanResidual = sumRes / N;
        stdResidual = N > 1
                ? Math.sqrt(Math.max(0.0, (sumRes2 - N * meanResidual * meanResidual) / (N - 1)))
                : 0.0;
        meanPsi = shiftedMeanToDouble(mean, shift);
        varPsi = var;
        cvPsi = (mean > 0) ? Math.sqrt(var) / mean : Double.MAX_VALUE;

        maybeDumpResidualSamples(estResults, xi, logZCorrected, meanResidual, stdResidual);

        double deltaPer = delta / 6.0; // 2 tails x 3 pfuncs
        double boundDelta = solveBernsteinDelta(N, var, range, deltaPer);
        if (!Double.isFinite(mean) || !Double.isFinite(var)
                || !Double.isFinite(range) || !Double.isFinite(boundDelta)) {
            setZeroBounds("two-stage: PAC statistics are non-finite");
            return;
        }
        double meanLower = Math.max(0.0, mean - boundDelta);
        double meanUpper = mean + boundDelta;
        if (!Double.isFinite(meanUpper) || meanUpper <= 0) {
            setZeroBounds("two-stage: PAC confidence interval is non-finite");
            return;
        }

        logZLowerPAC = meanLower > 0
                ? logZCorrected + shift + Math.log(meanLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = logZCorrected + shift + Math.log(meanUpper);
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);

        System.out.println("[PAC-2stage] Bernstein(residual): Delta=" + String.format("%.6f", boundDelta)
                + ", meanW=" + String.format("%.6f", mean)
                + ", rangeW=" + String.format("%.6f", range)
                + ", sampleRangeW=" + String.format("%.6f", sampleRange)
                + ", varW=" + String.format("%.6f", var)
                + ", shift=" + String.format("%.4f", shift)
                + ", logZcorr=" + String.format("%.4f", logZCorrected)
                + ", meanResidual=" + String.format("%.4f", meanResidual) + " kcal/mol"
                + ", stdResidual=" + String.format("%.4f", stdResidual) + " kcal/mol"
                + ", log10ZLower=" + formatLog10(logZLowerPAC)
                + ", log10ZUpper=" + formatLog10(logZUpperPAC)
                + ", epsilon=" + String.format("%.6f", epsilon)
                + ", N=" + N);
    }

    /**
     * Clipped residual-leg bound (step 4: drop the global residualBound assumption).
     *
     * The unbounded estimand is mu_w = E_{p_eta}[w], w = exp(-xi/RT). A rare
     * over-corrected conformation (xi << 0) makes w explode, so the empirical
     * Bernstein "range" term -- and any deterministic |xi|<=xi_max bound -- blows
     * up exponentially. Instead we clip at a per-pfunc threshold c = exp(logCap),
     * where logCap is a data-driven quantile of logW, and estimate E[min(w,c)].
     *
     * Working in units of c (wbar = min(w,c)/c in [0,1]):
     *  - Lower bound on mu_w: E[min(w,c)] <= E[w], and the (1-delta') Bernstein
     *    lower CI on E[min(w,c)] is therefore a rigorous lower bound on mu_w,
     *    WITH NO bounded-residual assumption. The range term is the deterministic
     *    constant 1 -- no exp explosion.
     *  - Upper bound: mu_w = E[min(w,c)] + E[(w-c)_+]. The clip bias E[(w-c)_+] is
     *    bounded one-sidedly: if xi >= -residualBound (one-sided), then w <= Wmax =
     *    exp(residualBound/RT), so (w-c)_+ <= (Wmax-c) on clipped samples and
     *    E[(w-c)_+] <= (Wmax/c - 1) * pUpper, with pUpper a (1-delta') Hoeffding
     *    upper CI on the clip probability P(w>c). This is additive (not exp-
     *    amplified through the whole bound) and only needs a far weaker one-sided
     *    tail assumption than the original two-sided |xi|<=xi_max.
     */
    private void computePACBoundResidualClipped(List<CCDResult> estResults,
                                                EnergyMatrix correctedEmat,
                                                double logZCorrected,
                                                double clipQuantile) {
        int N = estResults.size();
        if (N == 0) { setZeroBounds("two-stage(clip): no estimation samples"); return; }
        double[] logW = new double[N];
        double[] xi = new double[N];
        for (int i = 0; i < N; i++) {
            CCDResult r = estResults.get(i);
            xi[i] = r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat);
            logW[i] = -xi[i] / RT;
        }
        double sumRes = 0, sumRes2 = 0;
        for (int i = 0; i < N; i++) { sumRes += xi[i]; sumRes2 += xi[i] * xi[i]; }
        meanResidual = sumRes / N;
        stdResidual = N > 1
                ? Math.sqrt(Math.max(0.0, (sumRes2 - N * meanResidual * meanResidual) / (N - 1)))
                : 0.0;

        // per-pfunc clip threshold: a data-driven quantile of logW fixed on the
        // held-out pilot (PAC: chosen before S2 -> clipped weights i.i.d. given c).
        // Fall back to the est-set quantile only if the pilot cap is unavailable.
        double logCap = Double.isFinite(clipLogCap) ? clipLogCap : quantile(logW, clipQuantile);
        if (!Double.isFinite(logCap)) {
            setZeroBounds("two-stage(clip): non-finite clip threshold");
            return;
        }

        // clipped weights in units of c = exp(logCap): wbar = exp(min(logW,logCap)-logCap) in (0,1]
        double sum = 0, sum2 = 0;
        int nClipped = 0;
        for (int i = 0; i < N; i++) {
            if (logW[i] > logCap) nClipped++;
            double w = Math.exp(Math.min(logW[i], logCap) - logCap);
            sum += w; sum2 += w * w;
        }
        double mean = sum / N;
        double var = N > 1 ? Math.max(0.0, (sum2 - N * mean * mean) / (N - 1)) : 0.0;
        double range = 1.0; // deterministic: clipped weights lie in [0,1] in units of c

        meanPsi = shiftedMeanToDouble(mean, logCap);
        varPsi = var;
        cvPsi = (mean > 0) ? Math.sqrt(var) / mean : Double.MAX_VALUE;

        maybeDumpResidualSamples(estResults, xi, logZCorrected, meanResidual, stdResidual);

        double deltaPer = delta / 6.0; // 2 tails x 3 pfuncs
        double boundDelta = solveBernsteinDelta(N, var, range, deltaPer);
        if (!Double.isFinite(mean) || !Double.isFinite(var) || !Double.isFinite(boundDelta)) {
            setZeroBounds("two-stage(clip): PAC statistics are non-finite");
            return;
        }
        double meanLower = Math.max(0.0, mean - boundDelta);

        // ---- upper side: clipped mass (Bernstein) + clip-bias E[(w-c)+]/c ----
        double meanUpper = mean + boundDelta;
        // Empirical clip-excess in units of c: u_i = max(0, w_i/c - 1), nonzero only on clipped samples.
        double sumU = 0, sumU2 = 0;
        for (int i = 0; i < N; i++) {
            double u = Math.max(0.0, Math.exp(logW[i] - logCap) - 1.0);
            sumU += u; sumU2 += u * u;
        }
        double meanU = sumU / N;
        double varU = N > 1 ? Math.max(0.0, (sumU2 - N * meanU * meanU) / (N - 1)) : 0.0;

        double biasUpper;
        boolean biasCertified;
        String biasMode;
        if (Double.isFinite(residualBoundKcal) && residualBoundKcal / RT - logCap > 0) {
            // residualBound caps the (possibly unobserved) tail: u in [0, uRange], uRange = Wmax/c - 1.
            double uRange = Math.expm1(residualBoundKcal / RT - logCap);
            double pHat = (double) nClipped / N;
            double pUpper = Math.min(1.0, pHat + hoeffdingHalfWidth(N, deltaPer));
            double biasWorst = uRange * pUpper;                                      // every clipped sample at Wmax (old, loose)
            double biasEmp = meanU + solveBernsteinDelta(N, varU, uRange, deltaPer); // observed excess + concentration (tight)
            biasUpper = Math.min(biasWorst, biasEmp);                                // both valid upper bounds -> take tighter
            biasCertified = true;
            biasMode = (biasEmp <= biasWorst) ? "empirical" : "worstcase";
        } else if (Double.isFinite(residualBoundKcal)) {
            biasUpper = 0.0; biasCertified = true; biasMode = "capAboveWmax"; // logCap already >= Wmax: no excess
        } else {
            biasUpper = Double.POSITIVE_INFINITY; biasCertified = false; biasMode = "none"; // no tail bound -> rely on q_m
        }
        double meanUpperTotal = meanUpper + biasUpper;

        // lower bound (assumption-free: clipped weights in [0,1] by construction)
        logZLowerPAC = meanLower > 0
                ? logZCorrected + logCap + Math.log(meanLower)
                : Double.NEGATIVE_INFINITY;
        // upper bound: sampling clip bound, then cap by the deterministic q_m (E_m <= E_true => q_m >= q, zero assumptions)
        double logZUpperSampling = (Double.isFinite(meanUpperTotal) && meanUpperTotal > 0)
                ? logZCorrected + logCap + Math.log(meanUpperTotal)
                : Double.POSITIVE_INFINITY;
        boolean qmBinds = logZMinDet < logZUpperSampling;
        logZUpperPAC = Math.min(logZUpperSampling, logZMinDet);
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);

        System.out.println("[PAC-2stage] Bernstein(clipped): Delta=" + String.format("%.6f", boundDelta)
                + ", clipQ=" + String.format("%.3f", clipQuantile)
                + ", logCap=" + String.format("%.4f", logCap)
                + ", nClipped=" + nClipped + "/" + N
                + ", meanWc=" + String.format("%.6f", mean)
                + ", varWc=" + String.format("%.6f", var)
                + ", range=1.0(det)"
                + ", biasUpper=" + String.format("%.6f", biasUpper) + "(" + biasMode + ")"
                + (biasCertified ? "" : "[no residualBound -> upper via q_m only]")
                + ", upperVia=" + (qmBinds ? "q_m(det)" : "sampling")
                + ", log10ZmDet=" + formatLog10(logZMinDet)
                + ", meanResidual=" + String.format("%.4f", meanResidual) + " kcal/mol"
                + ", stdResidual=" + String.format("%.4f", stdResidual) + " kcal/mol"
                + ", log10ZLower=" + formatLog10(logZLowerPAC)
                + ", log10ZUpper=" + formatLog10(logZUpperPAC)
                + ", epsilon=" + String.format("%.6f", epsilon)
                + ", N=" + N);
    }

    // linear-interpolated quantile of a value array (does not mutate input)
    private static double quantile(double[] vals, double q) {
        int n = vals.length;
        if (n == 0) return Double.NaN;
        if (n == 1) return vals[0];
        double[] s = vals.clone();
        Arrays.sort(s);
        double pos = Math.max(0.0, Math.min(1.0, q)) * (n - 1);
        int lo = (int) Math.floor(pos);
        int hi = (int) Math.ceil(pos);
        if (lo == hi) return s[lo];
        double frac = pos - lo;
        return s[lo] * (1.0 - frac) + s[hi] * frac;
    }

    // Hoeffding upper-CI half-width for a [0,1]-bounded mean (the clip indicator):
    // P(p > pHat + sqrt(ln(1/delta)/(2N))) <= delta.
    private static double hoeffdingHalfWidth(int n, double delta) {
        if (n <= 0 || !(delta > 0)) return 1.0;
        return Math.sqrt(Math.log(1.0 / delta) / (2.0 * n));
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
            int edgeOrdinal = 0;
            for (RootedTreeEdge edge : topDownOrder) {
                if (!edge.getIsLambdaEdge()) continue;

                long edgeStart = System.currentTimeMillis();
                Map<Long, SampleGroup> groups = groupSamplesByMIdx(edge, confs);
                long edgeSeed = mix64(masterSeed + 0xD1B54A32D192ED03L * (long) edgeOrdinal++);
                sampleEdgeGroups(edge, confs, sampleRngs, groups, samplingPool, edgeSeed);
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
                                  ExecutorService samplingPool,
                                  long edgeSeed) {
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

        if (trySampleEdgeGroupsGpu(edge, groups, confs, edgeSeed)) {
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

    private boolean trySampleEdgeGroupsGpu(RootedTreeEdge edge,
                                           Map<Long, SampleGroup> groups,
                                           int[][] confs,
                                           long edgeSeed) {
        if (!gpuSampling || !edge.canUseGpuSampling()) {
            return false;
        }

        long[] mIdxPerSample = new long[confs.length];
        for (SampleGroup group : groups.values()) {
            for (int i = 0; i < group.count; i++) {
                mIdxPerSample[group.sampleIndices[i]] = group.mIdx;
            }
        }

        int[] lIdxPerSample = edge.sampleLambdaStatesGpu(mIdxPerSample, edgeSeed, samplingProgress);
        if (lIdxPerSample == null || lIdxPerSample.length != confs.length) {
            return false;
        }

        for (int sampleIndex = 0; sampleIndex < confs.length; sampleIndex++) {
            writeLambdaState(edge, lIdxPerSample[sampleIndex], confs[sampleIndex]);
        }
        return true;
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
        // The corrected-DP cache is structurally unhittable: its key embeds the raw
        // double bits of the per-edge eta corrections (appendLocalEtaSignature), and
        // eta is re-learned from random CCD samples on every pfunc / sequence /
        // sample-split phase, so the floats never repeat bit-for-bit. Caching here is
        // therefore pure overhead (key construction, hashing, storing tables that are
        // never reused). Always take the direct path.
        RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
        return null;
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

    /**
     * PAC epsilon from log-space Z bounds (design 1.3). Computes
     *   epsilon = 1 - Z_lower/Z_upper = 1 - exp(logZLower - logZUpper)
     * from the log ratio, so it stays finite when |logZ| ~ 1e6 and exp(logZ)
     * itself would underflow/overflow a double (the old direct-double path wrote
     * Z=0 -> fake epsilon=NaN/1, e.g. 5dc4). Tiny FP excursions just outside
     * [0,1] are snapped; non-finite bounds mean no usable bound -> epsilon = 1.
     * Unit-tested in TestPACLogSpaceBounds.
     */
    static double epsilonFromLogBounds(double logZLower, double logZUpper) {
        if (Double.isFinite(logZLower) && Double.isFinite(logZUpper)) {
            double eps = 1.0 - Math.exp(logZLower - logZUpper);
            if (eps < 0.0 && eps > -1e-12) eps = 0.0;
            if (eps > 1.0 && eps < 1.0 + 1e-12) eps = 1.0;
            return eps;
        }
        return 1.0;
    }

    /**
     * Empirical Bernstein bound (Maurer & Pontil, 2009) for the one-sided
     * deviation of the sample mean from E[X], at confidence deltaTarget:
     *
     *   E[X] - mean(X) <= sqrt(2 * s2 * ln(2/delta) / N)
     *                     + 7 * range * ln(2/delta) / (3 * (N - 1))
     *
     * This is the *fully empirical* form. Unlike plugging an empirical variance
     * into the variance-oracle Bernstein exponent, the second (range) term is
     * the sample-variance penalty that keeps the bound valid when s2 is itself
     * estimated from the same N samples. By symmetry the same Delta bounds the
     * opposite tail (-X has identical empirical variance and range); the caller
     * passes deltaTarget = delta/6 to cover both tails of all three partition
     * functions (2 tails x 3 pfuncs = 6 events).
     */
    private double solveBernsteinDelta(int N, double s2, double range, double deltaTarget) {
        if (N <= 1) return Double.POSITIVE_INFINITY;
        double logTerm = Math.log(2.0 / deltaTarget);
        double varTerm = Math.sqrt(Math.max(0.0, 2.0 * s2 * logTerm / N));
        double rangeTerm = 7.0 * range * logTerm / (3.0 * (N - 1));
        return varTerm + rangeTerm;
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
