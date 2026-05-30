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
import java.util.concurrent.atomic.AtomicInteger;

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
    private static final String DP_CACHE_ENABLED_PROPERTY = "branchmarkstar.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "branchmarkstar.dp.cache.maxEntries";
    private static final int DEFAULT_SAMPLES = 500;
    private static final double DEFAULT_CONFIDENCE = 0.05; // delta = 0.05 => 95% confidence
    private static final int DEFAULT_DP_CACHE_MAX_ENTRIES = 20000;

    private static final Object CORRECTED_DP_CACHE_LOCK = new Object();
    private static final LinkedHashMap<String, CachedDPTable> CORRECTED_DP_CACHE =
            new LinkedHashMap<>(1024, 0.75f, true);

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
    private final boolean dpCacheEnabled;
    private final int dpCacheMaxEntries;

    // Results
    private BigDecimal zLower;
    private BigDecimal zUpper;
    private double epsilon;
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

    public PACPartitionFunction(RootedTreeNode rootedRoot,
                                 RootedTreeEdge rootedRootEdge,
                                 EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace) {
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
        this.dpCacheEnabled = getConfigBoolean(DP_CACHE_ENABLED_PROPERTY, true);
        this.dpCacheMaxEntries = Math.max(0,
                getConfigInteger(DP_CACHE_MAX_ENTRIES_PROPERTY, DEFAULT_DP_CACHE_MAX_ENTRIES));
    }

    private static class CachedDPTable {
        final double[] lower;
        final double[] upper;

        CachedDPTable(double[] lower, double[] upper) {
            this.lower = lower;
            this.upper = upper;
        }
    }

    private static class DPCacheStats {
        int hits = 0;
        int misses = 0;
        int stores = 0;
    }

    private static class CorrectedDPResult {
        final double logZCorrected;
        final DPCacheStats cacheStats;

        CorrectedDPResult(double logZCorrected, DPCacheStats cacheStats) {
            this.logZCorrected = logZCorrected;
            this.cacheStats = cacheStats;
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

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        String value = getConfigProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    /**
     * Run the full PAC estimation pipeline.
     * Returns the epsilon achieved.
     */
    public double compute() {
        long startTime = System.currentTimeMillis();

        // Phase 0: Z_min from existing DP (already computed)
        double logZMin = rootedRootEdge.getLogZUpper()[0];
        double logZRigid = rootedRootEdge.getLogZLower()[0];
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

        // Phase 1: Sample conformations from p(c) ∝ exp(-E_min/kT)
        long phase1Start = System.currentTimeMillis();
        List<int[]> sampledConfs = sampleConformationsFromDP(numSamples);
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

        // Phase 3: Extract per-term eta corrections
        long phase3Start = System.currentTimeMillis();
        EtaCorrections eta = extractEtaCorrections(ccdResults);
        long phase3Time = System.currentTimeMillis() - phase3Start;
        System.out.println("[PAC] Phase 3: Eta extraction in " + phase3Time + " ms"
                + ", oneBodyTerms=" + eta.oneBodyCount
                + ", pairTerms=" + eta.pairCount);

        // Phase 4: Build corrected emat and recompute DP
        long phase4Start = System.currentTimeMillis();
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta);
        double logZCorrected = correctedDP.logZCorrected;
        long phase4Time = System.currentTimeMillis() - phase4Start;
        System.out.println("[PAC] Phase 4: logZ_corrected=" + String.format("%.4f", logZCorrected)
                + " (improvement over logZ_rigid: " + String.format("%.4f", logZCorrected - logZRigid) + ")"
                + " in " + phase4Time + " ms"
                + formatDPCacheStats(correctedDP.cacheStats));

        // Phase 5: Compute residuals and PAC bound
        long phase5Start = System.currentTimeMillis();
        computePACBound(ccdResults, correctedEmat, logZCorrected);
        long phase5Time = System.currentTimeMillis() - phase5Start;

        long totalTime = System.currentTimeMillis() - startTime;
        System.out.println("[PAC] Phase 5: epsilon=" + String.format("%.6f", epsilon)
                + ", meanPsi=" + String.format("%.6f", meanPsi)
                + ", cvPsi=" + String.format("%.4f", cvPsi)
                + ", meanResidual=" + String.format("%.4f", meanResidual) + " kcal/mol"
                + ", stdResidual=" + String.format("%.4f", stdResidual) + " kcal/mol"
                + " in " + phase5Time + " ms");
        System.out.println("[PAC] Total: " + totalTime + " ms, " + totalCCDCalls
                + " CCD calls, epsilon=" + String.format("%.6f", epsilon)
                + ", confidence=" + String.format("%.2f%%", (1.0 - delta) * 100));
        System.out.println("[PAC] Z bounds: lower=" + String.format("%.6e", zLower.doubleValue())
                + ", upper=" + String.format("%.6e", zUpper.doubleValue()));

        return epsilon;
    }

    // ========== Phase 1: Ancestral sampling from DP ==========

    /**
     * Sample conformations from p(c) ∝ exp(-E_min(c)/kT) using ancestral sampling
     * on the tree decomposition.
     *
     * Top-down traversal: for each edge, sample lambda-state proportional to
     * exp(-fullEnergyMin[mIdx][lIdx]/RT), conditioned on the M-state from parent.
     */
    private List<int[]> sampleConformationsFromDP(int n) {
        Random rng = new Random(42); // reproducible
        int numPos = rcs.getNumPos();
        List<int[]> conformations = new ArrayList<>(n);

        // Collect all lambda-edges in compact tree order (post-order gives leaves first)
        List<RootedTreeEdge> lambdaEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, lambdaEdges);

        // Build a top-down traversal order: root edge first, then children
        List<RootedTreeEdge> topDownOrder = buildTopDownOrder(rootedRootEdge);

        for (int s = 0; s < n; s++) {
            int[] conf = new int[numPos];
            Arrays.fill(conf, -1);

            // Top-down: sample each edge's lambda-state given its M-state (from parent)
            for (RootedTreeEdge edge : topDownOrder) {
                if (!edge.getIsLambdaEdge()) continue;

                int[] mPositions = edge.getMPositionsSorted();
                int[] lambdaPositions = edge.getLambdaPositionsSorted();

                // Determine M-state from already-sampled positions
                int[] mLocalRCs = new int[mPositions.length];
                for (int i = 0; i < mPositions.length; i++) {
                    int globalRC = conf[mPositions[i]];
                    // Convert global RC back to local index
                    mLocalRCs[i] = findLocalRCIndex(mPositions[i], globalRC);
                }
                int mIdx = edge.computeIndexInA(mLocalRCs);

                // Sample lambda-state proportional to exp(-fullEnergyMin/RT)
                int lIdx = sampleLambdaState(edge, mIdx, rng);

                // Decode lambda-state and write to conf
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

    /**
     * Find the local RC index for a position given its global RC.
     */
    private int findLocalRCIndex(int pos, int globalRC) {
        for (int localIdx = 0; localIdx < rcs.getNum(pos); localIdx++) {
            if (rcs.get(pos, localIdx) == globalRC) {
                return localIdx;
            }
        }
        throw new IllegalStateException("Global RC " + globalRC + " not found at position " + pos);
    }

    /**
     * Sample a lambda-state proportional to its Boltzmann weight exp(-fullEnergyMin/RT).
     * For leaf edges, uses pre-computed fullEnergyMin.
     * For non-leaf edges, uses local energy + child logZ.
     */
    private int sampleLambdaState(RootedTreeEdge edge, int mIdx, Random rng) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] logWeights = new double[totalLambda];

        double[][] fullEnergyMin = edge.getFullEnergyMin();

        if (!edge.hasFsetChildren()) {
            // Leaf edge: weight = exp(-fullEnergyMin[mIdx][lIdx]/RT)
            for (int lIdx = 0; lIdx < totalLambda; lIdx++) {
                logWeights[lIdx] = -fullEnergyMin[mIdx][lIdx] / RT;
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
                        int fIdx = fEdge.computeIndexInA(fM);
                        fSumUpper += fEdge.getLogZUpper()[fIdx];
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

        return new CorrectedDPResult(rootedRootEdge.getLogZUpper()[0], stats);
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

        String key = buildCorrectedDPCacheKey(namespace, edge, eta, edgeKeys);
        CachedDPTable cached = getCorrectedCachedDPTable(key, edge.getMArraySize());
        if (cached != null) {
            System.arraycopy(cached.lower, 0, edge.getLogZLower(), 0, cached.lower.length);
            System.arraycopy(cached.upper, 0, edge.getLogZUpper(), 0, cached.upper.length);
            stats.hits++;
        } else {
            edge.computeFullDP();
            putCorrectedCachedDPTable(key, edge.getLogZLower(), edge.getLogZUpper());
            stats.misses++;
            stats.stores++;
        }
        edgeKeys.put(edge, key);
    }

    private CachedDPTable getCorrectedCachedDPTable(String key, int expectedLength) {
        synchronized (CORRECTED_DP_CACHE_LOCK) {
            CachedDPTable cached = CORRECTED_DP_CACHE.get(key);
            if (cached == null) return null;
            if (cached.lower.length != expectedLength || cached.upper.length != expectedLength) {
                CORRECTED_DP_CACHE.remove(key);
                return null;
            }
            return cached;
        }
    }

    private void putCorrectedCachedDPTable(String key, double[] lower, double[] upper) {
        synchronized (CORRECTED_DP_CACHE_LOCK) {
            CORRECTED_DP_CACHE.put(key, new CachedDPTable(
                    Arrays.copyOf(lower, lower.length),
                    Arrays.copyOf(upper, upper.length)));
            while (CORRECTED_DP_CACHE.size() > dpCacheMaxEntries) {
                Iterator<String> it = CORRECTED_DP_CACHE.keySet().iterator();
                if (!it.hasNext()) break;
                it.next();
                it.remove();
            }
        }
    }

    private String formatDPCacheStats(DPCacheStats stats) {
        if (stats == null) return "";
        int size;
        synchronized (CORRECTED_DP_CACHE_LOCK) {
            size = CORRECTED_DP_CACHE.size();
        }
        return ", correctedDPCacheHits=" + stats.hits
                + ", correctedDPCacheMisses=" + stats.misses
                + ", correctedDPCacheSize=" + size;
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
        double[] phiValues = new double[N];    // phi(c) = exp(-g/kT), the basic IS weight
        double[] alphaValues = new double[N];  // alpha(c) = exp(-eta/kT), the control variate
        double[] residuals = new double[N];    // residual = g - eta (for diagnostics)

        for (int i = 0; i < N; i++) {
            CCDResult result = ccdResults.get(i);
            double g = result.eTrue - result.eMin; // g(c) = E_true - E_min
            double eta = computeFullConfPairwiseEnergy(result.conf, correctedEmat)
                    - result.eMin;  // eta(c) = E_corrected - E_min
            residuals[i] = g - eta;
            phiValues[i] = Math.exp(-g / RT);
            alphaValues[i] = Math.exp(-eta / RT);
        }

        // E_p[alpha] = Z_corrected / Z_min (exact from DP)
        // logZMinOriginal was saved before Phase 4 modified the DP tables.
        double logRatio = logZCorrected - logZMinOriginal;
        double expectedAlpha = Math.exp(logRatio);
        if (Double.isNaN(logRatio) || Double.isNaN(expectedAlpha)) {
            setZeroBounds("corrected/min DP ratio is NaN");
            return;
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
        covPhiAlpha /= (N - 1);
        varAlpha /= (N - 1);

        double beta = (varAlpha > 1e-30) ? covPhiAlpha / varAlpha : 0.0;

        // Compute control-variate-adjusted values: phi_cv = phi - beta*(alpha - E[alpha])
        double[] phiCV = new double[N];
        double sumCV = 0, sumCV2 = 0;
        double minCV = Double.MAX_VALUE, maxCV = Double.MIN_VALUE;
        for (int i = 0; i < N; i++) {
            phiCV[i] = phiValues[i] - beta * (alphaValues[i] - expectedAlpha);
            sumCV += phiCV[i];
            sumCV2 += phiCV[i] * phiCV[i];
            minCV = Math.min(minCV, phiCV[i]);
            maxCV = Math.max(maxCV, phiCV[i]);
        }

        double meanCV = sumCV / N;
        double varCV = Math.max(0, (sumCV2 - N * meanCV * meanCV) / (N - 1));
        double rangeCV = maxCV - minCV;

        // Residual diagnostics (for reporting)
        double sumRes = 0, sumRes2 = 0;
        for (int i = 0; i < N; i++) {
            sumRes += residuals[i];
            sumRes2 += residuals[i] * residuals[i];
        }
        meanResidual = sumRes / N;
        stdResidual = Math.sqrt(Math.max(0, (sumRes2 - N * meanResidual * meanResidual) / (N - 1)));

        // Also compute basic phi stats for reporting
        meanPsi = meanPhi; // report the basic IS mean as "meanPsi" for backward compatibility
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
        // Z_lower = Z_min * (meanCV - Delta)
        // Z_upper = Z_min * (meanCV + Delta)
        BigDecimal zMin = bigExpFromLog(logZMinOriginal);

        double cvBarLower = Math.max(0, meanCV - boundDelta);
        double cvBarUpper = meanCV + boundDelta;
        if (!Double.isFinite(cvBarLower) || !Double.isFinite(cvBarUpper) || cvBarUpper < 0) {
            setZeroBounds("PAC confidence interval is non-finite");
            return;
        }

        zLower = zMin.multiply(BigDecimal.valueOf(cvBarLower), PartitionFunction.decimalPrecision);
        zUpper = zMin.multiply(BigDecimal.valueOf(cvBarUpper), PartitionFunction.decimalPrecision);

        if (zUpper.compareTo(BigDecimal.ZERO) > 0) {
            epsilon = 1.0 - zLower.doubleValue() / zUpper.doubleValue();
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
                + ", E[alpha]=" + String.format("%.6f", expectedAlpha)
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
