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

package edu.duke.cs.osprey.packstar;

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
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.branchdp.RootedTreeEdge;
import edu.duke.cs.osprey.branchdp.RootedTreeNode;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.function.IntConsumer;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;

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
 *   The frequency/severity production path selects a count/context-regularized
 *   eta by cross-fitting one q_m training batch with alpha and triple-residual
 *   gamma both fixed at one, solves one corrected DP for the winner, and draws
 *   independent pilot, monitor, and final samples.  Its tail upper bound is
 *   conditional on a frozen external severity premise.
 *   The older conditional-B and unconditional bulk-tail paths remain for
 *   archived reproducibility.
 */
public class PackStarEstimator {

    // Configuration
    private static final String PAC_SAMPLES_PROPERTY = "packstar.pac.samples";
    private static final String PAC_CONFIDENCE_PROPERTY = "packstar.pac.confidence";
    private static final String PAC_TARGET_EPSILON_PROPERTY = "packstar.pac.targetEpsilon";
    private static final String PAC_RANDOM_SEED_PROPERTY = "packstar.pac.randomSeed";
    private static final String PAC_SAMPLING_BATCHED_PROPERTY = "packstar.pac.sampling.batched";
    private static final String PAC_SAMPLING_PARALLEL_PROPERTY = "packstar.pac.sampling.parallel";
    private static final String PAC_SAMPLING_THREADS_PROPERTY = "packstar.pac.sampling.threads";
    private static final String PAC_SAMPLING_LARGE_LAMBDA_PROPERTY = "packstar.pac.sampling.largeLambdaThreshold";
    private static final String PAC_SAMPLING_PROGRESS_PROPERTY = "packstar.pac.sampling.progress";
    private static final String PAC_SAMPLING_GPU_PROPERTY = "packstar.pac.sampling.gpu";
    // Train/adaptation/final sizing shared by the production conditional path
    // and the archived legacy conditional implementation.
    private static final String PAC_TRAIN_SAMPLES_PROPERTY = "packstar.pac.trainSamples";
    private static final String PAC_PILOT_SAMPLES_PROPERTY = "packstar.pac.pilotSamples";
    private static final String PAC_MAX_EST_SAMPLES_PROPERTY = "packstar.pac.maxEstSamples";
    private static final String PAC_NSTAR_INFLATE_PROPERTY = "packstar.pac.nstarInflate";
    // Legacy fixed conditional bound.  Production instead reads the capped
    // residualBoundGrid and freezes one grid value from adaptation data.
    private static final String PAC_RESIDUAL_BOUND_PROPERTY = "packstar.pac.residualBound";
    // Practical certified-or-abort path.  Eta and B are adaptation-only;
    // monitor/final samples are fresh.  The certificate is conditional on the
    // selected global one-sided residual bound, and any observed violation
    // aborts instead of falling back to q_m or deterministic search.
    private static final String PAC_CONDITIONAL_REPAIR_AUDIT_PROPERTY =
            "packstar.pac.conditionalRepairAudit";
    private static final String PAC_RESIDUAL_BOUND_GRID_PROPERTY =
            "packstar.pac.residualBoundGrid";
    private static final String PAC_MONITOR_SAMPLES_PROPERTY =
            "packstar.pac.monitorSamples";
    // Frequency/severity production: choose a count/context-shrunk eta without
    // candidate-specific DP or CCD, then use bounded on-policy refits to repair
    // proposal-support shift.  Every adaptation batch records its exact source
    // proposal; monitor/final samples remain fresh after the proposal freezes.
    private static final String PAC_FREQUENCY_SEVERITY_PRODUCTION_PROPERTY =
            "packstar.pac.frequencySeverityProduction";
    private static final String PAC_FREQUENCY_SEVERITY_RELATIVE_BOUND_PROPERTY =
            "packstar.pac.frequencySeverity.relativeBoundKcal";
    private static final String PAC_FREQUENCY_SEVERITY_CAP_PROPERTY =
            "packstar.pac.frequencySeverity.severityCap";
    private static final String PAC_FREQUENCY_SEVERITY_PREMISE_ID_PROPERTY =
            "packstar.pac.frequencySeverity.severityPremiseId";
    private static final String PAC_FREQUENCY_SEVERITY_SHRINK_GRID_PROPERTY =
            "packstar.pac.frequencySeverity.shrinkGrid";
    private static final String PAC_FREQUENCY_SEVERITY_ALPHA_GRID_PROPERTY =
            "packstar.pac.frequencySeverity.alphaGrid";
    private static final String PAC_FREQUENCY_SEVERITY_FOLDS_PROPERTY =
            "packstar.pac.frequencySeverity.folds";
    private static final String PAC_FREQUENCY_SEVERITY_MIN_SHIFT_ESS_FRACTION_PROPERTY =
            "packstar.pac.frequencySeverity.minShiftEssFraction";
    private static final String PAC_FREQUENCY_SEVERITY_DISCOVERY_MIN_SHIFT_ESS_FRACTION_PROPERTY =
            "packstar.pac.frequencySeverity.discoveryMinShiftEssFraction";
    private static final String PAC_FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY =
            "packstar.pac.frequencySeverity.maxRefits";
    private static final String PAC_FREQUENCY_SEVERITY_DISCOVERY_SAMPLES_PROPERTY =
            "packstar.pac.frequencySeverity.discoverySamples";
    private static final String PAC_FREQUENCY_SEVERITY_DISCOVERY_MAX_SAMPLES_PROPERTY =
            "packstar.pac.frequencySeverity.discoveryMaxSamples";
    private static final String PAC_FREQUENCY_SEVERITY_VALIDATION_SAMPLES_PROPERTY =
            "packstar.pac.frequencySeverity.validationSamples";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEta";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaScale";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_GRID_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaScaleGrid";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_ASSIGNMENTS_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaMaxAssignments";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_POSITION_TRIPLES_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaMaxPositionTriples";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_FILL_EDGES_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaMaxFillEdges";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MIN_CELL_CONTEXTS_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaMinCellContexts";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_PRIOR_STRENGTH_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaPriorStrength";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_LOCAL_CAP_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaLocalCapKcal";
    private static final String PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_RESIDUAL_CAP_PROPERTY =
            "packstar.pac.frequencySeverity.tripleEtaResidualCapKcal";
    private static final String PAC_FREQUENCY_SEVERITY_MIN_TRAIN_COUNT_PROPERTY =
            "packstar.pac.frequencySeverity.minTrainCount";
    private static final String PAC_FREQUENCY_SEVERITY_MAX_UNDERTRAINED_AMPLIFICATION_PROPERTY =
            "packstar.pac.frequencySeverity.maxUndertrainedAmplification";
    private static final String PAC_FREQUENCY_SEVERITY_SIZE_SAFETY_PROPERTY =
            "packstar.pac.frequencySeverity.sizeSafety";
    private static final String PAC_FREQUENCY_SEVERITY_TEST_ALPHA_PROPERTY =
            "packstar.pac.frequencySeverity.severityTestAlpha";
    private static final String PAC_FREQUENCY_SEVERITY_OUTPUT_DIR_PROPERTY =
            "packstar.pac.frequencySeverity.outputDir";
    // Reuse the same adaptation CCD pool to screen a pre-registered family of
    // eta-repair candidates.  Each candidate gets an exact DP normalizer; no
    // additional CCD is used for screening, and the winner is frozen before
    // monitor/final data exist.
    private static final String PAC_DP_RICH_PROPERTY =
            "packstar.pac.dpRich";
    private static final String PAC_DP_RICH_MIN_ESS_FRACTION_PROPERTY =
            "packstar.pac.dpRich.minEssFraction";
    // v3 keeps every proposal that actually generated adaptation samples and
    // selects an on-policy shape by predicted certificate cost / chi-square
    // proxy.  A fresh winner-calibration batch is still adaptation data; it is
    // drawn before monitor/final and may only shift the proposal energy by a
    // conformation-independent gauge (which leaves the proposal unchanged).
    private static final String PAC_ETA_V3_PROPERTY =
            "packstar.pac.etaV3";
    private static final String PAC_ETA_V3_CALIBRATION_SAMPLES_PROPERTY =
            "packstar.pac.etaV3.calibrationSamples";
    private static final String PAC_ETA_V3_PROOF_BOUND_PROPERTY =
            "packstar.pac.etaV3.proofBoundKcal";
    private static final String PAC_ETA_V3_HEADROOM_PROPERTY =
            "packstar.pac.etaV3.headroomKcal";
    // v4 learns proposal shape against the importance-weight second moment,
    // then freezes an IID mixture of exact-DP components.  One anchor retains
    // the conditional pointwise B premise and supplies the final range cap.
    private static final String PAC_ETA_V4_PROPERTY =
            "packstar.pac.etaV4";
    private static final String PAC_ETA_V4_CALIBRATION_SAMPLES_PROPERTY =
            "packstar.pac.etaV4.calibrationSamples";
    private static final String PAC_ETA_V4_PROOF_BOUND_PROPERTY =
            "packstar.pac.etaV4.proofBoundKcal";
    private static final String PAC_ETA_V4_HEADROOM_PROPERTY =
            "packstar.pac.etaV4.headroomKcal";
    private static final String PAC_ETA_V4_TAIL_POWERS_PROPERTY =
            "packstar.pac.etaV4.tailPowers";
    private static final String PAC_ETA_V4_TRUST_KCAL_PROPERTY =
            "packstar.pac.etaV4.trustKcal";
    private static final String PAC_ETA_V4_ANCHOR_WEIGHTS_PROPERTY =
            "packstar.pac.etaV4.anchorWeights";
    private static final String PAC_ETA_V4_MIN_TAIL_ESS_FRACTION_PROPERTY =
            "packstar.pac.etaV4.minTailEssFraction";
    private static final String PAC_ETA_V4_OAIS_MAX_ROUNDS_PROPERTY =
            "packstar.pac.etaV4.oaisMaxRounds";
    // v5 Phase A: intercept the final single OAIS proposal before the
    // conditional-B mixture certificate and exhaustively measure minimizing
    // mass above a frozen eta-threshold grid. This is an oracle-only
    // diagnostic, not yet the scalable certified threshold DP.
    private static final String PAC_ETA_V5_THRESHOLD_ORACLE_PROPERTY =
            "packstar.pac.etaV5.thresholdOracle";
    private static final String PAC_ETA_V5_THRESHOLDS_PROPERTY =
            "packstar.pac.etaV5.thresholdsKcal";
    private static final String PAC_ETA_V5_ORACLE_OUTPUT_PROPERTY =
            "packstar.pac.etaV5.oracleOutputTsv";
    private static final String PAC_ETA_V5_MAX_STATES_PROPERTY =
            "packstar.pac.etaV5.maxEnumeratedStates";
    private static final String PAC_ETA_V5_REBASE_INTERVAL_PROPERTY =
            "packstar.pac.etaV5.rebaseInterval";
    private static final String PAC_ETA_V5_LOGZ_TOLERANCE_PROPERTY =
            "packstar.pac.etaV5.logZTolerance";
    private static final String PAC_ETA_V5_MECHANISM_PILOT_PROPERTY =
            "packstar.pac.etaV5.mechanismPilot";
    private static final String PAC_ETA_V5_MECHANISM_SAMPLES_PROPERTY =
            "packstar.pac.etaV5.mechanismSamples";
    private static final String PAC_ETA_V5_MECHANISM_ASSIGNMENTS_PROPERTY =
            "packstar.pac.etaV5.mechanismAssignmentsTsv";
    private static final String PAC_ETA_V5_MECHANISM_SUMMARY_PROPERTY =
            "packstar.pac.etaV5.mechanismSummaryTsv";
    // Diagnostic-only fresh-IID audit of the frozen proposal's empirical
    // importance-weight second moment and tail.  It never changes proposal
    // selection, sample sizing, support decisions, or certificate output.
    private static final String PAC_MOMENT_TAIL_STAGE_PROPERTY =
            "packstar.pac.momentTailAudit.stage";
    private static final String PAC_MOMENT_TAIL_THRESHOLDS_PROPERTY =
            "packstar.pac.momentTailAudit.thresholdsKcal";
    private static final String PAC_MOMENT_TAIL_ASSIGNMENTS_PROPERTY =
            "packstar.pac.momentTailAudit.assignmentsTsv";
    private static final String PAC_MOMENT_TAIL_SUMMARY_PROPERTY =
            "packstar.pac.momentTailAudit.summaryTsv";
    // B-free certificate: estimate the clipped bulk under p_eta and the exact
    // omitted tail under a pairwise exponential tilt of p_m.
    private static final String PAC_UNCONDITIONAL_TAIL_PROPERTY =
            "packstar.pac.unconditionalTail";
    private static final String PAC_BULK_SAMPLES_PROPERTY =
            "packstar.pac.bulk.samples";
    private static final String PAC_TAIL_SAMPLES_PROPERTY =
            "packstar.pac.tail.samples";
    private static final String PAC_TAIL_TILT_LAMBDAS_PROPERTY =
            "packstar.pac.tail.tiltLambdas";
    private static final String PAC_ETA_REPAIR_ENABLED_PROPERTY =
            "packstar.pac.etaRepair";
    private static final String PAC_ETA_REPAIR_MAX_ROUNDS_PROPERTY =
            "packstar.pac.etaRepair.maxRounds";
    private static final String PAC_ETA_REPAIR_TRIGGER_PROPERTY =
            "packstar.pac.etaRepair.triggerKcal";
    private static final String PAC_ETA_REPAIR_TARGET_PROPERTY =
            "packstar.pac.etaRepair.targetKcal";
    private static final String PAC_ETA_REPAIR_SAFETY_PROPERTY =
            "packstar.pac.etaRepair.safety";
    private static final String PAC_ETA_REPAIR_MIN_ALPHA_PROPERTY =
            "packstar.pac.etaRepair.minAlpha";
    // Ablation switch (PACK* paper Table 3 / Fig 1 "no-eta" baseline): when false,
    // Phase 3 (eta correction learning) is skipped and eta is fixed to the zero
    // correction, so E_eta == E_m, q_eta == q_m, and Stage B reweights the full
    // gap g = E_true - E_m sampled from p_m instead of the residual xi = E_true -
    // E_eta sampled from p_eta (see SI Lemma 2, Remark "the residual factor is the
    // only estimated object": the two-stage identity holds for any eta, including
    // eta === 0). Does not change any other phase or the PAC sample-splitting logic.
    private static final String PAC_ETA_ENABLED_PROPERTY = "packstar.pac.etaEnabled";
    private static final String DP_PARALLEL_THREADS_PROPERTY = "packstar.dp.parallel.threads";
    private static final String DP_CACHE_ENABLED_PROPERTY = "packstar.dp.cache";
    private static final String DP_CACHE_MAX_ENTRIES_PROPERTY = "packstar.dp.cache.maxEntries";
    private static final String DP_CACHE_MAX_TABLE_BYTES_PROPERTY = "packstar.dp.cache.maxTableBytes";
    private static final String DP_CACHE_MAX_TOTAL_BYTES_PROPERTY = "packstar.dp.cache.maxTotalBytes";
    private static final String DP_CACHE_SKIP_IF_M_STATES_PROPERTY = "packstar.dp.cache.skipIfMStates";
    private static final int DEFAULT_SAMPLES = 1000;
    private static final double DEFAULT_CONFIDENCE = 0.05; // delta = 0.05 => 95% confidence
    private static final double DEFAULT_TARGET_EPSILON = 0.683;
    private static final long DEFAULT_RANDOM_SEED = 42L;
    private static final int DEFAULT_SAMPLING_LARGE_LAMBDA = 65_536;
    private static final double DEFAULT_TRAIN_FRACTION = 0.5;
    private static final double DEFAULT_PILOT_FRACTION = 0.1;
    private static final int DEFAULT_MAX_EST_SAMPLES = 4000;
    private static final double DEFAULT_NSTAR_INFLATE = 1.3;
    private static final double DEFAULT_RESIDUAL_BOUND = 1.0;
    private static final boolean DEFAULT_CONDITIONAL_REPAIR_AUDIT = true;
    private static final String DEFAULT_RESIDUAL_BOUND_GRID = "1,1.5,2";
    private static final int DEFAULT_MONITOR_SAMPLES = 100;
    private static final boolean DEFAULT_FREQUENCY_SEVERITY_PRODUCTION = false;
    private static final double DEFAULT_FREQUENCY_SEVERITY_RELATIVE_BOUND_KCAL = 1.0;
    private static final double DEFAULT_FREQUENCY_SEVERITY_CAP = 20.0;
    private static final String DEFAULT_FREQUENCY_SEVERITY_PREMISE_ID =
            "conditional-relative-gauge-S0-20-not-externally-recalibrated";
    private static final String DEFAULT_FREQUENCY_SEVERITY_SHRINK_GRID =
            "0:0,2:5,5:10,10:20";
    private static final String DEFAULT_FREQUENCY_SEVERITY_ALPHA_GRID = "1";
    private static final int DEFAULT_FREQUENCY_SEVERITY_FOLDS = 2;
    private static final double DEFAULT_FREQUENCY_SEVERITY_MIN_SHIFT_ESS_FRACTION = 0.25;
    private static final double DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MIN_SHIFT_ESS_FRACTION = 0.05;
    private static final int DEFAULT_FREQUENCY_SEVERITY_MAX_REFITS = 8;
    private static final int DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_SAMPLES = 100;
    private static final int DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MAX_SAMPLES = 400;
    private static final int DEFAULT_FREQUENCY_SEVERITY_VALIDATION_SAMPLES = 400;
    private static final boolean DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA = false;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE = 1.0;
    private static final String DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_GRID = "1";
    private static final long DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_ASSIGNMENTS = 1_000_000L;
    private static final int DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_POSITION_TRIPLES = 3;
    private static final int DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_FILL_EDGES = 3;
    private static final int DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MIN_CELL_CONTEXTS = 1;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_PRIOR_STRENGTH = 4.0;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_LOCAL_CAP_KCAL = 2.0;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_RESIDUAL_CAP_KCAL = 3.0;
    private static final int DEFAULT_FREQUENCY_SEVERITY_MIN_TRAIN_COUNT = 5;
    private static final double DEFAULT_FREQUENCY_SEVERITY_MAX_UNDERTRAINED_AMPLIFICATION = 1.25;
    private static final double DEFAULT_FREQUENCY_SEVERITY_SIZE_SAFETY = 0.9;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TEST_ALPHA = 0.05;
    // v1 remains the production default.  A dedicated v2 runner must opt in
    // after freezing the candidate-screening and reweighting policy.
    private static final boolean DEFAULT_DP_RICH = false;
    private static final double DEFAULT_DP_RICH_MIN_ESS_FRACTION = 0.25;
    private static final boolean DEFAULT_ETA_V3 = false;
    private static final int DEFAULT_ETA_V3_CALIBRATION_SAMPLES = 200;
    private static final double DEFAULT_ETA_V3_PROOF_BOUND_KCAL = 2.0;
    private static final double DEFAULT_ETA_V3_HEADROOM_KCAL = 0.25;
    private static final boolean DEFAULT_ETA_V4 = false;
    private static final int DEFAULT_ETA_V4_CALIBRATION_SAMPLES = 200;
    private static final double DEFAULT_ETA_V4_PROOF_BOUND_KCAL = 2.0;
    private static final double DEFAULT_ETA_V4_HEADROOM_KCAL = 0.25;
    private static final String DEFAULT_ETA_V4_TAIL_POWERS = "1,2";
    private static final String DEFAULT_ETA_V4_TRUST_KCAL = "0.125,0.25,0.5,1";
    private static final String DEFAULT_ETA_V4_ANCHOR_WEIGHTS = "0.5,0.8,0.9,0.95,1";
    private static final double DEFAULT_ETA_V4_MIN_TAIL_ESS_FRACTION = 0.02;
    private static final int DEFAULT_ETA_V4_OAIS_MAX_ROUNDS = 3;
    private static final boolean DEFAULT_ETA_V5_THRESHOLD_ORACLE = false;
    private static final String DEFAULT_ETA_V5_THRESHOLDS_KCAL =
            "-4,-3,-2,-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,"
                    + "2,2.25,2.5,2.75,3,3.5,4,5,6,8,12";
    private static final long DEFAULT_ETA_V5_MAX_ENUMERATED_STATES =
            50_000_000L;
    private static final int DEFAULT_ETA_V5_REBASE_INTERVAL = 1_048_576;
    private static final double DEFAULT_ETA_V5_LOGZ_TOLERANCE = 1.0e-6;
    private static final boolean DEFAULT_ETA_V5_MECHANISM_PILOT = false;
    private static final int DEFAULT_ETA_V5_MECHANISM_SAMPLES = 200;
    private static final boolean DEFAULT_UNCONDITIONAL_TAIL = false;
    // Provisional development defaults. Formal jobs override these with the
    // fixed budget frozen after production-scale micro-pilots.
    private static final int DEFAULT_BULK_SAMPLES = 128;
    private static final int DEFAULT_TAIL_SAMPLES = 128;
    private static final String DEFAULT_TAIL_TILT_LAMBDAS = "0,0.25,0.5,1,2,4";
    private static final boolean DEFAULT_ETA_REPAIR_ENABLED = true;
    private static final int DEFAULT_ETA_REPAIR_MAX_ROUNDS = 2;
    private static final double DEFAULT_ETA_REPAIR_TRIGGER_KCAL = 0.5;
    private static final double DEFAULT_ETA_REPAIR_TARGET_KCAL = 0.0;
    private static final double DEFAULT_ETA_REPAIR_SAFETY = 0.95;
    private static final double DEFAULT_ETA_REPAIR_MIN_ALPHA = 0.0;
    private static final boolean DEFAULT_ETA_ENABLED = true;
    private static final boolean DEFAULT_CLIP = true;
    private static final double DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE = 0.85;
    private static final boolean DEFAULT_ITERATE = true;
    private static final int DEFAULT_ITERATE_MAX_ROUNDS = 4;
    private static final double DEFAULT_SIZE_SAFETY = 0.9;
    private static final int DEFAULT_DP_CACHE_MAX_ENTRIES = 20000;
    private static final long DEFAULT_DP_CACHE_MAX_TABLE_BYTES = 256L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_MAX_TOTAL_BYTES = 4L * 1024L * 1024L * 1024L;
    private static final long DEFAULT_DP_CACHE_SKIP_IF_M_STATES = 8_000_000L;

    private static final Object CORRECTED_DP_CACHE_LOCK = new Object();
    private static final LinkedHashMap<String, CachedDPTable> CORRECTED_DP_CACHE =
            new LinkedHashMap<>(1024, 0.75f, true);
    private static long correctedDPCacheBytes = 0L;

    // Inputs
    private RootedTreeNode rootedRoot;
    private RootedTreeEdge rootedRootEdge;
    private final RootedTreeNode initialRootedRoot;
    private final RootedTreeEdge initialRootedRootEdge;
    private final EnergyMatrix branchMinimizingEmat;
    private final EnergyMatrix branchRigidEmat;
    private final InteractionGraph interactionGraph;
    private InteractionGraph proposalInteractionGraph;
    private String proposalInteractionGraphSignature;
    private final ConfEnergyCalculator minimizingEcalc;
    private final RCs rcs;
    private final SimpleConfSpace confSpace;
    private final double RT;
    /**
     * Energy matrix represented by the DP tables currently loaded in rootedRoot.
     * CPU ancestral sampling must use this exact matrix for non-leaf local terms;
     * otherwise the conditional weights disagree with the child DP normalizers.
     */
    private EnergyMatrix activeProposalEmat;

    // Configuration
    private final int numSamples;
    private final double delta; // confidence parameter
    private final double targetEpsilon;
    private final long randomSeed;
    private final String randomStreamIdentity;
    private final boolean batchedSampling;
    private final boolean gpuSampling;
    private final int samplingThreads;
    private final int samplingLargeLambdaThreshold;
    private final boolean samplingProgress;
    private final int trainSamples;
    private final int pilotSamples;
    private final int maxEstSamples;
    private final int sampleBudget;
    private double clipLogCap = Double.NaN; // dynamic B/RT, frozen independently of final sample values
    private int pilotResidualBoundViolations = 0;
    private double pilotMaxOverCorrectionKcal = Double.NEGATIVE_INFINITY;
    private double logZMinDet = Double.POSITIVE_INFINITY; // log q_m: deterministic, assumption-free upper bound on log q
    private final double nstarInflate;
    private final double residualBoundKcal;
    private final boolean conditionalRepairAudit;
    private final double[] residualBoundGridKcal;
    private final int monitorSamples;
    private final boolean frequencySeverityProduction;
    private final double frequencySeverityRelativeBoundKcal;
    private final double frequencySeverityCap;
    private final String frequencySeverityPremiseId;
    private final FrequencySeverityShrinkPair[] frequencySeverityShrinkGrid;
    private final double[] frequencySeverityAlphaGrid;
    private final int frequencySeverityFolds;
    private final double frequencySeverityMinShiftEssFraction;
    private final double frequencySeverityDiscoveryMinShiftEssFraction;
    private final int frequencySeverityMaxRefits;
    private final int frequencySeverityDiscoverySamples;
    private final int frequencySeverityDiscoveryMaxSamples;
    private final int frequencySeverityValidationSamples;
    private final boolean frequencySeverityTripleEtaEnabled;
    private final double frequencySeverityTripleEtaScale;
    private final double[] frequencySeverityTripleEtaScaleGrid;
    private final long frequencySeverityTripleEtaMaxAssignments;
    private final int frequencySeverityTripleEtaMaxPositionTriples;
    private final int frequencySeverityTripleEtaMaxFillEdges;
    private final int frequencySeverityTripleEtaMinCellContexts;
    private final double frequencySeverityTripleEtaPriorStrength;
    private final double frequencySeverityTripleEtaLocalCapKcal;
    private final double frequencySeverityTripleEtaResidualCapKcal;
    private PackStarTripleEtaCorrections frequencySeverityTripleEta = null;
    private String frequencySeverityTripleEtaFallbackReason = null;
    private final int frequencySeverityMinTrainCount;
    private final double frequencySeverityMaxUndertrainedAmplification;
    private final double frequencySeveritySizeSafety;
    private final double frequencySeverityTestAlpha;
    private final String frequencySeverityOutputDir;
    private final boolean dpRich;
    private final double dpRichMinEssFraction;
    private final boolean etaV3;
    private final int etaV3CalibrationSamples;
    private final double etaV3ProofBoundKcal;
    private final double etaV3HeadroomKcal;
    private final boolean etaV4;
    private final int etaV4CalibrationSamples;
    private final double etaV4ProofBoundKcal;
    private final double etaV4HeadroomKcal;
    private final double[] etaV4TailPowers;
    private final double[] etaV4TrustKcal;
    private final double[] etaV4AnchorWeights;
    private final double etaV4MinTailEssFraction;
    private final int etaV4OaisMaxRounds;
    private final boolean etaV5ThresholdOracle;
    private final double[] etaV5ThresholdsKcal;
    private final String etaV5OracleOutputTsv;
    private final long etaV5MaxEnumeratedStates;
    private final int etaV5RebaseInterval;
    private final double etaV5LogZTolerance;
    private final boolean etaV5MechanismPilot;
    private final int etaV5MechanismSamples;
    private final String etaV5MechanismAssignmentsTsv;
    private final String etaV5MechanismSummaryTsv;
    private double selectedResidualBoundKcal = Double.NaN;
    private final boolean unconditionalTail;
    private final int fixedBulkSamples;
    private final int fixedTailSamples;
    private final double[] tailTiltLambdas;
    private final boolean etaRepairEnabled;
    private final int etaRepairMaxRounds;
    private final double etaRepairTriggerKcal;
    private final double etaRepairTargetKcal;
    private final double etaRepairSafety;
    private final double etaRepairMinAlpha;
    private final boolean etaEnabled; // false => no-eta ablation (Table 3 / Fig 1 baseline)
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
    private boolean certificateValid = false;
    private String certificateFailureReason = "not computed";

    private static final BoltzmannCalculator LOG_SPACE_EXP =
            new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private final LocalRCMap[] localRCByGlobalRC;
    private PackStarSampleListener sampleListener = null;

    public PackStarEstimator(RootedTreeNode rootedRoot,
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

    public PackStarEstimator(RootedTreeNode rootedRoot,
                                RootedTreeEdge rootedRootEdge,
                                EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace,
                                 double requestedTargetEpsilon) {
        this(rootedRoot, rootedRootEdge,
                branchMinimizingEmat, branchRigidEmat,
                interactionGraph, minimizingEcalc,
                rcs, confSpace,
                requestedTargetEpsilon,
                Integer.MAX_VALUE);
    }

    public PackStarEstimator(RootedTreeNode rootedRoot,
                                RootedTreeEdge rootedRootEdge,
                                EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace,
                                 double requestedTargetEpsilon,
                                 int requestedSampleBudget) {
        this(rootedRoot, rootedRootEdge,
                branchMinimizingEmat, branchRigidEmat,
                interactionGraph, minimizingEcalc,
                rcs, confSpace,
                requestedTargetEpsilon,
                requestedSampleBudget,
                defaultRandomStreamIdentity(rcs));
    }

    /**
     * Construct an estimator with a stable logical random-stream identity.
     *
     * <p>The identity describes the calculation rather than its process-local
     * construction order. K-star supplies a state role plus the exact allowed-RC
     * signature, so shard count and rank-to-node mapping do not change the
     * estimator's random stream.</p>
     */
    public PackStarEstimator(RootedTreeNode rootedRoot,
                                RootedTreeEdge rootedRootEdge,
                                EnergyMatrix branchMinimizingEmat,
                                 EnergyMatrix branchRigidEmat,
                                 InteractionGraph interactionGraph,
                                 ConfEnergyCalculator minimizingEcalc,
                                 RCs rcs,
                                 SimpleConfSpace confSpace,
                                 double requestedTargetEpsilon,
                                 int requestedSampleBudget,
                                 String randomStreamIdentity) {
        this.rootedRoot = rootedRoot;
        this.rootedRootEdge = rootedRootEdge;
        this.initialRootedRoot = rootedRoot;
        this.initialRootedRootEdge = rootedRootEdge;
        this.branchMinimizingEmat = branchMinimizingEmat;
        this.branchRigidEmat = branchRigidEmat;
        this.interactionGraph = interactionGraph;
        this.proposalInteractionGraph = interactionGraph;
        this.proposalInteractionGraphSignature =
                interactionGraphSignature(interactionGraph);
        this.minimizingEcalc = minimizingEcalc;
        this.rcs = rcs;
        this.confSpace = confSpace;
        this.RT = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
        this.activeProposalEmat = branchMinimizingEmat;
        this.sampleBudget = sanitizeSampleBudget(requestedSampleBudget);

        this.numSamples = capSampleCount(
                Math.max(1, getConfigInteger(PAC_SAMPLES_PROPERTY, DEFAULT_SAMPLES)),
                1);
        String deltaStr = getConfigProperty(PAC_CONFIDENCE_PROPERTY, null);
        this.delta = (deltaStr != null) ? Double.parseDouble(deltaStr) : DEFAULT_CONFIDENCE;
        double defaultTarget = Double.isFinite(requestedTargetEpsilon) && requestedTargetEpsilon > 0.0
                ? requestedTargetEpsilon
                : DEFAULT_TARGET_EPSILON;
        this.targetEpsilon = Math.max(0.0,
                getConfigDouble(PAC_TARGET_EPSILON_PROPERTY, defaultTarget));
        long baseRandomSeed = getConfigLong(PAC_RANDOM_SEED_PROPERTY, DEFAULT_RANDOM_SEED);
        this.randomStreamIdentity = requireRandomStreamIdentity(randomStreamIdentity);
        this.randomSeed = deriveRandomSeed(baseRandomSeed, this.randomStreamIdentity);
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
        int defaultMaxEstSamples = Math.max(2, DEFAULT_MAX_EST_SAMPLES);
        this.trainSamples = capSampleCount(
                Math.max(2, getConfigInteger(PAC_TRAIN_SAMPLES_PROPERTY, defaultTrainSamples)),
                2);
        this.pilotSamples = capSampleCount(
                Math.max(2, getConfigInteger(PAC_PILOT_SAMPLES_PROPERTY, defaultPilotSamples)),
                2);
        this.maxEstSamples = capSampleCount(
                Math.max(2, getConfigInteger(PAC_MAX_EST_SAMPLES_PROPERTY, defaultMaxEstSamples)),
                2);
        this.nstarInflate = Math.max(1.0, getConfigDouble(PAC_NSTAR_INFLATE_PROPERTY, DEFAULT_NSTAR_INFLATE));
        double configuredResidualBound = getConfigDouble(PAC_RESIDUAL_BOUND_PROPERTY, DEFAULT_RESIDUAL_BOUND);
        this.residualBoundKcal = configuredResidualBound >= 0.0 ? configuredResidualBound : Double.NaN;
        this.conditionalRepairAudit = getConfigBoolean(
                PAC_CONDITIONAL_REPAIR_AUDIT_PROPERTY,
                DEFAULT_CONDITIONAL_REPAIR_AUDIT);
        this.residualBoundGridKcal = parsePositiveGrid(
                getConfigProperty(PAC_RESIDUAL_BOUND_GRID_PROPERTY,
                        DEFAULT_RESIDUAL_BOUND_GRID),
                "PACK* residual-bound grid");
        this.monitorSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_MONITOR_SAMPLES_PROPERTY,
                        DEFAULT_MONITOR_SAMPLES)), 2);
        this.frequencySeverityProduction = getConfigBoolean(
                PAC_FREQUENCY_SEVERITY_PRODUCTION_PROPERTY, DEFAULT_FREQUENCY_SEVERITY_PRODUCTION);
        this.frequencySeverityRelativeBoundKcal = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_RELATIVE_BOUND_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_RELATIVE_BOUND_KCAL);
        this.frequencySeverityCap = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_CAP_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_CAP);
        this.frequencySeverityPremiseId = getConfigProperty(
                PAC_FREQUENCY_SEVERITY_PREMISE_ID_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_PREMISE_ID).trim();
        this.frequencySeverityShrinkGrid = parseFrequencySeverityShrinkGrid(getConfigProperty(
                PAC_FREQUENCY_SEVERITY_SHRINK_GRID_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_SHRINK_GRID));
        this.frequencySeverityAlphaGrid = parseFixedOneGrid(getConfigProperty(
                PAC_FREQUENCY_SEVERITY_ALPHA_GRID_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_ALPHA_GRID),
                "PACK* frequency/severity alpha grid");
        this.frequencySeverityFolds = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_FOLDS_PROPERTY, DEFAULT_FREQUENCY_SEVERITY_FOLDS);
        this.frequencySeverityMinShiftEssFraction = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_MIN_SHIFT_ESS_FRACTION_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_MIN_SHIFT_ESS_FRACTION);
        this.frequencySeverityDiscoveryMinShiftEssFraction = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_DISCOVERY_MIN_SHIFT_ESS_FRACTION_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MIN_SHIFT_ESS_FRACTION);
        this.frequencySeverityMaxRefits = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_MAX_REFITS);
        this.frequencySeverityDiscoverySamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_FREQUENCY_SEVERITY_DISCOVERY_SAMPLES_PROPERTY,
                        DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_SAMPLES)), 2);
        this.frequencySeverityDiscoveryMaxSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_FREQUENCY_SEVERITY_DISCOVERY_MAX_SAMPLES_PROPERTY,
                        DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MAX_SAMPLES)), 2);
        this.frequencySeverityValidationSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_FREQUENCY_SEVERITY_VALIDATION_SAMPLES_PROPERTY,
                        DEFAULT_FREQUENCY_SEVERITY_VALIDATION_SAMPLES)), 2);
        this.frequencySeverityTripleEtaEnabled = getConfigBoolean(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA);
        this.frequencySeverityTripleEtaScale = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE);
        this.frequencySeverityTripleEtaScaleGrid = parseFixedOneGrid(
                getConfigProperty(
                        PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_GRID_PROPERTY,
                        DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_GRID),
                "PACK* frequency/severity triple eta scale grid");
        this.frequencySeverityTripleEtaMaxAssignments = getConfigLong(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_ASSIGNMENTS_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_ASSIGNMENTS);
        this.frequencySeverityTripleEtaMaxPositionTriples = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_POSITION_TRIPLES_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_POSITION_TRIPLES);
        this.frequencySeverityTripleEtaMaxFillEdges = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_FILL_EDGES_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_FILL_EDGES);
        this.frequencySeverityTripleEtaMinCellContexts = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_MIN_CELL_CONTEXTS_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MIN_CELL_CONTEXTS);
        this.frequencySeverityTripleEtaPriorStrength = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_PRIOR_STRENGTH_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_PRIOR_STRENGTH);
        this.frequencySeverityTripleEtaLocalCapKcal = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_LOCAL_CAP_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_LOCAL_CAP_KCAL);
        this.frequencySeverityTripleEtaResidualCapKcal = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_RESIDUAL_CAP_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_RESIDUAL_CAP_KCAL);
        this.frequencySeverityMinTrainCount = getConfigInteger(
                PAC_FREQUENCY_SEVERITY_MIN_TRAIN_COUNT_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_MIN_TRAIN_COUNT);
        this.frequencySeverityMaxUndertrainedAmplification = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_MAX_UNDERTRAINED_AMPLIFICATION_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_MAX_UNDERTRAINED_AMPLIFICATION);
        this.frequencySeveritySizeSafety = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_SIZE_SAFETY_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_SIZE_SAFETY);
        this.frequencySeverityTestAlpha = getConfigDouble(
                PAC_FREQUENCY_SEVERITY_TEST_ALPHA_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TEST_ALPHA);
        String configuredFrequencySeverityOutputDir = getConfigProperty(
                PAC_FREQUENCY_SEVERITY_OUTPUT_DIR_PROPERTY, null);
        this.frequencySeverityOutputDir = configuredFrequencySeverityOutputDir == null
                ? null : configuredFrequencySeverityOutputDir.trim();
        if (!Double.isFinite(frequencySeverityRelativeBoundKcal)
                || !(frequencySeverityRelativeBoundKcal > 0.0)) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity relative bound must be finite and positive: "
                            + frequencySeverityRelativeBoundKcal);
        }
        if (!Double.isFinite(frequencySeverityCap) || frequencySeverityCap < 0.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity severity cap must be finite and nonnegative: "
                            + frequencySeverityCap);
        }
        if (frequencySeverityPremiseId.isEmpty()) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity severity premise id must not be empty");
        }
        if (frequencySeverityFolds != 2) {
            throw new IllegalArgumentException(
                    "AdaptiveFrequencySeverityPAC currently freezes exactly two folds; got "
                            + frequencySeverityFolds);
        }
        if (!Double.isFinite(frequencySeverityMinShiftEssFraction)
                || !(frequencySeverityMinShiftEssFraction > 0.0)
                || frequencySeverityMinShiftEssFraction > 1.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity minimum shift ESS fraction must be in (0,1]: "
                            + frequencySeverityMinShiftEssFraction);
        }
        if (!Double.isFinite(frequencySeverityDiscoveryMinShiftEssFraction)
                || !(frequencySeverityDiscoveryMinShiftEssFraction > 0.0)
                || frequencySeverityDiscoveryMinShiftEssFraction
                > frequencySeverityMinShiftEssFraction) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity discovery shift ESS fraction must be in"
                            + " (0, minimum final shift ESS fraction]: "
                            + frequencySeverityDiscoveryMinShiftEssFraction);
        }
        if (frequencySeverityMaxRefits < 0 || frequencySeverityMaxRefits > 16) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity maximum refits must be in [0,16]: "
                            + frequencySeverityMaxRefits);
        }
        if (frequencySeverityDiscoveryMaxSamples
                < frequencySeverityDiscoverySamples) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity maximum discovery samples must"
                            + " be at least the initial discovery samples: "
                            + frequencySeverityDiscoveryMaxSamples + " < "
                            + frequencySeverityDiscoverySamples);
        }
        if (!Double.isFinite(frequencySeverityTripleEtaScale)
                || frequencySeverityTripleEtaScale != 1.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity triple eta scale is not tunable"
                            + " and must be fixed to 1: "
                            + frequencySeverityTripleEtaScale);
        }
        if (frequencySeverityTripleEtaMaxAssignments < 1L) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity triple eta assignment cap must"
                            + " be positive: "
                            + frequencySeverityTripleEtaMaxAssignments);
        }
        if (frequencySeverityTripleEtaMaxPositionTriples < 1
                || frequencySeverityTripleEtaMaxPositionTriples > 16) {
            throw new IllegalArgumentException(
                    "PACK* selected triple count must be in [1,16]: "
                            + frequencySeverityTripleEtaMaxPositionTriples);
        }
        if (frequencySeverityTripleEtaMaxFillEdges < 0
                || frequencySeverityTripleEtaMaxFillEdges > 16) {
            throw new IllegalArgumentException(
                    "PACK* selected triple fill-edge budget must be in [0,16]: "
                            + frequencySeverityTripleEtaMaxFillEdges);
        }
        if (frequencySeverityTripleEtaMinCellContexts < 1) {
            throw new IllegalArgumentException(
                    "PACK* selected triple minimum cell contexts must be positive: "
                            + frequencySeverityTripleEtaMinCellContexts);
        }
        if (!Double.isFinite(frequencySeverityTripleEtaPriorStrength)
                || frequencySeverityTripleEtaPriorStrength < 0.0
                || !Double.isFinite(frequencySeverityTripleEtaLocalCapKcal)
                || frequencySeverityTripleEtaLocalCapKcal < 0.0
                || !Double.isFinite(frequencySeverityTripleEtaResidualCapKcal)
                || !(frequencySeverityTripleEtaResidualCapKcal > 0.0)) {
            throw new IllegalArgumentException(
                    "PACK* selected triple shrink/cap controls are invalid");
        }
        if (frequencySeverityMinTrainCount < 1) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity minimum train count must be positive: "
                            + frequencySeverityMinTrainCount);
        }
        if (!Double.isFinite(frequencySeverityMaxUndertrainedAmplification)
                || frequencySeverityMaxUndertrainedAmplification < 1.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity maximum undertrained amplification must be finite and >=1: "
                            + frequencySeverityMaxUndertrainedAmplification);
        }
        if (!Double.isFinite(frequencySeveritySizeSafety)
                || !(frequencySeveritySizeSafety > 0.0) || frequencySeveritySizeSafety > 1.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity sizing safety must be in (0,1]: "
                            + frequencySeveritySizeSafety);
        }
        if (!Double.isFinite(frequencySeverityTestAlpha)
                || !(frequencySeverityTestAlpha > 0.0)
                || frequencySeverityTestAlpha >= 1.0) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity severity-test alpha must be in (0,1): "
                            + frequencySeverityTestAlpha);
        }
        if (frequencySeverityProduction && (frequencySeverityOutputDir == null
                || !isXtmpOutputPath(frequencySeverityOutputDir))) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity production requires an absolute output directory under"
                            + " /usr/xtmp/lz280: " + frequencySeverityOutputDir);
        }
        this.dpRich = getConfigBoolean(PAC_DP_RICH_PROPERTY,
                DEFAULT_DP_RICH);
        this.dpRichMinEssFraction = getConfigDouble(
                PAC_DP_RICH_MIN_ESS_FRACTION_PROPERTY,
                DEFAULT_DP_RICH_MIN_ESS_FRACTION);
        if (!Double.isFinite(dpRichMinEssFraction)
                || !(dpRichMinEssFraction > 0.0)
                || dpRichMinEssFraction > 1.0) {
            throw new IllegalArgumentException(
                    "PACK* DP-rich minimum ESS fraction must be in (0,1]: "
                            + dpRichMinEssFraction);
        }
        this.etaV3 = getConfigBoolean(PAC_ETA_V3_PROPERTY,
                DEFAULT_ETA_V3);
        this.etaV3CalibrationSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_ETA_V3_CALIBRATION_SAMPLES_PROPERTY,
                        DEFAULT_ETA_V3_CALIBRATION_SAMPLES)), 2);
        this.etaV3ProofBoundKcal = getConfigDouble(
                PAC_ETA_V3_PROOF_BOUND_PROPERTY,
                DEFAULT_ETA_V3_PROOF_BOUND_KCAL);
        this.etaV3HeadroomKcal = getConfigDouble(
                PAC_ETA_V3_HEADROOM_PROPERTY,
                DEFAULT_ETA_V3_HEADROOM_KCAL);
        if (etaV3 && !dpRich) {
            throw new IllegalArgumentException(
                    "PACK* eta-v3 requires packstar.pac.dpRich=true");
        }
        if (etaV3 && (!Double.isFinite(etaV3ProofBoundKcal)
                || !(etaV3ProofBoundKcal > 0.0)
                || !gridContains(residualBoundGridKcal,
                etaV3ProofBoundKcal))) {
            throw new IllegalArgumentException(
                    "PACK* eta-v3 proof bound must be a member of the residual-bound grid: "
                            + etaV3ProofBoundKcal);
        }
        if (etaV3 && (!Double.isFinite(etaV3HeadroomKcal)
                || etaV3HeadroomKcal < 0.0
                || etaV3HeadroomKcal >= etaV3ProofBoundKcal)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v3 headroom must be finite and in [0,B): "
                            + etaV3HeadroomKcal);
        }
        this.etaV4 = getConfigBoolean(PAC_ETA_V4_PROPERTY,
                DEFAULT_ETA_V4);
        this.etaV4CalibrationSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_ETA_V4_CALIBRATION_SAMPLES_PROPERTY,
                        DEFAULT_ETA_V4_CALIBRATION_SAMPLES)), 2);
        this.etaV4ProofBoundKcal = getConfigDouble(
                PAC_ETA_V4_PROOF_BOUND_PROPERTY,
                DEFAULT_ETA_V4_PROOF_BOUND_KCAL);
        this.etaV4HeadroomKcal = getConfigDouble(
                PAC_ETA_V4_HEADROOM_PROPERTY,
                DEFAULT_ETA_V4_HEADROOM_KCAL);
        this.etaV4TailPowers = parsePositiveGrid(getConfigProperty(
                PAC_ETA_V4_TAIL_POWERS_PROPERTY,
                DEFAULT_ETA_V4_TAIL_POWERS),
                "PACK* eta-v4 tail-power grid");
        this.etaV4TrustKcal = parsePositiveGrid(getConfigProperty(
                PAC_ETA_V4_TRUST_KCAL_PROPERTY,
                DEFAULT_ETA_V4_TRUST_KCAL),
                "PACK* eta-v4 trust-kcal grid");
        this.etaV4AnchorWeights = parsePositiveGrid(getConfigProperty(
                PAC_ETA_V4_ANCHOR_WEIGHTS_PROPERTY,
                DEFAULT_ETA_V4_ANCHOR_WEIGHTS),
                "PACK* eta-v4 anchor-weight grid");
        this.etaV4MinTailEssFraction = getConfigDouble(
                PAC_ETA_V4_MIN_TAIL_ESS_FRACTION_PROPERTY,
                DEFAULT_ETA_V4_MIN_TAIL_ESS_FRACTION);
        this.etaV4OaisMaxRounds = getConfigInteger(
                PAC_ETA_V4_OAIS_MAX_ROUNDS_PROPERTY,
                DEFAULT_ETA_V4_OAIS_MAX_ROUNDS);
        if (etaV4 && !dpRich) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 requires packstar.pac.dpRich=true");
        }
        if (etaV4 && etaV3) {
            throw new IllegalArgumentException(
                    "PACK* eta-v3 and eta-v4 modes are mutually exclusive");
        }
        if (etaV4 && (!Double.isFinite(etaV4ProofBoundKcal)
                || !(etaV4ProofBoundKcal > 0.0)
                || !gridContains(residualBoundGridKcal,
                etaV4ProofBoundKcal))) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 proof bound must be a member of the residual-bound grid: "
                            + etaV4ProofBoundKcal);
        }
        if (etaV4 && (!Double.isFinite(etaV4HeadroomKcal)
                || etaV4HeadroomKcal < 0.0
                || etaV4HeadroomKcal >= etaV4ProofBoundKcal)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 headroom must be finite and in [0,B): "
                            + etaV4HeadroomKcal);
        }
        for (double power : etaV4TailPowers) {
            if (power > 2.0) {
                throw new IllegalArgumentException(
                        "PACK* eta-v4 tail powers must be in (0,2]: "
                                + power);
            }
        }
        for (double weight : etaV4AnchorWeights) {
            if (weight > 1.0) {
                throw new IllegalArgumentException(
                        "PACK* eta-v4 anchor weights must be in (0,1]: "
                                + weight);
            }
        }
        if (!gridContains(etaV4AnchorWeights, 1.0)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 anchor-weight grid must contain 1");
        }
        if (!Double.isFinite(etaV4MinTailEssFraction)
                || !(etaV4MinTailEssFraction > 0.0)
                || etaV4MinTailEssFraction > 1.0) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 minimum tail ESS fraction must be in (0,1]: "
                            + etaV4MinTailEssFraction);
        }
        if (etaV4OaisMaxRounds < 1 || etaV4OaisMaxRounds > 8) {
            throw new IllegalArgumentException(
                    "PACK* eta-v4 OAIS max rounds must be in [1,8]: "
                            + etaV4OaisMaxRounds);
        }
        this.etaV5ThresholdOracle = getConfigBoolean(
                PAC_ETA_V5_THRESHOLD_ORACLE_PROPERTY,
                DEFAULT_ETA_V5_THRESHOLD_ORACLE);
        this.etaV5ThresholdsKcal = parseFiniteIncreasingGrid(
                getConfigProperty(PAC_ETA_V5_THRESHOLDS_PROPERTY,
                        DEFAULT_ETA_V5_THRESHOLDS_KCAL),
                "PACK* eta-v5 threshold-kcal grid");
        String configuredEtaV5Output = getConfigProperty(
                PAC_ETA_V5_ORACLE_OUTPUT_PROPERTY, null);
        this.etaV5OracleOutputTsv = configuredEtaV5Output == null
                ? null : configuredEtaV5Output.trim();
        this.etaV5MaxEnumeratedStates = getConfigLong(
                PAC_ETA_V5_MAX_STATES_PROPERTY,
                DEFAULT_ETA_V5_MAX_ENUMERATED_STATES);
        this.etaV5RebaseInterval = getConfigInteger(
                PAC_ETA_V5_REBASE_INTERVAL_PROPERTY,
                DEFAULT_ETA_V5_REBASE_INTERVAL);
        this.etaV5LogZTolerance = getConfigDouble(
                PAC_ETA_V5_LOGZ_TOLERANCE_PROPERTY,
                DEFAULT_ETA_V5_LOGZ_TOLERANCE);
        this.etaV5MechanismPilot = getConfigBoolean(
                PAC_ETA_V5_MECHANISM_PILOT_PROPERTY,
                DEFAULT_ETA_V5_MECHANISM_PILOT);
        this.etaV5MechanismSamples = getConfigInteger(
                PAC_ETA_V5_MECHANISM_SAMPLES_PROPERTY,
                DEFAULT_ETA_V5_MECHANISM_SAMPLES);
        String configuredEtaV5MechanismAssignments = getConfigProperty(
                PAC_ETA_V5_MECHANISM_ASSIGNMENTS_PROPERTY, null);
        this.etaV5MechanismAssignmentsTsv =
                configuredEtaV5MechanismAssignments == null
                        ? null
                        : configuredEtaV5MechanismAssignments.trim();
        String configuredEtaV5MechanismSummary = getConfigProperty(
                PAC_ETA_V5_MECHANISM_SUMMARY_PROPERTY, null);
        this.etaV5MechanismSummaryTsv =
                configuredEtaV5MechanismSummary == null
                        ? null : configuredEtaV5MechanismSummary.trim();
        if (etaV5ThresholdOracle && !etaV4) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 Phase-A oracle requires packstar.pac.etaV4=true"
                            + " so it can reuse the frozen OAIS proposal path");
        }
        if (etaV5ThresholdOracle && !conditionalRepairAudit) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 Phase-A oracle requires conditionalRepairAudit=true");
        }
        if (etaV5ThresholdOracle
                && (etaV5OracleOutputTsv == null
                || etaV5OracleOutputTsv.isEmpty())) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 Phase-A oracle requires an output TSV path");
        }
        if (etaV5ThresholdOracle
                && (!new File(etaV5OracleOutputTsv).isAbsolute()
                || !etaV5OracleOutputTsv.startsWith(
                "/usr/xtmp/lz280/"))) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 oracle output must be an absolute path under"
                            + " /usr/xtmp/lz280: " + etaV5OracleOutputTsv);
        }
        if (etaV5MaxEnumeratedStates < 1L) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 max enumerated states must be positive: "
                            + etaV5MaxEnumeratedStates);
        }
        if (etaV5RebaseInterval < 1) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 rebase interval must be positive: "
                            + etaV5RebaseInterval);
        }
        if (!Double.isFinite(etaV5LogZTolerance)
                || !(etaV5LogZTolerance > 0.0)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 logZ tolerance must be finite and positive: "
                            + etaV5LogZTolerance);
        }
        if (etaV5MechanismPilot && !etaV5ThresholdOracle) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 mechanism pilot requires the exact"
                            + " threshold oracle integrity gate");
        }
        if (etaV5MechanismPilot
                && etaV5MechanismSamples
                != DEFAULT_ETA_V5_MECHANISM_SAMPLES) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 frozen mechanism pilot requires exactly "
                            + DEFAULT_ETA_V5_MECHANISM_SAMPLES
                            + " p_m samples; got "
                            + etaV5MechanismSamples);
        }
        if (etaV5MechanismPilot
                && !isEtaV5XtmpOutputPath(
                etaV5MechanismAssignmentsTsv)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 mechanism assignment output must be an"
                            + " absolute path under /usr/xtmp/lz280: "
                            + etaV5MechanismAssignmentsTsv);
        }
        if (etaV5MechanismPilot
                && !isEtaV5XtmpOutputPath(
                etaV5MechanismSummaryTsv)) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 mechanism summary output must be an"
                            + " absolute path under /usr/xtmp/lz280: "
                            + etaV5MechanismSummaryTsv);
        }
        if (etaV5MechanismPilot
                && (etaV5MechanismAssignmentsTsv.equals(
                etaV5MechanismSummaryTsv)
                || etaV5MechanismAssignmentsTsv.equals(
                etaV5OracleOutputTsv)
                || etaV5MechanismSummaryTsv.equals(
                etaV5OracleOutputTsv))) {
            throw new IllegalArgumentException(
                    "PACK* eta-v5 oracle, mechanism assignment, and"
                            + " mechanism summary outputs must be distinct");
        }
        this.unconditionalTail = getConfigBoolean(
                PAC_UNCONDITIONAL_TAIL_PROPERTY, DEFAULT_UNCONDITIONAL_TAIL);
        if (conditionalRepairAudit && unconditionalTail) {
            throw new IllegalArgumentException(
                    "PACK* conditionalRepairAudit and unconditionalTail modes are mutually exclusive");
        }
        this.fixedBulkSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_BULK_SAMPLES_PROPERTY,
                        DEFAULT_BULK_SAMPLES)), 2);
        this.fixedTailSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_TAIL_SAMPLES_PROPERTY,
                        DEFAULT_TAIL_SAMPLES)), 2);
        this.tailTiltLambdas = parseTiltLambdas(getConfigProperty(
                PAC_TAIL_TILT_LAMBDAS_PROPERTY, DEFAULT_TAIL_TILT_LAMBDAS));
        this.etaRepairEnabled = getConfigBoolean(
                PAC_ETA_REPAIR_ENABLED_PROPERTY, DEFAULT_ETA_REPAIR_ENABLED);
        this.etaRepairMaxRounds = Math.max(0, getConfigInteger(
                PAC_ETA_REPAIR_MAX_ROUNDS_PROPERTY, DEFAULT_ETA_REPAIR_MAX_ROUNDS));
        this.etaRepairTriggerKcal = Math.max(0.0, getConfigDouble(
                PAC_ETA_REPAIR_TRIGGER_PROPERTY, DEFAULT_ETA_REPAIR_TRIGGER_KCAL));
        this.etaRepairTargetKcal = Math.max(0.0, getConfigDouble(
                PAC_ETA_REPAIR_TARGET_PROPERTY, DEFAULT_ETA_REPAIR_TARGET_KCAL));
        this.etaRepairSafety = clamp(getConfigDouble(
                PAC_ETA_REPAIR_SAFETY_PROPERTY, DEFAULT_ETA_REPAIR_SAFETY), 0.0, 1.0);
        this.etaRepairMinAlpha = clamp(getConfigDouble(
                PAC_ETA_REPAIR_MIN_ALPHA_PROPERTY, DEFAULT_ETA_REPAIR_MIN_ALPHA), 0.0, 1.0);
        this.etaEnabled = getConfigBoolean(PAC_ETA_ENABLED_PROPERTY, DEFAULT_ETA_ENABLED);
        if (frequencySeverityProduction && !conditionalRepairAudit) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity production requires conditionalRepairAudit=true");
        }
        if (frequencySeverityProduction && !etaEnabled) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity production requires etaEnabled=true");
        }
        if (frequencySeverityProduction && (dpRich || etaV3 || etaV4
                || etaV5ThresholdOracle || unconditionalTail)) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity production is mutually exclusive with legacy"
                            + " DP-rich/V3/V4/V5 and unconditional-tail modes");
        }
        if (frequencySeverityProduction && trainSamples < 4) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity two-fold selection requires at least four training samples");
        }
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

    private int sanitizeSampleBudget(int requestedSampleBudget) {
        return requestedSampleBudget > 0 ? requestedSampleBudget : Integer.MAX_VALUE;
    }

    private int capSampleCount(int configured, int minimum) {
        int capped = sampleBudget == Integer.MAX_VALUE
                ? configured
                : Math.min(configured, Math.max(minimum, sampleBudget));
        return Math.max(minimum, capped);
    }

    public void setSampleListener(PackStarSampleListener sampleListener) {
        this.sampleListener = sampleListener;
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

    static final class FrequencySeverityShrinkPair {
        final double unary;
        final double pair;

        FrequencySeverityShrinkPair(double unary, double pair) {
            this.unary = unary;
            this.pair = pair;
        }

        String id() {
            return String.format(Locale.ROOT, "ku-%g-kp-%g", unary, pair);
        }
    }

    private static class ConditionalDpRichDesign {
        final EtaCorrections eta;
        final EnergyMatrix emat;
        final double alpha;
        final double selectedBoundKcal;
        final double adaptationMaxOverCorrectionKcal;
        final double logZ;
        final double logRangeProxy;
        final String source;
        final int candidateDpSweeps;
        final double gaugeShiftDownKcal;
        final int predictedFinalSamples;
        final double predictedEpsilonAtMax;
        final double chiSquareProxy;
        final double screeningEss;

        ConditionalDpRichDesign(EtaCorrections eta, EnergyMatrix emat,
                                double alpha, double selectedBoundKcal,
                                double adaptationMaxOverCorrectionKcal,
                                double logZ, double logRangeProxy,
                                String source, int candidateDpSweeps) {
            this.eta = eta;
            this.emat = emat;
            this.alpha = alpha;
            this.selectedBoundKcal = selectedBoundKcal;
            this.adaptationMaxOverCorrectionKcal =
                    adaptationMaxOverCorrectionKcal;
            this.logZ = logZ;
            this.logRangeProxy = logRangeProxy;
            this.source = source;
            this.candidateDpSweeps = candidateDpSweeps;
            this.gaugeShiftDownKcal = 0.0;
            this.predictedFinalSamples = -1;
            this.predictedEpsilonAtMax = Double.NaN;
            this.chiSquareProxy = Double.NaN;
            this.screeningEss = Double.NaN;
        }

        ConditionalDpRichDesign(EtaCorrections eta, EnergyMatrix emat,
                                double alpha, double selectedBoundKcal,
                                double adaptationMaxOverCorrectionKcal,
                                double logZ, double logRangeProxy,
                                String source, int candidateDpSweeps,
                                double gaugeShiftDownKcal,
                                int predictedFinalSamples,
                                double predictedEpsilonAtMax,
                                double chiSquareProxy,
                                double screeningEss) {
            this.eta = eta;
            this.emat = emat;
            this.alpha = alpha;
            this.selectedBoundKcal = selectedBoundKcal;
            this.adaptationMaxOverCorrectionKcal =
                    adaptationMaxOverCorrectionKcal;
            this.logZ = logZ;
            this.logRangeProxy = logRangeProxy;
            this.source = source;
            this.candidateDpSweeps = candidateDpSweeps;
            this.gaugeShiftDownKcal = gaugeShiftDownKcal;
            this.predictedFinalSamples = predictedFinalSamples;
            this.predictedEpsilonAtMax = predictedEpsilonAtMax;
            this.chiSquareProxy = chiSquareProxy;
            this.screeningEss = screeningEss;
        }
    }

    private static class ConditionalProposalSnapshot {
        final String sourceProposalId;
        final EtaCorrections eta;
        final EnergyMatrix emat;
        final double logZ;
        final double reportedAlpha;

        ConditionalProposalSnapshot(String sourceProposalId,
                                    EtaCorrections eta,
                                    EnergyMatrix emat,
                                    double logZ,
                                    double reportedAlpha) {
            this.sourceProposalId = sourceProposalId;
            this.eta = eta;
            this.emat = emat;
            this.logZ = logZ;
            this.reportedAlpha = reportedAlpha;
        }
    }

    private static class ConditionalEtaV4Anchor {
        final ConditionalProposalSnapshot proposal;
        final String sourceProposalId;
        final List<CCDResult> onPolicySamples;
        final double adaptationMaxOverCorrectionKcal;
        final double chiSquareProxy;
        final double screeningEss;
        final double epsilonAtMax;
        final int predictedFinalSamples;
        final int candidateDpSweeps;

        ConditionalEtaV4Anchor(
                ConditionalProposalSnapshot proposal,
                String sourceProposalId,
                List<CCDResult> onPolicySamples,
                double adaptationMaxOverCorrectionKcal,
                double chiSquareProxy,
                double screeningEss,
                double epsilonAtMax,
                int predictedFinalSamples,
                int candidateDpSweeps) {
            this.proposal = proposal;
            this.sourceProposalId = sourceProposalId;
            this.onPolicySamples = onPolicySamples;
            this.adaptationMaxOverCorrectionKcal =
                    adaptationMaxOverCorrectionKcal;
            this.chiSquareProxy = chiSquareProxy;
            this.screeningEss = screeningEss;
            this.epsilonAtMax = epsilonAtMax;
            this.predictedFinalSamples = predictedFinalSamples;
            this.candidateDpSweeps = candidateDpSweeps;
        }
    }

    private static class ConditionalEtaV4MixtureDesign {
        final ConditionalProposalSnapshot[] components;
        final double[] weights;
        final int anchorIndex;
        final double anchorBoundKcal;
        final double effectiveBoundKcal;
        final double anchorLogRange;
        final boolean componentCover;
        final ConditionalAdaptationStats adaptation;
        final double clipLogCap;
        final ConditionalSizing sizing;
        final int candidateDpSweeps;

        ConditionalEtaV4MixtureDesign(
                ConditionalProposalSnapshot[] components,
                double[] weights, int anchorIndex,
                double anchorBoundKcal, double effectiveBoundKcal,
                double anchorLogRange, boolean componentCover,
                ConditionalAdaptationStats adaptation,
                double clipLogCap, ConditionalSizing sizing,
                int candidateDpSweeps) {
            this.components = components;
            this.weights = weights;
            this.anchorIndex = anchorIndex;
            this.anchorBoundKcal = anchorBoundKcal;
            this.effectiveBoundKcal = effectiveBoundKcal;
            this.anchorLogRange = anchorLogRange;
            this.componentCover = componentCover;
            this.adaptation = adaptation;
            this.clipLogCap = clipLogCap;
            this.sizing = sizing;
            this.candidateDpSweeps = candidateDpSweeps;
        }
    }

    private static class ConditionalEtaV4Tail {
        final ConditionalProposalSnapshot proposal;
        final List<CCDResult> baseScreeningSamples;
        final double tailPower;
        final double trustKcal;
        final double gradientTailEss;
        final double offPolicyProposalEss;
        final double logSecondMomentImprovement;
        final int candidateDpSweeps;
        final int roundIndex;

        ConditionalEtaV4Tail(ConditionalProposalSnapshot proposal,
                             List<CCDResult> baseScreeningSamples,
                             double tailPower, double trustKcal,
                             double gradientTailEss,
                             double offPolicyProposalEss,
                             double logSecondMomentImprovement,
                             int candidateDpSweeps,
                             int roundIndex) {
            this.proposal = proposal;
            this.baseScreeningSamples = baseScreeningSamples;
            this.tailPower = tailPower;
            this.trustKcal = trustKcal;
            this.gradientTailEss = gradientTailEss;
            this.offPolicyProposalEss = offPolicyProposalEss;
            this.logSecondMomentImprovement =
                    logSecondMomentImprovement;
            this.candidateDpSweeps = candidateDpSweeps;
            this.roundIndex = roundIndex;
        }
    }

    private static class EtaV4FeatureCell {
        final boolean pair;
        final int pos1;
        final int rc1;
        final int pos2;
        final int rc2;

        EtaV4FeatureCell(int pos1, int rc1) {
            this.pair = false;
            this.pos1 = pos1;
            this.rc1 = rc1;
            this.pos2 = -1;
            this.rc2 = -1;
        }

        EtaV4FeatureCell(int pos1, int rc1, int pos2, int rc2) {
            this.pair = true;
            this.pos1 = pos1;
            this.rc1 = rc1;
            this.pos2 = pos2;
            this.rc2 = rc2;
        }
    }

    private static class EtaV4SparseFeatures {
        final EtaV4FeatureCell[] cells;
        final int[][] sampleIndices;

        EtaV4SparseFeatures(EtaV4FeatureCell[] cells,
                            int[][] sampleIndices) {
            this.cells = cells;
            this.sampleIndices = sampleIndices;
        }
    }

    private static class EtaV4ProposalUpdate {
        final EtaCorrections eta;
        final double coefficientScaleKcal;
        final double maxObservedEnergyDeltaKcal;
        final double maxCellDeltaKcal;

        EtaV4ProposalUpdate(EtaCorrections eta,
                            double coefficientScaleKcal,
                            double maxObservedEnergyDeltaKcal,
                            double maxCellDeltaKcal) {
            this.eta = eta;
            this.coefficientScaleKcal = coefficientScaleKcal;
            this.maxObservedEnergyDeltaKcal =
                    maxObservedEnergyDeltaKcal;
            this.maxCellDeltaKcal = maxCellDeltaKcal;
        }
    }

    private static class EtaV4OffPolicyScore {
        final double logSecondMoment;
        final double proposalEffectiveSampleSize;
        final double minLogProposalRatio;
        final double maxLogProposalRatio;

        EtaV4OffPolicyScore(double logSecondMoment,
                            double proposalEffectiveSampleSize,
                            double minLogProposalRatio,
                            double maxLogProposalRatio) {
            this.logSecondMoment = logSecondMoment;
            this.proposalEffectiveSampleSize =
                    proposalEffectiveSampleSize;
            this.minLogProposalRatio = minLogProposalRatio;
            this.maxLogProposalRatio = maxLogProposalRatio;
        }
    }

    static class ConditionalImportanceWeights {
        final double[] normalized;
        final double effectiveSampleSize;
        final double minLogWeight;
        final double maxLogWeight;

        ConditionalImportanceWeights(double[] normalized,
                                     double effectiveSampleSize,
                                     double minLogWeight,
                                     double maxLogWeight) {
            this.normalized = normalized;
            this.effectiveSampleSize = effectiveSampleSize;
            this.minLogWeight = minLogWeight;
            this.maxLogWeight = maxLogWeight;
        }
    }

    private static class ConditionalAdaptationStats {
        final double[] logResidualWeights;
        final double[] normalizedProposalWeights;
        final double effectiveSampleSize;
        final double minLogProposalWeight;
        final double maxLogProposalWeight;
        final int sourceProposalCount;

        ConditionalAdaptationStats(double[] logResidualWeights,
                                   double[] normalizedProposalWeights,
                                   double effectiveSampleSize,
                                   double minLogProposalWeight,
                                   double maxLogProposalWeight,
                                   int sourceProposalCount) {
            this.logResidualWeights = logResidualWeights;
            this.normalizedProposalWeights = normalizedProposalWeights;
            this.effectiveSampleSize = effectiveSampleSize;
            this.minLogProposalWeight = minLogProposalWeight;
            this.maxLogProposalWeight = maxLogProposalWeight;
            this.sourceProposalCount = sourceProposalCount;
        }
    }

    private static class ConditionalSizing {
        final double meanClipped;
        final double varianceClipped;
        final double meanExcess;
        final double varianceExcess;
        final double excessRange;
        final double clippedProbability;
        final int finalSamples;
        final double epsilonAtMaxSamples;

        ConditionalSizing(double meanClipped, double varianceClipped,
                          double meanExcess, double varianceExcess,
                          double excessRange, double clippedProbability,
                          int finalSamples, double epsilonAtMaxSamples) {
            this.meanClipped = meanClipped;
            this.varianceClipped = varianceClipped;
            this.meanExcess = meanExcess;
            this.varianceExcess = varianceExcess;
            this.excessRange = excessRange;
            this.clippedProbability = clippedProbability;
            this.finalSamples = finalSamples;
            this.epsilonAtMaxSamples = epsilonAtMaxSamples;
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

    private static int getConfigInteger(String key, int defaultValue) {
        return PackStarConfig.getInteger(key, defaultValue, "[PACK*]");
    }

    private static long getConfigLong(String key, long defaultValue) {
        return PackStarConfig.getLong(key, defaultValue, "[PACK*]");
    }

    private static double getConfigDouble(String key, double defaultValue) {
        return PackStarConfig.getDouble(key, defaultValue, "[PACK*]");
    }

    private static long getConfigBytes(String key, long defaultValue) {
        return PackStarConfig.getBytes(key, defaultValue, "[PACK*]");
    }

    private static boolean getConfigBoolean(String key, boolean defaultValue) {
        return PackStarConfig.getBoolean(key, defaultValue);
    }

    private static String getConfigProperty(String key, String defaultValue) {
        return PackStarConfig.getProperty(key, defaultValue);
    }

    private static double clamp(double value, double lower, double upper) {
        if (Double.isNaN(value)) return lower;
        return Math.max(lower, Math.min(upper, value));
    }

    private static double[] parseTiltLambdas(String configured) {
        TreeSet<Double> values = new TreeSet<>();
        values.add(0.0);
        if (configured != null) {
            for (String token : configured.split(",")) {
                String trimmed = token.trim();
                if (trimmed.isEmpty()) continue;
                double value;
                try {
                    value = Double.parseDouble(trimmed);
                } catch (NumberFormatException ex) {
                    throw new IllegalArgumentException(
                            "Invalid PACK* tail tilt lambda: " + trimmed, ex);
                }
                if (!Double.isFinite(value) || value < 0.0) {
                    throw new IllegalArgumentException(
                            "PACK* tail tilt lambdas must be finite and nonnegative: " + trimmed);
                }
                values.add(value == 0.0 ? 0.0 : value);
            }
        }
        double[] parsed = new double[values.size()];
        int index = 0;
        for (double value : values) parsed[index++] = value;
        return parsed;
    }

    static double[] parsePositiveGrid(String configured, String label) {
        TreeSet<Double> values = new TreeSet<>();
        if (configured != null) {
            for (String token : configured.split(",")) {
                String trimmed = token.trim();
                if (trimmed.isEmpty()) continue;
                double value;
                try {
                    value = Double.parseDouble(trimmed);
                } catch (NumberFormatException ex) {
                    throw new IllegalArgumentException(
                            "Invalid " + label + " value: " + trimmed, ex);
                }
                if (!Double.isFinite(value) || !(value > 0.0)) {
                    throw new IllegalArgumentException(
                            label + " values must be finite and positive: " + trimmed);
                }
                values.add(value);
            }
        }
        if (values.isEmpty()) {
            throw new IllegalArgumentException(label + " must contain at least one value");
        }
        double[] parsed = new double[values.size()];
        int index = 0;
        for (double value : values) parsed[index++] = value;
        return parsed;
    }

    static double[] parseFixedOneGrid(String configured, String label) {
        if (configured == null || configured.trim().isEmpty()) {
            throw new IllegalArgumentException(label + " must be fixed to 1");
        }
        String[] tokens = configured.split(",", -1);
        if (tokens.length != 1) {
            throw new IllegalArgumentException(
                    label + " is not tunable and must contain only 1: "
                            + configured);
        }
        final double value;
        try {
            value = Double.parseDouble(tokens[0].trim());
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException(
                    "Invalid " + label + " value: " + tokens[0].trim(), ex);
        }
        if (!Double.isFinite(value) || value != 1.0) {
            throw new IllegalArgumentException(
                    label + " is not tunable and must be fixed to 1: "
                            + configured);
        }
        return new double[]{1.0};
    }

    static FrequencySeverityShrinkPair[] parseFrequencySeverityShrinkGrid(String configured) {
        if (configured == null || configured.trim().isEmpty()) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity shrink grid must not be empty");
        }
        List<FrequencySeverityShrinkPair> parsed = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        boolean containsRaw = false;
        for (String token : configured.split(",")) {
            String trimmed = token.trim();
            String[] parts = trimmed.split(":", -1);
            if (parts.length != 2) {
                throw new IllegalArgumentException(
                        "PACK* frequency/severity shrink values must have unary:pair form: "
                                + trimmed);
            }
            final double unary;
            final double pair;
            try {
                unary = Double.parseDouble(parts[0].trim());
                pair = Double.parseDouble(parts[1].trim());
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException(
                        "Invalid PACK* frequency/severity shrink value: " + trimmed, ex);
            }
            if (!Double.isFinite(unary) || unary < 0.0
                    || !Double.isFinite(pair) || pair < 0.0) {
                throw new IllegalArgumentException(
                        "PACK* frequency/severity shrink strengths must be finite and nonnegative: "
                                + trimmed);
            }
            String key = Long.toHexString(Double.doubleToLongBits(unary))
                    + ":" + Long.toHexString(Double.doubleToLongBits(pair));
            if (!seen.add(key)) continue;
            parsed.add(new FrequencySeverityShrinkPair(unary, pair));
            if (unary == 0.0 && pair == 0.0) containsRaw = true;
        }
        if (parsed.isEmpty() || !containsRaw) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity shrink grid must contain raw eta as 0:0");
        }
        return parsed.toArray(new FrequencySeverityShrinkPair[0]);
    }

    static double[] parseFiniteIncreasingGrid(
            String configured, String label) {
        if (configured == null || configured.trim().isEmpty()) {
            throw new IllegalArgumentException(
                    label + " must contain at least one value");
        }
        String[] tokens = configured.split(",");
        double[] parsed = new double[tokens.length];
        int count = 0;
        for (String token : tokens) {
            String trimmed = token.trim();
            if (trimmed.isEmpty()) {
                throw new IllegalArgumentException(
                        label + " contains an empty value: " + configured);
            }
            final double value;
            try {
                value = Double.parseDouble(trimmed);
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException(
                        "Invalid " + label + " value: " + trimmed, ex);
            }
            if (!Double.isFinite(value)) {
                throw new IllegalArgumentException(
                        label + " values must be finite: " + trimmed);
            }
            if (count > 0 && !(value > parsed[count - 1])) {
                throw new IllegalArgumentException(
                        label + " values must be strictly increasing: "
                                + configured);
            }
            parsed[count++] = value;
        }
        return count == parsed.length
                ? parsed : Arrays.copyOf(parsed, count);
    }

    static double selectSmallestResidualBound(double observedMaxKcal,
                                              double[] sortedGridKcal) {
        if (!Double.isFinite(observedMaxKcal)
                || sortedGridKcal == null || sortedGridKcal.length == 0) {
            return Double.NaN;
        }
        for (double candidate : sortedGridKcal) {
            if (!Double.isFinite(candidate) || !(candidate > 0.0)) {
                return Double.NaN;
            }
            if (observedMaxKcal <= candidate) return candidate;
        }
        return Double.NaN;
    }

    static boolean gridContains(double[] grid, double value) {
        if (grid == null || !Double.isFinite(value)) return false;
        for (double candidate : grid) {
            if (Double.isFinite(candidate)
                    && Math.abs(candidate - value) <= 1.0e-12) {
                return true;
            }
        }
        return false;
    }

    static double conditionalClippedNormalized(double logWeight,
                                               double logCap) {
        if (!Double.isFinite(logWeight) || !Double.isFinite(logCap)) {
            return Double.NaN;
        }
        return Math.exp(Math.min(logWeight, logCap) - logCap);
    }

    static double conditionalExcessNormalized(double logWeight,
                                              double logCap) {
        if (!Double.isFinite(logWeight) || !Double.isFinite(logCap)) {
            return Double.NaN;
        }
        return Math.max(0.0, Math.expm1(logWeight - logCap));
    }

    static double conditionalExcessRange(double boundKcal, double rt,
                                         double logCap) {
        if (!Double.isFinite(boundKcal) || boundKcal < 0.0
                || !Double.isFinite(rt) || !(rt > 0.0)
                || !Double.isFinite(logCap)) {
            return Double.NaN;
        }
        return Math.max(0.0, Math.expm1(boundKcal / rt - logCap));
    }

    static double conditionalLogRangeProxy(double logZ, double boundKcal,
                                           double rt) {
        if (!Double.isFinite(logZ) || !Double.isFinite(boundKcal)
                || boundKcal < 0.0 || !Double.isFinite(rt) || !(rt > 0.0)) {
            return Double.NaN;
        }
        return logZ + boundKcal / rt;
    }

    static boolean isBetterConditionalDpRichCandidate(
            double logRange, double boundKcal, double alpha,
            double bestLogRange, double bestBoundKcal, double bestAlpha) {
        if (!Double.isFinite(logRange) || !Double.isFinite(boundKcal)
                || !Double.isFinite(alpha)) return false;
        if (!Double.isFinite(bestLogRange)) return true;
        double tolerance = 1.0e-12;
        if (logRange < bestLogRange - tolerance) return true;
        if (Math.abs(logRange - bestLogRange) > tolerance) return false;
        if (boundKcal < bestBoundKcal - tolerance) return true;
        if (Math.abs(boundKcal - bestBoundKcal) > tolerance) return false;
        return alpha > bestAlpha + tolerance;
    }

    /**
     * Certificate-aware v3 ordering.  Ineligible candidates never win.  Among
     * candidates predicted to reach the requested epsilon at the frozen sample
     * cap, prefer lower sample cost, then the empirical chi-square proxy, then
     * the conservative absolute-range proxy.
     */
    static boolean isBetterConditionalEtaV3Candidate(
            boolean eligible, int finalSamples, double chiSquareProxy,
            double logRange, double epsilonAtMax,
            int bestFinalSamples, double bestChiSquareProxy,
            double bestLogRange, double bestEpsilonAtMax) {
        if (!eligible || finalSamples < 2
                || !Double.isFinite(chiSquareProxy)
                || chiSquareProxy < 0.0
                || !Double.isFinite(logRange)
                || !Double.isFinite(epsilonAtMax)) return false;
        if (bestFinalSamples < 2) return true;
        if (finalSamples != bestFinalSamples) {
            return finalSamples < bestFinalSamples;
        }
        double tolerance = 1.0e-12;
        if (chiSquareProxy < bestChiSquareProxy - tolerance) return true;
        if (Math.abs(chiSquareProxy - bestChiSquareProxy) > tolerance) {
            return false;
        }
        if (logRange < bestLogRange - tolerance) return true;
        if (Math.abs(logRange - bestLogRange) > tolerance) return false;
        return epsilonAtMax < bestEpsilonAtMax - tolerance;
    }

    static double conditionalImportanceLogRatio(
            double sourceEnergy, double targetEnergy,
            double sourceLogZ, double targetLogZ, double rt) {
        if (!Double.isFinite(sourceEnergy) || !Double.isFinite(targetEnergy)
                || !Double.isFinite(sourceLogZ)
                || !Double.isFinite(targetLogZ)
                || !Double.isFinite(rt) || !(rt > 0.0)) {
            return Double.NaN;
        }
        double logRatio = (sourceEnergy - targetEnergy) / rt
                + sourceLogZ - targetLogZ;
        return Double.isFinite(logRatio) ? logRatio : Double.NaN;
    }

    static ConditionalImportanceWeights normalizeConditionalImportanceWeights(
            double[] logWeights) {
        if (logWeights == null || logWeights.length == 0) return null;
        double minLogWeight = Double.POSITIVE_INFINITY;
        double maxLogWeight = Double.NEGATIVE_INFINITY;
        for (double logWeight : logWeights) {
            if (!Double.isFinite(logWeight)) return null;
            minLogWeight = Math.min(minLogWeight, logWeight);
            maxLogWeight = Math.max(maxLogWeight, logWeight);
        }

        double[] normalized = new double[logWeights.length];
        double sum = 0.0;
        for (int i = 0; i < logWeights.length; i++) {
            normalized[i] = Math.exp(logWeights[i] - maxLogWeight);
            sum += normalized[i];
        }
        if (!Double.isFinite(sum) || !(sum > 0.0)) return null;

        double sumSquares = 0.0;
        for (int i = 0; i < normalized.length; i++) {
            normalized[i] /= sum;
            sumSquares += normalized[i] * normalized[i];
        }
        if (!Double.isFinite(sumSquares) || !(sumSquares > 0.0)) return null;
        double ess = 1.0 / sumSquares;
        if (!Double.isFinite(ess)) return null;
        return new ConditionalImportanceWeights(
                normalized, ess, minLogWeight, maxLogWeight);
    }

    static boolean conditionalDpRichEssPasses(
            double effectiveSampleSize, int sampleCount,
            double minimumFraction) {
        if (!Double.isFinite(effectiveSampleSize) || sampleCount < 2
                || !Double.isFinite(minimumFraction)
                || !(minimumFraction > 0.0) || minimumFraction > 1.0) {
            return false;
        }
        double threshold = Math.max(2.0, minimumFraction * sampleCount);
        return effectiveSampleSize + 1.0e-12 >= threshold;
    }

    static double conditionalWeightedMean(double[] values, double[] weights) {
        if (values == null || weights == null
                || values.length == 0 || values.length != weights.length) {
            return Double.NaN;
        }
        double sumWeights = 0.0;
        double sum = 0.0;
        for (int i = 0; i < values.length; i++) {
            if (!Double.isFinite(values[i]) || !Double.isFinite(weights[i])
                    || weights[i] < 0.0) return Double.NaN;
            sumWeights += weights[i];
            sum += weights[i] * values[i];
        }
        if (!Double.isFinite(sumWeights) || !(sumWeights > 0.0)
                || !Double.isFinite(sum)) return Double.NaN;
        return sum / sumWeights;
    }

    /**
     * Empirical CV^2 of the unnormalised residual weights under the target
     * proposal.  For exact on-policy expectations this equals
     * chi^2(p_true || q_eta), and is invariant to a global energy gauge.
     */
    static double conditionalWeightChiSquareProxy(
            double[] logResidualWeights, double[] proposalWeights) {
        if (logResidualWeights == null || proposalWeights == null
                || logResidualWeights.length == 0
                || logResidualWeights.length != proposalWeights.length) {
            return Double.NaN;
        }
        double maxLogWeight = Double.NEGATIVE_INFINITY;
        for (double logWeight : logResidualWeights) {
            if (!Double.isFinite(logWeight)) return Double.NaN;
            maxLogWeight = Math.max(maxLogWeight, logWeight);
        }
        double[] scaled = new double[logResidualWeights.length];
        double[] scaledSquared = new double[logResidualWeights.length];
        for (int i = 0; i < logResidualWeights.length; i++) {
            scaled[i] = Math.exp(logResidualWeights[i] - maxLogWeight);
            scaledSquared[i] = scaled[i] * scaled[i];
        }
        double mean = conditionalWeightedMean(scaled, proposalWeights);
        double secondMoment = conditionalWeightedMean(
                scaledSquared, proposalWeights);
        if (!Double.isFinite(mean) || !(mean > 0.0)
                || !Double.isFinite(secondMoment)) return Double.NaN;
        double cvSquared = secondMoment / (mean * mean) - 1.0;
        if (!Double.isFinite(cvSquared)) return Double.NaN;
        return Math.max(0.0, cvSquared);
    }

    static double logMeanExp(double[] logValues) {
        if (logValues == null || logValues.length == 0) return Double.NaN;
        double max = Double.NEGATIVE_INFINITY;
        for (double value : logValues) {
            if (!Double.isFinite(value)) return Double.NaN;
            max = Math.max(max, value);
        }
        double sum = 0.0;
        for (double value : logValues) sum += Math.exp(value - max);
        if (!Double.isFinite(sum) || !(sum > 0.0)) return Double.NaN;
        return max + Math.log(sum) - Math.log(logValues.length);
    }

    /** Estimate a candidate's true IS second moment from a source-q batch. */
    private EtaV4OffPolicyScore scoreEtaV4CandidateOffPolicy(
            List<CCDResult> sourceSamples,
            EnergyMatrix sourceEmat, double sourceLogZ,
            EnergyMatrix candidateEmat, double candidateLogZ) {
        if (sourceSamples == null || sourceSamples.size() < 2
                || sourceEmat == null || candidateEmat == null
                || !Double.isFinite(sourceLogZ)
                || !Double.isFinite(candidateLogZ)) return null;
        double[] logSecondMomentTerms =
                new double[sourceSamples.size()];
        double[] logProposalRatios = new double[sourceSamples.size()];
        for (int i = 0; i < sourceSamples.size(); i++) {
            CCDResult sample = sourceSamples.get(i);
            if (!Double.isFinite(sample.eTrue)) return null;
            double sourceEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, sourceEmat);
            double candidateEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, candidateEmat);
            if (!Double.isFinite(sourceEnergy)
                    || !Double.isFinite(candidateEnergy)) return null;
            double logSourceQ = -sourceEnergy / RT - sourceLogZ;
            double logCandidateQ =
                    -candidateEnergy / RT - candidateLogZ;
            double logTarget = -sample.eTrue / RT;
            logSecondMomentTerms[i] = 2.0 * logTarget
                    - logCandidateQ - logSourceQ;
            logProposalRatios[i] = logCandidateQ - logSourceQ;
            if (!Double.isFinite(logSecondMomentTerms[i])
                    || !Double.isFinite(logProposalRatios[i])) return null;
        }
        ConditionalImportanceWeights normalized =
                normalizeConditionalImportanceWeights(logProposalRatios);
        double logSecondMoment = logMeanExp(logSecondMomentTerms);
        if (normalized == null || !Double.isFinite(logSecondMoment)) {
            return null;
        }
        return new EtaV4OffPolicyScore(logSecondMoment,
                normalized.effectiveSampleSize,
                normalized.minLogWeight, normalized.maxLogWeight);
    }

    static double conditionalWeightedSampleVariance(
            double[] values, double[] weights, double mean) {
        if (values == null || weights == null
                || values.length < 2 || values.length != weights.length
                || !Double.isFinite(mean)) return Double.NaN;
        double sumWeights = 0.0;
        double sumWeightSquares = 0.0;
        double sumSquaredDeviations = 0.0;
        for (int i = 0; i < values.length; i++) {
            double value = values[i];
            double weight = weights[i];
            if (!Double.isFinite(value) || !Double.isFinite(weight)
                    || weight < 0.0) return Double.NaN;
            double deviation = value - mean;
            sumWeights += weight;
            sumWeightSquares += weight * weight;
            sumSquaredDeviations += weight * deviation * deviation;
        }
        if (!Double.isFinite(sumWeights) || !(sumWeights > 0.0)
                || !Double.isFinite(sumWeightSquares)
                || !Double.isFinite(sumSquaredDeviations)) return Double.NaN;
        double denominator = sumWeights
                - sumWeightSquares / sumWeights;
        if (!(denominator > 0.0) || !Double.isFinite(denominator)) {
            return Double.NaN;
        }
        return Math.max(0.0, sumSquaredDeviations / denominator);
    }

    /**
     * Linear weighted quantile using weighted sample midpoints.  Equal weights
     * reduce exactly to the legacy order-statistic interpolation at i/(n-1).
     */
    static double conditionalWeightedQuantile(
            double[] values, double[] weights, double quantile) {
        if (values == null || weights == null
                || values.length == 0 || values.length != weights.length
                || !Double.isFinite(quantile)) return Double.NaN;
        int positive = 0;
        double totalWeight = 0.0;
        for (int i = 0; i < values.length; i++) {
            if (!Double.isFinite(values[i]) || !Double.isFinite(weights[i])
                    || weights[i] < 0.0) return Double.NaN;
            if (weights[i] > 0.0) positive++;
            totalWeight += weights[i];
        }
        if (positive == 0 || !Double.isFinite(totalWeight)
                || !(totalWeight > 0.0)) return Double.NaN;

        Integer[] order = new Integer[positive];
        int next = 0;
        for (int i = 0; i < weights.length; i++) {
            if (weights[i] > 0.0) order[next++] = i;
        }
        Arrays.sort(order, (left, right) -> {
            int comparison = Double.compare(values[left], values[right]);
            return comparison != 0 ? comparison : Integer.compare(left, right);
        });
        if (order.length == 1) return values[order[0]];

        double firstMidpoint = 0.5 * weights[order[0]] / totalWeight;
        double lastMidpoint = 1.0
                - 0.5 * weights[order[order.length - 1]] / totalWeight;
        if (!(lastMidpoint > firstMidpoint)) return values[order[0]];
        double target = clamp(quantile, 0.0, 1.0);
        double cumulative = 0.0;
        double previousPosition = 0.0;
        double previousValue = values[order[0]];
        for (int rank = 0; rank < order.length; rank++) {
            int index = order[rank];
            cumulative += weights[index];
            double midpoint = (cumulative - 0.5 * weights[index])
                    / totalWeight;
            double position = clamp(
                    (midpoint - firstMidpoint)
                            / (lastMidpoint - firstMidpoint),
                    0.0, 1.0);
            double value = values[index];
            if (target <= position) {
                if (rank == 0 || !(position > previousPosition)) {
                    return value;
                }
                double fraction = (target - previousPosition)
                        / (position - previousPosition);
                return previousValue * (1.0 - fraction) + value * fraction;
            }
            previousPosition = position;
            previousValue = value;
        }
        return values[order[order.length - 1]];
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

    static long deriveRandomSeed(long baseRandomSeed, String randomStreamIdentity) {
        String identity = requireRandomStreamIdentity(randomStreamIdentity);
        long identityHash = 0xcbf29ce484222325L;
        for (byte value : identity.getBytes(StandardCharsets.UTF_8)) {
            identityHash ^= value & 0xffL;
            identityHash *= 0x100000001b3L;
        }
        return mix64(baseRandomSeed ^ identityHash);
    }

    private Random stageRandom(String stage) {
        return new Random(deriveRandomSeed(
                randomSeed, randomStreamIdentity + "|unconditional-v1|" + stage));
    }

    private static String requireRandomStreamIdentity(String identity) {
        if (identity == null || identity.trim().isEmpty()) {
            throw new IllegalArgumentException(
                    "PACK* random-stream identity is required");
        }
        return identity;
    }

    private static String defaultRandomStreamIdentity(RCs rcs) {
        StringBuilder identity = new StringBuilder(
                "packstar|standalone|pac-v1");
        identity.append('|').append(rcs.getNumPos());
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            identity.append('|').append(pos).append(':');
            for (int rc : rcs.get(pos)) {
                identity.append(rc).append(',');
            }
        }
        return identity.toString();
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

        System.out.println("[PACK*] random stream: seed="
                + Long.toUnsignedString(randomSeed)
                + ", identityHash="
                + Long.toUnsignedString(
                deriveRandomSeed(0L, randomStreamIdentity), 16));

        // Phase 0: Z_min from existing DP (already computed)
        double logZMin = rootedRootEdge.getLogZUpper(0);
        double logZRigid = rootedRootEdge.getLogZLower(0);
        System.out.println("[PACK*] Phase 0: logZ_min=" + String.format("%.4f", logZMin)
                + ", logZ_rigid=" + String.format("%.4f", logZRigid)
                + ", gap=" + String.format("%.4f", logZMin - logZRigid));

        if (logZMin == Double.NEGATIVE_INFINITY) {
            writeFrequencySeverityInitialDpFailureQuietly(
                    logZMin, logZRigid, "initial DP upper bound is zero");
            setZeroBounds("initial DP upper bound is zero");
            System.out.println("[PACK*] Degenerate zero-mass state; skipping sampling and returning epsilon=1.0");
            return epsilon;
        }
        if (!Double.isFinite(logZMin)) {
            writeFrequencySeverityInitialDpFailureQuietly(
                    logZMin, logZRigid,
                    "initial DP upper bound is non-finite: " + logZMin);
            setZeroBounds("initial DP upper bound is non-finite: " + logZMin);
            System.out.println("[PACK*] Degenerate non-finite state; skipping sampling and returning epsilon=1.0");
            return epsilon;
        }

        System.out.println("[PACK*-2stage] Strict pilot holdout is enabled"
                + ", samples(train/pilot/maxEst)="
                + trainSamples + "/" + pilotSamples + "/" + maxEstSamples
                + ", fixedBulkSamples=" + fixedBulkSamples
                + ", fixedTailSamples=" + fixedTailSamples
                + ", sampleBudget="
                + (sampleBudget == Integer.MAX_VALUE ? "unbounded" : sampleBudget)
                + ", residualBound="
                + (Double.isFinite(residualBoundKcal)
                ? String.format("%.4f kcal/mol", residualBoundKcal)
                : "not configured")
                + ", conditionalRepairAudit=" + conditionalRepairAudit
                + (conditionalRepairAudit
                ? ", residualBoundGrid=" + Arrays.toString(residualBoundGridKcal)
                    + ", monitorSamples=" + monitorSamples
                    + ", dpRich=" + dpRich
                    + ", etaV5ThresholdOracle="
                    + etaV5ThresholdOracle
                    + ", etaV5MechanismPilot="
                    + etaV5MechanismPilot
                : "")
                + ", unconditionalTail=" + unconditionalTail
                + (unconditionalTail ? " (residualBound diagnostic only)" : ""));
        // Retained for diagnostics and archived estimators.  The production
        // conditional path never uses q_m to tighten or rescue an interval.
        this.logZMinDet = logZMin;
        runTwoStagePAC(startTime, logZRigid);

        if (!conditionalRepairAudit && !unconditionalTail
                && isValidCertificate(zLower, zUpper, epsilon)) {
            certificateValid = true;
            certificateFailureReason = "";
        }

        return epsilon;
    }

    private void printFinalSummary(long startTime) {
        long totalTime = System.currentTimeMillis() - startTime;
        System.out.println("[PACK*] Total: " + totalTime + " ms, " + totalCCDCalls
                + " CCD calls, epsilon=" + String.format("%.6f", epsilon)
                + ", confidence=" + String.format("%.2f%%", (1.0 - delta) * 100));
        System.out.println("[PACK*] Z bounds: lower="
                + String.format(Locale.ROOT, "%.6e", zLower)
                + ", upper="
                + String.format(Locale.ROOT, "%.6e", zUpper)
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
        if (conditionalRepairAudit) {
            runConditionalRepairAuditPAC(startTime, logZRigid);
        } else if (unconditionalTail) {
            runUnconditionalBulkTailPAC(startTime, logZRigid);
        } else {
            System.out.println("[PACK*-2stage] WARNING: legacy conditional residual-bound path enabled; "
                    + "the sampling upper certificate requires the configured global residual bound.");
            runConditionalTwoStagePAC(startTime, logZRigid);
        }
    }

    private static class TiltProposal {
        final double lambda;
        final double logZ;
        final double logRange;
        final double logCap;
        final double logRadiusProxy;
        final EnergyMatrix emat;
        final int evaluatedPositiveLambdas;
        final int dpSweeps;

        TiltProposal(double lambda, double logZ, double logRange,
                     double logCap, double logRadiusProxy,
                     EnergyMatrix emat, int evaluatedPositiveLambdas,
                     int dpSweeps) {
            this.lambda = lambda;
            this.logZ = logZ;
            this.logRange = logRange;
            this.logCap = logCap;
            this.logRadiusProxy = logRadiusProxy;
            this.emat = emat;
            this.evaluatedPositiveLambdas = evaluatedPositiveLambdas;
            this.dpSweeps = dpSweeps;
        }
    }

    private static class UnitIntervalCI {
        final double mean;
        final double variance;
        final double delta;
        final double lower;
        final double upper;

        UnitIntervalCI(double mean, double variance, double delta,
                       double lower, double upper) {
            this.mean = mean;
            this.variance = variance;
            this.delta = delta;
            this.lower = lower;
            this.upper = upper;
        }
    }

    /**
     * B-free PACK* certificate.
     *
     * <p>The clipped bulk and the omitted tail are two exact nonnegative
     * components.  Their normalized estimands both lie in [0,1].  Eta repair
     * changes only proposal efficiency; correctness does not depend on a
     * residual support bound.</p>
     */
    private void runUnconditionalBulkTailPAC(long startTime, double logZRigid) {
        System.out.println("[PACK*-unconditional] exact clipped-bulk + tilted-tail certificate enabled"
                + "; configured residualBound is diagnostic only"
                + "; tiltLambdas=" + Arrays.toString(tailTiltLambdas));

        // ---- Stage A: train eta from a fixed p_m stream ----
        long tTrain = System.currentTimeMillis();
        List<CCDResult> trainCCD = runParallelCCD(sampleConformationsFromDP(
                trainSamples, stageRandom("train")));
        System.out.println("[PACK*-unconditional] train: " + trainCCD.size()
                + " CCD from p_m in " + (System.currentTimeMillis() - tTrain) + " ms");
        if (trainCCD.isEmpty()
                || !validateObservedLowerBound(trainCCD, "train")) {
            failCertificate("unconditional: no valid train samples or observed E_m lower-bound violation",
                    startTime);
            return;
        }

        List<CCDResult> trainPool = new ArrayList<>(trainCCD);
        EtaCorrections eta = etaEnabled
                ? extractEtaCorrections(trainPool)
                : zeroEtaCorrections();
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        maybeDumpTrainingSamples(trainCCD, eta, correctedEmat);
        long tCorrected = System.currentTimeMillis();
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta);
        double logZCorrected = correctedDP.logZCorrected;
        if (!Double.isFinite(logZCorrected)) {
            failCertificate("unconditional: initial corrected DP normalizer is non-finite", startTime);
            return;
        }
        System.out.println("[PACK*-unconditional] corrected DP: logZ_eta="
                + String.format(Locale.ROOT, "%.6f", logZCorrected)
                + ", logZ_eta-logZ_rigid="
                + String.format(Locale.ROOT, "%.6f", logZCorrected - logZRigid)
                + ", oneBodyTerms=" + eta.oneBodyCount
                + ", pairTerms=" + eta.pairCount
                + ", elapsed=" + (System.currentTimeMillis() - tCorrected) + " ms");

        // ---- Stage B: proposal adaptation.  Every sample used here is excluded
        // from the final certificate. ----
        int adaptationCCDCalls = 0;
        boolean iterate = etaEnabled
                && getConfigBoolean("packstar.pac.iterate", DEFAULT_ITERATE);
        int iterateMaxRounds = iterate
                ? Math.max(0, getConfigInteger("packstar.pac.iterate.maxRounds",
                DEFAULT_ITERATE_MAX_ROUNDS)) : 0;
        int repairMaxRounds = etaEnabled && etaRepairEnabled
                ? etaRepairMaxRounds : 0;
        int adaptationMaxRounds = Math.max(iterateMaxRounds, repairMaxRounds);
        double collapseThresh = getConfigDouble(
                "packstar.pac.iterate.meanWThreshold", 0.3);
        double driftFracThresh = getConfigDouble(
                "packstar.pac.iterate.driftFraction", 0.2);
        int minTrainCount = getConfigInteger(
                "packstar.pac.iterate.minTrainCount", 5);
        double finalRepairAlpha = 1.0;
        boolean repairActive = false;
        int repairApplications = 0;
        List<CCDResult> repairWitnessPool = new ArrayList<>();

        if (adaptationMaxRounds > 0) {
            List<CCDResult> adaptationPilot = runParallelCCD(sampleConformationsFromDP(
                    Math.min(pilotSamples, maxEstSamples), stageRandom("adapt-0")));
            adaptationCCDCalls += adaptationPilot.size();
            if (adaptationPilot.isEmpty()
                    || !validateObservedLowerBound(adaptationPilot, "adapt-0")) {
                failCertificate("unconditional: invalid initial adaptation pilot", startTime);
                return;
            }

            for (int round = 1; round <= adaptationMaxRounds; round++) {
                double meanWPilot = estimateMeanW(adaptationPilot, correctedEmat);
                double driftFrac = pilotUndertrainedFraction(
                        adaptationPilot, eta, minTrainCount);
                double maxOverCorrection = maxOverCorrectionKcal(
                        adaptationPilot, correctedEmat);
                boolean refineRequested = round <= iterateMaxRounds
                        && ((Double.isFinite(meanWPilot) && meanWPilot < collapseThresh)
                        || driftFrac >= driftFracThresh);
                boolean repairRequested = repairApplications < repairMaxRounds
                        && Double.isFinite(maxOverCorrection)
                        && maxOverCorrection > etaRepairTriggerKcal;

                System.out.println("[PACK*-eta-repair] check round " + round + "/"
                        + adaptationMaxRounds
                        + ": meanW=" + String.format(Locale.ROOT, "%.6f", meanWPilot)
                        + ", driftFrac=" + String.format(Locale.ROOT, "%.4f", driftFrac)
                        + ", maxOverCorrectionKcal="
                        + String.format(Locale.ROOT, "%.6f", maxOverCorrection)
                        + ", refit=" + refineRequested
                        + ", repair=" + repairRequested
                        + ", repairApplications=" + repairApplications
                        + "/" + repairMaxRounds
                        + ", repairActive=" + repairActive);
                if (!(refineRequested || repairRequested)) break;

                trainPool.addAll(adaptationPilot);
                if (repairRequested) {
                    repairActive = true;
                    repairApplications++;
                }
                if (repairActive) {
                    // A later refit rebuilds eta from the enlarged training pool.
                    // Re-apply every held-out repair constraint so that refitting
                    // cannot silently erase a previously activated scalar repair.
                    repairWitnessPool.addAll(adaptationPilot);
                }
                EtaCorrections candidateEta = extractEtaCorrections(trainPool);
                EnergyMatrix candidateEmat = buildCorrectedEmat(candidateEta);
                double alpha = repairActive
                        ? chooseEtaRepairAlpha(repairWitnessPool, candidateEmat)
                        : 1.0;
                if (!Double.isFinite(alpha)) {
                    failCertificate("unconditional: eta repair encountered an observed E_m lower-bound violation",
                            startTime);
                    return;
                }
                if (alpha < 1.0) {
                    candidateEta = candidateEta.scaled(alpha);
                    candidateEmat = buildCorrectedEmat(candidateEta);
                }

                CorrectedDPResult candidateDP = recomputeDP(candidateEmat, candidateEta);
                if (!Double.isFinite(candidateDP.logZCorrected)) {
                    failCertificate("unconditional: eta repair produced non-finite corrected DP",
                            startTime);
                    return;
                }
                eta = candidateEta;
                correctedEmat = candidateEmat;
                logZCorrected = candidateDP.logZCorrected;
                finalRepairAlpha = alpha;

                adaptationPilot = runParallelCCD(sampleConformationsFromDP(
                        Math.min(pilotSamples, maxEstSamples),
                        stageRandom("adapt-" + round)));
                adaptationCCDCalls += adaptationPilot.size();
                if (adaptationPilot.isEmpty()
                        || !validateObservedLowerBound(adaptationPilot,
                        "adapt-" + round)) {
                    failCertificate("unconditional: invalid fresh adaptation pilot after repair",
                            startTime);
                    return;
                }
                System.out.println("[PACK*-eta-repair] applied round " + round
                        + ": alpha=" + String.format(Locale.ROOT, "%.6f", alpha)
                        + ", trainPool=" + trainPool.size()
                        + ", repairWitnesses=" + repairWitnessPool.size()
                        + ", repairApplications=" + repairApplications
                        + ", freshPilot=" + adaptationPilot.size()
                        + ", logZ_eta="
                        + String.format(Locale.ROOT, "%.6f", logZCorrected));
            }
        }

        // ---- Stage C: eta is frozen. Draw the fixed-budget final bulk set
        // while p_eta is still loaded. B/lambda selection below is a function
        // only of exact normalizers and the pre-registered sample counts; it
        // never inspects these final samples. ----
        int nBulk = fixedBulkSamples;
        List<CCDResult> bulkCCD = runParallelCCD(sampleConformationsFromDP(
                nBulk, stageRandom("bulk-final")));
        if (bulkCCD.isEmpty()
                || !validateObservedLowerBound(bulkCCD, "bulk-final")) {
            failCertificate("unconditional: invalid final bulk sample", startTime);
            return;
        }

        // ---- Stage D: jointly choose the normalizer-balanced B/lambda pair
        // and load p_lambda.  The final tail count is fixed before sampling. ----
        int nTail = fixedTailSamples;
        TiltProposal tilt = selectTiltProposal(
                eta, logZCorrected, nBulk, nTail);
        if (tilt == null || !Double.isFinite(tilt.logZ)
                || !Double.isFinite(tilt.logRange)
                || !Double.isFinite(tilt.logCap)) {
            failCertificate("unconditional: no finite dynamic-B/tilted-tail proposal", startTime);
            return;
        }
        clipLogCap = tilt.logCap;
        System.out.println("[PACK*-unconditional] frozen certificate design: repairAlpha="
                + String.format(Locale.ROOT, "%.6f", finalRepairAlpha)
                + ", B_kcal=" + String.format(Locale.ROOT, "%.6f", RT * clipLogCap)
                + ", clipLogCap=" + String.format(Locale.ROOT, "%.6f", clipLogCap)
                + ", lambda=" + String.format(Locale.ROOT, "%.6f", tilt.lambda)
                + ", logRadiusProxy="
                + String.format(Locale.ROOT, "%.8f", tilt.logRadiusProxy)
                + ", nBulk=" + nBulk
                + ", nTail=" + nTail);
        List<CCDResult> tailCCD = runParallelCCD(sampleConformationsFromDP(
                nTail, stageRandom("tail-final")));
        if (tailCCD.isEmpty()
                || !validateObservedLowerBound(tailCCD, "tail-final")) {
            failCertificate("unconditional: invalid final tilted-tail sample", startTime);
            return;
        }

        totalCCDCalls = trainCCD.size() + adaptationCCDCalls
                + bulkCCD.size() + tailCCD.size();
        computeUnconditionalBulkTailBound(
                bulkCCD, tailCCD, correctedEmat, logZCorrected,
                clipLogCap, tilt);
        printFinalSummary(startTime);
    }

    private static class FrequencySeverityEtaCoverage {
        final int[][] unaryDistinctContexts;
        final Map<Long, Integer> pairDistinctContexts;

        FrequencySeverityEtaCoverage(int numPos, int[] numRCs) {
            unaryDistinctContexts = new int[numPos][];
            for (int pos = 0; pos < numPos; pos++) {
                unaryDistinctContexts[pos] = new int[numRCs[pos]];
            }
            pairDistinctContexts = new HashMap<>();
        }

        int unary(int pos, int rc) {
            return unaryDistinctContexts[pos][rc];
        }

        int pair(int pos1, int rc1, int pos2, int rc2) {
            return pairDistinctContexts.getOrDefault(
                    EtaCorrections.packPairKey(pos1, rc1, pos2, rc2), 0);
        }
    }

    private static class FrequencySeverityEtaCellObservation {
        final int sampleIndex;
        final int fold;
        final boolean pair;
        final int pos1;
        final int rc1;
        final int pos2;
        final int rc2;
        final double correctionKcal;

        FrequencySeverityEtaCellObservation(
                int sampleIndex, int fold, boolean pair,
                int pos1, int rc1, int pos2, int rc2,
                double correctionKcal) {
            this.sampleIndex = sampleIndex;
            this.fold = fold;
            this.pair = pair;
            this.pos1 = pos1;
            this.rc1 = rc1;
            this.pos2 = pos2;
            this.rc2 = rc2;
            this.correctionKcal = correctionKcal;
        }
    }

    private static class FrequencySeverityEtaTraining {
        final EtaCorrections all;
        final EtaCorrections[] folds;
        final FrequencySeverityEtaCoverage allCoverage;
        final FrequencySeverityEtaCoverage[] foldCoverage;
        final List<FrequencySeverityEtaCellObservation> observations;
        final int allFitSampleCount;
        final int[] foldFitSampleCounts;
        final int historySampleCount;
        final List<CCDResult> allFitSamples;
        final List<CCDResult>[] foldFitSamples;

        FrequencySeverityEtaTraining(EtaCorrections all, EtaCorrections[] folds,
                      FrequencySeverityEtaCoverage allCoverage,
                      FrequencySeverityEtaCoverage[] foldCoverage,
                      List<FrequencySeverityEtaCellObservation> observations,
                      int allFitSampleCount, int[] foldFitSampleCounts,
                      int historySampleCount,
                      List<CCDResult> allFitSamples,
                      List<CCDResult>[] foldFitSamples) {
            this.all = all;
            this.folds = folds;
            this.allCoverage = allCoverage;
            this.foldCoverage = foldCoverage;
            this.observations = observations;
            this.allFitSampleCount = allFitSampleCount;
            this.foldFitSampleCounts = foldFitSampleCounts;
            this.historySampleCount = historySampleCount;
            this.allFitSamples = allFitSamples;
            this.foldFitSamples = foldFitSamples;
        }
    }

    private static class FrequencySeverityFoldScore {
        final int validationFold;
        final int sampleCount;
        final double logMu;
        final double shiftEss;
        final double shiftEssFraction;
        final double targetM2;
        final double targetEssFraction;
        final double bulkMean;
        final double bulkVariance;
        final double tailProbability;
        final double baselineUnaryUndertrainedMass;
        final double candidateUnaryUndertrainedMass;
        final double unaryUndertrainedAmplification;
        final double baselinePairUndertrainedMass;
        final double candidatePairUndertrainedMass;
        final double pairUndertrainedAmplification;
        final PackStarFrequencySeverityPAC.Sizing sizing;
        final boolean finite;
        final boolean shiftPass;
        final boolean coveragePass;
        final boolean bulkPass;
        final boolean reachabilityPass;

        FrequencySeverityFoldScore(int validationFold, int sampleCount,
                    double logMu, double shiftEss,
                    double shiftEssFraction, double targetM2,
                    double targetEssFraction, double bulkMean,
                    double bulkVariance, double tailProbability,
                    double baselineUnaryUndertrainedMass,
                    double candidateUnaryUndertrainedMass,
                    double unaryUndertrainedAmplification,
                    double baselinePairUndertrainedMass,
                    double candidatePairUndertrainedMass,
                    double pairUndertrainedAmplification,
                    PackStarFrequencySeverityPAC.Sizing sizing,
                    boolean finite, boolean shiftPass,
                    boolean coveragePass, boolean bulkPass,
                    boolean reachabilityPass) {
            this.validationFold = validationFold;
            this.sampleCount = sampleCount;
            this.logMu = logMu;
            this.shiftEss = shiftEss;
            this.shiftEssFraction = shiftEssFraction;
            this.targetM2 = targetM2;
            this.targetEssFraction = targetEssFraction;
            this.bulkMean = bulkMean;
            this.bulkVariance = bulkVariance;
            this.tailProbability = tailProbability;
            this.baselineUnaryUndertrainedMass =
                    baselineUnaryUndertrainedMass;
            this.candidateUnaryUndertrainedMass =
                    candidateUnaryUndertrainedMass;
            this.unaryUndertrainedAmplification =
                    unaryUndertrainedAmplification;
            this.baselinePairUndertrainedMass =
                    baselinePairUndertrainedMass;
            this.candidatePairUndertrainedMass =
                    candidatePairUndertrainedMass;
            this.pairUndertrainedAmplification =
                    pairUndertrainedAmplification;
            this.sizing = sizing;
            this.finite = finite;
            this.shiftPass = shiftPass;
            this.coveragePass = coveragePass;
            this.bulkPass = bulkPass;
            this.reachabilityPass = reachabilityPass;
        }

        boolean eligible() {
            return finite && shiftPass && coveragePass
                    && bulkPass && reachabilityPass;
        }
    }

    /**
     * Sample-size-weighted out-of-fold model-selection score.  This pools the
     * two independently validated eta fits for a less noisy adaptation
     * diagnostic; it is not a confidence interval and is always followed by a
     * full-refit audit and an independent on-policy pilot.
     */
    private static class FrequencySeverityCrossfitScore {
        final int sampleCount;
        final double logMu;
        final double shiftEss;
        final double shiftEssFraction;
        final double targetM2;
        final double targetEssFraction;
        final double bulkMean;
        final double bulkVariance;
        final double tailProbability;
        final double baselineUnaryUndertrainedMass;
        final double candidateUnaryUndertrainedMass;
        final double unaryUndertrainedAmplification;
        final double baselinePairUndertrainedMass;
        final double candidatePairUndertrainedMass;
        final double pairUndertrainedAmplification;
        final PackStarFrequencySeverityPAC.Sizing sizing;
        final boolean allFoldsFinite;
        final boolean allFoldsShiftPass;
        final boolean allFoldsCoveragePass;
        final boolean finite;
        final boolean shiftPass;
        final boolean coveragePass;
        final boolean bulkPass;
        final boolean reachabilityPass;

        FrequencySeverityCrossfitScore(
                int sampleCount, double logMu,
                double shiftEss, double shiftEssFraction,
                double targetM2, double targetEssFraction,
                double bulkMean, double bulkVariance,
                double tailProbability,
                double baselineUnaryUndertrainedMass,
                double candidateUnaryUndertrainedMass,
                double unaryUndertrainedAmplification,
                double baselinePairUndertrainedMass,
                double candidatePairUndertrainedMass,
                double pairUndertrainedAmplification,
                PackStarFrequencySeverityPAC.Sizing sizing,
                boolean allFoldsFinite,
                boolean allFoldsShiftPass,
                boolean allFoldsCoveragePass,
                boolean finite, boolean shiftPass,
                boolean coveragePass, boolean bulkPass,
                boolean reachabilityPass) {
            this.sampleCount = sampleCount;
            this.logMu = logMu;
            this.shiftEss = shiftEss;
            this.shiftEssFraction = shiftEssFraction;
            this.targetM2 = targetM2;
            this.targetEssFraction = targetEssFraction;
            this.bulkMean = bulkMean;
            this.bulkVariance = bulkVariance;
            this.tailProbability = tailProbability;
            this.baselineUnaryUndertrainedMass =
                    baselineUnaryUndertrainedMass;
            this.candidateUnaryUndertrainedMass =
                    candidateUnaryUndertrainedMass;
            this.unaryUndertrainedAmplification =
                    unaryUndertrainedAmplification;
            this.baselinePairUndertrainedMass =
                    baselinePairUndertrainedMass;
            this.candidatePairUndertrainedMass =
                    candidatePairUndertrainedMass;
            this.pairUndertrainedAmplification =
                    pairUndertrainedAmplification;
            this.sizing = sizing;
            this.allFoldsFinite = allFoldsFinite;
            this.allFoldsShiftPass = allFoldsShiftPass;
            this.allFoldsCoveragePass = allFoldsCoveragePass;
            this.finite = finite;
            this.shiftPass = shiftPass;
            this.coveragePass = coveragePass;
            this.bulkPass = bulkPass;
            this.reachabilityPass = reachabilityPass;
        }

        boolean eligible() {
            // Half-fold coverage is deliberately diagnostic: fitting on 250
            // observations must not be required to have the same sparse-cell
            // support as the actual all-500 refit.  The all-training coverage
            // audit remains a hard gate below.
            return allFoldsFinite && allFoldsShiftPass
                    && finite && shiftPass && bulkPass
                    && reachabilityPass;
        }
    }

    private static class FrequencySeverityCandidateScore {
        final String id;
        final FrequencySeverityShrinkPair shrink;
        final double alpha;
        final double tripleEtaScale;
        final EtaCorrections fullEta;
        final FrequencySeverityFoldScore[] folds;
        final FrequencySeverityCrossfitScore crossfit;
        final boolean crossfitEligible;
        final FrequencySeverityRefitAudit refitAudit;
        final boolean eligible;

        FrequencySeverityCandidateScore(String id, FrequencySeverityShrinkPair shrink,
                         double alpha, double tripleEtaScale,
                         EtaCorrections fullEta,
                         FrequencySeverityFoldScore[] folds,
                         FrequencySeverityCrossfitScore crossfit,
                         FrequencySeverityRefitAudit refitAudit) {
            this.id = id;
            this.shrink = shrink;
            this.alpha = alpha;
            this.tripleEtaScale = tripleEtaScale;
            this.fullEta = fullEta;
            this.folds = folds;
            this.crossfit = crossfit;
            this.refitAudit = refitAudit;
            this.crossfitEligible = crossfit != null
                    && crossfit.eligible();
            this.eligible = crossfitEligible && refitAudit != null
                    && refitAudit.eligible();
        }
    }

    private static class FrequencySeverityRefitAudit {
        final double shiftEssFraction;
        final double baselineUnaryUndertrainedMass;
        final double candidateUnaryUndertrainedMass;
        final double unaryUndertrainedAmplification;
        final double baselinePairUndertrainedMass;
        final double candidatePairUndertrainedMass;
        final double pairUndertrainedAmplification;
        final boolean finite;
        final boolean shiftPass;
        final boolean coveragePass;

        FrequencySeverityRefitAudit(
                double shiftEssFraction,
                double baselineUnaryUndertrainedMass,
                double candidateUnaryUndertrainedMass,
                double unaryUndertrainedAmplification,
                double baselinePairUndertrainedMass,
                double candidatePairUndertrainedMass,
                double pairUndertrainedAmplification,
                boolean finite, boolean shiftPass,
                boolean coveragePass) {
            this.shiftEssFraction = shiftEssFraction;
            this.baselineUnaryUndertrainedMass =
                    baselineUnaryUndertrainedMass;
            this.candidateUnaryUndertrainedMass =
                    candidateUnaryUndertrainedMass;
            this.unaryUndertrainedAmplification =
                    unaryUndertrainedAmplification;
            this.baselinePairUndertrainedMass =
                    baselinePairUndertrainedMass;
            this.candidatePairUndertrainedMass =
                    candidatePairUndertrainedMass;
            this.pairUndertrainedAmplification =
                    pairUndertrainedAmplification;
            this.finite = finite;
            this.shiftPass = shiftPass;
            this.coveragePass = coveragePass;
        }

        boolean eligible() {
            return finite && shiftPass && coveragePass;
        }
    }

    private static class FrequencySeverityPilotEvaluation {
        final double[] logRelativeWeights;
        final PackStarFrequencySeverityPAC.Interval interval;
        final PackStarFrequencySeverityPAC.Sizing sizing;

        FrequencySeverityPilotEvaluation(
                double[] logRelativeWeights,
                PackStarFrequencySeverityPAC.Interval interval,
                PackStarFrequencySeverityPAC.Sizing sizing) {
            this.logRelativeWeights = logRelativeWeights;
            this.interval = interval;
            this.sizing = sizing;
        }
    }

    private List<FrequencySeverityCandidateScore> scoreFrequencySeverityEtaCandidates(
            List<CCDResult> trainingSamples,
            FrequencySeverityEtaTraining training) {
        List<CCDResult>[] validationFolds = new List[frequencySeverityFolds];
        for (int fold = 0; fold < frequencySeverityFolds; fold++) {
            validationFolds[fold] = new ArrayList<>();
        }
        for (int index = 0; index < trainingSamples.size(); index++) {
            validationFolds[index % frequencySeverityFolds].add(trainingSamples.get(index));
        }

        List<FrequencySeverityCandidateScore> candidates = new ArrayList<>();
        // Alpha remains fixed at one.  Pair-only and one structural
        // selected-triple model are both cross-fitted; this is model-order
        // selection, not attenuation of a fitted correction.
        final double alpha = 1.0;
        for (FrequencySeverityShrinkPair shrink : frequencySeverityShrinkGrid) {
            EtaCorrections fullPairEta = shrinkEta(
                    training.all, training.allCoverage,
                    shrink.unary, shrink.pair).scaled(alpha);
            EtaCorrections[] foldPairEta =
                    new EtaCorrections[frequencySeverityFolds];
            for (int fitFold = 0;
                 fitFold < frequencySeverityFolds; fitFold++) {
                foldPairEta[fitFold] = shrinkEta(
                        training.folds[fitFold],
                        training.foldCoverage[fitFold],
                        shrink.unary, shrink.pair).scaled(alpha);
            }

            FrequencySeverityFoldScore[] pairFoldScores =
                    new FrequencySeverityFoldScore[frequencySeverityFolds];
            for (int validationFold = 0;
                 validationFold < frequencySeverityFolds; validationFold++) {
                int fitFold = 1 - validationFold;
                pairFoldScores[validationFold] =
                        scoreFrequencySeverityEtaFold(
                        validationFold,
                        validationFolds[validationFold],
                        foldPairEta[fitFold],
                        training.foldCoverage[fitFold]);
            }
            FrequencySeverityCrossfitScore pairCrossfit =
                    poolFrequencySeverityFoldScores(pairFoldScores);
            FrequencySeverityRefitAudit pairRefitAudit =
                    auditFrequencySeverityRefit(
                            trainingSamples, fullPairEta,
                            training.allCoverage);
            candidates.add(new FrequencySeverityCandidateScore(
                    frequencySeverityCandidateId(
                            shrink, alpha, 0.0),
                    shrink, alpha, 0.0, fullPairEta,
                    pairFoldScores, pairCrossfit, pairRefitAudit));

            if (frequencySeverityTripleEta == null) continue;
            PackStarTripleEtaCorrections fullSelectedTriple =
                    fitFrequencySeverityTripleEta(
                            training.allFitSamples, fullPairEta);
            if (!hasTripleEta(fullSelectedTriple)) continue;
            PackStarTripleEtaCorrections[] foldSelectedTriples =
                    new PackStarTripleEtaCorrections[
                    frequencySeverityFolds];
            EtaCorrections[] foldTripleEta =
                    new EtaCorrections[frequencySeverityFolds];
            for (int fitFold = 0;
                 fitFold < frequencySeverityFolds; fitFold++) {
                foldSelectedTriples[fitFold] =
                        fitFrequencySeverityTripleEta(
                                training.foldFitSamples[fitFold],
                                foldPairEta[fitFold]);
                foldTripleEta[fitFold] =
                        attachFrequencySeverityTripleEta(
                                foldPairEta[fitFold],
                                foldSelectedTriples[fitFold],
                                frequencySeverityTripleEtaScale);
            }
            EtaCorrections fullTripleEta =
                    attachFrequencySeverityTripleEta(
                            fullPairEta, fullSelectedTriple,
                            frequencySeverityTripleEtaScale);
            FrequencySeverityFoldScore[] tripleFoldScores =
                    new FrequencySeverityFoldScore[
                    frequencySeverityFolds];
            for (int validationFold = 0;
                 validationFold < frequencySeverityFolds;
                 validationFold++) {
                int fitFold = 1 - validationFold;
                tripleFoldScores[validationFold] =
                        scoreFrequencySeverityEtaFold(
                                validationFold,
                                validationFolds[validationFold],
                                foldTripleEta[fitFold],
                                training.foldCoverage[fitFold]);
            }
            FrequencySeverityCrossfitScore tripleCrossfit =
                    poolFrequencySeverityFoldScores(tripleFoldScores);
            FrequencySeverityRefitAudit tripleRefitAudit =
                    auditFrequencySeverityRefit(
                            trainingSamples, fullTripleEta,
                            training.allCoverage);
            candidates.add(new FrequencySeverityCandidateScore(
                    frequencySeverityCandidateId(
                            shrink, alpha, 1.0),
                    shrink, alpha, 1.0, fullTripleEta,
                    tripleFoldScores, tripleCrossfit,
                    tripleRefitAudit));
        }
        return candidates;
    }

    static String frequencySeverityCandidateId(
            FrequencySeverityShrinkPair shrink, double alpha,
            double tripleScale) {
        if (alpha != 1.0) {
            throw new IllegalArgumentException(
                    "frequency/severity alpha is fixed to 1: " + alpha);
        }
        if (tripleScale != 0.0 && tripleScale != 1.0) {
            throw new IllegalArgumentException(
                    "frequency/severity gamma is fixed to 1 when triple"
                            + " residuals are enabled: " + tripleScale);
        }
        String pair = shrink.id() + String.format(
                Locale.ROOT, "-alpha-%g", alpha);
        return tripleScale == 0.0
                ? pair + "-pair-only"
                : pair + String.format(Locale.ROOT,
                "-plus-triple-eta-gamma-%g", tripleScale);
    }

    private EtaCorrections attachFrequencySeverityTripleEta(
            EtaCorrections pairEta,
            PackStarTripleEtaCorrections selectedTriple,
            double tripleScale) {
        if (!hasTripleEta(selectedTriple)
                || tripleScale == 0.0) {
            return pairEta;
        }
        return pairEta.withTripleEta(
                selectedTriple,
                tripleScale);
    }

    private PackStarTripleEtaCorrections fitFrequencySeverityTripleEta(
            List<CCDResult> samples, EtaCorrections pairEta) {
        if (frequencySeverityTripleEta == null
                || samples == null || samples.size() < 2) {
            return null;
        }
        int[][] conformations = new int[samples.size()][];
        double[] residuals = new double[samples.size()];
        for (int sample = 0; sample < samples.size(); sample++) {
            CCDResult result = samples.get(sample);
            conformations[sample] = result.conf;
            residuals[sample] = result.eTrue - result.eMin
                    - computeEtaEnergy(result.conf, pairEta);
        }
        return frequencySeverityTripleEta.fitSelectedResidual(
                rcs, interactionGraph, conformations, residuals,
                frequencySeverityTripleEtaMaxPositionTriples,
                frequencySeverityTripleEtaMaxFillEdges,
                frequencySeverityTripleEtaMinCellContexts,
                frequencySeverityTripleEtaPriorStrength,
                frequencySeverityTripleEtaLocalCapKcal,
                frequencySeverityTripleEtaResidualCapKcal,
                frequencySeverityTripleEtaMaxAssignments);
    }

    private static boolean hasTripleEta(EtaCorrections eta) {
        return eta != null && eta.tripleEta != null
                && eta.tripleEtaScale != 0.0
                && eta.tripleEta.factorAssignments > 0L;
    }

    private static boolean hasTripleEta(
            PackStarTripleEtaCorrections tripleEta) {
        return tripleEta != null
                && tripleEta.factorAssignments > 0L;
    }

    private static PackStarTripleEtaCorrections.ResidualSummary
    summarizeFrequencySeverityTripleResidual(EtaCorrections eta) {
        if (!hasTripleEta(eta)) return null;
        return eta.tripleEta.summarizeResidual(eta::getPairEta);
    }

    private FrequencySeverityFoldScore scoreFrequencySeverityEtaFold(
            int validationFold,
            List<CCDResult> samples,
            EtaCorrections eta,
            FrequencySeverityEtaCoverage fitCoverage) {
        int n = samples.size();
        if (n <= 1) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }
        double[] logA = new double[n];
        double[] logV = new double[n];
        boolean[] unaryUndertrained = new boolean[n];
        boolean[] pairUndertrained = new boolean[n];
        for (int i = 0; i < n; i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + computeEtaEnergy(sample.conf, eta);
            logA[i] = frequencySeveritySourceLogWeight(
                    sample, targetEnergy);
            logV[i] = logA[i]
                    + (targetEnergy - sample.eTrue) / RT;
            unaryUndertrained[i] = touchesFrequencySeverityUndertrainedUnary(
                    sample.conf, fitCoverage);
            pairUndertrained[i] = touchesFrequencySeverityUndertrainedPair(
                    sample.conf, fitCoverage);
            if (!Double.isFinite(logA[i]) || !Double.isFinite(logV[i])) {
                return invalidFrequencySeverityFoldScore(validationFold, n);
            }
        }

        double logSumA = logSumExp(logA);
        double logSumV = logSumExp(logV);
        if (!Double.isFinite(logSumA) || !Double.isFinite(logSumV)) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }
        double logMu = logSumV - logSumA;
        double[] proposalWeights = new double[n];
        double[] logRelative = new double[n];
        double sumWeightSquares = 0.0;
        double bulkMean = 0.0;
        double tailProbability = 0.0;
        double candidateUnaryUndertrained = 0.0;
        double candidatePairUndertrained = 0.0;
        int baselineUnaryUndertrained = 0;
        int baselinePairUndertrained = 0;
        double logClip = frequencySeverityRelativeBoundKcal / RT;
        double[] bulk = new double[n];
        for (int i = 0; i < n; i++) {
            double weight = Math.exp(logA[i] - logSumA);
            double logR = logV[i] - logA[i] - logMu;
            if (!Double.isFinite(weight) || weight < 0.0
                    || !Double.isFinite(logR)) {
                return invalidFrequencySeverityFoldScore(validationFold, n);
            }
            proposalWeights[i] = weight;
            logRelative[i] = logR;
            sumWeightSquares += weight * weight;
            double value = logR >= logClip
                    ? 1.0 : Math.exp(logR - logClip);
            bulk[i] = value;
            bulkMean += weight * value;
            if (logR > logClip) tailProbability += weight;
            if (unaryUndertrained[i]) {
                baselineUnaryUndertrained++;
                candidateUnaryUndertrained += weight;
            }
            if (pairUndertrained[i]) {
                baselinePairUndertrained++;
                candidatePairUndertrained += weight;
            }
        }
        if (!(sumWeightSquares > 0.0)
                || !Double.isFinite(sumWeightSquares)) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }
        double shiftEss = 1.0 / sumWeightSquares;
        double shiftEssFraction = shiftEss / n;
        double varianceDenominator = 1.0 - sumWeightSquares;
        if (!(varianceDenominator > 0.0)) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }
        double bulkVariance = 0.0;
        for (int i = 0; i < n; i++) {
            double centered = bulk[i] - bulkMean;
            bulkVariance += proposalWeights[i] * centered * centered;
        }
        bulkVariance /= varianceDenominator;

        double logMeanR = weightedLogMeanExp(
                proposalWeights, logRelative, 1.0);
        double logSecondR = weightedLogMeanExp(
                proposalWeights, logRelative, 2.0);
        double logTargetM2 = logSecondR - 2.0 * logMeanR;
        double targetM2 = logTargetM2 < Math.log(Double.MAX_VALUE)
                ? Math.max(1.0, Math.exp(logTargetM2))
                : Double.POSITIVE_INFINITY;
        double targetEssFraction = Double.isFinite(targetM2)
                ? 1.0 / targetM2 : 0.0;
        double baselineUnaryMass =
                (double) baselineUnaryUndertrained / n;
        double baselinePairMass =
                (double) baselinePairUndertrained / n;
        double unaryAmplification = massAmplification(
                candidateUnaryUndertrained, baselineUnaryMass);
        double pairAmplification = massAmplification(
                candidatePairUndertrained, baselinePairMass);

        PackStarFrequencySeverityPAC.Moments moments;
        PackStarFrequencySeverityPAC.Sizing sizing;
        PackStarFrequencySeverityPAC.Interval atMax;
        try {
            moments = new PackStarFrequencySeverityPAC.Moments(
                    clamp(bulkMean, 0.0, 1.0),
                    Math.max(0.0, bulkVariance) * nstarInflate,
                    clamp(tailProbability, 0.0, 1.0));
            sizing = PackStarFrequencySeverityPAC.size(
                    moments, maxEstSamples,
                    configuredFrequencySeverityUnreachableSamples(),
                    targetEpsilon, frequencySeveritySizeSafety, frequencySeverityCap,
                    frequencySeverityEventDelta(), frequencySeverityEventDelta());
            atMax = PackStarFrequencySeverityPAC.project(
                    maxEstSamples, moments, frequencySeverityCap,
                    frequencySeverityEventDelta(), frequencySeverityEventDelta());
        } catch (IllegalArgumentException ex) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }

        boolean finite = Double.isFinite(logMu)
                && Double.isFinite(shiftEssFraction)
                && Double.isFinite(logTargetM2)
                && Double.isFinite(bulkMean)
                && Double.isFinite(bulkVariance)
                && Double.isFinite(tailProbability)
                && Double.isFinite(unaryAmplification)
                && Double.isFinite(pairAmplification);
        boolean shiftPass = shiftEssFraction + 1.0e-12
                >= frequencySeverityMinShiftEssFraction;
        boolean coveragePass = unaryAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12
                && pairAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12;
        boolean bulkPass = atMax.bulkLower > 0.0;
        boolean reachabilityPass = sizing.reachableAtMax
                && sizing.epsilonAtMaxSamples
                <= targetEpsilon + 1.0e-12;
        return new FrequencySeverityFoldScore(
                validationFold, n, logMu, shiftEss,
                shiftEssFraction, targetM2, targetEssFraction,
                bulkMean, bulkVariance, tailProbability,
                baselineUnaryMass, candidateUnaryUndertrained,
                unaryAmplification, baselinePairMass,
                candidatePairUndertrained, pairAmplification,
                sizing, finite, shiftPass, coveragePass,
                bulkPass, reachabilityPass);
    }

    private FrequencySeverityFoldScore invalidFrequencySeverityFoldScore(int fold, int sampleCount) {
        PackStarFrequencySeverityPAC.Moments fallback =
                new PackStarFrequencySeverityPAC.Moments(0.0, 0.0, 1.0);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        fallback, maxEstSamples,
                        configuredFrequencySeverityUnreachableSamples(),
                        targetEpsilon, frequencySeveritySizeSafety, frequencySeverityCap,
                        frequencySeverityEventDelta(), frequencySeverityEventDelta());
        return new FrequencySeverityFoldScore(
                fold, sampleCount, Double.NaN, 0.0, 0.0,
                Double.POSITIVE_INFINITY, 0.0, 0.0, 0.0, 1.0,
                1.0, 1.0, Double.POSITIVE_INFINITY,
                1.0, 1.0, Double.POSITIVE_INFINITY,
                sizing, false, false, false, false, false);
    }

    private FrequencySeverityCrossfitScore poolFrequencySeverityFoldScores(
            FrequencySeverityFoldScore[] folds) {
        if (folds == null || folds.length == 0) {
            return invalidFrequencySeverityCrossfitScore();
        }
        int[] sampleCounts = new int[folds.length];
        double[] shiftEss = new double[folds.length];
        double[] bulkMeans = new double[folds.length];
        double[] bulkVariances = new double[folds.length];
        double[] tailProbabilities = new double[folds.length];
        int total = 0;
        boolean allFoldsFinite = true;
        boolean allFoldsShiftPass = true;
        boolean allFoldsCoveragePass = true;
        for (int fold = 0; fold < folds.length; fold++) {
            FrequencySeverityFoldScore score = folds[fold];
            if (score == null) return invalidFrequencySeverityCrossfitScore();
            sampleCounts[fold] = score.sampleCount;
            shiftEss[fold] = score.shiftEss;
            bulkMeans[fold] = score.bulkMean;
            bulkVariances[fold] = score.bulkVariance;
            tailProbabilities[fold] = score.tailProbability;
            total += score.sampleCount;
            allFoldsFinite &= score.finite;
            allFoldsShiftPass &= score.shiftPass;
            allFoldsCoveragePass &= score.coveragePass;
        }
        if (!allFoldsFinite || total <= 1) {
            return invalidFrequencySeverityCrossfitScore();
        }

        PackStarFrequencySeverityPAC.PooledCrossfitMoments pooled;
        try {
            pooled = PackStarFrequencySeverityPAC.poolCrossfitMoments(
                    sampleCounts, shiftEss, bulkMeans,
                    bulkVariances, tailProbabilities);
        } catch (IllegalArgumentException ex) {
            return invalidFrequencySeverityCrossfitScore();
        }

        double logMu = 0.0;
        double targetM2 = 0.0;
        double baselineUnaryMass = 0.0;
        double candidateUnaryMass = 0.0;
        double baselinePairMass = 0.0;
        double candidatePairMass = 0.0;
        for (FrequencySeverityFoldScore fold : folds) {
            double foldMass = (double) fold.sampleCount / total;
            logMu += foldMass * fold.logMu;
            targetM2 += foldMass * fold.targetM2;
            baselineUnaryMass += foldMass
                    * fold.baselineUnaryUndertrainedMass;
            candidateUnaryMass += foldMass
                    * fold.candidateUnaryUndertrainedMass;
            baselinePairMass += foldMass
                    * fold.baselinePairUndertrainedMass;
            candidatePairMass += foldMass
                    * fold.candidatePairUndertrainedMass;
        }
        double unaryAmplification = massAmplification(
                candidateUnaryMass, baselineUnaryMass);
        double pairAmplification = massAmplification(
                candidatePairMass, baselinePairMass);
        double targetEssFraction = Double.isFinite(targetM2)
                && targetM2 > 0.0 ? 1.0 / targetM2 : 0.0;

        PackStarFrequencySeverityPAC.Sizing sizing;
        PackStarFrequencySeverityPAC.Interval atMax;
        try {
            PackStarFrequencySeverityPAC.Moments inflated =
                    new PackStarFrequencySeverityPAC.Moments(
                            pooled.moments.bulkMean,
                            pooled.moments.bulkVariance * nstarInflate,
                            pooled.moments.tailProbability);
            sizing = PackStarFrequencySeverityPAC.size(
                    inflated, maxEstSamples,
                    configuredFrequencySeverityUnreachableSamples(),
                    targetEpsilon, frequencySeveritySizeSafety,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            atMax = PackStarFrequencySeverityPAC.project(
                    maxEstSamples, inflated, frequencySeverityCap,
                    frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
        } catch (IllegalArgumentException ex) {
            return invalidFrequencySeverityCrossfitScore();
        }

        boolean finite = Double.isFinite(logMu)
                && Double.isFinite(pooled.effectiveSampleSize)
                && Double.isFinite(pooled.effectiveSampleFraction)
                && Double.isFinite(targetEssFraction)
                && Double.isFinite(pooled.moments.bulkMean)
                && Double.isFinite(pooled.moments.bulkVariance)
                && Double.isFinite(pooled.moments.tailProbability)
                && Double.isFinite(unaryAmplification)
                && Double.isFinite(pairAmplification);
        boolean shiftPass = finite
                && pooled.effectiveSampleFraction + 1.0e-12
                >= frequencySeverityMinShiftEssFraction;
        boolean coveragePass = finite
                && unaryAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12
                && pairAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12;
        boolean bulkPass = atMax.bulkLower > 0.0;
        boolean reachabilityPass = sizing.reachableAtMax
                && sizing.epsilonAtMaxSamples
                <= targetEpsilon + 1.0e-12;
        return new FrequencySeverityCrossfitScore(
                pooled.sampleCount, logMu,
                pooled.effectiveSampleSize,
                pooled.effectiveSampleFraction,
                targetM2, targetEssFraction,
                pooled.moments.bulkMean,
                pooled.moments.bulkVariance,
                pooled.moments.tailProbability,
                baselineUnaryMass, candidateUnaryMass,
                unaryAmplification,
                baselinePairMass, candidatePairMass,
                pairAmplification, sizing,
                allFoldsFinite, allFoldsShiftPass,
                allFoldsCoveragePass, finite, shiftPass,
                coveragePass, bulkPass, reachabilityPass);
    }

    private FrequencySeverityCrossfitScore invalidFrequencySeverityCrossfitScore() {
        PackStarFrequencySeverityPAC.Moments fallback =
                new PackStarFrequencySeverityPAC.Moments(0.0, 0.0, 1.0);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        fallback, maxEstSamples,
                        configuredFrequencySeverityUnreachableSamples(),
                        targetEpsilon, frequencySeveritySizeSafety,
                        frequencySeverityCap, frequencySeverityEventDelta(),
                        frequencySeverityEventDelta());
        return new FrequencySeverityCrossfitScore(
                0, Double.NaN, 0.0, 0.0,
                Double.POSITIVE_INFINITY, 0.0,
                0.0, 0.0, 1.0,
                1.0, 1.0, Double.POSITIVE_INFINITY,
                1.0, 1.0, Double.POSITIVE_INFINITY,
                sizing, false, false, false,
                false, false, false, false, false);
    }

    private FrequencySeverityCandidateScore selectBestFrequencySeverityCandidate(
            List<FrequencySeverityCandidateScore> candidates) {
        FrequencySeverityCandidateScore best = null;
        for (FrequencySeverityCandidateScore candidate : candidates) {
            if (!candidate.eligible) continue;
            if (best == null || isBetterFrequencySeverityCandidate(candidate, best)) {
                best = candidate;
            }
        }
        return best;
    }

    /**
     * Adaptation-only trust-region probe.  If no proposal is ready to issue a
     * bound, use the strongest unshrunk pair-eta step that still retains at
     * least the pre-registered discovery ESS in every fold and in the full
     * refit.  At that step, choose the pre-registered triple scale with the
     * best off-policy A3 projection; this keeps pair-only and higher-order
     * proposals in the same evidence-gated family.
     * Reachability and sparse-cell coverage are deliberately not required:
     * the probe exists to obtain the on-policy data those q_m diagnostics
     * cannot supply, and it can never reach monitor/final without a fresh
     * on-policy validation batch.
     */
    private FrequencySeverityCandidateScore selectFrequencySeverityDiscoveryProbe(
            List<FrequencySeverityCandidateScore> candidates) {
        FrequencySeverityCandidateScore strongest = null;
        for (FrequencySeverityCandidateScore candidate : candidates) {
            if (!(candidate.alpha > 0.0)
                    || candidate.shrink.unary != 0.0
                    || candidate.shrink.pair != 0.0
                    || candidate.crossfit == null
                    || candidate.refitAudit == null
                    || !candidate.crossfit.finite
                    || !candidate.crossfit.allFoldsFinite
                    || !candidate.refitAudit.finite
                    || candidate.crossfit.shiftEssFraction + 1.0e-12
                    < frequencySeverityDiscoveryMinShiftEssFraction
                    || candidate.refitAudit.shiftEssFraction + 1.0e-12
                    < frequencySeverityDiscoveryMinShiftEssFraction) {
                continue;
            }
            boolean everyFoldRetainsDiscoveryOverlap = true;
            for (FrequencySeverityFoldScore fold : candidate.folds) {
                if (!fold.finite
                        || fold.shiftEssFraction + 1.0e-12
                        < frequencySeverityDiscoveryMinShiftEssFraction) {
                    everyFoldRetainsDiscoveryOverlap = false;
                    break;
                }
            }
            if (!everyFoldRetainsDiscoveryOverlap) continue;
            if (strongest == null
                    || candidate.alpha > strongest.alpha
                    || (candidate.alpha == strongest.alpha
                    && isBetterFrequencySeverityDiscoveryStructure(
                    candidate, strongest))) {
                strongest = candidate;
            }
        }
        return strongest;
    }

    private static boolean isBetterFrequencySeverityDiscoveryStructure(
            FrequencySeverityCandidateScore candidate,
            FrequencySeverityCandidateScore incumbent) {
        int compare = Double.compare(
                candidate.crossfit.sizing.epsilonAtMaxSamples,
                incumbent.crossfit.sizing.epsilonAtMaxSamples);
        if (compare != 0) return compare < 0;
        compare = Double.compare(candidate.crossfit.targetM2,
                incumbent.crossfit.targetM2);
        if (compare != 0) return compare < 0;
        compare = Double.compare(
                candidate.crossfit.shiftEssFraction,
                incumbent.crossfit.shiftEssFraction);
        if (compare != 0) return compare > 0;
        compare = Double.compare(candidate.tripleEtaScale,
                incumbent.tripleEtaScale);
        if (compare != 0) return compare < 0;
        return candidate.id.compareTo(incumbent.id) < 0;
    }

    private boolean frequencySeverityPilotPass(
            PackStarFrequencySeverityPAC.Interval interval,
            PackStarFrequencySeverityPAC.Sizing sizing) {
        return interval != null && sizing != null
                && interval.hasPositiveBulkLower()
                && sizing.reachableAtMax
                && sizing.epsilonAtMaxSamples
                <= targetEpsilon + 1.0e-12;
    }

    /**
     * Decide whether a failed, adaptation-only discovery pilot is promising
     * enough to receive more IID samples from the exact same frozen proposal.
     * This rule never changes eta and never consumes validation or final data.
     */
    static boolean shouldExtendFrequencySeverityDiscovery(
            int currentSamples, int maximumSamples,
            PackStarFrequencySeverityPAC.Interval interval,
            PackStarFrequencySeverityPAC.Sizing sizing,
            double targetEpsilon) {
        return currentSamples >= 2
                && currentSamples < maximumSamples
                && interval != null
                && interval.isFinite()
                && !interval.hasPositiveBulkLower()
                && sizing != null
                && sizing.reachableAtMax
                && Double.isFinite(sizing.epsilonAtMaxSamples)
                && Double.isFinite(targetEpsilon)
                && sizing.epsilonAtMaxSamples
                <= targetEpsilon + 1.0e-12;
    }

    private FrequencySeverityPilotEvaluation evaluateFrequencySeverityPilot(
            List<CCDResult> samples, EtaCorrections eta,
            double logMuTrain) {
        double[] logRelative = computeFrequencySeverityLogRelativeWeights(
                samples, eta, logMuTrain);
        PackStarFrequencySeverityPAC.Interval interval =
                PackStarFrequencySeverityPAC.evaluate(
                        logRelative,
                        frequencySeverityRelativeBoundKcal / RT,
                        frequencySeverityCap,
                        frequencySeverityEventDelta(),
                        frequencySeverityEventDelta());
        PackStarFrequencySeverityPAC.Moments sizingMoments =
                new PackStarFrequencySeverityPAC.Moments(
                        interval.bulkMean,
                        interval.bulkVariance * nstarInflate,
                        interval.tailProbabilityEmpirical);
        PackStarFrequencySeverityPAC.Sizing sizing =
                PackStarFrequencySeverityPAC.size(
                        sizingMoments, maxEstSamples,
                        configuredFrequencySeverityUnreachableSamples(),
                        targetEpsilon, frequencySeveritySizeSafety,
                        frequencySeverityCap,
                        frequencySeverityEventDelta(),
                        frequencySeverityEventDelta());
        return new FrequencySeverityPilotEvaluation(
                logRelative, interval, sizing);
    }

    private static boolean isBetterFrequencySeverityCandidate(
            FrequencySeverityCandidateScore candidate,
            FrequencySeverityCandidateScore incumbent) {
        FrequencySeverityCrossfitScore candidateScore = candidate.crossfit;
        FrequencySeverityCrossfitScore incumbentScore = incumbent.crossfit;
        if (candidateScore.sizing.finalSamples
                != incumbentScore.sizing.finalSamples) {
            return candidateScore.sizing.finalSamples
                    < incumbentScore.sizing.finalSamples;
        }
        int compare = Double.compare(
                candidateScore.sizing.epsilonAtMaxSamples,
                incumbentScore.sizing.epsilonAtMaxSamples);
        if (compare != 0) return compare < 0;
        compare = Double.compare(candidateScore.targetM2,
                incumbentScore.targetM2);
        if (compare != 0) return compare < 0;
        compare = Double.compare(
                Math.max(candidateScore.unaryUndertrainedAmplification,
                        candidateScore.pairUndertrainedAmplification),
                Math.max(incumbentScore.unaryUndertrainedAmplification,
                        incumbentScore.pairUndertrainedAmplification));
        if (compare != 0) return compare < 0;
        compare = Double.compare(candidate.alpha, incumbent.alpha);
        if (compare != 0) return compare < 0;
        compare = Double.compare(candidate.tripleEtaScale,
                incumbent.tripleEtaScale);
        if (compare != 0) return compare < 0;
        return candidate.id.compareTo(incumbent.id) < 0;
    }

    private int configuredFrequencySeverityUnreachableSamples() {
        int configured = getConfigInteger(
                "packstar.pac.unreachableCap",
                Math.min(maxEstSamples, Math.max(2 * pilotSamples, 400)));
        return Math.max(2, Math.min(maxEstSamples, configured));
    }

    /** Two confidence events for each of the three K* partition functions. */
    private double frequencySeverityEventDelta() {
        return delta / 6.0;
    }

    private double massAmplification(double candidateMass,
                                     double baselineMass) {
        if (baselineMass <= 1.0e-15) {
            return candidateMass <= 1.0e-15
                    ? 1.0 : Double.POSITIVE_INFINITY;
        }
        return candidateMass / baselineMass;
    }

    private FrequencySeverityRefitAudit auditFrequencySeverityRefit(
            List<CCDResult> samples, EtaCorrections eta,
            FrequencySeverityEtaCoverage coverage) {
        int n = samples.size();
        double[] logA = new double[n];
        double logSumA = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + computeEtaEnergy(sample.conf, eta);
            logA[i] = frequencySeveritySourceLogWeight(
                    sample, targetEnergy);
            logSumA = logAddExp(logSumA, logA[i]);
        }
        if (!Double.isFinite(logSumA)) {
            return new FrequencySeverityRefitAudit(
                    Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, false, false, false);
        }

        int baselineUnaryCount = 0;
        int baselinePairCount = 0;
        double candidateUnaryMass = 0.0;
        double candidatePairMass = 0.0;
        double sumWeightSquares = 0.0;
        for (int i = 0; i < n; i++) {
            double weight = Math.exp(logA[i] - logSumA);
            sumWeightSquares += weight * weight;
            int[] conf = samples.get(i).conf;
            if (touchesFrequencySeverityUndertrainedUnary(conf, coverage)) {
                baselineUnaryCount++;
                candidateUnaryMass += weight;
            }
            if (touchesFrequencySeverityUndertrainedPair(conf, coverage)) {
                baselinePairCount++;
                candidatePairMass += weight;
            }
        }
        double shiftEssFraction = sumWeightSquares > 0.0
                && Double.isFinite(sumWeightSquares)
                ? 1.0 / (n * sumWeightSquares) : Double.NaN;
        double baselineUnaryMass = (double) baselineUnaryCount / n;
        double baselinePairMass = (double) baselinePairCount / n;
        double unaryAmplification = massAmplification(
                candidateUnaryMass, baselineUnaryMass);
        double pairAmplification = massAmplification(
                candidatePairMass, baselinePairMass);
        boolean finite = Double.isFinite(shiftEssFraction)
                && Double.isFinite(candidateUnaryMass)
                && Double.isFinite(candidatePairMass)
                && Double.isFinite(unaryAmplification)
                && Double.isFinite(pairAmplification);
        boolean shiftPass = finite && shiftEssFraction + 1.0e-12
                >= frequencySeverityMinShiftEssFraction;
        boolean coveragePass = finite
                && unaryAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12
                && pairAmplification
                <= frequencySeverityMaxUndertrainedAmplification + 1.0e-12;
        return new FrequencySeverityRefitAudit(
                shiftEssFraction,
                baselineUnaryMass, candidateUnaryMass,
                unaryAmplification,
                baselinePairMass, candidatePairMass,
                pairAmplification,
                finite, shiftPass, coveragePass);
    }

    private double weightedLogMeanExp(double[] normalizedWeights,
                                      double[] logValues,
                                      double power) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < normalizedWeights.length; i++) {
            if (!(normalizedWeights[i] > 0.0)) continue;
            max = Math.max(max,
                    Math.log(normalizedWeights[i]) + power * logValues[i]);
        }
        if (!Double.isFinite(max)) return Double.NaN;
        double sum = 0.0;
        for (int i = 0; i < normalizedWeights.length; i++) {
            if (!(normalizedWeights[i] > 0.0)) continue;
            sum += Math.exp(Math.log(normalizedWeights[i])
                    + power * logValues[i] - max);
        }
        return max + Math.log(sum);
    }

    private static double logSumExp(double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) max = Math.max(max, value);
        if (!Double.isFinite(max)) return max;
        double sum = 0.0;
        for (double value : values) sum += Math.exp(value - max);
        return max + Math.log(sum);
    }

    /**
     * Unnormalized log(q_target/q_source).  The omitted -logZ_target is
     * common to every adaptation sample and therefore cancels under
     * self-normalization.  Retaining logZ_source is essential when batches
     * came from different sequential proposals.
     */
    static double frequencySeveritySourceLogWeight(
            double targetEnergy, double sourceEnergy,
            double sourceLogZ, double rt) {
        if (!Double.isFinite(targetEnergy)
                || !Double.isFinite(sourceEnergy)
                || !Double.isFinite(sourceLogZ)
                || !Double.isFinite(rt) || !(rt > 0.0)) {
            return Double.NaN;
        }
        return -(targetEnergy - sourceEnergy) / rt + sourceLogZ;
    }

    private double frequencySeveritySourceLogWeight(
            CCDResult sample, double targetEnergy) {
        if (Double.isFinite(sample.sourceProposalEnergy)
                && Double.isFinite(sample.sourceProposalLogZ)
                && sample.sourceProposalId != null
                && !sample.sourceProposalId.isEmpty()) {
            return frequencySeveritySourceLogWeight(
                    targetEnergy, sample.sourceProposalEnergy,
                    sample.sourceProposalLogZ, RT);
        }
        // Compatibility for non-production diagnostics produced before
        // proposal provenance was recorded.  In that case the source is q_m,
        // and its common log normalizer may be omitted.
        return -(targetEnergy - sample.eMin) / RT;
    }

    private double computeEtaEnergy(int[] conf, EtaCorrections eta) {
        double energy = eta.globalOffsetKcal;
        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            if (rc < 0) continue;
            energy += eta.getOneBodyEta(pos, rc);
            for (int pos2 = pos + 1; pos2 < conf.length; pos2++) {
                int rc2 = conf[pos2];
                if (rc2 < 0 || !interactionGraph.hasEdge(pos, pos2)) continue;
                energy += eta.getPairEta(pos, rc, pos2, rc2);
            }
        }
        if (eta.tripleEta != null && eta.tripleEtaScale != 0.0) {
            energy += eta.tripleEtaScale
                    * eta.tripleEta.scoreResidual(
                    conf, eta::getPairEta);
        }
        return energy;
    }

    private boolean touchesFrequencySeverityUndertrainedUnary(
            int[] conf, FrequencySeverityEtaCoverage coverage) {
        return touchesFrequencySeverityUnaryBelowCount(
                conf, coverage, frequencySeverityMinTrainCount);
    }

    private boolean touchesFrequencySeverityUnaryBelowCount(
            int[] conf, FrequencySeverityEtaCoverage coverage,
            int minimumCount) {
        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            if (rc >= 0 && coverage.unary(pos, rc) < minimumCount) {
                return true;
            }
        }
        return false;
    }

    private boolean touchesFrequencySeverityUndertrainedPair(
            int[] conf, FrequencySeverityEtaCoverage coverage) {
        return touchesFrequencySeverityPairBelowCount(
                conf, coverage, frequencySeverityMinTrainCount);
    }

    private boolean touchesFrequencySeverityPairBelowCount(
            int[] conf, FrequencySeverityEtaCoverage coverage,
            int minimumCount) {
        for (int pos1 = 0; pos1 < conf.length; pos1++) {
            int rc1 = conf[pos1];
            if (rc1 < 0) continue;
            for (int pos2 = pos1 + 1; pos2 < conf.length; pos2++) {
                int rc2 = conf[pos2];
                if (rc2 < 0 || !interactionGraph.hasEdge(pos1, pos2)) {
                    continue;
                }
                if (coverage.pair(pos1, rc1, pos2, rc2)
                        < minimumCount) return true;
            }
        }
        return false;
    }

    private EtaCorrections shrinkEta(EtaCorrections raw,
                                     FrequencySeverityEtaCoverage coverage,
                                     double unaryShrink,
                                     double pairShrink) {
        int[] numRCs = new int[raw.oneBody.length];
        for (int pos = 0; pos < raw.oneBody.length; pos++) {
            numRCs[pos] = raw.oneBody[pos].length;
        }
        EtaCorrections shrunk = new EtaCorrections(
                raw.oneBody.length, numRCs);
        for (int pos = 0; pos < raw.oneBody.length; pos++) {
            for (int rc = 0; rc < raw.oneBody[pos].length; rc++) {
                int rawCount = raw.oneBodyCounts[pos][rc];
                if (rawCount <= 0) continue;
                int contexts = coverage.unary(pos, rc);
                double factor = unaryShrink == 0.0
                        ? 1.0 : contexts / (contexts + unaryShrink);
                double value = factor * raw.getOneBodyEta(pos, rc);
                shrunk.oneBody[pos][rc] = value * rawCount;
                shrunk.oneBodyCounts[pos][rc] = rawCount;
                shrunk.oneBodyCount++;
            }
        }
        for (Map.Entry<Long, double[]> entry : raw.pairSums.entrySet()) {
            double[] stored = entry.getValue();
            if (!(stored[1] > 0.0)) continue;
            int contexts = coverage.pairDistinctContexts.getOrDefault(
                    entry.getKey(), 0);
            double factor = pairShrink == 0.0
                    ? 1.0 : contexts / (contexts + pairShrink);
            double mean = stored[0] / stored[1];
            shrunk.pairSums.put(entry.getKey(),
                    new double[]{factor * mean * stored[1], stored[1]});
        }
        shrunk.pairCount = shrunk.pairSums.size();
        shrunk.globalOffsetKcal = raw.globalOffsetKcal;
        return shrunk;
    }

    private void runAdaptiveFrequencySeverityPAC(long startTime, double logZRigid) {
        File artifactDir;
        try {
            artifactDir = prepareFrequencySeverityArtifactDirectory();
            writeFrequencySeverityProtocolArtifact(artifactDir, logZRigid);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC artifact initialization failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }

        double initialProposalLogZ = rootedRootEdge.getLogZUpper(0);
        if (!Double.isFinite(initialProposalLogZ)) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir, "initial q_m proposal normalizer is non-finite",
                    0, 0, 0, 0);
            failCertificate("AdaptiveFrequencySeverityPAC: initial q_m proposal"
                    + " normalizer is non-finite", startTime);
            return;
        }

        frequencySeverityTripleEta = null;
        frequencySeverityTripleEtaFallbackReason = null;
        if (frequencySeverityTripleEtaEnabled
                && frequencySeverityTripleEtaScale > 0.0) {
            try {
                long workload =
                        PackStarTripleEtaCorrections.countCliqueAssignments(
                                rcs, interactionGraph);
                System.out.println(
                        "[PACK*-adaptive-frequency-severity] triple-eta preflight:"
                                + " cliqueAssignments=" + workload
                                + ", cap="
                                + frequencySeverityTripleEtaMaxAssignments
                                + ", maximumScale="
                                + frequencySeverityTripleEtaScale
                                + ", scaleGrid="
                                + Arrays.toString(
                                frequencySeverityTripleEtaScaleGrid));
                frequencySeverityTripleEta =
                        PackStarTripleEtaCorrections.compute(
                                rcs, minimizingEcalc, interactionGraph,
                                branchMinimizingEmat,
                                frequencySeverityTripleEtaMaxAssignments);
                frequencySeverityTripleEta.writeArtifacts(
                        new File(artifactDir, "triple_eta_table.tsv"),
                        new File(artifactDir, "triple_eta_summary.tsv"));
            } catch (PackStarTripleEtaCorrections.AssignmentCapExceededException ex) {
                frequencySeverityTripleEta = null;
                frequencySeverityTripleEtaFallbackReason = ex.getMessage();
                writeFrequencySeverityTripleEtaFallbackArtifactQuietly(
                        artifactDir, ex);
                System.out.println(
                        "[PACK*-adaptive-frequency-severity] triple eta"
                                + " assignment cap exceeded; continuing with"
                                + " pair-only proposal: " + ex.getMessage());
            } catch (Exception ex) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "triple eta construction failed: "
                                + sanitizeTsv(ex.getMessage()),
                        0, 0, 0, 0);
                failCertificate(
                        "AdaptiveFrequencySeverityPAC: triple eta"
                                + " construction failed: "
                                + sanitizeTsv(ex.getMessage()),
                        startTime);
                return;
            }
        }

        System.out.println("[PACK*-adaptive-frequency-severity] enabled"
                + "; qMTrain=" + trainSamples
                + "; discoveryPerRound=" + frequencySeverityDiscoverySamples
                + "->" + frequencySeverityDiscoveryMaxSamples
                + "; maxRefits=" + frequencySeverityMaxRefits
                + "; validation=" + frequencySeverityValidationSamples
                + "; folds=" + frequencySeverityFolds
                + "; finalMinShiftESS=" + frequencySeverityMinShiftEssFraction
                + "; discoveryMinShiftESS="
                + frequencySeverityDiscoveryMinShiftEssFraction
                + "; candidateDP=0; proposalDP<="
                + (frequencySeverityMaxRefits + 1)
                + "; tripleEta="
                + (frequencySeverityTripleEta == null
                ? (frequencySeverityTripleEtaFallbackReason == null
                ? "disabled" : "pair-only-cap-fallback")
                : "crossfit-selected-signed-residual")
                + "; tripleEtaMaximumScale="
                + frequencySeverityTripleEtaScale
                + "; tripleEtaScaleGrid="
                + Arrays.toString(frequencySeverityTripleEtaScaleGrid)
                + "; selectedTripleLimit="
                + frequencySeverityTripleEtaMaxPositionTriples
                + "; fillEdgeLimit="
                + frequencySeverityTripleEtaMaxFillEdges
                + "; alpha=1-fixed; triple-model=pair-only-vs-selected");

        String trainStage = "adaptive-frequency-severity-qm-train";
        List<CCDResult> trainCCD = runParallelCCD(
                sampleConformationsFromDP(
                        trainSamples, stageRandom(trainStage)),
                branchMinimizingEmat, initialProposalLogZ,
                "adaptive-frequency-severity-qm");
        totalCCDCalls = trainCCD.size();
        boolean qMLowerBoundValid = validateObservedLowerBound(
                trainCCD, trainStage);
        if (trainCCD.size() != trainSamples || !qMLowerBoundValid) {
            if (!qMLowerBoundValid) {
                writeFrequencySeverityLowerBoundViolationAuditQuietly(
                        artifactDir, trainStage, trainCCD);
            }
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir, "invalid or incomplete q_m training sample",
                    trainCCD.size(), 0, 0, 0);
            failCertificate("AdaptiveFrequencySeverityPAC: invalid or incomplete"
                    + " q_m training sample", startTime);
            return;
        }

        List<CCDResult> fitPool = new ArrayList<>(trainCCD);
        List<CCDResult> history = Collections.emptyList();
        List<CCDResult> scoringBatch = trainCCD;
        int discoveryCCDCalls = 0;
        int proposalDpSweeps = 0;

        FrequencySeverityCandidateScore selected = null;
        FrequencySeverityEtaTraining selectedTraining = null;
        EtaCorrections selectedEta = null;
        EnergyMatrix selectedEmat = null;
        double selectedLogZ = Double.NaN;
        double selectedLogMu = Double.NaN;
        FrequencySeverityPilotEvaluation preliminaryEvaluation = null;
        List<CCDResult> preliminarySamples = null;
        int selectedRound = -1;
        boolean selectedWasDiscoveryProbe = false;

        for (int round = 0; round <= frequencySeverityMaxRefits; round++) {
            String artifactPrefix = String.format(
                    Locale.ROOT, "adaptive_eta_round_%02d", round);
            FrequencySeverityEtaTraining training;
            List<FrequencySeverityCandidateScore> candidates;
            FrequencySeverityCandidateScore finalEligible;
            FrequencySeverityCandidateScore discoveryProbe;
            try {
                training = extractFrequencySeverityEtaTraining(
                        history, scoringBatch);
                candidates = scoreFrequencySeverityEtaCandidates(
                        scoringBatch, training);
                finalEligible = selectBestFrequencySeverityCandidate(candidates);
                discoveryProbe = finalEligible == null
                        ? selectFrequencySeverityDiscoveryProbe(candidates)
                        : null;
                FrequencySeverityCandidateScore roundSelection =
                        finalEligible != null ? finalEligible : discoveryProbe;
                writeFrequencySeverityTrainingArtifacts(
                        artifactDir, artifactPrefix, fitPool, training);
                writeFrequencySeverityCandidateArtifact(
                        artifactDir, artifactPrefix, candidates,
                        roundSelection == null ? null : roundSelection.id);
                writeFrequencySeverityCrossfitArtifact(
                        artifactDir, artifactPrefix, candidates,
                        roundSelection == null ? null : roundSelection.id);
                writeFrequencySeverityRefitArtifact(
                        artifactDir, artifactPrefix, candidates,
                        roundSelection == null ? null : roundSelection.id);
                if (roundSelection != null
                        && hasTripleEta(roundSelection.fullEta)) {
                    roundSelection.fullEta.tripleEta.writeArtifacts(
                            new File(artifactDir,
                                    artifactPrefix
                                            + "_selected_triple_eta_table.tsv"),
                            new File(artifactDir,
                                    artifactPrefix
                                            + "_selected_triple_eta_summary.tsv"));
                }
            } catch (Exception ex) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "round " + round + " eta fit/score/artifact failed: "
                                + sanitizeTsv(ex.getMessage()),
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC round " + round
                        + " eta fit/score/artifact failed: "
                        + sanitizeTsv(ex.getMessage()), startTime);
                return;
            }

            FrequencySeverityCandidateScore roundSelection =
                    finalEligible != null ? finalEligible : discoveryProbe;
            boolean roundUsesDiscoveryProbe = finalEligible == null;
            if (roundSelection == null) {
                writeFrequencySeverityCoverageArtifactQuietly(
                        artifactDir, training.all,
                        training.allCoverage, null);
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "no final candidate or overlap-safe unshrunk discovery probe",
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC: no final candidate"
                        + " and no discovery probe retained source-shift ESS >= "
                        + frequencySeverityDiscoveryMinShiftEssFraction,
                        startTime);
                return;
            }

            EtaCorrections roundEta = roundSelection.fullEta;
            double roundLogMu = computeFrequencySeverityTrainingLogMu(
                    fitPool, roundEta);
            if (!Double.isFinite(roundLogMu)) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir, "source-aware training gauge is non-finite",
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC: source-aware"
                        + " training gauge is non-finite", startTime);
                return;
            }

            EnergyMatrix roundEmat;
            double roundLogZ;
            if (round == 0 && roundSelection.alpha == 0.0
                    && !hasTripleEta(roundEta)) {
                roundEmat = branchMinimizingEmat;
                roundLogZ = initialProposalLogZ;
            } else {
                roundEmat = buildCorrectedEmat(roundEta);
                CorrectedDPResult correctedDP = recomputeDP(
                        roundEmat, roundEta);
                proposalDpSweeps++;
                roundLogZ = correctedDP.logZCorrected;
            }
            if (!Double.isFinite(roundLogZ)) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "round " + round + " proposal DP normalizer is non-finite",
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC round " + round
                        + ": proposal DP normalizer is non-finite", startTime);
                return;
            }

            String discoveryStage = String.format(
                    Locale.ROOT,
                    "adaptive-frequency-severity-discovery-%02d", round);
            List<CCDResult> discoverySamples = new ArrayList<>(runParallelCCD(
                    sampleConformationsFromDP(
                            frequencySeverityDiscoverySamples,
                            stageRandom(discoveryStage)),
                    roundEmat, roundLogZ,
                    "adaptive-frequency-severity-proposal-" + round));
            discoveryCCDCalls += discoverySamples.size();
            totalCCDCalls = trainCCD.size() + discoveryCCDCalls;
            if (discoverySamples.size() != frequencySeverityDiscoverySamples
                    || !validateObservedLowerBound(
                    discoverySamples, discoveryStage)) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "invalid or incomplete discovery sample at round " + round,
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC round " + round
                        + ": invalid or incomplete discovery sample", startTime);
                return;
            }

            int initialDiscoveryCount = discoverySamples.size();
            int extensionDiscoveryCount = 0;
            FrequencySeverityPilotEvaluation initialEvaluation;
            try {
                initialEvaluation = evaluateFrequencySeverityPilot(
                        discoverySamples, roundEta, roundLogMu);
            } catch (RuntimeException ex) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "round " + round + " initial discovery evaluation failed: "
                                + sanitizeTsv(ex.getMessage()),
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC round " + round
                        + " initial discovery evaluation failed: "
                        + sanitizeTsv(ex.getMessage()), startTime);
                return;
            }

            boolean discoveryExtended =
                    shouldExtendFrequencySeverityDiscovery(
                            discoverySamples.size(),
                            frequencySeverityDiscoveryMaxSamples,
                            initialEvaluation.interval,
                            initialEvaluation.sizing,
                            targetEpsilon);
            if (discoveryExtended) {
                String initialStage = discoveryStage + "-initial";
                try {
                    writeFrequencySeverityStageArtifact(
                            artifactDir, initialStage,
                            initialEvaluation.interval, null,
                            initialEvaluation.sizing, roundLogMu, roundLogZ);
                    writeFrequencySeveritySamplesArtifact(
                            artifactDir, initialStage,
                            discoverySamples, roundEta, roundLogMu);
                } catch (RuntimeException ex) {
                    writeFrequencySeverityFailureArtifactQuietly(
                            artifactDir,
                            "round " + round
                                    + " initial discovery artifact failed: "
                                    + sanitizeTsv(ex.getMessage()),
                            trainCCD.size(), discoveryCCDCalls, 0, 0);
                    failCertificate("AdaptiveFrequencySeverityPAC round " + round
                            + " initial discovery artifact failed: "
                            + sanitizeTsv(ex.getMessage()), startTime);
                    return;
                }

                int requestedExtension = frequencySeverityDiscoveryMaxSamples
                        - discoverySamples.size();
                String extensionStage = discoveryStage + "-extension";
                System.out.println("[PACK*-adaptive-frequency-severity] round="
                        + round
                        + ", initialDiscovery=" + discoverySamples.size()
                        + ", initialTail=" + initialEvaluation.interval.tailCount
                        + ", initialBulkLower="
                        + String.format(Locale.ROOT, "%.9g",
                        initialEvaluation.interval.bulkLower)
                        + ", initialEpsilonAtMax="
                        + String.format(Locale.ROOT, "%.9f",
                        initialEvaluation.sizing.epsilonAtMaxSamples)
                        + ", action=extend-same-frozen-proposal-to-"
                        + frequencySeverityDiscoveryMaxSamples);
                List<CCDResult> extensionSamples = runParallelCCD(
                        sampleConformationsFromDP(
                                requestedExtension,
                                stageRandom(extensionStage)),
                        roundEmat, roundLogZ,
                        "adaptive-frequency-severity-proposal-" + round);
                extensionDiscoveryCount = extensionSamples.size();
                discoveryCCDCalls += extensionDiscoveryCount;
                totalCCDCalls = trainCCD.size() + discoveryCCDCalls;
                if (extensionDiscoveryCount != requestedExtension
                        || !validateObservedLowerBound(
                        extensionSamples, extensionStage)) {
                    writeFrequencySeverityFailureArtifactQuietly(
                            artifactDir,
                            "invalid or incomplete same-proposal discovery"
                                    + " extension at round " + round,
                            trainCCD.size(), discoveryCCDCalls, 0, 0);
                    failCertificate("AdaptiveFrequencySeverityPAC round " + round
                            + ": invalid or incomplete same-proposal discovery"
                            + " extension", startTime);
                    return;
                }
                try {
                    writeFrequencySeveritySamplesArtifact(
                            artifactDir, extensionStage,
                            extensionSamples, roundEta, roundLogMu);
                } catch (RuntimeException ex) {
                    writeFrequencySeverityFailureArtifactQuietly(
                            artifactDir,
                            "round " + round
                                    + " discovery extension artifact failed: "
                                    + sanitizeTsv(ex.getMessage()),
                            trainCCD.size(), discoveryCCDCalls, 0, 0);
                    failCertificate("AdaptiveFrequencySeverityPAC round " + round
                            + " discovery extension artifact failed: "
                            + sanitizeTsv(ex.getMessage()), startTime);
                    return;
                }
                discoverySamples.addAll(extensionSamples);
            }

            FrequencySeverityPilotEvaluation roundEvaluation;
            boolean preliminaryPass;
            try {
                roundEvaluation = discoveryExtended
                        ? evaluateFrequencySeverityPilot(
                        discoverySamples, roundEta, roundLogMu)
                        : initialEvaluation;
                preliminaryPass = frequencySeverityPilotPass(
                        roundEvaluation.interval, roundEvaluation.sizing);
                writeFrequencySeverityStageArtifact(
                        artifactDir, discoveryStage,
                        roundEvaluation.interval, null,
                        roundEvaluation.sizing, roundLogMu, roundLogZ);
                writeFrequencySeveritySamplesArtifact(
                        artifactDir, discoveryStage,
                        discoverySamples, roundEta, roundLogMu);
                writeAdaptiveFrequencySeverityRoundArtifact(
                        artifactDir, artifactPrefix, round,
                        history.size(), scoringBatch.size(), fitPool.size(),
                        roundSelection, roundUsesDiscoveryProbe,
                        roundLogMu, roundLogZ, proposalDpSweeps,
                        discoverySamples, initialDiscoveryCount,
                        extensionDiscoveryCount,
                        frequencySeverityDiscoveryMaxSamples,
                        training.allCoverage,
                        roundEvaluation, preliminaryPass);
            } catch (RuntimeException ex) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "round " + round + " discovery evaluation/artifact failed: "
                                + sanitizeTsv(ex.getMessage()),
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC round " + round
                        + " discovery evaluation/artifact failed: "
                        + sanitizeTsv(ex.getMessage()), startTime);
                return;
            }

            System.out.println("[PACK*-adaptive-frequency-severity] round=" + round
                    + ", role=" + (roundUsesDiscoveryProbe
                    ? "discovery-probe" : "final-eligible")
                    + ", candidate=" + roundSelection.id
                    + ", fitPool=" + fitPool.size()
                    + ", sameProposalExtended=" + discoveryExtended
                    + ", tail=" + roundEvaluation.interval.tailCount
                    + "/" + roundEvaluation.interval.sampleCount
                    + ", bulkLower="
                    + String.format(Locale.ROOT, "%.9g",
                    roundEvaluation.interval.bulkLower)
                    + ", reachable=" + roundEvaluation.sizing.reachableAtMax
                    + ", epsilonAtMax="
                    + String.format(Locale.ROOT, "%.9f",
                    roundEvaluation.sizing.epsilonAtMaxSamples)
                    + ", action=" + (preliminaryPass
                    ? "validate" : round < frequencySeverityMaxRefits
                    ? "merge-and-refit" : "abort"));

            if (preliminaryPass) {
                selected = roundSelection;
                selectedTraining = training;
                selectedEta = roundEta;
                selectedEmat = roundEmat;
                selectedLogZ = roundLogZ;
                selectedLogMu = roundLogMu;
                preliminaryEvaluation = roundEvaluation;
                preliminarySamples = discoverySamples;
                selectedRound = round;
                selectedWasDiscoveryProbe = roundUsesDiscoveryProbe;
                break;
            }

            if (round == frequencySeverityMaxRefits) {
                writeFrequencySeverityFailureArtifactQuietly(
                        artifactDir,
                        "discovery pilot remained unreachable or bulk-collapsed at refit cap",
                        trainCCD.size(), discoveryCCDCalls, 0, 0);
                failCertificate("AdaptiveFrequencySeverityPAC: discovery pilot"
                        + " remained unreachable or bulk-collapsed after "
                        + frequencySeverityMaxRefits + " refits", startTime);
                return;
            }

            history = new ArrayList<>(fitPool);
            scoringBatch = discoverySamples;
            fitPool.addAll(discoverySamples);
        }

        if (selected == null || selectedTraining == null
                || selectedEta == null || selectedEmat == null
                || preliminaryEvaluation == null
                || preliminarySamples == null) {
            failCertificate("AdaptiveFrequencySeverityPAC: internal proposal-selection"
                    + " state is incomplete", startTime);
            return;
        }

        String validationStage = "adaptive-frequency-severity-validation";
        List<CCDResult> validationCCD = runParallelCCD(
                sampleConformationsFromDP(
                        frequencySeverityValidationSamples,
                        stageRandom(validationStage)),
                selectedEmat, selectedLogZ,
                "adaptive-frequency-severity-proposal-" + selectedRound);
        totalCCDCalls = trainCCD.size() + discoveryCCDCalls
                + validationCCD.size();
        if (validationCCD.size() != frequencySeverityValidationSamples
                || !validateObservedLowerBound(
                validationCCD, validationStage)) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir, "invalid or incomplete independent validation",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), 0, 0);
            failCertificate("AdaptiveFrequencySeverityPAC: invalid or incomplete"
                    + " independent validation", startTime);
            return;
        }

        FrequencySeverityPilotEvaluation validation;
        try {
            validation = evaluateFrequencySeverityPilot(
                    validationCCD, selectedEta, selectedLogMu);
            writeFrequencySeverityStageArtifact(
                    artifactDir, "validation", validation.interval,
                    null, validation.sizing,
                    selectedLogMu, selectedLogZ);
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "validation", validationCCD,
                    selectedEta, selectedLogMu);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC validation statistics/artifact"
                    + " failed: " + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (!frequencySeverityPilotPass(
                validation.interval, validation.sizing)) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir,
                    "independent 400-sample validation was unreachable or bulk-collapsed",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), 0, 0);
            failCertificate("AdaptiveFrequencySeverityPAC: independent validation"
                    + " did not confirm a reachable positive-bulk proposal", startTime);
            return;
        }

        try {
            writeFrequencySeverityCoverageArtifact(
                    artifactDir, selectedTraining.all,
                    selectedTraining.allCoverage, selectedEta);
            writeFrequencySeverityWinnerArtifact(
                    artifactDir, selected, selected.refitAudit,
                    selectedLogMu, proposalDpSweeps);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC selected-proposal artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }

        System.out.println("[PACK*-adaptive-frequency-severity] frozen proposal:"
                + " round=" + selectedRound
                + ", candidate=" + selected.id
                + ", selectedVia=" + (selectedWasDiscoveryProbe
                ? "on-policy-validated-discovery-probe" : "final-eligible-grid")
                + ", fitPool=" + fitPool.size()
                + ", discoveryCCD=" + discoveryCCDCalls
                + ", validationTail=" + validation.interval.tailCount
                + "/" + validation.interval.sampleCount
                + ", validationEpsilonAtMax="
                + String.format(Locale.ROOT, "%.9f",
                validation.sizing.epsilonAtMaxSamples)
                + ", finalN=" + validation.sizing.finalSamples
                + ", proposalDpSweeps=" + proposalDpSweeps);

        String monitorStage = "adaptive-frequency-severity-monitor";
        List<CCDResult> monitorCCD = runParallelCCD(
                sampleConformationsFromDP(
                        monitorSamples, stageRandom(monitorStage)),
                selectedEmat, selectedLogZ, monitorStage);
        totalCCDCalls = trainCCD.size() + discoveryCCDCalls
                + validationCCD.size() + monitorCCD.size();
        if (monitorCCD.size() != monitorSamples
                || !validateObservedLowerBound(monitorCCD, monitorStage)) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir, "invalid or incomplete independent monitor",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), monitorCCD.size(), 0);
            failCertificate("AdaptiveFrequencySeverityPAC: independent monitor"
                    + " observed an implementation-integrity failure", startTime);
            return;
        }
        double[] monitorLogR = computeFrequencySeverityLogRelativeWeights(
                monitorCCD, selectedEta, selectedLogMu);
        PackStarFrequencySeverityPAC.Interval monitorInterval;
        PackStarFrequencySeverityPAC.SeverityTest severityTest;
        try {
            monitorInterval = PackStarFrequencySeverityPAC.evaluate(
                    monitorLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            severityTest = PackStarFrequencySeverityPAC.testConditionalSeverity(
                    monitorLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityTestAlpha);
            writeFrequencySeverityStageArtifact(
                    artifactDir, "monitor", monitorInterval,
                    severityTest, null, selectedLogMu, selectedLogZ);
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "monitor", monitorCCD,
                    selectedEta, selectedLogMu);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC monitor statistics/artifact"
                    + " failed: " + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (severityTest.rejected) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir,
                    "independent monitor rejected the external severity premise",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), monitorCCD.size(), 0);
            failCertificate("AdaptiveFrequencySeverityPAC: independent monitor"
                    + " rejected severity premise "
                    + frequencySeverityPremiseId, startTime);
            return;
        }

        String finalStage = "adaptive-frequency-severity-final";
        List<CCDResult> finalCCD = runParallelCCD(
                sampleConformationsFromDP(
                        validation.sizing.finalSamples,
                        stageRandom(finalStage)),
                selectedEmat, selectedLogZ, finalStage);
        totalCCDCalls = trainCCD.size() + discoveryCCDCalls
                + validationCCD.size() + monitorCCD.size()
                + finalCCD.size();
        if (finalCCD.size() != validation.sizing.finalSamples
                || !validateObservedLowerBound(finalCCD, finalStage)) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir, "invalid or incomplete fresh final sample",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), monitorCCD.size(),
                    finalCCD.size());
            failCertificate("AdaptiveFrequencySeverityPAC: invalid or incomplete"
                    + " fresh final sample", startTime);
            return;
        }
        double[] finalLogR = computeFrequencySeverityLogRelativeWeights(
                finalCCD, selectedEta, selectedLogMu);
        PackStarFrequencySeverityPAC.Interval finalInterval;
        try {
            finalInterval = PackStarFrequencySeverityPAC.evaluate(
                    finalLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "final", finalCCD,
                    selectedEta, selectedLogMu);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC final statistics failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (frequencySeverityCap == 0.0 && finalInterval.tailCount > 0) {
            writeFrequencySeverityFailureArtifactQuietly(
                    artifactDir,
                    "final sample logically violates zero severity premise",
                    trainCCD.size(), discoveryCCDCalls
                            + validationCCD.size(), monitorCCD.size(),
                    finalCCD.size());
            failCertificate("AdaptiveFrequencySeverityPAC: final sample logically"
                    + " violates the zero severity premise", startTime);
            return;
        }

        double logScale = selectedLogZ + selectedLogMu
                + frequencySeverityRelativeBoundKcal / RT;
        logZLowerPAC = finalInterval.normalizedMeanLower > 0.0
                ? logScale + Math.log(finalInterval.normalizedMeanLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = finalInterval.normalizedMeanUpper > 0.0
                ? logScale + Math.log(finalInterval.normalizedMeanUpper)
                : Double.NaN;
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);
        setFrequencySeverityDiagnosticStatistics(
                finalCCD, selectedEta, finalLogR, finalInterval);
        certificateValid = finalInterval.isFinite()
                && isValidCertificate(zLower, zUpper, epsilon);
        certificateFailureReason = certificateValid ? ""
                : "AdaptiveFrequencySeverityPAC final interval failed validation";
        try {
            writeFrequencySeverityStageArtifact(
                    artifactDir, "final", finalInterval,
                    null, validation.sizing,
                    selectedLogMu, selectedLogZ);
            writeFrequencySeverityRunSummaryArtifact(
                    artifactDir, selected, validation.sizing,
                    finalInterval, trainCCD.size(),
                    discoveryCCDCalls, validationCCD.size(),
                    monitorCCD.size(), finalCCD.size(),
                    proposalDpSweeps, certificateValid);
        } catch (RuntimeException ex) {
            failCertificate("AdaptiveFrequencySeverityPAC final artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (!certificateValid) {
            printFinalSummary(startTime);
            return;
        }
        System.out.println("[PACK*-adaptive-frequency-severity] issued"
                + " assumption-conditional interval: premiseId="
                + frequencySeverityPremiseId
                + ", tail=" + finalInterval.tailCount
                + "/" + finalInterval.sampleCount
                + ", epsilon="
                + String.format(Locale.ROOT, "%.9f", epsilon)
                + ", targetReached="
                + (epsilon <= targetEpsilon + 1.0e-12));
        printFinalSummary(startTime);
    }

    /** Retained temporarily for artifact comparison with the failed one-shot replay. */
    private void runOneShotFrequencySeverityPACLegacy(long startTime, double logZRigid) {
        File artifactDir;
        try {
            artifactDir = prepareFrequencySeverityArtifactDirectory();
            writeFrequencySeverityProtocolArtifact(artifactDir, logZRigid);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 artifact initialization failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }

        System.out.println("[PACK*-frequency-severity] production frequency/severity route enabled"
                + "; train=" + trainSamples
                + "; folds=" + frequencySeverityFolds
                + "; shrinkCandidates=" + frequencySeverityShrinkGrid.length
                + "; alphaGrid=" + Arrays.toString(frequencySeverityAlphaGrid)
                + "; minShiftESSFraction=" + frequencySeverityMinShiftEssFraction
                + "; B_rel_kcal=" + frequencySeverityRelativeBoundKcal
                + "; S0=" + frequencySeverityCap
                + "; premiseId=" + frequencySeverityPremiseId
                + "; candidateDP=0; candidateCCD=0; winnerDP<=1"
                + "; legacyRefitBypassed=true");

        // Stage A: the only proposal-training CCD batch, drawn from q_m.
        List<CCDResult> trainCCD = runParallelCCD(
                sampleConformationsFromDP(
                        trainSamples, stageRandom("frequency-severity-v21-train")));
        totalCCDCalls = trainCCD.size();
        if (trainCCD.size() != trainSamples
                || !validateObservedLowerBound(
                trainCCD, "frequency-severity-v21-train")) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "invalid or incomplete q_m training sample",
                    trainCCD.size(), 0, 0, 0);
            failCertificate("frequency/severity PAC V2.1: invalid or incomplete q_m training sample",
                    startTime);
            return;
        }

        FrequencySeverityEtaTraining training;
        try {
            training = extractFrequencySeverityEtaTraining(trainCCD);
            writeFrequencySeverityTrainingArtifacts(
                    artifactDir, trainCCD, training);
        } catch (RuntimeException ex) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "eta training extraction/artifact failed: "
                            + sanitizeTsv(ex.getMessage()),
                    trainCCD.size(), 0, 0, 0);
            failCertificate("frequency/severity PAC V2.1 eta training extraction/artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        List<FrequencySeverityCandidateScore> candidates =
                scoreFrequencySeverityEtaCandidates(trainCCD, training);
        FrequencySeverityCandidateScore winner = selectBestFrequencySeverityCandidate(candidates);
        try {
            writeFrequencySeverityCandidateArtifact(artifactDir, candidates,
                    winner == null ? null : winner.id);
            writeFrequencySeverityCrossfitArtifact(artifactDir, candidates,
                    winner == null ? null : winner.id);
            writeFrequencySeverityRefitArtifact(artifactDir, candidates,
                    winner == null ? null : winner.id);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 candidate artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (winner == null) {
            writeFrequencySeverityCoverageArtifactQuietly(
                    artifactDir, training.all, training.allCoverage, null);
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "no eta-only candidate passed cross-fit and full-refit gates",
                    trainCCD.size(), 0, 0, 0);
            failCertificate("frequency/severity PAC V2.1: no eta-only candidate passed"
                    + " cross-fit shift/bulk/reachability and full-refit"
                    + " shift/sparse-cell gates;"
                    + " defensive active tier not enabled", startTime);
            return;
        }

        EtaCorrections eta = winner.fullEta;
        FrequencySeverityRefitAudit refitAudit = winner.refitAudit;
        int winnerDpSweeps = winner.alpha == 0.0
                && !hasTripleEta(winner.fullEta) ? 0 : 1;
        double logMuTrain = computeFrequencySeverityTrainingLogMu(trainCCD, eta);
        if (!Double.isFinite(logMuTrain)) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "training-relative gauge is non-finite",
                    trainCCD.size(), 0, 0, 0);
            failCertificate("frequency/severity PAC V2.1: training-relative gauge is non-finite",
                    startTime);
            return;
        }
        try {
            writeFrequencySeverityCoverageArtifact(
                    artifactDir, training.all,
                    training.allCoverage, eta);
            writeFrequencySeverityWinnerArtifact(
                    artifactDir, winner,
                    refitAudit, logMuTrain, winnerDpSweeps);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 coverage artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }

        // At most one candidate-corrected DP solve.  The eta-zero winner reuses
        // the still-loaded initial q_m tables.  The relative gauge is an
        // analytic scale and never triggers another DP sweep.
        EnergyMatrix correctedEmat;
        double logZCorrected;
        if (winnerDpSweeps == 0) {
            correctedEmat = branchMinimizingEmat;
            logZCorrected = rootedRootEdge.getLogZUpper(0);
        } else {
            correctedEmat = buildCorrectedEmat(eta);
            CorrectedDPResult correctedDP = recomputeDP(
                    correctedEmat, eta);
            logZCorrected = correctedDP.logZCorrected;
        }
        if (!Double.isFinite(logZCorrected)) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "winner corrected DP normalizer is non-finite",
                    trainCCD.size(), 0, 0, 0);
            failCertificate("frequency/severity PAC V2.1: winner corrected DP normalizer is non-finite",
                    startTime);
            return;
        }

        // Stage B: one on-policy pilot.  It sizes/vetoes only; it never refits
        // eta, changes the gauge, or competes with another proposal.
        String pilotStage = "frequency-severity-v21-pilot";
        List<CCDResult> pilotCCD = runParallelCCD(
                sampleConformationsFromDP(
                        pilotSamples, stageRandom(pilotStage)),
                correctedEmat, logZCorrected, pilotStage);
        totalCCDCalls = trainCCD.size() + pilotCCD.size();
        if (pilotCCD.size() != pilotSamples
                || !validateObservedLowerBound(pilotCCD, pilotStage)) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "invalid or incomplete winner pilot",
                    trainCCD.size(), pilotCCD.size(), 0, 0);
            failCertificate("frequency/severity PAC V2.1: invalid or incomplete winner pilot",
                    startTime);
            return;
        }
        double[] pilotLogR = computeFrequencySeverityLogRelativeWeights(
                pilotCCD, eta, logMuTrain);
        PackStarFrequencySeverityPAC.Interval pilotInterval;
        PackStarFrequencySeverityPAC.Sizing sizing;
        try {
            pilotInterval = PackStarFrequencySeverityPAC.evaluate(
                    pilotLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            PackStarFrequencySeverityPAC.Moments sizingMoments =
                    new PackStarFrequencySeverityPAC.Moments(
                            pilotInterval.bulkMean,
                            pilotInterval.bulkVariance * nstarInflate,
                            pilotInterval.tailProbabilityEmpirical);
            sizing = PackStarFrequencySeverityPAC.size(
                    sizingMoments, maxEstSamples,
                    configuredFrequencySeverityUnreachableSamples(),
                    targetEpsilon, frequencySeveritySizeSafety,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            writeFrequencySeverityStageArtifact(artifactDir, "pilot", pilotInterval,
                    null, sizing, logMuTrain, logZCorrected);
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "pilot", pilotCCD,
                    eta, logMuTrain);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 pilot statistics/artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (!pilotInterval.hasPositiveBulkLower()) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "pilot bounded-bulk lower endpoint is zero",
                    trainCCD.size(), pilotCCD.size(), 0, 0);
            failCertificate("frequency/severity PAC V2.1: pilot bounded-bulk lower endpoint is zero;"
                    + " severity cannot repair bulk collapse", startTime);
            return;
        }

        System.out.println("[PACK*-frequency-severity] frozen winner: id=" + winner.id
                + ", logZ_eta="
                + String.format(Locale.ROOT, "%.9f", logZCorrected)
                + ", logMuTrain="
                + String.format(Locale.ROOT, "%.9f", logMuTrain)
                + ", pilotTail=" + pilotInterval.tailCount
                + "/" + pilotInterval.sampleCount
                + ", pilotBulkLower="
                + String.format(Locale.ROOT, "%.9g",
                pilotInterval.bulkLower)
                + ", reachableAtMax=" + sizing.reachableAtMax
                + ", epsilonAtMax="
                + String.format(Locale.ROOT, "%.9f",
                sizing.epsilonAtMaxSamples)
                + ", finalN=" + sizing.finalSamples);

        // Stage C: independent falsification monitor.  A tail observation is
        // expected frequency/severity tail data, not a pointwise-B abort.  Only a valid e-test
        // rejection of S0 or an implementation-integrity failure vetoes.
        String monitorStage = "frequency-severity-v21-monitor";
        List<CCDResult> monitorCCD = runParallelCCD(
                sampleConformationsFromDP(
                        monitorSamples, stageRandom(monitorStage)),
                correctedEmat, logZCorrected, monitorStage);
        totalCCDCalls = trainCCD.size() + pilotCCD.size()
                + monitorCCD.size();
        if (monitorCCD.size() != monitorSamples
                || !validateObservedLowerBound(monitorCCD, monitorStage)) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "invalid or incomplete independent monitor",
                    trainCCD.size(), pilotCCD.size(),
                    monitorCCD.size(), 0);
            failCertificate("frequency/severity PAC V2.1: independent monitor observed an"
                    + " implementation-integrity failure", startTime);
            return;
        }
        double[] monitorLogR = computeFrequencySeverityLogRelativeWeights(
                monitorCCD, eta, logMuTrain);
        PackStarFrequencySeverityPAC.Interval monitorInterval;
        PackStarFrequencySeverityPAC.SeverityTest severityTest;
        try {
            monitorInterval = PackStarFrequencySeverityPAC.evaluate(
                    monitorLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            severityTest =
                    PackStarFrequencySeverityPAC.testConditionalSeverity(
                            monitorLogR, frequencySeverityRelativeBoundKcal / RT,
                            frequencySeverityCap, frequencySeverityTestAlpha);
            writeFrequencySeverityStageArtifact(artifactDir, "monitor", monitorInterval,
                    severityTest, null, logMuTrain, logZCorrected);
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "monitor", monitorCCD,
                    eta, logMuTrain);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 monitor statistics/artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (severityTest.rejected) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "independent monitor rejected the external severity premise",
                    trainCCD.size(), pilotCCD.size(),
                    monitorCCD.size(), 0);
            failCertificate("frequency/severity PAC V2.1: independent monitor rejected severity"
                    + " premise " + frequencySeverityPremiseId, startTime);
            return;
        }

        // Stage D: the only sample entering the issued frequency/severity confidence interval.
        String finalStage = "frequency-severity-v21-final";
        List<CCDResult> finalCCD = runParallelCCD(
                sampleConformationsFromDP(
                        sizing.finalSamples, stageRandom(finalStage)),
                correctedEmat, logZCorrected, finalStage);
        totalCCDCalls = trainCCD.size() + pilotCCD.size()
                + monitorCCD.size() + finalCCD.size();
        if (finalCCD.size() != sizing.finalSamples
                || !validateObservedLowerBound(finalCCD, finalStage)) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "invalid or incomplete fresh final sample",
                    trainCCD.size(), pilotCCD.size(),
                    monitorCCD.size(), finalCCD.size());
            failCertificate("frequency/severity PAC V2.1: invalid or incomplete fresh final sample",
                    startTime);
            return;
        }
        double[] finalLogR = computeFrequencySeverityLogRelativeWeights(
                finalCCD, eta, logMuTrain);
        PackStarFrequencySeverityPAC.Interval finalInterval;
        try {
            finalInterval = PackStarFrequencySeverityPAC.evaluate(
                    finalLogR, frequencySeverityRelativeBoundKcal / RT,
                    frequencySeverityCap, frequencySeverityEventDelta(),
                    frequencySeverityEventDelta());
            writeFrequencySeveritySamplesArtifact(
                    artifactDir, "final", finalCCD,
                    eta, logMuTrain);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 final statistics failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (frequencySeverityCap == 0.0 && finalInterval.tailCount > 0) {
            writeFrequencySeverityFailureArtifactQuietly(artifactDir,
                    "final sample logically violates zero severity premise",
                    trainCCD.size(), pilotCCD.size(),
                    monitorCCD.size(), finalCCD.size());
            failCertificate("frequency/severity PAC V2.1: final sample logically violates the"
                    + " zero severity premise", startTime);
            return;
        }

        double logScale = logZCorrected + logMuTrain
                + frequencySeverityRelativeBoundKcal / RT;
        logZLowerPAC = finalInterval.normalizedMeanLower > 0.0
                ? logScale + Math.log(finalInterval.normalizedMeanLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = finalInterval.normalizedMeanUpper > 0.0
                ? logScale + Math.log(finalInterval.normalizedMeanUpper)
                : Double.NaN;
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);
        setFrequencySeverityDiagnosticStatistics(finalCCD, eta, finalLogR,
                finalInterval);
        certificateValid = finalInterval.isFinite()
                && isValidCertificate(zLower, zUpper, epsilon);
        certificateFailureReason = certificateValid
                ? "" : "frequency/severity PAC V2.1 final interval failed validation";
        try {
            writeFrequencySeverityStageArtifact(artifactDir, "final", finalInterval,
                    null, sizing, logMuTrain, logZCorrected);
            writeFrequencySeverityRunSummaryArtifact(artifactDir, winner, sizing,
                    finalInterval, trainCCD.size(), pilotCCD.size(),
                    monitorCCD.size(), finalCCD.size(),
                    certificateValid);
        } catch (RuntimeException ex) {
            failCertificate("frequency/severity PAC V2.1 final artifact failed: "
                    + sanitizeTsv(ex.getMessage()), startTime);
            return;
        }
        if (!certificateValid) {
            printFinalSummary(startTime);
            return;
        }
        System.out.println("[PACK*-frequency-severity] issued assumption-conditional frequency/severity interval:"
                + " premiseId=" + frequencySeverityPremiseId
                + ", tail=" + finalInterval.tailCount
                + "/" + finalInterval.sampleCount
                + ", pUpper="
                + String.format(Locale.ROOT, "%.9g",
                finalInterval.tailProbabilityUpper)
                + ", S0=" + frequencySeverityCap
                + ", epsilon="
                + String.format(Locale.ROOT, "%.9f", epsilon)
                + ", targetReached="
                + (epsilon <= targetEpsilon + 1.0e-12));
        printFinalSummary(startTime);
    }

    private double computeFrequencySeverityTrainingLogMu(
            List<CCDResult> samples, EtaCorrections eta) {
        double[] logA = new double[samples.size()];
        double[] logV = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + computeEtaEnergy(sample.conf, eta);
            logA[i] = frequencySeveritySourceLogWeight(
                    sample, targetEnergy);
            logV[i] = logA[i]
                    + (targetEnergy - sample.eTrue) / RT;
        }
        return logSumExp(logV) - logSumExp(logA);
    }

    private double[] computeFrequencySeverityLogRelativeWeights(
            List<CCDResult> samples, EtaCorrections eta,
            double logMuTrain) {
        double[] logRelative = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            CCDResult sample = samples.get(i);
            double etaEnergy = computeEtaEnergy(sample.conf, eta);
            logRelative[i] = (sample.eMin + etaEnergy
                    - sample.eTrue) / RT - logMuTrain;
        }
        return logRelative;
    }

    private void setFrequencySeverityDiagnosticStatistics(
            List<CCDResult> samples, EtaCorrections eta,
            double[] logRelative,
            PackStarFrequencySeverityPAC.Interval interval) {
        double sumResidual = 0.0;
        double sumResidual2 = 0.0;
        for (CCDResult sample : samples) {
            double eEta = sample.eMin
                    + computeEtaEnergy(sample.conf, eta);
            double residual = sample.eTrue - eEta;
            sumResidual += residual;
            sumResidual2 += residual * residual;
        }
        int n = samples.size();
        meanResidual = sumResidual / n;
        stdResidual = n > 1
                ? Math.sqrt(Math.max(0.0,
                (sumResidual2 - n * meanResidual * meanResidual)
                        / (n - 1))) : 0.0;
        double logMeanR = logSumExp(logRelative) - Math.log(n);
        double[] twice = new double[n];
        for (int i = 0; i < n; i++) twice[i] = 2.0 * logRelative[i];
        double logM2 = logSumExp(twice) - Math.log(n)
                - 2.0 * logMeanR;
        double m2 = logM2 < Math.log(Double.MAX_VALUE)
                ? Math.max(1.0, Math.exp(logM2))
                : Double.POSITIVE_INFINITY;
        meanPsi = logMeanR < Math.log(Double.MAX_VALUE)
                ? Math.exp(logMeanR) : Double.MAX_VALUE;
        varPsi = interval.bulkVariance;
        cvPsi = Double.isFinite(m2)
                ? Math.sqrt(Math.max(0.0, m2 - 1.0))
                : Double.MAX_VALUE;
    }

    private File prepareFrequencySeverityArtifactDirectory() {
        File root = new File(frequencySeverityOutputDir);
        if ((!root.exists() && !root.mkdirs()) || !root.isDirectory()) {
            throw new IllegalStateException(
                    "cannot create frequency/severity output root " + root);
        }
        String identityHash = Long.toUnsignedString(
                deriveRandomSeed(0L, randomStreamIdentity), 16);
        File stateDir = new File(root, "state-" + identityHash);
        if (stateDir.exists() || !stateDir.mkdir()) {
            throw new IllegalStateException(
                    "frequency/severity route refuses to overwrite state artifact directory "
                            + stateDir);
        }
        return stateDir;
    }

    /** Record why eta training cannot start when the initial q_m DP is empty. */
    private void writeFrequencySeverityInitialDpFailureQuietly(
            double logZUpper, double logZLower, String reason) {
        if (!frequencySeverityProduction
                || frequencySeverityOutputDir == null) return;
        try {
            File dir = prepareFrequencySeverityArtifactDirectory();
            writeFrequencySeverityProtocolArtifact(dir, logZLower);
            long unaryCells = 0L;
            long finiteRigidUnaryCells = 0L;
            long finiteMinimizingUnaryCells = 0L;
            int positionsWithoutFiniteRigidUnary = 0;
            int positionsWithoutFiniteMinimizingUnary = 0;
            int numPos = rcs.getNumPos();
            for (int pos = 0; pos < numPos; pos++) {
                boolean anyFiniteRigid = false;
                boolean anyFiniteMinimizing = false;
                int numRCs = rcs.getNum(pos);
                for (int localRc = 0; localRc < numRCs; localRc++) {
                    int globalRc = rcs.get(pos, localRc);
                    unaryCells++;
                    if (Double.isFinite(branchRigidEmat.getOneBody(
                            pos, globalRc))) {
                        finiteRigidUnaryCells++;
                        anyFiniteRigid = true;
                    }
                    if (Double.isFinite(branchMinimizingEmat.getOneBody(
                            pos, globalRc))) {
                        finiteMinimizingUnaryCells++;
                        anyFiniteMinimizing = true;
                    }
                }
                if (!anyFiniteRigid) positionsWithoutFiniteRigidUnary++;
                if (!anyFiniteMinimizing) {
                    positionsWithoutFiniteMinimizingUnary++;
                }
            }
            int graphEdges = 0;
            int edgesWithoutFiniteRigidPair = 0;
            int edgesWithoutFiniteMinimizingPair = 0;
            long pairCells = 0L;
            long finiteRigidPairCells = 0L;
            long finiteMinimizingPairCells = 0L;
            for (int pos1 = 0; pos1 < numPos; pos1++) {
                int numRCs1 = rcs.getNum(pos1);
                for (int pos2 = pos1 + 1; pos2 < numPos; pos2++) {
                    if (!interactionGraph.hasEdge(pos1, pos2)) continue;
                    graphEdges++;
                    boolean anyFiniteRigid = false;
                    boolean anyFiniteMinimizing = false;
                    int numRCs2 = rcs.getNum(pos2);
                    for (int localRc1 = 0;
                         localRc1 < numRCs1; localRc1++) {
                        int globalRc1 = rcs.get(pos1, localRc1);
                        for (int localRc2 = 0;
                             localRc2 < numRCs2; localRc2++) {
                            int globalRc2 = rcs.get(pos2, localRc2);
                            pairCells++;
                            if (Double.isFinite(branchRigidEmat.getPairwise(
                                    pos1, globalRc1,
                                    pos2, globalRc2))) {
                                finiteRigidPairCells++;
                                anyFiniteRigid = true;
                            }
                            if (Double.isFinite(
                                    branchMinimizingEmat.getPairwise(
                                            pos1, globalRc1,
                                            pos2, globalRc2))) {
                                finiteMinimizingPairCells++;
                                anyFiniteMinimizing = true;
                            }
                        }
                    }
                    if (!anyFiniteRigid) edgesWithoutFiniteRigidPair++;
                    if (!anyFiniteMinimizing) {
                        edgesWithoutFiniteMinimizingPair++;
                    }
                }
            }
            File output = newFrequencySeverityArtifact(
                    dir, "initial_proposal_support_preflight.tsv");
            try (PrintWriter writer = new PrintWriter(
                    output, StandardCharsets.UTF_8.name())) {
                writer.println("key\tvalue");
                writeFrequencySeverityKey(writer, "schema",
                        "packstar-adaptive-frequency-severity-support-preflight-v1");
                writeFrequencySeverityKey(writer, "classification",
                        logZUpper == Double.NEGATIVE_INFINITY
                                ? "initial-dp-zero-mass"
                                : "initial-dp-nonfinite");
                writeFrequencySeverityKey(writer, "reason", reason);
                writeFrequencySeverityKey(writer, "initialLogZUpper", logZUpper);
                writeFrequencySeverityKey(writer, "initialLogZLower", logZLower);
                writeFrequencySeverityKey(writer, "positions", numPos);
                writeFrequencySeverityKey(writer, "unaryCells", unaryCells);
                writeFrequencySeverityKey(writer, "finiteRigidUnaryCells",
                        finiteRigidUnaryCells);
                writeFrequencySeverityKey(writer,
                        "finiteMinimizingUnaryCells",
                        finiteMinimizingUnaryCells);
                writeFrequencySeverityKey(writer,
                        "positionsWithoutFiniteRigidUnary",
                        positionsWithoutFiniteRigidUnary);
                writeFrequencySeverityKey(writer,
                        "positionsWithoutFiniteMinimizingUnary",
                        positionsWithoutFiniteMinimizingUnary);
                writeFrequencySeverityKey(writer, "interactionGraphEdges",
                        graphEdges);
                writeFrequencySeverityKey(writer, "pairCells", pairCells);
                writeFrequencySeverityKey(writer, "finiteRigidPairCells",
                        finiteRigidPairCells);
                writeFrequencySeverityKey(writer,
                        "finiteMinimizingPairCells",
                        finiteMinimizingPairCells);
                writeFrequencySeverityKey(writer,
                        "edgesWithoutFiniteRigidPair",
                        edgesWithoutFiniteRigidPair);
                writeFrequencySeverityKey(writer,
                        "edgesWithoutFiniteMinimizingPair",
                        edgesWithoutFiniteMinimizingPair);
                writeFrequencySeverityKey(writer, "qMSamplingPossible", false);
                writeFrequencySeverityKey(writer, "etaTrainingPossible", false);
                writeFrequencySeverityKey(writer,
                        "ordinaryEtaTrainingCanRepair", false);
                writeFrequencySeverityKey(writer, "requiredNextTier",
                        "DP-support-or-numerical-audit-before-eta");
            }
            writeFrequencySeverityInitialSelectedCellAudit(dir);
            writeFrequencySeverityInitialDpEdgeAudit(dir);
            writeFrequencySeverityFailureArtifactQuietly(
                    dir, reason, 0, 0, 0, 0);
        } catch (IOException | RuntimeException ex) {
            System.out.println("[PACK*-frequency-severity] WARNING: cannot write initial DP support artifact: "
                    + sanitizeTsv(ex.getMessage()));
        }
    }

    /**
     * Emit the exact sequence-filtered unary and pair cells used by q_m.
     * Local RC indices are DP mixed-radix digits; global RC indices address
     * the energy matrices. Keeping both distinguishes a support failure from
     * a local/global RC indexing defect.
     */
    private void writeFrequencySeverityInitialSelectedCellAudit(File dir)
            throws IOException {
        File output = newFrequencySeverityArtifact(
                dir, "initial_proposal_selected_cells.tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tkind\tpos1\tpos1Name\tlocalRc1"
                    + "\tglobalRc1\tconfType1\tpos2\tpos2Name"
                    + "\tlocalRc2\tglobalRc2\tconfType2"
                    + "\tinteractionGraphEdge"
                    + "\trigidEnergy\trigidStatus"
                    + "\tminimizingEnergy\tminimizingStatus");
            int numPos = rcs.getNumPos();
            for (int pos = 0; pos < numPos; pos++) {
                for (int localRc = 0; localRc < rcs.getNum(pos); localRc++) {
                    int globalRc = rcs.get(pos, localRc);
                    double rigid = branchRigidEmat.getOneBody(pos, globalRc);
                    double minimizing = branchMinimizingEmat.getOneBody(
                            pos, globalRc);
                    writer.printf(Locale.ROOT,
                            "packstar-adaptive-frequency-severity-selected-support-v1"
                                    + "\tunary\t%d\t%s\t%d\t%d\t%s"
                                    + "\t-1\tNA\t-1\t-1\tNA\tfalse"
                                    + "\t%.17g\t%s\t%.17g\t%s%n",
                            pos, sanitizeTsv(confSpace.name(pos)),
                            localRc, globalRc,
                            sanitizeTsv(confSpace.confType(pos, globalRc)),
                            rigid, frequencySeverityNumericStatus(rigid),
                            minimizing,
                            frequencySeverityNumericStatus(minimizing));
                }
            }
            for (int pos1 = 0; pos1 < numPos; pos1++) {
                for (int pos2 = pos1 + 1; pos2 < numPos; pos2++) {
                    boolean graphEdge = interactionGraph.hasEdge(pos1, pos2);
                    for (int localRc1 = 0;
                         localRc1 < rcs.getNum(pos1); localRc1++) {
                        int globalRc1 = rcs.get(pos1, localRc1);
                        for (int localRc2 = 0;
                             localRc2 < rcs.getNum(pos2); localRc2++) {
                            int globalRc2 = rcs.get(pos2, localRc2);
                            double rigid = branchRigidEmat.getPairwise(
                                    pos1, globalRc1, pos2, globalRc2);
                            double minimizing =
                                    branchMinimizingEmat.getPairwise(
                                            pos1, globalRc1,
                                            pos2, globalRc2);
                            writer.printf(Locale.ROOT,
                                    "packstar-adaptive-frequency-severity-selected-support-v1"
                                            + "\tpair\t%d\t%s\t%d\t%d\t%s"
                                            + "\t%d\t%s\t%d\t%d\t%s\t%s"
                                            + "\t%.17g\t%s\t%.17g\t%s%n",
                                    pos1, sanitizeTsv(confSpace.name(pos1)),
                                    localRc1, globalRc1,
                                    sanitizeTsv(confSpace.confType(
                                            pos1, globalRc1)),
                                    pos2, sanitizeTsv(confSpace.name(pos2)),
                                    localRc2, globalRc2,
                                    sanitizeTsv(confSpace.confType(
                                            pos2, globalRc2)),
                                    Boolean.toString(graphEdge),
                                    rigid,
                                    frequencySeverityNumericStatus(rigid),
                                    minimizing,
                                    frequencySeverityNumericStatus(minimizing));
                        }
                    }
                }
            }
            if (writer.checkError()) {
                throw new IOException("write failed for " + output);
            }
        }
    }

    private static String frequencySeverityNumericStatus(double value) {
        if (Double.isFinite(value)) return "finite";
        if (Double.isNaN(value)) return "nan";
        return value < 0.0 ? "negative-infinity" : "positive-infinity";
    }

    private void writeFrequencySeverityInitialDpEdgeAudit(File dir) {
        File output = newFrequencySeverityArtifact(
                dir, "initial_dp_edge_table_audit.tsv");
        List<RootedTreeEdge> edges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, edges);
        long remainingAuditEntries = 10_000_000L;
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tedgeIndex\trootEdge\tmPositions"
                    + "\tlambdaPositions\tmStateCount\tscannedEntries"
                    + "\tcompleteScan\tlowerFinite\tlowerNegativeInfinity"
                    + "\tlowerPositiveInfinity\tlowerNaN\tlowerMinimumFinite"
                    + "\tlowerMaximumFinite\tupperFinite"
                    + "\tupperNegativeInfinity\tupperPositiveInfinity"
                    + "\tupperNaN\tupperMinimumFinite\tupperMaximumFinite");
            for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++) {
                RootedTreeEdge edge = edges.get(edgeIndex);
                long stateCount = edge.getMStateCount();
                long scanCount = Math.min(stateCount,
                        Math.max(0L, remainingAuditEntries));
                remainingAuditEntries -= scanCount;
                long[] finite = new long[2];
                long[] negativeInfinity = new long[2];
                long[] positiveInfinity = new long[2];
                long[] nan = new long[2];
                double[] minimum = {
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY};
                double[] maximum = {
                        Double.NEGATIVE_INFINITY,
                        Double.NEGATIVE_INFINITY};
                for (long state = 0; state < scanCount; state++) {
                    double[] values = {
                            edge.getLogZLower(state),
                            edge.getLogZUpper(state)};
                    for (int bound = 0; bound < values.length; bound++) {
                        double value = values[bound];
                        if (Double.isFinite(value)) {
                            finite[bound]++;
                            minimum[bound] = Math.min(minimum[bound], value);
                            maximum[bound] = Math.max(maximum[bound], value);
                        } else if (Double.isNaN(value)) {
                            nan[bound]++;
                        } else if (value < 0.0) {
                            negativeInfinity[bound]++;
                        } else {
                            positiveInfinity[bound]++;
                        }
                    }
                }
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-dp-edge-audit-v1"
                                + "\t%d\t%s\t%s\t%s\t%d\t%d\t%s"
                                + "\t%d\t%d\t%d\t%d\t%.17g\t%.17g"
                                + "\t%d\t%d\t%d\t%d\t%.17g\t%.17g%n",
                        edgeIndex,
                        Boolean.toString(edge == rootedRootEdge),
                        formatFrequencySeverityConf(
                                edge.getMPositionsSorted()),
                        formatFrequencySeverityConf(
                                edge.getLambdaPositionsSorted()),
                        stateCount, scanCount,
                        Boolean.toString(scanCount == stateCount),
                        finite[0], negativeInfinity[0],
                        positiveInfinity[0], nan[0],
                        finite[0] == 0 ? Double.NaN : minimum[0],
                        finite[0] == 0 ? Double.NaN : maximum[0],
                        finite[1], negativeInfinity[1],
                        positiveInfinity[1], nan[1],
                        finite[1] == 0 ? Double.NaN : minimum[1],
                        finite[1] == 0 ? Double.NaN : maximum[1]);
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity initial DP edge audit", ex);
        }
    }

    private void writeFrequencySeverityProtocolArtifact(File dir, double logZRigid) {
        File output = newFrequencySeverityArtifact(dir, "frequency_severity_protocol.tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("key\tvalue");
            writeFrequencySeverityKey(writer, "schema",
                    "packstar-adaptive-frequency-severity-protocol-v7");
            writeFrequencySeverityKey(writer, "algorithm",
                    "AdaptiveFrequencySeverityPAC-crossfit-selected-signed-triple-residual-source-aware-online-refit");
            writeFrequencySeverityKey(writer, "proposalModel",
                    "fixed-alpha-one-single-pair-eta-plus-zero-or-one-model-containing-up-to-K-greedily-selected-sparse-signed-triple-factors");
            writeFrequencySeverityKey(writer, "tripleEtaEnabled",
                    frequencySeverityTripleEtaEnabled);
            writeFrequencySeverityKey(writer, "tripleEtaMaximumScale",
                    frequencySeverityTripleEtaScale);
            writeFrequencySeverityKey(writer, "tripleEtaScaleGrid",
                    Arrays.toString(frequencySeverityTripleEtaScaleGrid));
            writeFrequencySeverityKey(writer, "tripleEtaScaleSelection",
                    "fixed-one-when-included-not-continuously-tuned");
            writeFrequencySeverityKey(writer, "tripleEtaModelOrderSelection",
                    "pair-only-or-selected-triple-by-crossfit");
            writeFrequencySeverityKey(writer,
                    "tripleEtaAssignmentCapFallback",
                    "pair-only-when-clique-workload-exceeds-cap");
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumPartialAssignments",
                    frequencySeverityTripleEtaMaxAssignments);
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumSelectedPositionTriples",
                    frequencySeverityTripleEtaMaxPositionTriples);
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumFillEdges",
                    frequencySeverityTripleEtaMaxFillEdges);
            writeFrequencySeverityKey(writer,
                    "tripleEtaMinimumCellContexts",
                    frequencySeverityTripleEtaMinCellContexts);
            writeFrequencySeverityKey(writer,
                    "tripleEtaPriorStrength",
                    frequencySeverityTripleEtaPriorStrength);
            writeFrequencySeverityKey(writer,
                    "tripleEtaLocalPriorCapKcal",
                    frequencySeverityTripleEtaLocalCapKcal);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualCapKcal",
                    frequencySeverityTripleEtaResidualCapKcal);
            writeFrequencySeverityKey(writer, "tripleEtaFactorDomain",
                    "all-position-triples-with-a-bounded-zero-energy-fill-edge-budget");
            writeFrequencySeverityKey(writer,
                    "tripleEtaLocalJointCorrection",
                    "signed-shrinkage-prior-from-shared-three-pair-minimum-minus-independent-three-pair-minima");
            writeFrequencySeverityKey(writer,
                    "tripleEtaProposalFactor",
                    "crossfit-signed-Etrue-minus-Em-minus-single-pair-eta-residual-table");
            writeFrequencySeverityKey(writer,
                    "tripleEtaPairAllocation",
                    "none-residual-is-defined-after-single-pair-eta");
            writeFrequencySeverityKey(writer, "tripleEtaNormalizer",
                    "exact-higher-order-branch-DP-with-zero-energy-fill-edge-redecomposition");
            writeFrequencySeverityKey(writer,
                    "tripleEtaAdditionalFullConformationCcdCalls", 0);
            writeFrequencySeverityKey(writer,
                    "tripleEtaFitUsesExistingSampledFullCcd", true);
            writeFrequencySeverityKey(writer, "candidateObjective",
                    "select-pair-only-or-selected-triple-and-shrink-by-min-pooled-out-of-fold-predicted-final-N-then-epsilon-at-cap-then-target-M2");
            writeFrequencySeverityKey(writer, "crossfitAggregation",
                    "history-plus-complementary-half-of-newest-on-policy-batch-fits-with-source-aware-held-out-scoring");
            writeFrequencySeverityKey(writer, "crossfitHardGates",
                    "all-fold-finite-and-shift-ESS-plus-pooled-bulk-and-reachability");
            writeFrequencySeverityKey(writer, "crossfitCoverageRole",
                    "diagnostic-only");
            writeFrequencySeverityKey(writer, "fullRefitCoverageRole",
                    "hard-gate-before-winner-DP");
            writeFrequencySeverityKey(writer, "guarantee",
                    "assumption-conditional-frequency-severity-PAC");
            writeFrequencySeverityKey(writer, "severityPremiseId", frequencySeverityPremiseId);
            writeFrequencySeverityKey(writer, "severityPremiseExternallyValidated", "false");
            writeFrequencySeverityKey(writer, "severityCapS0", frequencySeverityCap);
            writeFrequencySeverityKey(writer, "relativeBoundKcal", frequencySeverityRelativeBoundKcal);
            writeFrequencySeverityKey(writer, "relativeGauge",
                    "self-normalized-source-aware-importance-estimate-of-E_qeta[exp((Eeta-Etrue)/RT)]");
            writeFrequencySeverityKey(writer, "confidenceDelta", delta);
            writeFrequencySeverityKey(writer, "confidenceFamily",
                    "two-events-times-three-Kstar-partition-functions");
            writeFrequencySeverityKey(writer, "bulkDelta",
                    frequencySeverityEventDelta());
            writeFrequencySeverityKey(writer, "frequencyDelta",
                    frequencySeverityEventDelta());
            writeFrequencySeverityKey(writer, "targetEpsilon", targetEpsilon);
            writeFrequencySeverityKey(writer, "trainSamples", trainSamples);
            writeFrequencySeverityKey(writer, "discoverySamplesPerRound",
                    frequencySeverityDiscoverySamples);
            writeFrequencySeverityKey(writer, "maximumDiscoverySamplesPerRound",
                    frequencySeverityDiscoveryMaxSamples);
            writeFrequencySeverityKey(writer, "discoveryExtensionPolicy",
                    "same-frozen-proposal-when-bulk-lower-is-zero-and-point-sizing-is-reachable-at-target");
            writeFrequencySeverityKey(writer, "maximumRefits",
                    frequencySeverityMaxRefits);
            writeFrequencySeverityKey(writer, "validationSamples",
                    frequencySeverityValidationSamples);
            writeFrequencySeverityKey(writer, "monitorSamples", monitorSamples);
            writeFrequencySeverityKey(writer, "maxFinalSamples", maxEstSamples);
            writeFrequencySeverityKey(writer, "unreachableFinalSamples",
                    configuredFrequencySeverityUnreachableSamples());
            writeFrequencySeverityKey(writer, "folds", frequencySeverityFolds);
            writeFrequencySeverityKey(writer, "shrinkGrid", formatFrequencySeverityShrinkGrid());
            writeFrequencySeverityKey(writer, "alphaGrid", Arrays.toString(frequencySeverityAlphaGrid));
            writeFrequencySeverityKey(writer, "alphaSelection",
                    "fixed-one-not-tuned");
            writeFrequencySeverityKey(writer, "minimumShiftEssFraction",
                    frequencySeverityMinShiftEssFraction);
            writeFrequencySeverityKey(writer,
                    "discoveryMinimumShiftEssFraction",
                    frequencySeverityDiscoveryMinShiftEssFraction);
            writeFrequencySeverityKey(writer, "minimumDistinctContextCount",
                    frequencySeverityMinTrainCount);
            writeFrequencySeverityKey(writer, "maximumUndertrainedMassAmplification",
                    frequencySeverityMaxUndertrainedAmplification);
            writeFrequencySeverityKey(writer, "sizingSafety", frequencySeveritySizeSafety);
            writeFrequencySeverityKey(writer, "severityTestAlpha", frequencySeverityTestAlpha);
            writeFrequencySeverityKey(writer, "severityTest",
                    "equal-mixture-of-sample-mean-and-fixed-betting-products");
            writeFrequencySeverityKey(writer, "severityBettingLambdas",
                    Arrays.toString(
                            PackStarFrequencySeverityPAC.severityBettingLambdas()));
            writeFrequencySeverityKey(writer, "candidateSpecificCcd", 0);
            writeFrequencySeverityKey(writer, "candidateSpecificDp", 0);
            writeFrequencySeverityKey(writer, "proposalCorrectedDpMaximum",
                    frequencySeverityMaxRefits + 1);
            writeFrequencySeverityKey(writer, "etaZeroReusesInitialDpOnlyAtRoundZero",
                    false);
            writeFrequencySeverityKey(writer, "legacyIterativeRepairBypassed", "true");
            writeFrequencySeverityKey(writer, "activeOnPolicyRefitEnabled", "true");
            writeFrequencySeverityKey(writer, "failedDiscoveryBatchReusedForEtaFit", "true");
            writeFrequencySeverityKey(writer, "successfulDiscoveryBatchReusedForEtaFit", "false");
            writeFrequencySeverityKey(writer, "validationBatchReusedForEtaFit", "false");
            writeFrequencySeverityKey(writer, "defensiveMixtureTierEnabled", "false");
            writeFrequencySeverityKey(writer, "initialLogZRigid", logZRigid);
            writeFrequencySeverityKey(writer, "randomStreamIdentityHash",
                    Long.toUnsignedString(
                            deriveRandomSeed(0L, randomStreamIdentity), 16));
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity protocol artifact", ex);
        }
    }

    private void writeFrequencySeverityTrainingArtifacts(
            File dir, List<CCDResult> samples,
            FrequencySeverityEtaTraining training) {
        writeFrequencySeverityTrainingArtifacts(
                dir, "", samples, training);
    }

    private void writeFrequencySeverityTrainingArtifacts(
            File dir, String artifactPrefix,
            List<CCDResult> samples,
            FrequencySeverityEtaTraining training) {
        File sampleOutput = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "frequency_severity_training_samples.tsv"));
        try (PrintWriter writer = new PrintWriter(
                sampleOutput, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tindex\tnewestBatchFold\tconf\teTrueKcal"
                    + "\teMinKcal\tlogBaselineResidualWeight"
                    + "\tsourceProposalId\tsourceProposalEnergyKcal"
                    + "\tsourceProposalLogZ");
            for (int index = 0; index < samples.size(); index++) {
                CCDResult sample = samples.get(index);
                int newestBatchFold = index < training.historySampleCount
                        ? -1 : (index - training.historySampleCount)
                        % frequencySeverityFolds;
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-training-sample-v1"
                                + "\t%d\t%d\t%s\t%.17g\t%.17g\t%.17g"
                                + "\t%s\t%.17g\t%.17g%n",
                        index, newestBatchFold,
                        formatFrequencySeverityConf(sample.conf),
                        sample.eTrue, sample.eMin,
                        (sample.eMin - sample.eTrue) / RT,
                        sanitizeTsv(sample.sourceProposalId == null
                                ? "NA" : sample.sourceProposalId),
                        sample.sourceProposalEnergy,
                        sample.sourceProposalLogZ);
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity training samples", ex);
        }

        File observationOutput = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_training_cell_observations.tsv"));
        try (PrintWriter writer = new PrintWriter(
                observationOutput, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tsampleIndex\tfold\tcellType"
                    + "\tpos1\trc1\tpos2\trc2\tcorrectionKcal");
            for (FrequencySeverityEtaCellObservation observation
                    : training.observations) {
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-eta-cell-observation-v1"
                                + "\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%.17g%n",
                        observation.sampleIndex, observation.fold,
                        observation.pair ? "pair" : "unary",
                        observation.pos1, observation.rc1,
                        observation.pos2, observation.rc2,
                        observation.correctionKcal);
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity eta cell observations", ex);
        }

        File statisticsOutput = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_fit_cell_statistics.tsv"));
        try (PrintWriter writer = new PrintWriter(
                statisticsOutput, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tfitScope\tfitSampleCount\tcellType"
                    + "\tpos1\trc1\tpos2\trc2\trawIidCount"
                    + "\tdistinctContextCount\trawMeanCorrectionKcal");
            writeFrequencySeverityEtaStatisticsScope(
                    writer, "all", training.allFitSampleCount,
                    training.all, training.allCoverage);
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                writeFrequencySeverityEtaStatisticsScope(
                        writer, "fold-" + fold,
                        training.foldFitSampleCounts[fold],
                        training.folds[fold],
                        training.foldCoverage[fold]);
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity eta fit statistics", ex);
        }
    }

    private void writeFrequencySeverityEtaStatisticsScope(
            PrintWriter writer, String scope, int fitSampleCount,
            EtaCorrections raw, FrequencySeverityEtaCoverage coverage) {
        for (int pos = 0; pos < raw.oneBody.length; pos++) {
            for (int rc = 0; rc < raw.oneBody[pos].length; rc++) {
                int count = raw.oneBodyCounts[pos][rc];
                if (count <= 0) continue;
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-eta-fit-cell-v1"
                                + "\t%s\t%d\tunary\t%d\t%d\t-1\t-1"
                                + "\t%d\t%d\t%.17g%n",
                        scope, fitSampleCount, pos, rc, count,
                        coverage.unary(pos, rc),
                        raw.getOneBodyEta(pos, rc));
            }
        }
        List<Long> pairKeys = new ArrayList<>(raw.pairSums.keySet());
        Collections.sort(pairKeys);
        for (long key : pairKeys) {
            double[] stored = raw.pairSums.get(key);
            int count = (int) Math.round(stored[1]);
            if (count <= 0) continue;
            int pos1 = (int) ((key >>> 48) & 0xffffL);
            int rc1 = (int) ((key >>> 32) & 0xffffL);
            int pos2 = (int) ((key >>> 16) & 0xffffL);
            int rc2 = (int) (key & 0xffffL);
            writer.printf(Locale.ROOT,
                    "packstar-adaptive-frequency-severity-eta-fit-cell-v1"
                            + "\t%s\t%d\tpair\t%d\t%d\t%d\t%d"
                            + "\t%d\t%d\t%.17g%n",
                    scope, fitSampleCount,
                    pos1, rc1, pos2, rc2, count,
                    coverage.pair(pos1, rc1, pos2, rc2),
                    raw.getPairEta(pos1, rc1, pos2, rc2));
        }
    }

    private void writeFrequencySeverityCandidateArtifact(
            File dir, List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        writeFrequencySeverityCandidateArtifact(
                dir, "", candidates, selectedId);
    }

    private void writeFrequencySeverityCandidateArtifact(
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(dir,
                frequencySeverityPrefixedArtifactName(
                        artifactPrefix, "eta_candidate_scores.tsv"));
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tcandidate\tselected\teligible\tvalidationFold"
                    + "\tunaryShrink\tpairShrink\talpha\ttripleEtaScale\tn"
                    + "\tlogMuOutOfFold\tshiftEss\tshiftEssFraction"
                    + "\tshiftPass\ttargetM2Diagnostic"
                    + "\ttargetEssFractionDiagnostic\tbulkMean"
                    + "\tbulkVariance\ttailProbability\tpredictedFinalN"
                    + "\tpredictedEpsilonAtFinalN\tpredictedEpsilonAtMax"
                    + "\treachableAtMax\tbulkPass"
                    + "\tbaselineUnaryUndertrainedMass"
                    + "\tcandidateUnaryUndertrainedMass"
                    + "\tunaryUndertrainedAmplification"
                    + "\tbaselinePairUndertrainedMass"
                    + "\tcandidatePairUndertrainedMass"
                    + "\tpairUndertrainedAmplification\tcoveragePass"
                    + "\tfinite");
            for (FrequencySeverityCandidateScore candidate : candidates) {
                for (FrequencySeverityFoldScore fold : candidate.folds) {
                    writer.printf(Locale.ROOT,
                            "packstar-adaptive-frequency-severity-candidate-score-v2\t%s\t%s\t%s"
                                    + "\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%d"
                                    + "\t%.17g\t%.17g\t%.17g\t%s"
                                    + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                    + "\t%.17g\t%d\t%.17g\t%.17g\t%s\t%s"
                                    + "\t%.17g\t%.17g\t%.17g"
                                    + "\t%.17g\t%.17g\t%.17g\t%s\t%s%n",
                            candidate.id,
                            Boolean.toString(candidate.id.equals(selectedId)),
                            Boolean.toString(candidate.eligible),
                            fold.validationFold,
                            candidate.shrink.unary,
                            candidate.shrink.pair,
                            candidate.alpha,
                            candidate.tripleEtaScale,
                            fold.sampleCount,
                            fold.logMu,
                            fold.shiftEss,
                            fold.shiftEssFraction,
                            Boolean.toString(fold.shiftPass),
                            fold.targetM2,
                            fold.targetEssFraction,
                            fold.bulkMean,
                            fold.bulkVariance,
                            fold.tailProbability,
                            fold.sizing.finalSamples,
                            fold.sizing.epsilonAtFinalSamples,
                            fold.sizing.epsilonAtMaxSamples,
                            Boolean.toString(fold.sizing.reachableAtMax),
                            Boolean.toString(fold.bulkPass),
                            fold.baselineUnaryUndertrainedMass,
                            fold.candidateUnaryUndertrainedMass,
                            fold.unaryUndertrainedAmplification,
                            fold.baselinePairUndertrainedMass,
                            fold.candidatePairUndertrainedMass,
                            fold.pairUndertrainedAmplification,
                            Boolean.toString(fold.coveragePass),
                            Boolean.toString(fold.finite));
                }
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity candidate scores", ex);
        }
    }

    private void writeFrequencySeverityCrossfitArtifact(
            File dir, List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        writeFrequencySeverityCrossfitArtifact(
                dir, "", candidates, selectedId);
    }

    private void writeFrequencySeverityCrossfitArtifact(
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_candidate_crossfit_pooled.tsv"));
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tcandidate\tselected\teligible"
                    + "\tunaryShrink\tpairShrink\talpha\ttripleEtaScale\tn"
                    + "\tlogMuPooledDiagnostic\tshiftEss"
                    + "\tshiftEssFraction\tshiftPass"
                    + "\tallFoldsFinite\tallFoldsShiftPass"
                    + "\tallFoldsCoveragePassDiagnostic"
                    + "\ttargetM2Diagnostic"
                    + "\ttargetEssFractionDiagnostic"
                    + "\tbulkMean\tbulkVariance\ttailProbability"
                    + "\tpredictedFinalN\tpredictedEpsilonAtFinalN"
                    + "\tpredictedEpsilonAtMax\treachableAtMax"
                    + "\tbulkPass"
                    + "\tbaselineUnaryUndertrainedMass"
                    + "\tcandidateUnaryUndertrainedMass"
                    + "\tunaryUndertrainedAmplification"
                    + "\tbaselinePairUndertrainedMass"
                    + "\tcandidatePairUndertrainedMass"
                    + "\tpairUndertrainedAmplification"
                    + "\tcoveragePassDiagnostic\tfinite");
            for (FrequencySeverityCandidateScore candidate : candidates) {
                FrequencySeverityCrossfitScore score = candidate.crossfit;
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-crossfit-pooled-v2"
                                + "\t%s\t%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%d"
                                + "\t%.17g\t%.17g\t%.17g\t%s\t%s\t%s\t%s"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%d\t%.17g\t%.17g\t%s\t%s"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%s\t%s%n",
                        candidate.id,
                        Boolean.toString(candidate.id.equals(selectedId)),
                        Boolean.toString(candidate.eligible),
                        candidate.shrink.unary,
                        candidate.shrink.pair,
                        candidate.alpha,
                        candidate.tripleEtaScale,
                        score.sampleCount,
                        score.logMu,
                        score.shiftEss,
                        score.shiftEssFraction,
                        Boolean.toString(score.shiftPass),
                        Boolean.toString(score.allFoldsFinite),
                        Boolean.toString(score.allFoldsShiftPass),
                        Boolean.toString(score.allFoldsCoveragePass),
                        score.targetM2,
                        score.targetEssFraction,
                        score.bulkMean,
                        score.bulkVariance,
                        score.tailProbability,
                        score.sizing.finalSamples,
                        score.sizing.epsilonAtFinalSamples,
                        score.sizing.epsilonAtMaxSamples,
                        Boolean.toString(score.sizing.reachableAtMax),
                        Boolean.toString(score.bulkPass),
                        score.baselineUnaryUndertrainedMass,
                        score.candidateUnaryUndertrainedMass,
                        score.unaryUndertrainedAmplification,
                        score.baselinePairUndertrainedMass,
                        score.candidatePairUndertrainedMass,
                        score.pairUndertrainedAmplification,
                        Boolean.toString(score.coveragePass),
                        Boolean.toString(score.finite));
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity pooled cross-fit scores", ex);
        }
    }

    private void writeFrequencySeverityRefitArtifact(
            File dir, List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        writeFrequencySeverityRefitArtifact(
                dir, "", candidates, selectedId);
    }

    private void writeFrequencySeverityRefitArtifact(
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_candidate_refit_audit.tsv"));
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tcandidate\ttripleEtaScale\tselected\tcrossfitEligible"
                    + "\tfullRefitEligible\tshiftEssFraction\tshiftPass"
                    + "\tbaselineUnaryUndertrainedMass"
                    + "\tcandidateUnaryUndertrainedMass"
                    + "\tunaryUndertrainedAmplification"
                    + "\tbaselinePairUndertrainedMass"
                    + "\tcandidatePairUndertrainedMass"
                    + "\tpairUndertrainedAmplification"
                    + "\tcoveragePass\tfinite");
            for (FrequencySeverityCandidateScore candidate : candidates) {
                FrequencySeverityRefitAudit audit = candidate.refitAudit;
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-refit-audit-v2"
                                + "\t%s\t%.17g\t%s\t%s\t%s\t%.17g\t%s"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%s\t%s%n",
                        candidate.id,
                        candidate.tripleEtaScale,
                        Boolean.toString(candidate.id.equals(selectedId)),
                        Boolean.toString(candidate.crossfitEligible),
                        Boolean.toString(audit.eligible()),
                        audit.shiftEssFraction,
                        Boolean.toString(audit.shiftPass),
                        audit.baselineUnaryUndertrainedMass,
                        audit.candidateUnaryUndertrainedMass,
                        audit.unaryUndertrainedAmplification,
                        audit.baselinePairUndertrainedMass,
                        audit.candidatePairUndertrainedMass,
                        audit.pairUndertrainedAmplification,
                        Boolean.toString(audit.coveragePass),
                        Boolean.toString(audit.finite));
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity full-refit audit", ex);
        }
    }

    private void writeFrequencySeverityWinnerArtifact(
            File dir, FrequencySeverityCandidateScore winner,
            FrequencySeverityRefitAudit refitAudit,
            double logMuTrain, int winnerDpSweeps) {
        PackStarTripleEtaCorrections.ResidualSummary tripleSummary =
                summarizeFrequencySeverityTripleResidual(winner.fullEta);
        File output = newFrequencySeverityArtifact(
                dir, "eta_selected.tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("key\tvalue");
            writeFrequencySeverityKey(writer, "schema",
                    "packstar-adaptive-frequency-severity-selected-eta-v3");
            writeFrequencySeverityKey(writer, "candidate", winner.id);
            writeFrequencySeverityKey(writer, "unaryShrink",
                    winner.shrink.unary);
            writeFrequencySeverityKey(writer, "pairShrink",
                    winner.shrink.pair);
            writeFrequencySeverityKey(writer, "alpha", winner.alpha);
            writeFrequencySeverityKey(writer, "tripleEtaScale",
                    winner.fullEta.tripleEtaScale);
            writeFrequencySeverityKey(writer, "tripleEtaSignatureSha256",
                    winner.fullEta.tripleEta == null ? "NA"
                            : winner.fullEta.tripleEta.signatureSha256);
            writeFrequencySeverityKey(writer, "tripleEtaFactorDefinition",
                    "crossfit-signed-Etrue-minus-Em-minus-single-pair-eta-residual-table");
            writeFrequencySeverityKey(writer,
                    "tripleEtaSelectedPositionTriples",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.positionTripleCount : 0L);
            writeFrequencySeverityKey(writer,
                    "tripleEtaSelectedPositionScopes",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.positionTriples()
                            : "[]");
            writeFrequencySeverityKey(writer,
                    "tripleEtaSelectedFillEdges",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.requiredFillEdges(
                            interactionGraph).size() : 0);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualPositiveFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.positiveFactors);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualNegativeFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.negativeFactors);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualZeroFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.zeroFactors);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualMinimumKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.minimumCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualMaximumKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.maximumCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualMaximumAbsoluteKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.maximumAbsoluteCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "pooledOutOfFoldPredictedFinalSamples",
                    winner.crossfit.sizing.finalSamples);
            writeFrequencySeverityKey(writer,
                    "pooledOutOfFoldEpsilonAtMaximumSamples",
                    winner.crossfit.sizing.epsilonAtMaxSamples);
            writeFrequencySeverityKey(writer,
                    "pooledOutOfFoldShiftEssFraction",
                    winner.crossfit.shiftEssFraction);
            writeFrequencySeverityKey(writer,
                    "allOutOfFoldShiftGatesPassed",
                    winner.crossfit.allFoldsShiftPass);
            writeFrequencySeverityKey(writer,
                    "outOfFoldCoverageDiagnosticPassed",
                    winner.crossfit.coveragePass
                            && winner.crossfit.allFoldsCoveragePass);
            writeFrequencySeverityKey(writer,
                    "fullTrainingShiftEssFraction",
                    refitAudit.shiftEssFraction);
            writeFrequencySeverityKey(writer,
                    "minimumShiftEssFraction",
                    frequencySeverityMinShiftEssFraction);
            writeFrequencySeverityKey(writer,
                    "fullTrainingBaselineUnaryUndertrainedMass",
                    refitAudit.baselineUnaryUndertrainedMass);
            writeFrequencySeverityKey(writer,
                    "fullTrainingCandidateUnaryUndertrainedMass",
                    refitAudit.candidateUnaryUndertrainedMass);
            writeFrequencySeverityKey(writer,
                    "fullTrainingUnaryUndertrainedAmplification",
                    refitAudit.unaryUndertrainedAmplification);
            writeFrequencySeverityKey(writer,
                    "fullTrainingBaselinePairUndertrainedMass",
                    refitAudit.baselinePairUndertrainedMass);
            writeFrequencySeverityKey(writer,
                    "fullTrainingCandidatePairUndertrainedMass",
                    refitAudit.candidatePairUndertrainedMass);
            writeFrequencySeverityKey(writer,
                    "fullTrainingPairUndertrainedAmplification",
                    refitAudit.pairUndertrainedAmplification);
            writeFrequencySeverityKey(writer,
                    "maximumUndertrainedMassAmplification",
                    frequencySeverityMaxUndertrainedAmplification);
            writeFrequencySeverityKey(writer, "logMuTrain", logMuTrain);
            writeFrequencySeverityKey(writer, "proposalDpSweeps",
                    winnerDpSweeps);
            writeFrequencySeverityKey(writer, "candidateDpSweeps", 0);
            writeFrequencySeverityKey(writer, "candidateCcdCalls", 0);
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity selected eta", ex);
        }
    }

    private void writeAdaptiveFrequencySeverityRoundArtifact(
            File dir, String artifactPrefix, int round,
            int historyCount, int scoringBatchCount, int fitPoolCount,
            FrequencySeverityCandidateScore selected,
            boolean discoveryProbe,
            double logMuTrain, double logZProposal,
            int cumulativeProposalDpSweeps,
            List<CCDResult> discoverySamples,
            int initialDiscoveryCount,
            int extensionDiscoveryCount,
            int maximumDiscoveryCount,
            FrequencySeverityEtaCoverage fitCoverage,
            FrequencySeverityPilotEvaluation evaluation,
            boolean passed) {
        PackStarTripleEtaCorrections.ResidualSummary tripleSummary =
                summarizeFrequencySeverityTripleResidual(
                        selected.fullEta);
        File output = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix, "selection.tsv"));
        int unseenUnary = 0;
        int unseenPair = 0;
        int undertrainedUnary = 0;
        int undertrainedPair = 0;
        int unseenUnaryTail = 0;
        int unseenPairTail = 0;
        double logClip = frequencySeverityRelativeBoundKcal / RT;
        for (int i = 0; i < discoverySamples.size(); i++) {
            int[] conf = discoverySamples.get(i).conf;
            boolean sampleUnseenUnary = touchesFrequencySeverityUnaryBelowCount(
                    conf, fitCoverage, 1);
            boolean sampleUnseenPair = touchesFrequencySeverityPairBelowCount(
                    conf, fitCoverage, 1);
            if (sampleUnseenUnary) unseenUnary++;
            if (sampleUnseenPair) unseenPair++;
            if (touchesFrequencySeverityUndertrainedUnary(conf, fitCoverage)) {
                undertrainedUnary++;
            }
            if (touchesFrequencySeverityUndertrainedPair(conf, fitCoverage)) {
                undertrainedPair++;
            }
            boolean tail = evaluation.logRelativeWeights[i] > logClip;
            if (tail && sampleUnseenUnary) unseenUnaryTail++;
            if (tail && sampleUnseenPair) unseenPairTail++;
        }
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("key\tvalue");
            writeFrequencySeverityKey(writer, "schema",
                    "packstar-adaptive-frequency-severity-round-v3");
            writeFrequencySeverityKey(writer, "round", round);
            writeFrequencySeverityKey(writer, "historyCcd", historyCount);
            writeFrequencySeverityKey(writer, "scoringBatchCcd", scoringBatchCount);
            writeFrequencySeverityKey(writer, "fitPoolCcd", fitPoolCount);
            writeFrequencySeverityKey(writer, "selectedCandidate", selected.id);
            writeFrequencySeverityKey(writer, "selectionRole",
                    discoveryProbe ? "discovery-probe" : "final-eligible-candidate");
            writeFrequencySeverityKey(writer, "selectedInitiallyEligible",
                    selected.eligible);
            writeFrequencySeverityKey(writer, "selectedAlpha", selected.alpha);
            writeFrequencySeverityKey(writer, "selectedTripleEtaScale",
                    selected.tripleEtaScale);
            writeFrequencySeverityKey(writer, "selectedUnaryShrink",
                    selected.shrink.unary);
            writeFrequencySeverityKey(writer, "selectedPairShrink",
                    selected.shrink.pair);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaResidualPositiveFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.positiveFactors);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaResidualNegativeFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.negativeFactors);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaResidualMinimumKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.minimumCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaResidualMaximumKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.maximumCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaResidualMaximumAbsoluteKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.maximumAbsoluteCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaPositionTriples",
                    hasTripleEta(selected.fullEta)
                            ? selected.fullEta.tripleEta.positionTripleCount
                            : 0L);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaPositionScopes",
                    hasTripleEta(selected.fullEta)
                            ? selected.fullEta.tripleEta.positionTriples()
                            : "[]");
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaFillEdges",
                    hasTripleEta(selected.fullEta)
                            ? selected.fullEta.tripleEta.requiredFillEdges(
                            interactionGraph).size() : 0);
            writeFrequencySeverityKey(writer, "logMuTrain", logMuTrain);
            writeFrequencySeverityKey(writer, "logZProposal", logZProposal);
            writeFrequencySeverityKey(writer, "cumulativeProposalDpSweeps",
                    cumulativeProposalDpSweeps);
            writeFrequencySeverityKey(writer, "discoveryCcd",
                    discoverySamples.size());
            writeFrequencySeverityKey(writer, "initialDiscoveryCcd",
                    initialDiscoveryCount);
            writeFrequencySeverityKey(writer, "extensionDiscoveryCcd",
                    extensionDiscoveryCount);
            writeFrequencySeverityKey(writer, "maximumDiscoveryCcd",
                    maximumDiscoveryCount);
            writeFrequencySeverityKey(writer,
                    "sameFrozenProposalExtensionTriggered",
                    extensionDiscoveryCount > 0);
            writeFrequencySeverityKey(writer,
                    "sameFrozenProposalExtensionPolicy",
                    "bulk-lower-zero-and-point-sizing-reachable-at-target");
            writeFrequencySeverityKey(writer, "discoveryTailCount",
                    evaluation.interval.tailCount);
            writeFrequencySeverityKey(writer, "discoveryBulkLower",
                    evaluation.interval.bulkLower);
            writeFrequencySeverityKey(writer, "discoveryReachableAtMax",
                    evaluation.sizing.reachableAtMax);
            writeFrequencySeverityKey(writer, "discoveryEpsilonAtMax",
                    evaluation.sizing.epsilonAtMaxSamples);
            writeFrequencySeverityKey(writer, "samplesTouchingUnseenUnary",
                    unseenUnary);
            writeFrequencySeverityKey(writer, "samplesTouchingUnseenPair",
                    unseenPair);
            writeFrequencySeverityKey(writer, "samplesTouchingUndertrainedUnary",
                    undertrainedUnary);
            writeFrequencySeverityKey(writer, "samplesTouchingUndertrainedPair",
                    undertrainedPair);
            writeFrequencySeverityKey(writer, "tailSamplesTouchingUnseenUnary",
                    unseenUnaryTail);
            writeFrequencySeverityKey(writer, "tailSamplesTouchingUnseenPair",
                    unseenPairTail);
            writeFrequencySeverityKey(writer, "preliminaryPass", passed);
            writeFrequencySeverityKey(writer, "nextAction",
                    passed ? "draw-independent-validation"
                            : round < frequencySeverityMaxRefits
                            ? "merge-discovery-and-refit" : "abort-at-refit-cap");
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write adaptive frequency/severity round artifact", ex);
        }
    }

    private void writeFrequencySeverityCoverageArtifact(
            File dir, EtaCorrections raw,
            FrequencySeverityEtaCoverage coverage, EtaCorrections selected) {
        File output = newFrequencySeverityArtifact(dir, "eta_cell_coverage.tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tcellType\tpos1\trc1\tpos2\trc2"
                    + "\trawIidCount\tdistinctContextCount"
                    + "\trawEtaKcal\tselectedEtaKcal\tundertrained");
            for (int pos = 0; pos < raw.oneBody.length; pos++) {
                for (int rc = 0; rc < raw.oneBody[pos].length; rc++) {
                    int count = raw.oneBodyCounts[pos][rc];
                    if (count == 0) continue;
                    int contexts = coverage.unary(pos, rc);
                    writer.printf(Locale.ROOT,
                            "packstar-adaptive-frequency-severity-cell-coverage-v1\tunary"
                                    + "\t%d\t%d\t-1\t-1\t%d\t%d"
                                    + "\t%.17g\t%.17g\t%s%n",
                            pos, rc, count, contexts,
                            raw.getOneBodyEta(pos, rc),
                            selected == null ? Double.NaN
                                    : selected.getOneBodyEta(pos, rc),
                            Boolean.toString(contexts < frequencySeverityMinTrainCount));
                }
            }
            for (Map.Entry<Long, double[]> entry : raw.pairSums.entrySet()) {
                long key = entry.getKey();
                int pos1 = (int) ((key >>> 48) & 0xffffL);
                int rc1 = (int) ((key >>> 32) & 0xffffL);
                int pos2 = (int) ((key >>> 16) & 0xffffL);
                int rc2 = (int) (key & 0xffffL);
                int count = (int) Math.round(entry.getValue()[1]);
                int contexts = coverage.pairDistinctContexts
                        .getOrDefault(key, 0);
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-cell-coverage-v1\tpair"
                                + "\t%d\t%d\t%d\t%d\t%d\t%d"
                                + "\t%.17g\t%.17g\t%s%n",
                        pos1, rc1, pos2, rc2, count, contexts,
                        raw.getPairEta(pos1, rc1, pos2, rc2),
                        selected == null ? Double.NaN
                                : selected.getPairEta(
                                pos1, rc1, pos2, rc2),
                        Boolean.toString(contexts < frequencySeverityMinTrainCount));
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity cell coverage", ex);
        }
    }

    private void writeFrequencySeverityCoverageArtifactQuietly(
            File dir, EtaCorrections raw,
            FrequencySeverityEtaCoverage coverage, EtaCorrections selected) {
        try {
            writeFrequencySeverityCoverageArtifact(dir, raw, coverage, selected);
        } catch (RuntimeException ex) {
            System.out.println("[PACK*-frequency-severity] WARNING: "
                    + sanitizeTsv(ex.getMessage()));
        }
    }

    private void writeFrequencySeveritySamplesArtifact(
            File dir, String stage, List<CCDResult> samples,
            EtaCorrections eta, double logMuTrain) {
        File output = newFrequencySeverityArtifact(
                dir, "frequency_severity_" + stage + "_samples.tsv");
        double logClip = frequencySeverityRelativeBoundKcal / RT;
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tstage\tindex\tconf\teTrueKcal\teMinKcal"
                    + "\tetaKcal\teProposalKcal\tlogRawWeight"
                    + "\tlogMuTrain\tlogRelativeWeight\tbulkY\ttailEvent");
            for (int index = 0; index < samples.size(); index++) {
                CCDResult sample = samples.get(index);
                double etaEnergy = computeEtaEnergy(sample.conf, eta);
                double eProposal = sample.eMin + etaEnergy;
                double logRawWeight = (eProposal - sample.eTrue) / RT;
                double logRelative = logRawWeight - logMuTrain;
                double bulk = logRelative >= logClip
                        ? 1.0 : Math.exp(logRelative - logClip);
                writer.printf(Locale.ROOT,
                        "packstar-adaptive-frequency-severity-sample-v1"
                                + "\t%s\t%d\t%s\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%s%n",
                        stage, index, formatFrequencySeverityConf(sample.conf),
                        sample.eTrue, sample.eMin, etaEnergy,
                        eProposal, logRawWeight, logMuTrain,
                        logRelative, bulk,
                        Boolean.toString(logRelative > logClip));
            }
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity " + stage
                            + " sample artifact", ex);
        }
    }

    private String formatFrequencySeverityConf(int[] conf) {
        StringBuilder text = new StringBuilder();
        for (int index = 0; index < conf.length; index++) {
            if (index > 0) text.append(',');
            text.append(conf[index]);
        }
        return text.toString();
    }

    private void writeFrequencySeverityStageArtifact(
            File dir, String stage,
            PackStarFrequencySeverityPAC.Interval interval,
            PackStarFrequencySeverityPAC.SeverityTest severityTest,
            PackStarFrequencySeverityPAC.Sizing sizing,
            double logMuTrain, double logZCorrected) {
        File output = newFrequencySeverityArtifact(dir, "frequency_severity_" + stage + ".tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tstage\tn\ttailCount\tbulkMean"
                    + "\tbulkVariance\tbulkRadius\tbulkLower\tbulkUpper"
                    + "\ttailProbabilityEmpirical\ttailProbabilityUpper"
                    + "\tseverityCapS0\tinducedTailMeanUpper"
                    + "\tnormalizedMeanLower\tnormalizedMeanUpper\tepsilon"
                    + "\tempiricalTailMean\tempiricalConditionalSeverity"
                    + "\tobservedMaxConditionalSeverity"
                    + "\tseverityTestTailCount\tseverityTestSufficient"
                    + "\tseverityTestRejected\tseverityTestLogE"
                    + "\tseverityTestPUpper\tpredictedFinalN"
                    + "\treachableAtMax\tpredictedEpsilonAtMax"
                    + "\tlogMuTrain\tlogZProposal\tseverityPremiseId");
            writer.printf(Locale.ROOT,
                    "packstar-adaptive-frequency-severity-stage-v1\t%s\t%d\t%d"
                            + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
                            + "\t%.17g\t%.17g\t%.17g\t%.17g"
                            + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
                            + "\t%.17g\t%d\t%s\t%s\t%.17g\t%.17g"
                            + "\t%d\t%s\t%.17g\t%.17g\t%.17g\t%s%n",
                    stage, interval.sampleCount, interval.tailCount,
                    interval.bulkMean, interval.bulkVariance,
                    interval.bulkRadius, interval.bulkLower,
                    interval.bulkUpper,
                    interval.tailProbabilityEmpirical,
                    interval.tailProbabilityUpper,
                    interval.conditionalSeverityCap,
                    interval.inducedTailMeanUpper,
                    interval.normalizedMeanLower,
                    interval.normalizedMeanUpper,
                    interval.epsilon,
                    interval.empiricalTailMean,
                    interval.empiricalConditionalSeverity,
                    interval.observedMaxConditionalSeverity,
                    severityTest == null ? -1 : severityTest.tailCount,
                    severityTest == null ? "NA"
                            : Boolean.toString(
                            severityTest.sufficientTailSamples),
                    severityTest == null ? "NA"
                            : Boolean.toString(severityTest.rejected),
                    severityTest == null ? Double.NaN
                            : severityTest.logEValue,
                    severityTest == null ? Double.NaN
                            : severityTest.pValueUpper,
                    sizing == null ? -1 : sizing.finalSamples,
                    sizing == null ? "NA"
                            : Boolean.toString(sizing.reachableAtMax),
                    sizing == null ? Double.NaN
                            : sizing.epsilonAtMaxSamples,
                    logMuTrain, logZCorrected,
                    sanitizeTsv(frequencySeverityPremiseId));
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity " + stage + " artifact", ex);
        }
    }

    private void writeFrequencySeverityRunSummaryArtifact(
            File dir, FrequencySeverityCandidateScore winner,
            PackStarFrequencySeverityPAC.Sizing sizing,
            PackStarFrequencySeverityPAC.Interval interval,
            int trainCount, int pilotCount, int monitorCount,
            int finalCount, boolean valid) {
        writeFrequencySeverityRunSummaryArtifact(
                dir, winner, sizing, interval,
                trainCount, 0, pilotCount, monitorCount,
                finalCount, winner.alpha == 0.0
                        && !hasTripleEta(winner.fullEta) ? 0 : 1,
                valid);
    }

    private void writeFrequencySeverityRunSummaryArtifact(
            File dir, FrequencySeverityCandidateScore winner,
            PackStarFrequencySeverityPAC.Sizing sizing,
            PackStarFrequencySeverityPAC.Interval interval,
            int trainCount, int discoveryCount, int validationCount,
            int monitorCount, int finalCount,
            int proposalDpSweeps, boolean valid) {
        File output = newFrequencySeverityArtifact(dir, "frequency_severity_final_interval.tsv");
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("key\tvalue");
            writeFrequencySeverityKey(writer, "schema",
                    "packstar-adaptive-frequency-severity-final-v3");
            writeFrequencySeverityKey(writer, "selectedCandidate", winner.id);
            writeFrequencySeverityKey(writer, "certificateValid", valid);
            writeFrequencySeverityKey(writer, "assumptionConditional", "true");
            writeFrequencySeverityKey(writer, "severityPremiseId", frequencySeverityPremiseId);
            writeFrequencySeverityKey(writer, "severityCapS0", frequencySeverityCap);
            writeFrequencySeverityKey(writer, "relativeBoundKcal", frequencySeverityRelativeBoundKcal);
            writeFrequencySeverityKey(writer, "trainCcd", trainCount);
            writeFrequencySeverityKey(writer, "discoveryCcd", discoveryCount);
            writeFrequencySeverityKey(writer, "validationCcd", validationCount);
            writeFrequencySeverityKey(writer, "monitorCcd", monitorCount);
            writeFrequencySeverityKey(writer, "finalCcd", finalCount);
            writeFrequencySeverityKey(writer, "totalCcd", totalCCDCalls);
            writeFrequencySeverityKey(writer, "candidateDpSweeps", 0);
            writeFrequencySeverityKey(writer, "proposalDpSweeps",
                    proposalDpSweeps);
            writeFrequencySeverityKey(writer, "tripleEtaEnabled",
                    frequencySeverityTripleEta != null);
            writeFrequencySeverityKey(writer,
                    "tripleEtaCapacityFallback",
                    frequencySeverityTripleEtaFallbackReason != null);
            writeFrequencySeverityKey(writer,
                    "tripleEtaCapacityFallbackReason",
                    frequencySeverityTripleEtaFallbackReason == null
                            ? "NA" : frequencySeverityTripleEtaFallbackReason);
            writeFrequencySeverityKey(writer, "selectedTripleEtaActive",
                    hasTripleEta(winner.fullEta));
            writeFrequencySeverityKey(writer, "tripleEtaScale",
                    winner.fullEta.tripleEtaScale);
            writeFrequencySeverityKey(writer,
                    "tripleEtaPartialCcdAssignments",
                    frequencySeverityTripleEta == null ? 0L
                            : frequencySeverityTripleEta.factorAssignments);
            writeFrequencySeverityKey(writer,
                    "tripleEtaPositiveJointFactors",
                    frequencySeverityTripleEta == null ? 0L
                            : frequencySeverityTripleEta.positiveFactors);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaPositionTriples",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.positionTripleCount : 0L);
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaPositionScopes",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.positionTriples()
                            : "[]");
            writeFrequencySeverityKey(writer,
                    "selectedTripleEtaFillEdges",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.requiredFillEdges(
                            interactionGraph).size() : 0);
            PackStarTripleEtaCorrections.ResidualSummary tripleSummary =
                    summarizeFrequencySeverityTripleResidual(
                            winner.fullEta);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualPositiveFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.positiveFactors);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualNegativeFactors",
                    tripleSummary == null ? 0L
                            : tripleSummary.negativeFactors);
            writeFrequencySeverityKey(writer,
                    "tripleEtaResidualMaximumAbsoluteKcal",
                    tripleSummary == null ? 0.0
                            : tripleSummary.maximumAbsoluteCorrectionKcal);
            writeFrequencySeverityKey(writer,
                    "tripleEtaSignatureSha256",
                    hasTripleEta(winner.fullEta)
                            ? winner.fullEta.tripleEta.signatureSha256 : "NA");
            writeFrequencySeverityKey(writer, "predictedReachableAtMax",
                    sizing.reachableAtMax);
            writeFrequencySeverityKey(writer, "predictedEpsilonAtMax",
                    sizing.epsilonAtMaxSamples);
            writeFrequencySeverityKey(writer, "tailCount", interval.tailCount);
            writeFrequencySeverityKey(writer, "tailProbabilityUpper",
                    interval.tailProbabilityUpper);
            writeFrequencySeverityKey(writer, "bulkLower", interval.bulkLower);
            writeFrequencySeverityKey(writer, "bulkUpper", interval.bulkUpper);
            writeFrequencySeverityKey(writer, "epsilon", epsilon);
            writeFrequencySeverityKey(writer, "targetEpsilon", targetEpsilon);
            writeFrequencySeverityKey(writer, "targetReached",
                    epsilon <= targetEpsilon + 1.0e-12);
            writeFrequencySeverityKey(writer, "logZLower", logZLowerPAC);
            writeFrequencySeverityKey(writer, "logZUpper", logZUpperPAC);
        } catch (Exception ex) {
            throw new IllegalStateException(
                    "cannot write frequency/severity final interval artifact", ex);
        }
    }

    private void writeFrequencySeverityFailureArtifactQuietly(
            File dir, String reason, int trainCount,
            int pilotCount, int monitorCount, int finalCount) {
        try {
            File output = new File(dir, "frequency_severity_failure.tsv");
            if (output.exists()) return;
            try (PrintWriter writer = new PrintWriter(
                    output, StandardCharsets.UTF_8.name())) {
                writer.println("key\tvalue");
                writeFrequencySeverityKey(writer, "schema", "packstar-adaptive-frequency-severity-failure-v1");
                writeFrequencySeverityKey(writer, "reason", sanitizeTsv(reason));
                writeFrequencySeverityKey(writer, "trainCcd", trainCount);
                writeFrequencySeverityKey(writer, "pilotCcd", pilotCount);
                writeFrequencySeverityKey(writer, "monitorCcd", monitorCount);
                writeFrequencySeverityKey(writer, "finalCcd", finalCount);
                writeFrequencySeverityKey(writer, "totalCcd", totalCCDCalls);
            }
        } catch (Exception ex) {
            System.out.println("[PACK*-frequency-severity] WARNING: cannot write failure artifact: "
                    + sanitizeTsv(ex.getMessage()));
        }
    }

    /** Record a typed capacity fallback while preserving the pair-only path. */
    private void writeFrequencySeverityTripleEtaFallbackArtifactQuietly(
            File dir,
            PackStarTripleEtaCorrections.AssignmentCapExceededException ex) {
        try {
            File output = newFrequencySeverityArtifact(
                    dir, "triple_eta_capacity_fallback.tsv");
            try (PrintWriter writer = new PrintWriter(
                    output, StandardCharsets.UTF_8.name())) {
                writer.println("key\tvalue");
                writeFrequencySeverityKey(writer, "schema",
                        "packstar-adaptive-frequency-severity-triple-cap-fallback-v1");
                writeFrequencySeverityKey(writer, "fallbackModel", "pair-only");
                writeFrequencySeverityKey(writer, "reason", ex.getMessage());
                writeFrequencySeverityKey(writer, "assignments", ex.assignments);
                writeFrequencySeverityKey(writer,
                        "maximumAssignments", ex.maximumAssignments);
                writeFrequencySeverityKey(writer, "certificatePolicy",
                        "optional-triple-model-rejected; pair-only-candidate-remains-eligible");
            }
        } catch (Exception artifactError) {
            System.out.println(
                    "[PACK*-frequency-severity] WARNING: cannot write"
                            + " triple eta capacity fallback artifact: "
                            + sanitizeTsv(artifactError.getMessage()));
        }
    }

    /**
     * Write a forensic decomposition for an observed q_m lower-bound failure.
     *
     * <p>This is deliberately diagnostic only.  The caller has already made
     * the certificate decision, and this method must never relax it.  The
     * summary compares the sparse minimizing-matrix sum with the energy of the
     * same minimized pose.  The term table then exposes the one-body and pair
     * contributions in both representations, together with local/global RC
     * identifiers.  It is emitted only on an invalid q_m audit so normal
     * production runs do not pay for a second force-field decomposition.</p>
     */
    private void writeFrequencySeverityLowerBoundViolationAuditQuietly(
            File dir, String stage, List<CCDResult> samples) {
        if (dir == null || samples == null) return;
        try {
            File summary = newFrequencySeverityArtifact(
                    dir, "q_m_lower_bound_violation_audit.tsv");
            File terms = newFrequencySeverityArtifact(
                    dir, "q_m_lower_bound_violation_terms.tsv");

            try (PrintWriter summaryWriter = new PrintWriter(
                    summary, StandardCharsets.UTF_8.name());
                 PrintWriter termWriter = new PrintWriter(
                         terms, StandardCharsets.UTF_8.name())) {

                summaryWriter.println("schema\tstage\tsampleIndex\tstatus"
                        + "\teTrue\teMin\tgapTrueMinusMin\tmatrixConst"
                        + "\tmatrixOneBodySum\tmatrixGraphPairSum"
                        + "\tmatrixTotal\tposeOneBodySum\tposePairSum"
                        + "\tposeTotal\ttrueMinusPoseTotal\tgraphEdges"
                        + "\tcutPairs\tinteractionCount\tconf");
                termWriter.println("schema\tstage\tsampleIndex\ttermType"
                        + "\tpos1\tpos1Name\tlocalRc1\tglobalRc1"
                        + "\tconfType1\tpos2\tpos2Name\tlocalRc2"
                        + "\tglobalRc2\tconfType2\tgraphEdge"
                        + "\tmatrixEnergy\tposeEnergySparse"
                        + "\tposeEnergyFullPair\tdifference");

                for (int sampleIndex = 0;
                     sampleIndex < samples.size(); sampleIndex++) {
                    CCDResult sample = samples.get(sampleIndex);
                    double gap = sample.eTrue - sample.eMin;
                    boolean nonFinite = !Double.isFinite(sample.eTrue)
                            || !Double.isFinite(sample.eMin);
                    if (!nonFinite && !(gap < 0.0)) continue;

                    EnergyMatrix poseBreakdown = null;
                    String status = nonFinite ? "nonfinite" : "violation";
                    try {
                        if (sample.epmol != null) {
                            poseBreakdown = new ResidueForcefieldBreakdown.ByPosition(
                                    minimizingEcalc, sample.conf, sample.epmol)
                                    .breakdownForcefield(
                                            ResidueForcefieldBreakdown.Type.All);
                        } else {
                            status += ":missing-pose";
                        }
                    } catch (RuntimeException ex) {
                        status += ":breakdown-error-"
                                + ex.getClass().getSimpleName();
                    }

                    double matrixConst = branchMinimizingEmat.getConstTerm();
                    double matrixOneBodySum = 0.0;
                    double matrixGraphPairSum = 0.0;
                    double poseOneBodySum = 0.0;
                    double posePairSum = 0.0;
                    int graphEdges = 0;
                    int cutPairs = 0;

                    for (int pos = 0; pos < sample.conf.length; pos++) {
                        int globalRc = sample.conf[pos];
                        if (globalRc < 0) continue;
                        int localRc = localRCByGlobalRC[pos].get(globalRc);
                        double matrixEnergy = branchMinimizingEmat
                                .getOneBody(pos, globalRc);
                        double poseEnergy = poseBreakdown == null
                                ? Double.NaN
                                : poseBreakdown.getOneBody(pos, 0);
                        matrixOneBodySum += matrixEnergy;
                        if (Double.isFinite(poseEnergy)) {
                            poseOneBodySum += poseEnergy;
                        }
                        termWriter.printf(Locale.ROOT,
                                "packstar-qm-lower-bound-audit-v1\t%s\t%d\tsingle"
                                        + "\t%d\t%s\t%d\t%d\t%s\t-1\tNA\t-1\t-1\tNA"
                                        + "\tfalse\t%.17g\t%.17g\t%.17g\t%.17g%n",
                                sanitizeTsv(stage), sampleIndex, pos,
                                sanitizeTsv(confSpace.name(pos)), localRc,
                                globalRc,
                                sanitizeTsv(confSpace.confType(pos, globalRc)),
                                matrixEnergy, poseEnergy, poseEnergy,
                                poseEnergy - matrixEnergy);
                    }

                    for (int pos1 = 0; pos1 < sample.conf.length; pos1++) {
                        int globalRc1 = sample.conf[pos1];
                        if (globalRc1 < 0) continue;
                        int localRc1 = localRCByGlobalRC[pos1].get(globalRc1);
                        for (int pos2 = pos1 + 1;
                             pos2 < sample.conf.length; pos2++) {
                            int globalRc2 = sample.conf[pos2];
                            if (globalRc2 < 0) continue;
                            boolean graphEdge = interactionGraph.hasEdge(
                                    pos1, pos2);
                            if (graphEdge) graphEdges++;
                            else cutPairs++;
                            int localRc2 = localRCByGlobalRC[pos2]
                                    .get(globalRc2);
                            double matrixEnergy = graphEdge
                                    ? branchMinimizingEmat.getPairwise(
                                    pos1, globalRc1, pos2, globalRc2)
                                    : 0.0;
                            double poseEnergyFullPair = poseBreakdown == null
                                    ? Double.NaN
                                    : poseBreakdown.getPairwise(
                                    pos1, 0, pos2, 0);
                            double poseEnergy = graphEdge
                                    ? poseEnergyFullPair : 0.0;
                            if (graphEdge) {
                                matrixGraphPairSum += matrixEnergy;
                            }
                            if (Double.isFinite(poseEnergy)) {
                                posePairSum += poseEnergy;
                            }
                            termWriter.printf(Locale.ROOT,
                                    "packstar-qm-lower-bound-audit-v1\t%s\t%d\tpair"
                                            + "\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s"
                                            + "\t%s\t%.17g\t%.17g\t%.17g%n",
                                    sanitizeTsv(stage), sampleIndex,
                                    pos1, sanitizeTsv(confSpace.name(pos1)),
                                    localRc1, globalRc1,
                                    sanitizeTsv(confSpace.confType(pos1,
                                            globalRc1)),
                                    pos2, sanitizeTsv(confSpace.name(pos2)),
                                    localRc2, globalRc2,
                                    sanitizeTsv(confSpace.confType(pos2,
                                            globalRc2)),
                                    Boolean.toString(graphEdge), matrixEnergy,
                                    poseEnergy, poseEnergyFullPair,
                                    poseEnergy - matrixEnergy);
                        }
                    }

                    double matrixTotal = matrixConst + matrixOneBodySum
                            + matrixGraphPairSum;
                    double poseTotal = poseOneBodySum + posePairSum;
                    int interactionCount = sample.epmol == null
                            || sample.epmol.inters == null
                            ? -1 : sample.epmol.inters.size();
                    summaryWriter.printf(Locale.ROOT,
                            "packstar-qm-lower-bound-audit-v1\t%s\t%d\t%s"
                                    + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
                                    + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
                                    + "\t%.17g\t%d\t%d\t%d\t%s%n",
                            sanitizeTsv(stage), sampleIndex, status,
                            sample.eTrue, sample.eMin, gap, matrixConst,
                            matrixOneBodySum, matrixGraphPairSum, matrixTotal,
                            poseOneBodySum, posePairSum, poseTotal,
                            sample.eTrue - poseTotal, graphEdges, cutPairs,
                            interactionCount,
                            sanitizeTsv(Arrays.toString(sample.conf)));
                }

                if (summaryWriter.checkError() || termWriter.checkError()) {
                    throw new IOException("write failed for q_m lower-bound audit");
                }
            }
        } catch (Exception ex) {
            System.out.println("[PACK*-frequency-severity] WARNING: cannot write q_m lower-bound audit: "
                    + sanitizeTsv(ex.getMessage()));
        }
    }

    private File newFrequencySeverityArtifact(File dir, String name) {
        File output = new File(dir, name);
        if (output.exists()) {
            throw new IllegalStateException(
                    "frequency/severity route refuses to overwrite artifact " + output);
        }
        return output;
    }

    private String frequencySeverityPrefixedArtifactName(
            String prefix, String baseName) {
        if (prefix == null || prefix.isEmpty()) return baseName;
        return sanitizeTsv(prefix).replace('\t', '_') + "_" + baseName;
    }

    private void writeFrequencySeverityKey(PrintWriter writer, String key, Object value) {
        writer.println(sanitizeTsv(key) + "\t" + sanitizeTsv(
                String.valueOf(value)));
    }

    private String formatFrequencySeverityShrinkGrid() {
        StringBuilder text = new StringBuilder();
        for (FrequencySeverityShrinkPair shrink : frequencySeverityShrinkGrid) {
            if (text.length() > 0) text.append(',');
            text.append(String.format(Locale.ROOT, "%g:%g",
                    shrink.unary, shrink.pair));
        }
        return text.toString();
    }

    /**
     * Practical conditional certificate selected for production.
     *
     * <p>All proposal adaptation, eta repair, B selection, clip selection, and
     * final sample sizing are functions of train/adaptation data only.  Eta,
     * p_eta, B, the clip threshold, and N are frozen before an independent
     * monitor and the fresh final estimation set are drawn.  The monitor can
     * falsify but cannot prove the global support premise.  Conditional on the
     * selected pointwise bound holding, the final bounded-weight interval has
     * the advertised PAC coverage.</p>
     */
    private void runConditionalRepairAuditPAC(long startTime, double logZRigid) {
        if (frequencySeverityProduction) {
            runAdaptiveFrequencySeverityPAC(startTime, logZRigid);
            return;
        }
        System.out.println("[PACK*-conditional] eta-repair + capped dynamic-B certificate enabled"
                + "; B_grid=" + Arrays.toString(residualBoundGridKcal)
                + "; monitorSamples=" + monitorSamples
                + "; dpRich=" + dpRich
                + "; etaV3=" + etaV3
                + "; etaV4=" + etaV4
                + "; observed violation policy=abort"
                + "; q_m/search fallback=false");

        int adaptationCCDCalls = 0;
        int monitorCCDCalls = 0;
        int finalCCDCalls = 0;

        // ---- Stage A: p_m training ----
        List<CCDResult> trainCCD = runParallelCCD(sampleConformationsFromDP(
                trainSamples, stageRandom("conditional-train")));
        totalCCDCalls = trainCCD.size();
        if (trainCCD.isEmpty()
                || !validateObservedLowerBound(trainCCD, "conditional-train")) {
            failCertificate("conditional: invalid train sample or observed E_m lower-bound violation",
                    startTime);
            return;
        }

        List<CCDResult> trainPool = new ArrayList<>(trainCCD);
        EtaCorrections eta = etaEnabled
                ? extractEtaCorrections(trainPool)
                : zeroEtaCorrections();
        EtaCorrections dpRichBaseEta = eta;
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        maybeDumpTrainingSamples(trainCCD, eta, correctedEmat);
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta);
        double logZCorrected = correctedDP.logZCorrected;
        if (!Double.isFinite(logZCorrected)) {
            failCertificate("conditional: initial corrected DP normalizer is non-finite",
                    startTime);
            return;
        }
        System.out.println("[PACK*-conditional] corrected DP: logZ_eta="
                + String.format(Locale.ROOT, "%.6f", logZCorrected)
                + ", logZ_eta-logZ_rigid="
                + String.format(Locale.ROOT, "%.6f", logZCorrected - logZRigid)
                + ", oneBodyTerms=" + eta.oneBodyCount
                + ", pairTerms=" + eta.pairCount);

        // ---- Stage B: adaptation/refit/repair only ----
        boolean iterate = etaEnabled
                && getConfigBoolean("packstar.pac.iterate", DEFAULT_ITERATE);
        int iterateMaxRounds = iterate
                ? Math.max(0, getConfigInteger("packstar.pac.iterate.maxRounds",
                DEFAULT_ITERATE_MAX_ROUNDS)) : 0;
        int repairMaxRounds = etaEnabled && etaRepairEnabled
                ? etaRepairMaxRounds : 0;
        int adaptationMaxRounds = Math.max(iterateMaxRounds, repairMaxRounds);
        double collapseThresh = getConfigDouble(
                "packstar.pac.iterate.meanWThreshold", 0.3);
        double driftFracThresh = getConfigDouble(
                "packstar.pac.iterate.driftFraction", 0.2);
        int minTrainCount = getConfigInteger(
                "packstar.pac.iterate.minTrainCount", 5);
        int repairApplications = 0;
        boolean repairActive = false;
        double finalRepairAlpha = 1.0;
        List<CCDResult> repairWitnessPool = new ArrayList<>();
        List<ConditionalProposalSnapshot> proposalArchive = new ArrayList<>();

        int adaptationCount = Math.min(pilotSamples, maxEstSamples);
        String adaptationSourceId = "conditional-adapt-0";
        proposalArchive.add(new ConditionalProposalSnapshot(
                adaptationSourceId, eta, correctedEmat,
                logZCorrected, finalRepairAlpha));
        List<CCDResult> adaptationPilot = runParallelCCD(sampleConformationsFromDP(
                adaptationCount, stageRandom(adaptationSourceId)),
                correctedEmat, logZCorrected, adaptationSourceId);
        List<CCDResult> adaptationPool = new ArrayList<>(adaptationPilot);
        adaptationCCDCalls += adaptationPilot.size();
        totalCCDCalls = trainCCD.size() + adaptationCCDCalls;
        if (adaptationPilot.isEmpty()
                || !validateObservedLowerBound(adaptationPilot,
                "conditional-adapt-0")) {
            failCertificate("conditional: invalid initial adaptation pilot",
                    startTime);
            return;
        }

        for (int round = 1; round <= adaptationMaxRounds; round++) {
            double meanWPilot = estimateMeanW(adaptationPilot, correctedEmat);
            double driftFrac = pilotUndertrainedFraction(
                    adaptationPilot, eta, minTrainCount);
            double maxOverCorrection = maxOverCorrectionKcal(
                    adaptationPilot, correctedEmat);
            boolean refineRequested = round <= iterateMaxRounds
                    && ((Double.isFinite(meanWPilot)
                    && meanWPilot < collapseThresh)
                    || driftFrac >= driftFracThresh);
            boolean repairRequested = repairApplications < repairMaxRounds
                    && Double.isFinite(maxOverCorrection)
                    && maxOverCorrection > etaRepairTriggerKcal;

            System.out.println("[PACK*-conditional-repair] check round "
                    + round + "/" + adaptationMaxRounds
                    + ": meanW="
                    + String.format(Locale.ROOT, "%.6f", meanWPilot)
                    + ", driftFrac="
                    + String.format(Locale.ROOT, "%.4f", driftFrac)
                    + ", maxOverCorrectionKcal="
                    + String.format(Locale.ROOT, "%.6f", maxOverCorrection)
                    + ", refit=" + refineRequested
                    + ", repair=" + repairRequested
                    + ", repairApplications=" + repairApplications
                    + "/" + repairMaxRounds);
            if (!(refineRequested || repairRequested)) break;

            trainPool.addAll(adaptationPilot);
            if (repairRequested) {
                repairActive = true;
                repairApplications++;
            }
            if (repairActive) repairWitnessPool.addAll(adaptationPilot);

            EtaCorrections candidateEta = extractEtaCorrections(trainPool);
            dpRichBaseEta = candidateEta;
            EnergyMatrix candidateEmat = buildCorrectedEmat(candidateEta);
            double alpha = repairActive
                    ? chooseEtaRepairAlpha(repairWitnessPool, candidateEmat)
                    : 1.0;
            if (!Double.isFinite(alpha)) {
                failCertificate("conditional: eta repair encountered an observed E_m lower-bound violation",
                        startTime);
                return;
            }
            if (alpha < 1.0) {
                candidateEta = candidateEta.scaled(alpha);
                candidateEmat = buildCorrectedEmat(candidateEta);
            }
            CorrectedDPResult candidateDP = recomputeDP(
                    candidateEmat, candidateEta);
            if (!Double.isFinite(candidateDP.logZCorrected)) {
                failCertificate("conditional: eta repair produced a non-finite corrected DP",
                        startTime);
                return;
            }
            eta = candidateEta;
            correctedEmat = candidateEmat;
            logZCorrected = candidateDP.logZCorrected;
            finalRepairAlpha = alpha;

            adaptationSourceId = "conditional-adapt-" + round;
            proposalArchive.add(new ConditionalProposalSnapshot(
                    adaptationSourceId, eta, correctedEmat,
                    logZCorrected, finalRepairAlpha));
            adaptationPilot = runParallelCCD(sampleConformationsFromDP(
                    adaptationCount,
                    stageRandom(adaptationSourceId)),
                    correctedEmat, logZCorrected,
                    adaptationSourceId);
            adaptationPool.addAll(adaptationPilot);
            adaptationCCDCalls += adaptationPilot.size();
            totalCCDCalls = trainCCD.size() + adaptationCCDCalls;
            if (adaptationPilot.isEmpty()
                    || !validateObservedLowerBound(adaptationPilot,
                    "conditional-adapt-" + round)) {
                failCertificate("conditional: invalid fresh adaptation pilot after repair",
                        startTime);
                return;
            }
            System.out.println("[PACK*-conditional-repair] applied round "
                    + round + ": alpha="
                    + String.format(Locale.ROOT, "%.6f", alpha)
                    + ", trainPool=" + trainPool.size()
                    + ", repairWitnesses=" + repairWitnessPool.size()
                    + ", freshPilot=" + adaptationPilot.size()
                    + ", logZ_eta="
                    + String.format(Locale.ROOT, "%.6f", logZCorrected));
        }

        if (etaV4) {
            runConditionalEtaV4MixturePAC(
                    startTime, trainCCD.size(), adaptationCCDCalls,
                    proposalArchive, adaptationPool);
            return;
        }

        // ---- Stage C: adaptation-only DP-rich candidate screen, then freeze ----
        double adaptationMaxOverCorrection;
        int dpRichSweeps = 0;
        int dpRichCalibrationCCD = 0;
        String dpRichSource = "disabled-current-eta";
        ConditionalAdaptationStats dpRichAdaptation = null;
        if (dpRich && etaEnabled) {
            ConditionalDpRichDesign design = etaV3
                    ? selectConditionalEtaV3Design(
                    proposalArchive, adaptationPool)
                    : selectConditionalDpRichDesign(
                    dpRichBaseEta, adaptationPool, finalRepairAlpha);
            if (design == null) {
                failCertificate("conditional: no eligible DP-rich eta/B candidate survived adaptation",
                        startTime);
                return;
            }
            eta = design.eta;
            correctedEmat = design.emat;
            logZCorrected = design.logZ;
            finalRepairAlpha = design.alpha;
            selectedResidualBoundKcal = design.selectedBoundKcal;
            adaptationMaxOverCorrection =
                    design.adaptationMaxOverCorrectionKcal;
            dpRichSweeps = design.candidateDpSweeps;
            dpRichSource = design.source;
            List<CCDResult> winnerAdaptation = adaptationPool;
            if (etaV3) {
                String calibrationSourceId = "conditional-v3-calibration";
                List<CCDResult> calibration = runParallelCCD(
                        sampleConformationsFromDP(
                                etaV3CalibrationSamples,
                                stageRandom(calibrationSourceId)),
                        correctedEmat, logZCorrected,
                        calibrationSourceId);
                dpRichCalibrationCCD = calibration.size();
                adaptationCCDCalls += calibration.size();
                totalCCDCalls = trainCCD.size() + adaptationCCDCalls;
                if (calibration.isEmpty()
                        || !validateObservedLowerBound(calibration,
                        calibrationSourceId)) {
                    failCertificate("conditional: invalid eta-v3 winner calibration batch",
                            startTime);
                    return;
                }
                adaptationPool.addAll(calibration);
                double calibratedMax = maxOverCorrectionKcal(
                        adaptationPool, correctedEmat);
                double targetMax = selectedResidualBoundKcal
                        - etaV3HeadroomKcal;
                double additionalShiftDown = Math.max(0.0,
                        calibratedMax - targetMax);
                if (additionalShiftDown > 0.0) {
                    eta = eta.shiftedConstant(-additionalShiftDown);
                    correctedEmat = buildCorrectedEmat(eta);
                    CorrectedDPResult calibratedDP = recomputeDP(
                            correctedEmat, eta);
                    dpRichSweeps++;
                    if (!Double.isFinite(calibratedDP.logZCorrected)) {
                        failCertificate("conditional: eta-v3 gauge calibration produced a non-finite DP",
                                startTime);
                        return;
                    }
                    logZCorrected = calibratedDP.logZCorrected;
                    calibratedMax -= additionalShiftDown;
                    System.out.println("[PACK*-conditional-eta-v3] calibration gauge:"
                            + " additionalShiftDownKcal="
                            + String.format(Locale.ROOT, "%.6f",
                            additionalShiftDown)
                            + ", calibratedMaxOverCorrectionKcal="
                            + String.format(Locale.ROOT, "%.6f",
                            calibratedMax)
                            + ", targetMaxKcal="
                            + String.format(Locale.ROOT, "%.6f",
                            targetMax)
                            + ", proposalDistributionChanged=false");
                }
                adaptationMaxOverCorrection = calibratedMax;
                // This is an exact on-policy batch from the frozen winner.
                // Earlier samples remain support witnesses but are not needed
                // for winner sizing, so no cross-proposal ESS approximation is
                // introduced after calibration.
                winnerAdaptation = calibration;
            }
            // Frequency/severity validation is deliberately separated from the legacy
            // pointwise-B certificate gates below.  The adaptation-only
            // proposal is frozen at this point and its exact DP tables are
            // active, so this batch is fresh IID from that proposal.  Full
            // CCD supplies the true energy used by the offline frequency/severity analysis.
            // The mode writes diagnostics and then stops without claiming a
            // partition-function certificate.
            if (shouldWriteMomentTailAudit("frequency-severity-validation")) {
                String validationSourceId = "frequency-severity-validation";
                List<CCDResult> validationCCD = runParallelCCD(
                        sampleConformationsFromDP(
                                monitorSamples,
                                stageRandom(validationSourceId)),
                        correctedEmat, logZCorrected,
                        validationSourceId);
                monitorCCDCalls = validationCCD.size();
                totalCCDCalls = trainCCD.size() + adaptationCCDCalls
                        + monitorCCDCalls;
                maybeWritePairwiseMomentTailAudit(
                        validationCCD, correctedEmat,
                        validationSourceId);
                boolean complete = validationCCD.size() == monitorSamples;
                boolean lowerBoundValid = !validationCCD.isEmpty()
                        && validateObservedLowerBound(
                        validationCCD, validationSourceId);
                String reason = complete && lowerBoundValid
                        ? "diagnostic-only frequency/severity validation complete; no certificate requested"
                        : "diagnostic-only frequency/severity validation incomplete or observed a minimizing lower-bound violation";
                failCertificate(reason, startTime);
                return;
            }
            dpRichAdaptation = buildConditionalAdaptationStats(
                    winnerAdaptation, correctedEmat, logZCorrected,
                    selectedResidualBoundKcal);
            if (dpRichAdaptation == null) {
                failCertificate("conditional: DP-rich winner adaptation reweighting is invalid",
                        startTime);
                return;
            }
            double essThreshold = Math.max(2.0,
                    dpRichMinEssFraction
                            * dpRichAdaptation.logResidualWeights.length);
            boolean essPass = conditionalDpRichEssPasses(
                    dpRichAdaptation.effectiveSampleSize,
                    dpRichAdaptation.logResidualWeights.length,
                    dpRichMinEssFraction);
            System.out.println("[PACK*-conditional-dp-rich] winner reweighting:"
                    + " adaptation="
                    + dpRichAdaptation.logResidualWeights.length
                    + ", sourceProposals="
                    + dpRichAdaptation.sourceProposalCount
                    + ", ESS="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichAdaptation.effectiveSampleSize)
                    + ", ESSFraction="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichAdaptation.effectiveSampleSize
                            / dpRichAdaptation.logResidualWeights.length)
                    + ", minESSFraction="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichMinEssFraction)
                    + ", ESSThreshold="
                    + String.format(Locale.ROOT, "%.6f", essThreshold)
                    + ", minLogRho="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichAdaptation.minLogProposalWeight)
                    + ", maxLogRho="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichAdaptation.maxLogProposalWeight)
                    + ", pass=" + essPass
                    + ", extraCandidateCCD="
                    + dpRichCalibrationCCD);
            if (!essPass) {
                failCertificate("conditional: DP-rich winner adaptation reweighting ESS is below the pre-registered threshold",
                        startTime);
                return;
            }
        } else {
            adaptationMaxOverCorrection = maxOverCorrectionKcal(
                    adaptationPilot, correctedEmat);
            selectedResidualBoundKcal = selectSmallestResidualBound(
                    adaptationMaxOverCorrection, residualBoundGridKcal);
            if (!Double.isFinite(selectedResidualBoundKcal)) {
                failCertificate("conditional: adaptation max over-correction "
                        + String.format(Locale.ROOT, "%.6f", adaptationMaxOverCorrection)
                        + " kcal/mol exceeds capped B grid "
                        + Arrays.toString(residualBoundGridKcal), startTime);
                return;
            }
        }
        clipLogCap = dpRichAdaptation != null
                ? selectConditionalClipLogCap(
                dpRichAdaptation, selectedResidualBoundKcal)
                : selectConditionalClipLogCap(
                adaptationPilot, correctedEmat,
                selectedResidualBoundKcal);
        if (!Double.isFinite(clipLogCap)) {
            failCertificate("conditional: failed to freeze a finite clip threshold",
                    startTime);
            return;
        }
        if (etaV3 && dpRichAdaptation != null) {
            ConditionalSizing calibratedSizing = computeConditionalSizing(
                    dpRichAdaptation, selectedResidualBoundKcal,
                    clipLogCap);
            if (calibratedSizing == null
                    || calibratedSizing.epsilonAtMaxSamples
                    > targetEpsilon + 1.0e-12) {
                failCertificate("conditional: eta-v3 calibrated winner is predicted unreachable at the final sample cap",
                        startTime);
                return;
            }
        }
        int nStar = dpRichAdaptation != null
                ? computeRequiredConditionalEstSamples(
                dpRichAdaptation, selectedResidualBoundKcal, clipLogCap)
                : computeRequiredConditionalEstSamples(
                adaptationPilot, correctedEmat,
                selectedResidualBoundKcal, clipLogCap);
        if (nStar < 2) {
            failCertificate("conditional: failed to select a valid final sample count",
                    startTime);
            return;
        }
        System.out.println("[PACK*-conditional] frozen design: repairAlpha="
                + String.format(Locale.ROOT, "%.6f", finalRepairAlpha)
                + ", repairApplications=" + repairApplications
                + ", dpRichSource=" + dpRichSource
                + ", dpRichSweeps=" + dpRichSweeps
                + ", extraCandidateCCD=" + dpRichCalibrationCCD
                + ", adaptationPool=" + adaptationPool.size()
                + (dpRichAdaptation != null
                ? ", adaptationESS="
                    + String.format(Locale.ROOT, "%.6f",
                    dpRichAdaptation.effectiveSampleSize)
                    + ", adaptationSourceProposals="
                    + dpRichAdaptation.sourceProposalCount
                : "")
                + ", adaptationMaxOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f",
                adaptationMaxOverCorrection)
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                selectedResidualBoundKcal)
                + ", clipLogCap="
                + String.format(Locale.ROOT, "%.6f", clipLogCap)
                + ", monitorN=" + monitorSamples
                + ", finalN=" + nStar
                + ", q_mFallback=false, deterministicFallback=false");

        // ---- Stage D: fresh independent monitor; never used to alter eta/B/N ----
        List<CCDResult> monitorCCD = runParallelCCD(sampleConformationsFromDP(
                monitorSamples, stageRandom("conditional-monitor")));
        monitorCCDCalls = monitorCCD.size();
        totalCCDCalls = trainCCD.size() + adaptationCCDCalls
                + monitorCCDCalls;
        if (!auditConditionalSupport(monitorCCD, correctedEmat,
                selectedResidualBoundKcal, "monitor")) {
            failCertificate("conditional: independent monitor observed a support or lower-bound violation",
                    startTime);
            return;
        }

        // ---- Stage E: wholly fresh final estimation ----
        List<CCDResult> finalCCD = runParallelCCD(sampleConformationsFromDP(
                nStar, stageRandom("conditional-final")));
        finalCCDCalls = finalCCD.size();
        totalCCDCalls = trainCCD.size() + adaptationCCDCalls
                + monitorCCDCalls + finalCCDCalls;
        if (!auditConditionalSupport(finalCCD, correctedEmat,
                selectedResidualBoundKcal, "final")) {
            failCertificate("conditional: final estimation observed a support or lower-bound violation",
                    startTime);
            return;
        }

        computeConditionalRepairedBound(finalCCD, correctedEmat,
                logZCorrected, selectedResidualBoundKcal);
        if (!certificateValid) {
            if ("not computed".equals(certificateFailureReason)) {
                certificateFailureReason = "conditional repaired interval failed validation";
            }
            printFinalSummary(startTime);
            return;
        }
        printFinalSummary(startTime);
    }

    private void runConditionalEtaV4MixturePAC(
            long startTime, int trainCCDCalls,
            int priorAdaptationCCDCalls,
            List<ConditionalProposalSnapshot> proposalArchive,
            List<CCDResult> adaptationPool) {
        int adaptationCCDCalls = priorAdaptationCCDCalls;
        int monitorCCDCalls = 0;
        int finalCCDCalls = 0;

        ConditionalEtaV4Anchor anchor =
                selectConditionalEtaV4Anchor(
                        proposalArchive, adaptationPool);
        if (anchor == null) {
            failCertificate("conditional eta-v4: no valid conditional anchor",
                    startTime);
            return;
        }
        ConditionalProposalSnapshot[] components;
        List<CCDResult> mixtureScreenPool = new ArrayList<>();
        double[] referenceWeights;
        int candidateDpSweeps = anchor.candidateDpSweeps;
        ConditionalProposalSnapshot currentOaisProposal =
                anchor.proposal;
        List<CCDResult> currentOaisSamples =
                anchor.onPolicySamples;
        List<CCDResult> originalAnchorScreening = null;
        List<CCDResult> allOaisCalibration = new ArrayList<>();
        List<CCDResult> finalOaisCalibration = null;
        int acceptedOaisRounds = 0;
        double cumulativeHeldOutImprovement = 0.0;

        for (int roundIndex = 1;
             roundIndex <= etaV4OaisMaxRounds; roundIndex++) {
            ConditionalEtaV4Tail round =
                    selectConditionalEtaV4Tail(
                            currentOaisProposal,
                            currentOaisSamples, roundIndex);
            if (round == null) break;
            candidateDpSweeps += round.candidateDpSweeps;
            if (originalAnchorScreening == null) {
                originalAnchorScreening = new ArrayList<>(
                        round.baseScreeningSamples);
            }

            String calibrationStage =
                    "conditional-v4-oais-round-"
                            + roundIndex + "-calibration";
            List<CCDResult> roundCalibration = runParallelCCD(
                    sampleConformationsFromDP(
                            etaV4CalibrationSamples,
                            stageRandom(calibrationStage)),
                    round.proposal.emat, round.proposal.logZ,
                    round.proposal.sourceProposalId);
            adaptationCCDCalls += roundCalibration.size();
            totalCCDCalls = trainCCDCalls + adaptationCCDCalls;
            if (roundCalibration.isEmpty()
                    || !validateObservedLowerBound(
                    roundCalibration, calibrationStage)) {
                failCertificate("conditional eta-v4: invalid fresh OAIS on-policy calibration batch",
                        startTime);
                return;
            }
            allOaisCalibration.addAll(roundCalibration);
            finalOaisCalibration = roundCalibration;
            currentOaisProposal = round.proposal;
            currentOaisSamples = roundCalibration;
            acceptedOaisRounds++;
            cumulativeHeldOutImprovement +=
                    round.logSecondMomentImprovement;
            System.out.println("[PACK*-conditional-eta-v4] OAIS round accepted:"
                    + " round=" + roundIndex
                    + ", freshOnPolicySamples="
                    + roundCalibration.size()
                    + ", roundLogSecondMomentImprovement="
                    + String.format(Locale.ROOT, "%.9f",
                    round.logSecondMomentImprovement)
                    + ", cumulativeHeldOutImprovement="
                    + String.format(Locale.ROOT, "%.9f",
                    cumulativeHeldOutImprovement)
                    + ", adaptationCCD=" + adaptationCCDCalls);
        }

        if (etaV5ThresholdOracle) {
            System.out.println("[PACK*-eta-v5] Phase-A proposal frozen:"
                    + " source="
                    + currentOaisProposal.sourceProposalId
                    + ", acceptedOaisRounds=" + acceptedOaisRounds
                    + ", proposalLogZ="
                    + String.format(Locale.ROOT, "%.9f",
                    currentOaisProposal.logZ)
                    + ", entering exhaustive threshold oracle before"
                    + " conditional-B mixture calibration");
            runEtaV5ThresholdOracle(
                    startTime, currentOaisProposal);
            return;
        }

        if (acceptedOaisRounds == 0) {
            components = new ConditionalProposalSnapshot[]{
                    anchor.proposal};
            mixtureScreenPool.addAll(anchor.onPolicySamples);
            referenceWeights = new double[]{1.0};
        } else {
            ConditionalProposalSnapshot preCalibrationAnchor =
                    anchor.proposal;
            ConditionalProposalSnapshot calibratedAnchor =
                    calibrateConditionalEtaV4AnchorGauge(
                            preCalibrationAnchor, allOaisCalibration,
                            "all-oais-calibration-anchor");
            if (calibratedAnchor == null) {
                failCertificate("conditional eta-v4: OAIS calibration could not calibrate the original anchor gauge",
                        startTime);
                return;
            }
            double additionalShiftDown = RT
                    * (calibratedAnchor.logZ
                    - preCalibrationAnchor.logZ);
            double calibratedAdaptationMax = Math.max(
                    anchor.adaptationMaxOverCorrectionKcal
                            - additionalShiftDown,
                    maxOverCorrectionKcal(
                            allOaisCalibration,
                            calibratedAnchor.emat));
            anchor = new ConditionalEtaV4Anchor(
                    calibratedAnchor, anchor.sourceProposalId,
                    anchor.onPolicySamples, calibratedAdaptationMax,
                    anchor.chiSquareProxy, anchor.screeningEss,
                    anchor.epsilonAtMax, anchor.predictedFinalSamples,
                    anchor.candidateDpSweeps);
            if (!auditConditionalSupport(
                    allOaisCalibration, anchor.proposal.emat,
                    etaV4ProofBoundKcal,
                    "eta-v4-all-oais-calibration-anchor")) {
                failCertificate("conditional eta-v4: calibrated anchor failed the OAIS support audit",
                        startTime);
                return;
            }
            List<CCDResult> componentGaugeWitnesses =
                    new ArrayList<>(adaptationPool);
            componentGaugeWitnesses.addAll(allOaisCalibration);
            ConditionalProposalSnapshot calibratedTail =
                    calibrateConditionalEtaV4AnchorGauge(
                            currentOaisProposal,
                            componentGaugeWitnesses,
                            "final-oais-component");
            if (calibratedTail == null
                    || !auditConditionalSupport(
                    componentGaugeWitnesses, calibratedTail.emat,
                    etaV4ProofBoundKcal,
                    "eta-v4-final-oais-component-anchor-candidate")) {
                failCertificate("conditional eta-v4: final OAIS component could not be calibrated as an anchor candidate",
                        startTime);
                return;
            }
            components = new ConditionalProposalSnapshot[]{
                    anchor.proposal, calibratedTail};
            mixtureScreenPool.addAll(originalAnchorScreening);
            mixtureScreenPool.addAll(finalOaisCalibration);
            referenceWeights = new double[]{
                    (double) originalAnchorScreening.size()
                            / mixtureScreenPool.size(),
                    (double) finalOaisCalibration.size()
                            / mixtureScreenPool.size()
            };
        }
        System.out.println("[PACK*-conditional-eta-v4] OAIS summary:"
                + " acceptedRounds=" + acceptedOaisRounds
                + ", maxRounds=" + etaV4OaisMaxRounds
                + ", freshOnPolicyCCD="
                + allOaisCalibration.size()
                + ", cumulativeHeldOutImprovement="
                + String.format(Locale.ROOT, "%.9f",
                cumulativeHeldOutImprovement)
                + ", finalComponents=" + components.length);

        ConditionalEtaV4MixtureDesign design =
                selectConditionalEtaV4Mixture(
                        components, mixtureScreenPool,
                        referenceWeights, candidateDpSweeps);
        if (design == null) {
            failCertificate("conditional eta-v4: no valid anchor/second-moment mixture survived adaptation",
                    startTime);
            return;
        }

        List<CCDResult> mixtureCalibration =
                sampleConditionalEtaV4Mixture(
                        design, etaV4CalibrationSamples,
                        "conditional-v4-mixture-calibration");
        adaptationCCDCalls += mixtureCalibration.size();
        totalCCDCalls = trainCCDCalls + adaptationCCDCalls;
        if (mixtureCalibration.isEmpty()
                || !validateObservedLowerBound(
                mixtureCalibration,
                "conditional-v4-mixture-calibration")) {
            failCertificate("conditional eta-v4: invalid fresh IID mixture calibration batch",
                    startTime);
            return;
        }
        ConditionalProposalSnapshot[] calibratedComponents =
                design.components.clone();
        if (design.componentCover) {
            for (int k = 0; k < calibratedComponents.length; k++) {
                if (design.weights[k] == 0.0) continue;
                calibratedComponents[k] =
                        calibrateConditionalEtaV4AnchorGauge(
                                calibratedComponents[k],
                                mixtureCalibration,
                                "mixture-calibration-component-" + k);
                if (calibratedComponents[k] == null) {
                    failCertificate("conditional eta-v4: fresh IID mixture calibration could not calibrate a cover component gauge",
                            startTime);
                    return;
                }
            }
        } else {
            calibratedComponents[design.anchorIndex] =
                    calibrateConditionalEtaV4AnchorGauge(
                            calibratedComponents[design.anchorIndex],
                            mixtureCalibration,
                            "mixture-calibration");
            if (calibratedComponents[design.anchorIndex] == null) {
                failCertificate("conditional eta-v4: fresh IID mixture calibration could not calibrate the anchor gauge",
                        startTime);
                return;
            }
        }
        double[] calibratedLogZ =
                new double[calibratedComponents.length];
        double[] calibratedBounds =
                new double[calibratedComponents.length];
        for (int k = 0; k < calibratedComponents.length; k++) {
            calibratedLogZ[k] = calibratedComponents[k].logZ;
            calibratedBounds[k] = design.anchorBoundKcal;
        }
        double calibratedAnchorLogRange =
                design.componentCover
                ? PackStarEtaV4Importance.componentCoverLogRange(
                calibratedLogZ, calibratedBounds, RT,
                design.weights)
                : PackStarEtaV4Importance.anchorLogRange(
                        calibratedComponents[design.anchorIndex].logZ,
                        design.anchorBoundKcal, RT,
                        design.weights[design.anchorIndex]);
        if (!Double.isFinite(calibratedAnchorLogRange)) {
            failCertificate("conditional eta-v4: calibrated certificate has a non-finite mixture range",
                    startTime);
            return;
        }
        double calibratedEffectiveBound = RT
                * (calibratedAnchorLogRange
                - calibratedComponents[design.anchorIndex].logZ);
        if (!Double.isFinite(calibratedEffectiveBound)
                || !(calibratedEffectiveBound > 0.0)) {
            failCertificate("conditional eta-v4: calibrated certificate has an invalid effective range",
                    startTime);
            return;
        }
        design = new ConditionalEtaV4MixtureDesign(
                calibratedComponents, design.weights,
                design.anchorIndex, design.anchorBoundKcal,
                calibratedEffectiveBound,
                calibratedAnchorLogRange, design.componentCover,
                design.adaptation, design.clipLogCap,
                design.sizing, design.candidateDpSweeps);
        if (!auditConditionalEtaV4Support(
                mixtureCalibration, design,
                "mixture-calibration")) {
            failCertificate("conditional eta-v4: calibrated IID mixture failed its support or range audit",
                    startTime);
            return;
        }
        ConditionalAdaptationStats calibratedStats =
                buildConditionalEtaV4Stats(
                        mixtureCalibration, design.components,
                        design.weights, null,
                        design.anchorIndex,
                        design.effectiveBoundKcal);
        double calibratedClip = calibratedStats == null
                ? Double.NaN : selectConditionalClipLogCap(
                calibratedStats, design.effectiveBoundKcal);
        ConditionalSizing calibratedSizing = calibratedStats == null
                || !Double.isFinite(calibratedClip) ? null
                : computeConditionalSizing(
                calibratedStats, design.effectiveBoundKcal,
                calibratedClip);
        if (calibratedSizing == null
                || calibratedSizing.epsilonAtMaxSamples
                > targetEpsilon + 1.0e-12
                || calibratedSizing.finalSamples < 2) {
            failCertificate("conditional eta-v4: calibrated IID mixture is predicted unreachable at the final sample cap",
                    startTime);
            return;
        }
        design = new ConditionalEtaV4MixtureDesign(
                design.components, design.weights,
                design.anchorIndex, design.anchorBoundKcal,
                design.effectiveBoundKcal,
                design.anchorLogRange, design.componentCover,
                calibratedStats, calibratedClip,
                calibratedSizing, design.candidateDpSweeps);
        selectedResidualBoundKcal = design.anchorBoundKcal;
        clipLogCap = calibratedClip;
        int nStar = calibratedSizing.finalSamples;
        System.out.println("[PACK*-conditional-eta-v4] frozen design:"
                + " components=" + design.components.length
                + ", anchorIndex=" + design.anchorIndex
                + ", anchorSource="
                + design.components[design.anchorIndex].sourceProposalId
                + ", certificateMode="
                + (design.componentCover
                ? "component-cover" : "single-anchor")
                + ", weights=" + Arrays.toString(design.weights)
                + ", anchorB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                design.anchorBoundKcal)
                + ", effectiveB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                design.effectiveBoundKcal)
                + ", adaptationCCD=" + adaptationCCDCalls
                + ", candidateDpSweeps="
                + design.candidateDpSweeps
                + ", calibrationESS="
                + String.format(Locale.ROOT, "%.6f",
                calibratedStats.effectiveSampleSize)
                + ", clipLogCap="
                + String.format(Locale.ROOT, "%.6f", clipLogCap)
                + ", monitorN=" + monitorSamples
                + ", finalN=" + nStar
                + ", iidComponentLabels=true"
                + ", balanceDenominator=true"
                + ", q_mFallback=false, deterministicFallback=false");

        List<CCDResult> monitorCCD = sampleConditionalEtaV4Mixture(
                design, monitorSamples,
                "conditional-v4-monitor");
        monitorCCDCalls = monitorCCD.size();
        totalCCDCalls = trainCCDCalls + adaptationCCDCalls
                + monitorCCDCalls;
        if (!auditConditionalEtaV4Support(
                monitorCCD, design, "monitor")) {
            failCertificate("conditional eta-v4: independent IID mixture monitor observed a support or range violation",
                    startTime);
            return;
        }

        List<CCDResult> finalCCD = sampleConditionalEtaV4Mixture(
                design, nStar, "conditional-v4-final");
        finalCCDCalls = finalCCD.size();
        totalCCDCalls = trainCCDCalls + adaptationCCDCalls
                + monitorCCDCalls + finalCCDCalls;
        if (!auditConditionalEtaV4Support(
                finalCCD, design, "final")) {
            failCertificate("conditional eta-v4: fresh IID mixture final set observed a support or range violation",
                    startTime);
            return;
        }
        computeConditionalEtaV4Bound(finalCCD, design);
        if (!certificateValid
                && "not computed".equals(certificateFailureReason)) {
            certificateFailureReason =
                    "conditional eta-v4 mixture interval failed validation";
        }
        printFinalSummary(startTime);
    }

    /**
     * Eta-v5 Phase A deliberately ends after an exhaustive, no-CCD threshold
     * mass oracle. Even a passing gate remains an Aborted PACK* pfunc because
     * no fresh two-stratum PAC estimate has been constructed yet.
     */
    private void runEtaV5ThresholdOracle(
            long startTime,
            ConditionalProposalSnapshot proposal) {
        File output = new File(etaV5OracleOutputTsv);
        File parent = output.getParentFile();
        if (parent == null || !parent.isDirectory()) {
            failCertificate(
                    "eta-v5 Phase-A oracle output parent does not exist: "
                            + (parent == null ? "null" : parent),
                    startTime);
            return;
        }
        if (output.exists()) {
            failCertificate(
                    "eta-v5 Phase-A oracle refuses to overwrite output: "
                            + output,
                    startTime);
            return;
        }

        long oracleStart = System.currentTimeMillis();
        final double canonicalShiftKcal;
        final PackStarThresholdTail.OracleResult result;
        try {
            canonicalShiftKcal =
                    PackStarThresholdTail.canonicalShiftKcal(
                            RT, logZMinDet, proposal.logZ);
            result = PackStarThresholdTail.enumerate(
                    rcs,
                    branchMinimizingEmat,
                    proposal.emat,
                    interactionGraph,
                    RT,
                    canonicalShiftKcal,
                    etaV5ThresholdsKcal,
                    etaV5MaxEnumeratedStates,
                    etaV5RebaseInterval,
                    (visited, total) -> System.out.println(
                            "[PACK*-eta-v5] Phase-A enumeration progress: "
                                    + visited + "/" + total));
        } catch (RuntimeException ex) {
            writeEtaV5OracleError(output, proposal, ex);
            failCertificate(
                    "eta-v5 Phase-A enumeration failed: "
                            + ex.getClass().getSimpleName() + ": "
                            + ex.getMessage(),
                    startTime);
            return;
        }
        long oracleElapsedMs = System.currentTimeMillis() - oracleStart;

        double minimizingLogZDelta = Math.abs(
                result.logZMinEnumerated - logZMinDet);
        double proposalLogZDelta = Math.abs(
                result.logZProposalEnumerated - proposal.logZ);
        double canonicalLogZDelta = Math.abs(
                result.logZProposalCanonicalEnumerated
                        - result.logZMinEnumerated);
        boolean stateAccountingPass =
                result.visitedStates == result.cartesianStates
                        && result.finiteMinimizingStates
                        + result.zeroMassStates
                        == result.cartesianStates;
        boolean monotonePass = result.hasMonotoneTail(
                etaV5LogZTolerance);
        boolean minimizingNormalizerPass =
                Double.isFinite(minimizingLogZDelta)
                        && minimizingLogZDelta
                        <= etaV5LogZTolerance;
        boolean proposalNormalizerPass =
                Double.isFinite(proposalLogZDelta)
                        && proposalLogZDelta
                        <= etaV5LogZTolerance;
        boolean canonicalNormalizerPass =
                Double.isFinite(canonicalLogZDelta)
                        && canonicalLogZDelta
                        <= etaV5LogZTolerance;
        boolean endpointPass = etaV5EndpointChecksPass(result);
        boolean gatePassed = stateAccountingPass
                && monotonePass
                && minimizingNormalizerPass
                && proposalNormalizerPass
                && canonicalNormalizerPass
                && endpointPass;

        String gateReason = gatePassed
                ? "pass"
                : "stateAccounting=" + stateAccountingPass
                + ",monotone=" + monotonePass
                + ",minimizingNormalizer="
                + minimizingNormalizerPass
                + ",proposalNormalizer="
                + proposalNormalizerPass
                + ",canonicalNormalizer="
                + canonicalNormalizerPass
                + ",endpoints=" + endpointPass;
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tproposalSource\tintegrityGatePassed"
                    + "\tintegrityGateReason\tthresholdGauge"
                    + "\tcanonicalThresholdKcal\trawThresholdKcal"
                    + "\tcanonicalCap\ttailStateCount"
                    + "\tlogTailMass\tlogTailFraction\ttailMassFraction"
                    + "\tcartesianStates\tfiniteMinimizingStates"
                    + "\tzeroMassStates\tinfiniteEtaStates"
                    + "\tlogZMinEnumerated\tlogZMinDp\tlogZMinAbsDelta"
                    + "\tlogZProposalRawEnumerated\tlogZProposalRawDp"
                    + "\tlogZProposalRawAbsDelta\tcanonicalShiftKcal"
                    + "\tlogZProposalCanonicalEnumerated"
                    + "\tlogZCanonicalAbsDelta"
                    + "\trawEtaMinKcal\trawEtaMaxKcal"
                    + "\tcanonicalEtaMinKcal\tcanonicalEtaMaxKcal"
                    + "\tmaxMinimizingRebaseDriftKcal"
                    + "\tmaxProposalRebaseDriftKcal\telapsedMs");
            for (int i = 0; i < result.thresholdsKcal.length; i++) {
                double logFraction = result.logTailFraction(i);
                double fraction = logFraction == Double.NEGATIVE_INFINITY
                        ? 0.0 : Math.exp(logFraction);
                writer.printf(Locale.ROOT,
                        "packstar-eta-v5-threshold-oracle-v2\t%s\t%s\t%s"
                                + "\tcanonical-normalizer"
                                + "\t%.17g\t%.17g\t%.17g\t%d"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%d\t%d\t%d\t%d"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%d%n",
                        sanitizeTsv(proposal.sourceProposalId),
                        Boolean.toString(gatePassed),
                        gateReason,
                        result.thresholdsKcal[i],
                        result.thresholdsKcal[i]
                                + result.canonicalShiftKcal,
                        Math.exp(result.thresholdsKcal[i] / RT),
                        result.tailStateCounts[i],
                        result.logTailMass[i],
                        logFraction,
                        fraction,
                        result.cartesianStates,
                        result.finiteMinimizingStates,
                        result.zeroMassStates,
                        result.infiniteEtaStates,
                        result.logZMinEnumerated,
                        logZMinDet,
                        minimizingLogZDelta,
                        result.logZProposalEnumerated,
                        proposal.logZ,
                        proposalLogZDelta,
                        result.canonicalShiftKcal,
                        result.logZProposalCanonicalEnumerated,
                        canonicalLogZDelta,
                        result.rawEtaMinKcal,
                        result.rawEtaMaxKcal,
                        result.canonicalEtaMinKcal,
                        result.canonicalEtaMaxKcal,
                        result.maxMinimizingRebaseDriftKcal,
                        result.maxProposalRebaseDriftKcal,
                        oracleElapsedMs);
            }
        } catch (Exception ex) {
            failCertificate(
                    "eta-v5 Phase-A oracle could not write TSV: "
                            + ex.getMessage(),
                    startTime);
            return;
        }

        System.out.println("[PACK*-eta-v5] Phase-A oracle integrity gate: "
                + (gatePassed ? "PASS" : "FAIL")
                + ", states=" + result.cartesianStates
                + ", finiteMassStates="
                + result.finiteMinimizingStates
                + ", zeroMassStates=" + result.zeroMassStates
                + ", canonicalShiftKcal="
                + String.format(Locale.ROOT, "%.9f",
                result.canonicalShiftKcal)
                + ", rawEtaRangeKcal=["
                + String.format(Locale.ROOT, "%.9f",
                result.rawEtaMinKcal)
                + ","
                + String.format(Locale.ROOT, "%.9f",
                result.rawEtaMaxKcal) + "]"
                + ", canonicalEtaRangeKcal=["
                + String.format(Locale.ROOT, "%.9f",
                result.canonicalEtaMinKcal)
                + ","
                + String.format(Locale.ROOT, "%.9f",
                result.canonicalEtaMaxKcal) + "]"
                + ", logZMinAbsDelta="
                + String.format(Locale.ROOT, "%.3e",
                minimizingLogZDelta)
                + ", logZProposalAbsDelta="
                + String.format(Locale.ROOT, "%.3e",
                proposalLogZDelta)
                + ", logZCanonicalAbsDelta="
                + String.format(Locale.ROOT, "%.3e",
                canonicalLogZDelta)
                + ", elapsedMs=" + oracleElapsedMs
                + ", output=" + output);
        if (gatePassed && etaV5MechanismPilot) {
            runEtaV5MechanismPilot(
                    startTime, proposal, result);
            return;
        }
        failCertificate(
                gatePassed
                        ? "eta-v5 Phase-A exact canonical-threshold oracle"
                        + " passed its integrity gate;"
                        + " diagnostic-only stop before any tail CCD"
                        : "eta-v5 Phase-A exact canonical-threshold oracle"
                        + " integrity gate failed: "
                        + gateReason,
                startTime);
    }

    /**
     * One frozen diagnostic batch from p_m.  This does not estimate a PACK*
     * pfunc or issue a certificate; it only tests whether the Phase-A safe
     * superset is dominated by false positives after direct CCD reveals g.
     */
    private void runEtaV5MechanismPilot(
            long startTime,
            ConditionalProposalSnapshot proposal,
            PackStarThresholdTail.OracleResult oracle) {
        File assignmentsOutput = new File(
                etaV5MechanismAssignmentsTsv);
        File summaryOutput = new File(
                etaV5MechanismSummaryTsv);
        for (File output : new File[]{
                assignmentsOutput, summaryOutput}) {
            File parent = output.getParentFile();
            if (parent == null || !parent.isDirectory()) {
                failCertificate(
                        "eta-v5 mechanism output parent does not exist: "
                                + (parent == null ? "null" : parent),
                        startTime);
                return;
            }
            if (output.exists()) {
                failCertificate(
                        "eta-v5 mechanism pilot refuses to overwrite output: "
                                + output,
                        startTime);
                return;
            }
        }

        final double minimizingPilotLogZ;
        final List<int[]> conformations;
        try {
            CorrectedDPResult minimizingDP = recomputeDP(
                    branchMinimizingEmat, zeroEtaCorrections());
            minimizingPilotLogZ = minimizingDP.logZCorrected;
            if (!Double.isFinite(minimizingPilotLogZ)
                    || Math.abs(minimizingPilotLogZ - logZMinDet)
                    > etaV5LogZTolerance) {
                failCertificate(
                        "eta-v5 mechanism p_m DP normalizer mismatch: pilot="
                                + minimizingPilotLogZ + " initial="
                                + logZMinDet,
                        startTime);
                return;
            }
            conformations = sampleConformationsFromDP(
                    etaV5MechanismSamples,
                    stageRandom("eta-v5-pm-mechanism-pilot"));
        } catch (RuntimeException ex) {
            failCertificate(
                    "eta-v5 mechanism p_m sampling failed: "
                            + ex.getClass().getSimpleName() + ": "
                            + ex.getMessage(),
                    startTime);
            return;
        }
        if (conformations.size() != etaV5MechanismSamples) {
            failCertificate(
                    "eta-v5 mechanism p_m sampler returned "
                            + conformations.size() + " of "
                            + etaV5MechanismSamples + " assignments",
                    startTime);
            return;
        }

        System.out.println("[PACK*-eta-v5] mechanism pilot: drawing "
                + etaV5MechanismSamples
                + " IID p_m assignments after normalizer audit"
                + ", logZMin="
                + String.format(Locale.ROOT, "%.12f",
                minimizingPilotLogZ)
                + ", proposalSource="
                + proposal.sourceProposalId);
        List<CCDResult> ccd = runParallelCCD(
                conformations, proposal.emat, proposal.logZ,
                proposal.sourceProposalId);
        totalCCDCalls += ccd.size();
        if (ccd.size() != etaV5MechanismSamples
                || !validateObservedLowerBound(
                ccd, "eta-v5-pm-mechanism-pilot")) {
            failCertificate(
                    "eta-v5 mechanism pilot produced an incomplete CCD"
                            + " batch or lower-bound violation",
                    startTime);
            return;
        }

        List<PackStarEtaV5MechanismPilot.Sample> samples =
                new ArrayList<>(ccd.size());
        for (CCDResult value : ccd) {
            if (!Double.isFinite(value.sourceProposalEnergy)) {
                failCertificate(
                        "eta-v5 mechanism sample lacks frozen proposal"
                                + " energy provenance",
                        startTime);
                return;
            }
            samples.add(new PackStarEtaV5MechanismPilot.Sample(
                    value.conf, value.eMin,
                    value.sourceProposalEnergy, value.eTrue));
        }

        String pilotStream = randomStreamIdentity
                + "|unconditional-v1|eta-v5-pm-mechanism-pilot";
        final PackStarEtaV5MechanismPilot.Result result;
        try {
            result = PackStarEtaV5MechanismPilot.analyze(
                    samples, oracle.canonicalShiftKcal,
                    etaV5ThresholdsKcal, RT);
            PackStarEtaV5MechanismPilot.writeAssignmentTsv(
                    assignmentsOutput, result,
                    proposal.sourceProposalId, pilotStream);
            PackStarEtaV5MechanismPilot.writeSummaryTsv(
                    summaryOutput, result,
                    proposal.sourceProposalId, pilotStream);
        } catch (Exception ex) {
            failCertificate(
                    "eta-v5 mechanism pilot analysis/output failed: "
                            + ex.getClass().getSimpleName() + ": "
                            + ex.getMessage(),
                    startTime);
            return;
        }

        int primaryIndex = -1;
        for (int i = 0; i < result.thresholdsKcal.length; i++) {
            if (Double.compare(result.thresholdsKcal[i], 0.0) == 0) {
                primaryIndex = i;
                break;
            }
        }
        if (primaryIndex < 0) {
            failCertificate(
                    "eta-v5 mechanism pilot is missing frozen primary"
                            + " threshold tau=0",
                    startTime);
            return;
        }
        System.out.println("[PACK*-eta-v5] mechanism primary endpoint:"
                + " tau=0, actualExceedCount="
                + result.exceedCounts[primaryIndex] + "/"
                + result.samples.size()
                + ", meanH="
                + String.format(Locale.ROOT, "%.12g",
                result.meanH[primaryIndex])
                + ", maxH="
                + String.format(Locale.ROOT, "%.12g",
                result.maxH[primaryIndex])
                + ", pearson(rawEta,g)="
                + String.format(Locale.ROOT, "%.9f",
                result.pearsonRawEtaGap)
                + ", spearman(rawEta,g)="
                + String.format(Locale.ROOT, "%.9f",
                result.spearmanRawEtaGap)
                + ", assignments=" + assignmentsOutput
                + ", summary=" + summaryOutput);
        failCertificate(
                "eta-v5 200-IID-p_m direct-CCD mechanism pilot completed;"
                        + " diagnostic-only stop with no certificate",
                startTime);
    }

    private boolean etaV5EndpointChecksPass(
            PackStarThresholdTail.OracleResult result) {
        for (int i = 0; i < result.thresholdsKcal.length; i++) {
            double threshold = result.thresholdsKcal[i];
            if (threshold < result.canonicalEtaMinKcal) {
                double deltaFromFullMass = Math.abs(
                        result.logTailMass[i]
                                - result.logZMinEnumerated);
                if (!Double.isFinite(deltaFromFullMass)
                        || deltaFromFullMass
                        > etaV5LogZTolerance) {
                    return false;
                }
            }
            if (Double.isFinite(result.canonicalEtaMaxKcal)
                    && threshold >= result.canonicalEtaMaxKcal
                    && (result.tailStateCounts[i] != 0L
                    || result.logTailMass[i]
                    != Double.NEGATIVE_INFINITY)) {
                return false;
            }
        }
        return true;
    }

    private void writeEtaV5OracleError(
            File output,
            ConditionalProposalSnapshot proposal,
            RuntimeException error) {
        try (PrintWriter writer = new PrintWriter(
                output, StandardCharsets.UTF_8.name())) {
            writer.println("schema\tproposalSource"
                    + "\tintegrityGatePassed\terrorClass\terrorMessage");
            writer.println("packstar-eta-v5-threshold-oracle-error-v2\t"
                    + sanitizeTsv(proposal.sourceProposalId)
                    + "\tfalse\t"
                    + sanitizeTsv(error.getClass().getName())
                    + "\t" + sanitizeTsv(error.getMessage()));
        } catch (Exception ignored) {
            System.out.println("[PACK*-eta-v5] WARNING: could not write"
                    + " the Phase-A oracle error TSV: " + ignored);
        }
    }

    private static String sanitizeTsv(String value) {
        if (value == null) return "";
        return value.replace('\t', ' ')
                .replace('\n', ' ')
                .replace('\r', ' ');
    }

    private static boolean isEtaV5XtmpOutputPath(String value) {
        return isXtmpOutputPath(value);
    }

    private static boolean isXtmpOutputPath(String value) {
        return value != null && !value.isEmpty()
                && new File(value).isAbsolute()
                && value.startsWith("/usr/xtmp/lz280/");
    }

    private double selectConditionalClipLogCap(
            List<CCDResult> adaptation,
            EnergyMatrix correctedEmat,
            double boundKcal) {
        if (adaptation.isEmpty() || !Double.isFinite(boundKcal)
                || !(boundKcal > 0.0)) return Double.NaN;
        if (!getConfigBoolean("packstar.pac.clip", DEFAULT_CLIP)) {
            return boundKcal / RT;
        }
        double[] logW = new double[adaptation.size()];
        for (int i = 0; i < adaptation.size(); i++) {
            CCDResult sample = adaptation.get(i);
            double eEta = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat);
            logW[i] = (eEta - sample.eTrue) / RT;
            if (!Double.isFinite(logW[i])
                    || logW[i] > boundKcal / RT + 1.0e-12) {
                return Double.NaN;
            }
        }
        double clipQuantile = getConfigDouble(
                "packstar.pac.clipQuantile",
                DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE);
        // C >= 1 avoids manufacturing an excess term for ordinary weights
        // below one.  The upper endpoint remains the proof-side cap exp(B/RT).
        return Math.max(0.0, Math.min(boundKcal / RT,
                quantile(logW, clipQuantile)));
    }

    private int computeRequiredConditionalEstSamples(
            List<CCDResult> adaptation,
            EnergyMatrix correctedEmat,
            double boundKcal,
            double logCap) {
        int n = adaptation.size();
        if (n <= 1 || !Double.isFinite(boundKcal)
                || !Double.isFinite(logCap)) return maxEstSamples;
        double uRange = conditionalExcessRange(boundKcal, RT, logCap);
        if (!Double.isFinite(uRange)) return maxEstSamples;

        double sum = 0.0;
        double sum2 = 0.0;
        double sumU = 0.0;
        double sumU2 = 0.0;
        int nClip = 0;
        for (CCDResult sample : adaptation) {
            double eEta = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat);
            double logW = (eEta - sample.eTrue) / RT;
            if (!Double.isFinite(logW)
                    || logW > boundKcal / RT + 1.0e-12) {
                return maxEstSamples;
            }
            double clipped = conditionalClippedNormalized(logW, logCap);
            double excess = conditionalExcessNormalized(logW, logCap);
            sum += clipped;
            sum2 += clipped * clipped;
            sumU += excess;
            sumU2 += excess * excess;
            if (logW > logCap) nClip++;
        }
        double mean = sum / n;
        double meanU = sumU / n;
        double var = Math.max(0.0,
                (sum2 - n * mean * mean) / (n - 1)) * nstarInflate;
        double varU = Math.max(0.0,
                (sumU2 - n * meanU * meanU) / (n - 1)) * nstarInflate;
        double pHat = (double) nClip / n;
        double deltaPer = perEventDelta(true);
        int nStar = solveClipSizing(mean, var, meanU, varU,
                uRange, pHat, deltaPer);
        System.out.println("[PACK*-conditional] N* sizing: pilot=" + n
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f", boundKcal)
                + ", meanWc="
                + String.format(Locale.ROOT, "%.6f", mean)
                + ", varWcInflated="
                + String.format(Locale.ROOT, "%.6f", var)
                + ", meanExcess="
                + String.format(Locale.ROOT, "%.6f", meanU)
                + ", uRange="
                + String.format(Locale.ROOT, "%.6f", uRange)
                + ", pClip="
                + String.format(Locale.ROOT, "%.4f", pHat)
                + ", finalN=" + nStar
                + ", cap=" + maxEstSamples);
        return nStar;
    }

    private ConditionalAdaptationStats buildConditionalAdaptationStats(
            List<CCDResult> adaptation,
            EnergyMatrix targetEmat,
            double targetLogZ,
            double boundKcal) {
        if (adaptation == null || adaptation.size() < 2
                || targetEmat == null || !Double.isFinite(targetLogZ)
                || !Double.isFinite(boundKcal) || !(boundKcal > 0.0)) {
            return null;
        }
        double[] logResidualWeights = new double[adaptation.size()];
        double[] logProposalWeights = new double[adaptation.size()];
        Set<String> sourceProposalIds = new LinkedHashSet<>();
        for (int i = 0; i < adaptation.size(); i++) {
            CCDResult sample = adaptation.get(i);
            if (!Double.isFinite(sample.eTrue)
                    || !Double.isFinite(sample.sourceProposalEnergy)
                    || !Double.isFinite(sample.sourceProposalLogZ)
                    || sample.sourceProposalId == null
                    || sample.sourceProposalId.isEmpty()) {
                return null;
            }
            double targetEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, targetEmat);
            if (!Double.isFinite(targetEnergy)) return null;
            double logResidualWeight = (targetEnergy - sample.eTrue) / RT;
            if (!Double.isFinite(logResidualWeight)
                    || logResidualWeight > boundKcal / RT + 1.0e-12) {
                return null;
            }
            double logProposalWeight = conditionalImportanceLogRatio(
                    sample.sourceProposalEnergy, targetEnergy,
                    sample.sourceProposalLogZ, targetLogZ, RT);
            if (!Double.isFinite(logProposalWeight)) return null;
            logResidualWeights[i] = logResidualWeight;
            logProposalWeights[i] = logProposalWeight;
            sourceProposalIds.add(sample.sourceProposalId);
        }
        ConditionalImportanceWeights normalized =
                normalizeConditionalImportanceWeights(logProposalWeights);
        if (normalized == null) return null;
        return new ConditionalAdaptationStats(
                logResidualWeights, normalized.normalized,
                normalized.effectiveSampleSize,
                normalized.minLogWeight, normalized.maxLogWeight,
                sourceProposalIds.size());
    }

    private double selectConditionalClipLogCap(
            ConditionalAdaptationStats adaptation,
            double boundKcal) {
        if (adaptation == null || !Double.isFinite(boundKcal)
                || !(boundKcal > 0.0)) return Double.NaN;
        if (!getConfigBoolean("packstar.pac.clip", DEFAULT_CLIP)) {
            return boundKcal / RT;
        }
        double clipQuantile = getConfigDouble(
                "packstar.pac.clipQuantile",
                DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE);
        double weightedQuantile = conditionalWeightedQuantile(
                adaptation.logResidualWeights,
                adaptation.normalizedProposalWeights,
                clipQuantile);
        if (!Double.isFinite(weightedQuantile)) return Double.NaN;
        return Math.max(0.0, Math.min(boundKcal / RT,
                weightedQuantile));
    }

    private ConditionalSizing computeConditionalSizing(
            ConditionalAdaptationStats adaptation,
            double boundKcal,
            double logCap) {
        if (adaptation == null
                || adaptation.logResidualWeights.length < 2
                || !Double.isFinite(boundKcal)
                || !Double.isFinite(logCap)) return null;
        double uRange = conditionalExcessRange(boundKcal, RT, logCap);
        if (!Double.isFinite(uRange)) return null;

        int n = adaptation.logResidualWeights.length;
        double[] clipped = new double[n];
        double[] excess = new double[n];
        double[] clippedIndicator = new double[n];
        for (int i = 0; i < n; i++) {
            double logWeight = adaptation.logResidualWeights[i];
            if (!Double.isFinite(logWeight)
                    || logWeight > boundKcal / RT + 1.0e-12) return null;
            clipped[i] = conditionalClippedNormalized(logWeight, logCap);
            excess[i] = conditionalExcessNormalized(logWeight, logCap);
            clippedIndicator[i] = logWeight > logCap ? 1.0 : 0.0;
        }
        double[] proposalWeights = adaptation.normalizedProposalWeights;
        double mean = conditionalWeightedMean(clipped, proposalWeights);
        double meanU = conditionalWeightedMean(excess, proposalWeights);
        double variance = conditionalWeightedSampleVariance(
                clipped, proposalWeights, mean);
        double varianceU = conditionalWeightedSampleVariance(
                excess, proposalWeights, meanU);
        double pHat = conditionalWeightedMean(
                clippedIndicator, proposalWeights);
        if (!Double.isFinite(mean) || !Double.isFinite(meanU)
                || !Double.isFinite(variance)
                || !Double.isFinite(varianceU)
                || !Double.isFinite(pHat)) return null;
        variance *= nstarInflate;
        varianceU *= nstarInflate;
        double deltaPer = perEventDelta(true);
        int nStar = solveClipSizing(mean, variance, meanU, varianceU,
                uRange, pHat, deltaPer);
        double epsilonAtMax = predictClipEps(
                maxEstSamples, mean, variance, meanU, varianceU,
                uRange, pHat, deltaPer);
        if (!Double.isFinite(epsilonAtMax)) return null;
        return new ConditionalSizing(
                mean, variance, meanU, varianceU,
                uRange, pHat, nStar, epsilonAtMax);
    }

    private int computeRequiredConditionalEstSamples(
            ConditionalAdaptationStats adaptation,
            double boundKcal,
            double logCap) {
        ConditionalSizing sizing = computeConditionalSizing(
                adaptation, boundKcal, logCap);
        if (sizing == null) return maxEstSamples;
        System.out.println("[PACK*-conditional-dp-rich] N* sizing: adaptation="
                + adaptation.logResidualWeights.length
                + ", sourceProposals=" + adaptation.sourceProposalCount
                + ", ESS="
                + String.format(Locale.ROOT, "%.3f",
                adaptation.effectiveSampleSize)
                + ", ESSFraction="
                + String.format(Locale.ROOT, "%.4f",
                adaptation.effectiveSampleSize
                        / adaptation.logResidualWeights.length)
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f", boundKcal)
                + ", meanWc="
                + String.format(Locale.ROOT, "%.6f", sizing.meanClipped)
                + ", varWcInflated="
                + String.format(Locale.ROOT, "%.6f",
                sizing.varianceClipped)
                + ", meanExcess="
                + String.format(Locale.ROOT, "%.6f", sizing.meanExcess)
                + ", uRange="
                + String.format(Locale.ROOT, "%.6f", sizing.excessRange)
                + ", pClip="
                + String.format(Locale.ROOT, "%.4f",
                sizing.clippedProbability)
                + ", epsilonAtMax="
                + String.format(Locale.ROOT, "%.6f",
                sizing.epsilonAtMaxSamples)
                + ", finalN=" + sizing.finalSamples
                + ", cap=" + maxEstSamples);
        return sizing.finalSamples;
    }

    private boolean auditConditionalSupport(
            List<CCDResult> samples,
            EnergyMatrix correctedEmat,
            double boundKcal,
            String stage) {
        if (samples.isEmpty() || !Double.isFinite(boundKcal)
                || !(boundKcal > 0.0)) return false;
        int lowerViolations = 0;
        int residualViolations = 0;
        double minGap = Double.POSITIVE_INFINITY;
        double maxOverCorrection = Double.NEGATIVE_INFINITY;
        for (CCDResult sample : samples) {
            if (!Double.isFinite(sample.eTrue)
                    || !Double.isFinite(sample.eMin)) return false;
            double eEta = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat);
            if (!Double.isFinite(eEta)) return false;
            double gap = sample.eTrue - sample.eMin;
            double overCorrection = eEta - sample.eTrue;
            minGap = Math.min(minGap, gap);
            maxOverCorrection = Math.max(maxOverCorrection,
                    overCorrection);
            if (gap < 0.0) lowerViolations++;
            if (overCorrection > boundKcal) residualViolations++;
        }
        System.out.println("[PACK*-conditional-audit] " + stage
                + ": samples=" + samples.size()
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f", boundKcal)
                + ", maxOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f", maxOverCorrection)
                + ", residualBoundViolations=" + residualViolations
                + ", min(E_true-E_m)="
                + String.format(Locale.ROOT, "%.9f", minGap)
                + ", lowerBoundViolations=" + lowerViolations);
        maybeWritePairwiseMomentTailAudit(
                samples, correctedEmat, stage);
        return lowerViolations == 0 && residualViolations == 0;
    }

    private void computeConditionalRepairedBound(
            List<CCDResult> finalSamples,
            EnergyMatrix correctedEmat,
            double logZCorrected,
            double boundKcal) {
        int n = finalSamples.size();
        if (n <= 1 || !Double.isFinite(logZCorrected)
                || !Double.isFinite(boundKcal)
                || !Double.isFinite(clipLogCap)) {
            failCertificate("conditional: invalid final interval inputs", -1L);
            return;
        }
        double uRange = conditionalExcessRange(
                boundKcal, RT, clipLogCap);
        if (!Double.isFinite(uRange)) {
            failCertificate("conditional: non-finite deterministic excess range", -1L);
            return;
        }

        double[] xi = new double[n];
        double[] logW = new double[n];
        double sum = 0.0;
        double sum2 = 0.0;
        double sumU = 0.0;
        double sumU2 = 0.0;
        double sumRes = 0.0;
        double sumRes2 = 0.0;
        int nClipped = 0;
        for (int i = 0; i < n; i++) {
            CCDResult sample = finalSamples.get(i);
            double eEta = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat);
            xi[i] = sample.eTrue - eEta;
            logW[i] = -xi[i] / RT;
            if (!Double.isFinite(logW[i])
                    || logW[i] > boundKcal / RT + 1.0e-12) {
                failCertificate("conditional: final support violation reached interval construction",
                        -1L);
                return;
            }
            double clipped = conditionalClippedNormalized(
                    logW[i], clipLogCap);
            double excess = conditionalExcessNormalized(
                    logW[i], clipLogCap);
            if (!Double.isFinite(clipped) || !Double.isFinite(excess)
                    || clipped < 0.0 || clipped > 1.0 + 1.0e-12
                    || excess < 0.0 || excess > uRange + 1.0e-10) {
                failCertificate("conditional: bounded final statistic is invalid",
                        -1L);
                return;
            }
            sum += clipped;
            sum2 += clipped * clipped;
            sumU += excess;
            sumU2 += excess * excess;
            sumRes += xi[i];
            sumRes2 += xi[i] * xi[i];
            if (logW[i] > clipLogCap) nClipped++;
        }

        double mean = sum / n;
        double meanU = sumU / n;
        double var = Math.max(0.0,
                (sum2 - n * mean * mean) / (n - 1));
        double varU = Math.max(0.0,
                (sumU2 - n * meanU * meanU) / (n - 1));
        double deltaPer = perEventDelta(true);
        double clippedDelta = solveBernsteinDelta(
                n, var, 1.0, deltaPer);
        double excessDelta = uRange > 0.0
                ? solveBernsteinDelta(n, varU, uRange, deltaPer)
                : 0.0;
        double pHat = (double) nClipped / n;
        double pUpper = Math.min(1.0,
                pHat + hoeffdingHalfWidth(n, deltaPer));
        double biasWorst = uRange * pUpper;
        double biasEmpirical = meanU + excessDelta;
        double biasUpper = Math.min(biasWorst, biasEmpirical);
        if (!Double.isFinite(mean) || !Double.isFinite(var)
                || !Double.isFinite(clippedDelta)
                || !Double.isFinite(biasUpper)) {
            failCertificate("conditional: non-finite final PAC statistics", -1L);
            return;
        }

        double meanLower = Math.max(0.0, mean - clippedDelta);
        double meanUpper = mean + clippedDelta + biasUpper;
        logZLowerPAC = meanLower > 0.0
                ? logZCorrected + clipLogCap + Math.log(meanLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = meanUpper > 0.0
                ? logZCorrected + clipLogCap + Math.log(meanUpper)
                : Double.NaN;
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);

        meanResidual = sumRes / n;
        stdResidual = Math.sqrt(Math.max(0.0,
                (sumRes2 - n * meanResidual * meanResidual) / (n - 1)));
        meanPsi = shiftedMeanToDouble(mean + meanU, clipLogCap);
        varPsi = var;
        cvPsi = mean > 0.0 ? Math.sqrt(var) / mean : Double.MAX_VALUE;
        emitResidualSampleTraces(finalSamples, xi, logW,
                logZCorrected, clipLogCap);
        maybeDumpResidualSamples(finalSamples, xi,
                logZCorrected, meanResidual, stdResidual);

        certificateValid = isValidCertificate(zLower, zUpper, epsilon);
        certificateFailureReason = certificateValid
                ? "" : "conditional repaired interval failed final validation";
        System.out.println("[PACK*-conditional] final certificate: selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f", boundKcal)
                + ", logCap="
                + String.format(Locale.ROOT, "%.6f", clipLogCap)
                + ", clippedDelta="
                + String.format(Locale.ROOT, "%.6f", clippedDelta)
                + ", excessMean="
                + String.format(Locale.ROOT, "%.6f", meanU)
                + ", excessRange="
                + String.format(Locale.ROOT, "%.6f", uRange)
                + ", excessUpper="
                + String.format(Locale.ROOT, "%.6f", biasUpper)
                + ", nClipped=" + nClipped + "/" + n
                + ", upperVia=conditional-sampling-only"
                + ", q_mFallback=false"
                + ", log10Lower=" + formatLog10(logZLowerPAC)
                + ", log10Upper=" + formatLog10(logZUpperPAC)
                + ", epsilon="
                + String.format(Locale.ROOT, "%.8f", epsilon)
                + ", certified=" + certificateValid);
    }

    private void computeConditionalEtaV4Bound(
            List<CCDResult> finalSamples,
            ConditionalEtaV4MixtureDesign design) {
        int n = finalSamples == null ? 0 : finalSamples.size();
        if (n <= 1 || design == null
                || !Double.isFinite(design.effectiveBoundKcal)
                || !Double.isFinite(clipLogCap)) {
            failCertificate("conditional eta-v4: invalid final interval inputs",
                    -1L);
            return;
        }
        double logZReference =
                design.components[design.anchorIndex].logZ;
        if (!Double.isFinite(logZReference)) {
            failCertificate("conditional eta-v4: non-finite anchor normalizer",
                    -1L);
            return;
        }
        double uRange = conditionalExcessRange(
                design.effectiveBoundKcal, RT, clipLogCap);
        if (!Double.isFinite(uRange)) {
            failCertificate("conditional eta-v4: non-finite deterministic excess range",
                    -1L);
            return;
        }

        double[] xi = new double[n];
        double[] logW = new double[n];
        double sum = 0.0;
        double sum2 = 0.0;
        double sumU = 0.0;
        double sumU2 = 0.0;
        double sumRes = 0.0;
        double sumRes2 = 0.0;
        int nClipped = 0;
        for (int i = 0; i < n; i++) {
            CCDResult sample = finalSamples.get(i);
            double[] logComponents =
                    conditionalEtaV4LogComponentProbabilities(
                            sample.conf, design.components);
            double logQMix = logComponents == null ? Double.NaN
                    : PackStarEtaV4Importance.logMixtureProbability(
                            logComponents, design.weights);
            logW[i] = -sample.eTrue / RT
                    - logQMix - logZReference;
            xi[i] = -RT * logW[i];
            if (!Double.isFinite(logW[i])
                    || logW[i] > design.effectiveBoundKcal / RT
                    + 1.0e-12) {
                failCertificate("conditional eta-v4: final support violation reached interval construction",
                        -1L);
                return;
            }
            double clipped = conditionalClippedNormalized(
                    logW[i], clipLogCap);
            double excess = conditionalExcessNormalized(
                    logW[i], clipLogCap);
            if (!Double.isFinite(clipped) || !Double.isFinite(excess)
                    || clipped < 0.0 || clipped > 1.0 + 1.0e-12
                    || excess < 0.0 || excess > uRange + 1.0e-10) {
                failCertificate("conditional eta-v4: bounded final statistic is invalid",
                        -1L);
                return;
            }
            sum += clipped;
            sum2 += clipped * clipped;
            sumU += excess;
            sumU2 += excess * excess;
            sumRes += xi[i];
            sumRes2 += xi[i] * xi[i];
            if (logW[i] > clipLogCap) nClipped++;
        }

        double mean = sum / n;
        double meanU = sumU / n;
        double var = Math.max(0.0,
                (sum2 - n * mean * mean) / (n - 1));
        double varU = Math.max(0.0,
                (sumU2 - n * meanU * meanU) / (n - 1));
        double deltaPer = perEventDelta(true);
        double clippedDelta = solveBernsteinDelta(
                n, var, 1.0, deltaPer);
        double excessDelta = uRange > 0.0
                ? solveBernsteinDelta(n, varU, uRange, deltaPer)
                : 0.0;
        double pHat = (double) nClipped / n;
        double pUpper = Math.min(1.0,
                pHat + hoeffdingHalfWidth(n, deltaPer));
        double biasWorst = uRange * pUpper;
        double biasEmpirical = meanU + excessDelta;
        double biasUpper = Math.min(biasWorst, biasEmpirical);
        if (!Double.isFinite(mean) || !Double.isFinite(var)
                || !Double.isFinite(clippedDelta)
                || !Double.isFinite(biasUpper)) {
            failCertificate("conditional eta-v4: non-finite final PAC statistics",
                    -1L);
            return;
        }

        double meanLower = Math.max(0.0, mean - clippedDelta);
        double meanUpper = mean + clippedDelta + biasUpper;
        logZLowerPAC = meanLower > 0.0
                ? logZReference + clipLogCap + Math.log(meanLower)
                : Double.NEGATIVE_INFINITY;
        logZUpperPAC = meanUpper > 0.0
                ? logZReference + clipLogCap + Math.log(meanUpper)
                : Double.NaN;
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);

        meanResidual = sumRes / n;
        stdResidual = Math.sqrt(Math.max(0.0,
                (sumRes2 - n * meanResidual * meanResidual) / (n - 1)));
        meanPsi = shiftedMeanToDouble(mean + meanU, clipLogCap);
        varPsi = var;
        cvPsi = mean > 0.0 ? Math.sqrt(var) / mean : Double.MAX_VALUE;
        emitResidualSampleTraces(finalSamples, xi, logW,
                logZReference, clipLogCap);

        certificateValid = isValidCertificate(zLower, zUpper, epsilon);
        certificateFailureReason = certificateValid
                ? "" : "conditional eta-v4 mixture interval failed final validation";
        System.out.println("[PACK*-conditional-eta-v4] final certificate:"
                + " certificateMode="
                + (design.componentCover
                ? "component-cover" : "single-anchor")
                + ", anchorB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                design.anchorBoundKcal)
                + ", anchorWeight="
                + String.format(Locale.ROOT, "%.6f",
                design.weights[design.anchorIndex])
                + ", effectiveB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                design.effectiveBoundKcal)
                + ", logCap="
                + String.format(Locale.ROOT, "%.6f", clipLogCap)
                + ", clippedDelta="
                + String.format(Locale.ROOT, "%.6f", clippedDelta)
                + ", excessMean="
                + String.format(Locale.ROOT, "%.6f", meanU)
                + ", excessRange="
                + String.format(Locale.ROOT, "%.6f", uRange)
                + ", excessUpper="
                + String.format(Locale.ROOT, "%.6f", biasUpper)
                + ", nClipped=" + nClipped + "/" + n
                + ", upperVia=conditional-iid-mixture-sampling"
                + ", q_mFallback=false"
                + ", log10Lower=" + formatLog10(logZLowerPAC)
                + ", log10Upper=" + formatLog10(logZUpperPAC)
                + ", epsilon="
                + String.format(Locale.ROOT, "%.8f", epsilon)
                + ", certified=" + certificateValid);
    }

    private void runConditionalTwoStagePAC(long startTime, double logZRigid) {
        Random rng = new Random(randomSeed);

        // ---- Stage A: train from p_m -> eta, corrected DP (p_eta, q_eta) ----
        long tA = System.currentTimeMillis();
        List<int[]> trainConfs = sampleConformationsFromDP(trainSamples, rng);
        List<CCDResult> trainCCD = runParallelCCD(trainConfs);
        System.out.println("[PACK*-2stage] Stage A train: " + trainCCD.size()
                + " CCD from p_m in " + (System.currentTimeMillis() - tA) + " ms");
        if (trainCCD.isEmpty()) {
            setZeroBounds("two-stage: no valid train samples");
            printFinalSummary(startTime);
            return;
        }

        EtaCorrections eta = etaEnabled ? extractEtaCorrections(trainCCD) : zeroEtaCorrections();
        if (!etaEnabled) {
            System.out.println("[PACK*-2stage] no-eta ablation active (packstar.pac.etaEnabled=false):"
                    + " eta fixed to 0, E_eta=E_m, q_eta=q_m; Stage B reweights the full gap g=E_true-E_m.");
        }
        EnergyMatrix correctedEmat = buildCorrectedEmat(eta);
        maybeDumpTrainingSamples(trainCCD, eta, correctedEmat);
        long tDP = System.currentTimeMillis();
        CorrectedDPResult correctedDP = recomputeDP(correctedEmat, eta); // loads p_eta into the tree
        double logZCorrected = correctedDP.logZCorrected;
        System.out.println("[PACK*-2stage] corrected DP solved ONCE: logZ_corrected="
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
        System.out.println("[PACK*-2stage] Stage B pilot: " + pilotCCD.size()
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
        boolean iterate = getConfigBoolean("packstar.pac.iterate", DEFAULT_ITERATE);
        double collapseThresh = getConfigDouble("packstar.pac.iterate.meanWThreshold", 0.3);
        double driftFracThresh = getConfigDouble("packstar.pac.iterate.driftFraction", 0.2);
        int minTrainCount = getConfigInteger("packstar.pac.iterate.minTrainCount", 5);
        int maxRounds = getConfigInteger("packstar.pac.iterate.maxRounds", DEFAULT_ITERATE_MAX_ROUNDS);
        // No-eta ablation must stay at eta===0 in every round; distribution-shift
        // refinement re-learns a real eta, so it is force-disabled here rather than
        // relying on callers to also pass packstar.pac.iterate=false.
        if (iterate && etaEnabled) {
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
                System.out.println("[PACK*-2stage] refinement check (round " + round + "/" + maxRounds
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
                    System.out.println("[PACK*-2stage] refinement round " + round
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
                System.out.println("[PACK*-2stage] distribution-shift refinement round " + round
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
        System.out.println("[PACK*-2stage] Stage B fresh estimation: " + estCCD.size() + "/" + nStar
                + " CCD from p_eta in " + (System.currentTimeMillis() - tEst) + " ms"
                + " (pilot held out)");
        totalCCDCalls = trainCCD.size() + pilotCCD.size() + estCCD.size() + extraRefineCCD;

        // ---- clip threshold (if enabled) is fixed on the held-out pilot, BEFORE the
        // estimation weights S2 are observed, so the clipped weights stay i.i.d. given c ----
        if (getConfigBoolean("packstar.pac.clip", DEFAULT_CLIP)) {
            double clipQuantile = getConfigDouble("packstar.pac.clipQuantile",
                    DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE);
            double[] plogW = new double[pilotCCD.size()];
            pilotResidualBoundViolations = 0;
            pilotMaxOverCorrectionKcal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < plogW.length; i++) {
                CCDResult r = pilotCCD.get(i);
                plogW[i] = -(r.eTrue - computeFullConfPairwiseEnergy(r.conf, correctedEmat)) / RT;
            }
            clipLogCap = quantile(plogW, clipQuantile);
            if (Double.isFinite(residualBoundKcal)) {
                for (double pilotLogW : plogW) {
                    if (pilotLogW > clipLogCap) {
                        double overCorrectionKcal = RT * pilotLogW;
                        pilotMaxOverCorrectionKcal = Math.max(
                                pilotMaxOverCorrectionKcal,
                                overCorrectionKcal
                        );
                        if (overCorrectionKcal > residualBoundKcal) {
                            pilotResidualBoundViolations++;
                        }
                    }
                }
            }
            System.out.println("[PACK*-2stage] clip threshold from pilot: logCap="
                    + String.format("%.4f", clipLogCap) + " (clipQ=" + clipQuantile
                    + ", nPilot=" + plogW.length + ")");
        }

        // ---- residual-leg empirical-Bernstein bound: Z = q_eta * E_{p_eta}[exp(-xi/RT)] ----
        computePACBoundResidual(estCCD, correctedEmat, logZCorrected);
        printFinalSummary(startTime);
    }

    private boolean validateObservedLowerBound(List<CCDResult> samples, String stage) {
        int violations = 0;
        double minGap = Double.POSITIVE_INFINITY;
        for (CCDResult sample : samples) {
            if (!Double.isFinite(sample.eTrue) || !Double.isFinite(sample.eMin)) {
                System.out.println("[PACK*-lower-bound] ERROR: non-finite energy at " + stage);
                return false;
            }
            double gap = sample.eTrue - sample.eMin;
            minGap = Math.min(minGap, gap);
            if (gap < 0.0) violations++;
        }
        System.out.println("[PACK*-lower-bound] audit " + stage
                + ": samples=" + samples.size()
                + ", min(E_true-E_m)="
                + String.format(Locale.ROOT, "%.9f", minGap)
                + " kcal/mol, strictThreshold=0"
                + ", violations=" + violations);
        if (violations > 0) {
            System.out.println("[PACK*-lower-bound] ERROR: observed E_m > E_true. "
                    + "This falsifies the structural lower-bound premise for the configured "
                    + "Hamiltonian/interaction set; refusing certification.");
        }
        return violations == 0;
    }

    private double maxOverCorrectionKcal(List<CCDResult> samples,
                                         EnergyMatrix correctedEmat) {
        double max = Double.NEGATIVE_INFINITY;
        for (CCDResult sample : samples) {
            double eEta = computeFullConfPairwiseEnergy(sample.conf, correctedEmat);
            max = Math.max(max, eEta - sample.eTrue);
        }
        return max;
    }

    /**
     * Calibrate only the energy gauge of an eta-v4 anchor on adaptation data.
     * A constant downward shift changes neither the exact-DP proposal nor any
     * IID mixture labels: E' = E-d and log Z' = log Z+d/RT cancel exactly in
     * log q.  It only moves newly observed residuals below the pre-registered
     * proof headroom.  This method must never be called from monitor/final.
     */
    private ConditionalProposalSnapshot calibrateConditionalEtaV4AnchorGauge(
            ConditionalProposalSnapshot anchor,
            List<CCDResult> adaptation, String stage) {
        if (anchor == null || anchor.eta == null || anchor.emat == null
                || !Double.isFinite(anchor.logZ)
                || adaptation == null || adaptation.isEmpty()
                || stage == null || stage.isEmpty()) {
            return null;
        }
        double observedMax = maxOverCorrectionKcal(
                adaptation, anchor.emat);
        double targetMax = etaV4ProofBoundKcal
                - etaV4HeadroomKcal;
        if (!Double.isFinite(observedMax)
                || !Double.isFinite(targetMax)) return null;
        double additionalShiftDown = Math.max(
                0.0, observedMax - targetMax);
        EtaCorrections calibratedEta = anchor.eta;
        EnergyMatrix calibratedEmat = anchor.emat;
        double calibratedLogZ = anchor.logZ;
        if (additionalShiftDown > 0.0) {
            calibratedEta = anchor.eta.shiftedConstant(
                    -additionalShiftDown);
            calibratedEmat = buildCorrectedEmat(calibratedEta);
            calibratedLogZ = anchor.logZ
                    + additionalShiftDown / RT;
        }
        double calibratedMax = observedMax - additionalShiftDown;
        boolean valid = Double.isFinite(calibratedLogZ)
                && Double.isFinite(calibratedMax)
                && calibratedMax <= targetMax + 1.0e-10;
        System.out.println("[PACK*-conditional-eta-v4] adaptation-only anchor gauge:"
                + " stage=" + stage
                + ", samples=" + adaptation.size()
                + ", observedMaxBeforeKcal="
                + String.format(Locale.ROOT, "%.6f", observedMax)
                + ", additionalShiftDownKcal="
                + String.format(Locale.ROOT, "%.6f",
                additionalShiftDown)
                + ", calibratedMaxKcal="
                + String.format(Locale.ROOT, "%.6f", calibratedMax)
                + ", targetMaxKcal="
                + String.format(Locale.ROOT, "%.6f", targetMax)
                + ", proposalDistributionChanged=false"
                + ", logZUpdatedAnalytically=true"
                + ", valid=" + valid);
        if (!valid) return null;
        return new ConditionalProposalSnapshot(
                anchor.sourceProposalId,
                calibratedEta, calibratedEmat,
                calibratedLogZ, anchor.reportedAlpha);
    }

    private static List<CCDResult> conditionalSamplesFromSource(
            List<CCDResult> samples, String sourceProposalId) {
        List<CCDResult> selected = new ArrayList<>();
        if (samples == null || sourceProposalId == null) return selected;
        for (CCDResult sample : samples) {
            if (sourceProposalId.equals(sample.sourceProposalId)) {
                selected.add(sample);
            }
        }
        return selected;
    }

    /**
     * Choose a conditional anchor without requiring it to be independently
     * reachable.  Eta-v4 may rescue its width with a second-moment component;
     * the anchor's role is to retain a valid B premise and a usable on-policy
     * batch while minimizing the standalone predicted epsilon.
     */
    private ConditionalEtaV4Anchor selectConditionalEtaV4Anchor(
            List<ConditionalProposalSnapshot> proposalArchive,
            List<CCDResult> adaptation) {
        if (proposalArchive == null || proposalArchive.isEmpty()
                || adaptation == null || adaptation.isEmpty()) return null;

        ConditionalProposalSnapshot bestProposal = null;
        String bestSourceId = null;
        List<CCDResult> bestOnPolicy = null;
        double bestMax = Double.NaN;
        double bestChiSquare = Double.NaN;
        double bestEss = Double.NaN;
        double bestEpsilon = Double.POSITIVE_INFINITY;
        int bestFinalSamples = Integer.MAX_VALUE;
        double bestLogRange = Double.POSITIVE_INFINITY;
        double targetMax = etaV4ProofBoundKcal - etaV4HeadroomKcal;
        int candidateIndex = 0;

        for (ConditionalProposalSnapshot snapshot : proposalArchive) {
            candidateIndex++;
            List<CCDResult> onPolicy = conditionalSamplesFromSource(
                    adaptation, snapshot.sourceProposalId);
            // Every adaptation CCD observation is a valid falsification
            // witness for the anchor's global B premise, regardless of which
            // archived proposal generated it.  Only variance statistics below
            // are restricted to the genuinely on-policy batch.
            double unshiftedMax = onPolicy.isEmpty()
                    ? Double.NaN
                    : maxOverCorrectionKcal(adaptation, snapshot.emat);
            double shiftDown = Double.isFinite(unshiftedMax)
                    ? Math.max(0.0, unshiftedMax - targetMax)
                    : Double.NaN;
            EtaCorrections candidateEta = Double.isFinite(shiftDown)
                    ? snapshot.eta.shiftedConstant(-shiftDown) : null;
            EnergyMatrix candidateEmat = candidateEta == null
                    ? null : buildCorrectedEmat(candidateEta);
            double candidateLogZ = Double.isFinite(shiftDown)
                    ? snapshot.logZ + shiftDown / RT : Double.NaN;
            double candidateMax = Double.isFinite(shiftDown)
                    ? unshiftedMax - shiftDown : Double.NaN;
            ConditionalAdaptationStats stats = candidateEmat == null
                    ? null : buildConditionalAdaptationStats(
                    onPolicy, candidateEmat, candidateLogZ,
                    etaV4ProofBoundKcal);
            boolean essPass = stats != null
                    && conditionalDpRichEssPasses(
                    stats.effectiveSampleSize, onPolicy.size(),
                    dpRichMinEssFraction);
            double candidateClip = stats == null ? Double.NaN
                    : selectConditionalClipLogCap(
                    stats, etaV4ProofBoundKcal);
            ConditionalSizing sizing = stats != null
                    && Double.isFinite(candidateClip)
                    ? computeConditionalSizing(stats,
                    etaV4ProofBoundKcal, candidateClip) : null;
            double chiSquare = stats == null ? Double.NaN
                    : conditionalWeightChiSquareProxy(
                    stats.logResidualWeights,
                    stats.normalizedProposalWeights);
            double logRange = conditionalLogRangeProxy(
                    candidateLogZ, etaV4ProofBoundKcal, RT);
            boolean eligible = essPass && sizing != null
                    && Double.isFinite(candidateMax)
                    && candidateMax
                    <= etaV4ProofBoundKcal + 1.0e-12
                    && Double.isFinite(chiSquare)
                    && Double.isFinite(logRange);
            boolean better = false;
            if (eligible) {
                double tolerance = 1.0e-12;
                if (sizing.epsilonAtMaxSamples
                        < bestEpsilon - tolerance) {
                    better = true;
                } else if (Math.abs(sizing.epsilonAtMaxSamples
                        - bestEpsilon) <= tolerance
                        && sizing.finalSamples < bestFinalSamples) {
                    better = true;
                } else if (Math.abs(sizing.epsilonAtMaxSamples
                        - bestEpsilon) <= tolerance
                        && sizing.finalSamples == bestFinalSamples
                        && (!Double.isFinite(bestChiSquare)
                        || chiSquare < bestChiSquare - tolerance)) {
                    better = true;
                } else if (Math.abs(sizing.epsilonAtMaxSamples
                        - bestEpsilon) <= tolerance
                        && sizing.finalSamples == bestFinalSamples
                        && Math.abs(chiSquare - bestChiSquare) <= tolerance
                        && logRange < bestLogRange - tolerance) {
                    better = true;
                }
            }

            System.out.println("[PACK*-conditional-eta-v4] anchor candidate "
                    + candidateIndex + "/" + proposalArchive.size()
                    + ": source=" + snapshot.sourceProposalId
                    + ", onPolicySamples=" + onPolicy.size()
                    + ", gaugeShiftDownKcal="
                    + String.format(Locale.ROOT, "%.6f", shiftDown)
                    + ", adaptationMaxOverCorrectionKcal="
                    + String.format(Locale.ROOT, "%.6f", candidateMax)
                    + ", B_kcal="
                    + String.format(Locale.ROOT, "%.6f",
                    etaV4ProofBoundKcal)
                    + ", ESS=" + (stats == null ? "NaN"
                    : String.format(Locale.ROOT, "%.6f",
                    stats.effectiveSampleSize))
                    + ", chiSquareProxy="
                    + String.format(Locale.ROOT, "%.6f", chiSquare)
                    + ", predictedFinalN=" + (sizing == null ? "NA"
                    : Integer.toString(sizing.finalSamples))
                    + ", predictedEpsilonAtMax=" + (sizing == null
                    ? "NaN" : String.format(Locale.ROOT, "%.6f",
                    sizing.epsilonAtMaxSamples))
                    + ", eligible=" + eligible
                    + ", bestSoFar=" + better);
            if (!better) continue;

            bestProposal = new ConditionalProposalSnapshot(
                    "conditional-v4-anchor",
                    candidateEta, candidateEmat,
                    candidateLogZ, snapshot.reportedAlpha);
            bestSourceId = snapshot.sourceProposalId;
            bestOnPolicy = onPolicy;
            bestMax = candidateMax;
            bestChiSquare = chiSquare;
            bestEss = stats.effectiveSampleSize;
            bestEpsilon = sizing.epsilonAtMaxSamples;
            bestFinalSamples = sizing.finalSamples;
            bestLogRange = logRange;
        }

        if (bestProposal == null || bestOnPolicy == null) return null;
        CorrectedDPResult winnerDP = recomputeDP(
                bestProposal.emat, bestProposal.eta);
        if (!Double.isFinite(winnerDP.logZCorrected)) return null;
        bestProposal = new ConditionalProposalSnapshot(
                bestProposal.sourceProposalId,
                bestProposal.eta, bestProposal.emat,
                winnerDP.logZCorrected,
                bestProposal.reportedAlpha);
        System.out.println("[PACK*-conditional-eta-v4] anchor winner: source="
                + bestSourceId
                + ", adaptationMaxOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f", bestMax)
                + ", B_kcal=" + String.format(Locale.ROOT, "%.6f",
                etaV4ProofBoundKcal)
                + ", chiSquareProxy="
                + String.format(Locale.ROOT, "%.6f", bestChiSquare)
                + ", ESS=" + String.format(Locale.ROOT, "%.6f", bestEss)
                + ", predictedFinalN=" + bestFinalSamples
                + ", predictedEpsilonAtMax="
                + String.format(Locale.ROOT, "%.6f", bestEpsilon));
        return new ConditionalEtaV4Anchor(
                bestProposal, bestSourceId, bestOnPolicy,
                bestMax, bestChiSquare, bestEss,
                bestEpsilon, bestFinalSamples, 1);
    }

    /** Build one OAIS round by direct second-moment descent. */
    private ConditionalEtaV4Tail selectConditionalEtaV4Tail(
            ConditionalProposalSnapshot baseProposal,
            List<CCDResult> onPolicySamples, int roundIndex) {
        if (baseProposal == null || onPolicySamples == null
                || onPolicySamples.size() < 2 || roundIndex < 1) return null;
        List<CCDResult> gradientSamples = new ArrayList<>();
        List<CCDResult> screeningSamples = new ArrayList<>();
        for (int i = 0; i < onPolicySamples.size(); i++) {
            if ((i & 1) == 0) gradientSamples.add(
                    onPolicySamples.get(i));
            else screeningSamples.add(onPolicySamples.get(i));
        }
        if (gradientSamples.size() < 2 || screeningSamples.size() < 2) {
            return null;
        }
        EtaV4SparseFeatures features = buildEtaV4SparseFeatures(
                gradientSamples);
        if (features == null || features.cells.length == 0) return null;

        int n = gradientSamples.size();
        double[] baseWeights = new double[n];
        Arrays.fill(baseWeights, 1.0);
        double[] logImportance = new double[n];
        for (int i = 0; i < n; i++) {
            CCDResult sample = gradientSamples.get(i);
            double anchorEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, baseProposal.emat);
            logImportance[i] = (anchorEnergy - sample.eTrue) / RT;
            if (!Double.isFinite(logImportance[i])) return null;
        }
        EtaV4OffPolicyScore baseline = scoreEtaV4CandidateOffPolicy(
                screeningSamples,
                baseProposal.emat, baseProposal.logZ,
                baseProposal.emat, baseProposal.logZ);
        if (baseline == null) return null;

        ConditionalProposalSnapshot bestProposal = null;
        double bestPower = Double.NaN;
        double bestTrust = Double.NaN;
        double bestGradientEss = Double.NaN;
        double bestProposalEss = Double.NaN;
        double bestImprovement = 0.0;
        int dpSweeps = 0;
        double minTailEss = Math.max(2.0,
                etaV4MinTailEssFraction * n);
        // ESS is a floating-point reduction of normalized exponentials.  Do
        // not reject the intended absolute ESS=2 boundary because of a few
        // ulps; the independent held-out objective remains the acceptance
        // gate for every resulting proposal.
        double minTailEssTolerance = Math.max(
                1.0e-6, 1.0e-4 * minTailEss);

        for (int powerIndex = etaV4TailPowers.length - 1;
             powerIndex >= 0; powerIndex--) {
            double tailPower = etaV4TailPowers[powerIndex];
            PackStarEtaV4Importance.SparseGradient sparse =
                    PackStarEtaV4Importance.sparseChiSquareGradient(
                            features.cells.length,
                            features.sampleIndices,
                            baseWeights, logImportance, tailPower);
            boolean gradientEligible = sparse != null
                    && PackStarEtaV4Importance.tailEssPasses(
                    sparse.tailEffectiveSampleSize, minTailEss);
            System.out.println("[PACK*-conditional-eta-v4] gradient:"
                    + " oaisRound=" + roundIndex
                    + ", tailPower="
                    + String.format(Locale.ROOT, "%.6f", tailPower)
                    + ", features=" + features.cells.length
                    + ", fitSamples=" + n
                    + ", heldOutScreenSamples="
                    + screeningSamples.size()
                    + ", tailESS=" + (sparse == null ? "NaN"
                    : String.format(Locale.ROOT, "%.6f",
                    sparse.tailEffectiveSampleSize))
                    + ", threshold="
                    + String.format(Locale.ROOT, "%.6f", minTailEss)
                    + ", tolerance="
                    + String.format(Locale.ROOT, "%.6g",
                    minTailEssTolerance)
                    + ", eligible=" + gradientEligible);
            if (!gradientEligible) continue;

            for (double trustKcal : etaV4TrustKcal) {
                EtaV4ProposalUpdate update = applyEtaV4EnergyGradient(
                        baseProposal.eta, features,
                        sparse.gradient, trustKcal);
                if (update == null) continue;
                EnergyMatrix candidateEmat = buildCorrectedEmat(update.eta);
                CorrectedDPResult candidateDP = recomputeDP(
                        candidateEmat, update.eta);
                dpSweeps++;
                EtaV4OffPolicyScore score =
                        scoreEtaV4CandidateOffPolicy(
                                screeningSamples,
                                baseProposal.emat,
                                baseProposal.logZ,
                                candidateEmat,
                                candidateDP.logZCorrected);
                double improvement = score == null ? Double.NaN
                        : baseline.logSecondMoment
                        - score.logSecondMoment;
                boolean proposalEssPass = score != null
                        && conditionalDpRichEssPasses(
                        score.proposalEffectiveSampleSize,
                        screeningSamples.size(),
                        dpRichMinEssFraction);
                boolean better = proposalEssPass
                        && Double.isFinite(improvement)
                        && improvement > bestImprovement + 1.0e-12;
                System.out.println("[PACK*-conditional-eta-v4] tail candidate:"
                        + " oaisRound=" + roundIndex
                        + ", tailPower="
                        + String.format(Locale.ROOT, "%.6f", tailPower)
                        + ", trustKcal="
                        + String.format(Locale.ROOT, "%.6f", trustKcal)
                        + ", maxObservedDeltaKcal="
                        + String.format(Locale.ROOT, "%.6f",
                        update.maxObservedEnergyDeltaKcal)
                        + ", maxCellDeltaKcal="
                        + String.format(Locale.ROOT, "%.6f",
                        update.maxCellDeltaKcal)
                        + ", offPolicyESS=" + (score == null ? "NaN"
                        : String.format(Locale.ROOT, "%.6f",
                        score.proposalEffectiveSampleSize))
                        + ", logSecondMomentImprovement="
                        + String.format(Locale.ROOT, "%.9f", improvement)
                        + ", eligible=" + proposalEssPass
                        + ", bestSoFar=" + better);
                if (!better) continue;
                bestProposal = new ConditionalProposalSnapshot(
                        "conditional-v4-oais-" + roundIndex,
                        update.eta, candidateEmat,
                        candidateDP.logZCorrected, 1.0);
                bestPower = tailPower;
                bestTrust = trustKcal;
                bestGradientEss = sparse.tailEffectiveSampleSize;
                bestProposalEss = score.proposalEffectiveSampleSize;
                bestImprovement = improvement;
            }
        }

        if (bestProposal == null) {
            recomputeDP(baseProposal.emat, baseProposal.eta);
            System.out.println("[PACK*-conditional-eta-v4] no tail candidate"
                    + " improved the direct off-policy second moment"
                    + " at OAIS round " + roundIndex + ";"
                    + " retaining the current proposal and stopping OAIS");
            return null;
        }
        CorrectedDPResult winnerDP = recomputeDP(
                bestProposal.emat, bestProposal.eta);
        dpSweeps++;
        if (!Double.isFinite(winnerDP.logZCorrected)) return null;
        bestProposal = new ConditionalProposalSnapshot(
                bestProposal.sourceProposalId,
                bestProposal.eta, bestProposal.emat,
                winnerDP.logZCorrected,
                bestProposal.reportedAlpha);
        System.out.println("[PACK*-conditional-eta-v4] tail winner:"
                + " oaisRound=" + roundIndex
                + ", tailPower="
                + String.format(Locale.ROOT, "%.6f", bestPower)
                + ", trustKcal="
                + String.format(Locale.ROOT, "%.6f", bestTrust)
                + ", gradientTailESS="
                + String.format(Locale.ROOT, "%.6f", bestGradientEss)
                + ", offPolicyProposalESS="
                + String.format(Locale.ROOT, "%.6f", bestProposalEss)
                + ", logSecondMomentImprovement="
                + String.format(Locale.ROOT, "%.9f", bestImprovement)
                + ", candidateDpSweeps=" + dpSweeps);
        return new ConditionalEtaV4Tail(bestProposal,
                Collections.unmodifiableList(
                        new ArrayList<>(screeningSamples)),
                bestPower, bestTrust, bestGradientEss,
                bestProposalEss, bestImprovement, dpSweeps,
                roundIndex);
    }

    private double[] conditionalEtaV4LogComponentProbabilities(
            int[] conf, ConditionalProposalSnapshot[] components) {
        if (conf == null || components == null
                || components.length == 0) return null;
        double[] logProbabilities = new double[components.length];
        for (int k = 0; k < components.length; k++) {
            ConditionalProposalSnapshot component = components[k];
            if (component == null || component.emat == null
                    || !Double.isFinite(component.logZ)) return null;
            double energy = computeFullConfPairwiseEnergy(
                    conf, component.emat);
            logProbabilities[k] = -energy / RT - component.logZ;
            if (!Double.isFinite(logProbabilities[k])) return null;
        }
        return logProbabilities;
    }

    /**
     * Build candidate-q_mix adaptation statistics from a pooled balance-MIS
     * reference.  Null reference weights mean the samples are already IID
     * from candidate q_mix and therefore receive uniform weight.
     */
    private ConditionalAdaptationStats buildConditionalEtaV4Stats(
            List<CCDResult> samples,
            ConditionalProposalSnapshot[] components,
            double[] candidateWeights,
            double[] referenceWeights,
            int anchorIndex, double effectiveBoundKcal) {
        if (samples == null || samples.size() < 2
                || components == null || components.length == 0
                || candidateWeights == null
                || candidateWeights.length != components.length
                || anchorIndex < 0 || anchorIndex >= components.length
                || !Double.isFinite(effectiveBoundKcal)
                || !(effectiveBoundKcal > 0.0)
                || (referenceWeights != null
                && referenceWeights.length != components.length)) {
            return null;
        }
        double anchorLogZ = components[anchorIndex].logZ;
        if (!Double.isFinite(anchorLogZ)) return null;
        double[] logResidualWeights = new double[samples.size()];
        double[] logProposalRatios = referenceWeights == null
                ? null : new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            CCDResult sample = samples.get(i);
            if (!Double.isFinite(sample.eTrue)) return null;
            double[] logComponents =
                    conditionalEtaV4LogComponentProbabilities(
                            sample.conf, components);
            if (logComponents == null) return null;
            double logCandidate = PackStarEtaV4Importance
                    .logMixtureProbability(
                            logComponents, candidateWeights);
            if (!Double.isFinite(logCandidate)) return null;
            double logResidual = -sample.eTrue / RT
                    - logCandidate - anchorLogZ;
            if (!Double.isFinite(logResidual)
                    || logResidual
                    > effectiveBoundKcal / RT + 1.0e-12) {
                return null;
            }
            logResidualWeights[i] = logResidual;
            if (referenceWeights != null) {
                double logReference = PackStarEtaV4Importance
                        .logMixtureProbability(
                                logComponents, referenceWeights);
                if (!Double.isFinite(logReference)) return null;
                logProposalRatios[i] = logCandidate - logReference;
            }
        }

        if (referenceWeights == null) {
            double[] uniform = new double[samples.size()];
            Arrays.fill(uniform, 1.0 / samples.size());
            return new ConditionalAdaptationStats(
                    logResidualWeights, uniform,
                    samples.size(), 0.0, 0.0,
                    components.length);
        }
        ConditionalImportanceWeights normalized =
                normalizeConditionalImportanceWeights(
                        logProposalRatios);
        if (normalized == null) return null;
        return new ConditionalAdaptationStats(
                logResidualWeights, normalized.normalized,
                normalized.effectiveSampleSize,
                normalized.minLogWeight,
                normalized.maxLogWeight,
                components.length);
    }

    private ConditionalEtaV4MixtureDesign selectConditionalEtaV4Mixture(
            ConditionalProposalSnapshot[] components,
            List<CCDResult> pooledAdaptation,
            double[] referenceWeights,
            int candidateDpSweeps) {
        if (components == null || components.length == 0
                || pooledAdaptation == null
                || pooledAdaptation.size() < 2
                || referenceWeights == null
                || referenceWeights.length != components.length) {
            return null;
        }
        double bestEpsilon = Double.POSITIVE_INFINITY;
        int bestFinalSamples = Integer.MAX_VALUE;
        double bestChiSquare = Double.POSITIVE_INFINITY;
        double bestAnchorWeight = Double.NaN;
        int bestAnchorIndex = -1;
        boolean bestComponentCover = false;
        boolean bestPreferredCover = false;
        double[] bestWeights = null;
        double bestEffectiveBound = Double.NaN;
        double bestAnchorLogRange = Double.NaN;
        ConditionalAdaptationStats bestStats = null;
        double bestClip = Double.NaN;
        ConditionalSizing bestSizing = null;

        double[] candidateAnchorWeights = components.length == 1
                ? new double[]{1.0} : etaV4AnchorWeights;
        double[] componentLogZ = new double[components.length];
        double[] componentBounds = new double[components.length];
        for (int k = 0; k < components.length; k++) {
            componentLogZ[k] = components[k].logZ;
            componentBounds[k] = etaV4ProofBoundKcal;
        }
        double[] rangeBalancedWeights = components.length > 1
                ? PackStarEtaV4Importance.rangeBalancedCoverWeights(
                componentLogZ, componentBounds, RT) : null;
        int choicesPerReference = candidateAnchorWeights.length
                + (rangeBalancedWeights == null ? 0 : 1);
        int candidateIndex = 0;
        int totalCandidates = components.length
                * choicesPerReference;
        for (int anchorIndex = 0;
             anchorIndex < components.length; anchorIndex++) {
            for (int choice = 0; choice < choicesPerReference; choice++) {
                candidateIndex++;
                boolean rangeBalanced = choice
                        == candidateAnchorWeights.length;
                double[] weights;
                double anchorWeight;
                if (rangeBalanced) {
                    weights = rangeBalancedWeights.clone();
                    anchorWeight = weights[anchorIndex];
                } else {
                    int weightIndex = candidateAnchorWeights.length
                            - 1 - choice;
                    anchorWeight = candidateAnchorWeights[weightIndex];
                    weights = new double[components.length];
                    weights[anchorIndex] = anchorWeight;
                    if (components.length == 2) {
                        weights[1 - anchorIndex] = 1.0 - anchorWeight;
                    } else if (components.length > 2) {
                        double remainder = (1.0 - anchorWeight)
                                / (components.length - 1.0);
                        for (int k = 0; k < components.length; k++) {
                            if (k != anchorIndex) weights[k] = remainder;
                        }
                    }
                }
                boolean componentCover = components.length > 1;
                for (double weight : weights) {
                    componentCover &= weight > 0.0;
                }
                double anchorLogRange =
                        componentCover
                        ? PackStarEtaV4Importance.componentCoverLogRange(
                        componentLogZ, componentBounds, RT, weights)
                        : PackStarEtaV4Importance.anchorLogRange(
                                components[anchorIndex].logZ,
                                etaV4ProofBoundKcal, RT,
                                anchorWeight);
                double effectiveBound = Double.isFinite(anchorLogRange)
                        ? RT * (anchorLogRange
                        - components[anchorIndex].logZ)
                        : Double.NaN;
                int coverViolations = componentCover
                        ? conditionalEtaV4CoverViolations(
                        pooledAdaptation, components, weights,
                        etaV4ProofBoundKcal)
                        : 0;
                ConditionalAdaptationStats stats =
                        coverViolations == 0
                        && Double.isFinite(effectiveBound)
                        && effectiveBound > 0.0
                        ? buildConditionalEtaV4Stats(
                                pooledAdaptation, components,
                                weights, referenceWeights,
                                anchorIndex, effectiveBound) : null;
                boolean essPass = stats != null
                        && conditionalDpRichEssPasses(
                        stats.effectiveSampleSize,
                        pooledAdaptation.size(),
                        dpRichMinEssFraction);
                double candidateClip = stats == null ? Double.NaN
                        : selectConditionalClipLogCap(
                        stats, effectiveBound);
                ConditionalSizing sizing = stats != null
                        && Double.isFinite(candidateClip)
                        ? computeConditionalSizing(
                        stats, effectiveBound, candidateClip) : null;
                double chiSquare = stats == null ? Double.NaN
                        : conditionalWeightChiSquareProxy(
                        stats.logResidualWeights,
                        stats.normalizedProposalWeights);
                boolean valid = essPass && sizing != null
                        && Double.isFinite(chiSquare)
                        && Double.isFinite(anchorLogRange)
                        && coverViolations == 0;
                boolean preferredCover = valid && componentCover
                        && sizing.epsilonAtMaxSamples
                        <= targetEpsilon + 1.0e-12;
                boolean better = false;
                if (valid) {
                    double tolerance = 1.0e-12;
                    if (preferredCover && !bestPreferredCover) {
                        better = true;
                    } else if (!preferredCover && bestPreferredCover) {
                        better = false;
                    } else if (sizing.epsilonAtMaxSamples
                            < bestEpsilon - tolerance) {
                        better = true;
                    } else if (Math.abs(sizing.epsilonAtMaxSamples
                            - bestEpsilon) <= tolerance
                            && sizing.finalSamples
                            < bestFinalSamples) {
                        better = true;
                    } else if (Math.abs(sizing.epsilonAtMaxSamples
                            - bestEpsilon) <= tolerance
                            && sizing.finalSamples
                            == bestFinalSamples
                            && chiSquare
                            < bestChiSquare - tolerance) {
                        better = true;
                    } else if (Math.abs(sizing.epsilonAtMaxSamples
                            - bestEpsilon) <= tolerance
                            && sizing.finalSamples
                            == bestFinalSamples
                            && Math.abs(chiSquare - bestChiSquare)
                            <= tolerance
                            && (!Double.isFinite(bestAnchorWeight)
                            || anchorWeight
                            > bestAnchorWeight + tolerance)) {
                        better = true;
                    } else if (Math.abs(sizing.epsilonAtMaxSamples
                            - bestEpsilon) <= tolerance
                            && sizing.finalSamples
                            == bestFinalSamples
                            && Math.abs(chiSquare - bestChiSquare)
                            <= tolerance
                            && Math.abs(anchorWeight
                            - bestAnchorWeight) <= tolerance
                            && (bestAnchorIndex < 0
                            || anchorIndex < bestAnchorIndex)) {
                        better = true;
                    }
                }
                System.out.println("[PACK*-conditional-eta-v4] mixture candidate "
                        + candidateIndex + "/" + totalCandidates
                        + ": anchorIndex=" + anchorIndex
                        + ", anchorSource="
                        + components[anchorIndex].sourceProposalId
                        + ", anchorWeight="
                        + String.format(Locale.ROOT, "%.6f",
                        anchorWeight)
                        + ", weights=" + Arrays.toString(weights)
                        + ", certificateMode="
                        + (componentCover
                        ? "component-cover" : "single-anchor")
                        + ", rangeBalanced=" + rangeBalanced
                        + ", coverViolations=" + coverViolations
                        + ", effectiveB_kcal="
                        + String.format(Locale.ROOT, "%.6f",
                        effectiveBound)
                        + ", ESS=" + (stats == null ? "NaN"
                        : String.format(Locale.ROOT, "%.6f",
                        stats.effectiveSampleSize))
                        + ", chiSquareProxy="
                        + String.format(Locale.ROOT, "%.6f",
                        chiSquare)
                        + ", predictedFinalN="
                        + (sizing == null ? "NA"
                        : Integer.toString(sizing.finalSamples))
                        + ", predictedEpsilonAtMax="
                        + (sizing == null ? "NaN"
                        : String.format(Locale.ROOT, "%.6f",
                        sizing.epsilonAtMaxSamples))
                        + ", valid=" + valid
                        + ", bestSoFar=" + better);
                if (!better) continue;
                bestWeights = weights;
                bestAnchorIndex = anchorIndex;
                bestAnchorWeight = anchorWeight;
                bestComponentCover = componentCover;
                bestPreferredCover = preferredCover;
                bestEffectiveBound = effectiveBound;
                bestAnchorLogRange = anchorLogRange;
                bestStats = stats;
                bestClip = candidateClip;
                bestSizing = sizing;
                bestEpsilon = sizing.epsilonAtMaxSamples;
                bestFinalSamples = sizing.finalSamples;
                bestChiSquare = chiSquare;
            }
        }
        if (bestWeights == null || bestStats == null
                || bestSizing == null || bestAnchorIndex < 0) return null;
        System.out.println("[PACK*-conditional-eta-v4] mixture winner:"
                + " anchorIndex=" + bestAnchorIndex
                + ", anchorSource="
                + components[bestAnchorIndex].sourceProposalId
                + ", anchorWeight="
                + String.format(Locale.ROOT, "%.6f", bestAnchorWeight)
                + ", weights=" + Arrays.toString(bestWeights)
                + ", certificateMode="
                + (bestComponentCover
                ? "component-cover" : "single-anchor")
                + ", effectiveB_kcal="
                + String.format(Locale.ROOT, "%.6f", bestEffectiveBound)
                + ", chiSquareProxy="
                + String.format(Locale.ROOT, "%.6f", bestChiSquare)
                + ", predictedFinalN=" + bestFinalSamples
                + ", predictedEpsilonAtMax="
                + String.format(Locale.ROOT, "%.6f", bestEpsilon));
        return new ConditionalEtaV4MixtureDesign(
                components, bestWeights, bestAnchorIndex,
                etaV4ProofBoundKcal, bestEffectiveBound,
                bestAnchorLogRange, bestComponentCover, bestStats,
                bestClip, bestSizing, candidateDpSweeps);
    }

    private int conditionalEtaV4CoverViolations(
            List<CCDResult> samples,
            ConditionalProposalSnapshot[] components,
            double[] weights, double boundKcal) {
        if (samples == null || components == null || weights == null
                || components.length == 0
                || components.length != weights.length
                || !Double.isFinite(boundKcal) || boundKcal < 0.0) {
            return Integer.MAX_VALUE;
        }
        double[] bounds = new double[components.length];
        Arrays.fill(bounds, boundKcal);
        double[] residuals = new double[components.length];
        int violations = 0;
        for (CCDResult sample : samples) {
            if (sample == null || !Double.isFinite(sample.eTrue)) {
                return Integer.MAX_VALUE;
            }
            for (int k = 0; k < components.length; k++) {
                residuals[k] = computeFullConfPairwiseEnergy(
                        sample.conf, components[k].emat)
                        - sample.eTrue;
            }
            if (!PackStarEtaV4Importance.componentCoverContains(
                    residuals, bounds, weights)) violations++;
        }
        return violations;
    }

    /** Draw IID component labels, then batch exact-DP sampling by component. */
    private List<CCDResult> sampleConditionalEtaV4Mixture(
            ConditionalEtaV4MixtureDesign design,
            int sampleCount, String stage) {
        if (design == null || design.components == null
                || design.weights == null || sampleCount < 1
                || stage == null || stage.isEmpty()) {
            return Collections.emptyList();
        }
        int componentCount = design.components.length;
        List<List<Integer>> positionsByComponent = new ArrayList<>();
        for (int k = 0; k < componentCount; k++) {
            positionsByComponent.add(new ArrayList<>());
        }
        Random labelRandom = stageRandom(stage + "-labels");
        for (int sample = 0; sample < sampleCount; sample++) {
            double draw = labelRandom.nextDouble();
            double cumulative = 0.0;
            int selected = componentCount - 1;
            for (int k = 0; k < componentCount; k++) {
                cumulative += design.weights[k];
                if (draw < cumulative) {
                    selected = k;
                    break;
                }
            }
            positionsByComponent.get(selected).add(sample);
        }

        CCDResult[] ordered = new CCDResult[sampleCount];
        int dpSweeps = 0;
        for (int k = 0; k < componentCount; k++) {
            List<Integer> positions = positionsByComponent.get(k);
            if (positions.isEmpty()) continue;
            ConditionalProposalSnapshot component = design.components[k];
            CorrectedDPResult loaded = recomputeDP(
                    component.emat, component.eta);
            dpSweeps++;
            if (!Double.isFinite(loaded.logZCorrected)
                    || Math.abs(loaded.logZCorrected - component.logZ)
                    > 1.0e-7) {
                System.out.println("[PACK*-conditional-eta-v4] ERROR: "
                        + stage + " component " + k
                        + " frozen/reloaded logZ mismatch="
                        + String.format(Locale.ROOT, "%.9g",
                        loaded.logZCorrected - component.logZ));
                return Collections.emptyList();
            }
            List<int[]> conformations = sampleConformationsFromDP(
                    positions.size(),
                    stageRandom(stage + "-component-" + k));
            List<CCDResult> componentResults = runParallelCCD(
                    conformations, component.emat, component.logZ,
                    stage + "-component-" + k);
            if (componentResults.size() != positions.size()) {
                return Collections.emptyList();
            }
            for (int i = 0; i < positions.size(); i++) {
                ordered[positions.get(i)] = componentResults.get(i);
            }
        }
        List<CCDResult> result = new ArrayList<>(sampleCount);
        for (CCDResult sample : ordered) {
            if (sample == null) return Collections.emptyList();
            result.add(sample);
        }
        StringBuilder counts = new StringBuilder();
        for (int k = 0; k < componentCount; k++) {
            if (k > 0) counts.append(',');
            counts.append(k).append(':')
                    .append(positionsByComponent.get(k).size());
        }
        System.out.println("[PACK*-conditional-eta-v4] " + stage
                + " IID mixture labels: N=" + sampleCount
                + ", counts=" + counts
                + ", dpSweeps=" + dpSweeps);
        return result;
    }

    private boolean auditConditionalEtaV4Support(
            List<CCDResult> samples,
            ConditionalEtaV4MixtureDesign design,
            String stage) {
        if (samples == null || samples.isEmpty() || design == null) {
            return false;
        }
        ConditionalProposalSnapshot anchor =
                design.components[design.anchorIndex];
        int lowerViolations = 0;
        int anchorViolations = 0;
        int allComponentsViolate = 0;
        int mixtureRangeViolations = 0;
        double minGap = Double.POSITIVE_INFINITY;
        double maxAnchorOverCorrection = Double.NEGATIVE_INFINITY;
        double[] maxComponentOverCorrection =
                new double[design.components.length];
        Arrays.fill(maxComponentOverCorrection,
                Double.NEGATIVE_INFINITY);
        double[] componentBounds =
                new double[design.components.length];
        Arrays.fill(componentBounds, design.anchorBoundKcal);
        double[] componentResiduals =
                new double[design.components.length];
        double maxNormalizedWeight = 0.0;
        for (CCDResult sample : samples) {
            if (!Double.isFinite(sample.eTrue)
                    || !Double.isFinite(sample.eMin)) return false;
            double gap = sample.eTrue - sample.eMin;
            double anchorEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, anchor.emat);
            double overCorrection = anchorEnergy - sample.eTrue;
            for (int k = 0; k < design.components.length; k++) {
                double componentEnergy = computeFullConfPairwiseEnergy(
                        sample.conf, design.components[k].emat);
                componentResiduals[k] = componentEnergy - sample.eTrue;
                maxComponentOverCorrection[k] = Math.max(
                        maxComponentOverCorrection[k],
                        componentResiduals[k]);
            }
            minGap = Math.min(minGap, gap);
            maxAnchorOverCorrection = Math.max(
                    maxAnchorOverCorrection, overCorrection);
            if (gap < 0.0) lowerViolations++;
            if (overCorrection > design.anchorBoundKcal) {
                anchorViolations++;
            }
            if (!PackStarEtaV4Importance.componentCoverContains(
                    componentResiduals, componentBounds,
                    design.weights)) {
                allComponentsViolate++;
            }
            double[] logComponents =
                    conditionalEtaV4LogComponentProbabilities(
                            sample.conf, design.components);
            double normalized = logComponents == null ? Double.NaN
                    : PackStarEtaV4Importance
                    .boundedNormalizedMixtureWeight(
                            -sample.eTrue / RT,
                            logComponents, design.weights,
                            design.anchorLogRange);
            if (!Double.isFinite(normalized)) {
                mixtureRangeViolations++;
            } else {
                maxNormalizedWeight = Math.max(
                        maxNormalizedWeight, normalized);
            }
        }
        System.out.println("[PACK*-conditional-eta-v4-audit] " + stage
                + ": samples=" + samples.size()
                + ", anchorIndex=" + design.anchorIndex
                + ", anchorSource=" + anchor.sourceProposalId
                + ", certificateMode="
                + (design.componentCover
                ? "component-cover" : "single-anchor")
                + ", anchorB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                design.anchorBoundKcal)
                + ", maxAnchorOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f",
                maxAnchorOverCorrection)
                + ", maxComponentOverCorrectionKcal="
                + Arrays.toString(maxComponentOverCorrection)
                + ", anchorBoundViolations=" + anchorViolations
                + ", allComponentsViolate="
                + allComponentsViolate
                + ", mixtureRangeViolations=" + mixtureRangeViolations
                + ", maxNormalizedWeight="
                + String.format(Locale.ROOT, "%.9f",
                maxNormalizedWeight)
                + ", min(E_true-E_m)="
                + String.format(Locale.ROOT, "%.9f", minGap)
                + ", lowerBoundViolations=" + lowerViolations);
        maybeWriteEtaV4MomentTailAudit(samples, design, stage);
        int certificateViolations = design.componentCover
                ? allComponentsViolate : anchorViolations;
        return lowerViolations == 0 && certificateViolations == 0
                && mixtureRangeViolations == 0;
    }

    /**
     * Eta-v3 treats the refit/repair trajectory as a trust-region proposal
     * archive instead of assuming the last refit is best.  Every archived
     * proposal has a genuinely on-policy adaptation batch.  Its energy gauge
     * may be shifted downward to preserve fixed proof headroom; this leaves
     * q_eta unchanged exactly.  Eligibility requires on-policy ESS and a
     * predicted certificate at maxEstSamples, then the winner minimizes the
     * predicted final sample count and empirical chi-square proxy.
     */
    private ConditionalDpRichDesign selectConditionalEtaV3Design(
            List<ConditionalProposalSnapshot> proposalArchive,
            List<CCDResult> adaptation) {
        if (proposalArchive == null || proposalArchive.isEmpty()
                || adaptation == null || adaptation.isEmpty()) {
            return null;
        }

        EtaCorrections bestEta = null;
        EnergyMatrix bestEmat = null;
        double bestAlpha = Double.NaN;
        double bestMaxOverCorrection = Double.NaN;
        double bestLogZ = Double.NaN;
        double bestLogRange = Double.NaN;
        double bestGaugeShift = Double.NaN;
        double bestEpsilonAtMax = Double.NaN;
        double bestChiSquare = Double.NaN;
        double bestEss = Double.NaN;
        int bestFinalSamples = -1;
        String bestSource = null;
        int candidateIndex = 0;
        double targetMaxKcal = etaV3ProofBoundKcal
                - etaV3HeadroomKcal;

        for (ConditionalProposalSnapshot snapshot : proposalArchive) {
            candidateIndex++;
            List<CCDResult> onPolicy = conditionalSamplesFromSource(
                    adaptation, snapshot.sourceProposalId);
            double unshiftedMax = maxOverCorrectionKcal(
                    adaptation, snapshot.emat);
            double shiftDown = Double.isFinite(unshiftedMax)
                    ? Math.max(0.0, unshiftedMax - targetMaxKcal)
                    : Double.NaN;
            EtaCorrections candidateEta = Double.isFinite(shiftDown)
                    ? snapshot.eta.shiftedConstant(-shiftDown) : null;
            EnergyMatrix candidateEmat = candidateEta == null
                    ? null : buildCorrectedEmat(candidateEta);
            double candidateMax = Double.isFinite(shiftDown)
                    ? unshiftedMax - shiftDown : Double.NaN;
            double logZ = Double.isFinite(shiftDown)
                    ? snapshot.logZ + shiftDown / RT : Double.NaN;
            double logRange = conditionalLogRangeProxy(
                    logZ, etaV3ProofBoundKcal, RT);
            ConditionalAdaptationStats onPolicyStats =
                    candidateEmat == null ? null
                            : buildConditionalAdaptationStats(
                            onPolicy, candidateEmat, logZ,
                            etaV3ProofBoundKcal);
            boolean reweightingValid = onPolicyStats != null;
            boolean essPass = reweightingValid
                    && conditionalDpRichEssPasses(
                    onPolicyStats.effectiveSampleSize,
                    onPolicy.size(), dpRichMinEssFraction);
            double candidateClip = reweightingValid
                    ? selectConditionalClipLogCap(
                    onPolicyStats, etaV3ProofBoundKcal)
                    : Double.NaN;
            ConditionalSizing candidateSizing = reweightingValid
                    && Double.isFinite(candidateClip)
                    ? computeConditionalSizing(onPolicyStats,
                    etaV3ProofBoundKcal, candidateClip)
                    : null;
            double chiSquare = reweightingValid
                    ? conditionalWeightChiSquareProxy(
                    onPolicyStats.logResidualWeights,
                    onPolicyStats.normalizedProposalWeights)
                    : Double.NaN;
            boolean reachable = candidateSizing != null
                    && candidateSizing.epsilonAtMaxSamples
                    <= targetEpsilon + 1.0e-12;
            boolean eligible = essPass && reachable
                    && Double.isFinite(candidateMax)
                    && candidateMax
                    <= etaV3ProofBoundKcal + 1.0e-12;
            boolean better = isBetterConditionalEtaV3Candidate(
                    eligible,
                    candidateSizing == null ? -1
                            : candidateSizing.finalSamples,
                    chiSquare, logRange,
                    candidateSizing == null ? Double.NaN
                            : candidateSizing.epsilonAtMaxSamples,
                    bestFinalSamples, bestChiSquare,
                    bestLogRange, bestEpsilonAtMax);

            System.out.println("[PACK*-conditional-eta-v3] candidate "
                    + candidateIndex + "/" + proposalArchive.size()
                    + ": source=" + snapshot.sourceProposalId
                    + ", alpha=" + String.format(Locale.ROOT,
                    "%.6f", snapshot.reportedAlpha)
                    + ", onPolicySamples=" + onPolicy.size()
                    + ", unshiftedAdaptationMaxKcal="
                    + String.format(Locale.ROOT, "%.6f",
                    unshiftedMax)
                    + ", gaugeShiftDownKcal="
                    + String.format(Locale.ROOT, "%.6f", shiftDown)
                    + ", adaptationMaxOverCorrectionKcal="
                    + String.format(Locale.ROOT, "%.6f",
                    candidateMax)
                    + ", selectedB_kcal="
                    + String.format(Locale.ROOT, "%.6f",
                    etaV3ProofBoundKcal)
                    + ", headroomKcal="
                    + String.format(Locale.ROOT, "%.6f",
                    etaV3HeadroomKcal)
                    + ", logZ_eta="
                    + String.format(Locale.ROOT, "%.6f", logZ)
                    + ", logRangeProxy="
                    + String.format(Locale.ROOT, "%.6f", logRange)
                    + ", ESS="
                    + (reweightingValid
                    ? String.format(Locale.ROOT, "%.6f",
                    onPolicyStats.effectiveSampleSize) : "NaN")
                    + ", ESSPass=" + essPass
                    + ", chiSquareProxy="
                    + String.format(Locale.ROOT, "%.6f", chiSquare)
                    + ", weightedClipLogCap="
                    + String.format(Locale.ROOT, "%.6f",
                    candidateClip)
                    + ", predictedFinalN="
                    + (candidateSizing == null ? "NA"
                    : Integer.toString(candidateSizing.finalSamples))
                    + ", predictedEpsilonAtMax="
                    + (candidateSizing == null ? "NaN"
                    : String.format(Locale.ROOT, "%.6f",
                    candidateSizing.epsilonAtMaxSamples))
                    + ", reachable=" + reachable
                    + ", eligible=" + eligible
                    + ", bestSoFar=" + better);
            if (!better) continue;

            bestEta = candidateEta;
            bestEmat = candidateEmat;
            bestAlpha = snapshot.reportedAlpha;
            bestMaxOverCorrection = candidateMax;
            bestLogZ = logZ;
            bestLogRange = logRange;
            bestGaugeShift = shiftDown;
            bestFinalSamples = candidateSizing.finalSamples;
            bestEpsilonAtMax = candidateSizing.epsilonAtMaxSamples;
            bestChiSquare = chiSquare;
            bestEss = onPolicyStats.effectiveSampleSize;
            bestSource = "archive:" + snapshot.sourceProposalId;
        }

        if (bestEta == null || bestEmat == null
                || bestFinalSamples < 2) return null;

        CorrectedDPResult winnerDP = recomputeDP(bestEmat, bestEta);
        double analyticLogZ = bestLogZ;
        bestLogZ = winnerDP.logZCorrected;
        bestLogRange = conditionalLogRangeProxy(
                bestLogZ, etaV3ProofBoundKcal, RT);
        if (!Double.isFinite(bestLogRange)) return null;
        double logZGaugeError = Math.abs(bestLogZ - analyticLogZ);
        if (!Double.isFinite(logZGaugeError)
                || logZGaugeError > 1.0e-7) {
            System.out.println("[PACK*-conditional-eta-v3] WARNING:"
                    + " analytic/recomputed gauge logZ mismatch="
                    + String.format(Locale.ROOT, "%.9g",
                    logZGaugeError));
        }

        System.out.println("[PACK*-conditional-eta-v3] winner: source="
                + bestSource
                + ", alpha=" + String.format(Locale.ROOT,
                "%.6f", bestAlpha)
                + ", gaugeShiftDownKcal="
                + String.format(Locale.ROOT, "%.6f",
                bestGaugeShift)
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f",
                etaV3ProofBoundKcal)
                + ", headroomKcal="
                + String.format(Locale.ROOT, "%.6f",
                etaV3HeadroomKcal)
                + ", adaptationMaxOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f",
                bestMaxOverCorrection)
                + ", logZ_eta=" + String.format(Locale.ROOT,
                "%.6f", bestLogZ)
                + ", logRangeProxy=" + String.format(Locale.ROOT,
                "%.6f", bestLogRange)
                + ", predictedFinalN=" + bestFinalSamples
                + ", predictedEpsilonAtMax="
                + String.format(Locale.ROOT, "%.6f",
                bestEpsilonAtMax)
                + ", chiSquareProxy="
                + String.format(Locale.ROOT, "%.6f",
                bestChiSquare)
                + ", screeningESS="
                + String.format(Locale.ROOT, "%.6f", bestEss)
                + ", candidateDpSweeps=1"
                + ", extraCandidateCCD="
                + etaV3CalibrationSamples);
        return new ConditionalDpRichDesign(
                bestEta, bestEmat, bestAlpha,
                etaV3ProofBoundKcal, bestMaxOverCorrection,
                bestLogZ, bestLogRange, bestSource, 1,
                bestGaugeShift, bestFinalSamples,
                bestEpsilonAtMax, bestChiSquare, bestEss);
    }

    /**
     * Screen eta-repair candidates using adaptation CCD only.  Candidate
     * alphas are the raw eta, the current witness-repair alpha, one analytic
     * repair alpha for every pre-registered B, and zero eta as a fail-safe.
     * Every feasible candidate receives an exact DP normalizer.  Since q is
     * fixed, minimizing log Z_eta + B/RT minimizes the deterministic absolute
     * range of the corresponding bounded importance contribution.
     */
    private ConditionalDpRichDesign selectConditionalDpRichDesign(
            EtaCorrections baseEta,
            List<CCDResult> adaptation,
            double currentRepairAlpha) {
        if (baseEta == null || adaptation == null || adaptation.isEmpty()) {
            return null;
        }

        EnergyMatrix baseEmat = buildCorrectedEmat(baseEta);
        LinkedHashMap<Double, String> candidates = new LinkedHashMap<>();
        addConditionalDpRichAlpha(candidates, 1.0, "raw-eta");
        addConditionalDpRichAlpha(candidates, currentRepairAlpha,
                "current-witness-repair");
        for (double proofBound : residualBoundGridKcal) {
            double alpha = chooseEtaRepairAlpha(adaptation, baseEmat,
                    proofBound, null);
            addConditionalDpRichAlpha(candidates, alpha,
                    "repair-for-B=" + String.format(Locale.ROOT,
                            "%.6f", proofBound));
        }
        addConditionalDpRichAlpha(candidates, 0.0, "zero-eta");

        EtaCorrections bestEta = null;
        EnergyMatrix bestEmat = null;
        double bestAlpha = Double.NaN;
        double bestBound = Double.NaN;
        double bestMaxOverCorrection = Double.NaN;
        double bestLogZ = Double.NaN;
        double bestLogRange = Double.NaN;
        String bestSource = null;
        int dpSweeps = 0;
        int candidateIndex = 0;

        for (Map.Entry<Double, String> entry : candidates.entrySet()) {
            candidateIndex++;
            double alpha = entry.getKey();
            EtaCorrections candidateEta = baseEta.scaled(alpha);
            EnergyMatrix candidateEmat = buildCorrectedEmat(candidateEta);
            double maxOverCorrection = maxOverCorrectionKcal(
                    adaptation, candidateEmat);
            double candidateBound = selectSmallestResidualBound(
                    maxOverCorrection, residualBoundGridKcal);
            if (!Double.isFinite(candidateBound)) {
                System.out.println("[PACK*-conditional-dp-rich] candidate "
                        + candidateIndex + "/" + candidates.size()
                        + ": source=" + entry.getValue()
                        + ", alpha=" + String.format(Locale.ROOT,
                        "%.6f", alpha)
                        + ", adaptationMaxOverCorrectionKcal="
                        + String.format(Locale.ROOT, "%.6f",
                        maxOverCorrection)
                        + ", rejected=exceeds-B-grid");
                continue;
            }

            CorrectedDPResult candidateDP = recomputeDP(
                    candidateEmat, candidateEta);
            dpSweeps++;
            double logZ = candidateDP.logZCorrected;
            double logRange = conditionalLogRangeProxy(
                    logZ, candidateBound, RT);
            ConditionalAdaptationStats reweighted =
                    buildConditionalAdaptationStats(
                            adaptation, candidateEmat, logZ,
                            candidateBound);
            boolean reweightingValid = reweighted != null;
            boolean essPass = reweightingValid
                    && conditionalDpRichEssPasses(
                    reweighted.effectiveSampleSize,
                    adaptation.size(), dpRichMinEssFraction);
            double candidateClip = reweightingValid
                    ? selectConditionalClipLogCap(
                    reweighted, candidateBound)
                    : Double.NaN;
            ConditionalSizing candidateSizing =
                    reweightingValid && Double.isFinite(candidateClip)
                            ? computeConditionalSizing(
                            reweighted, candidateBound, candidateClip)
                            : null;
            boolean better = isBetterConditionalDpRichCandidate(
                    logRange, candidateBound, alpha,
                    bestLogRange, bestBound, bestAlpha);
            System.out.println("[PACK*-conditional-dp-rich] candidate "
                    + candidateIndex + "/" + candidates.size()
                    + ": source=" + entry.getValue()
                    + ", alpha=" + String.format(Locale.ROOT,
                    "%.6f", alpha)
                    + ", adaptationMaxOverCorrectionKcal="
                    + String.format(Locale.ROOT, "%.6f",
                    maxOverCorrection)
                    + ", selectedB_kcal="
                    + String.format(Locale.ROOT, "%.6f", candidateBound)
                    + ", logZ_eta="
                    + String.format(Locale.ROOT, "%.6f", logZ)
                    + ", logRangeProxy="
                    + String.format(Locale.ROOT, "%.6f", logRange)
                    + ", reweighting="
                    + (reweightingValid ? "valid" : "invalid")
                    + ", sourceProposals="
                    + (reweightingValid
                    ? Integer.toString(reweighted.sourceProposalCount)
                    : "NA")
                    + ", ESS="
                    + (reweightingValid
                    ? String.format(Locale.ROOT, "%.6f",
                    reweighted.effectiveSampleSize)
                    : "NaN")
                    + ", ESSPass=" + essPass
                    + ", weightedClipLogCap="
                    + String.format(Locale.ROOT, "%.6f", candidateClip)
                    + ", predictedFinalN="
                    + (candidateSizing != null
                    ? Integer.toString(candidateSizing.finalSamples)
                    : "NA")
                    + ", predictedEpsilonAtMax="
                    + (candidateSizing != null
                    ? String.format(Locale.ROOT, "%.6f",
                    candidateSizing.epsilonAtMaxSamples)
                    : "NaN")
                    + ", bestSoFar=" + better);
            if (!better) continue;
            bestEta = candidateEta;
            bestEmat = candidateEmat;
            bestAlpha = alpha;
            bestBound = candidateBound;
            bestMaxOverCorrection = maxOverCorrection;
            bestLogZ = logZ;
            bestLogRange = logRange;
            bestSource = entry.getValue();
        }

        if (bestEta == null || bestEmat == null
                || !Double.isFinite(bestLogRange)) {
            return null;
        }

        // Candidate screening leaves the final candidate's tables active, not
        // necessarily the winner's. Reload the winner before monitor/final.
        CorrectedDPResult winnerDP = recomputeDP(bestEmat, bestEta);
        dpSweeps++;
        bestLogZ = winnerDP.logZCorrected;
        bestLogRange = conditionalLogRangeProxy(
                bestLogZ, bestBound, RT);
        if (!Double.isFinite(bestLogRange)) return null;

        System.out.println("[PACK*-conditional-dp-rich] winner: source="
                + bestSource
                + ", alpha=" + String.format(Locale.ROOT, "%.6f", bestAlpha)
                + ", selectedB_kcal="
                + String.format(Locale.ROOT, "%.6f", bestBound)
                + ", adaptationMaxOverCorrectionKcal="
                + String.format(Locale.ROOT, "%.6f",
                bestMaxOverCorrection)
                + ", logZ_eta="
                + String.format(Locale.ROOT, "%.6f", bestLogZ)
                + ", logRangeProxy="
                + String.format(Locale.ROOT, "%.6f", bestLogRange)
                + ", candidateDpSweeps=" + dpSweeps
                + ", extraCandidateCCD=0");
        return new ConditionalDpRichDesign(
                bestEta, bestEmat, bestAlpha, bestBound,
                bestMaxOverCorrection, bestLogZ, bestLogRange,
                bestSource, dpSweeps);
    }

    private static void addConditionalDpRichAlpha(
            LinkedHashMap<Double, String> candidates,
            double alpha, String source) {
        if (!Double.isFinite(alpha) || alpha < 0.0 || alpha > 1.0) return;
        double canonical = alpha == 0.0 ? 0.0 : alpha == 1.0 ? 1.0 : alpha;
        candidates.putIfAbsent(canonical, source);
    }

    /**
     * Choose a scalar eta shrink from adaptation data only.  This is an
     * efficiency repair, not a support certificate.  In production, unseen
     * conformations remain covered only conditionally on the frozen global B
     * premise; monitor/final observations can falsify but not prove it.
     */
    private double chooseEtaRepairAlpha(List<CCDResult> adaptation,
                                        EnergyMatrix candidateEmat) {
        return chooseEtaRepairAlpha(adaptation, candidateEmat,
                etaRepairTargetKcal, "scalar shrink selection");
    }

    private double chooseEtaRepairAlpha(List<CCDResult> adaptation,
                                        EnergyMatrix candidateEmat,
                                        double targetOverCorrectionKcal,
                                        String logContext) {
        double alpha = 1.0;
        int constraints = 0;
        for (CCDResult sample : adaptation) {
            double gap = sample.eTrue - sample.eMin;
            if (!Double.isFinite(gap) || gap < 0.0) {
                return Double.NaN;
            }
            double etaValue = computeFullConfPairwiseEnergy(
                    sample.conf, candidateEmat) - sample.eMin;
            if (!Double.isFinite(etaValue) || etaValue <= 0.0) continue;
            double overCorrection = etaValue - gap;
            if (overCorrection <= targetOverCorrectionKcal) continue;
            double admissible = (gap + targetOverCorrectionKcal) / etaValue;
            alpha = Math.min(alpha, clamp(admissible, 0.0, 1.0));
            constraints++;
        }
        if (alpha < 1.0) {
            alpha = Math.max(etaRepairMinAlpha,
                    clamp(alpha * etaRepairSafety, 0.0, 1.0));
        }
        if (logContext != null) {
            System.out.println("[PACK*-eta-repair] " + logContext
                    + ": constraints=" + constraints
                    + ", targetOverCorrection="
                    + String.format(Locale.ROOT, "%.6f",
                    targetOverCorrectionKcal)
                    + " kcal/mol, safety="
                    + String.format(Locale.ROOT, "%.4f", etaRepairSafety)
                    + ", alpha="
                    + String.format(Locale.ROOT, "%.6f", alpha));
        }
        return alpha;
    }

    /** A valid upper bound on unbiased sample variance for N values in [0,1]. */
    static double maximumUnitIntervalSampleVariance(int sampleCount) {
        if (sampleCount <= 1) return Double.POSITIVE_INFINITY;
        return sampleCount / (4.0 * (sampleCount - 1.0));
    }

    static final class RadiusBalancedThreshold {
        final double logCap;
        final double logRadius;

        RadiusBalancedThreshold(double logCap, double logRadius) {
            this.logCap = logCap;
            this.logRadius = logRadius;
        }
    }

    /**
     * Minimize the deterministic confidence-radius proxy
     *
     *   q_eta C Delta_b + q_lambda C^{-lambda} Delta_t,  C=exp(logCap)>=1.
     *
     * For lambda>0 the unconstrained stationary point is available in closed
     * form; lambda=0 has a C-independent tail radius and therefore chooses
     * logCap=0.  This is an efficiency choice only and never enters coverage.
     */
    static RadiusBalancedThreshold balanceLogCap(
            double logZEta, double logZLambda, double lambda,
            double bulkDelta, double tailDelta) {
        if (!Double.isFinite(logZEta) || !Double.isFinite(logZLambda)
                || !Double.isFinite(lambda) || lambda < 0.0
                || !Double.isFinite(bulkDelta) || !(bulkDelta > 0.0)
                || !Double.isFinite(tailDelta) || !(tailDelta > 0.0)) {
            return null;
        }
        double logA = logZEta + Math.log(bulkDelta);
        double logD = logZLambda + Math.log(tailDelta);
        double logCap = 0.0;
        if (lambda > 0.0) {
            logCap = Math.max(0.0,
                    (Math.log(lambda) + logD - logA) / (lambda + 1.0));
        }
        double logBulkRadius = logA + logCap;
        double logTailRadius = logD - lambda * logCap;
        double logRadius = logAddExp(logBulkRadius, logTailRadius);
        if (!Double.isFinite(logCap) || !Double.isFinite(logRadius)) {
            return null;
        }
        return new RadiusBalancedThreshold(logCap, logRadius);
    }

    private TiltProposal selectTiltProposal(
            EtaCorrections eta, double logZCorrected,
            int nBulk, int nTail) {
        double deltaPer = perEventDelta(true);
        double bulkDelta = solveBernsteinDelta(
                nBulk, maximumUnitIntervalSampleVariance(nBulk),
                1.0, deltaPer);
        double tailDelta = solveBernsteinDelta(
                nTail, maximumUnitIntervalSampleVariance(nTail),
                1.0, deltaPer);
        if (!Double.isFinite(bulkDelta) || !(bulkDelta > 0.0)
                || !Double.isFinite(tailDelta) || !(tailDelta > 0.0)) {
            return null;
        }

        double bestLambda = 0.0;
        double bestLogZ = logZMinDet;
        RadiusBalancedThreshold bestDesign = balanceLogCap(
                logZCorrected, bestLogZ, bestLambda,
                bulkDelta, tailDelta);
        if (bestDesign == null) return null;
        double bestLogCap = bestDesign.logCap;
        double bestLogRange = bestLogZ;
        double bestLogRadius = bestDesign.logRadius;
        EnergyMatrix bestEmat = branchMinimizingEmat;
        double loadedLambda = Double.NaN;
        int evaluatedPositive = 0;
        int dpSweeps = 0;

        System.out.println("[PACK*-tail] fixed-budget design: N_bulk="
                + nBulk + ", N_tail=" + nTail
                + ", Delta_bulk(worst)="
                + String.format(Locale.ROOT, "%.8f", bulkDelta)
                + ", Delta_tail(worst)="
                + String.format(Locale.ROOT, "%.8f", tailDelta));
        System.out.println("[PACK*-tail] candidate lambda=0.000000, logZ="
                + String.format(Locale.ROOT, "%.8f", logZMinDet)
                + ", B_kcal=0.000000"
                + ", logRange="
                + String.format(Locale.ROOT, "%.8f", logZMinDet)
                + ", logRadiusProxy="
                + String.format(Locale.ROOT, "%.8f", bestLogRadius)
                + " (known initial q_m normalizer)");
        for (double lambda : tailTiltLambdas) {
            if (!(lambda > 0.0)) continue;
            EnergyMatrix tiltedEmat = buildEtaAdjustedEmat(eta, -lambda);
            CorrectedDPResult tiltedDP = recomputeDP(tiltedEmat, eta);
            loadedLambda = lambda;
            evaluatedPositive++;
            dpSweeps++;
            double logZ = tiltedDP.logZCorrected;
            RadiusBalancedThreshold design = balanceLogCap(
                    logZCorrected, logZ, lambda,
                    bulkDelta, tailDelta);
            if (design == null) continue;
            double logRange = logZ - lambda * design.logCap;
            System.out.println("[PACK*-tail] candidate lambda="
                    + String.format(Locale.ROOT, "%.6f", lambda)
                    + ", logZ=" + String.format(Locale.ROOT, "%.8f", logZ)
                    + ", B_kcal="
                    + String.format(Locale.ROOT, "%.6f", RT * design.logCap)
                    + ", logRange="
                    + String.format(Locale.ROOT, "%.8f", logRange)
                    + ", logRadiusProxy="
                    + String.format(Locale.ROOT, "%.8f", design.logRadius));
            if (Double.isFinite(logRange)
                    && design.logRadius < bestLogRadius) {
                bestLambda = lambda;
                bestLogZ = logZ;
                bestLogRange = logRange;
                bestLogCap = design.logCap;
                bestLogRadius = design.logRadius;
                bestEmat = tiltedEmat;
            }
        }

        if (!Double.isFinite(bestLogRange)) return null;
        if (Double.compare(loadedLambda, bestLambda) != 0) {
            CorrectedDPResult selectedDP = recomputeDP(bestEmat, eta);
            dpSweeps++;
            bestLogZ = selectedDP.logZCorrected;
            RadiusBalancedThreshold selectedDesign = balanceLogCap(
                    logZCorrected, bestLogZ, bestLambda,
                    bulkDelta, tailDelta);
            if (selectedDesign == null) return null;
            bestLogCap = selectedDesign.logCap;
            bestLogRadius = selectedDesign.logRadius;
            bestLogRange = bestLogZ - bestLambda * bestLogCap;
        }
        System.out.println("[PACK*-tail] selected lambda="
                + String.format(Locale.ROOT, "%.6f", bestLambda)
                + ", B_kcal="
                + String.format(Locale.ROOT, "%.6f", RT * bestLogCap)
                + ", logZ_lambda="
                + String.format(Locale.ROOT, "%.8f", bestLogZ)
                + ", logR_lambda="
                + String.format(Locale.ROOT, "%.8f", bestLogRange)
                + ", logRadiusProxy="
                + String.format(Locale.ROOT, "%.8f", bestLogRadius)
                + ", positiveCandidates=" + evaluatedPositive
                + ", tailDpSweeps=" + dpSweeps);
        return new TiltProposal(bestLambda, bestLogZ, bestLogRange,
                bestLogCap, bestLogRadius,
                bestEmat, evaluatedPositive, dpSweeps);
    }

    static double clippedBulkNormalized(double gap, double eta,
                                        double rt, double logCap) {
        if (!Double.isFinite(gap) || !Double.isFinite(eta)
                || !Double.isFinite(rt) || !(rt > 0.0)
                || !Double.isFinite(logCap)) return Double.NaN;
        double logWeight = (eta - gap) / rt;
        return Math.exp(Math.min(logWeight, logCap) - logCap);
    }

    /** h = max(0, exp(-g/RT) - C exp(-eta/RT)), evaluated stably. */
    static double tailRemainderH(double gap, double eta,
                                 double rt, double logCap) {
        if (!Double.isFinite(gap) || !Double.isFinite(eta)
                || !Double.isFinite(rt) || !(rt > 0.0)
                || !Double.isFinite(logCap)) return Double.NaN;
        double logA = -gap / rt;
        double logB = logCap - eta / rt;
        if (!(logA > logB)) return 0.0;
        double h = Math.exp(logA) * (-Math.expm1(logB - logA));
        if (h < 0.0 && h > -1.0e-15) h = 0.0;
        if (h > 1.0 && h < 1.0 + 1.0e-12) h = 1.0;
        return h;
    }

    /** Tail sample normalized by R_lambda, hence deterministically in [0,1]. */
    static double tiltedTailNormalized(double gap, double eta,
                                       double rt, double logCap,
                                       double lambda) {
        if (!Double.isFinite(lambda) || lambda < 0.0) return Double.NaN;
        double h = tailRemainderH(gap, eta, rt, logCap);
        if (!(h > 0.0)) return h;
        double etaAboveThreshold = eta / rt - logCap;
        if (etaAboveThreshold < -1.0e-12) return Double.NaN;
        etaAboveThreshold = Math.max(0.0, etaAboveThreshold);
        double value = h * Math.exp(-lambda * etaAboveThreshold);
        if (value < 0.0 && value > -1.0e-15) value = 0.0;
        if (value > 1.0 && value < 1.0 + 1.0e-12) value = 1.0;
        return value;
    }

    private UnitIntervalCI unitIntervalCI(double[] values, double deltaTarget) {
        if (values.length <= 1) return null;
        double sum = 0.0;
        double sum2 = 0.0;
        for (double raw : values) {
            if (!Double.isFinite(raw) || raw < -1.0e-12
                    || raw > 1.0 + 1.0e-12) return null;
            double value = clamp(raw, 0.0, 1.0);
            sum += value;
            sum2 += value * value;
        }
        double mean = sum / values.length;
        double variance = Math.max(0.0,
                (sum2 - values.length * mean * mean) / (values.length - 1));
        double bound = solveBernsteinDelta(
                values.length, variance, 1.0, deltaTarget);
        if (!Double.isFinite(bound)) return null;
        return new UnitIntervalCI(mean, variance, bound,
                Math.max(0.0, mean - bound),
                Math.min(1.0, mean + bound));
    }

    private void computeUnconditionalBulkTailBound(
            List<CCDResult> bulkSamples,
            List<CCDResult> tailSamples,
            EnergyMatrix correctedEmat,
            double logZCorrected,
            double logCap,
            TiltProposal tilt) {

        double[] bulkValues = new double[bulkSamples.size()];
        double[] bulkXi = new double[bulkSamples.size()];
        double[] bulkLogW = new double[bulkSamples.size()];
        double sumResidual = 0.0;
        double sumResidual2 = 0.0;
        for (int i = 0; i < bulkSamples.size(); i++) {
            CCDResult sample = bulkSamples.get(i);
            double gap = sample.eTrue - sample.eMin;
            double etaValue = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat) - sample.eMin;
            double xi = gap - etaValue;
            bulkXi[i] = xi;
            bulkLogW[i] = -xi / RT;
            bulkValues[i] = clippedBulkNormalized(
                    gap, etaValue, RT, logCap);
            sumResidual += xi;
            sumResidual2 += xi * xi;
        }

        double[] tailValues = new double[tailSamples.size()];
        int positiveTailSamples = 0;
        double maxTailValue = 0.0;
        for (int i = 0; i < tailSamples.size(); i++) {
            CCDResult sample = tailSamples.get(i);
            double gap = sample.eTrue - sample.eMin;
            double etaValue = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat) - sample.eMin;
            double value = tiltedTailNormalized(
                    gap, etaValue, RT, logCap, tilt.lambda);
            tailValues[i] = value;
            if (value > 0.0) positiveTailSamples++;
            if (Double.isFinite(value)) maxTailValue = Math.max(maxTailValue, value);
        }

        double deltaPer = perEventDelta(true);
        UnitIntervalCI bulkCI = unitIntervalCI(bulkValues, deltaPer);
        UnitIntervalCI tailCI = unitIntervalCI(tailValues, deltaPer);
        if (bulkCI == null || tailCI == null) {
            failCertificate("unconditional: non-finite or out-of-range bulk/tail statistic", -1L);
            return;
        }

        double logBulkScale = logZCorrected + logCap;
        double logTailScale = tilt.logRange;
        double logBulkLower = scaledMeanLog(logBulkScale, bulkCI.lower);
        double logBulkUpper = scaledMeanLog(logBulkScale, bulkCI.upper);
        double logTailLower = scaledMeanLog(logTailScale, tailCI.lower);
        double logTailUpper = scaledMeanLog(logTailScale, tailCI.upper);
        double candidateLower = logAddExp(logBulkLower, logTailLower);
        double samplingUpper = logAddExp(logBulkUpper, logTailUpper);
        double candidateUpper = Math.min(samplingUpper, logZMinDet);

        if (Double.isNaN(candidateLower) || !Double.isFinite(candidateUpper)
                || candidateLower > candidateUpper + 1.0e-12) {
            failCertificate("unconditional: combined log interval is invalid", -1L);
            return;
        }

        logZLowerPAC = candidateLower;
        logZUpperPAC = candidateUpper;
        zLower = bigExpFromLog(logZLowerPAC);
        zUpper = bigExpFromLog(logZUpperPAC);
        epsilon = epsilonFromLogBounds(logZLowerPAC, logZUpperPAC);
        certificateValid = isValidCertificate(zLower, zUpper, epsilon);
        certificateFailureReason = certificateValid
                ? "" : "combined bulk-tail interval failed final validation";
        if (!certificateValid) {
            System.out.println("[PACK*-unconditional] ERROR: " + certificateFailureReason);
            return;
        }

        int nBulk = bulkSamples.size();
        meanResidual = sumResidual / nBulk;
        stdResidual = nBulk > 1
                ? Math.sqrt(Math.max(0.0,
                (sumResidual2 - nBulk * meanResidual * meanResidual) / (nBulk - 1)))
                : 0.0;
        meanPsi = shiftedMeanToDouble(bulkCI.mean, logCap);
        varPsi = bulkCI.variance;
        cvPsi = bulkCI.mean > 0.0
                ? Math.sqrt(bulkCI.variance) / bulkCI.mean
                : Double.MAX_VALUE;

        emitResidualSampleTraces(
                bulkSamples, bulkXi, bulkLogW, logZCorrected, logCap);
        maybeDumpResidualSamples(
                bulkSamples, bulkXi, logZCorrected, meanResidual, stdResidual);

        boolean qmBinds = logZMinDet < samplingUpper;
        System.out.println("[PACK*-unconditional] bulk CI: mean="
                + String.format(Locale.ROOT, "%.8f", bulkCI.mean)
                + ", var=" + String.format(Locale.ROOT, "%.8f", bulkCI.variance)
                + ", Delta=" + String.format(Locale.ROOT, "%.8f", bulkCI.delta)
                + ", normalized=["
                + String.format(Locale.ROOT, "%.8f", bulkCI.lower) + ","
                + String.format(Locale.ROOT, "%.8f", bulkCI.upper) + "]"
                + ", N=" + bulkValues.length);
        System.out.println("[PACK*-unconditional] tail CI: lambda="
                + String.format(Locale.ROOT, "%.6f", tilt.lambda)
                + ", logR=" + String.format(Locale.ROOT, "%.8f", tilt.logRange)
                + ", mean=" + String.format(Locale.ROOT, "%.8f", tailCI.mean)
                + ", var=" + String.format(Locale.ROOT, "%.8f", tailCI.variance)
                + ", Delta=" + String.format(Locale.ROOT, "%.8f", tailCI.delta)
                + ", normalized=["
                + String.format(Locale.ROOT, "%.8f", tailCI.lower) + ","
                + String.format(Locale.ROOT, "%.8f", tailCI.upper) + "]"
                + ", positive=" + positiveTailSamples + "/" + tailValues.length
                + ", maxY=" + String.format(Locale.ROOT, "%.8f", maxTailValue));
        System.out.println("[PACK*-unconditional] combined certificate: upperVia="
                + (qmBinds ? "q_m(det)" : "bulk+tilted-tail")
                + ", log10Lower=" + formatLog10(logZLowerPAC)
                + ", log10Upper=" + formatLog10(logZUpperPAC)
                + ", epsilon=" + String.format(Locale.ROOT, "%.8f", epsilon)
                + ", certified=true");
    }

    private static double scaledMeanLog(double logScale, double mean) {
        if (mean <= 0.0) return Double.NEGATIVE_INFINITY;
        if (!Double.isFinite(logScale) || !Double.isFinite(mean)) return Double.NaN;
        return logScale + Math.log(mean);
    }

    /** Predicted/final relative width for a clipped-bulk plus tilted-tail sum. */
    static double combinedBulkTailEpsilon(
            double logBulkScale,
            double bulkLower,
            double bulkUpper,
            double logTailScale,
            double tailMean,
            double tailDelta,
            double deterministicLogUpper) {
        if (!Double.isFinite(logBulkScale)
                || !Double.isFinite(logTailScale)
                || !Double.isFinite(deterministicLogUpper)
                || !Double.isFinite(bulkLower)
                || !Double.isFinite(bulkUpper)
                || !Double.isFinite(tailMean)
                || !Double.isFinite(tailDelta)
                || bulkLower < 0.0 || bulkUpper < bulkLower
                || tailMean < 0.0 || tailMean > 1.0
                || tailDelta < 0.0) return 1.0;
        double tailLower = Math.max(0.0, tailMean - tailDelta);
        double tailUpper = Math.min(1.0, tailMean + tailDelta);
        double candidateLower = logAddExp(
                scaledMeanLog(logBulkScale, bulkLower),
                scaledMeanLog(logTailScale, tailLower));
        double samplingUpper = logAddExp(
                scaledMeanLog(logBulkScale, bulkUpper),
                scaledMeanLog(logTailScale, tailUpper));
        double candidateUpper = Math.min(samplingUpper, deterministicLogUpper);
        return epsilonFromLogBounds(candidateLower, candidateUpper);
    }

    static double logAddExp(double first, double second) {
        if (first == Double.NEGATIVE_INFINITY) return second;
        if (second == Double.NEGATIVE_INFINITY) return first;
        if (!Double.isFinite(first) || !Double.isFinite(second)) {
            return first == Double.POSITIVE_INFINITY || second == Double.POSITIVE_INFINITY
                    ? Double.POSITIVE_INFINITY : Double.NaN;
        }
        double max = Math.max(first, second);
        return max + Math.log1p(Math.exp(Math.min(first, second) - max));
    }

    static boolean isValidCertificate(BigDecimal lower, BigDecimal upper,
                                      double epsilon) {
        return lower != null && upper != null
                && MathTools.isFinite(lower) && MathTools.isFinite(upper)
                && lower.signum() >= 0 && upper.signum() > 0
                && lower.compareTo(upper) <= 0
                && Double.isFinite(epsilon)
                && epsilon >= 0.0 && epsilon <= 1.0;
    }

    private void failCertificate(String reason, long startTime) {
        certificateValid = false;
        certificateFailureReason = reason;
        setZeroBounds(reason);
        if (startTime >= 0L) printFinalSummary(startTime);
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
        boolean clip = getConfigBoolean("packstar.pac.clip", DEFAULT_CLIP);
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
            double clipQuantile = getConfigDouble("packstar.pac.clipQuantile",
                    DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE);
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
        double deltaPer = perEventDelta(clip); // 6 plain / 12 clipped (see perEventDelta)

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
            System.out.println("[PACK*-2stage] N* calc(clip): pilot=" + n
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
        System.out.println("[PACK*-2stage] N* calc: pilot=" + n
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
        double sizeSafety = getConfigDouble("packstar.pac.sizeSafety", DEFAULT_SIZE_SAFETY);
        double sizeTarget = targetEpsilon * sizeSafety;
        int unreachableCap = getConfigInteger("packstar.pac.unreachableCap",
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
            System.out.println("[PACK*-2stage] WARNING: packstar.pac.residualBound is smaller than "
                    + "the observed residual-weight range; " + context
                    + " uses the sample range and is empirical, not a strict PAC certificate.");
            return sampleRange;
        }
        System.out.println("[PACK*-2stage] WARNING: no finite packstar.pac.residualBound is configured; "
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

    private void emitResidualSampleTraces(List<CCDResult> samples,
                                          double[] xi,
                                          double[] logW,
                                          double logZCorrected,
                                          double logCap) {
        if (sampleListener == null) {
            return;
        }
        for (int i = 0; i < samples.size(); i++) {
            CCDResult r = samples.get(i);
            double eProposal = r.eTrue - xi[i];
            boolean clipped = Double.isFinite(logCap) && logW[i] > logCap;
            sampleListener.onSample(new PackStarSampleTrace(
                    PackStarSampleTrace.Stage.ESTIMATION,
                    i,
                    r.conf,
                    r.eTrue,
                    r.eMin,
                    eProposal,
                    xi[i],
                    logW[i],
                    logZCorrected,
                    logZMinDet,
                    logCap,
                    clipped));
        }
    }

    private boolean shouldWriteMomentTailAudit(String stage) {
        String assignments = getConfigProperty(
                PAC_MOMENT_TAIL_ASSIGNMENTS_PROPERTY, null);
        String summary = getConfigProperty(
                PAC_MOMENT_TAIL_SUMMARY_PROPERTY, null);
        if (assignments == null || assignments.trim().isEmpty()
                || summary == null || summary.trim().isEmpty()) {
            return false;
        }
        String selectedStage = getConfigProperty(
                PAC_MOMENT_TAIL_STAGE_PROPERTY, "monitor");
        return selectedStage != null
                && selectedStage.trim().equals(stage);
    }

    private double[] momentTailThresholds() {
        return parseFiniteIncreasingGrid(getConfigProperty(
                        PAC_MOMENT_TAIL_THRESHOLDS_PROPERTY,
                        "-4,-3,-2,-1,0,0.25,0.5,0.75,1,1.25,1.5,"
                                + "1.75,2,2.25,2.5,2.75,3,3.5,4,5,6,8,12"),
                "PACK* moment/tail thresholds");
    }

    private void maybeWritePairwiseMomentTailAudit(
            List<CCDResult> samples, EnergyMatrix correctedEmat,
            String stage) {
        if (!shouldWriteMomentTailAudit(stage)) return;
        String scheme = etaV3 ? "eta-v3-on-policy"
                : dpRich ? "dp-rich-v2" : "eta-v1";
        List<PackStarMomentTailAudit.Sample> audit =
                new ArrayList<>(samples.size());
        for (CCDResult sample : samples) {
            double eProposal = computeFullConfPairwiseEnergy(
                    sample.conf, correctedEmat);
            double overCorrection = eProposal - sample.eTrue;
            audit.add(new PackStarMomentTailAudit.Sample(
                    sample.conf, sample.eMin, sample.eTrue,
                    overCorrection / RT, overCorrection));
        }
        writeMomentTailAudit(audit, scheme, stage);
    }

    private void maybeWriteEtaV4MomentTailAudit(
            List<CCDResult> samples,
            ConditionalEtaV4MixtureDesign design,
            String stage) {
        if (!shouldWriteMomentTailAudit(stage)) return;
        ConditionalProposalSnapshot anchor =
                design.components[design.anchorIndex];
        double logZReference = anchor.logZ;
        List<PackStarMomentTailAudit.Sample> audit =
                new ArrayList<>(samples.size());
        for (CCDResult sample : samples) {
            double[] logComponents =
                    conditionalEtaV4LogComponentProbabilities(
                            sample.conf, design.components);
            double logQMix = logComponents == null ? Double.NaN
                    : PackStarEtaV4Importance.logMixtureProbability(
                            logComponents, design.weights);
            double logWeight = -sample.eTrue / RT
                    - logQMix - logZReference;
            double anchorEnergy = computeFullConfPairwiseEnergy(
                    sample.conf, anchor.emat);
            audit.add(new PackStarMomentTailAudit.Sample(
                    sample.conf, sample.eMin, sample.eTrue,
                    logWeight, anchorEnergy - sample.eTrue));
        }
        writeMomentTailAudit(audit, "eta-v4.2-balance-mixture", stage);
    }

    private void writeMomentTailAudit(
            List<PackStarMomentTailAudit.Sample> samples,
            String scheme, String stage) {
        try {
            PackStarMomentTailAudit.Result result =
                    PackStarMomentTailAudit.analyze(
                            samples, momentTailThresholds(), RT);
            File assignments = new File(getConfigProperty(
                    PAC_MOMENT_TAIL_ASSIGNMENTS_PROPERTY, "").trim());
            File summary = new File(getConfigProperty(
                    PAC_MOMENT_TAIL_SUMMARY_PROPERTY, "").trim());
            PackStarMomentTailAudit.writeAssignmentTsv(
                    assignments, result, scheme, stage);
            PackStarMomentTailAudit.writeSummaryTsv(
                    summary, result, scheme, stage);
            System.out.println("[PACK*-moment-tail] " + scheme
                    + "/" + stage + ": N=" + samples.size()
                    + ", normalizedM2="
                    + String.format(Locale.ROOT, "%.9f",
                    result.normalizedSecondMoment)
                    + ", ESSFraction="
                    + String.format(Locale.ROOT, "%.9f",
                    result.effectiveSampleSizeFraction)
                    + ", maxSecondMomentShare="
                    + String.format(Locale.ROOT, "%.9f",
                    result.maxSecondMomentShare)
                    + ", summary=" + summary);
        } catch (Exception ex) {
            System.err.println("[PACK*-moment-tail] diagnostic output failed: "
                    + ex.getMessage());
        }
    }

    // ===== Diagnostic dump for multi-body residual attribution =====
    // Enabled by -Dpackstar.pac.dumpDir=<dir>. Writes, per pfunc, the
    // per-sample residual xi = E_true - E_eta (pairwise part already removed) plus
    // the branch-decomposition bags, so offline analysis can locate the 3-body
    // interactions that the pairwise eta cannot absorb.
    private static final AtomicInteger PAC_DUMP_COUNTER = new AtomicInteger(0);

    private void maybeDumpResidualSamples(List<CCDResult> est, double[] xi,
                                          double logZCorrected,
                                          double meanXi, double stdXi) {
        String dir = getConfigProperty("packstar.pac.dumpDir", null);
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
            System.out.println("[PACK*-2stage] residual dump -> " + rf.getName()
                    + " (N=" + est.size() + ", logZcorr=" + String.format("%.4f", logZCorrected)
                    + ", meanXi=" + String.format("%.4f", meanXi) + ")");
        } catch (Exception ex) {
            System.err.println("[PACK*-2stage] residual dump failed: " + ex.getMessage());
        }
    }

    private static final AtomicInteger PAC_TRAIN_COUNTER = new AtomicInteger(0);

    private void maybeDumpTrainingSamples(List<CCDResult> train, EtaCorrections eta,
                                          EnergyMatrix correctedEmat) {
        String dir = getConfigProperty("packstar.pac.dumpDir", null);
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
            System.out.println("[PACK*-2stage] training dump -> " + tf.getName()
                    + " (trainN=" + train.size() + ")");
        } catch (Exception ex) {
            System.err.println("[PACK*-2stage] training dump failed: " + ex.getMessage());
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
        if (getConfigBoolean("packstar.pac.clip", DEFAULT_CLIP)) {
            double clipQuantile = getConfigDouble("packstar.pac.clipQuantile",
                    DEFAULT_LEGACY_CONDITIONAL_CLIP_QUANTILE);
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

        emitResidualSampleTraces(estResults, xi, logW, logZCorrected, Double.NaN);
        maybeDumpResidualSamples(estResults, xi, logZCorrected, meanResidual, stdResidual);

        double deltaPer = perEventDelta(false); // non-clip: 2 tails x 3 pfuncs = 6 events
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

        System.out.println("[PACK*-2stage] Bernstein(residual): Delta=" + String.format("%.6f", boundDelta)
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

        emitResidualSampleTraces(estResults, xi, logW, logZCorrected, logCap);
        maybeDumpResidualSamples(estResults, xi, logZCorrected, meanResidual, stdResidual);

        double deltaPer = perEventDelta(true); // clip: lower+upper+Hoeffding+excess = 4 x 3 = 12 events
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
        int estimationResidualBoundViolations = 0;
        double estimationMaxOverCorrectionKcal = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < N; i++) {
            double u = Math.max(0.0, Math.exp(logW[i] - logCap) - 1.0);
            sumU += u; sumU2 += u * u;
            if (logW[i] > logCap) {
                double overCorrectionKcal = -xi[i];
                estimationMaxOverCorrectionKcal = Math.max(
                        estimationMaxOverCorrectionKcal,
                        overCorrectionKcal
                );
                if (Double.isFinite(residualBoundKcal) && overCorrectionKcal > residualBoundKcal) {
                    estimationResidualBoundViolations++;
                }
            }
        }
        double meanU = sumU / N;
        double varU = N > 1 ? Math.max(0.0, (sumU2 - N * meanU * meanU) / (N - 1)) : 0.0;
        int residualBoundViolations = pilotResidualBoundViolations + estimationResidualBoundViolations;
        double maxOverCorrectionKcal = Math.max(
                pilotMaxOverCorrectionKcal,
                estimationMaxOverCorrectionKcal
        );
        boolean residualBoundViolated = Double.isFinite(residualBoundKcal)
                && residualBoundViolations > 0;

        double biasUpper;
        boolean biasCertified;
        String biasMode;
        if (residualBoundViolated) {
            // A sample from the fixed corrected proposal directly falsifies the assumed
            // one-sided support bound. A sample maximum cannot certify unseen tail mass,
            // so fail closed and let the deterministic q_m ceiling supply the upper bound.
            biasUpper = Double.POSITIVE_INFINITY;
            biasCertified = false;
            biasMode = "assumptionViolated";
            System.out.println("[PACK*-2stage] WARNING: one-sided PAC residual bound violated by "
                    + residualBoundViolations + " pilot/estimation sample(s): max(E_eta-E_true)="
                    + String.format("%.6f", maxOverCorrectionKcal) + " kcal/mol > configured "
                    + String.format("%.6f", residualBoundKcal)
                    + "; conditional sampling upper disabled, falling back to q_m.");
        } else if (Double.isFinite(residualBoundKcal) && residualBoundKcal / RT - logCap > 0) {
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

        System.out.println("[PACK*-2stage] Bernstein(clipped): Delta=" + String.format("%.6f", boundDelta)
                + ", clipQ=" + String.format("%.3f", clipQuantile)
                + ", logCap=" + String.format("%.4f", logCap)
                + ", nClipped=" + nClipped + "/" + N
                + ", meanWc=" + String.format("%.6f", mean)
                + ", varWc=" + String.format("%.6f", var)
                + ", range=1.0(det)"
                + ", biasUpper=" + String.format("%.6f", biasUpper) + "(" + biasMode + ")"
                + (biasCertified ? "" : "[tail uncertified -> upper via q_m only]")
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
                        daemonThreadFactory("packstar-sample"))
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
                    throw new IllegalStateException("PACK* sampling reached edge before M position "
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
        if (!gpuSampling || activeProposalEmat.hasHigherOrderTerms()
                || !edge.canUseGpuSampling()) {
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
        double localEnergy = computeLocalEnergy(edge, mRCs, lambdaRCs, activeProposalEmat);
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
            throw new RuntimeException("PACK* sampling interrupted", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new RuntimeException("PACK* sampling failed", cause);
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

        System.out.println("[PACK*] Phase 1 edge: lambdaStates=" + totalLambda
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
                double localEnergy = computeLocalEnergy(edge, mRCs, lambdaRCs, activeProposalEmat);
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
        return edge.computeLocalEnergyPublic(
                mRCs, lambdaRCs, emat, interactionGraph);
    }

    // ========== Phase 2: Parallel CCD minimization ==========

    private static class CCDResult {
        final int[] conf;
        final double eTrue;
        final double eMin;
        final EnergyCalculator.EnergiedParametricMolecule epmol;
        final double sourceProposalEnergy;
        final double sourceProposalLogZ;
        final String sourceProposalId;

        CCDResult(int[] conf, double eTrue, double eMin,
                  EnergyCalculator.EnergiedParametricMolecule epmol) {
            this(conf, eTrue, eMin, epmol,
                    Double.NaN, Double.NaN, null);
        }

        CCDResult(int[] conf, double eTrue, double eMin,
                  EnergyCalculator.EnergiedParametricMolecule epmol,
                  double sourceProposalEnergy,
                  double sourceProposalLogZ,
                  String sourceProposalId) {
            this.conf = conf;
            this.eTrue = eTrue;
            this.eMin = eMin;
            this.epmol = epmol;
            this.sourceProposalEnergy = sourceProposalEnergy;
            this.sourceProposalLogZ = sourceProposalLogZ;
            this.sourceProposalId = sourceProposalId;
        }
    }

    /**
     * Run CCD minimizations in parallel via minimizingEcalc.tasks.
     */
    private List<CCDResult> runParallelCCD(List<int[]> conformations) {
        return runParallelCCD(conformations, null, Double.NaN, null);
    }

    /**
     * Run CCD while retaining the exact proposal provenance of this batch.
     * The source energy and log normalizer are later used only by adaptation
     * to reweight samples to a DP-rich winner; monitor/final remain fresh.
     */
    private List<CCDResult> runParallelCCD(
            List<int[]> conformations,
            EnergyMatrix sourceProposalEmat,
            double sourceProposalLogZ,
            String sourceProposalId) {
        boolean recordSourceProposal = sourceProposalEmat != null
                && Double.isFinite(sourceProposalLogZ)
                && sourceProposalId != null
                && !sourceProposalId.isEmpty();
        int total = conformations.size();
        AtomicReferenceArray<CCDResult> results =
                new AtomicReferenceArray<>(total);
        AtomicInteger completed = new AtomicInteger(0);

        // ---- parallelism instrumentation ----
        long batchStartNs = System.nanoTime();
        long[] ccdWallUs = new long[total];
        java.util.Set<String> workerThreads = java.util.concurrent.ConcurrentHashMap.newKeySet();
        AtomicInteger inFlight = new AtomicInteger(0);
        AtomicInteger peakConcurrency = new AtomicInteger(0);
        // --------------------------------------

        for (int idx = 0; idx < total; idx++) {
            int[] conf = conformations.get(idx);
            final int ci = idx;
            final long submitNs = System.nanoTime();
            ResidueInteractions inters = makeSparseFullConfInters(conf);
            RCTuple tuple = new RCTuple(conf);

            minimizingEcalc.calcEnergyAsync(tuple, inters, epmol -> {
                int cur = inFlight.incrementAndGet();
                peakConcurrency.accumulateAndGet(cur, Math::max);
                workerThreads.add(Thread.currentThread().getName());

                double eTrue = epmol.energy;
                double eMin = computeFullConfPairwiseEnergy(conf, branchMinimizingEmat);
                double sourceEnergy = recordSourceProposal
                        ? computeFullConfProposalEnergy(
                        conf, sourceProposalEmat)
                        : Double.NaN;
                // Callback completion order depends on task scheduling and GPU
                // load.  Preserve the seeded sample order so eta accumulation,
                // pilot statistics, and all later floating-point reductions do
                // not inherit that nondeterminism.
                results.set(ci, new CCDResult(
                        conf, eTrue, eMin, epmol,
                        sourceEnergy,
                        recordSourceProposal
                                ? sourceProposalLogZ : Double.NaN,
                        recordSourceProposal ? sourceProposalId : null));
                ccdWallUs[ci] = (System.nanoTime() - submitNs) / 1000L;

                inFlight.decrementAndGet();
                int done = completed.incrementAndGet();
                if (done % 100 == 0 || done == total) {
                    System.out.println("[PACK*] CCD progress: " + done + "/" + total);
                }
            });
        }

        minimizingEcalc.tasks.waitForFinish();

        // ---- parallelism summary ----
        long batchWallMs = (System.nanoTime() - batchStartNs) / 1_000_000L;
        long[] sorted = java.util.Arrays.copyOf(ccdWallUs, total);
        java.util.Arrays.sort(sorted);
        long sumUs = 0L;
        for (long w : ccdWallUs) sumUs += w;
        double p50 = total > 0 ? sorted[Math.min(total - 1, total / 2)] / 1000.0 : 0.0;
        double p95 = total > 0 ? sorted[Math.min(total - 1, (int) (total * 0.95))] / 1000.0 : 0.0;
        double maxMs = total > 0 ? sorted[total - 1] / 1000.0 : 0.0;
        double sumMs = sumUs / 1000.0;
        double effParallel = batchWallMs > 0 ? sumMs / batchWallMs : 0.0;
        System.out.printf("[PACK*-INSTR] CCD batch: n=%d batchWall=%dms distinctThreads=%d peakConcurrency=%d "
                        + "perCCD(ms) p50=%.1f p95=%.1f max=%.1f sumCCD=%.0fms effParallel=%.1fx%n",
                total, batchWallMs, workerThreads.size(), peakConcurrency.get(),
                p50, p95, maxMs, sumMs, effParallel);
        // -----------------------------

        List<CCDResult> orderedResults = new ArrayList<>(total);
        for (int i = 0; i < total; i++) {
            CCDResult result = results.get(i);
            if (result != null) orderedResults.add(result);
        }
        return orderedResults;
    }

    private double computeFullConfProposalEnergy(
            int[] conf, EnergyMatrix proposalEmat) {
        if (!proposalEmat.hasHigherOrderTerms()) {
            return computeFullConfPairwiseEnergy(conf, proposalEmat);
        }
        return proposalEmat.getInternalEnergy(new RCTuple(conf));
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
        // Added once per complete conformation.  This is an energy gauge: it
        // changes Z_eta and residual weights by reciprocal constants but leaves
        // the normalized proposal q_eta exactly unchanged.
        double globalOffsetKcal = 0.0;
        PackStarTripleEtaCorrections tripleEta = null;
        double tripleEtaScale = 0.0;

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

        EtaCorrections scaled(double alpha) {
            if (!Double.isFinite(alpha) || alpha < 0.0 || alpha > 1.0) {
                throw new IllegalArgumentException(
                        "eta scale must be finite and in [0,1]: " + alpha);
            }
            if (tripleEta != null && alpha != 1.0) {
                throw new IllegalStateException(
                        "scale pair eta before attaching the eta-relative"
                                + " triple residual");
            }
            int[] numRCs = new int[oneBody.length];
            for (int pos = 0; pos < oneBody.length; pos++) {
                numRCs[pos] = oneBody[pos].length;
            }
            EtaCorrections scaled = new EtaCorrections(oneBody.length, numRCs);
            for (int pos = 0; pos < oneBody.length; pos++) {
                for (int rc = 0; rc < oneBody[pos].length; rc++) {
                    scaled.oneBody[pos][rc] = alpha * oneBody[pos][rc];
                    scaled.oneBodyCounts[pos][rc] = oneBodyCounts[pos][rc];
                }
            }
            for (Map.Entry<Long, double[]> entry : pairSums.entrySet()) {
                double[] value = entry.getValue();
                scaled.pairSums.put(entry.getKey(),
                        new double[]{alpha * value[0], value[1]});
            }
            scaled.oneBodyCount = oneBodyCount;
            scaled.pairCount = pairCount;
            scaled.globalOffsetKcal = alpha * globalOffsetKcal;
            scaled.tripleEta = tripleEta;
            scaled.tripleEtaScale = tripleEtaScale;
            return scaled;
        }

        EtaCorrections withTripleEta(
                PackStarTripleEtaCorrections prior, double scale) {
            if (prior == null) {
                throw new IllegalArgumentException(
                        "triple eta prior must be non-null");
            }
            if (!Double.isFinite(scale) || scale != 1.0) {
                throw new IllegalArgumentException(
                        "triple eta scale is not tunable and must be 1: "
                                + scale);
            }
            EtaCorrections combined = scaled(1.0);
            combined.tripleEta = prior;
            combined.tripleEtaScale = scale;
            return combined;
        }

        EtaCorrections shiftedConstant(double deltaKcal) {
            if (!Double.isFinite(deltaKcal)) {
                throw new IllegalArgumentException(
                        "eta global energy shift must be finite: " + deltaKcal);
            }
            EtaCorrections shifted = scaled(1.0);
            shifted.globalOffsetKcal += deltaKcal;
            return shifted;
        }
    }

    private EtaV4SparseFeatures buildEtaV4SparseFeatures(
            List<CCDResult> samples) {
        if (samples == null || samples.isEmpty()) return null;
        LinkedHashMap<Long, Integer> indices = new LinkedHashMap<>();
        List<EtaV4FeatureCell> cells = new ArrayList<>();
        int[][] sampleIndices = new int[samples.size()][];

        for (int sampleIndex = 0;
             sampleIndex < samples.size(); sampleIndex++) {
            int[] conf = samples.get(sampleIndex).conf;
            if (conf == null || conf.length != rcs.getNumPos()) return null;
            List<Integer> active = new ArrayList<>();
            for (int pos = 0; pos < conf.length; pos++) {
                int rc = conf[pos];
                if (rc < 0) return null;
                long key = Long.MIN_VALUE
                        | ((long) pos << 32)
                        | (rc & 0xffffffffL);
                Integer index = indices.get(key);
                if (index == null) {
                    index = cells.size();
                    indices.put(key, index);
                    cells.add(new EtaV4FeatureCell(pos, rc));
                }
                active.add(index);
            }
            for (int pos1 = 0; pos1 < conf.length; pos1++) {
                for (int pos2 = pos1 + 1;
                     pos2 < conf.length; pos2++) {
                    if (!interactionGraph.hasEdge(pos1, pos2)) continue;
                    int rc1 = conf[pos1];
                    int rc2 = conf[pos2];
                    long key = EtaCorrections.packPairKey(
                            pos1, rc1, pos2, rc2);
                    Integer index = indices.get(key);
                    if (index == null) {
                        index = cells.size();
                        indices.put(key, index);
                        cells.add(new EtaV4FeatureCell(
                                pos1, rc1, pos2, rc2));
                    }
                    active.add(index);
                }
            }
            sampleIndices[sampleIndex] = new int[active.size()];
            for (int i = 0; i < active.size(); i++) {
                sampleIndices[sampleIndex][i] = active.get(i);
            }
        }
        return new EtaV4SparseFeatures(
                cells.toArray(new EtaV4FeatureCell[0]),
                sampleIndices);
    }

    /**
     * Convert a natural-parameter descent direction into an energy update.
     * Both every observed conformation's total change and every individual
     * cell change are capped by the pre-registered trust radius.
     */
    private EtaV4ProposalUpdate applyEtaV4EnergyGradient(
            EtaCorrections baseEta, EtaV4SparseFeatures features,
            double[] gradient, double trustKcal) {
        if (baseEta == null || features == null || gradient == null
                || gradient.length != features.cells.length
                || !Double.isFinite(trustKcal) || !(trustKcal > 0.0)) {
            return null;
        }
        double maxCell = 0.0;
        for (double value : gradient) {
            if (!Double.isFinite(value)) return null;
            maxCell = Math.max(maxCell, Math.abs(value));
        }
        double maxObserved = 0.0;
        for (int[] active : features.sampleIndices) {
            double total = 0.0;
            for (int index : active) total += gradient[index];
            if (!Double.isFinite(total)) return null;
            maxObserved = Math.max(maxObserved, Math.abs(total));
        }
        double rawScale = Math.max(maxCell, maxObserved);
        if (!(rawScale > 0.0) || !Double.isFinite(rawScale)) return null;
        double scaleKcal = trustKcal / rawScale;

        EtaCorrections updated = baseEta.scaled(1.0);
        for (int index = 0; index < features.cells.length; index++) {
            double deltaKcal = scaleKcal * gradient[index];
            EtaV4FeatureCell cell = features.cells[index];
            if (cell.pair) {
                long key = EtaCorrections.packPairKey(
                        cell.pos1, cell.rc1, cell.pos2, cell.rc2);
                double[] stored = updated.pairSums.get(key);
                double count = stored == null ? 1.0 : stored[1];
                double value = updated.getPairEta(
                        cell.pos1, cell.rc1, cell.pos2, cell.rc2)
                        + deltaKcal;
                updated.pairSums.put(key,
                        new double[]{value * count, count});
                if (stored == null) updated.pairCount++;
            } else {
                int oldCount = updated.oneBodyCounts[
                        cell.pos1][cell.rc1];
                int count = Math.max(1, oldCount);
                double value = updated.getOneBodyEta(
                        cell.pos1, cell.rc1) + deltaKcal;
                updated.oneBody[cell.pos1][cell.rc1] = value * count;
                updated.oneBodyCounts[cell.pos1][cell.rc1] = count;
                if (oldCount == 0) updated.oneBodyCount++;
            }
        }
        return new EtaV4ProposalUpdate(updated, scaleKcal,
                scaleKcal * maxObserved,
                scaleKcal * maxCell);
    }

    /**
     * No-eta ablation baseline (packstar.pac.etaEnabled=false): an EtaCorrections
     * with no samples ever added. getOneBodyEta/getPairEta both return 0.0 for
     * every (pos,rc) with a zero count, so this is exactly eta === 0, and
     * buildCorrectedEmat(eta) below reduces to E_eta == E_m (SI Lemma 2 remark).
     */
    private EtaCorrections zeroEtaCorrections() {
        int numPos = rcs.getNumPos();
        int[] numRCs = new int[numPos];
        for (int p = 0; p < numPos; p++) {
            numRCs[p] = branchMinimizingEmat.getNumConfAtPos(p);
        }
        return new EtaCorrections(numPos, numRCs);
    }

    /**
     * Extract the full and two complementary fold sufficient statistics in a
     * single forcefield-breakdown pass.  Distinct full assignments are counted
     * separately from raw IID multiplicities so duplicate proposal draws do
     * not masquerade as new correction contexts.
     */
    private FrequencySeverityEtaTraining extractFrequencySeverityEtaTraining(
            List<CCDResult> ccdResults) {
        return extractFrequencySeverityEtaTraining(
                Collections.emptyList(), ccdResults);
    }

    /**
     * Online two-fold fit.  Historical adaptation samples are fixed before the
     * newest on-policy batch exists, so both fold fits may use all history.
     * Each newest-batch observation is used by exactly one fold fit and scored
     * by the complementary fit.  This prevents a failed q_eta batch from being
     * evaluated as if it were an IID extension of the original q_m batch.
     */
    private FrequencySeverityEtaTraining extractFrequencySeverityEtaTraining(
            List<CCDResult> history,
            List<CCDResult> newestOnPolicyBatch) {
        int numPos = rcs.getNumPos();
        int[] numRCs = new int[numPos];
        for (int pos = 0; pos < numPos; pos++) {
            numRCs[pos] = branchMinimizingEmat.getNumConfAtPos(pos);
        }
        EtaCorrections all = new EtaCorrections(numPos, numRCs);
        EtaCorrections[] folds = new EtaCorrections[frequencySeverityFolds];
        FrequencySeverityEtaCoverage[] foldCoverage = new FrequencySeverityEtaCoverage[frequencySeverityFolds];
        Set<String>[] foldAssignments = new Set[frequencySeverityFolds];
        List<CCDResult>[] foldFitSamples = new List[frequencySeverityFolds];
        for (int fold = 0; fold < frequencySeverityFolds; fold++) {
            folds[fold] = new EtaCorrections(numPos, numRCs);
            foldCoverage[fold] = new FrequencySeverityEtaCoverage(numPos, numRCs);
            foldAssignments[fold] = new HashSet<>();
            foldFitSamples[fold] = new ArrayList<>();
        }
        List<CCDResult> allFitSamples = new ArrayList<>(
                history.size() + newestOnPolicyBatch.size());
        FrequencySeverityEtaCoverage allCoverage = new FrequencySeverityEtaCoverage(numPos, numRCs);
        Set<String> allAssignments = new HashSet<>();
        int totalSamples = history.size() + newestOnPolicyBatch.size();
        int[] foldFitSampleCounts = new int[frequencySeverityFolds];
        List<FrequencySeverityEtaCellObservation> observations =
                new ArrayList<>(totalSamples * Math.max(1, numPos));

        for (int sampleIndex = 0; sampleIndex < totalSamples; sampleIndex++) {
            boolean historical = sampleIndex < history.size();
            int newestIndex = sampleIndex - history.size();
            CCDResult result = historical
                    ? history.get(sampleIndex)
                    : newestOnPolicyBatch.get(newestIndex);
            int newestFold = historical ? -1
                    : newestIndex % frequencySeverityFolds;
            boolean[] includeInFoldFit = new boolean[frequencySeverityFolds];
            allFitSamples.add(result);
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                includeInFoldFit[fold] = historical || fold == newestFold;
                if (includeInFoldFit[fold]) {
                    foldFitSampleCounts[fold]++;
                    foldFitSamples[fold].add(result);
                }
            }
            int[] conf = result.conf;
            ResidueForcefieldEnergy efunc =
                    (ResidueForcefieldEnergy) minimizingEcalc.ecalc
                            .makeEnergyFunction(result.epmol);

            for (int pos = 0; pos < numPos; pos++) {
                int rc = conf[pos];
                if (rc < 0) continue;
                ResidueInteractions singleInters =
                        minimizingEcalc.makeSingleInters(pos, rc);
                double correction = efunc.makeSubset(singleInters).getEnergy()
                        - branchMinimizingEmat.getOneBody(pos, rc);
                all.addOneBodySample(pos, rc, correction);
                for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                    if (includeInFoldFit[fold]) {
                        folds[fold].addOneBodySample(pos, rc, correction);
                    }
                }
                observations.add(new FrequencySeverityEtaCellObservation(
                        sampleIndex, newestFold, false,
                        pos, rc, -1, -1, correction));
            }
            for (int pos1 = 0; pos1 < numPos; pos1++) {
                int rc1 = conf[pos1];
                if (rc1 < 0) continue;
                for (int pos2 = pos1 + 1; pos2 < numPos; pos2++) {
                    int rc2 = conf[pos2];
                    if (rc2 < 0 || !interactionGraph.hasEdge(pos1, pos2)) {
                        continue;
                    }
                    ResidueInteractions pairInters =
                            minimizingEcalc.makePairInters(
                                    pos1, rc1, pos2, rc2);
                    double correction = efunc.makeSubset(pairInters).getEnergy()
                            - branchMinimizingEmat.getPairwise(
                            pos1, rc1, pos2, rc2);
                    all.addPairSample(pos1, rc1, pos2, rc2, correction);
                    for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                        if (includeInFoldFit[fold]) {
                            folds[fold].addPairSample(
                                    pos1, rc1, pos2, rc2, correction);
                        }
                    }
                    observations.add(new FrequencySeverityEtaCellObservation(
                            sampleIndex, newestFold, true,
                            pos1, rc1, pos2, rc2, correction));
                }
            }

            String assignment = Arrays.toString(conf);
            if (allAssignments.add(assignment)) {
                addFrequencySeverityCoverage(conf, allCoverage);
            }
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                if (includeInFoldFit[fold]
                        && foldAssignments[fold].add(assignment)) {
                    addFrequencySeverityCoverage(conf, foldCoverage[fold]);
                }
            }
        }
        finalizeEtaTermCounts(all);
        for (EtaCorrections fold : folds) finalizeEtaTermCounts(fold);
        return new FrequencySeverityEtaTraining(
                all, folds, allCoverage, foldCoverage, observations,
                totalSamples, foldFitSampleCounts, history.size(),
                List.copyOf(allFitSamples), foldFitSamples);
    }

    private void addFrequencySeverityCoverage(int[] conf, FrequencySeverityEtaCoverage coverage) {
        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            if (rc >= 0) coverage.unaryDistinctContexts[pos][rc]++;
        }
        for (int pos1 = 0; pos1 < conf.length; pos1++) {
            int rc1 = conf[pos1];
            if (rc1 < 0) continue;
            for (int pos2 = pos1 + 1; pos2 < conf.length; pos2++) {
                int rc2 = conf[pos2];
                if (rc2 < 0 || !interactionGraph.hasEdge(pos1, pos2)) {
                    continue;
                }
                long key = EtaCorrections.packPairKey(
                        pos1, rc1, pos2, rc2);
                coverage.pairDistinctContexts.merge(key, 1, Integer::sum);
            }
        }
    }

    private void finalizeEtaTermCounts(EtaCorrections eta) {
        eta.oneBodyCount = 0;
        for (int pos = 0; pos < eta.oneBody.length; pos++) {
            for (int rc = 0; rc < eta.oneBody[pos].length; rc++) {
                if (eta.oneBodyCounts[pos][rc] > 0) eta.oneBodyCount++;
            }
        }
        eta.pairCount = eta.pairSums.size();
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
        return buildEtaAdjustedEmat(eta, 1.0);
    }

    /** Build E_factor = E_m + factor * eta. */
    private EnergyMatrix buildEtaAdjustedEmat(EtaCorrections eta, double factor) {
        if (!Double.isFinite(factor)) {
            throw new IllegalArgumentException("eta energy factor must be finite");
        }
        ensureProposalInteractionGraph(eta);
        EnergyMatrix corrected = new EnergyMatrix(branchMinimizingEmat);
        int numPos = rcs.getNumPos();

        // Add eta corrections to one-body terms
        for (int pos = 0; pos < numPos; pos++) {
            int numRC = corrected.getNumConfAtPos(pos);
            for (int rc = 0; rc < numRC; rc++) {
                double original = corrected.getOneBody(pos, rc);
                double correction = factor * eta.getOneBodyEta(pos, rc);
                if (pos == 0) {
                    correction += factor * eta.globalOffsetKcal;
                }
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
                        double correction = factor * eta.getPairEta(pos1, rc1, pos2, rc2);
                        corrected.setPairwise(pos1, rc1, pos2, rc2, original + correction);
                    }
                }
            }
        }

        if (eta.tripleEta != null && eta.tripleEtaScale != 0.0) {
            eta.tripleEta.applyResidualTo(
                    corrected,
                    factor * eta.tripleEtaScale,
                    eta::getPairEta);
        }

        zeroProposalFillPairEnergies(corrected);

        return corrected;
    }

    /**
     * Make every selected triple a clique in a proposal-only primal graph.
     * Added edges carry exactly zero pair energy; they change only the branch
     * decomposition needed to evaluate the higher-order factors exactly.
     */
    private void ensureProposalInteractionGraph(EtaCorrections eta) {
        List<int[]> fillEdges = eta != null && eta.tripleEta != null
                && eta.tripleEtaScale != 0.0
                ? eta.tripleEta.requiredFillEdges(interactionGraph)
                : List.of();
        List<int[]> desiredEdges = new ArrayList<>(
                interactionGraph.getEdgeList().size()
                        + fillEdges.size());
        desiredEdges.addAll(interactionGraph.getEdgeList());
        desiredEdges.addAll(fillEdges);
        InteractionGraph desired = fillEdges.isEmpty()
                ? interactionGraph
                : InteractionGraph.buildFromEdges(
                interactionGraph.getNumPositions(), desiredEdges);
        String signature = interactionGraphSignature(desired);
        if (signature.equals(proposalInteractionGraphSignature)) return;

        if (fillEdges.isEmpty()) {
            rootedRoot = initialRootedRoot;
            rootedRootEdge = initialRootedRootEdge;
            proposalInteractionGraph = interactionGraph;
            proposalInteractionGraphSignature = signature;
            System.out.println("[PACK*-triple-eta-dp] restored base proposal"
                    + " graph: edges=" + interactionGraph.getNumEdges());
            return;
        }

        int[] stateCounts = new int[rcs.getNumPos()];
        for (int pos = 0; pos < stateCounts.length; pos++) {
            stateCounts[pos] = rcs.getNum(pos);
        }
        BranchDecomposition decomposition = new BranchDecomposition(
                desired, BranchDecomposition.Strategy.WEIGHTED_HICKS,
                stateCounts);
        decomposition.compute();
        RootedTreeNode rebuiltRoot = decomposition.rootBranchTree(rcs);
        if (rebuiltRoot == null) {
            throw new IllegalStateException(
                    "selected triple proposal graph produced no branch tree");
        }
        RootedTreeEdge.postOrderCompLlambda(rebuiltRoot, true);
        RootedTreeEdge rebuiltRootEdge =
                rebuiltRoot.getLeftChild().getChildOfEdge();
        rebuiltRootEdge.compactTree();
        rootedRoot = rebuiltRoot;
        rootedRootEdge = rebuiltRootEdge;
        proposalInteractionGraph = desired;
        proposalInteractionGraphSignature = signature;
        System.out.println("[PACK*-triple-eta-dp] rebuilt exact proposal"
                + " graph: baseEdges=" + interactionGraph.getNumEdges()
                + ", fillEdges=" + fillEdges.size()
                + ", proposalEdges=" + desired.getNumEdges()
                + ", branchwidth=" + decomposition.getBranchwidth()
                + ", triples="
                + (eta.tripleEta == null ? "[]"
                : eta.tripleEta.positionTriples()));
    }

    private void zeroProposalFillPairEnergies(EnergyMatrix emat) {
        if (proposalInteractionGraph == interactionGraph) return;
        for (int[] edge : proposalInteractionGraph.getEdgeList()) {
            int pos1 = edge[0];
            int pos2 = edge[1];
            if (interactionGraph.hasEdge(pos1, pos2)) continue;
            for (int rc1 : rcs.get(pos1)) {
                for (int rc2 : rcs.get(pos2)) {
                    emat.setPairwise(pos1, rc1, pos2, rc2, 0.0);
                }
            }
        }
    }

    private static String interactionGraphSignature(
            InteractionGraph graph) {
        List<Long> edges = new ArrayList<>();
        for (int[] edge : graph.getEdgeList()) {
            int lower = Math.min(edge[0], edge[1]);
            int upper = Math.max(edge[0], edge[1]);
            edges.add(((long) lower << 32)
                    | (upper & 0xffffffffL));
        }
        Collections.sort(edges);
        return graph.getNumPositions() + ":" + edges;
    }

    /**
     * Recompute DP with corrected emat.
     * The corrected emat replaces the minimizing emat (upper bound side).
     * Returns logZ_corrected.
     */
    private CorrectedDPResult recomputeDP(EnergyMatrix correctedEmat, EtaCorrections eta) {
        // Re-initialize DP with corrected emat as the "min" side
        // Keep rigid emat the same (lower bound)
        EnergyMatrix proposalRigidEmat =
                new EnergyMatrix(branchRigidEmat);
        zeroProposalFillPairEnergies(proposalRigidEmat);
        RootedTreeEdge.postOrderInitIncremental(rootedRoot,
                proposalRigidEmat, correctedEmat,
                proposalInteractionGraph, RT);
        DPCacheStats stats = computeCorrectedDPTables(eta);
        activeProposalEmat = correctedEmat;

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

        sb.append("global=");
        appendDoubleBits(sb, eta.globalOffsetKcal);
        sb.append("|triple=");
        if (eta.tripleEta == null || eta.tripleEtaScale == 0.0) {
            sb.append("none");
        } else {
            sb.append(eta.tripleEta.signatureSha256).append('@');
            appendDoubleBits(sb, eta.tripleEtaScale);
        }
        sb.append("|ob=");
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
     * Unit-tested in TestPackStarLogSpaceBounds.
     */
    public static double epsilonFromLogBounds(double logZLower, double logZUpper) {
        if (Double.isFinite(logZLower) && Double.isFinite(logZUpper)) {
            double orderingError = logZLower - logZUpper;
            if (orderingError > 1.0e-12) return 1.0;
            double eps = 1.0 - Math.exp(logZLower - logZUpper);
            if (eps < 0.0 && eps > -1e-12) eps = 0.0;
            if (eps > 1.0 && eps < 1.0 + 1e-12) eps = 1.0;
            return Double.isFinite(eps) && eps >= 0.0 && eps <= 1.0
                    ? eps : 1.0;
        }
        return 1.0;
    }

    /**
     * Per-event confidence for the union bound across the two interval endpoints
     * of all three K* partition functions.
     * Non-clipped path: 2 tails x 3 pfuncs = 6 events.
     * Clipped path: lower + upper clipped-mass Bernstein + Hoeffding tail
     * probability + empirical Bernstein on the observed clip-excess
     * = 4 events x 3 pfuncs = 12 events.
     */
    private double perEventDelta(boolean clip) {
        return clip ? delta / 12.0 : delta / 6.0;
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
     * passes deltaTarget = perEventDelta(clip): the non-clipped path uses delta/6
     * (2 tails x 3 pfuncs = 6 events); the clipped path uses delta/12 (lower +
     * upper clipped-mass Bernstein + Hoeffding tail + excess Bernstein = 4 x 3 = 12).
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
    static BigDecimal bigExpFromLog(double logVal) {
        if (logVal == Double.NEGATIVE_INFINITY) return BigDecimal.ZERO;
        if (!Double.isFinite(logVal)) return MathTools.BigNaN;

        // Do not route negative log values through BoltzmannCalculator.calc().
        // That legacy path uses ExpFunction.exp(), whose MathContext precision
        // is incorrectly interpreted as a fixed decimal scale and therefore
        // rounds values below roughly 1e-64 to zero.  The precise exp() path
        // retains MathContext significant digits and BigDecimal's arbitrary
        // decimal exponent.
        return LOG_SPACE_EXP.exp(logVal);
    }

    private void setZeroBounds(String reason) {
        System.out.println("[PACK*] WARNING: " + reason + "; returning zero Z bounds");
        certificateValid = false;
        certificateFailureReason = reason;
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
    public boolean hasValidCertificate() { return certificateValid; }
    public String getCertificateFailureReason() { return certificateFailureReason; }
}
