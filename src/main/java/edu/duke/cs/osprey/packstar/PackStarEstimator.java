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
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.function.Consumer;
import java.util.function.IntConsumer;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

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
 *   The sole production path selects a count/context-regularized eta by
 *   cross-fitting one q_m training batch with alpha and triple-residual gamma
 *   both fixed at one, repairs proposal support with bounded on-policy refits,
 *   and draws independent pilot, monitor, and final samples.  Its tail upper
 *   bound is conditional on a frozen external severity premise.
 */
public class PackStarEstimator {

    private static final long MISSING_PAIR_KEY = Long.MIN_VALUE;

    // Configuration
    private static final String PAC_SAMPLES_PROPERTY = "packstar.pac.samples";
    private static final String PAC_CONFIDENCE_PROPERTY = "packstar.pac.confidence";
    private static final String PAC_TARGET_EPSILON_PROPERTY = "packstar.pac.targetEpsilon";
    private static final String PAC_RANDOM_SEED_PROPERTY = "packstar.pac.randomSeed";
    private static final String PAC_SAMPLING_PARALLEL_PROPERTY = "packstar.pac.sampling.parallel";
    private static final String PAC_SAMPLING_THREADS_PROPERTY = "packstar.pac.sampling.threads";
    private static final String PAC_SAMPLING_LARGE_LAMBDA_PROPERTY = "packstar.pac.sampling.largeLambdaThreshold";
    private static final String PAC_SAMPLING_PROGRESS_PROPERTY = "packstar.pac.sampling.progress";
    private static final String PAC_SAMPLING_GPU_PROPERTY = "packstar.pac.sampling.gpu";
    private static final String PAC_CCD_SUBMISSION_BATCH_SIZE_PROPERTY =
            "packstar.pac.ccd.submissionBatchSize";
    private static final String PAC_CCD_INSTRUMENTATION_PROPERTY =
            "packstar.pac.ccd.instrumentation";
    // Train/adaptation/final sizing for the adaptive frequency/severity route.
    private static final String PAC_TRAIN_SAMPLES_PROPERTY = "packstar.pac.trainSamples";
    private static final String PAC_PILOT_SAMPLES_PROPERTY = "packstar.pac.pilotSamples";
    private static final String PAC_MAX_EST_SAMPLES_PROPERTY = "packstar.pac.maxEstSamples";
    private static final String PAC_NSTAR_INFLATE_PROPERTY = "packstar.pac.nstarInflate";
    private static final String PAC_MONITOR_SAMPLES_PROPERTY =
            "packstar.pac.monitorSamples";
    // Adaptive frequency/severity: choose a count/context-shrunk eta, repair
    // proposal-support shift with bounded on-policy refits, then freeze the
    // proposal before monitor/final samples.
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
    private static final String DP_PARALLEL_THREADS_PROPERTY = "packstar.dp.parallel.threads";
    private static final int DEFAULT_SAMPLES = 1000;
    private static final double DEFAULT_CONFIDENCE = 0.05; // delta = 0.05 => 95% confidence
    private static final double DEFAULT_TARGET_EPSILON = 0.683;
    private static final long DEFAULT_RANDOM_SEED = 42L;
    private static final int DEFAULT_SAMPLING_LARGE_LAMBDA = 65_536;
    private static final double DEFAULT_TRAIN_FRACTION = 0.5;
    private static final double DEFAULT_PILOT_FRACTION = 0.1;
    private static final int DEFAULT_MAX_EST_SAMPLES = 4000;
    private static final double DEFAULT_NSTAR_INFLATE = 1.3;
    private static final int DEFAULT_MONITOR_SAMPLES = 100;
    private static final int DEFAULT_PAC_CCD_SUBMISSION_BATCH_SIZE = 512;
    private static final boolean DEFAULT_PAC_CCD_INSTRUMENTATION = false;
    private static final String DEFAULT_FREQUENCY_SEVERITY_OUTPUT_DIR =
            "/usr/xtmp/lz280/packstar_adaptive_frequency_severity/default";
    private static final double DEFAULT_FREQUENCY_SEVERITY_RELATIVE_BOUND_KCAL = 1.0;
    private static final double DEFAULT_FREQUENCY_SEVERITY_CAP = 20.0;
    private static final String DEFAULT_FREQUENCY_SEVERITY_PREMISE_ID =
            "conditional-relative-gauge-S0-20-not-externally-recalibrated";
    private static final String DEFAULT_FREQUENCY_SEVERITY_SHRINK_GRID =
            "0:0,2:5,5:10,10:20";
    private static final String DEFAULT_FREQUENCY_SEVERITY_ALPHA_GRID = "1";
    private static final int DEFAULT_FREQUENCY_SEVERITY_FOLDS = 2;
    // Defaults copied from the successful 42-case frequency/severity run;
    // its source-shift diagnostics used permissive 1e-12 ESS thresholds.
    private static final double DEFAULT_FREQUENCY_SEVERITY_MIN_SHIFT_ESS_FRACTION = 1.0e-12;
    private static final double DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MIN_SHIFT_ESS_FRACTION = 1.0e-12;
    private static final int DEFAULT_FREQUENCY_SEVERITY_MAX_REFITS = 8;
    private static final int DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_SAMPLES = 100;
    private static final int DEFAULT_FREQUENCY_SEVERITY_DISCOVERY_MAX_SAMPLES = 400;
    private static final int DEFAULT_FREQUENCY_SEVERITY_VALIDATION_SAMPLES = 400;
    private static final boolean DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA = true;
    private static final double DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE = 1.0;
    private static final String DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_SCALE_GRID = "1";
    private static final long DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA_MAX_ASSIGNMENTS = 500_000L;
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
    /** Sorted immutable base-graph edges used by all hot conformation loops. */
    private final int[][] interactionEdges;
    /** Edge indices grouped by their lower position, preserving old pair-sum order. */
    private final int[][] interactionEdgeIndicesByFirstPosition;
    /** Lazily materialized interaction templates; templates are read-only after creation. */
    private final ResidueInteractions[][] singleInteractionTemplates;
    private final Map<Long, ResidueInteractions> pairInteractionTemplates = new HashMap<>();
    private ResidueInteractions shellInteractionTemplate;
    /** Rigid proposal matrix with proposal-only fill edges zeroed, cached by graph signature. */
    private EnergyMatrix proposalRigidEmat;
    private String proposalRigidEmatSignature;
    /** Minimizing proposal matrix with the same fill-edge gauge, cached by graph signature. */
    private EnergyMatrix proposalMinimizingEmat;
    private String proposalMinimizingEmatSignature;
    /**
     * Reusable pair-only corrected matrix.  Pair eta refits mutate only the
     * observed cells; the touched-cell lists let the next refit restore the
     * exact base values without cloning the full EMAT again.
     */
    private EnergyMatrix correctedPairwiseEmat;
    private String correctedPairwiseEmatSignature;
    private long[] correctedPairwiseOneBodyKeys = new long[0];
    private int correctedPairwiseOneBodyCount;
    private long[] correctedPairwisePairKeys = new long[0];
    private int correctedPairwisePairCount;
    private double correctedPairwiseGlobalOffsetKcal;
    /** Sampling structure/resources are reused until the proposal tree changes or compute ends. */
    private List<RootedTreeEdge> samplingTopDownOrder;
    private RootedTreeEdge samplingOrderRoot;
    private ExecutorService samplingPool;
    private boolean samplingResourcesClosed;
    /** Reuse conditional-CDF storage on each sampling worker between groups. */
    private final ThreadLocal<double[]> samplingCdfScratch =
            ThreadLocal.withInitial(() -> new double[0]);
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
    private final boolean gpuSampling;
    private final int samplingThreads;
    private final int samplingLargeLambdaThreshold;
    private final boolean samplingProgress;
    private final int ccdSubmissionBatchSize;
    private final boolean ccdInstrumentation;
    private final int trainSamples;
    private final int pilotSamples;
    private final int maxEstSamples;
    private final int sampleBudget;
    private final double nstarInflate;
    private final int monitorSamples;
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
    // Controlled ablation: keep q_m fixed while retaining calibration and PAC stages.
    private final boolean frequencySeverityProposalLearning;
    private final boolean frequencySeverityJointMomentLearning;
    private final boolean frequencySeverityBudgetForward;
    private final boolean frequencySeverityDecompositionCostSelection;
    private final PackStarTripleDecompositionCosts.Limits frequencySeverityDecompositionLimits;
    private PackStarTripleDecompositionCosts.Cache frequencySeverityDecompositionCache;
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

    // Results
    private BigDecimal zLower;
    private BigDecimal zUpper;
    private double epsilon;
    private double logZLowerPAC;
    private double logZUpperPAC;
    private int totalCCDCalls;
    /** Number of sample records presented to CCD, including duplicates. */
    private int totalCCDSampleRecords;

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
    private String functionalEventName = null;
    private PackStarFunctionalEvent functionalEvent = null;
    private PackStarFunctionalObservableResult functionalObservableResult =
            PackStarFunctionalObservableResult.notConfigured();

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
        this.interactionEdges = copySortedInteractionEdges(interactionGraph);
        this.interactionEdgeIndicesByFirstPosition =
                indexInteractionEdgesByFirstPosition(
                        rcs.getNumPos(), interactionEdges);
        this.singleInteractionTemplates = new ResidueInteractions[rcs.getNumPos()][];
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            int maxGlobalRC = -1;
            for (int local = 0; local < rcs.getNum(pos); local++) {
                maxGlobalRC = Math.max(maxGlobalRC, rcs.get(pos, local));
            }
            this.singleInteractionTemplates[pos] =
                    new ResidueInteractions[Math.max(
                            branchMinimizingEmat.getNumConfAtPos(pos),
                            maxGlobalRC + 1)];
        }
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
        this.gpuSampling = getConfigBoolean(PAC_SAMPLING_GPU_PROPERTY, false);
        this.samplingThreads = resolveConfiguredSamplingThreads();
        this.samplingLargeLambdaThreshold = Math.max(1,
                getConfigInteger(PAC_SAMPLING_LARGE_LAMBDA_PROPERTY,
                        DEFAULT_SAMPLING_LARGE_LAMBDA));
        this.samplingProgress = getConfigBoolean(PAC_SAMPLING_PROGRESS_PROPERTY, true);
        this.ccdSubmissionBatchSize = Math.max(1,
                getConfigInteger(PAC_CCD_SUBMISSION_BATCH_SIZE_PROPERTY,
                        DEFAULT_PAC_CCD_SUBMISSION_BATCH_SIZE));
        this.ccdInstrumentation = getConfigBoolean(
                PAC_CCD_INSTRUMENTATION_PROPERTY,
                DEFAULT_PAC_CCD_INSTRUMENTATION);
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
        this.monitorSamples = capSampleCount(Math.max(2,
                getConfigInteger(PAC_MONITOR_SAMPLES_PROPERTY,
                        DEFAULT_MONITOR_SAMPLES)), 2);
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
        this.frequencySeverityAlphaGrid = parseFixedOneGrid(
                getConfigProperty(PAC_FREQUENCY_SEVERITY_ALPHA_GRID_PROPERTY,
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
        this.frequencySeverityProposalLearning = getConfigBoolean(
                "packstar.pac.frequencySeverity.proposalLearning", true);
        this.frequencySeverityTripleEtaEnabled = frequencySeverityProposalLearning && getConfigBoolean(
                PAC_FREQUENCY_SEVERITY_TRIPLE_ETA_PROPERTY,
                DEFAULT_FREQUENCY_SEVERITY_TRIPLE_ETA);
        this.frequencySeverityJointMomentLearning = frequencySeverityProposalLearning && getConfigBoolean(
                "packstar.pac.frequencySeverity.jointMomentLearning", true);
        this.frequencySeverityDecompositionCostSelection = frequencySeverityJointMomentLearning
                && PackStarTripleDecompositionCosts.configuredEnabled();
        this.frequencySeverityBudgetForward = frequencySeverityJointMomentLearning
                && PackStarTripleDecompositionCosts.configuredStrategy().equals("budget-forward");
        this.frequencySeverityDecompositionLimits = frequencySeverityDecompositionCostSelection
                ? PackStarTripleDecompositionCosts.Limits.configured() : null;
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
        if (configuredFrequencySeverityOutputDir == null) {
            String benchOutputDir = getConfigProperty(
                    "osprey.bench.outputDir", null);
            String trimmedBenchOutputDir = benchOutputDir == null
                    ? null : benchOutputDir.trim();
            configuredFrequencySeverityOutputDir = isXtmpOutputPath(
                    trimmedBenchOutputDir)
                    ? new File(trimmedBenchOutputDir,
                            "adaptive_frequency_severity").getPath()
                    : DEFAULT_FREQUENCY_SEVERITY_OUTPUT_DIR;
        }
        this.frequencySeverityOutputDir =
                configuredFrequencySeverityOutputDir.trim();
        if (!Double.isFinite(frequencySeverityRelativeBoundKcal)
                || !(frequencySeverityRelativeBoundKcal > 0.0)) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity relative bound must be finite and positive: "
                            + frequencySeverityRelativeBoundKcal);
        }
        if (!Double.isFinite(frequencySeverityCap)
                || frequencySeverityCap < 0.0) {
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
        if (frequencySeverityMaxRefits < 0
                || frequencySeverityMaxRefits > 16) {
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
        if (frequencySeverityBudgetForward && frequencySeverityTripleEtaMaxPositionTriples > 3)
            throw new IllegalArgumentException("budget-forward final protocol supports K=0,1,2,3 only");
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
                || !(frequencySeveritySizeSafety > 0.0)
                || frequencySeveritySizeSafety > 1.0) {
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
        if (frequencySeverityOutputDir == null
                || !isXtmpOutputPath(frequencySeverityOutputDir)) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity production requires an absolute output directory under"
                            + " /usr/xtmp/lz280: " + frequencySeverityOutputDir);
        }
        if (trainSamples < 4) {
            throw new IllegalArgumentException(
                    "PACK* frequency/severity two-fold selection requires at least four training samples");
        }
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

    /**
     * Configure one event to measure on the independent final sample.
     *
     * <p>The event is evaluated after proposal learning and after the final
     * proposal is frozen.  Passing null for both arguments disables the
     * observable.</p>
     */
    public void setFunctionalEvent(String name, PackStarFunctionalEvent event) {
        if (event == null) {
            if (name != null && !name.trim().isEmpty()) {
                throw new IllegalArgumentException(
                        "a functional event name requires a non-null event");
            }
            functionalEventName = null;
            functionalEvent = null;
            functionalObservableResult =
                    PackStarFunctionalObservableResult.notConfigured();
            return;
        }
        if (name == null || name.trim().isEmpty()) {
            throw new IllegalArgumentException(
                    "a functional event requires a non-empty name");
        }
        functionalEventName = name.trim();
        functionalEvent = event;
        functionalObservableResult =
                PackStarFunctionalObservableResult.notComputed(functionalEventName);
    }

    /** Return the most recent functional-event estimate and diagnostics. */
    public PackStarFunctionalObservableResult getFunctionalObservableResult() {
        return functionalObservableResult;
    }

    private static class CorrectedDPResult {
        final double logZCorrected;

        CorrectedDPResult(double logZCorrected) {
            this.logZCorrected = logZCorrected;
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
        /** Local M digits captured while grouping; avoids a second decode. */
        int[] mRCs;
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

    private static int[][] copySortedInteractionEdges(InteractionGraph graph) {
        List<int[]> edges = new ArrayList<>(graph.getEdgeList().size());
        for (int[] edge : graph.getEdgeList()) {
            if (edge == null || edge.length < 2) {
                throw new IllegalArgumentException(
                        "interaction graph contains a malformed edge");
            }
            int pos1 = Math.min(edge[0], edge[1]);
            int pos2 = Math.max(edge[0], edge[1]);
            edges.add(new int[]{pos1, pos2});
        }
        edges.sort((a, b) -> {
            int compare = Integer.compare(a[0], b[0]);
            return compare != 0 ? compare : Integer.compare(a[1], b[1]);
        });
        // InteractionGraph.hasEdge() is set-like.  Keep the same semantics if
        // a caller supplied a duplicate edge, otherwise the cached edge list
        // would silently double-count every pair term.
        List<int[]> unique = new ArrayList<>(edges.size());
        int previousPos1 = -1;
        int previousPos2 = -1;
        for (int[] edge : edges) {
            if (edge[0] == previousPos1 && edge[1] == previousPos2) {
                continue;
            }
            unique.add(edge);
            previousPos1 = edge[0];
            previousPos2 = edge[1];
        }
        return unique.toArray(new int[0][]);
    }

    private static int[][] indexInteractionEdgesByFirstPosition(
            int numPositions, int[][] edges) {
        List<Integer>[] byPosition = new List[numPositions];
        for (int pos = 0; pos < numPositions; pos++) {
            byPosition[pos] = new ArrayList<>();
        }
        for (int edgeIndex = 0; edgeIndex < edges.length; edgeIndex++) {
            byPosition[edges[edgeIndex][0]].add(edgeIndex);
        }
        int[][] result = new int[numPositions][];
        for (int pos = 0; pos < numPositions; pos++) {
            result[pos] = byPosition[pos].stream()
                    .mapToInt(Integer::intValue).toArray();
        }
        return result;
    }

    private ResidueInteractions getSingleInteractionTemplate(int pos, int rc) {
        if (pos < 0 || pos >= singleInteractionTemplates.length
                || rc < 0 || rc >= singleInteractionTemplates[pos].length) {
            throw new IllegalArgumentException(
                    "single interaction template index out of range: pos="
                            + pos + " rc=" + rc);
        }
        ResidueInteractions template = singleInteractionTemplates[pos][rc];
        if (template == null) {
            template = minimizingEcalc.makeSingleInters(pos, rc);
            singleInteractionTemplates[pos][rc] = template;
        }
        return template;
    }

    private static long packInteractionTemplateKey(
            int pos1, int rc1, int pos2, int rc2) {
        if (pos1 > pos2) {
            int tmp = pos1;
            pos1 = pos2;
            pos2 = tmp;
            tmp = rc1;
            rc1 = rc2;
            rc2 = tmp;
        }
        return ((long) pos1 << 48) | ((long) rc1 << 32)
                | ((long) pos2 << 16) | (rc2 & 0xffffL);
    }

    private ResidueInteractions getPairInteractionTemplate(
            int pos1, int rc1, int pos2, int rc2) {
        long key = packInteractionTemplateKey(pos1, rc1, pos2, rc2);
        ResidueInteractions template = pairInteractionTemplates.get(key);
        if (template == null) {
            template = minimizingEcalc.makePairInters(pos1, rc1, pos2, rc2);
            pairInteractionTemplates.put(key, template);
        }
        return template;
    }

    private ResidueInteractions getShellInteractionTemplate() {
        if (shellInteractionTemplate == null) {
            shellInteractionTemplate = minimizingEcalc.makeShellInters();
        }
        return shellInteractionTemplate;
    }

    /**
     * Run the full PAC estimation pipeline.
     * Returns the epsilon achieved.
     */
    public double compute() {
        long startTime = System.currentTimeMillis();

        // A PackStarEstimator is normally used once, but resetting these
        // counters makes repeated compute() calls well-defined and keeps the
        // reported CCD savings tied to this run only.
        totalCCDCalls = 0;
        totalCCDSampleRecords = 0;
        functionalObservableResult = functionalEvent == null
                ? PackStarFunctionalObservableResult.notConfigured()
                : PackStarFunctionalObservableResult.notComputed(functionalEventName);

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

        System.out.println("[PACK*-adaptive-frequency-severity] fixed production route"
                + ", samples(train/pilot/maxEst)="
                + trainSamples + "/" + pilotSamples + "/" + maxEstSamples
                + ", monitorSamples=" + monitorSamples
                + ", sampleBudget="
                + (sampleBudget == Integer.MAX_VALUE ? "unbounded" : sampleBudget)
                + ", B_rel=" + frequencySeverityRelativeBoundKcal
                + ", S0=" + frequencySeverityCap
                + ", maxRefits=" + frequencySeverityMaxRefits);
        try {
            runAdaptiveFrequencySeverityPAC(startTime, logZRigid);
        } finally {
            shutdownSamplingResources();
        }

        return epsilon;
    }

    private void printFinalSummary(long startTime) {
        long totalTime = System.currentTimeMillis() - startTime;
        System.out.println("[PACK*] Total: " + totalTime + " ms, " + totalCCDCalls
                + " unique CCD calls for " + totalCCDSampleRecords
                + " sample records, epsilon=" + String.format("%.6f", epsilon)
                + ", confidence=" + String.format("%.2f%%", (1.0 - delta) * 100));
        System.out.println("[PACK*] Z bounds: lower="
                + String.format(Locale.ROOT, "%.6e", zLower)
                + ", upper="
                + String.format(Locale.ROOT, "%.6e", zUpper)
                + ", log10Lower=" + formatLog10(logZLowerPAC)
                + ", log10Upper=" + formatLog10(logZUpperPAC));
    }

    // ========== Adaptive frequency/severity estimator ==========

    /**
     * Adaptive frequency/severity PAC estimator.
     *
     *   Z = q_eta * E_{p_eta}[ exp(-xi/RT) ],   xi = E_true - E_eta   (exact identity)
     *
     * The q_m training batch is cross-fitted to choose a count/context-shrunk
     * eta.  Bounded on-policy discovery refits repair proposal-support shift;
     * the selected proposal is then frozen before validation, monitor, and
     * final samples.  The pilot sizes the fresh final sample once the proposal
     * is frozen; there is no alternate conditional/unconditional estimator and
     * no proposal restore after the final selection.
     *
     * Sample-splitting: eta and p_eta depend only on training/adaptation data,
     * while monitor/final samples are fresh.  Conditional on the frozen
     * severity premise and the pilot, the final sample size is fixed for the
     * empirical-Bernstein PAC calculation.
     */
    private static class FrequencySeverityEtaCoverage {
        final int[][] unaryDistinctContexts;
        final EtaCorrections.PairCountMap pairDistinctContexts;

        FrequencySeverityEtaCoverage(int numPos, int[] numRCs) {
            unaryDistinctContexts = new int[numPos][];
            for (int pos = 0; pos < numPos; pos++) {
                unaryDistinctContexts[pos] = new int[numRCs[pos]];
            }
            pairDistinctContexts = new EtaCorrections.PairCountMap();
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
        int requestedTripleCount;
        int[] foldTripleCounts = new int[2];
        String pathStopReason = "legacy-or-pair-only";
        String[] foldPathStopReasons = {"legacy-or-pair-only", "legacy-or-pair-only"};
        double innerLogRhoGain = Double.NaN;
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
            this.requestedTripleCount = hasTripleEta(fullEta) ? (int) fullEta.tripleEta.positionTripleCount : 0;
            if (this.requestedTripleCount == 0) this.innerLogRhoGain = 0;
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

    /** Scratch arrays reused across sequential candidate/fold scoring. */
    private static class FrequencySeverityFoldScratch {
        double[] logA = new double[0];
        double[] logV = new double[0];
        double[] proposalWeights = new double[0];
        double[] logRelative = new double[0];
        double[] bulk = new double[0];

        void ensureCapacity(int n) {
            if (logA.length >= n) return;
            logA = new double[n];
            logV = new double[n];
            proposalWeights = new double[n];
            logRelative = new double[n];
            bulk = new double[n];
        }
    }

    /**
     * Cache eta energies during one candidate-scoring pass.
     *
     * The pair part is keyed separately from the selected triple residual so
     * that the pair-only candidate and its corresponding triple candidate can
     * share the same full/fold pair energies.  The cache is deliberately
     * scoped to one scoring pass: CCDResult features remain eta-independent,
     * while these energies depend on the fitted candidate model.
     */
    private class FrequencySeverityEtaEnergyCache {

        private static class CachedEnergyVector {
            final double[] values;
            final boolean[] computed;

            CachedEnergyVector(int size) {
                values = new double[size];
                computed = new boolean[size];
            }
        }

        private final IdentityHashMap<CCDResult, Integer> sampleIndices =
                new IdentityHashMap<>();
        private final IdentityHashMap<EtaCorrections, CachedEnergyVector>
                pairEnergies = new IdentityHashMap<>();
        private final IdentityHashMap<EtaCorrections, CachedEnergyVector>
                tripleResidualEnergies = new IdentityHashMap<>();

        FrequencySeverityEtaEnergyCache(
                List<CCDResult> trainingSamples,
                List<CCDResult>[] validationFolds,
                FrequencySeverityEtaTraining training) {
            addSamples(trainingSamples);
            addSamples(training.allFitSamples);
            for (List<CCDResult> foldSamples : training.foldFitSamples) {
                addSamples(foldSamples);
            }
            for (List<CCDResult> foldSamples : validationFolds) {
                addSamples(foldSamples);
            }
        }

        private void addSamples(List<CCDResult> samples) {
            if (samples == null) return;
            for (CCDResult sample : samples) {
                if (!sampleIndices.containsKey(sample)) {
                    sampleIndices.put(sample, sampleIndices.size());
                }
            }
        }

        private int sampleIndex(CCDResult sample) {
            Integer index = sampleIndices.get(sample);
            if (index == null) {
                throw new IllegalArgumentException(
                        "eta energy cache does not contain CCD sample");
            }
            return index;
        }

        private CachedEnergyVector getPairEnergyVector(
                EtaCorrections pairEta) {
            CachedEnergyVector vector = pairEnergies.get(pairEta);
            if (vector == null) {
                vector = new CachedEnergyVector(sampleIndices.size());
                pairEnergies.put(pairEta, vector);
            }
            return vector;
        }

        private CachedEnergyVector getTripleResidualVector(
                EtaCorrections eta) {
            CachedEnergyVector vector = tripleResidualEnergies.get(eta);
            if (vector == null) {
                vector = new CachedEnergyVector(sampleIndices.size());
                tripleResidualEnergies.put(eta, vector);
            }
            return vector;
        }

        double getPairEnergy(CCDResult sample, EtaCorrections pairEta) {
            if (pairEta.tripleEta != null
                    && pairEta.tripleEtaScale != 0.0) {
                throw new IllegalArgumentException(
                        "pair eta energy cache received triple eta model");
            }
            int index = sampleIndex(sample);
            CachedEnergyVector vector = getPairEnergyVector(pairEta);
            if (!vector.computed[index]) {
                vector.values[index] = computeEtaEnergy(sample, pairEta);
                // Cache non-finite values as well.  The caller must see the
                // same invalid-candidate behavior on every reuse.
                vector.computed[index] = true;
            }
            return vector.values[index];
        }

        private double getTripleResidualEnergy(
                CCDResult sample, EtaCorrections eta) {
            int index = sampleIndex(sample);
            CachedEnergyVector vector = getTripleResidualVector(eta);
            if (!vector.computed[index]) {
                vector.values[index] = eta.tripleEta.scoreResidual(
                        sample.conf, eta::getPairEta);
                vector.computed[index] = true;
            }
            return vector.values[index];
        }

        double getEnergy(
                CCDResult sample, EtaCorrections eta,
                EtaCorrections pairEta) {
            double energy = getPairEnergy(sample, pairEta);
            if (eta.tripleEta != null && eta.tripleEtaScale != 0.0) {
                energy += eta.tripleEtaScale
                        * getTripleResidualEnergy(sample, eta);
            }
            return energy;
        }
    }

    private static class FrequencySeverityCoverageFlags {
        final boolean[] unaryUndertrained;
        final boolean[] pairUndertrained;

        FrequencySeverityCoverageFlags(int sampleCount) {
            unaryUndertrained = new boolean[sampleCount];
            pairUndertrained = new boolean[sampleCount];
        }
    }

    private FrequencySeverityCoverageFlags buildFrequencySeverityCoverageFlags(
            List<CCDResult> samples, FrequencySeverityEtaCoverage coverage) {
        FrequencySeverityCoverageFlags flags =
                new FrequencySeverityCoverageFlags(samples.size());
        for (int i = 0; i < samples.size(); i++) {
            int[] conf = samples.get(i).conf;
            flags.unaryUndertrained[i] =
                    touchesFrequencySeverityUndertrainedUnary(conf, coverage);
            flags.pairUndertrained[i] =
                    touchesFrequencySeverityUndertrainedPair(conf, coverage);
        }
        return flags;
    }

    private List<FrequencySeverityCandidateScore> scoreFrequencySeverityEtaCandidates(
            List<CCDResult> trainingSamples,
            FrequencySeverityEtaTraining training) {
        List<CCDResult>[] validationFolds = new List[frequencySeverityFolds];
        FrequencySeverityCoverageFlags[] validationCoverageFlags =
                new FrequencySeverityCoverageFlags[frequencySeverityFolds];
        for (int fold = 0; fold < frequencySeverityFolds; fold++) {
            validationFolds[fold] = new ArrayList<>();
        }
        for (int index = 0; index < trainingSamples.size(); index++) {
            validationFolds[index % frequencySeverityFolds].add(trainingSamples.get(index));
        }
        for (int validationFold = 0;
             validationFold < frequencySeverityFolds; validationFold++) {
            int fitFold = 1 - validationFold;
            validationCoverageFlags[validationFold] =
                    buildFrequencySeverityCoverageFlags(
                            validationFolds[validationFold],
                            training.foldCoverage[fitFold]);
        }
        FrequencySeverityCoverageFlags allCoverageFlags =
                buildFrequencySeverityCoverageFlags(
                        trainingSamples, training.allCoverage);

        List<FrequencySeverityCandidateScore> candidates = new ArrayList<>();
        FrequencySeverityFoldScratch foldScratch =
                new FrequencySeverityFoldScratch();
        FrequencySeverityEtaEnergyCache energyCache =
                new FrequencySeverityEtaEnergyCache(
                        trainingSamples, validationFolds, training);
        if (!frequencySeverityProposalLearning) {
            // Extraction retains support/provenance but leaves every eta cell zero.
            // Evaluate one fixed proposal, without fitting or selecting corrections.
            FrequencySeverityFoldScore[] scores =
                    new FrequencySeverityFoldScore[frequencySeverityFolds];
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                scores[fold] = scoreFrequencySeverityEtaFold(
                        fold, validationFolds[fold], training.folds[1 - fold],
                        training.folds[1 - fold], foldScratch,
                        validationCoverageFlags[fold], energyCache);
            }
            FrequencySeverityCandidateScore fixed = new FrequencySeverityCandidateScore(
                    "fixed-qm-no-learning", new FrequencySeverityShrinkPair(0.0, 0.0),
                    0.0, 0.0, training.all, scores,
                    poolFrequencySeverityFoldScores(scores),
                    auditFrequencySeverityRefit(trainingSamples, training.all,
                            training.all, allCoverageFlags, foldScratch, energyCache));
            fixed.pathStopReason = "proposal-learning-disabled";
            Arrays.fill(fixed.foldPathStopReasons, "proposal-learning-disabled");
            candidates.add(fixed);
            return candidates;
        }
        // Alpha remains fixed at one. Crossfit pair-only and the retained
        // triple prefixes; legacy comparison strategies retain only their
        // terminal set. No fitted correction is continuously attenuated.
        final double alpha = 1.0;
        for (FrequencySeverityShrinkPair shrink : frequencySeverityShrinkGrid) {
            EtaCorrections fullPairEta = shrinkEta(
                    training.all, training.allCoverage,
                    shrink.unary, shrink.pair);
            if (frequencySeverityJointMomentLearning) {
                fullPairEta = refitFrequencySeverityPairEta(training.allFitSamples,
                        fullPairEta, training.allCoverage, shrink);
            }
            EtaCorrections[] foldPairEta =
                    new EtaCorrections[frequencySeverityFolds];
            for (int fitFold = 0;
                 fitFold < frequencySeverityFolds; fitFold++) {
                foldPairEta[fitFold] = shrinkEta(
                        training.folds[fitFold],
                        training.foldCoverage[fitFold],
                        shrink.unary, shrink.pair);
                if (frequencySeverityJointMomentLearning) {
                    foldPairEta[fitFold] = refitFrequencySeverityPairEta(
                            training.foldFitSamples[fitFold], foldPairEta[fitFold],
                            training.foldCoverage[fitFold], shrink);
                }
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
                        foldPairEta[fitFold],
                        foldScratch,
                        validationCoverageFlags[validationFold],
                        energyCache);
            }
            FrequencySeverityCrossfitScore pairCrossfit =
                    poolFrequencySeverityFoldScores(pairFoldScores);
            FrequencySeverityRefitAudit pairRefitAudit =
                    auditFrequencySeverityRefit(
                    trainingSamples, fullPairEta,
                            fullPairEta, allCoverageFlags, foldScratch,
                            energyCache);
            candidates.add(new FrequencySeverityCandidateScore(
                    frequencySeverityCandidateId(
                            shrink, alpha, 0.0),
                    shrink, alpha, 0.0, fullPairEta,
                    pairFoldScores, pairCrossfit, pairRefitAudit));

            if (frequencySeverityTripleEta == null) continue;
            if (frequencySeverityBudgetForward) {
                scoreFrequencySeverityTriplePath(candidates, shrink, training, trainingSamples,
                        validationFolds, fullPairEta, foldPairEta, validationCoverageFlags,
                        allCoverageFlags, foldScratch, energyCache);
                continue;
            }
            PackStarTripleEtaCorrections fullSelectedTriple =
                    fitFrequencySeverityTripleEta(
                            training.allFitSamples, fullPairEta,
                            energyCache, shrink);
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
                                foldPairEta[fitFold], energyCache, shrink);
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
                        foldPairEta[fitFold],
                        foldScratch,
                        validationCoverageFlags[validationFold],
                        energyCache);
            }
            FrequencySeverityCrossfitScore tripleCrossfit =
                    poolFrequencySeverityFoldScores(tripleFoldScores);
            FrequencySeverityRefitAudit tripleRefitAudit =
                    auditFrequencySeverityRefit(
                    trainingSamples, fullTripleEta, fullPairEta,
                            allCoverageFlags, foldScratch, energyCache);
            FrequencySeverityCandidateScore legacyCandidate = new FrequencySeverityCandidateScore(
                    frequencySeverityCandidateId(shrink, alpha, 1.0),
                    shrink, alpha, 1.0, fullTripleEta,
                    tripleFoldScores, tripleCrossfit, tripleRefitAudit);
            for (int f = 0; f < frequencySeverityFolds; f++)
                legacyCandidate.foldTripleCounts[f] = hasTripleEta(foldSelectedTriples[f])
                        ? (int) foldSelectedTriples[f].positionTripleCount : 0;
            candidates.add(legacyCandidate);
        }
        return candidates;
    }

    private void scoreFrequencySeverityTriplePath(
            List<FrequencySeverityCandidateScore> candidates, FrequencySeverityShrinkPair shrink,
            FrequencySeverityEtaTraining training, List<CCDResult> trainingSamples,
            List<CCDResult>[] validationFolds, EtaCorrections fullPairEta, EtaCorrections[] foldPairEta,
            FrequencySeverityCoverageFlags[] validationCoverageFlags,
            FrequencySeverityCoverageFlags allCoverageFlags, FrequencySeverityFoldScratch scratch,
            FrequencySeverityEtaEnergyCache cache) {
        PackStarTripleEtaCorrections.MomentPath full = fitFrequencySeverityTriplePath(
                training.allFitSamples, fullPairEta, shrink, frequencySeverityTripleEtaMaxPositionTriples);
        PackStarTripleEtaCorrections.MomentPath[] paths = new PackStarTripleEtaCorrections.MomentPath[2];
        // Only K values present in the full path become scoring candidates.
        // The at-most-K fold procedure has the same prefixes regardless of
        // its eventual upper limit. Deeper fold fits cannot affect any score.
        int requiredK = full.models.size() - 1;
        for (int f = 0; f < 2; f++) paths[f] = fitFrequencySeverityTriplePath(
                training.foldFitSamples[f], foldPairEta[f], shrink, requiredK);
        FrequencySeverityCandidateScore pairCandidate = candidates.get(candidates.size() - 1);
        pairCandidate.pathStopReason = full.stopReason;
        for (int f = 0; f < 2; f++) pairCandidate.foldPathStopReasons[f] = paths[f].stopReason;
        for (int k = 1; k < full.models.size(); k++) {
            EtaCorrections fullEta = attachFrequencySeverityTripleEta(fullPairEta, full.models.get(k), 1);
            FrequencySeverityFoldScore[] scores = new FrequencySeverityFoldScore[2];
            for (int validationFold = 0; validationFold < 2; validationFold++) {
                int fitFold = 1 - validationFold;
                // Crossfit the at-most-K procedure. A fold may stop earlier;
                // never force a factor without evidence or reuse full-fit scopes.
                EtaCorrections eta = attachFrequencySeverityTripleEta(foldPairEta[fitFold], paths[fitFold].atMost(k), 1);
                scores[validationFold] = scoreFrequencySeverityEtaFold(validationFold,
                        validationFolds[validationFold], eta, foldPairEta[fitFold], scratch,
                        validationCoverageFlags[validationFold], cache);
            }
            FrequencySeverityCandidateScore candidate = new FrequencySeverityCandidateScore(
                    frequencySeverityCandidateId(shrink, 1, 1) + "-K-" + k,
                    shrink, 1, 1, fullEta, scores, poolFrequencySeverityFoldScores(scores),
                    auditFrequencySeverityRefit(trainingSamples, fullEta, fullPairEta,
                            allCoverageFlags, scratch, cache));
            candidate.requestedTripleCount = k;
            candidate.pathStopReason = full.stopReason;
            candidate.innerLogRhoGain = full.gains.get(k);
            for (int f = 0; f < 2; f++) {
                candidate.foldTripleCounts[f] = Math.min(k, paths[f].models.size() - 1);
                candidate.foldPathStopReasons[f] = paths[f].stopReason;
            }
            candidates.add(candidate);
        }
    }

    private PackStarTripleEtaCorrections.MomentPath fitFrequencySeverityTriplePath(
            List<CCDResult> samples, EtaCorrections pairEta, FrequencySeverityShrinkPair shrink,
            int requiredK) {
        PackStarTripleEtaCorrections.MomentPath empty = new PackStarTripleEtaCorrections.MomentPath(rcs.getNumPos());
        if (requiredK == 0) { empty.stopReason = "required-prefix-limit"; return empty; }
        PackStarProposalLearning.Data[] inner = new PackStarProposalLearning.Data[2];
        for (int heldOut = 0; heldOut < 2; heldOut++) {
            List<CCDResult> fit = new ArrayList<>();
            for (CCDResult sample : samples)
                if (PackStarProposalLearning.fold(sample.conf) != heldOut) fit.add(sample);
            if (fit.size() < 2) { empty.stopReason = "insufficient-inner-data"; return empty; }
            FrequencySeverityEtaTraining training = extractFrequencySeverityEtaTraining(fit);
            EtaCorrections prior = shrinkEta(training.all, training.allCoverage, shrink.unary, shrink.pair);
            EtaCorrections pair = refitFrequencySeverityPairEta(fit, prior, training.allCoverage, shrink);
            inner[heldOut] = proposalLearningData(samples, pair);
        }
        if (frequencySeverityDecompositionCache == null)
            frequencySeverityDecompositionCache = new PackStarTripleDecompositionCosts.Cache(
                    rcs, interactionGraph, initialRootedRoot,
                    getConfigInteger(PackStarTripleDecompositionCosts.PREFIX + "PreviewCacheEntries", 1024), confSpace);
        PackStarTripleEtaCorrections.MomentPath path = frequencySeverityTripleEta.fitSecondMomentPath(rcs, interactionGraph,
                proposalLearningData(samples, pairEta), inner, RT,
                requiredK, frequencySeverityTripleEtaMaxFillEdges,
                frequencySeverityTripleEtaMinCellContexts, frequencySeverityTripleEtaPriorStrength,
                frequencySeverityTripleEtaResidualCapKcal, frequencySeverityTripleEtaMaxAssignments,
                frequencySeverityDecompositionCache, frequencySeverityDecompositionLimits);
        if (requiredK < frequencySeverityTripleEtaMaxPositionTriples && path.stopReason.equals("maximum-K"))
            path.stopReason = "required-prefix-limit";
        return path;
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
            List<CCDResult> samples, EtaCorrections pairEta,
            FrequencySeverityEtaEnergyCache energyCache,
            FrequencySeverityShrinkPair shrink) {
        if (frequencySeverityTripleEta == null
                || samples == null || samples.size() < 2) {
            return null;
        }
        if (frequencySeverityJointMomentLearning) {
            PackStarProposalLearning.Data[] inner = new PackStarProposalLearning.Data[2];
            for (int heldOut = 0; heldOut < 2; heldOut++) {
                List<CCDResult> fit = new ArrayList<>();
                for (CCDResult sample : samples)
                    if (PackStarProposalLearning.fold(sample.conf) != heldOut) fit.add(sample);
                if (fit.size() < 2) return null;
                // Refit the entire pair pipeline without this inner fold; a
                // frozen pair fit on all samples would leak held-out CCD labels.
                FrequencySeverityEtaTraining training = extractFrequencySeverityEtaTraining(fit);
                EtaCorrections prior = shrinkEta(training.all, training.allCoverage,
                        shrink.unary, shrink.pair);
                EtaCorrections pair = refitFrequencySeverityPairEta(fit, prior,
                        training.allCoverage, shrink);
                inner[heldOut] = proposalLearningData(samples, pair);
            }
            if (frequencySeverityDecompositionCostSelection && frequencySeverityDecompositionCache == null) {
                frequencySeverityDecompositionCache = new PackStarTripleDecompositionCosts.Cache(
                        rcs, interactionGraph, initialRootedRoot,
                        getConfigInteger(PackStarTripleDecompositionCosts.PREFIX + "PreviewCacheEntries", 1024), confSpace);
            }
            return frequencySeverityTripleEta.fitSelectedSecondMoment(
                    rcs, interactionGraph, proposalLearningData(samples, pairEta), inner, RT,
                    frequencySeverityTripleEtaMaxPositionTriples,
                    frequencySeverityTripleEtaMaxFillEdges,
                    frequencySeverityTripleEtaMinCellContexts,
                    frequencySeverityTripleEtaPriorStrength,
                    frequencySeverityTripleEtaResidualCapKcal,
                    frequencySeverityTripleEtaMaxAssignments,
                    frequencySeverityDecompositionCache, frequencySeverityDecompositionLimits);
        }
        int[][] conformations = new int[samples.size()][];
        double[] residuals = new double[samples.size()];
        for (int sample = 0; sample < samples.size(); sample++) {
            CCDResult result = samples.get(sample);
            conformations[sample] = result.conf;
            residuals[sample] = result.eTrue - result.eMin
                    - energyCache.getPairEnergy(result, pairEta);
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

    private PackStarProposalLearning.Data proposalLearningData(
            List<CCDResult> samples, EtaCorrections pair) {
        int[][] conf = new int[samples.size()][];
        double[] residual = new double[samples.size()], logWeight = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            CCDResult sample = samples.get(i);
            double energy = sample.eMin + computeEtaEnergy(sample, pair);
            conf[i] = sample.conf;
            residual[i] = sample.eTrue - energy;
            logWeight[i] = frequencySeveritySourceLogWeight(sample, energy);
        }
        return new PackStarProposalLearning.Data(conf, residual, logWeight);
    }

    /** Jointly fit full CCD residual using only cells of the existing graph. */
    private EtaCorrections refitFrequencySeverityPairEta(List<CCDResult> samples,
            EtaCorrections prior, FrequencySeverityEtaCoverage coverage,
            FrequencySeverityShrinkPair shrink) {
        if (samples.size() < 2) return prior;
        int[][] unaryIds = new int[prior.oneBody.length][];
        List<Double> ratios = new ArrayList<>();
        for (int p = 0; p < unaryIds.length; p++) {
            unaryIds[p] = new int[prior.oneBody[p].length];
            Arrays.fill(unaryIds[p], -1);
            for (int rc = 0; rc < unaryIds[p].length; rc++) {
                if (prior.oneBodyCounts[p][rc] == 0) continue;
                unaryIds[p][rc] = ratios.size();
                ratios.add(shrink.unary / Math.max(1, coverage.unary(p, rc)));
            }
        }
        long[] pairKeys = prior.pairSums.keys();
        Arrays.sort(pairKeys);
        Map<Long, Integer> pairIds = new HashMap<>();
        for (long key : pairKeys) {
            pairIds.put(key, ratios.size());
            ratios.add(shrink.pair / Math.max(1,
                    coverage.pairDistinctContexts.getOrDefault(key, 0)));
        }
        int[][] features = new int[samples.size()][];
        for (int i = 0; i < samples.size(); i++) {
            int[] conf = samples.get(i).conf;
            int[] active = new int[conf.length + interactionEdges.length];
            int count = 0;
            for (int p = 0; p < conf.length; p++) {
                if (conf[p] >= 0 && unaryIds[p][conf[p]] >= 0)
                    active[count++] = unaryIds[p][conf[p]];
            }
            for (int[] edge : interactionEdges) {
                int a = edge[0], b = edge[1];
                if (conf[a] < 0 || conf[b] < 0) continue;
                Integer id = pairIds.get(EtaCorrections.packPairKey(a, conf[a], b, conf[b]));
                if (id != null) active[count++] = id;
            }
            features[i] = Arrays.copyOf(active, count);
        }
        PackStarProposalLearning.Data data = proposalLearningData(samples, prior);
        double[] regularization = new double[ratios.size()];
        for (int j = 0; j < ratios.size(); j++) regularization[j] = ratios.get(j);
        PackStarProposalLearning.JointFit fit = PackStarProposalLearning.jointFit(
                features, data.residual, data.logWeight, regularization,
                frequencySeverityTripleEtaResidualCapKcal, 40);
        for (int p = 0; p < unaryIds.length; p++) {
            for (int rc = 0; rc < unaryIds[p].length; rc++) {
                int id = unaryIds[p][rc];
                if (id >= 0) prior.oneBody[p][rc] += fit.delta[id] * prior.oneBodyCounts[p][rc];
            }
        }
        for (long key : pairKeys) {
            int count = prior.pairSums.count(key);
            prior.pairSums.put(key, prior.pairSums.sum(key) + fit.delta[pairIds.get(key)] * count, count);
        }
        prior.globalOffsetKcal += fit.offset;
        System.out.println("[PACK*-joint-pair-fit] samples=" + samples.size()
                + ", observedCells=" + ratios.size() + ", offsetKcal=" + fit.offset);
        return prior;
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
            EtaCorrections pairEta,
            FrequencySeverityFoldScratch scratch,
            FrequencySeverityCoverageFlags coverageFlags,
            FrequencySeverityEtaEnergyCache energyCache) {
        int n = samples.size();
        if (n <= 1) {
            return invalidFrequencySeverityFoldScore(validationFold, n);
        }
        scratch.ensureCapacity(n);
        double[] logA = scratch.logA;
        double[] logV = scratch.logV;
        for (int i = 0; i < n; i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + energyCache.getEnergy(sample, eta, pairEta);
            logA[i] = frequencySeveritySourceLogWeight(
                    sample, targetEnergy);
            logV[i] = logA[i]
                    + (targetEnergy - sample.eTrue) / RT;
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
        double[] proposalWeights = scratch.proposalWeights;
        double[] logRelative = scratch.logRelative;
        double sumWeightSquares = 0.0;
        double bulkMean = 0.0;
        double tailProbability = 0.0;
        double candidateUnaryUndertrained = 0.0;
        double candidatePairUndertrained = 0.0;
        int baselineUnaryUndertrained = 0;
        int baselinePairUndertrained = 0;
        double logClip = frequencySeverityRelativeBoundKcal / RT;
        double[] bulk = scratch.bulk;
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
            if (coverageFlags.unaryUndertrained[i]) {
                baselineUnaryUndertrained++;
                candidateUnaryUndertrained += weight;
            }
            if (coverageFlags.pairUndertrained[i]) {
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
        if (frequencySeverityBudgetForward) {
            List<FrequencySeverityCandidateScore> eligible = new ArrayList<>();
            for (FrequencySeverityCandidateScore candidate : candidates)
                if (candidate.eligible) eligible.add(candidate);
            return PackStarProposalLearning.chooseWithinFivePercent(eligible,
                    candidate -> candidate.crossfit.sizing.finalSamples,
                    (a, b) -> {
                        int c = Integer.compare(a.requestedTripleCount, b.requestedTripleCount);
                        if (c == 0 && a.requestedTripleCount > 0)
                            c = PackStarTripleDecompositionCosts.compare(candidateStructure(a), candidateStructure(b));
                        if (c != 0) return c;
                        if (isBetterFrequencySeverityCandidate(a, b)) return -1;
                        if (isBetterFrequencySeverityCandidate(b, a)) return 1;
                        return 0;
                    });
        }
        FrequencySeverityCandidateScore best = null;
        for (FrequencySeverityCandidateScore candidate : candidates) {
            if (!candidate.eligible) continue;
            if (best == null || isBetterFrequencySeverityCandidate(candidate, best)) {
                best = candidate;
            }
        }
        return best;
    }

    private PackStarTripleDecompositionCosts.Cost candidateStructure(FrequencySeverityCandidateScore candidate) {
        Set<Long> fill = new LinkedHashSet<>();
        for (int[] edge : candidate.fullEta.tripleEta.requiredFillEdges(interactionGraph))
            fill.add(((long) edge[0] << 32) | edge[1]);
        return frequencySeverityDecompositionCache.preview(fill);
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
            if ((frequencySeverityProposalLearning && !(candidate.alpha > 0.0))
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
            EtaCorrections pairEta,
            FrequencySeverityCoverageFlags coverageFlags,
            FrequencySeverityFoldScratch scratch,
            FrequencySeverityEtaEnergyCache energyCache) {
        int n = samples.size();
        scratch.ensureCapacity(n);
        double[] logA = scratch.logA;
        double logSumA = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + energyCache.getEnergy(sample, eta, pairEta);
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
            if (coverageFlags.unaryUndertrained[i]) {
                baselineUnaryCount++;
                candidateUnaryMass += weight;
            }
            if (coverageFlags.pairUndertrained[i]) {
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
        // Every sample entering candidate fitting must carry its exact source
        // proposal.  Returning NaN fails the candidate closed instead of
        // silently mixing the production route with an obsolete q_m fallback.
        return Double.NaN;
    }

    private double computeEtaEnergy(CCDResult sample, EtaCorrections eta) {
        int[] conf = sample.conf;
        long[] pairKeys = sample.features == null
                ? null : sample.features.pairKeys;
        double energy = eta.globalOffsetKcal;
        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            if (rc < 0) continue;
            energy += eta.getOneBodyEta(pos, rc);
            for (int edgeIndex : interactionEdgeIndicesByFirstPosition[pos]) {
                int pos2 = interactionEdges[edgeIndex][1];
                int rc2 = conf[pos2];
                if (rc2 < 0) continue;
                long pairKey = pairKeys == null
                        ? EtaCorrections.packPairKey(pos, rc, pos2, rc2)
                        : pairKeys[edgeIndex];
                energy += eta.getPairEta(pairKey);
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
        for (int[] edge : interactionEdges) {
            int pos1 = edge[0];
            int pos2 = edge[1];
            int rc1 = conf[pos1];
            int rc2 = conf[pos2];
            if (rc1 < 0 || rc2 < 0) continue;
            if (coverage.pair(pos1, rc1, pos2, rc2)
                    < minimumCount) return true;
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
                raw.oneBody.length, numRCs, raw.pairSums.size());
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
        raw.pairSums.forEach((key, sum, count) -> {
            if (count <= 0) return;
            int contexts = coverage.pairDistinctContexts.getOrDefault(
                    key, 0);
            double factor = pairShrink == 0.0
                    ? 1.0 : contexts / (contexts + pairShrink);
            shrunk.pairSums.put(key, factor * sum, count);
        });
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
                long workload = frequencySeverityJointMomentLearning ? 0L
                        : PackStarTripleEtaCorrections.countCliqueAssignments(
                                rcs, interactionGraph);
                System.out.println(
                        "[PACK*-adaptive-frequency-severity] triple-eta preflight:"
                                + " cliqueAssignments=" + workload
                                + ", jointMomentLearning=" + frequencySeverityJointMomentLearning
                                + ", cap="
                                + frequencySeverityTripleEtaMaxAssignments
                                + ", maximumScale="
                                + frequencySeverityTripleEtaScale
                                + ", scaleGrid="
                                + Arrays.toString(
                                frequencySeverityTripleEtaScaleGrid));
                frequencySeverityTripleEta = frequencySeverityJointMomentLearning
                        ? PackStarTripleEtaCorrections.emptyFitted(rcs.getNumPos())
                        : PackStarTripleEtaCorrections.compute(
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
                + "; proposalLearning=" + frequencySeverityProposalLearning
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
                + "; alpha=" + (frequencySeverityProposalLearning ? "1-fixed" : "0-fixed-qm")
                + "; triple-model=" + (!frequencySeverityProposalLearning ? "disabled"
                : frequencySeverityBudgetForward ? "K-prefix-path" : "pair-only-vs-selected"));

        String trainStage = "adaptive-frequency-severity-qm-train";
        List<CCDResult> trainCCD = runParallelCCD(
                sampleConformationsFromDP(
                        trainSamples, stageRandom(trainStage)),
                branchMinimizingEmat, initialProposalLogZ,
                "adaptive-frequency-severity-qm");
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
                writeFrequencySeverityModelPathArtifact(artifactDir, candidates, artifactPrefix, roundSelection == null ? "NA" : roundSelection.id);
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
            if ((!frequencySeverityProposalLearning || round == 0) && roundSelection.alpha == 0.0
                    && !hasTripleEta(roundEta)) {
                roundEmat = branchMinimizingEmat;
                roundLogZ = initialProposalLogZ;
            } else {
                roundEmat = buildCorrectedEmat(roundEta);
                CorrectedDPResult correctedDP = recomputeDP(
                        roundEmat);
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
                null, Double.NaN, null);
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
                null, Double.NaN, null);
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
                null, Double.NaN, null);
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
        computeFunctionalObservable(finalCCD, finalLogR);
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

    private double computeFrequencySeverityTrainingLogMu(
            List<CCDResult> samples, EtaCorrections eta) {
        double[] logA = new double[samples.size()];
        double[] logV = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            CCDResult sample = samples.get(i);
            double targetEnergy = sample.eMin
                    + computeEtaEnergy(sample, eta);
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
            double etaEnergy = computeEtaEnergy(sample, eta);
            logRelative[i] = (sample.eMin + etaEnergy
                    - sample.eTrue) / RT - logMuTrain;
        }
        return logRelative;
    }

    /**
     * Evaluate the configured functional event on the fresh final sample.
     * The proposal normalizer is deliberately absent from logRelative: it is
     * common to every sample and cancels in the self-normalized event ratio.
     */
    private void computeFunctionalObservable(
            List<CCDResult> samples, double[] logRelative) {
        if (functionalEvent == null) {
            return;
        }
        if (samples.size() != logRelative.length) {
            throw new IllegalArgumentException(
                    "functional observable sample and weight counts differ");
        }

        functionalObservableResult = PackStarFunctionalObservableEvaluator.evaluate(
                functionalEventName, functionalEvent, RT, logRelative, i -> {
                    CCDResult sample = samples.get(i);
                    return new PackStarFunctionalEvent.Sample(
                            i, sample.conf, sample.eTrue, sample.eMin, sample.epmol);
                });

        PackStarFunctionalObservableResult result = functionalObservableResult;
        System.out.println("[PACK*-observable] event="
                + sanitizeTsv(result.getEventName())
                + ", status=" + result.getStatus()
                + ", P_C=" + formatObservableValue(result.getProbability())
                + ", G_restrict="
                + formatObservableValue(result.getRestrictionFreeEnergy())
                + ", samples=" + result.getSampleCount()
                + ", weightedSamples=" + result.getWeightedSampleCount()
                + ", hits=" + result.getHitCount()
                + ", ESS=" + formatObservableValue(
                result.getEffectiveSampleSize())
                + ", maxNormalizedWeight=" + formatObservableValue(
                result.getMaxNormalizedWeight())
                + ", diagnostic=" + sanitizeTsv(result.getDiagnostic()));
    }

    private static String formatObservableValue(double value) {
        return Double.isFinite(value)
                ? String.format(Locale.ROOT, "%.9g", value)
                : Double.toString(value);
    }

    private void setFrequencySeverityDiagnosticStatistics(
            List<CCDResult> samples, EtaCorrections eta,
            double[] logRelative,
            PackStarFrequencySeverityPAC.Interval interval) {
        double sumResidual = 0.0;
        double sumResidual2 = 0.0;
        for (CCDResult sample : samples) {
            double eEta = sample.eMin
                    + computeEtaEnergy(sample, eta);
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
        if (frequencySeverityOutputDir == null) return;
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
            try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
            writer.println("key\tvalue");
            writeFrequencySeverityKey(writer, "schema",
                    "packstar-adaptive-frequency-severity-protocol-v10");
            writeFrequencySeverityKey(writer, "proposalLearningEnabled",
                    frequencySeverityProposalLearning);
            writeFrequencySeverityKey(writer, "tripleFitThreads",
                    frequencySeverityBudgetForward ? PackStarTripleEtaCorrections.configuredFitThreads() : 1);
            writeFrequencySeverityKey(writer, "tripleMomentReduction",
                    "FP64-grouped-feature-masses-max-shifted-Kahan-cached-gradient-masses;pair-fit-sequential-logadd");
            writeFrequencySeverityKey(writer, "fixedProposalCalibration",
                    frequencySeverityProposalLearning ? "not-applicable"
                            : "same-train-discovery-validation-monitor-final-budgets;gauge-only-adaptation;zero-corrected-DP-sweeps");
            writeFrequencySeverityKey(writer, "proposalLearningObjective",
                    !frequencySeverityProposalLearning ? "none-fixed-qm"
                            : frequencySeverityJointMomentLearning
                            ? "source-weighted-joint-residual+nested-heldout-log-rho"
                            : "legacy-termwise+deduplicated-residual-mse");
            writeFrequencySeverityKey(writer, "algorithm",
                    "AdaptiveFrequencySeverityPAC-crossfit-selected-signed-triple-residual-source-aware-online-refit");
            writeFrequencySeverityKey(writer, "proposalModel",
                    !frequencySeverityProposalLearning ? "original-DP-proposal-no-learned-correction"
                            : frequencySeverityBudgetForward ? "fixed-pair-plus-jointly-refitted-forward-triple-prefixes-K-0-1-2-3"
                            : "pair-only-versus-legacy-terminal-triple-set");
            writeFrequencySeverityKey(writer, "tripleEtaEnabled",
                    frequencySeverityTripleEtaEnabled);
            writeFrequencySeverityKey(writer, "tripleEtaMaximumScale",
                    frequencySeverityTripleEtaScale);
            writeFrequencySeverityKey(writer, "tripleEtaScaleGrid",
                    Arrays.toString(frequencySeverityTripleEtaScaleGrid));
            writeFrequencySeverityKey(writer, "tripleEtaScaleSelection",
                    "fixed-one-when-included-not-continuously-tuned");
            writeFrequencySeverityKey(writer, "tripleEtaModelOrderSelection",
                    frequencySeverityBudgetForward ? "all-achieved-prefixes;crossfit-at-most-K-with-fold-early-stop"
                            : "pair-only-or-selected-triple-by-crossfit");
            writeFrequencySeverityKey(writer,
                    "tripleEtaAssignmentCapFallback",
                    frequencySeverityJointMomentLearning
                            ? "selected-table-budget;no-eager-clique-prior"
                            : "pair-only-when-clique-workload-exceeds-cap");
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumPartialAssignments",
                    frequencySeverityTripleEtaMaxAssignments);
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumSelectedPositionTriples",
                    frequencySeverityTripleEtaMaxPositionTriples);
            writeFrequencySeverityKey(writer,
                    "tripleEtaMaximumFillEdges",
                    frequencySeverityTripleEtaMaxFillEdges);
            writeFrequencySeverityKey(writer, "tripleEtaSelectionStrategy",
                    frequencySeverityBudgetForward ? "budget-forward"
                            : frequencySeverityDecompositionCostSelection ? "decomposition-cost" : "fill-edge");
            writeFrequencySeverityKey(writer, "tripleEtaSelectionCostMetric",
                    frequencySeverityBudgetForward ? "hard-resource-gates;heldout-log-rho-gain;absolute-1e-6-ties-by-work-host-file-table"
                            : frequencySeverityDecompositionCostSelection
                            ? "held-out-log-rho-gain/(1+max-relative-work-state-table-host-file-growth)"
                            : "gain/(1+new-fill-edges)");
            if (frequencySeverityDecompositionCostSelection) {
                writeFrequencySeverityKey(writer, "tripleEtaDecompositionLimits",
                        frequencySeverityDecompositionLimits.toString());
                writeFrequencySeverityKey(writer, "tripleEtaDecompositionPreview",
                        "weighted-hicks-configured-root;shared-with-proposal-execution;no-enumeration-arrays-or-dp-tables;per-pfunc-LRU-cache");
            }
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
                    frequencySeverityJointMomentLearning
                            ? "conditional-second-moment-shrunk-toward-global-moment;unseen-zero"
                            : "signed-shrinkage-prior-from-shared-three-pair-minimum-minus-independent-three-pair-minima");
            writeFrequencySeverityKey(writer,
                    "tripleEtaProposalFactor",
                    frequencySeverityJointMomentLearning
                            ? "nested-crossfit-conditional-second-moment-table"
                            : "crossfit-signed-Etrue-minus-Em-minus-single-pair-eta-residual-table");
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
                    frequencySeverityBudgetForward ? "minimum-projected-final-N;within-1.05-of-global-minimum-prefer-smaller-K-then-structure"
                            : "legacy-min-projected-final-N-then-epsilon-then-M2");
            writeFrequencySeverityKey(writer, "modelOrderFinalNTolerance", frequencySeverityBudgetForward ? 0.05 : 0.0);
            writeFrequencySeverityKey(writer, "tripleJointRefit", frequencySeverityBudgetForward);
            writeFrequencySeverityKey(writer, "tripleJointPenalty", "0.5*sum(mass*strength/min(contexts,sourceESS)*(h/RT)^2);unsupported-zero");
            writeFrequencySeverityKey(writer, "tripleJointMaxIterations", PackStarTripleEtaCorrections.JOINT_MOMENT_ITERATIONS);
            writeFrequencySeverityKey(writer, "tripleEarlyStopRule", "both-inner-fold-gains>1e-8;converged-refits;otherwise-stop");
            writeFrequencySeverityKey(writer, "tripleSelectionDataBoundary", "adaptation-only;post-freeze-validation-and-final-never-select-scopes");
            writeFrequencySeverityKey(writer, "proposalImprovementChangesSeverityPremise", false);
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
            File dir, String artifactPrefix,
            List<CCDResult> samples,
            FrequencySeverityEtaTraining training) {
        File sampleOutput = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "frequency_severity_training_samples.tsv"));
        try (PrintWriter writer = openFrequencySeverityWriter(sampleOutput)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(observationOutput)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(statisticsOutput)) {
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
        long[] pairKeys = raw.pairSums.keys();
        Arrays.sort(pairKeys);
        for (long key : pairKeys) {
            int count = raw.pairSums.count(key);
            if (count <= 0) continue;
            int pos1 = EtaCorrections.pairPos1(key);
            int rc1 = EtaCorrections.pairRc1(key);
            int pos2 = EtaCorrections.pairPos2(key);
            int rc2 = EtaCorrections.pairRc2(key);
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

    private void writeFrequencySeverityModelPathArtifact(File dir,
            List<FrequencySeverityCandidateScore> candidates, String prefix, String selectedId) {
        File output = newFrequencySeverityArtifact(dir,
                frequencySeverityPrefixedArtifactName(prefix, "eta_candidate_model_path.tsv"));
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
            writer.println("schema\tcandidate\tselected\teligible\tK\tfitFold0K\tfitFold1K"
                    + "\tinnerLogRhoGain\tpredictedFinalN\tstopReason\tfitFold0Stop\tfitFold1Stop\tscopes\tsignature\tstructure");
            for (FrequencySeverityCandidateScore c : candidates) {
                boolean active = hasTripleEta(c.fullEta);
                writer.printf(Locale.ROOT, "packstar-model-path-v1\t%s\t%s\t%s\t%d\t%d\t%d\t%.17g\t%d\t%s\t%s\t%s\t%s\t%s\t%s%n",
                        c.id, c.id.equals(selectedId), c.eligible, c.requestedTripleCount,
                        c.foldTripleCounts[0], c.foldTripleCounts[1], c.innerLogRhoGain,
                        c.crossfit.sizing.finalSamples, c.pathStopReason,
                        c.foldPathStopReasons[0], c.foldPathStopReasons[1],
                        active ? c.fullEta.tripleEta.positionTriples() : "[]",
                        active ? c.fullEta.tripleEta.signatureSha256 : "NA",
                        frequencySeverityBudgetForward && active ? candidateStructure(c) : "NA");
            }
        } catch (Exception ex) {
            throw new IllegalStateException("cannot write candidate model path", ex);
        }
    }

    private void writeFrequencySeverityCandidateArtifact(
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(dir,
                frequencySeverityPrefixedArtifactName(
                        artifactPrefix, "eta_candidate_scores.tsv"));
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_candidate_crossfit_pooled.tsv"));
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
            File dir, String artifactPrefix,
            List<FrequencySeverityCandidateScore> candidates,
            String selectedId) {
        File output = newFrequencySeverityArtifact(
                dir, frequencySeverityPrefixedArtifactName(
                        artifactPrefix,
                        "eta_candidate_refit_audit.tsv"));
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
            for (long key : raw.pairSums.keys()) {
                int pos1 = EtaCorrections.pairPos1(key);
                int rc1 = EtaCorrections.pairRc1(key);
                int pos2 = EtaCorrections.pairPos2(key);
                int rc2 = EtaCorrections.pairRc2(key);
                int count = raw.pairSums.count(key);
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
            writer.println("schema\tstage\tindex\tconf\teTrueKcal\teMinKcal"
                    + "\tetaKcal\teProposalKcal\tlogRawWeight"
                    + "\tlogMuTrain\tlogRelativeWeight\tbulkY\ttailEvent");
            for (int index = 0; index < samples.size(); index++) {
                CCDResult sample = samples.get(index);
                double etaEnergy = computeEtaEnergy(sample, eta);
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
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
            int trainCount, int discoveryCount, int validationCount,
            int monitorCount, int finalCount,
            int proposalDpSweeps, boolean valid) {
        File output = newFrequencySeverityArtifact(dir, "frequency_severity_final_interval.tsv");
        try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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
            writeFrequencySeverityKey(writer, "functionalEvent",
                    functionalEvent == null ? "NA" : functionalEventName);
            writeFrequencySeverityKey(writer, "functionalEventStatus",
                    functionalObservableResult.getStatus());
            writeFrequencySeverityKey(writer, "functionalEventProbability",
                    functionalObservableResult.getProbability());
            writeFrequencySeverityKey(writer, "functionalEventRestrictionFreeEnergy",
                    functionalObservableResult.getRestrictionFreeEnergy());
            writeFrequencySeverityKey(writer, "functionalEventHitCount",
                    functionalObservableResult.getHitCount());
            writeFrequencySeverityKey(writer, "functionalEventEffectiveSampleSize",
                    functionalObservableResult.getEffectiveSampleSize());
            writeFrequencySeverityKey(writer, "functionalEventMaxNormalizedWeight",
                    functionalObservableResult.getMaxNormalizedWeight());
            writeFrequencySeverityKey(writer, "totalCcd", totalCCDCalls);
            writeFrequencySeverityKey(writer, "totalSampleRecords",
                    totalCCDSampleRecords);
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
            try (PrintWriter writer = openFrequencySeverityWriter(output)) {
                writer.println("key\tvalue");
                writeFrequencySeverityKey(writer, "schema", "packstar-adaptive-frequency-severity-failure-v1");
                writeFrequencySeverityKey(writer, "reason", sanitizeTsv(reason));
                writeFrequencySeverityKey(writer, "trainCcd", trainCount);
                writeFrequencySeverityKey(writer, "pilotCcd", pilotCount);
                writeFrequencySeverityKey(writer, "monitorCcd", monitorCount);
                writeFrequencySeverityKey(writer, "finalCcd", finalCount);
                writeFrequencySeverityKey(writer, "totalCcd", totalCCDCalls);
                writeFrequencySeverityKey(writer, "totalSampleRecords",
                        totalCCDSampleRecords);
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
            try (PrintWriter writer = openFrequencySeverityWriter(output)) {
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

            try (PrintWriter summaryWriter = openFrequencySeverityWriter(summary);
                 PrintWriter termWriter = openFrequencySeverityWriter(terms)) {

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

    /** Open artifact output with a real byte buffer and create-new semantics. */
    private static PrintWriter openFrequencySeverityWriter(File output)
            throws IOException {
        return new PrintWriter(Files.newBufferedWriter(
                output.toPath(), StandardCharsets.UTF_8,
                StandardOpenOption.CREATE_NEW,
                StandardOpenOption.WRITE));
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

    /** Replace tabs and line breaks before writing a TSV artifact field. */
    private static String sanitizeTsv(String value) {
        if (value == null) return "";
        return value.replace('\t', ' ')
                .replace('\n', ' ')
                .replace('\r', ' ');
    }

    private static boolean isXtmpOutputPath(String value) {
        return value != null && !value.isEmpty()
                && new File(value).isAbsolute()
                && value.startsWith("/usr/xtmp/lz280/");
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
    private List<int[]> sampleConformationsFromDP(int n, Random rng) {
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

        List<RootedTreeEdge> topDownOrder = getSamplingTopDownOrder();
        ExecutorService samplingPool = getSamplingPool();

        int edgeOrdinal = 0;
        for (RootedTreeEdge edge : topDownOrder) {
            if (!edge.getIsLambdaEdge()) continue;

            long edgeStart = System.currentTimeMillis();
            Map<Long, SampleGroup> groups = groupSamplesByMIdx(edge, confs);
            long edgeSeed = mix64(masterSeed + 0xD1B54A32D192ED03L * (long) edgeOrdinal++);
            sampleEdgeGroups(edge, confs, sampleRngs, groups, samplingPool, edgeSeed);
            logSamplingEdgeProgress(edge, n, groups, System.currentTimeMillis() - edgeStart);
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
            SampleGroup group = groups.get(mIdx);
            if (group == null) {
                group = new SampleGroup(mIdx);
                group.mRCs = Arrays.copyOf(mLocalRCs, mLocalRCs.length);
                groups.put(mIdx, group);
            }
            group.add(s);
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
        List<SampleGroup> singletonGroups = new ArrayList<>();
        for (SampleGroup group : groups) {
            if (group.count == 1) singletonGroups.add(group);
        }
        runSamplingGroupsParallel(singletonGroups, group -> {
                int sampleIndex = group.sampleIndices[0];
                int lIdx = sampleLambdaStateStreaming(
                        edge, group.mIdx, group.mRCs,
                        sampleRngs[sampleIndex]);
                writeLambdaState(edge, lIdx, confs[sampleIndex]);
            }, samplingPool);
    }

    private void processSmallGroupsParallel(RootedTreeEdge edge,
                                            Collection<SampleGroup> groups,
                                            int[][] confs,
                                            SplittableRandom[] sampleRngs,
                                            ExecutorService samplingPool) {
        runSamplingGroupsParallel(groups,
                group -> sampleGroupSerial(edge, group, confs, sampleRngs),
                samplingPool);
    }

    /**
     * Submit at most one task per sampling worker.  The old implementation
     * submitted one future per M-group, which made a high-cardinality edge
     * build a large executor queue even though groups are independent.
     */
    private void runSamplingGroupsParallel(Collection<SampleGroup> groups,
        Consumer<SampleGroup> task,
                                           ExecutorService samplingPool) {
        if (groups.isEmpty()) return;
        List<SampleGroup> work = new ArrayList<>(groups);
        int workers = Math.min(samplingThreads, work.size());
        AtomicInteger next = new AtomicInteger(0);
        List<Future<?>> futures = new ArrayList<>(workers);
        for (int worker = 0; worker < workers; worker++) {
            futures.add(samplingPool.submit(() -> {
                int index;
                while ((index = next.getAndIncrement()) < work.size()) {
                    task.accept(work.get(index));
                }
            }));
        }
        waitForSamplingFutures(futures);
    }

    private void sampleGroupSerial(RootedTreeEdge edge, SampleGroup group,
                                   int[][] confs, SplittableRandom[] sampleRngs) {
        if (group.count == 1) {
            int sampleIndex = group.sampleIndices[0];
            int lIdx = sampleLambdaStateStreaming(
                    edge, group.mIdx, group.mRCs,
                    sampleRngs[sampleIndex]);
            writeLambdaState(edge, lIdx, confs[sampleIndex]);
            return;
        }

        double[] cdf = buildConditionalCDFSerial(
                edge, group.mIdx, group.mRCs);
        sampleGroupFromCDF(edge, group, confs, sampleRngs, cdf);
    }

    private void sampleGroupWithParallelCDF(RootedTreeEdge edge, SampleGroup group,
                                            int[][] confs, SplittableRandom[] sampleRngs,
                                            ExecutorService samplingPool) {
        double[] cdf = buildConditionalCDFParallel(
                edge, group.mIdx, group.mRCs, samplingPool);
        sampleGroupFromCDF(edge, group, confs, sampleRngs, cdf);
    }

    private void sampleGroupFromCDF(RootedTreeEdge edge, SampleGroup group,
                                    int[][] confs, SplittableRandom[] sampleRngs,
                                    double[] cdf) {
        double total = cdf[cdf.length - 1];
        if (!(total > 0.0) || !Double.isFinite(total)) {
            for (int i = 0; i < group.count; i++) {
                int sampleIndex = group.sampleIndices[i];
                int lIdx = sampleLambdaStateStreaming(
                        edge, group.mIdx, group.mRCs,
                        sampleRngs[sampleIndex]);
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

    private double[] buildConditionalCDFSerial(
            RootedTreeEdge edge, long mIdx, int[] mRCs) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] cdf = acquireSamplingCdf(totalLambda);

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

    private double[] buildConditionalCDFParallel(
            RootedTreeEdge edge, long mIdx, int[] mRCs,
            ExecutorService samplingPool) {
        int totalLambda = edge.getTotalLambdaStates();
        double[] cdf = acquireSamplingCdf(totalLambda);
        int chunks = Math.min(samplingThreads, totalLambda);
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

    private double[] acquireSamplingCdf(int length) {
        double[] cdf = samplingCdfScratch.get();
        // Callers use the array length as the lambda-state count, so retaining
        // a larger buffer across differently sized edges would be incorrect.
        if (cdf.length != length) {
            cdf = new double[length];
            samplingCdfScratch.set(cdf);
        }
        return cdf;
    }

    private int sampleLambdaStateStreaming(
            RootedTreeEdge edge, long mIdx, int[] mRCs,
            SplittableRandom rng) {
        int totalLambda = edge.getTotalLambdaStates();
        if (totalLambda <= 1) {
            return 0;
        }

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
        return edge.computeLambdaLogWeightPublic(
                mIdx, mRCs, lIdx, activeProposalEmat,
                proposalInteractionGraph, RT);
    }

    private void writeLambdaState(RootedTreeEdge edge, int lIdx, int[] conf) {
        edge.writeLambdaStateGlobal(lIdx, conf);
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
     * Build top-down traversal order for the compact tree.
     */
    private List<RootedTreeEdge> buildTopDownOrder(RootedTreeEdge rootEdge) {
        List<RootedTreeEdge> order = new ArrayList<>();
        Queue<RootedTreeEdge> queue = new ArrayDeque<>();
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

    private List<RootedTreeEdge> getSamplingTopDownOrder() {
        if (samplingTopDownOrder == null || samplingOrderRoot != rootedRootEdge) {
            samplingTopDownOrder = List.copyOf(buildTopDownOrder(rootedRootEdge));
            samplingOrderRoot = rootedRootEdge;
        }
        return samplingTopDownOrder;
    }

    private ExecutorService getSamplingPool() {
        if (samplingThreads <= 1) return null;
        if (samplingResourcesClosed) {
            throw new IllegalStateException("PACK* sampling resources already closed");
        }
        if (samplingPool == null) {
            samplingPool = Executors.newFixedThreadPool(samplingThreads,
                    daemonThreadFactory("packstar-sample"));
        }
        return samplingPool;
    }

    private void shutdownSamplingResources() {
        if (samplingResourcesClosed) return;
        samplingResourcesClosed = true;
        if (samplingPool != null) {
            samplingPool.shutdownNow();
            samplingPool = null;
        }
        samplingCdfScratch.remove();
    }

    // ========== Phase 2: Parallel CCD minimization ==========

    /**
     * Immutable key for an exact discrete conformation assignment.  Using a
     * value key rather than {@code int[]} identity is important here because
     * DP sampling allocates a fresh array for every draw.  The assignment is
     * copied once into the key, so the map remains safe even if a caller later
     * reuses or mutates its input array.
     */
    private static final class ConformationKey {
        private final int[] values;
        private final int hash;

        ConformationKey(int[] values) {
            if (values == null) {
                throw new IllegalArgumentException(
                        "CCD conformation cannot be null");
            }
            this.values = values.clone();
            this.hash = Arrays.hashCode(this.values);
        }

        @Override
        public boolean equals(Object other) {
            return other instanceof ConformationKey
                    && Arrays.equals(values,
                    ((ConformationKey) other).values);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }

    /**
     * Unique CCD work plus the inverse map needed to scatter results back to
     * the original IID sample order.  The returned sample list deliberately
     * keeps duplicate references to one CCDResult: all later reductions still
     * visit every draw, while feature caches can also avoid recomputing the
     * same forcefield decomposition for a repeated assignment.
     */
    static final class UniqueConformationBatch {
        final List<int[]> uniqueConformations;
        final int[] uniqueIndexBySample;

        UniqueConformationBatch(List<int[]> uniqueConformations,
                                int[] uniqueIndexBySample) {
            this.uniqueConformations = uniqueConformations;
            this.uniqueIndexBySample = uniqueIndexBySample;
        }
    }

    static UniqueConformationBatch deduplicateConformations(
            List<int[]> conformations) {
        if (conformations == null) {
            throw new IllegalArgumentException(
                    "CCD conformation list cannot be null");
        }

        Map<ConformationKey, Integer> indexByKey =
                new LinkedHashMap<>(Math.max(16, conformations.size()));
        List<int[]> unique = new ArrayList<>(conformations.size());
        int[] uniqueIndexBySample = new int[conformations.size()];

        for (int sample = 0; sample < conformations.size(); sample++) {
            int[] conformation = conformations.get(sample);
            ConformationKey key = new ConformationKey(conformation);
            Integer uniqueIndex = indexByKey.get(key);
            if (uniqueIndex == null) {
                uniqueIndex = unique.size();
                indexByKey.put(key, uniqueIndex);
                unique.add(conformation.clone());
            }
            uniqueIndexBySample[sample] = uniqueIndex;
        }

        return new UniqueConformationBatch(
                List.copyOf(unique), uniqueIndexBySample);
    }

    private static class CCDResult {
        final int[] conf;
        final double eTrue;
        final double eMin;
        final EnergyCalculator.EnergiedParametricMolecule epmol;
        final double sourceProposalEnergy;
        final double sourceProposalLogZ;
        final String sourceProposalId;
        SampleFeatures features;
        String assignmentKey;

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
     * Eta-training features for one minimized conformation.  These are raw
     * forcefield-minus-base-EMAT term corrections and do not depend on the
     * fitted eta, shrinkage choice, or fold.  Keeping them on CCDResult makes
     * every later refit an aggregation pass instead of another forcefield
     * decomposition pass.
     */
    private static class SampleFeatures {
        final double[] oneBodyCorrections;
        final double[] pairCorrections;
        /** Pair-cell keys in the estimator's stable edge order. */
        final long[] pairKeys;

        SampleFeatures(double[] oneBodyCorrections,
                       double[] pairCorrections,
                       long[] pairKeys) {
            this.oneBodyCorrections = oneBodyCorrections;
            this.pairCorrections = pairCorrections;
            this.pairKeys = pairKeys;
        }
    }

    /**
     * Run CCD while retaining the exact proposal provenance of this batch.
     * The source energy and log normalizer are later used only by adaptation
     * to reweight samples to a selected proposal; monitor/final remain fresh.
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
        UniqueConformationBatch uniqueBatch =
                deduplicateConformations(conformations);
        int total = conformations.size();
        int uniqueTotal = uniqueBatch.uniqueConformations.size();

        // Keep the estimator's external CCD count honest: it is the number of
        // actual minimizer invocations, not the number of IID records that
        // consume their results.  The latter is retained separately for audit
        // files and for checking the requested sample count.
        totalCCDSampleRecords += total;
        totalCCDCalls += uniqueTotal;

        if (uniqueTotal < total) {
            System.out.println("[PACK*-CCD] exact conformation dedup: samples="
                    + total + ", unique=" + uniqueTotal
                    + ", reused=" + (total - uniqueTotal)
                    + ", saved=" + String.format(Locale.ROOT, "%.2f%%",
                    total > 0 ? 100.0 * (total - uniqueTotal) / total : 0.0));
        }

        AtomicReferenceArray<CCDResult> results =
                new AtomicReferenceArray<>(uniqueTotal);
        AtomicInteger completed = new AtomicInteger(0);

        // Instrumentation is useful for tuning, but sorting a timing array for
        // every CCD batch is production overhead.  Keep it opt-in.
        long batchStartNs = ccdInstrumentation ? System.nanoTime() : 0L;
        long[] ccdWallUs = ccdInstrumentation ? new long[uniqueTotal] : null;
        Set<String> workerThreads = ccdInstrumentation
                ? java.util.concurrent.ConcurrentHashMap.newKeySet() : null;
        AtomicInteger inFlight = ccdInstrumentation
                ? new AtomicInteger(0) : null;
        AtomicInteger peakConcurrency = ccdInstrumentation
                ? new AtomicInteger(0) : null;
        int pendingTasks = 0;

        for (int idx = 0; idx < uniqueTotal; idx++) {
            int[] conf = uniqueBatch.uniqueConformations.get(idx);
            final int ci = idx;
            final long submitNs = System.nanoTime();
            ResidueInteractions inters = makeSparseFullConfInters(conf);
            RCTuple tuple = new RCTuple(conf);

            minimizingEcalc.calcEnergyAsync(tuple, inters, epmol -> {
                if (ccdInstrumentation) {
                    int cur = inFlight.incrementAndGet();
                    peakConcurrency.accumulateAndGet(cur, Math::max);
                    workerThreads.add(Thread.currentThread().getName());
                }

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
                if (ccdInstrumentation) {
                    ccdWallUs[ci] = (System.nanoTime() - submitNs) / 1000L;
                }

                if (ccdInstrumentation) {
                    inFlight.decrementAndGet();
                }
                int done = completed.incrementAndGet();
                if (done % 100 == 0 || done == uniqueTotal) {
                    System.out.println("[PACK*] CCD progress: " + done + "/"
                            + uniqueTotal + " unique (logical samples="
                            + total + ")");
                }
            });

            if (++pendingTasks >= ccdSubmissionBatchSize) {
                minimizingEcalc.tasks.waitForFinish();
                pendingTasks = 0;
            }
        }

        if (pendingTasks > 0) {
            minimizingEcalc.tasks.waitForFinish();
        }

        if (ccdInstrumentation) {
            long batchWallMs = (System.nanoTime() - batchStartNs) / 1_000_000L;
            long[] sorted = java.util.Arrays.copyOf(ccdWallUs, uniqueTotal);
            java.util.Arrays.sort(sorted);
            long sumUs = 0L;
            for (long w : ccdWallUs) sumUs += w;
            double p50 = uniqueTotal > 0
                    ? sorted[Math.min(uniqueTotal - 1, uniqueTotal / 2)] / 1000.0 : 0.0;
            double p95 = uniqueTotal > 0
                    ? sorted[Math.min(uniqueTotal - 1, (int) (uniqueTotal * 0.95))]
                    / 1000.0 : 0.0;
            double maxMs = uniqueTotal > 0
                    ? sorted[uniqueTotal - 1] / 1000.0 : 0.0;
            double sumMs = sumUs / 1000.0;
            double effParallel = batchWallMs > 0 ? sumMs / batchWallMs : 0.0;
            System.out.printf(
                    "[PACK*-INSTR] CCD batch: n=%d batchWall=%dms "
                            + "logicalSamples=%d "
                            + "distinctThreads=%d peakConcurrency=%d "
                            + "perCCD(ms) p50=%.1f p95=%.1f max=%.1f "
                            + "sumCCD=%.0fms effParallel=%.1fx%n",
                    uniqueTotal, batchWallMs, total, workerThreads.size(),
                    peakConcurrency.get(), p50, p95, maxMs, sumMs,
                    effParallel);
        }

        List<CCDResult> orderedResults = new ArrayList<>(total);
        for (int i = 0; i < total; i++) {
            CCDResult result = results.get(uniqueBatch.uniqueIndexBySample[i]);
            if (result == null) {
                throw new IllegalStateException(
                        "CCD result missing for unique conformation "
                                + uniqueBatch.uniqueIndexBySample[i]
                                + " (sample " + i + ")");
            }
            orderedResults.add(result);
        }
        return orderedResults;
    }

    private SampleFeatures getSampleFeatures(CCDResult result) {
        if (result.features == null) {
            result.features = buildSampleFeatures(result);
        }
        return result.features;
    }

    private static String getAssignmentKey(CCDResult result) {
        if (result.assignmentKey == null) {
            result.assignmentKey = Arrays.toString(result.conf);
        }
        return result.assignmentKey;
    }

    private SampleFeatures buildSampleFeatures(CCDResult result) {
        int[] conf = result.conf;
        ResidueForcefieldEnergy efunc =
                (ResidueForcefieldEnergy) minimizingEcalc.ecalc
                        .makeEnergyFunction(result.epmol);
        double[] oneBody = new double[conf.length];
        Arrays.fill(oneBody, Double.NaN);
        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            if (rc < 0) continue;
            oneBody[pos] = efunc.makeSubset(
                    getSingleInteractionTemplate(pos, rc)).getEnergy()
                    - branchMinimizingEmat.getOneBody(pos, rc);
        }

        double[] pair = new double[interactionEdges.length];
        long[] pairKeys = new long[interactionEdges.length];
        Arrays.fill(pairKeys, MISSING_PAIR_KEY);
        Arrays.fill(pair, Double.NaN);
        for (int edgeIndex = 0; edgeIndex < interactionEdges.length;
             edgeIndex++) {
            int pos1 = interactionEdges[edgeIndex][0];
            int pos2 = interactionEdges[edgeIndex][1];
            int rc1 = conf[pos1];
            int rc2 = conf[pos2];
            if (rc1 < 0 || rc2 < 0) continue;
            pairKeys[edgeIndex] = EtaCorrections.packPairKey(
                    pos1, rc1, pos2, rc2);
            pair[edgeIndex] = efunc.makeSubset(
                    getPairInteractionTemplate(pos1, rc1, pos2, rc2))
                    .getEnergy()
                    - branchMinimizingEmat.getPairwise(
                    pos1, rc1, pos2, rc2);
        }
        return new SampleFeatures(oneBody, pair, pairKeys);
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
        @FunctionalInterface
        interface PairStatsConsumer {
            void accept(long key, double sum, int count);
        }

        /** Primitive long-to-int map for distinct-context coverage counts. */
        static class PairCountMap {
            private static final double LOAD_FACTOR = 0.70;
            private long[] keys;
            private int[] values;
            private boolean[] used;
            private int mask;
            private int threshold;
            private int size;

            PairCountMap() {
                allocate(16);
            }

            void increment(long key) {
                int slot = findSlot(key);
                if (!used[slot]) {
                    ensureCapacityForInsert();
                    slot = findSlot(key);
                    used[slot] = true;
                    keys[slot] = key;
                    size++;
                }
                if (values[slot] == Integer.MAX_VALUE) {
                    throw new ArithmeticException(
                            "eta context count overflow");
                }
                values[slot]++;
            }

            int getOrDefault(long key, int defaultValue) {
                int slot = findSlot(key);
                return used[slot] ? values[slot] : defaultValue;
            }

            private void allocate(int capacity) {
                keys = new long[capacity];
                values = new int[capacity];
                used = new boolean[capacity];
                mask = capacity - 1;
                threshold = Math.max(1,
                        (int) Math.floor(capacity * LOAD_FACTOR));
                size = 0;
            }

            private void ensureCapacityForInsert() {
                if (size + 1 <= threshold) return;
                long[] oldKeys = keys;
                int[] oldValues = values;
                boolean[] oldUsed = used;
                if (oldKeys.length > (1 << 29)) {
                    throw new IllegalStateException(
                            "eta context table exceeds supported capacity");
                }
                allocate(oldKeys.length << 1);
                for (int oldSlot = 0; oldSlot < oldUsed.length; oldSlot++) {
                    if (!oldUsed[oldSlot]) continue;
                    int slot = findSlot(oldKeys[oldSlot]);
                    used[slot] = true;
                    keys[slot] = oldKeys[oldSlot];
                    values[slot] = oldValues[oldSlot];
                    size++;
                }
            }

            private int findSlot(long key) {
                int slot = PairStatsMap.mix(key) & mask;
                while (used[slot] && keys[slot] != key) {
                    slot = (slot + 1) & mask;
                }
                return slot;
            }
        }

        /**
         * Small open-addressed map for the sparse pair cells.  Pair eta is
         * updated in a single estimator thread, so a primitive table is enough
         * and avoids both boxed Long keys and one double[2] allocation per
         * observed cell.  Iteration order is intentionally unspecified; all
         * serialized views sort their keys when order matters.
         */
        static class PairStatsMap {
            private static final double LOAD_FACTOR = 0.70;
            private long[] keys;
            private double[] sums;
            private int[] counts;
            private boolean[] used;
            private int mask;
            private int threshold;
            private int size;

            PairStatsMap() {
                this(16);
            }

            PairStatsMap(int requestedCapacity) {
                int capacity = 1;
                while (capacity < Math.max(4, requestedCapacity)) {
                    capacity <<= 1;
                }
                allocate(capacity);
            }

            int size() {
                return size;
            }

            void add(long key, double value) {
                int slot = findSlot(key);
                if (!used[slot]) {
                    ensureCapacityForInsert();
                    slot = findSlot(key);
                    used[slot] = true;
                    keys[slot] = key;
                    size++;
                }
                sums[slot] += value;
                if (counts[slot] == Integer.MAX_VALUE) {
                    throw new ArithmeticException(
                            "eta pair observation count overflow");
                }
                counts[slot]++;
            }

            void put(long key, double sum, int count) {
                if (count < 0) {
                    throw new IllegalArgumentException(
                            "eta pair observation count must be nonnegative");
                }
                int slot = findSlot(key);
                if (!used[slot]) {
                    ensureCapacityForInsert();
                    slot = findSlot(key);
                    used[slot] = true;
                    keys[slot] = key;
                    size++;
                }
                sums[slot] = sum;
                counts[slot] = count;
            }

            int count(long key) {
                int slot = findSlot(key);
                return used[slot] ? counts[slot] : 0;
            }

            double sum(long key) {
                int slot = findSlot(key);
                return used[slot] ? sums[slot] : 0.0;
            }

            double mean(long key) {
                int slot = findSlot(key);
                if (!used[slot] || counts[slot] == 0) return 0.0;
                return sums[slot] / counts[slot];
            }

            long[] keys() {
                long[] result = new long[size];
                int out = 0;
                for (int slot = 0; slot < used.length; slot++) {
                    if (used[slot]) result[out++] = keys[slot];
                }
                return result;
            }

            void forEach(PairStatsConsumer consumer) {
                for (int slot = 0; slot < used.length; slot++) {
                    if (used[slot]) {
                        consumer.accept(keys[slot], sums[slot], counts[slot]);
                    }
                }
            }

            private void allocate(int capacity) {
                keys = new long[capacity];
                sums = new double[capacity];
                counts = new int[capacity];
                used = new boolean[capacity];
                mask = capacity - 1;
                threshold = Math.max(1,
                        (int) Math.floor(capacity * LOAD_FACTOR));
                size = 0;
            }

            private void ensureCapacityForInsert() {
                if (size + 1 <= threshold) return;
                long[] oldKeys = keys;
                double[] oldSums = sums;
                int[] oldCounts = counts;
                boolean[] oldUsed = used;
                if (oldKeys.length > (1 << 29)) {
                    throw new IllegalStateException(
                            "eta pair table exceeds supported capacity");
                }
                allocate(oldKeys.length << 1);
                for (int oldSlot = 0; oldSlot < oldUsed.length; oldSlot++) {
                    if (!oldUsed[oldSlot]) continue;
                    int slot = findSlot(oldKeys[oldSlot]);
                    used[slot] = true;
                    keys[slot] = oldKeys[oldSlot];
                    sums[slot] = oldSums[oldSlot];
                    counts[slot] = oldCounts[oldSlot];
                    size++;
                }
            }

            private int findSlot(long key) {
                int slot = mix(key) & mask;
                while (used[slot] && keys[slot] != key) {
                    slot = (slot + 1) & mask;
                }
                return slot;
            }

            private static int mix(long value) {
                value ^= value >>> 33;
                value *= 0xff51afd7ed558ccdl;
                value ^= value >>> 33;
                value *= 0xc4ceb9fe1a85ec53l;
                value ^= value >>> 33;
                return (int) value;
            }
        }

        // eta_oneBody[pos][rc] = mean correction for one-body term
        final double[][] oneBody;
        final int[][] oneBodyCounts; // per-(pos,rc) counts, NOT per-pos
        // eta_pair[pos1][rc1][pos2][rc2] = mean correction for pair term
        // Stored as map to avoid sparse array issues
        final PairStatsMap pairSums;
        int oneBodyCount = 0;
        int pairCount = 0;
        // Added once per complete conformation.  This is an energy gauge: it
        // changes Z_eta and residual weights by reciprocal constants but leaves
        // the normalized proposal q_eta exactly unchanged.
        double globalOffsetKcal = 0.0;
        PackStarTripleEtaCorrections tripleEta = null;
        double tripleEtaScale = 0.0;

        EtaCorrections(int numPos, int[] numRCsPerPos) {
            this(numPos, numRCsPerPos, 16);
        }

        EtaCorrections(int numPos, int[] numRCsPerPos,
                       int expectedPairCells) {
            oneBody = new double[numPos][];
            oneBodyCounts = new int[numPos][];
            for (int p = 0; p < numPos; p++) {
                oneBody[p] = new double[numRCsPerPos[p]];
                oneBodyCounts[p] = new int[numRCsPerPos[p]];
            }
            pairSums = new PairStatsMap(expectedPairCells);
        }

        private EtaCorrections(
                EtaCorrections source,
                PackStarTripleEtaCorrections tripleEta,
                double tripleEtaScale) {
            oneBody = source.oneBody;
            oneBodyCounts = source.oneBodyCounts;
            pairSums = source.pairSums;
            oneBodyCount = source.oneBodyCount;
            pairCount = source.pairCount;
            globalOffsetKcal = source.globalOffsetKcal;
            this.tripleEta = tripleEta;
            this.tripleEtaScale = tripleEtaScale;
        }

        static long packPairKey(int pos1, int rc1, int pos2, int rc2) {
            // Ensure pos1 < pos2
            if (pos1 > pos2) {
                int tmp = pos1; pos1 = pos2; pos2 = tmp;
                tmp = rc1; rc1 = rc2; rc2 = tmp;
            }
            return ((long) pos1 << 48) | ((long) rc1 << 32) | ((long) pos2 << 16) | rc2;
        }

        static int pairPos1(long key) {
            return (int) ((key >>> 48) & 0xffffL);
        }

        static int pairRc1(long key) {
            return (int) ((key >>> 32) & 0xffffL);
        }

        static int pairPos2(long key) {
            return (int) ((key >>> 16) & 0xffffL);
        }

        static int pairRc2(long key) {
            return (int) (key & 0xffffL);
        }

        void addOneBodySample(int pos, int rc, double correction) {
            oneBody[pos][rc] += correction;
            oneBodyCounts[pos][rc]++;
        }

        void addPairSample(int pos1, int rc1, int pos2, int rc2, double correction) {
            long key = packPairKey(pos1, rc1, pos2, rc2);
            pairSums.add(key, correction);
        }

        double getOneBodyEta(int pos, int rc) {
            if (oneBodyCounts[pos][rc] == 0) return 0.0;
            return oneBody[pos][rc] / oneBodyCounts[pos][rc];
        }

        double getPairEta(int pos1, int rc1, int pos2, int rc2) {
            long key = packPairKey(pos1, rc1, pos2, rc2);
            return getPairEta(key);
        }

        double getPairEta(long key) {
            return pairSums.mean(key);
        }

        /**
         * Make a read-only model-order view over a fitted pair correction.
         * Candidate scoring never mutates an EtaCorrections after fitting, so
         * sharing the dense unary arrays and primitive pair table avoids a
         * second full copy just to attach a selected triple residual.
         */
        private EtaCorrections withTripleModel(
                PackStarTripleEtaCorrections prior, double scale) {
            return new EtaCorrections(this, prior, scale);
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
            return withTripleModel(prior, scale);
        }
    }

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
            allFitSamples.add(result);
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                if (historical || fold == newestFold) {
                    foldFitSampleCounts[fold]++;
                    foldFitSamples[fold].add(result);
                }
            }
            int[] conf = result.conf;
            if (frequencySeverityProposalLearning) {
            SampleFeatures features = getSampleFeatures(result);

            for (int pos = 0; pos < numPos; pos++) {
                int rc = conf[pos];
                if (rc < 0) continue;
                double correction = features.oneBodyCorrections[pos];
                all.addOneBodySample(pos, rc, correction);
                for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                    if (historical || fold == newestFold) {
                        folds[fold].addOneBodySample(pos, rc, correction);
                    }
                }
                observations.add(new FrequencySeverityEtaCellObservation(
                        sampleIndex, newestFold, false,
                        pos, rc, -1, -1, correction));
            }
            for (int edgeIndex = 0; edgeIndex < interactionEdges.length;
                 edgeIndex++) {
                int pos1 = interactionEdges[edgeIndex][0];
                int pos2 = interactionEdges[edgeIndex][1];
                int rc1 = conf[pos1];
                int rc2 = conf[pos2];
                if (rc1 < 0 || rc2 < 0) continue;
                double correction = features.pairCorrections[edgeIndex];
                all.addPairSample(pos1, rc1, pos2, rc2, correction);
                for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                    if (historical || fold == newestFold) {
                        folds[fold].addPairSample(
                                pos1, rc1, pos2, rc2, correction);
                    }
                }
                observations.add(new FrequencySeverityEtaCellObservation(
                        sampleIndex, newestFold, true,
                        pos1, rc1, pos2, rc2, correction));
            }
            }

            String assignment = getAssignmentKey(result);
            if (allAssignments.add(assignment)) {
                addFrequencySeverityCoverage(conf, allCoverage);
            }
            for (int fold = 0; fold < frequencySeverityFolds; fold++) {
                if ((historical || fold == newestFold)
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
        for (int[] edge : interactionEdges) {
            int pos1 = edge[0];
            int pos2 = edge[1];
            int rc1 = conf[pos1];
            int rc2 = conf[pos2];
            if (rc1 < 0 || rc2 < 0) continue;
            long key = EtaCorrections.packPairKey(pos1, rc1, pos2, rc2);
            coverage.pairDistinctContexts.increment(key);
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
    // ========== Phase 4: Build corrected emat and recompute DP ==========

    /** Build E_corrected = E_m + eta. */
    private EnergyMatrix buildCorrectedEmat(EtaCorrections eta) {
        ensureProposalInteractionGraph(eta);
        boolean hasTriple = eta != null
                && eta.tripleEta != null
                && eta.tripleEtaScale != 0.0;
        EnergyMatrix corrected = hasTriple
                ? new EnergyMatrix(getProposalMinimizingEmat())
                : beginPairwiseCorrectedEmat();
        int numPos = rcs.getNumPos();

        // Add the energy gauge to position zero.  This is the only dense
        // update: a constant must reach every conformation in order to keep
        // the proposal normalization explicit.
        if (eta.globalOffsetKcal != 0.0) {
            int numRC = corrected.getNumConfAtPos(0);
            for (int rc = 0; rc < numRC; rc++) {
                corrected.setOneBody(
                        0, rc,
                        corrected.getOneBody(0, rc)
                                + eta.globalOffsetKcal);
            }
        }

        // Eta is observed only on a sparse subset of RC cells.  Updating all
        // one-body cells and all graph RC pairs made every refit pay for the
        // full EMAT even when only a small training support changed.
        if (!hasTriple) {
            ensureCorrectedOneBodyKeyCapacity(eta.oneBodyCount);
        }
        for (int pos = 0; pos < numPos; pos++) {
            for (int rc = 0; rc < eta.oneBody[pos].length; rc++) {
                if (eta.oneBodyCounts[pos][rc] == 0) continue;
                if (!hasTriple) {
                    correctedPairwiseOneBodyKeys[
                            correctedPairwiseOneBodyCount++] =
                            ((long) pos << 32) | (rc & 0xffffffffL);
                }
                corrected.setOneBody(
                        pos, rc,
                        corrected.getOneBody(pos, rc)
                                + eta.getOneBodyEta(pos, rc));
            }
        }

        if (!hasTriple) {
            ensureCorrectedPairKeyCapacity(eta.pairSums.size());
        }
        eta.pairSums.forEach((key, sum, count) -> {
            if (count <= 0) return;
            int pos1 = EtaCorrections.pairPos1(key);
            int rc1 = EtaCorrections.pairRc1(key);
            int pos2 = EtaCorrections.pairPos2(key);
            int rc2 = EtaCorrections.pairRc2(key);
            if (!hasTriple) {
                correctedPairwisePairKeys[
                        correctedPairwisePairCount++] = key;
            }
            corrected.setPairwise(
                    pos1, rc1, pos2, rc2,
                    corrected.getPairwise(pos1, rc1, pos2, rc2)
                            + sum / count);
        });

        if (hasTriple) {
            eta.tripleEta.applyResidualTo(
                    corrected, eta.tripleEtaScale,
                    eta::getPairEta);
        }

        // Proposal-only fill pairs were zeroed when the graph-specific base
        // minimizing matrix was first built.  Eta updates touch only the
        // original interaction graph and triple residuals are higher-order
        // terms, so repeating the full fill-cell sweep here is redundant.

        if (!hasTriple) {
            correctedPairwiseGlobalOffsetKcal = eta.globalOffsetKcal;
        }

        return corrected;
    }

    /**
     * Start a pair-only refit on the reusable corrected matrix.  Restoring
     * touched cells from the immutable base matrix avoids round-off drift from
     * subtracting the previous correction and keeps the refit numerically
     * equivalent to a fresh {@code new EnergyMatrix(branchMinimizingEmat)}.
     */
    private EnergyMatrix beginPairwiseCorrectedEmat() {
        if (correctedPairwiseEmat == null
                || !Objects.equals(
                correctedPairwiseEmatSignature,
                proposalInteractionGraphSignature)) {
            correctedPairwiseEmat = new EnergyMatrix(branchMinimizingEmat);
            correctedPairwiseEmatSignature = proposalInteractionGraphSignature;
            correctedPairwiseOneBodyCount = 0;
            correctedPairwisePairCount = 0;
            correctedPairwiseGlobalOffsetKcal = 0.0;
            return correctedPairwiseEmat;
        }

        if (correctedPairwiseGlobalOffsetKcal != 0.0) {
            int numRC = correctedPairwiseEmat.getNumConfAtPos(0);
            for (int rc = 0; rc < numRC; rc++) {
                correctedPairwiseEmat.setOneBody(
                        0, rc, branchMinimizingEmat.getOneBody(0, rc));
            }
        }
        for (int i = 0; i < correctedPairwiseOneBodyCount; i++) {
            long key = correctedPairwiseOneBodyKeys[i];
            int pos = (int) (key >>> 32);
            int rc = (int) key;
            correctedPairwiseEmat.setOneBody(
                    pos, rc, branchMinimizingEmat.getOneBody(pos, rc));
        }
        for (int i = 0; i < correctedPairwisePairCount; i++) {
            long key = correctedPairwisePairKeys[i];
            int pos1 = EtaCorrections.pairPos1(key);
            int rc1 = EtaCorrections.pairRc1(key);
            int pos2 = EtaCorrections.pairPos2(key);
            int rc2 = EtaCorrections.pairRc2(key);
            correctedPairwiseEmat.setPairwise(
                    pos1, rc1, pos2, rc2,
                    branchMinimizingEmat.getPairwise(
                            pos1, rc1, pos2, rc2));
        }
        correctedPairwiseOneBodyCount = 0;
        correctedPairwisePairCount = 0;
        correctedPairwiseGlobalOffsetKcal = 0.0;
        return correctedPairwiseEmat;
    }

    private void ensureCorrectedOneBodyKeyCapacity(int capacity) {
        if (correctedPairwiseOneBodyKeys.length < capacity) {
            correctedPairwiseOneBodyKeys = Arrays.copyOf(
                    correctedPairwiseOneBodyKeys,
                    Math.max(capacity,
                            Math.max(4, correctedPairwiseOneBodyKeys.length * 2)));
        }
    }

    private void ensureCorrectedPairKeyCapacity(int capacity) {
        if (correctedPairwisePairKeys.length < capacity) {
            correctedPairwisePairKeys = Arrays.copyOf(
                    correctedPairwisePairKeys,
                    Math.max(capacity,
                            Math.max(4, correctedPairwisePairKeys.length * 2)));
        }
    }

    private void invalidateProposalRigidEmat() {
        proposalRigidEmat = null;
        proposalRigidEmatSignature = null;
        proposalMinimizingEmat = null;
        proposalMinimizingEmatSignature = null;
        correctedPairwiseEmat = null;
        correctedPairwiseEmatSignature = null;
        correctedPairwiseOneBodyCount = 0;
        correctedPairwisePairCount = 0;
        correctedPairwiseGlobalOffsetKcal = 0.0;
    }

    private EnergyMatrix getProposalRigidEmat() {
        if (proposalRigidEmat == null
                || !Objects.equals(
                proposalRigidEmatSignature,
                proposalInteractionGraphSignature)) {
            proposalRigidEmat = new EnergyMatrix(branchRigidEmat);
            zeroProposalFillPairEnergies(proposalRigidEmat);
            proposalRigidEmatSignature = proposalInteractionGraphSignature;
        }
        return proposalRigidEmat;
    }

    private EnergyMatrix getProposalMinimizingEmat() {
        if (proposalMinimizingEmat == null
                || !Objects.equals(
                proposalMinimizingEmatSignature,
                proposalInteractionGraphSignature)) {
            proposalMinimizingEmat = new EnergyMatrix(branchMinimizingEmat);
            zeroProposalFillPairEnergies(proposalMinimizingEmat);
            proposalMinimizingEmatSignature = proposalInteractionGraphSignature;
        }
        return proposalMinimizingEmat;
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

        invalidateProposalRigidEmat();
        samplingTopDownOrder = null;

        if (fillEdges.isEmpty()) {
            rootedRoot = initialRootedRoot;
            rootedRootEdge = initialRootedRootEdge;
            proposalInteractionGraph = interactionGraph;
            proposalInteractionGraphSignature = signature;
            System.out.println("[PACK*-triple-eta-dp] restored base proposal"
                    + " graph: edges=" + interactionGraph.getNumEdges());
            return;
        }

        PackStarTripleDecompositionCosts.ProposalRoot proposal =
                PackStarTripleDecompositionCosts.rootProposalGraph(desired, rcs, confSpace, true);
        RootedTreeNode rebuiltRoot = proposal.selected.root;
        RootedTreeEdge rebuiltRootEdge =
                rebuiltRoot.getLeftChild().getChildOfEdge();
        rootedRoot = rebuiltRoot;
        rootedRootEdge = rebuiltRootEdge;
        proposalInteractionGraph = desired;
        proposalInteractionGraphSignature = signature;
        System.out.println("[PACK*-triple-eta-dp] rebuilt exact proposal"
                + " graph: baseEdges=" + interactionGraph.getNumEdges()
                + ", fillEdges=" + fillEdges.size()
                + ", proposalEdges=" + desired.getNumEdges()
                + ", branchwidth=" + proposal.branchwidth
                + ", rootSplit=" + proposal.selected.splitEdgeIndex
                + ", rootSelectionMs=" + proposal.rootSelectionNanos / 1e6
                + ", gpuWork=" + proposal.selected.gpuWork
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
    private CorrectedDPResult recomputeDP(EnergyMatrix correctedEmat) {
        // Keep the rigid lower-bound matrix unchanged and recompute the
        // corrected proposal tables directly. The only retained route never
        // reuses corrected DP tables across independently learned proposals.
        EnergyMatrix proposalRigidEmat = getProposalRigidEmat();
        RootedTreeEdge.postOrderInitIncremental(
                rootedRoot, proposalRigidEmat, correctedEmat,
                proposalInteractionGraph, RT);
        RootedTreeEdge.postOrderComputeFullDP(rootedRoot);
        activeProposalEmat = correctedEmat;
        return new CorrectedDPResult(rootedRootEdge.getLogZUpper(0));
    }


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
        for (int i = 0; i < conf.length; i++) {
            if (conf[i] < 0) continue;
            energy += emat.getOneBody(i, conf[i]);
            for (int edgeIndex : interactionEdgeIndicesByFirstPosition[i]) {
                int pos2 = interactionEdges[edgeIndex][1];
                if (conf[pos2] < 0) continue;
                energy += emat.getPairwise(
                        i, conf[i], pos2, conf[pos2]);
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
            inters.addAll(getSingleInteractionTemplate(pos, conf[pos]));
        }
        for (int[] edge : interactionEdges) {
            int pos1 = edge[0];
            int pos2 = edge[1];
            if (conf[pos1] < 0 || conf[pos2] < 0) continue;
            inters.addAll(getPairInteractionTemplate(
                    pos1, conf[pos1], pos2, conf[pos2]));
        }
        if (minimizingEcalc.addShellInters) {
            inters.addAll(getShellInteractionTemplate());
        }
        return inters;
    }

    /**
     * Convert a log-space value to BigDecimal.
     */
    static BigDecimal bigExpFromLog(double logVal) {
        if (logVal == Double.NEGATIVE_INFINITY) return BigDecimal.ZERO;
        if (!Double.isFinite(logVal)) return MathTools.BigNaN;

        // Use the arbitrary-precision exponential path for very negative logs;
        // ordinary double-space exponentiation would underflow there.
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
    /** Number of logical CCD sample records, including deduplicated repeats. */
    public int getTotalCCDSampleRecords() { return totalCCDSampleRecords; }
    public double getMeanPsi() { return meanPsi; }
    public double getCvPsi() { return cvPsi; }
    public double getMeanResidual() { return meanResidual; }
    public double getStdResidual() { return stdResidual; }
    public boolean hasValidCertificate() { return certificateValid; }
    public String getCertificateFailureReason() { return certificateFailureReason; }
}
