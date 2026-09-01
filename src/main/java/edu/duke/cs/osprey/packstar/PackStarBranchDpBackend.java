package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.branchdp.BranchDpBackend;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;

/**
 * PACK* branch-DP backend.
 *
 * <p>The DP and sampling primitives live in the shared branch-DP engine, but
 * PACK* owns the runtime policy here: packstar.* config aliases,
 * no deterministic MARK* fallback, no branch lookahead/region-atom controls,
 * and a fixed-budget certified-or-abort estimator path.</p>
 */
final class PackStarBranchDpBackend extends BranchDpBackend implements PackStarBackend {

    private static final String ETA_ENABLED_PROPERTY = "packstar.pac.etaEnabled";
    private static final String ITERATE_PROPERTY = "packstar.pac.iterate";
    private static final String ITERATE_MAX_ROUNDS_PROPERTY =
            "packstar.pac.iterate.maxRounds";
    private static final String UNCONDITIONAL_TAIL_PROPERTY =
            "packstar.pac.unconditionalTail";
    private static final String DP_RICH_PROPERTY = "packstar.pac.dpRich";
    private static final String ETA_V4_PROPERTY = "packstar.pac.etaV4";
    private static final String ETA_V4_TAIL_POWERS_PROPERTY =
            "packstar.pac.etaV4.tailPowers";
    private static final String ETA_V4_TRUST_KCAL_PROPERTY =
            "packstar.pac.etaV4.trustKcal";
    private static final String CONDITIONAL_REPAIR_AUDIT_PROPERTY =
            "packstar.pac.conditionalRepairAudit";
    private static final String FREQUENCY_SEVERITY_PRODUCTION_PROPERTY =
            "packstar.pac.frequencySeverityProduction";
    private static final String FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY =
            "packstar.pac.frequencySeverity.maxRefits";
    private static final String TILT_LAMBDAS_PROPERTY =
            "packstar.pac.tail.tiltLambdas";
    private static final String ETA_REPAIR_PROPERTY =
            "packstar.pac.etaRepair";
    private static final String ETA_REPAIR_MAX_ROUNDS_PROPERTY =
            "packstar.pac.etaRepair.maxRounds";
    private static final int DEFAULT_ITERATE_MAX_ROUNDS = 4;
    private static final int DEFAULT_ETA_REPAIR_MAX_ROUNDS = 2;
    private static final String DEFAULT_TILT_LAMBDAS = "0,0.25,0.5,1,2,4";
    private static final String DEFAULT_ETA_V4_TAIL_POWERS = "1,2";
    private static final String DEFAULT_ETA_V4_TRUST_KCAL =
            "0.125,0.25,0.5,1";
    // raw eta, current repair, one repair per B={1,1.5,2}, zero eta;
    // candidate alphas can deduplicate, but admission must budget the maximum,
    // followed by one winner reload before monitor/final sampling.
    private static final int MAX_DP_RICH_CANDIDATES = 6;

    private PackStarSampleListener sampleListener = null;
    private final String configuredSeedStateRole;
    private Integer calculationInstanceId = null;

    PackStarBranchDpBackend(SimpleConfSpace confSpace,
                            EnergyMatrix rigidEmat,
                            EnergyMatrix minimizingEmat,
                            ConfEnergyCalculator minimizingConfEcalc,
                            RCs rcs,
                            Parallelism parallelism,
                            String stateNameOverride) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
                rcs, parallelism, stateNameOverride);
        this.configuredSeedStateRole = stateNameOverride == null
                || stateNameOverride.trim().isEmpty()
                ? "standalone" : stateNameOverride.trim();
    }

    @Override
    protected String getConfigProperty(String key, String defaultValue) {
        return PackStarConfig.getProperty(key, defaultValue);
    }

    @Override
    protected EdgeSelectionStrategy getEdgeSelectionStrategy() {
        return EdgeSelectionStrategy.LAMBDA_STATES;
    }

    @Override
    protected int getEdgeLookaheadMaxStates() {
        return 1;
    }

    @Override
    protected int getEdgeLookaheadMaxPendingEdges() {
        return 2;
    }

    @Override
    protected boolean getEdgeLookaheadParallel() {
        return false;
    }

    @Override
    protected boolean getUseHigherOrderCorrections() {
        return false;
    }

    @Override
    protected int getAdmissionDpSweeps() {
        int explicit = getConfigInteger(
                BranchDpAdmission.DP_SWEEPS_PROPERTY, 0);
        if (explicit > 0) {
            return explicit;
        }

        if (getConfigBoolean(
                FREQUENCY_SEVERITY_PRODUCTION_PROPERTY, false)) {
            int maxRefits = Math.max(0, getConfigInteger(
                    FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY, 8));
            // One initial q_m sweep plus one proposal sweep for each round
            // 0..maxRefits.  A round-zero eta=0 proposal can reuse q_m only
            // when no triple eta is active, so budget the non-reuse case.
            return 2 + maxRefits;
        }

        // Fail-safe upper bound for the current PACK* estimator: one initial
        // p_m DP, one corrected p_eta DP, and up to maxRounds corrected
        // refinement DPs.  Archived unconditional mode also needs tilted
        // proposal sweeps. Certificate miss aborts and never restores the
        // original DP for deterministic search.
        boolean etaEnabled = getConfigBoolean(ETA_ENABLED_PROPERTY, true);
        boolean iterate = getConfigBoolean(ITERATE_PROPERTY, true);
        boolean etaRepair = etaEnabled
                && getConfigBoolean(ETA_REPAIR_PROPERTY, true);
        int maxRounds = Math.max(0, getConfigInteger(
                ITERATE_MAX_ROUNDS_PROPERTY,
                DEFAULT_ITERATE_MAX_ROUNDS));
        if (etaRepair) {
            maxRounds = Math.max(maxRounds, Math.max(0, getConfigInteger(
                    ETA_REPAIR_MAX_ROUNDS_PROPERTY,
                    DEFAULT_ETA_REPAIR_MAX_ROUNDS)));
        }
        int additionalDpSweeps = 0;
        boolean conditionalRepairAudit = getConfigBoolean(
                CONDITIONAL_REPAIR_AUDIT_PROPERTY, true);
        if (!conditionalRepairAudit
                && getConfigBoolean(UNCONDITIONAL_TAIL_PROPERTY, false)) {
            additionalDpSweeps += countPositiveTiltLambdas(getConfigProperty(
                    TILT_LAMBDAS_PROPERTY, DEFAULT_TILT_LAMBDAS)) + 1;
        }
        if (conditionalRepairAudit && etaEnabled
                && getConfigBoolean(DP_RICH_PROPERTY, false)) {
            if (getConfigBoolean(ETA_V4_PROPERTY, false)) {
                int tailPowers = countPositiveGridValues(
                        getConfigProperty(
                                ETA_V4_TAIL_POWERS_PROPERTY,
                                DEFAULT_ETA_V4_TAIL_POWERS));
                int trustRadii = countPositiveGridValues(
                        getConfigProperty(
                                ETA_V4_TRUST_KCAL_PROPERTY,
                                DEFAULT_ETA_V4_TRUST_KCAL));
                additionalDpSweeps +=
                        conservativeEtaV4AdditionalDpSweeps(
                                tailPowers, trustRadii);
            } else {
                additionalDpSweeps += MAX_DP_RICH_CANDIDATES + 1;
            }
        }
        return conservativeDpSweeps(
                etaEnabled, iterate || etaRepair,
                maxRounds, additionalDpSweeps, false);
    }

    static int conservativeDpSweeps(boolean etaEnabled, boolean iterate,
                                    int maxRounds, boolean restore) {
        return conservativeDpSweeps(
                etaEnabled, iterate, maxRounds, 0, restore);
    }

    static int conservativeDpSweeps(boolean etaEnabled, boolean iterate,
                                    int maxRounds, int additionalDpSweeps,
                                    boolean restore) {
        int sweeps = 2;
        if (etaEnabled && iterate) {
            sweeps += Math.max(0, maxRounds);
        }
        sweeps += Math.max(0, additionalDpSweeps);
        return restore ? sweeps + 1 : sweeps;
    }

    /**
     * One anchor reload, all tail-power/trust candidates, one tail winner
     * reload, and up to two component reloads for each of mixture calibration,
     * monitor, and final sampling.
     */
    static int conservativeEtaV4AdditionalDpSweeps(
            int tailPowers, int trustRadii) {
        return Math.max(0, tailPowers) * Math.max(0, trustRadii) + 8;
    }

    private static int countPositiveTiltLambdas(String configured) {
        java.util.Set<Double> unique = new java.util.HashSet<>();
        if (configured != null) {
            for (String token : configured.split(",")) {
                String trimmed = token.trim();
                if (trimmed.isEmpty()) continue;
                try {
                    double value = Double.parseDouble(trimmed);
                    if (Double.isFinite(value) && value > 0.0) unique.add(value);
                } catch (NumberFormatException ignored) {
                    // The estimator reports the configuration error. Admission
                    // remains conservative for all parsable candidates.
                }
            }
        }
        return unique.size();
    }

    private static int countPositiveGridValues(String configured) {
        java.util.Set<Double> unique = new java.util.HashSet<>();
        if (configured != null) {
            for (String token : configured.split(",")) {
                String trimmed = token.trim();
                if (trimmed.isEmpty()) continue;
                try {
                    double value = Double.parseDouble(trimmed);
                    if (Double.isFinite(value) && value > 0.0) {
                        unique.add(value);
                    }
                } catch (NumberFormatException ignored) {
                    // The estimator owns detailed configuration validation.
                }
            }
        }
        return unique.size();
    }

    @Override
    protected void logBackendControlOverrides() {
        System.out.println("PACK*: edge-lookahead and triple-correction controls are disabled for the PACK* backend.");
    }

    @Override
    protected void computeWithoutBranchDecomposition(int maxNumConfs) {
        if (rootedRoot == null
                && interactionGraph != null
                && interactionGraph.getNumPositions() <= 1) {
            System.out.println("PACK*: singleton/empty branch graph; using the exact small-state base case (no estimator miss or search fallback).");
            super.computeWithoutBranchDecomposition(maxNumConfs);
            return;
        }
        abortPackStar("branch decomposition is unavailable; refusing to fall back to deterministic MARK*");
    }

    @Override
    protected boolean beforeBranchDecompositionSearch(int maxNumConfs) {
        if (!dpTablesReady) {
            abortPackStar("DP tables are unavailable; PACK* requires the DP/sampling backend");
            return true;
        }

        if (epsilonBound <= targetEpsilon) {
            finishWithCurrentDPBounds("PACK*: sampling skipped; first-round DP bounds already meet target, "
                    + "so CCD sampling and corrected DP are unnecessary. "
                    + "epsilon=" + String.format("%.6f", epsilonBound)
                    + " <= target=" + String.format("%.6f", targetEpsilon));
            return true;
        }

        computeWithPackStarEstimator(maxNumConfs);
        return true;
    }

    @Override
    public void setReduceMinimizations(boolean enabled) {
        reduceMinimizations = enabled;
    }

    @Override
    public void setSampleListener(PackStarSampleListener listener) {
        this.sampleListener = listener;
    }

    @Override
    public void setInstanceId(int val) {
        this.calculationInstanceId = val;
    }

    @Override
    public void close() {
        releaseLargeMemory();
    }

    private void abortPackStar(String reason) {
        System.out.println("PACK*: aborted: " + reason + ".");
        values = new MARKStarBound.Values();
        values.qprime = MathTools.BigPositiveInfinity;
        setStatus(PartitionFunction.Status.Aborted);
    }

    private void computeWithPackStarEstimator(int sampleBudget) {
        System.out.println("PACK*: estimator activated for state=" + stateName
                + ", sampleBudget=" + (sampleBudget == Integer.MAX_VALUE ? "unbounded" : sampleBudget));

        PackStarEstimator estimator = new PackStarEstimator(
                rootedRoot, rootedRootEdge,
                branchMinimizingEmat, branchRigidEmat,
                interactionGraph, getMinimizingEcalc(),
                searchRCs, confSpace,
                targetEpsilon,
                sampleBudget,
                randomStreamIdentity());
        estimator.setSampleListener(sampleListener);

        double estimatorEpsilon = estimator.compute();

        if (!estimator.hasValidCertificate()
                || !PackStarEstimator.isValidCertificate(
                estimator.getZLower(), estimator.getZUpper(), estimatorEpsilon)) {
            totalMinimizations = estimator.getTotalCCDCalls();
            abortPackStar("estimator certificate is invalid: "
                    + estimator.getCertificateFailureReason());
            return;
        }

        PartitionFunction.Values vals = getValues();
        vals.qstar = estimator.getZLower();
        vals.pstar = estimator.getZUpper();
        vals.qprime = vals.pstar.subtract(vals.qstar);

        flatSumZLower = estimator.getZLower();
        flatSumZUpper = estimator.getZUpper();
        epsilonBound = estimatorEpsilon;
        totalMinimizations = estimator.getTotalCCDCalls();

        System.out.println("PACK*: estimator finished. epsilon=" + String.format("%.6f", estimatorEpsilon)
                + ", CCD calls=" + estimator.getTotalCCDCalls()
                + ", cvPsi=" + String.format("%.4f", estimator.getCvPsi())
                + ", meanResidual=" + String.format("%.4f", estimator.getMeanResidual()) + " kcal/mol"
                + ", stdResidual=" + String.format("%.4f", estimator.getStdResidual()) + " kcal/mol");

        if (Double.isFinite(estimatorEpsilon)
                && estimatorEpsilon >= 0.0
                && estimatorEpsilon <= targetEpsilon) {
            setStatus(PartitionFunction.Status.Estimated);
        } else {
            System.out.println("PACK*: aborted: valid PAC interval missed target; epsilon="
                    + String.format("%.6f", estimatorEpsilon)
                    + " > target=" + String.format("%.6f", targetEpsilon)
                    + ", CCD calls=" + estimator.getTotalCCDCalls()
                    + ". No deterministic fallback or branch search was started.");
            setStatus(PartitionFunction.Status.Aborted);
        }
    }

    private String randomStreamIdentity() {
        String stateRole = calculationInstanceId == null
                ? configuredSeedStateRole
                : "kstar-state-" + calculationInstanceId;
        return BranchDpAdmission.stateKey(
                "packstar|" + stateRole + "|pac-v1", searchRCs);
    }
}
