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
 * no deterministic MARK* fallback, no branch lookahead/region-atom controls, and
 * the PACK* estimator path instead of the deterministic search loop.</p>
 */
final class PackStarBranchDpBackend extends BranchDpBackend implements PackStarBackend {

    private static final String RESTORE_DP_PROPERTY = "packstar.pac.restoreDP";
    private static final String ETA_ENABLED_PROPERTY = "packstar.pac.etaEnabled";
    private static final String ITERATE_PROPERTY = "packstar.pac.iterate";
    private static final String ITERATE_MAX_ROUNDS_PROPERTY =
            "packstar.pac.iterate.maxRounds";
    private static final int DEFAULT_ITERATE_MAX_ROUNDS = 4;

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

        // Fail-safe upper bound for the current PACK* estimator: one initial
        // p_m DP, one corrected p_eta DP, up to maxRounds corrected refinement
        // DPs, and an optional compatibility restore of the initial DP.
        boolean etaEnabled = getConfigBoolean(ETA_ENABLED_PROPERTY, true);
        boolean iterate = getConfigBoolean(ITERATE_PROPERTY, true);
        int maxRounds = Math.max(0, getConfigInteger(
                ITERATE_MAX_ROUNDS_PROPERTY,
                DEFAULT_ITERATE_MAX_ROUNDS));
        boolean restore = getConfigBoolean(RESTORE_DP_PROPERTY, false);
        return conservativeDpSweeps(etaEnabled, iterate, maxRounds, restore);
    }

    static int conservativeDpSweeps(boolean etaEnabled, boolean iterate,
                                    int maxRounds, boolean restore) {
        int sweeps = 2;
        if (etaEnabled && iterate) {
            sweeps += Math.max(0, maxRounds);
        }
        return restore ? sweeps + 1 : sweeps;
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
            System.out.println("PACK*: singleton/empty branch graph; using standard MARK* fallback for this exact small state.");
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

        if (getConfigBoolean(RESTORE_DP_PROPERTY, false)) {
            // Optional compatibility behavior; normal completion only needs
            // PACK*'s returned q bounds, so this remains disabled by default.
            restoreInitialDPTables("PACK*:");
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

        if (estimatorEpsilon <= targetEpsilon) {
            setStatus(PartitionFunction.Status.Estimated);
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
