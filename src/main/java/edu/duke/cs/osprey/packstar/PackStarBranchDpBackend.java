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

    private static final String FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY =
            "packstar.pac.frequencySeverity.maxRefits";
    private PackStarSampleListener sampleListener = null;
    private String functionalEventName = null;
    private PackStarFunctionalEvent functionalEvent = null;
    private PackStarFunctionalObservableResult functionalObservableResult =
            PackStarFunctionalObservableResult.notConfigured();
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
        int maxRefits = Math.max(0, getConfigInteger(
                FREQUENCY_SEVERITY_MAX_REFITS_PROPERTY, 8));
        return frequencySeverityDpSweeps(maxRefits);
    }

    static int frequencySeverityDpSweeps(int maxRefits) {
        // One initial q_m sweep plus one proposal sweep for each round
        // 0..maxRefits. Round zero may reuse q_m only in a special no-triple
        // case, so admission reserves the non-reuse upper bound.
        return 2 + Math.max(0, maxRefits);
    }

    @Override
    protected void computeWithoutBranchDecomposition(int maxNumConfs) {
        if (rootedRoot == null
                && interactionGraph != null
                && interactionGraph.getNumPositions() <= 1) {
            markFunctionalSampleSkipped("exact singleton/empty graph path has no final event sample");
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
            markFunctionalSampleSkipped("initial DP bounds met the Z target; no final event sample");
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

    @Override
    public PackStarFunctionalObservableResult getFunctionalObservableResult() {
        return functionalObservableResult;
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
        markFunctionalSampleSkipped("calculation aborted before final event evaluation: " + reason);
        System.out.println("PACK*: aborted: " + reason + ".");
        values = new MARKStarBound.Values();
        values.qprime = MathTools.BigPositiveInfinity;
        setStatus(PartitionFunction.Status.Aborted);
    }

    private void markFunctionalSampleSkipped(String reason) {
        if (functionalEvent != null && functionalObservableResult.getStatus()
                == PackStarFunctionalObservableResult.Status.NOT_COMPUTED) {
            functionalObservableResult = PackStarFunctionalObservableResult.notComputed(
                    functionalEventName, reason);
        }
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
        estimator.setFunctionalEvent(functionalEventName, functionalEvent);

        double estimatorEpsilon = estimator.compute();
        functionalObservableResult = estimator.getFunctionalObservableResult();

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
