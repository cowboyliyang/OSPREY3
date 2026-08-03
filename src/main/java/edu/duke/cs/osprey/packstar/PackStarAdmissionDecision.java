package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.function.Consumer;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Allocation-free decision core for whole-case PACK* admission.
 *
 * <p>This class deliberately knows nothing about K*, PDB files, CUDA, mmap, or
 * DP-table construction.  Each state supplies an allocation-free preview
 * callback.  The core sums the initial predictions, spends a bounded exact
 * retry budget on the largest contributors only when necessary, and either
 * returns an admitted case or throws before the caller can start the formal
 * materialization path.</p>
 */
final class PackStarAdmissionDecision {

    private static final String LOG_PREFIX = "PACK* preflight:";

    @FunctionalInterface
    interface Previewer {
        BranchDpAdmission.Prediction preview();
    }

    static final class State {

        final String sequence;
        final Previewer previewer;
        private BranchDpAdmission.Prediction prediction;

        State(String sequence, Previewer previewer) {
            if (previewer == null) {
                throw new IllegalArgumentException("previewer is required");
            }
            this.sequence = sequence == null ? "unknown" : sequence;
            this.previewer = previewer;
        }

        BranchDpAdmission.Prediction prediction() {
            return prediction;
        }

        private BranchDpAdmission.Prediction preview() {
            BranchDpAdmission.Prediction result = previewer.preview();
            if (result == null) {
                throw new IllegalStateException(LOG_PREFIX
                        + " state preview returned no admission prediction for sequence="
                        + sequence);
            }
            return result;
        }
    }

    static final class Settings {

        final double slaHours;
        final int finalMaxStates;
        final long finalExactMaxMillis;
        final int maxRounds;
        final int previewThreads;
        final boolean lockLoadedPolicies;
        /** Run exact deepening even when this shard is below its local SLA. */
        final boolean forceFinalPass;
        final BranchDecomposition.ExactImproveOptions initialExactOptions;

        Settings(double slaHours, int finalMaxStates,
                 long finalExactMaxMillis,
                 BranchDecomposition.ExactImproveOptions initialExactOptions) {
            this(slaHours, finalMaxStates, finalExactMaxMillis, 1, 1,
                    false, false, initialExactOptions);
        }

        Settings(double slaHours, int finalMaxStates,
                 long finalExactMaxMillis, int maxRounds, int previewThreads,
                 BranchDecomposition.ExactImproveOptions initialExactOptions) {
            this(slaHours, finalMaxStates, finalExactMaxMillis, maxRounds,
                    previewThreads, false, false, initialExactOptions);
        }

        Settings(double slaHours, int finalMaxStates,
                 long finalExactMaxMillis, int maxRounds, int previewThreads,
                 boolean lockLoadedPolicies,
                 BranchDecomposition.ExactImproveOptions initialExactOptions) {
            this(slaHours, finalMaxStates, finalExactMaxMillis, maxRounds,
                    previewThreads, lockLoadedPolicies, false,
                    initialExactOptions);
        }

        Settings(double slaHours, int finalMaxStates,
                 long finalExactMaxMillis, int maxRounds, int previewThreads,
                 boolean lockLoadedPolicies, boolean forceFinalPass,
                 BranchDecomposition.ExactImproveOptions initialExactOptions) {
            if (!(slaHours > 0.0) || !Double.isFinite(slaHours)) {
                throw new IllegalArgumentException(
                        "whole-case SLA must be finite and positive");
            }
            if (initialExactOptions == null) {
                throw new IllegalArgumentException(
                        "initial exact-improvement options are required");
            }
            this.slaHours = slaHours;
            this.finalMaxStates = Math.max(0, finalMaxStates);
            this.finalExactMaxMillis = Math.max(1L, finalExactMaxMillis);
            this.maxRounds = Math.max(1, maxRounds);
            this.previewThreads = Math.max(1, previewThreads);
            this.lockLoadedPolicies = lockLoadedPolicies;
            this.forceFinalPass = forceFinalPass;
            this.initialExactOptions = initialExactOptions;
        }
    }

    static final class Outcome {

        final BranchDpAdmission.CaseSummary initialSummary;
        final BranchDpAdmission.CaseSummary finalSummary;
        final int attemptedFinalStates;
        final int acceptedFinalStates;

        Outcome(BranchDpAdmission.CaseSummary initialSummary,
                BranchDpAdmission.CaseSummary finalSummary,
                int attemptedFinalStates, int acceptedFinalStates) {
            this.initialSummary = initialSummary;
            this.finalSummary = finalSummary;
            this.attemptedFinalStates = attemptedFinalStates;
            this.acceptedFinalStates = acceptedFinalStates;
        }

        boolean usedFinalPass() {
            return attemptedFinalStates > 0;
        }
    }

    private PackStarAdmissionDecision() {
    }

    static Outcome decide(List<State> states, Settings settings) {
        if (states == null || states.isEmpty()) {
            throw new IllegalArgumentException(
                    "whole-case admission needs at least one state");
        }
        if (settings == null) {
            throw new IllegalArgumentException("admission settings are required");
        }

        ExecutorService pool = settings.previewThreads > 1
                ? Executors.newFixedThreadPool(settings.previewThreads)
                : null;
        try {
            // ---- initial pass: preview every unique state (optionally parallel) ----
            runPreviews(states, state -> {
                state.prediction = state.preview();
                if (!settings.lockLoadedPolicies) {
                    retainAcceptedInitialPlan(state.prediction,
                            settings.initialExactOptions);
                }
            }, pool);
            for (int i = 0; i < states.size(); i++) {
                State state = states.get(i);
                System.out.println(LOG_PREFIX + " phase=initial state="
                        + (i + 1) + "/" + states.size()
                        + " sequence=" + state.sequence
                        + " " + state.prediction.format());
            }

            BranchDpAdmission.CaseSummary initial = summarize(
                    states, settings.slaHours);
            logSummary("initial", initial, states);
            if (initial.withinSla() && !settings.forceFinalPass) {
                System.out.println(LOG_PREFIX
                        + " admitted after initial pass; retainedExactPolicies="
                        + BranchDpAdmission.exactPolicyCount() + ".");
                return new Outcome(initial, initial, 0, 0);
            }

            if (settings.lockLoadedPolicies) {
                // Loaded plan is authoritative; do not re-search in production.
                reject(initial, states,
                        "loaded exact-policy plan still exceeds SLA"
                                + " (lock mode: no re-search)");
            }

            if (settings.finalMaxStates <= 0) {
                reject(initial, states,
                        "whole case exceeds SLA and finalMaxStates=0");
            }

            // ---- bounded multi-round exact deepening on largest contributors ----
            Set<String> stalled = ConcurrentHashMap.newKeySet();
            int attempted = 0;
            int accepted = 0;
            int round = 0;
            while (round < settings.maxRounds) {
                round++;
                if (!settings.forceFinalPass
                        && summarize(states, settings.slaHours).withinSla()) {
                    break;
                }
                final int currentRound = round;
                List<State> sorted = new ArrayList<>(states);
                sorted.sort((a, b) -> Double.compare(
                        b.prediction.predictedSeconds,
                        a.prediction.predictedSeconds));
                List<State> targets = new ArrayList<>();
                for (State state : sorted) {
                    if (targets.size() >= settings.finalMaxStates) break;
                    if (!(state.prediction.predictedSeconds > 0.0)) continue;
                    if (stalled.contains(state.prediction.stateKey)) continue;
                    targets.add(state);
                }
                if (targets.isEmpty()) {
                    break;
                }

                int roundAttempted;
                int roundAccepted;
                if (pool == null) {
                    // sequential: early-stop as soon as the case fits the SLA
                    // (behavior-identical to the original single final pass).
                    int a = 0;
                    int acc = 0;
                    for (State state : targets) {
                        if (!settings.forceFinalPass
                                && summarize(states, settings.slaHours).withinSla()) {
                            break;
                        }
                        a++;
                        if (attemptDeepen(state, settings, currentRound)) {
                            acc++;
                        } else {
                            stalled.add(state.prediction.stateKey);
                        }
                    }
                    roundAttempted = a;
                    roundAccepted = acc;
                } else {
                    // parallel: attempt every selected target this round.
                    AtomicInteger acc = new AtomicInteger(0);
                    runPreviews(targets, state -> {
                        if (attemptDeepen(state, settings, currentRound)) {
                            acc.incrementAndGet();
                        } else {
                            stalled.add(state.prediction.stateKey);
                        }
                    }, pool);
                    roundAttempted = targets.size();
                    roundAccepted = acc.get();
                }
                attempted += roundAttempted;
                accepted += roundAccepted;
                System.out.println(LOG_PREFIX + " phase=final round=" + round
                        + "/" + settings.maxRounds + " attempted=" + roundAttempted
                        + " roundAccepted=" + roundAccepted
                        + " predictedCaseHours=" + String.format(Locale.ROOT,
                        "%.4f",
                        summarize(states, settings.slaHours).totalHours()));
                if (roundAccepted == 0) {
                    break;
                }
            }

            BranchDpAdmission.CaseSummary result = summarize(
                    states, settings.slaHours);
            logSummary("final", result, states);
            if (!result.withinSla()) {
                reject(result, states,
                        "whole case still exceeds SLA after bounded multi-round"
                                + " pass (rounds=" + round + ", attempted="
                                + attempted + ", accepted=" + accepted + ")");
            }

            System.out.println(LOG_PREFIX
                    + " admitted after multi-round pass; rounds=" + round
                    + ", attempted=" + attempted + ", accepted=" + accepted
                    + ", retainedExactPolicies="
                    + BranchDpAdmission.exactPolicyCount()
                    + ", allocationPhase=before-materialization.");
            return new Outcome(initial, result, attempted, accepted);
        } finally {
            if (pool != null) {
                pool.shutdown();
            }
        }
    }

    /** Run an allocation-free preview action over states, optionally in parallel. */
    private static void runPreviews(List<State> states,
                                    Consumer<State> action,
                                    ExecutorService pool) {
        if (pool == null) {
            for (State state : states) {
                action.accept(state);
            }
            return;
        }
        List<Future<?>> futures = new ArrayList<>();
        for (State state : states) {
            futures.add(pool.submit(() -> action.accept(state)));
        }
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
                throw new RuntimeException("preflight preview interrupted", ex);
            } catch (ExecutionException ex) {
                Throwable cause = ex.getCause();
                if (cause instanceof RuntimeException) {
                    throw (RuntimeException) cause;
                }
                throw new RuntimeException("preflight preview failed", cause);
            }
        }
    }

    /**
     * One bounded exact-improvement attempt for a single state, escalating the
     * branchwidth-drop target with the round index.  Returns true iff a strictly
     * faster plan was accepted (its exact policy is then retained for the real
     * run via the prediction ceiling).  Distinct states have distinct state
     * keys and distinct State objects, so concurrent attempts do not share
     * mutable fields; the exact-policy store is a ConcurrentHashMap.
     */
    private static boolean attemptDeepen(State state, Settings settings,
                                         int round) {
        BranchDpAdmission.Prediction previous = state.prediction;
        BranchDpAdmission.ExactPolicy previousPolicy =
                BranchDpAdmission.getExactPolicy(previous.stateKey);
        int base = previous.adaptiveAccepted
                ? previous.acceptedBranchwidthDrop() : 0;
        int nextDrop = Math.max(base + 1, round);
        BranchDpAdmission.ExactPolicy trialPolicy =
                new BranchDpAdmission.ExactPolicy(
                        nextDrop, nextDrop, settings.finalExactMaxMillis);
        BranchDpAdmission.putExactPolicy(previous.stateKey, trialPolicy);

        System.out.println(LOG_PREFIX + " phase=final round=" + round
                + " state=" + previous.stateName
                + " sequence=" + state.sequence
                + " targetDrop=" + nextDrop
                + " maxMillis=" + settings.finalExactMaxMillis
                + " previousHours=" + String.format(Locale.ROOT,
                "%.4f", previous.predictedHours()));

        BranchDpAdmission.Prediction candidate;
        try {
            candidate = state.preview();
        } catch (RuntimeException ex) {
            BranchDpAdmission.putExactPolicy(previous.stateKey, previousPolicy);
            System.err.println(LOG_PREFIX
                    + " deepen candidate failed before materialization for state="
                    + previous.stateName + ": " + ex.getMessage()
                    + "; keeping previous plan.");
            return false;
        }

        if (!sameStateIdentity(candidate, previous)) {
            BranchDpAdmission.putExactPolicy(previous.stateKey, previousPolicy);
            System.err.println(LOG_PREFIX
                    + " deepen candidate changed state identity for state="
                    + previous.stateName + "; keeping previous plan.");
            return false;
        }

        if (isStrictImprovement(candidate, previous)) {
            state.prediction = candidate;
            BranchDpAdmission.putExactPolicy(previous.stateKey,
                    trialPolicy.withPredictionCeiling(
                            candidate.predictedSeconds));
            System.out.println(LOG_PREFIX
                    + " phase=final accepted state=" + candidate.stateName
                    + " branchwidth=" + previous.branchwidth
                    + "->" + candidate.branchwidth
                    + " predictedHours=" + String.format(Locale.ROOT,
                    "%.4f->%.4f", previous.predictedHours(),
                    candidate.predictedHours()));
            return true;
        }
        BranchDpAdmission.putExactPolicy(previous.stateKey, previousPolicy);
        System.out.println(LOG_PREFIX
                + " phase=final rejected state=" + previous.stateName
                + " candidateHours=" + String.format(Locale.ROOT,
                "%.4f", candidate.predictedHours())
                + " previousHours=" + String.format(Locale.ROOT,
                "%.4f", previous.predictedHours()));
        return false;
    }

    private static boolean sameStateIdentity(
            BranchDpAdmission.Prediction candidate,
            BranchDpAdmission.Prediction previous) {
        return candidate.stateKey != null
                && candidate.stateKey.equals(previous.stateKey);
    }

    private static boolean isStrictImprovement(
            BranchDpAdmission.Prediction candidate,
            BranchDpAdmission.Prediction previous) {
        if (!candidate.hasFinitePrediction()) return false;
        if (!previous.hasFinitePrediction()) return true;
        double tolerance = Math.max(1.0e-9,
                previous.predictedSeconds * 1.0e-12);
        return candidate.predictedSeconds
                < previous.predictedSeconds - tolerance;
    }

    private static void retainAcceptedInitialPlan(
            BranchDpAdmission.Prediction prediction,
            BranchDecomposition.ExactImproveOptions options) {
        if (!prediction.adaptiveAccepted) return;
        BranchDpAdmission.putExactPolicy(prediction.stateKey,
                new BranchDpAdmission.ExactPolicy(
                        options.minDrop, options.maxDrop, options.maxMillis,
                        prediction.predictedSeconds));
    }

    private static BranchDpAdmission.CaseSummary summarize(
            List<State> states, double slaHours) {
        List<BranchDpAdmission.Prediction> predictions = new ArrayList<>();
        for (State state : states) {
            predictions.add(state.prediction);
        }
        return new BranchDpAdmission.CaseSummary(predictions, slaHours);
    }

    private static void logSummary(String phase,
                                   BranchDpAdmission.CaseSummary summary,
                                   List<State> states) {
        int finite = 0;
        for (State state : states) {
            if (state.prediction.hasFinitePrediction()) finite++;
        }
        System.out.println(LOG_PREFIX + " phase=" + phase
                + " uniqueStates=" + states.size()
                + " finiteStates=" + finite
                + " predictedCaseHours=" + String.format(Locale.ROOT,
                "%.4f", summary.totalHours())
                + " caseSlaHours=" + String.format(Locale.ROOT,
                "%.4f", summary.slaHours)
                + " admitted=" + summary.withinSla()
                + ", allocationPhase=before-materialization");
    }

    private static void reject(BranchDpAdmission.CaseSummary summary,
                               List<State> states, String reason) {
        List<State> contributors = new ArrayList<>(states);
        contributors.sort((a, b) -> Double.compare(
                b.prediction.predictedSeconds,
                a.prediction.predictedSeconds));
        StringBuilder top = new StringBuilder();
        int limit = Math.min(5, contributors.size());
        for (int i = 0; i < limit; i++) {
            State state = contributors.get(i);
            if (i > 0) top.append("; ");
            top.append(state.prediction.stateName)
                    .append('/')
                    .append(state.sequence)
                    .append('=')
                    .append(String.format(Locale.ROOT, "%.4fh",
                            state.prediction.predictedHours()));
        }
        throw new IllegalStateException(LOG_PREFIX + " rejected before DP-table"
                + " materialization: " + reason
                + "; predictedCaseHours="
                + String.format(Locale.ROOT, "%.4f", summary.totalHours())
                + "; caseSlaHours="
                + String.format(Locale.ROOT, "%.4f", summary.slaHours)
                + "; topContributors=[" + top + "].");
    }
}
