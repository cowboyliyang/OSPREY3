package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDpAdmission;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Supplier;

/** Hardware-independent, bounded structural search. This is not SLA admission. */
final class PackStarAutoPreflight {

    static final class Settings {
        final int patience, maxRounds, threads;
        final double minGain;
        final long initialMillis, maxMillis;

        Settings(int patience, int maxRounds, int threads, double minGain,
                 long initialMillis, long maxMillis) {
            if (patience < 1 || maxRounds < 1 || threads < 1
                    || !Double.isFinite(minGain) || minGain <= 0 || minGain >= 1
                    || initialMillis < 1 || maxMillis < initialMillis) {
                throw new IllegalArgumentException("invalid automatic preflight settings");
            }
            this.patience = patience;
            this.maxRounds = maxRounds;
            this.threads = threads;
            this.minGain = minGain;
            this.initialMillis = initialMillis;
            this.maxMillis = maxMillis;
        }

        long budget(int round) {
            long budget = initialMillis;
            for (int i = 1; i < round && budget < maxMillis; i++) {
                budget = budget > maxMillis / 4 ? maxMillis : budget * 4;
            }
            return budget;
        }
    }

    static final class State {
        final Supplier<BranchDpAdmission.Prediction> preview;
        BranchDpAdmission.Prediction best;
        int quietRounds;
        int failures;

        State(Supplier<BranchDpAdmission.Prediction> preview) {
            this.preview = preview;
        }

        BranchDpAdmission.Prediction read() {
            BranchDpAdmission.Prediction p = preview.get();
            if (p == null || p.stateKey == null || p.gpuWork.signum() < 0) {
                throw new IllegalStateException("automatic preflight has no valid structural preview");
            }
            return p;
        }
    }

    static final class Result {
        final BranchDpAdmission.CaseSummary summary;
        final BigInteger initialWork, finalWork;
        final String stopReason;
        final int rounds;

        Result(List<State> states, BigInteger initialWork, String reason, int rounds) {
            List<BranchDpAdmission.Prediction> predictions = new ArrayList<>();
            for (State state : states) predictions.add(state.best);
            // No fabricated throughput or hours: structural search is not admission.
            summary = BranchDpAdmission.CaseSummary.structural(predictions);
            this.initialWork = initialWork;
            finalWork = totalWork(states);
            stopReason = reason;
            this.rounds = rounds;
        }
    }

    static Result optimize(List<State> states, Settings settings,
                           boolean locked, Runnable checkpoint) {
        if (states.isEmpty()) throw new IllegalArgumentException("no preflight states");
        ExecutorService pool = Executors.newFixedThreadPool(settings.threads);
        try {
            parallel(states, pool, state -> state.best = state.read());
            BigInteger initial = totalWork(states);
            checkpoint.run();
            if (locked) return new Result(states, initial, "loaded-policy", 0);
            int round = 0;
            String reason = "search-limit";
            while (round < settings.maxRounds) {
                List<State> active = new ArrayList<>();
                for (State state : states) {
                    if (state.quietRounds < settings.patience
                            && state.best.totalGpuWork().signum() > 0) active.add(state);
                }
                if (active.isEmpty()) {
                    reason = "plateau";
                    break;
                }
                active.sort((a, b) -> b.best.totalGpuWork().compareTo(a.best.totalGpuWork()));
                final long millis = settings.budget(++round);
                parallel(active, pool, state -> deepen(state, settings, millis));
                checkpoint.run();
                System.out.println("PACK* auto preflight: round=" + round
                        + " activeStates=" + active.size() + " totalDpWork=" + totalWork(states)
                        + " exactMillisPerState=" + millis);
            }
            // Distinguish a plateau on the final allowed round from a budget stop.
            boolean quiet = states.stream().allMatch(s -> s.quietRounds >= settings.patience
                    || s.best.totalGpuWork().signum() == 0);
            if (quiet) reason = "plateau";
            if (states.stream().anyMatch(s -> s.failures > 0)) reason += "-with-preview-errors";
            Result result = new Result(states, initial, reason, round);
            System.out.println("PACK* auto preflight: stopReason=" + reason
                    + " rounds=" + round + " initialDpWork=" + initial
                    + " finalDpWork=" + result.finalWork
                    + " reduction=" + relativeGain(initial, result.finalWork)
                    + " (structural optimization only; no runtime admission)");
            return result;
        } finally {
            pool.shutdownNow();
        }
    }

    private static void deepen(State state, Settings settings, long millis) {
        BranchDpAdmission.Prediction previous = state.best;
        BranchDpAdmission.ExactPolicy saved = BranchDpAdmission.getExactPolicy(previous.stateKey);
        int drop = previous.acceptedBranchwidthDrop() + 1;
        BranchDpAdmission.ExactPolicy trial = new BranchDpAdmission.ExactPolicy(drop, drop, millis);
        BranchDpAdmission.putExactPolicy(previous.stateKey, trial);
        boolean keep = false;
        try {
            BranchDpAdmission.Prediction candidate = state.read();
            if (!previous.stateKey.equals(candidate.stateKey)) {
                throw new IllegalStateException("preview changed state identity");
            }
            // Never exchange lower work for loss of GPU structural support.
            boolean improvement = candidate.gpuUnsupportedEdges <= previous.gpuUnsupportedEdges
                    && candidate.totalGpuWork().compareTo(previous.totalGpuWork()) < 0;
            if (improvement) {
                double gain = relativeGain(previous.totalGpuWork(), candidate.totalGpuWork());
                state.best = candidate;
                BranchDpAdmission.putExactPolicy(previous.stateKey,
                        trial.withWorkCeiling(candidate.totalGpuWork()));
                keep = true;
                state.quietRounds = gain >= settings.minGain ? 0 : state.quietRounds + 1;
                System.out.println("PACK* auto preflight: accepted state=" + candidate.stateKey
                        + " bw=" + previous.branchwidth + "->" + candidate.branchwidth
                        + " root=" + candidate.rootSplitEdge
                        + " dpWork=" + previous.totalGpuWork() + "->" + candidate.totalGpuWork());
            } else {
                state.quietRounds++;
            }
        } catch (RuntimeException ex) {
            state.failures++;
            state.quietRounds++;
            System.err.println("PACK* auto preflight: keeping previous plan for "
                    + previous.stateKey + ": " + ex.getMessage());
        } finally {
            if (!keep) BranchDpAdmission.putExactPolicy(previous.stateKey, saved);
        }
    }

    static double relativeGain(BigInteger before, BigInteger after) {
        if (before.signum() == 0) return 0.0;
        return new BigDecimal(before.subtract(after))
                .divide(new BigDecimal(before), MathContext.DECIMAL64).doubleValue();
    }

    private static BigInteger totalWork(List<State> states) {
        BigInteger total = BigInteger.ZERO;
        for (State state : states) total = total.add(state.best.totalGpuWork());
        return total;
    }

    private static void parallel(List<State> states, ExecutorService pool,
                                 java.util.function.Consumer<State> action) {
        List<Future<?>> futures = new ArrayList<>();
        for (State state : states) futures.add(pool.submit(() -> action.accept(state)));
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException("automatic preflight interrupted", ex);
            } catch (ExecutionException ex) {
                throw new IllegalStateException("automatic preflight preview failed", ex.getCause());
            }
        }
    }
}
