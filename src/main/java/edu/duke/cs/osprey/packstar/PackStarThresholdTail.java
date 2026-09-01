package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.math.BigInteger;
import java.util.Arrays;

/**
 * Exact small-domain diagnostics for the eta-v5 threshold-stratified
 * estimator.
 *
 * <p>This class deliberately does no CCD work and stores no conformation
 * table. It streams the complete RC Cartesian product, accumulates the
 * minimizing-matrix mass in canonical-eta threshold bins, and returns the
 * exact combinatorial oracle evaluated in double precision. Canonical eta
 * removes the arbitrary additive constant in a proposal energy matrix by
 * matching its exact normalizer to the minimizing-matrix normalizer. A
 * scalable certified threshold DP is a later gate; this enumeration is its
 * ground-truth diagnostic on protein/10.</p>
 */
public final class PackStarThresholdTail {

    private PackStarThresholdTail() {}

    public interface ProgressListener {
        void onProgress(long visitedStates, long totalStates);
    }

    public static final class OracleResult {
        public final long cartesianStates;
        public final long visitedStates;
        public final long finiteMinimizingStates;
        public final long zeroMassStates;
        public final long infiniteEtaStates;
        public final double logZMinEnumerated;
        public final double logZProposalEnumerated;
        public final double canonicalShiftKcal;
        public final double logZProposalCanonicalEnumerated;
        public final double rawEtaMinKcal;
        public final double rawEtaMaxKcal;
        public final double canonicalEtaMinKcal;
        public final double canonicalEtaMaxKcal;
        public final double maxMinimizingRebaseDriftKcal;
        public final double maxProposalRebaseDriftKcal;
        public final double[] thresholdsKcal;
        public final long[] tailStateCounts;
        public final double[] logTailMass;

        private OracleResult(
                long cartesianStates,
                long visitedStates,
                long finiteMinimizingStates,
                long zeroMassStates,
                long infiniteEtaStates,
                double logZMinEnumerated,
                double logZProposalEnumerated,
                double canonicalShiftKcal,
                double logZProposalCanonicalEnumerated,
                double rawEtaMinKcal,
                double rawEtaMaxKcal,
                double canonicalEtaMinKcal,
                double canonicalEtaMaxKcal,
                double maxMinimizingRebaseDriftKcal,
                double maxProposalRebaseDriftKcal,
                double[] thresholdsKcal,
                long[] tailStateCounts,
                double[] logTailMass) {
            this.cartesianStates = cartesianStates;
            this.visitedStates = visitedStates;
            this.finiteMinimizingStates = finiteMinimizingStates;
            this.zeroMassStates = zeroMassStates;
            this.infiniteEtaStates = infiniteEtaStates;
            this.logZMinEnumerated = logZMinEnumerated;
            this.logZProposalEnumerated = logZProposalEnumerated;
            this.canonicalShiftKcal = canonicalShiftKcal;
            this.logZProposalCanonicalEnumerated =
                    logZProposalCanonicalEnumerated;
            this.rawEtaMinKcal = rawEtaMinKcal;
            this.rawEtaMaxKcal = rawEtaMaxKcal;
            this.canonicalEtaMinKcal = canonicalEtaMinKcal;
            this.canonicalEtaMaxKcal = canonicalEtaMaxKcal;
            this.maxMinimizingRebaseDriftKcal =
                    maxMinimizingRebaseDriftKcal;
            this.maxProposalRebaseDriftKcal =
                    maxProposalRebaseDriftKcal;
            this.thresholdsKcal = thresholdsKcal.clone();
            this.tailStateCounts = tailStateCounts.clone();
            this.logTailMass = logTailMass.clone();
        }

        public double logTailFraction(int thresholdIndex) {
            double logMass = logTailMass[thresholdIndex];
            if (logMass == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
            return logMass - logZMinEnumerated;
        }

        public boolean hasMonotoneTail(double tolerance) {
            for (int i = 1; i < thresholdsKcal.length; i++) {
                if (tailStateCounts[i] > tailStateCounts[i - 1]) {
                    return false;
                }
                double previous = logTailMass[i - 1];
                double current = logTailMass[i];
                if (current == Double.NEGATIVE_INFINITY) {
                    continue;
                }
                if (previous == Double.NEGATIVE_INFINITY
                        || current > previous + tolerance) {
                    return false;
                }
            }
            return true;
        }
    }

    /**
     * Stream the complete RC Cartesian product and exactly stratify the
     * minimizing-matrix Boltzmann mass by canonical-eta threshold.
     */
    public static OracleResult enumerate(
            RCs rcs,
            EnergyMatrix minimizingEmat,
            EnergyMatrix proposalEmat,
            InteractionGraph interactionGraph,
            double rt,
            double canonicalShiftKcal,
            double[] thresholdsKcal,
            long maxStates,
            int rebaseInterval,
            ProgressListener progressListener) {

        if (rcs == null || minimizingEmat == null || proposalEmat == null
                || interactionGraph == null) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle inputs must be non-null");
        }
        if (!Double.isFinite(rt) || !(rt > 0.0)) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle RT must be finite and positive: " + rt);
        }
        if (!Double.isFinite(canonicalShiftKcal)) {
            throw new IllegalArgumentException(
                    "eta-v5 canonical shift must be finite: "
                            + canonicalShiftKcal);
        }
        double[] thresholds = validateThresholds(thresholdsKcal);
        if (maxStates < 1L) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle maxStates must be positive: " + maxStates);
        }
        if (rebaseInterval < 1) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle rebaseInterval must be positive: "
                            + rebaseInterval);
        }
        int numPos = rcs.getNumPos();
        if (interactionGraph.getNumPositions() != numPos) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle graph/RC position mismatch: graph="
                            + interactionGraph.getNumPositions()
                            + " rcs=" + numPos);
        }

        BigInteger cartesianBig = rcs.getNumConformations();
        if (cartesianBig.signum() <= 0) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle requires a nonempty conformation space");
        }
        if (cartesianBig.compareTo(BigInteger.valueOf(maxStates)) > 0) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle Cartesian product exceeds frozen cap: states="
                            + cartesianBig + " cap=" + maxStates);
        }
        final long cartesianStates;
        try {
            cartesianStates = cartesianBig.longValueExact();
        } catch (ArithmeticException ex) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle Cartesian product does not fit in long: "
                            + cartesianBig, ex);
        }

        int[][] allowedRCs = new int[numPos][];
        int[] digits = new int[numPos];
        int[] conf = new int[numPos];
        for (int pos = 0; pos < numPos; pos++) {
            allowedRCs[pos] = rcs.get(pos).clone();
            if (allowedRCs[pos].length == 0) {
                throw new IllegalArgumentException(
                        "eta-v5 oracle position has no allowed RCs: " + pos);
            }
            conf[pos] = allowedRCs[pos][0];
        }
        int[][] neighbors = buildNeighbors(interactionGraph, numPos);

        MatrixAccumulator minimizing = new MatrixAccumulator(
                minimizingEmat, interactionGraph, neighbors, conf);
        MatrixAccumulator proposal = new MatrixAccumulator(
                proposalEmat, interactionGraph, neighbors, conf);
        LogSumExpAccumulator totalMass = new LogSumExpAccumulator();
        LogSumExpAccumulator totalProposalMass =
                new LogSumExpAccumulator();
        LogSumExpAccumulator[] massByExceededThresholdCount =
                new LogSumExpAccumulator[thresholds.length + 1];
        long[] statesByExceededThresholdCount =
                new long[thresholds.length + 1];
        for (int i = 0; i < massByExceededThresholdCount.length; i++) {
            massByExceededThresholdCount[i] =
                    new LogSumExpAccumulator();
        }

        long finiteStates = 0L;
        long zeroMassStates = 0L;
        long infiniteEtaStates = 0L;
        double rawEtaMin = Double.POSITIVE_INFINITY;
        double rawEtaMax = Double.NEGATIVE_INFINITY;
        double canonicalEtaMin = Double.POSITIVE_INFINITY;
        double canonicalEtaMax = Double.NEGATIVE_INFINITY;
        double maxMinimizingDrift = 0.0;
        double maxProposalDrift = 0.0;

        for (long stateIndex = 0L;
             stateIndex < cartesianStates; stateIndex++) {
            double minimizingEnergy = minimizing.value();
            double proposalEnergy = proposal.value();
            if (Double.isNaN(minimizingEnergy)
                    || minimizingEnergy == Double.NEGATIVE_INFINITY) {
                throw new IllegalStateException(
                        "eta-v5 oracle encountered invalid minimizing energy"
                                + " at state " + stateIndex + ": "
                                + minimizingEnergy);
            }

            if (minimizingEnergy == Double.POSITIVE_INFINITY) {
                if (proposalEnergy != Double.POSITIVE_INFINITY) {
                    throw new IllegalStateException(
                            "eta-v5 proposal resurrected a zero-mass minimizing"
                                    + " state " + stateIndex
                                    + ": proposalEnergy=" + proposalEnergy);
                }
                zeroMassStates = checkedIncrement(
                        zeroMassStates, "zero-mass state count");
            } else {
                if (Double.isNaN(proposalEnergy)
                        || proposalEnergy == Double.NEGATIVE_INFINITY) {
                    throw new IllegalStateException(
                            "eta-v5 oracle encountered invalid proposal energy"
                                    + " at finite-mass state " + stateIndex
                                    + ": " + proposalEnergy);
                }
                double rawEta = proposalEnergy == Double.POSITIVE_INFINITY
                        ? Double.POSITIVE_INFINITY
                        : proposalEnergy - minimizingEnergy;
                if (Double.isNaN(rawEta)
                        || rawEta == Double.NEGATIVE_INFINITY) {
                    throw new IllegalStateException(
                            "eta-v5 oracle encountered invalid raw eta at state "
                                    + stateIndex + ": " + rawEta);
                }
                double canonicalEta = rawEta == Double.POSITIVE_INFINITY
                        ? Double.POSITIVE_INFINITY
                        : rawEta - canonicalShiftKcal;
                if (Double.isNaN(canonicalEta)
                        || canonicalEta == Double.NEGATIVE_INFINITY) {
                    throw new IllegalStateException(
                            "eta-v5 oracle encountered invalid canonical eta"
                                    + " at state " + stateIndex + ": "
                                    + canonicalEta);
                }
                if (rawEta == Double.POSITIVE_INFINITY) {
                    infiniteEtaStates = checkedIncrement(
                            infiniteEtaStates, "infinite-eta state count");
                } else {
                    double proposalLogWeight = -proposalEnergy / rt;
                    if (!Double.isFinite(proposalLogWeight)) {
                        throw new IllegalStateException(
                                "eta-v5 proposal Boltzmann log weight is"
                                        + " non-finite at state "
                                        + stateIndex + ": "
                                        + proposalLogWeight);
                    }
                    totalProposalMass.add(proposalLogWeight);
                }
                rawEtaMin = Math.min(rawEtaMin, rawEta);
                rawEtaMax = Math.max(rawEtaMax, rawEta);
                canonicalEtaMin = Math.min(
                        canonicalEtaMin, canonicalEta);
                canonicalEtaMax = Math.max(
                        canonicalEtaMax, canonicalEta);

                double logWeight = -minimizingEnergy / rt;
                if (!Double.isFinite(logWeight)) {
                    throw new IllegalStateException(
                            "eta-v5 oracle minimizing Boltzmann log weight"
                                    + " is non-finite at state " + stateIndex
                                    + ": " + logWeight);
                }
                int exceeded = upperBoundStrict(
                        thresholds, canonicalEta);
                totalMass.add(logWeight);
                massByExceededThresholdCount[exceeded].add(logWeight);
                statesByExceededThresholdCount[exceeded] = checkedIncrement(
                        statesByExceededThresholdCount[exceeded],
                        "eta-bin state count");
                finiteStates = checkedIncrement(
                        finiteStates, "finite-mass state count");
            }

            long visited = stateIndex + 1L;
            if (visited % rebaseInterval == 0L
                    && visited < cartesianStates) {
                maxMinimizingDrift = Math.max(maxMinimizingDrift,
                        minimizing.rebase(conf));
                maxProposalDrift = Math.max(maxProposalDrift,
                        proposal.rebase(conf));
            }
            if (progressListener != null
                    && (visited == cartesianStates
                    || visited % 1_000_000L == 0L)) {
                progressListener.onProgress(visited, cartesianStates);
            }
            if (visited < cartesianStates) {
                advance(allowedRCs, digits, conf,
                        minimizing, proposal);
            }
        }

        if (finiteStates + zeroMassStates != cartesianStates) {
            throw new IllegalStateException(
                    "eta-v5 oracle state accounting mismatch: finite="
                            + finiteStates + " zeroMass=" + zeroMassStates
                            + " cartesian=" + cartesianStates);
        }
        if (finiteStates == 0L) {
            throw new IllegalStateException(
                    "eta-v5 oracle found no finite minimizing-matrix mass");
        }

        long[] tailStateCounts = new long[thresholds.length];
        double[] logTailMass = new double[thresholds.length];
        LogSumExpAccumulator suffixMass = new LogSumExpAccumulator();
        long suffixStates = 0L;
        for (int thresholdIndex = thresholds.length - 1;
             thresholdIndex >= 0; thresholdIndex--) {
            int bin = thresholdIndex + 1;
            suffixMass.addLogMass(
                    massByExceededThresholdCount[bin].logValue());
            suffixStates = checkedAdd(
                    suffixStates,
                    statesByExceededThresholdCount[bin],
                    "tail state count");
            tailStateCounts[thresholdIndex] = suffixStates;
            logTailMass[thresholdIndex] = suffixMass.logValue();
        }

        double logZProposalCanonical =
                totalProposalMass.logValue() + canonicalShiftKcal / rt;
        if (!Double.isFinite(logZProposalCanonical)) {
            throw new IllegalStateException(
                    "eta-v5 canonical proposal normalizer is non-finite: "
                            + logZProposalCanonical);
        }

        return new OracleResult(
                cartesianStates,
                cartesianStates,
                finiteStates,
                zeroMassStates,
                infiniteEtaStates,
                totalMass.logValue(),
                totalProposalMass.logValue(),
                canonicalShiftKcal,
                logZProposalCanonical,
                rawEtaMin,
                rawEtaMax,
                canonicalEtaMin,
                canonicalEtaMax,
                maxMinimizingDrift,
                maxProposalDrift,
                thresholds,
                tailStateCounts,
                logTailMass);
    }

    /**
     * Unique proposal-energy gauge that makes the canonical proposal
     * normalizer equal to the minimizing-matrix normalizer.
     */
    public static double canonicalShiftKcal(
            double rt, double logZMin, double logZProposalRaw) {
        if (!Double.isFinite(rt) || !(rt > 0.0)
                || !Double.isFinite(logZMin)
                || !Double.isFinite(logZProposalRaw)) {
            throw new IllegalArgumentException(
                    "eta-v5 canonical shift inputs must be finite and RT"
                            + " must be positive");
        }
        double shift = rt * (logZMin - logZProposalRaw);
        if (!Double.isFinite(shift)) {
            throw new IllegalArgumentException(
                    "eta-v5 canonical shift is non-finite: " + shift);
        }
        return shift;
    }

    /** Bulk integrand for the threshold decomposition. */
    public static double bulkValue(
            double etaKcal, double gapKcal,
            double tauKcal, double rt) {
        requireDecompositionInputs(etaKcal, gapKcal, tauKcal, rt);
        return Math.exp(Math.min(
                (etaKcal - gapKcal) / rt,
                tauKcal / rt));
    }

    /** Tail integrand, bounded in [0,1] whenever gapKcal is nonnegative. */
    public static double tailValue(
            double etaKcal, double gapKcal,
            double tauKcal, double rt) {
        requireDecompositionInputs(etaKcal, gapKcal, tauKcal, rt);
        double targetFactor = Math.exp(-gapKcal / rt);
        double clippedProposalFactor = etaKcal == Double.POSITIVE_INFINITY
                ? 0.0 : Math.exp((tauKcal - etaKcal) / rt);
        return Math.max(0.0, targetFactor - clippedProposalFactor);
    }

    private static void requireDecompositionInputs(
            double etaKcal, double gapKcal,
            double tauKcal, double rt) {
        if (Double.isNaN(etaKcal)
                || etaKcal == Double.NEGATIVE_INFINITY
                || !Double.isFinite(gapKcal) || gapKcal < 0.0
                || !Double.isFinite(tauKcal)
                || !Double.isFinite(rt) || !(rt > 0.0)) {
            throw new IllegalArgumentException(
                    "invalid eta-v5 threshold decomposition inputs");
        }
    }

    private static double[] validateThresholds(double[] configured) {
        if (configured == null || configured.length == 0) {
            throw new IllegalArgumentException(
                    "eta-v5 oracle requires at least one threshold");
        }
        double[] thresholds = configured.clone();
        for (int i = 0; i < thresholds.length; i++) {
            if (!Double.isFinite(thresholds[i])) {
                throw new IllegalArgumentException(
                        "eta-v5 thresholds must be finite: "
                                + Arrays.toString(configured));
            }
            if (i > 0 && !(thresholds[i] > thresholds[i - 1])) {
                throw new IllegalArgumentException(
                        "eta-v5 thresholds must be strictly increasing: "
                                + Arrays.toString(configured));
            }
        }
        return thresholds;
    }

    private static int[][] buildNeighbors(
            InteractionGraph graph, int numPos) {
        int[] counts = new int[numPos];
        for (int[] edge : graph.getEdgeList()) {
            counts[edge[0]]++;
            counts[edge[1]]++;
        }
        int[][] neighbors = new int[numPos][];
        for (int pos = 0; pos < numPos; pos++) {
            neighbors[pos] = new int[counts[pos]];
        }
        Arrays.fill(counts, 0);
        for (int[] edge : graph.getEdgeList()) {
            neighbors[edge[0]][counts[edge[0]]++] = edge[1];
            neighbors[edge[1]][counts[edge[1]]++] = edge[0];
        }
        return neighbors;
    }

    private static int upperBoundStrict(
            double[] sortedThresholds, double eta) {
        int lower = 0;
        int upper = sortedThresholds.length;
        while (lower < upper) {
            int middle = (lower + upper) >>> 1;
            if (sortedThresholds[middle] < eta) {
                lower = middle + 1;
            } else {
                upper = middle;
            }
        }
        return lower;
    }

    private static void advance(
            int[][] allowedRCs,
            int[] digits,
            int[] conf,
            MatrixAccumulator minimizing,
            MatrixAccumulator proposal) {
        for (int pos = digits.length - 1; pos >= 0; pos--) {
            int nextDigit = digits[pos] + 1;
            if (nextDigit < allowedRCs[pos].length) {
                changePosition(pos, allowedRCs[pos][nextDigit], conf,
                        minimizing, proposal);
                digits[pos] = nextDigit;
                return;
            }
            changePosition(pos, allowedRCs[pos][0], conf,
                    minimizing, proposal);
            digits[pos] = 0;
        }
        throw new IllegalStateException(
                "eta-v5 oracle mixed-radix counter overflowed early");
    }

    private static void changePosition(
            int pos, int newRC, int[] conf,
            MatrixAccumulator minimizing,
            MatrixAccumulator proposal) {
        minimizing.removePosition(pos, conf);
        proposal.removePosition(pos, conf);
        conf[pos] = newRC;
        minimizing.addPosition(pos, conf);
        proposal.addPosition(pos, conf);
    }

    private static long checkedIncrement(long value, String label) {
        if (value == Long.MAX_VALUE) {
            throw new ArithmeticException(label + " overflow");
        }
        return value + 1L;
    }

    private static long checkedAdd(long left, long right, String label) {
        if (right > 0L && left > Long.MAX_VALUE - right) {
            throw new ArithmeticException(label + " overflow");
        }
        return left + right;
    }

    private static final class MatrixAccumulator {
        private final EnergyMatrix emat;
        private final InteractionGraph graph;
        private final int[][] neighbors;
        private FactorAccumulator sum;

        MatrixAccumulator(
                EnergyMatrix emat,
                InteractionGraph graph,
                int[][] neighbors,
                int[] conf) {
            this.emat = emat;
            this.graph = graph;
            this.neighbors = neighbors;
            this.sum = build(conf);
        }

        double value() {
            return sum.value();
        }

        void removePosition(int pos, int[] conf) {
            sum.remove(emat.getOneBody(pos, conf[pos]));
            for (int neighbor : neighbors[pos]) {
                sum.remove(emat.getPairwise(
                        pos, conf[pos], neighbor, conf[neighbor]));
            }
        }

        void addPosition(int pos, int[] conf) {
            sum.add(emat.getOneBody(pos, conf[pos]));
            for (int neighbor : neighbors[pos]) {
                sum.add(emat.getPairwise(
                        pos, conf[pos], neighbor, conf[neighbor]));
            }
        }

        double rebase(int[] conf) {
            double before = value();
            FactorAccumulator rebuilt = build(conf);
            double after = rebuilt.value();
            double drift;
            if (Double.isFinite(before) && Double.isFinite(after)) {
                drift = Math.abs(before - after);
            } else if (Double.doubleToLongBits(before)
                    == Double.doubleToLongBits(after)) {
                drift = 0.0;
            } else {
                throw new IllegalStateException(
                        "eta-v5 incremental energy changed finiteness during rebase: before="
                                + before + " after=" + after);
            }
            sum = rebuilt;
            return drift;
        }

        private FactorAccumulator build(int[] conf) {
            FactorAccumulator rebuilt = new FactorAccumulator();
            rebuilt.add(emat.getConstTerm());
            for (int pos = 0; pos < conf.length; pos++) {
                rebuilt.add(emat.getOneBody(pos, conf[pos]));
            }
            for (int[] edge : graph.getEdgeList()) {
                rebuilt.add(emat.getPairwise(
                        edge[0], conf[edge[0]],
                        edge[1], conf[edge[1]]));
            }
            return rebuilt;
        }
    }

    private static final class FactorAccumulator {
        private double sum = 0.0;
        private double correction = 0.0;
        private long positiveInfinityCount = 0L;
        private long negativeInfinityCount = 0L;
        private long nanCount = 0L;

        void add(double value) {
            update(value, 1L);
        }

        void remove(double value) {
            update(value, -1L);
        }

        double value() {
            if (nanCount > 0L
                    || (positiveInfinityCount > 0L
                    && negativeInfinityCount > 0L)) {
                return Double.NaN;
            }
            if (positiveInfinityCount > 0L) {
                return Double.POSITIVE_INFINITY;
            }
            if (negativeInfinityCount > 0L) {
                return Double.NEGATIVE_INFINITY;
            }
            return sum + correction;
        }

        private void update(double value, long direction) {
            if (Double.isNaN(value)) {
                nanCount = checkedCounterUpdate(
                        nanCount, direction, "NaN factor count");
            } else if (value == Double.POSITIVE_INFINITY) {
                positiveInfinityCount = checkedCounterUpdate(
                        positiveInfinityCount, direction,
                        "positive-infinity factor count");
            } else if (value == Double.NEGATIVE_INFINITY) {
                negativeInfinityCount = checkedCounterUpdate(
                        negativeInfinityCount, direction,
                        "negative-infinity factor count");
            } else {
                addFinite(direction > 0L ? value : -value);
            }
        }

        private void addFinite(double value) {
            double updated = sum + value;
            if (Math.abs(sum) >= Math.abs(value)) {
                correction += (sum - updated) + value;
            } else {
                correction += (value - updated) + sum;
            }
            sum = updated;
        }

        private static long checkedCounterUpdate(
                long current, long direction, String label) {
            if (direction > 0L) {
                return checkedIncrement(current, label);
            }
            if (current <= 0L) {
                throw new IllegalStateException(label + " underflow");
            }
            return current - 1L;
        }
    }

    private static final class LogSumExpAccumulator {
        private double maxLog = Double.NEGATIVE_INFINITY;
        private double scaledSum = 0.0;
        private double correction = 0.0;

        void add(double logValue) {
            if (!Double.isFinite(logValue)) {
                throw new IllegalArgumentException(
                        "log-sum-exp input must be finite: " + logValue);
            }
            if (maxLog == Double.NEGATIVE_INFINITY) {
                maxLog = logValue;
                scaledSum = 1.0;
                correction = 0.0;
            } else if (logValue > maxLog) {
                double oldTotal = scaledSum + correction;
                scaledSum = oldTotal * Math.exp(maxLog - logValue) + 1.0;
                correction = 0.0;
                maxLog = logValue;
            } else {
                addScaled(Math.exp(logValue - maxLog));
            }
        }

        void addLogMass(double logMass) {
            if (logMass == Double.NEGATIVE_INFINITY) {
                return;
            }
            add(logMass);
        }

        double logValue() {
            if (maxLog == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
            double total = scaledSum + correction;
            if (!Double.isFinite(total) || !(total > 0.0)) {
                throw new IllegalStateException(
                        "invalid log-sum-exp scaled total: " + total);
            }
            return maxLog + Math.log(total);
        }

        private void addScaled(double value) {
            double updated = scaledSum + value;
            if (Math.abs(scaledSum) >= Math.abs(value)) {
                correction += (scaledSum - updated) + value;
            } else {
                correction += (value - updated) + scaledSum;
            }
            scaledSum = updated;
        }
    }
}
