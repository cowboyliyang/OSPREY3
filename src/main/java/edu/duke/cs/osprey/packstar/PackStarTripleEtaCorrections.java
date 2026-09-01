package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Sparse signed three-position factors used only to define a PACK* proposal
 * energy.
 *
 * <p>The initial local prior partially minimizes each retained graph-clique's
 * three pair terms with the historical {@code 1/(numPos-2)} scale.  It makes
 * zero full-conformation CCD calls.  The prior is not added wholesale.  Each
 * adaptive round instead fits a small number of signed RC-triple tables to
 * {@code E_true-E_m-eta_single_pair}, with the local value used only as a
 * shrinkage target.  Consequently these factors are residuals by definition:
 * no pair-incidence allocation is needed.</p>
 *
 * <p>Selected factors may be non-cliques.  The estimator must add zero-energy
 * fill edges and rebuild the proposal-only branch decomposition before using
 * them.  {@link #requiredFillEdges(InteractionGraph)} exposes that exact
 * requirement.</p>
 *
 * <p>This is a proposal model, not a claimed lower bound.  Final PACK* weights
 * still use full-conformation CCD energies.</p>
 */
final class PackStarTripleEtaCorrections {

    /**
     * Signals that the optional triple table was rejected only because its
     * enumerated workload exceeded the configured capacity.  The estimator
     * may safely fall back to its already-supported pair-only proposal for
     * this case; other construction failures must remain fatal.
     */
    static final class AssignmentCapExceededException
            extends IllegalArgumentException {
        private static final long serialVersionUID = 1L;

        final long assignments;
        final long maximumAssignments;

        AssignmentCapExceededException(
                long assignments, long maximumAssignments) {
            super("triple eta clique workload exceeds cap: assignments="
                    + assignments + " cap=" + maximumAssignments);
            this.assignments = assignments;
            this.maximumAssignments = maximumAssignments;
        }
    }

    @FunctionalInterface
    interface PairEtaLookup {
        double getPairEta(
                int pos1, int rc1, int pos2, int rc2);
    }

    static final class ResidualSummary {
        final long positiveFactors;
        final long negativeFactors;
        final long zeroFactors;
        final double minimumCorrectionKcal;
        final double maximumCorrectionKcal;
        final double maximumAbsoluteCorrectionKcal;

        ResidualSummary(
                long positiveFactors,
                long negativeFactors,
                long zeroFactors,
                double minimumCorrectionKcal,
                double maximumCorrectionKcal,
                double maximumAbsoluteCorrectionKcal) {
            this.positiveFactors = positiveFactors;
            this.negativeFactors = negativeFactors;
            this.zeroFactors = zeroFactors;
            this.minimumCorrectionKcal = minimumCorrectionKcal;
            this.maximumCorrectionKcal = maximumCorrectionKcal;
            this.maximumAbsoluteCorrectionKcal =
                    maximumAbsoluteCorrectionKcal;
        }

        boolean hasNonzeroFactors() {
            return positiveFactors > 0L || negativeFactors > 0L;
        }
    }

    static final String TABLE_SCHEMA =
            "packstar-selected-triple-eta-table-v3";
    static final String SUMMARY_SCHEMA =
            "packstar-selected-triple-eta-summary-v3";

    static final class Entry {
        final long index;
        final int pos1;
        final int rc1;
        final int pos2;
        final int rc2;
        final int pos3;
        final int rc3;
        final double pairwiseOffsetKcal;
        volatile double partialMinimumKcal = Double.NaN;
        volatile double rawCorrectionKcal = Double.NaN;
        volatile double storedCorrectionKcal = 0.0;
        volatile boolean finite = false;

        Entry(long index,
              int pos1, int rc1,
              int pos2, int rc2,
              int pos3, int rc3,
              double pairwiseOffsetKcal) {
            this.index = index;
            this.pos1 = pos1;
            this.rc1 = rc1;
            this.pos2 = pos2;
            this.rc2 = rc2;
            this.pos3 = pos3;
            this.rc3 = rc3;
            this.pairwiseOffsetKcal = pairwiseOffsetKcal;
        }

        void complete(double partialMinimumKcal) {
            this.partialMinimumKcal = partialMinimumKcal;
            this.rawCorrectionKcal = partialMinimumKcal
                    - pairwiseOffsetKcal;
            this.finite = Double.isFinite(pairwiseOffsetKcal)
                    && Double.isFinite(partialMinimumKcal)
                    && Double.isFinite(rawCorrectionKcal);
            // Keep the sign.  A tiny negative local value can arise from
            // numerical minimization noise, and fitted residual tables must
            // be able to repair pair over-correction with a negative factor.
            this.storedCorrectionKcal = finite
                    ? rawCorrectionKcal : 0.0;
        }

        static Entry fitted(
                long index,
                int pos1, int rc1,
                int pos2, int rc2,
                int pos3, int rc3,
                double correctionKcal) {
            Entry entry = new Entry(index,
                    pos1, rc1, pos2, rc2, pos3, rc3, 0.0);
            entry.partialMinimumKcal = correctionKcal;
            entry.rawCorrectionKcal = correctionKcal;
            entry.storedCorrectionKcal = correctionKcal;
            entry.finite = Double.isFinite(correctionKcal);
            return entry;
        }

        TripleKey key() {
            return new TripleKey(
                    pos1, rc1, pos2, rc2, pos3, rc3);
        }

        RCTuple tuple() {
            return new RCTuple(
                    pos3, rc3, pos2, rc2, pos1, rc1);
        }
    }

    static final class PositionTriple
            implements Comparable<PositionTriple> {
        final int pos1;
        final int pos2;
        final int pos3;

        PositionTriple(int a, int b, int c) {
            int[] sorted = {a, b, c};
            Arrays.sort(sorted);
            if (sorted[0] == sorted[1]
                    || sorted[1] == sorted[2]) {
                throw new IllegalArgumentException(
                        "triple positions must be distinct");
            }
            this.pos1 = sorted[2];
            this.pos2 = sorted[1];
            this.pos3 = sorted[0];
        }

        long assignmentCount(RCs rcs) {
            long count = Math.multiplyExact(
                    (long) rcs.getNum(pos1),
                    (long) rcs.getNum(pos2));
            return Math.multiplyExact(count,
                    (long) rcs.getNum(pos3));
        }

        List<Long> requiredFillKeys(InteractionGraph graph) {
            List<Long> keys = new ArrayList<>(3);
            addFillKey(keys, graph, pos1, pos2);
            addFillKey(keys, graph, pos1, pos3);
            addFillKey(keys, graph, pos2, pos3);
            return keys;
        }

        private static void addFillKey(
                List<Long> keys, InteractionGraph graph,
                int a, int b) {
            if (!graph.hasEdge(a, b)) {
                keys.add(packPositionPair(a, b));
            }
        }

        @Override
        public int compareTo(PositionTriple other) {
            int compare = Integer.compare(pos1, other.pos1);
            if (compare != 0) return compare;
            compare = Integer.compare(pos2, other.pos2);
            if (compare != 0) return compare;
            return Integer.compare(pos3, other.pos3);
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (!(other instanceof PositionTriple)) return false;
            PositionTriple rhs = (PositionTriple) other;
            return pos1 == rhs.pos1 && pos2 == rhs.pos2
                    && pos3 == rhs.pos3;
        }

        @Override
        public int hashCode() {
            return Objects.hash(pos1, pos2, pos3);
        }

        @Override
        public String toString() {
            return pos3 + "," + pos2 + "," + pos1;
        }
    }

    private static final class CellStats {
        double sum;
        int contexts;

        void add(double value) {
            sum += value;
            contexts++;
        }

        double mean() {
            return contexts > 0 ? sum / contexts : 0.0;
        }
    }

    private static final class UniqueSample {
        final int[] conf;
        double residualSum;
        int multiplicity;

        UniqueSample(int[] conf, double residual) {
            this.conf = conf.clone();
            this.residualSum = residual;
            this.multiplicity = 1;
        }

        void add(double residual) {
            residualSum += residual;
            multiplicity++;
        }

        double residual() {
            return residualSum / multiplicity;
        }
    }

    private static final class CandidateFit {
        final PositionTriple triple;
        final Map<TripleKey, Double> values;
        final double[] predictions;
        final double meanSquaredGain;
        final double priority;
        final Set<Long> newFillEdges;
        final long assignments;

        CandidateFit(PositionTriple triple,
                     Map<TripleKey, Double> values,
                     double[] predictions,
                     double meanSquaredGain,
                     double priority,
                     Set<Long> newFillEdges,
                     long assignments) {
            this.triple = triple;
            this.values = values;
            this.predictions = predictions;
            this.meanSquaredGain = meanSquaredGain;
            this.priority = priority;
            this.newFillEdges = newFillEdges;
            this.assignments = assignments;
        }
    }

    final int numPositions;
    final double pairWeight;
    final long positionTripleCount;
    final long factorAssignments;
    final long positiveFactors;
    final long nonpositiveFactors;
    final long nonfiniteFactors;
    final double maximumStoredCorrectionKcal;
    final double wallSeconds;
    final String signatureSha256;
    final List<Entry> entries;
    private final Map<TripleKey, Entry> entriesByKey;
    private final List<PositionTriple> positionTriples;

    PackStarTripleEtaCorrections(
            int numPositions,
            double pairWeight,
            long positionTripleCount,
            List<Entry> entries,
            double wallSeconds) {
        this.numPositions = numPositions;
        this.pairWeight = pairWeight;
        this.positionTripleCount = positionTripleCount;
        this.factorAssignments = entries.size();
        this.entries = List.copyOf(entries);
        this.wallSeconds = wallSeconds;
        this.entriesByKey = new HashMap<>();
        Set<PositionTriple> uniquePositionTriples =
                new LinkedHashSet<>();
        long positive = 0L;
        long nonpositive = 0L;
        long nonfinite = 0L;
        double maximum = 0.0;
        for (Entry entry : entries) {
            if (!entry.finite) {
                nonfinite++;
            } else {
                entriesByKey.put(entry.key(), entry);
                uniquePositionTriples.add(new PositionTriple(
                        entry.pos1, entry.pos2, entry.pos3));
                if (entry.storedCorrectionKcal > 0.0) {
                    positive++;
                    maximum = Math.max(maximum,
                            entry.storedCorrectionKcal);
                } else {
                    nonpositive++;
                }
            }
        }
        this.positiveFactors = positive;
        this.nonpositiveFactors = nonpositive;
        this.nonfiniteFactors = nonfinite;
        this.maximumStoredCorrectionKcal = maximum;
        List<PositionTriple> sortedTriples =
                new ArrayList<>(uniquePositionTriples);
        Collections.sort(sortedTriples);
        this.positionTriples = List.copyOf(sortedTriples);
        if (positionTripleCount != this.positionTriples.size()) {
            throw new IllegalArgumentException(
                    "triple eta position-scope count mismatch: declared="
                            + positionTripleCount + " observed="
                            + this.positionTriples.size());
        }
        this.signatureSha256 = computeSignature(entries);
    }

    static long countCliquePositionTriples(
            int numPositions, InteractionGraph graph) {
        long count = 0L;
        for (int pos1 = 2; pos1 < numPositions; pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    if (isClique(graph, pos1, pos2, pos3)) {
                        count = Math.addExact(count, 1L);
                    }
                }
            }
        }
        return count;
    }

    static long countCliqueAssignments(RCs rcs, InteractionGraph graph) {
        validateShape(rcs, graph);
        long count = 0L;
        for (int pos1 = 2; pos1 < rcs.getNumPos(); pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    if (!isClique(graph, pos1, pos2, pos3)) continue;
                    long assignments = Math.multiplyExact(
                            (long) rcs.getNum(pos1),
                            (long) rcs.getNum(pos2));
                    assignments = Math.multiplyExact(assignments,
                            (long) rcs.getNum(pos3));
                    count = Math.addExact(count, assignments);
                }
            }
        }
        return count;
    }

    static void ensureAssignmentCapacity(
            long assignments, long maximumAssignments) {
        if (maximumAssignments < 1L) {
            throw new IllegalArgumentException(
                    "triple eta assignment cap must be positive");
        }
        if (assignments > maximumAssignments) {
            throw new AssignmentCapExceededException(
                    assignments, maximumAssignments);
        }
    }

    static PackStarTripleEtaCorrections compute(
            RCs rcs,
            ConfEnergyCalculator confEcalc,
            InteractionGraph graph,
            EnergyMatrix minimizingEmat,
            long maximumAssignments) {
        long startNanos = System.nanoTime();
        if (confEcalc == null || minimizingEmat == null) {
            throw new IllegalArgumentException(
                    "triple eta inputs must be non-null");
        }
        validateShape(rcs, graph);
        if (minimizingEmat.getNumPos() != rcs.getNumPos()) {
            throw new IllegalArgumentException(
                    "triple eta matrix/RC position mismatch");
        }
        if (minimizingEmat.hasHigherOrderTerms()) {
            throw new IllegalArgumentException(
                    "triple eta requires a pairwise base energy matrix");
        }
        if (rcs.getNumPos() < 3) {
            return new PackStarTripleEtaCorrections(
                    rcs.getNumPos(), 0.0, 0L,
                    List.of(), 0.0);
        }
        if (maximumAssignments < 1L) {
            throw new IllegalArgumentException(
                    "triple eta assignment cap must be positive");
        }
        long expected = countCliqueAssignments(rcs, graph);
        ensureAssignmentCapacity(expected, maximumAssignments);
        if (expected > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(
                    "triple eta clique workload exceeds in-memory table limit: "
                            + expected);
        }
        long cliqueTriples = countCliquePositionTriples(
                rcs.getNumPos(), graph);
        double pairWeight = 1.0 / (rcs.getNumPos() - 2.0);
        List<Entry> entries = new ArrayList<>((int) expected);
        AtomicLong completed = new AtomicLong();

        System.out.println("[PACK*-triple-eta] cliquePositionTriples="
                + cliqueTriples + ", factorAssignments=" + expected
                + ", pairWeight="
                + String.format(Locale.ROOT, "%.9g", pairWeight)
                + ", additionalFullConformationCcdCalls=0");

        for (int pos1 = 2; pos1 < rcs.getNumPos(); pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    if (!isClique(graph, pos1, pos2, pos3)) continue;
                    for (int rc1 : rcs.get(pos1)) {
                        for (int rc2 : rcs.get(pos2)) {
                            for (int rc3 : rcs.get(pos3)) {
                                double offset = pairWeight * (
                                        minimizingEmat.getPairwise(
                                                pos1, rc1, pos2, rc2)
                                                + minimizingEmat.getPairwise(
                                                pos1, rc1, pos3, rc3)
                                                + minimizingEmat.getPairwise(
                                                pos2, rc2, pos3, rc3));
                                Entry entry = new Entry(
                                        entries.size(),
                                        pos1, rc1, pos2, rc2,
                                        pos3, rc3, offset);
                                entries.add(entry);

                                ResInterGen generator = ResInterGen.of(
                                        confEcalc.confSpace);
                                generator.addInter(pos1, pos2,
                                        pairWeight, 0.0);
                                generator.addInter(pos1, pos3,
                                        pairWeight, 0.0);
                                generator.addInter(pos2, pos3,
                                        pairWeight, 0.0);
                                ResidueInteractions interactions =
                                        generator.make();
                                RCTuple tuple = entry.tuple();
                                confEcalc.tasks.submit(
                                        () -> confEcalc.calcEnergy(
                                                tuple, interactions).energy,
                                        energy -> {
                                            entry.complete(energy);
                                            long done = completed.incrementAndGet();
                                            if (done % 5000L == 0L
                                                    || done == expected) {
                                                System.out.println(
                                                        "[PACK*-triple-eta] partialCcdCompleted="
                                                                + done + "/" + expected);
                                            }
                                        });
                            }
                        }
                    }
                }
            }
        }
        if (entries.size() != expected) {
            throw new IllegalStateException(
                    "triple eta enumeration mismatch: expected="
                            + expected + " observed=" + entries.size());
        }
        confEcalc.tasks.waitForFinish();
        if (completed.get() != expected) {
            throw new IllegalStateException(
                    "triple eta partial-CCD completion mismatch: expected="
                            + expected + " observed=" + completed.get());
        }
        double wallSeconds = (System.nanoTime() - startNanos) / 1.0e9;
        PackStarTripleEtaCorrections result =
                new PackStarTripleEtaCorrections(
                        rcs.getNumPos(), pairWeight, cliqueTriples,
                        entries, wallSeconds);
        if (result.nonfiniteFactors > 0L) {
            throw new IllegalStateException(
                    "triple eta produced non-finite factors: "
                            + result.nonfiniteFactors + "/"
                            + result.factorAssignments);
        }
        System.out.println("[PACK*-triple-eta] complete: positiveFactors="
                + result.positiveFactors + "/" + result.factorAssignments
                + ", maxCorrectionKcal="
                + String.format(Locale.ROOT, "%.9g",
                result.maximumStoredCorrectionKcal)
                + ", wallSeconds="
                + String.format(Locale.ROOT, "%.3f", wallSeconds)
                + ", signature=" + result.signatureSha256);
        return result;
    }

    /**
     * Fit a small signed residual expansion.  Exact duplicate full
     * assignments are collapsed before fitting, so repeated IID draws do not
     * masquerade as independent contexts for one sparse RC-triple cell.
     */
    PackStarTripleEtaCorrections fitSelectedResidual(
            RCs rcs,
            InteractionGraph baseGraph,
            int[][] conformations,
            double[] residualsKcal,
            int maximumPositionTriples,
            int maximumFillEdges,
            int minimumCellContexts,
            double priorStrength,
            double localPriorCapKcal,
            double residualCapKcal,
            long maximumAssignments) {
        validateShape(rcs, baseGraph);
        if (conformations == null || residualsKcal == null
                || conformations.length != residualsKcal.length) {
            throw new IllegalArgumentException(
                    "triple residual samples have inconsistent shapes");
        }
        if (maximumPositionTriples < 0 || maximumFillEdges < 0
                || minimumCellContexts < 1
                || !Double.isFinite(priorStrength)
                || priorStrength < 0.0
                || !Double.isFinite(localPriorCapKcal)
                || localPriorCapKcal < 0.0
                || !Double.isFinite(residualCapKcal)
                || !(residualCapKcal > 0.0)
                || maximumAssignments < 1L) {
            throw new IllegalArgumentException(
                    "invalid selected-triple fitting controls");
        }
        if (maximumPositionTriples == 0
                || conformations.length < 2) {
            return emptyFitted(rcs.getNumPos());
        }

        LinkedHashMap<String, UniqueSample> uniqueByAssignment =
                new LinkedHashMap<>();
        for (int sample = 0; sample < conformations.length; sample++) {
            int[] conf = conformations[sample];
            double residual = residualsKcal[sample];
            if (conf == null || conf.length != rcs.getNumPos()
                    || !Double.isFinite(residual)) {
                throw new IllegalArgumentException(
                        "invalid triple residual sample " + sample);
            }
            String key = Arrays.toString(conf);
            UniqueSample prior = uniqueByAssignment.get(key);
            if (prior == null) {
                uniqueByAssignment.put(key,
                        new UniqueSample(conf, residual));
            } else {
                prior.add(residual);
            }
        }
        List<UniqueSample> samples =
                new ArrayList<>(uniqueByAssignment.values());
        if (samples.size() < 2) {
            return emptyFitted(rcs.getNumPos());
        }

        double residualMean = 0.0;
        for (UniqueSample sample : samples) {
            residualMean += sample.residual();
        }
        residualMean /= samples.size();
        double[] workingResidual = new double[samples.size()];
        for (int sample = 0; sample < samples.size(); sample++) {
            workingResidual[sample] = samples.get(sample).residual()
                    - residualMean;
        }

        List<PositionTriple> candidates = new ArrayList<>();
        for (int pos1 = 2; pos1 < rcs.getNumPos(); pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    candidates.add(new PositionTriple(
                            pos1, pos2, pos3));
                }
            }
        }

        Set<PositionTriple> selected = new LinkedHashSet<>();
        Set<Long> selectedFillEdges = new LinkedHashSet<>();
        List<Entry> fittedEntries = new ArrayList<>();
        long selectedAssignments = 0L;
        for (int step = 0;
             step < maximumPositionTriples; step++) {
            CandidateFit best = null;
            for (PositionTriple candidate : candidates) {
                if (selected.contains(candidate)) continue;
                long assignments = candidate.assignmentCount(rcs);
                if (assignments > maximumAssignments
                        - selectedAssignments) continue;
                Set<Long> newFillEdges = new LinkedHashSet<>(
                        candidate.requiredFillKeys(baseGraph));
                newFillEdges.removeAll(selectedFillEdges);
                if (selectedFillEdges.size() + newFillEdges.size()
                        > maximumFillEdges) continue;
                CandidateFit fit = fitCandidate(
                        rcs, candidate, samples, workingResidual,
                        minimumCellContexts, priorStrength,
                        localPriorCapKcal, residualCapKcal,
                        newFillEdges, assignments);
                if (fit == null) continue;
                if (best == null
                        || fit.priority > best.priority + 1.0e-12
                        || (Math.abs(fit.priority - best.priority)
                        <= 1.0e-12
                        && fit.triple.compareTo(best.triple) < 0)) {
                    best = fit;
                }
            }
            if (best == null || !(best.meanSquaredGain > 1.0e-8)) {
                break;
            }

            selected.add(best.triple);
            selectedFillEdges.addAll(best.newFillEdges);
            selectedAssignments = Math.addExact(
                    selectedAssignments, best.assignments);
            for (int sample = 0; sample < workingResidual.length;
                 sample++) {
                workingResidual[sample] -= best.predictions[sample];
            }
            appendFittedEntries(
                    fittedEntries, rcs, best.triple, best.values);
            System.out.println("[PACK*-triple-eta-fit] selected="
                    + best.triple
                    + ", step=" + (step + 1)
                    + ", meanSquaredGain="
                    + String.format(Locale.ROOT, "%.9g",
                    best.meanSquaredGain)
                    + ", newFillEdges=" + best.newFillEdges.size()
                    + ", totalFillEdges=" + selectedFillEdges.size()
                    + ", assignments=" + best.assignments);
        }

        return new PackStarTripleEtaCorrections(
                rcs.getNumPos(), 0.0, selected.size(),
                fittedEntries, 0.0);
    }

    private CandidateFit fitCandidate(
            RCs rcs,
            PositionTriple triple,
            List<UniqueSample> samples,
            double[] workingResidual,
            int minimumCellContexts,
            double priorStrength,
            double localPriorCapKcal,
            double residualCapKcal,
            Set<Long> newFillEdges,
            long assignments) {
        Map<TripleKey, CellStats> stats = new HashMap<>();
        Set<TripleKey> distinctObservedCells = new HashSet<>();
        for (int sample = 0; sample < samples.size(); sample++) {
            TripleKey key = keyFor(
                    triple, samples.get(sample).conf);
            stats.computeIfAbsent(key, ignored -> new CellStats())
                    .add(workingResidual[sample]);
            distinctObservedCells.add(key);
        }
        if (distinctObservedCells.size() < 2) return null;

        Map<TripleKey, Double> values = new HashMap<>();
        for (int rc1 : rcs.get(triple.pos1)) {
            for (int rc2 : rcs.get(triple.pos2)) {
                for (int rc3 : rcs.get(triple.pos3)) {
                    TripleKey key = new TripleKey(
                            triple.pos1, rc1,
                            triple.pos2, rc2,
                            triple.pos3, rc3);
                    Entry localEntry = entriesByKey.get(key);
                    double localPrior = localEntry == null
                            ? 0.0 : clamp(
                            localEntry.storedCorrectionKcal,
                            -localPriorCapKcal,
                            localPriorCapKcal);
                    CellStats cell = stats.get(key);
                    double value = localPrior;
                    if (cell != null
                            && cell.contexts >= minimumCellContexts) {
                        value = (cell.contexts * cell.mean()
                                + priorStrength * localPrior)
                                / (cell.contexts + priorStrength);
                    }
                    values.put(key, clamp(value,
                            -residualCapKcal, residualCapKcal));
                }
            }
        }

        double[] predictions = predictions(
                triple, samples, values);
        centerAndClip(values, triple, samples, predictions,
                residualCapKcal);
        predictions = predictions(triple, samples, values);

        double sseBefore = 0.0;
        double sseAfter = 0.0;
        for (int sample = 0; sample < samples.size(); sample++) {
            double before = workingResidual[sample];
            double after = before - predictions[sample];
            sseBefore += before * before;
            sseAfter += after * after;
        }
        double gain = (sseBefore - sseAfter) / samples.size();
        if (!Double.isFinite(gain) || !(gain > 0.0)) return null;
        double priority = gain / (1.0 + newFillEdges.size());
        return new CandidateFit(
                triple, values, predictions, gain, priority,
                Set.copyOf(newFillEdges), assignments);
    }

    private static double[] predictions(
            PositionTriple triple,
            List<UniqueSample> samples,
            Map<TripleKey, Double> values) {
        double[] predictions = new double[samples.size()];
        for (int sample = 0; sample < samples.size(); sample++) {
            predictions[sample] = values.getOrDefault(
                    keyFor(triple, samples.get(sample).conf), 0.0);
        }
        return predictions;
    }

    private static void centerAndClip(
            Map<TripleKey, Double> values,
            PositionTriple triple,
            List<UniqueSample> samples,
            double[] predictions,
            double cap) {
        // Two passes are sufficient to make the empirical gauge nearly zero
        // even when the first subtraction hits the trust-region cap.
        for (int pass = 0; pass < 2; pass++) {
            double center = 0.0;
            for (double prediction : predictions) center += prediction;
            center /= predictions.length;
            for (Map.Entry<TripleKey, Double> entry : values.entrySet()) {
                entry.setValue(clamp(entry.getValue() - center,
                        -cap, cap));
            }
            predictions = predictions(triple, samples, values);
        }
    }

    private static void appendFittedEntries(
            List<Entry> output,
            RCs rcs,
            PositionTriple triple,
            Map<TripleKey, Double> values) {
        for (int rc1 : rcs.get(triple.pos1)) {
            for (int rc2 : rcs.get(triple.pos2)) {
                for (int rc3 : rcs.get(triple.pos3)) {
                    TripleKey key = new TripleKey(
                            triple.pos1, rc1,
                            triple.pos2, rc2,
                            triple.pos3, rc3);
                    double value = values.getOrDefault(key, 0.0);
                    output.add(Entry.fitted(output.size(),
                            triple.pos1, rc1,
                            triple.pos2, rc2,
                            triple.pos3, rc3,
                            value));
                }
            }
        }
    }

    private static TripleKey keyFor(
            PositionTriple triple, int[] conf) {
        return new TripleKey(
                triple.pos1, conf[triple.pos1],
                triple.pos2, conf[triple.pos2],
                triple.pos3, conf[triple.pos3]);
    }

    private static PackStarTripleEtaCorrections emptyFitted(
            int numPositions) {
        return new PackStarTripleEtaCorrections(
                numPositions, 0.0, 0L,
                List.of(), 0.0);
    }

    List<PositionTriple> positionTriples() {
        return positionTriples;
    }

    List<int[]> requiredFillEdges(InteractionGraph graph) {
        if (graph == null
                || graph.getNumPositions() != numPositions) {
            throw new IllegalArgumentException(
                    "triple eta fill-edge graph mismatch");
        }
        Set<Long> packed = new LinkedHashSet<>();
        for (PositionTriple triple : positionTriples) {
            packed.addAll(triple.requiredFillKeys(graph));
        }
        List<Long> sorted = new ArrayList<>(packed);
        Collections.sort(sorted);
        List<int[]> edges = new ArrayList<>(sorted.size());
        for (long key : sorted) {
            edges.add(new int[]{
                    (int) (key >>> 32), (int) key});
        }
        return edges;
    }

    private static long packPositionPair(int pos1, int pos2) {
        int lower = Math.min(pos1, pos2);
        int upper = Math.max(pos1, pos2);
        return ((long) lower << 32) | (upper & 0xffffffffL);
    }

    private static double clamp(double value, double lower, double upper) {
        return Math.max(lower, Math.min(upper, value));
    }

    double scoreJointCorrection(int[] conf) {
        if (conf == null || conf.length != numPositions) {
            throw new IllegalArgumentException(
                    "triple eta conformation length mismatch");
        }
        double score = 0.0;
        for (int pos1 = 2; pos1 < conf.length; pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    Entry entry = entriesByKey.get(
                            new TripleKey(
                                    pos1, conf[pos1],
                                    pos2, conf[pos2],
                                    pos3, conf[pos3]));
                    if (entry != null) {
                        score += entry.storedCorrectionKcal;
                    }
                }
            }
        }
        return score;
    }

    double scoreResidual(int[] conf, PairEtaLookup pairEta) {
        if (conf == null || conf.length != numPositions) {
            throw new IllegalArgumentException(
                    "triple eta conformation length mismatch");
        }
        Objects.requireNonNull(pairEta, "triple eta pair lookup");
        double score = 0.0;
        for (int pos1 = 2; pos1 < conf.length; pos1++) {
            for (int pos2 = 1; pos2 < pos1; pos2++) {
                for (int pos3 = 0; pos3 < pos2; pos3++) {
                    Entry entry = entriesByKey.get(
                            new TripleKey(
                                    pos1, conf[pos1],
                                    pos2, conf[pos2],
                                    pos3, conf[pos3]));
                    if (entry != null) {
                        score += residualCorrection(entry, pairEta);
                    }
                }
            }
        }
        return score;
    }

    double residualCorrection(Entry entry, PairEtaLookup pairEta) {
        Objects.requireNonNull(entry, "triple eta entry");
        Objects.requireNonNull(pairEta, "triple eta pair lookup");
        // The selected table was fitted after single/pair eta was removed.
        // It is already a residual and must not subtract any allocated share
        // of pair eta here.  The old 1/(n-2) subtraction was only exact for a
        // complete triangle family and caused a large sparse-incidence bias.
        double residual = entry.storedCorrectionKcal;
        if (!Double.isFinite(residual)) {
            throw new IllegalArgumentException(
                    "triple eta residual is non-finite at assignment "
                            + entry.index + ": correction=" + residual);
        }
        return residual;
    }

    ResidualSummary summarizeResidual(PairEtaLookup pairEta) {
        Objects.requireNonNull(pairEta, "triple eta pair lookup");
        long positive = 0L;
        long negative = 0L;
        long zero = 0L;
        double minimum = Double.POSITIVE_INFINITY;
        double maximum = Double.NEGATIVE_INFINITY;
        double maximumAbsolute = 0.0;
        for (Entry entry : entries) {
            if (!entry.finite) continue;
            double residual = residualCorrection(entry, pairEta);
            if (residual > 0.0) {
                positive++;
            } else if (residual < 0.0) {
                negative++;
            } else {
                zero++;
            }
            minimum = Math.min(minimum, residual);
            maximum = Math.max(maximum, residual);
            maximumAbsolute = Math.max(
                    maximumAbsolute, Math.abs(residual));
        }
        if (positive + negative + zero == 0L) {
            minimum = 0.0;
            maximum = 0.0;
        }
        return new ResidualSummary(
                positive, negative, zero,
                minimum, maximum, maximumAbsolute);
    }

    void applyResidualTo(
            EnergyMatrix target,
            double scale,
            PairEtaLookup pairEta) {
        if (target == null || target.getNumPos() != numPositions) {
            throw new IllegalArgumentException(
                    "triple eta target matrix shape mismatch");
        }
        if (target.hasHigherOrderTerms()) {
            throw new IllegalArgumentException(
                    "triple eta target already has higher-order terms");
        }
        if (!Double.isFinite(scale)) {
            throw new IllegalArgumentException(
                    "triple eta matrix scale must be finite: "
                            + scale);
        }
        Objects.requireNonNull(pairEta, "triple eta pair lookup");
        if (scale == 0.0) return;
        for (Entry entry : entries) {
            if (entry.finite) {
                double correction = scale
                        * residualCorrection(entry, pairEta);
                if (correction == 0.0) continue;
                if (!Double.isFinite(correction)) {
                    throw new IllegalArgumentException(
                            "scaled triple eta residual is non-finite at"
                                    + " assignment " + entry.index);
                }
                target.setHigherOrder(entry.tuple(),
                        correction);
            }
        }
    }

    void writeArtifacts(File tableOutput, File summaryOutput)
            throws IOException {
        requireOutput(tableOutput);
        requireOutput(summaryOutput);
        if (tableOutput.getCanonicalFile().equals(
                summaryOutput.getCanonicalFile())) {
            throw new IllegalArgumentException(
                    "triple eta table and summary outputs must differ");
        }
        writeAtomically(tableOutput.toPath(), writer -> {
            writer.write("schema\tindex\tpos1\trc1\tpos2\trc2"
                    + "\tpos3\trc3\tpairWeight"
                    + "\tpairwiseOffsetKcal\tpartialMinimumKcal"
                    + "\trawCorrectionKcal\tstoredCorrectionKcal\tfinite");
            writer.newLine();
            for (Entry entry : entries) {
                writer.write(String.format(Locale.ROOT,
                        "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
                                + "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%s%n",
                        TABLE_SCHEMA, entry.index,
                        entry.pos1, entry.rc1,
                        entry.pos2, entry.rc2,
                        entry.pos3, entry.rc3,
                        pairWeight, entry.pairwiseOffsetKcal,
                        entry.partialMinimumKcal,
                        entry.rawCorrectionKcal,
                        entry.storedCorrectionKcal,
                        Boolean.toString(entry.finite)));
            }
        });
        writeAtomically(summaryOutput.toPath(), writer -> {
            writer.write("schema\tnumPositions\tpositionTriples"
                    + "\tfactorAssignments\tpositiveFactors"
                    + "\tnonpositiveFactors\tnonfiniteFactors"
                    + "\tpairWeight\tmaximumStoredCorrectionKcal"
                    + "\tadditionalFullConformationCcdCalls\twallSeconds"
                    + "\tsignatureSha256");
            writer.newLine();
            writer.write(String.format(Locale.ROOT,
                    "%s\t%d\t%d\t%d\t%d\t%d\t%d"
                            + "\t%.17g\t%.17g\t0\t%.17g\t%s%n",
                    SUMMARY_SCHEMA, numPositions,
                    positionTripleCount, factorAssignments,
                    positiveFactors, nonpositiveFactors,
                    nonfiniteFactors, pairWeight,
                    maximumStoredCorrectionKcal, wallSeconds,
                    signatureSha256));
        });
    }

    private static boolean isClique(
            InteractionGraph graph, int pos1, int pos2, int pos3) {
        return graph.hasEdge(pos1, pos2)
                && graph.hasEdge(pos1, pos3)
                && graph.hasEdge(pos2, pos3);
    }

    private static void validateShape(RCs rcs, InteractionGraph graph) {
        if (rcs == null || graph == null) {
            throw new IllegalArgumentException(
                    "triple eta RCs and graph must be non-null");
        }
        if (rcs.getNumPos() != graph.getNumPositions()) {
            throw new IllegalArgumentException(
                    "triple eta graph/RC position mismatch");
        }
    }

    private static String computeSignature(List<Entry> entries) {
        final MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException ex) {
            throw new IllegalStateException("SHA-256 is unavailable", ex);
        }
        ByteBuffer buffer = ByteBuffer.allocate(7 * Long.BYTES);
        for (Entry entry : entries) {
            buffer.clear();
            buffer.putLong(entry.pos1);
            buffer.putLong(entry.rc1);
            buffer.putLong(entry.pos2);
            buffer.putLong(entry.rc2);
            buffer.putLong(entry.pos3);
            buffer.putLong(entry.rc3);
            buffer.putLong(Double.doubleToLongBits(
                    entry.storedCorrectionKcal));
            digest.update(buffer.array());
        }
        StringBuilder hex = new StringBuilder(64);
        for (byte value : digest.digest()) {
            hex.append(String.format(Locale.ROOT, "%02x", value & 0xff));
        }
        return hex.toString();
    }

    private static void requireOutput(File output) {
        if (output == null) {
            throw new IllegalArgumentException(
                    "triple eta output must be non-null");
        }
        File parent = output.getAbsoluteFile().getParentFile();
        if (parent == null || !parent.isDirectory()) {
            throw new IllegalArgumentException(
                    "triple eta output parent does not exist: " + output);
        }
        if (output.exists()) {
            throw new IllegalArgumentException(
                    "triple eta refuses to overwrite output: " + output);
        }
    }

    @FunctionalInterface
    private interface WriterAction {
        void write(BufferedWriter writer) throws IOException;
    }

    private static void writeAtomically(Path output, WriterAction action)
            throws IOException {
        Path parent = output.toAbsolutePath().getParent();
        Path temporary = Files.createTempFile(parent,
                output.getFileName().toString() + ".", ".tmp");
        boolean moved = false;
        try {
            try (BufferedWriter writer = Files.newBufferedWriter(
                    temporary, StandardCharsets.UTF_8)) {
                action.write(writer);
            }
            Files.move(temporary, output,
                    StandardCopyOption.ATOMIC_MOVE);
            moved = true;
        } finally {
            if (!moved) Files.deleteIfExists(temporary);
        }
    }

    private static final class TripleKey {
        final int pos1;
        final int rc1;
        final int pos2;
        final int rc2;
        final int pos3;
        final int rc3;

        TripleKey(int pos1, int rc1,
                  int pos2, int rc2,
                  int pos3, int rc3) {
            this.pos1 = pos1;
            this.rc1 = rc1;
            this.pos2 = pos2;
            this.rc2 = rc2;
            this.pos3 = pos3;
            this.rc3 = rc3;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (!(other instanceof TripleKey)) return false;
            TripleKey rhs = (TripleKey) other;
            return pos1 == rhs.pos1 && rc1 == rhs.rc1
                    && pos2 == rhs.pos2 && rc2 == rhs.rc2
                    && pos3 == rhs.pos3 && rc3 == rhs.rc3;
        }

        @Override
        public int hashCode() {
            return Objects.hash(pos1, rc1, pos2, rc2, pos3, rc3);
        }
    }
}
