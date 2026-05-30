package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;

/**
 * Strategy 11: leaf GNN replacement plus auditable post hoc correction, with
 * optional subtree-GNN navigation that only reorders internal nodes.
 */
public class MARKStarBoundGNNS11 extends MARKStarBoundFastQueues {

    private static final double KT = 0.5922;

    private GNNSubtreeEnergyCalculator subtreeGNN;
    private AuditLeafLogger auditLogger;

    private boolean subtreeNavigatorEnabled = Boolean.parseBoolean(
            System.getProperty("osprey.gnn.s11.subtreeNavigator", "true"));
    private String landscapeMode = System.getProperty("osprey.gnn.s11.landscapeMode", "mix")
            .trim().toLowerCase();
    private double landscapeMix = clamp01(Double.parseDouble(
            System.getProperty("osprey.gnn.s11.landscapeMix", "0.25")));
    private double originalQuota = clamp01(Double.parseDouble(
            System.getProperty("osprey.gnn.s11.originalQuota",
                    Double.toString(1.0 - landscapeMix))));
    private double landscapeEpsilon = clamp01(Double.parseDouble(
            System.getProperty("osprey.gnn.s11.landscapeEpsilon", "0.25")));
    private double landscapeClamp = Double.parseDouble(
            System.getProperty("osprey.gnn.s11.landscapeClamp", "10.0"));
    private int gpuBatchSize = Integer.getInteger("osprey.gnn.s11.landscapeBatchSize",
            Integer.getInteger("osprey.gnn.gpuBatchSize", 1000));
    private String subtreeOutputMode = System.getProperty("osprey.gnn.s11.subtreeOutput", "deltaF")
            .trim().toLowerCase();

    private final IdentityHashMap<MARKStarNode, Double> landscapeScores = new IdentityHashMap<>();
    private final Random epsilonRng = new Random(Long.getLong("osprey.gnn.s11.landscapeSeed", 20260529L));

    private int s11RoundCounter = 0;
    private int s11LeafGNNReplacements = 0;
    private int s11LeafRejectedAboveUpper = 0;
    private int s11LandscapeReranks = 0;
    private int s11LandscapeCalls = 0;
    private int s11LandscapeEvals = 0;
    private int s11LandscapeSkipped = 0;
    private double s11LandscapeTimeMs = 0;

    public MARKStarBoundGNNS11(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat,
                               ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        System.out.println("[S11] Initialized. subtreeNavigator=" + subtreeNavigatorEnabled
                + ", landscapeMode=" + landscapeMode
                + ", originalQuota=" + originalQuota
                + ", landscapeClamp=" + landscapeClamp
                + ", auditMode=" + System.getProperty("osprey.gnn.s11.auditMode", "full"));
    }

    public void setSubtreeGNN(GNNSubtreeEnergyCalculator calc) {
        this.subtreeGNN = calc;
        System.out.println("[S11] Subtree GNN attached as navigator only"
                + (subtreeNavigatorEnabled ? "" : " (currently disabled)"));
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = Math.max(1, size);
    }

    public void setAuditContext(String designId, String state, int sequenceIndex, String sequence) {
        this.auditLogger = new AuditLeafLogger(designId, state, sequenceIndex, sequence);
    }

    @Override
    public void compute(int maxNumConfs) {
        try {
            super.compute(maxNumConfs);
        } finally {
            closeAuditLogger();
            System.out.println("[S11] Done. eps=" + String.format("%.6f", epsilonBound)
                    + ", leafGNN=" + s11LeafGNNReplacements
                    + ", rejectedAboveUpper=" + s11LeafRejectedAboveUpper
                    + ", navigator=" + (navigatorActive() ? "on" : "off")
                    + ", landscapeCalls=" + s11LandscapeCalls
                    + ", landscapeEvals=" + s11LandscapeEvals
                    + ", reranks=" + s11LandscapeReranks
                    + ", auditRecords=" + getAuditRecordedCount()
                    + ", auditSelected=" + getAuditSelectedCount());
        }
    }

    @Override
    protected void attachExtraStats(PartitionFunction.Result result) {
        result.setStat("s11LeafGNNReplacements", s11LeafGNNReplacements);
        result.setStat("s11LeafRejectedAboveUpper", s11LeafRejectedAboveUpper);
        result.setStat("s11SubtreeNavigatorEnabled", subtreeNavigatorEnabled ? 1 : 0);
        result.setStat("s11LandscapeMode", landscapeModeCode());
        result.setStat("s11LandscapeReranks", s11LandscapeReranks);
        result.setStat("s11LandscapeCalls", s11LandscapeCalls);
        result.setStat("s11LandscapeEvals", s11LandscapeEvals);
        result.setStat("s11LandscapeSkipped", s11LandscapeSkipped);
        result.setStat("s11LandscapeTimeMs", (long) s11LandscapeTimeMs);
        result.setStat("s11AuditMode", auditLogger == null ? 0 : auditLogger.getModeCode());
        result.setStat("s11AuditRecords", getAuditRecordedCount());
        result.setStat("s11AuditSelected", getAuditSelectedCount());
        result.setStat("s11AuditTopSelected", getAuditTopSelectedCount());
        result.setStat("s11AuditTailSampleSelected", getAuditTailSampleSelectedCount());
    }

    @Override
    protected void onGNNLeafReplacement(
            MARKStarNode markNode,
            Node confNode,
            double oldLower,
            double oldUpper,
            double gnnEnergy,
            double newLower,
            double newUpper,
            boolean rejectedAboveUpper) {

        s11LeafGNNReplacements++;
        if (rejectedAboveUpper) s11LeafRejectedAboveUpper++;
        if (auditLogger != null && auditLogger.isEnabled()) {
            auditLogger.record(confNode.assignments, confNode.getLevel(),
                    oldLower, oldUpper, gnnEnergy, newLower, newUpper, rejectedAboveUpper);
        }
    }

    @Override
    protected void tightenBoundInPhases() {
        s11RoundCounter++;

        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] zSums = new BigDecimal[]{internalZ, leafZ};

        Stopwatch loopWatch = new Stopwatch().start();
        Stopwatch phaseWatch = new Stopwatch();
        double epsilonBeforePhase = epsilonBound;

        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, zSums);
        updateBound();
        internalZ = zSums[0];
        leafZ = zSums[1];

        boolean isLeafRound = !leafNodes.isEmpty()
                && (internalNodes.isEmpty() || MathTools.isLessThan(internalZ, leafZ));
        int numNodes = isLeafRound ? leafNodes.size() : internalNodes.size();

        if (isLeafRound) {
            phaseWatch.start();
            processLeafRound(leafNodes, newNodes);
            loopTasks.waitForFinish();
            phaseWatch.stop();
            leafTimeAverage = phaseWatch.getTimeS();
            totalLeafTime += leafTimeAverage;
            totalLeafRounds++;
            queue.addAll(internalNodes);
            if (maxMinimizations < parallelism.numThreads) maxMinimizations++;
        } else {
            phaseWatch.start();
            rerankInternalNodesWithLandscape(internalNodes);
            processInternalRound(internalNodes, newNodes);
            loopTasks.waitForFinish();
            phaseWatch.stop();
            double internalTime = phaseWatch.getTimeS();
            internalTimeAverage = internalTime / Math.max(1, internalNodes.size());
            totalInternalTime += internalTime;
            totalInternalRounds++;
            queue.addAll(leafNodes);
            numInternalNodesProcessed += internalNodes.size();
        }

        updateBound();
        double epsilonReduction = Math.max(0, epsilonBeforePhase - epsilonBound);
        if (epsilonBound <= targetEpsilon) {
            if (isLeafRound) totalLeafEpsilonReduction += epsilonReduction;
            else totalInternalEpsilonReduction += epsilonReduction;
            return;
        }

        loopCleanup(newNodes, loopWatch, numNodes);
        if (isLeafRound) totalLeafEpsilonReduction += epsilonReduction;
        else totalInternalEpsilonReduction += epsilonReduction;

        if (s11RoundCounter % 50 == 0) {
            System.out.println(String.format(
                    "[S11 r%d] eps=%.6f leaf=%d internal=%d leafGNN=%d navigator=%s calls=%d evals=%d audit=%d",
                    s11RoundCounter, epsilonBound, leafNodes.size(), internalNodes.size(),
                    s11LeafGNNReplacements, navigatorActive() ? landscapeMode : "off",
                    s11LandscapeCalls, s11LandscapeEvals, getAuditRecordedCount()));
        }
    }

    private void processLeafRound(List<MARKStarNode> leafNodes, List<MARKStarNode> newNodes) {
        for (MARKStarNode leafNode : leafNodes) {
            processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
            leafNode.markUpdated();
            debugPrint("Processing Node: " + leafNode.getConfSearchNode().toString());
        }
    }

    private void processInternalRound(List<MARKStarNode> internalNodes, List<MARKStarNode> newNodes) {
        for (MARKStarNode internalNode : internalNodes) {
            Node node = internalNode.getConfSearchNode();
            if (!MathTools.isGreaterThan(internalNode.getLowerBound(), BigDecimal.ONE)
                    && MathTools.isGreaterThan(
                    MathTools.bigDivide(internalNode.getUpperBound(), rootNode.getUpperBound(),
                            PartitionFunction.decimalPrecision),
                    new BigDecimal(1 - targetEpsilon))) {
                loopTasks.submit(() -> {
                    boundLowestBoundConfUnderNode(internalNode, newNodes);
                    return null;
                }, ignored -> { });
            } else {
                processPartialConfNode(newNodes, internalNode, node);
            }
            internalNode.markUpdated();
        }
    }

    private void rerankInternalNodesWithLandscape(List<MARKStarNode> internalNodes) {
        if (internalNodes.size() < 2) return;
        if (!navigatorActive()) {
            s11LandscapeSkipped++;
            return;
        }

        scoreLandscape(internalNodes);
        if ("rerank".equals(landscapeMode)) {
            internalNodes.sort(Comparator.comparingDouble(this::landscapeScoreOrDefault).reversed());
            s11LandscapeReranks++;
        } else if ("mix".equals(landscapeMode)) {
            mixInternalOrder(internalNodes);
            s11LandscapeReranks++;
        } else if ("epsilon".equals(landscapeMode)) {
            epsilonInternalOrder(internalNodes);
            s11LandscapeReranks++;
        } else {
            s11LandscapeSkipped++;
        }
    }

    private void scoreLandscape(List<MARKStarNode> nodes) {
        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : nodes) {
            if (landscapeScores.containsKey(node)) continue;
            if (node.getConfSearchNode().isLeaf()) continue;
            candidates.add(node);
        }
        if (candidates.isEmpty()) return;

        for (int start = 0; start < candidates.size(); start += gpuBatchSize) {
            int end = Math.min(candidates.size(), start + gpuBatchSize);
            int[][] assignments = new int[end - start][];
            for (int i = start; i < end; i++) {
                assignments[i - start] = candidates.get(i).getConfSearchNode().assignments.clone();
            }

            long t = System.nanoTime();
            double[] rawOutputs = subtreeGNN.predictSubtreeModelOutputs(assignments);
            s11LandscapeTimeMs += (System.nanoTime() - t) / 1e6;
            s11LandscapeCalls++;
            s11LandscapeEvals += rawOutputs.length;

            for (int i = 0; i < rawOutputs.length; i++) {
                MARKStarNode node = candidates.get(start + i);
                double logZResidual = subtreeRawToLogZResidual(rawOutputs[i]);
                double logUpper = logBigDecimal(node.getUpperBound());
                double score = logUpper + clamp(logZResidual, -landscapeClamp, landscapeClamp);
                if (!Double.isFinite(score)) score = logUpper;
                landscapeScores.put(node, score);
            }
        }
    }

    private void mixInternalOrder(List<MARKStarNode> internalNodes) {
        int n = internalNodes.size();
        int nOrig = Math.max(0, Math.min(n, (int)Math.ceil(originalQuota * n)));
        int nGnn = n - nOrig;
        if (nGnn <= 0) return;

        List<MARKStarNode> original = new ArrayList<>(internalNodes);
        List<MARKStarNode> byGnn = new ArrayList<>(internalNodes);
        byGnn.sort(Comparator.comparingDouble(this::landscapeScoreOrDefault).reversed());

        LinkedHashSet<MARKStarNode> mixed = new LinkedHashSet<>();
        for (int i = 0; i < nOrig; i++) mixed.add(original.get(i));
        for (int i = 0; i < nGnn; i++) mixed.add(byGnn.get(i));
        mixed.addAll(original);

        internalNodes.clear();
        for (MARKStarNode node : mixed) {
            internalNodes.add(node);
            if (internalNodes.size() >= n) break;
        }
    }

    private void epsilonInternalOrder(List<MARKStarNode> internalNodes) {
        if (epsilonRng.nextDouble() >= landscapeEpsilon) return;
        MARKStarNode best = Collections.max(internalNodes,
                Comparator.comparingDouble(this::landscapeScoreOrDefault));
        internalNodes.remove(best);
        internalNodes.add(0, best);
    }

    private boolean navigatorActive() {
        return subtreeNavigatorEnabled
                && subtreeGNN != null
                && !"off".equals(landscapeMode)
                && !"none".equals(landscapeMode)
                && !"false".equals(landscapeMode);
    }

    private double landscapeScoreOrDefault(MARKStarNode node) {
        Double score = landscapeScores.get(node);
        if (score != null) return score;
        return logBigDecimal(node.getUpperBound());
    }

    private double subtreeRawToLogZResidual(double raw) {
        if ("logz".equals(subtreeOutputMode)
                || "logz_residual".equals(subtreeOutputMode)
                || "residual".equals(subtreeOutputMode)) {
            return raw;
        }
        return -raw / KT;
    }

    private void closeAuditLogger() {
        if (auditLogger == null) return;
        auditLogger.close();
    }

    private int getAuditRecordedCount() {
        return auditLogger == null ? 0 : auditLogger.getRecordedCount();
    }

    private int getAuditSelectedCount() {
        return auditLogger == null ? 0 : auditLogger.getSelectedCount();
    }

    private int getAuditTopSelectedCount() {
        return auditLogger == null ? 0 : auditLogger.getTopSelectedCount();
    }

    private int getAuditTailSampleSelectedCount() {
        return auditLogger == null ? 0 : auditLogger.getTailSampleSelectedCount();
    }

    private int landscapeModeCode() {
        if ("off".equals(landscapeMode) || "none".equals(landscapeMode) || "false".equals(landscapeMode)) return 0;
        if ("rerank".equals(landscapeMode)) return 1;
        if ("mix".equals(landscapeMode)) return 2;
        if ("epsilon".equals(landscapeMode)) return 3;
        return -1;
    }

    private static double logBigDecimal(BigDecimal z) {
        if (z == null || z.signum() <= 0) return Double.NEGATIVE_INFINITY;
        int adjustedExponent = z.precision() - z.scale() - 1;
        BigDecimal mantissa = z.movePointLeft(adjustedExponent);
        return Math.log(mantissa.doubleValue()) + adjustedExponent * Math.log(10.0);
    }

    private static double clamp(double value, double min, double max) {
        return Math.max(min, Math.min(max, value));
    }

    private static double clamp01(double value) {
        if (Double.isNaN(value)) return 0.0;
        return clamp(value, 0.0, 1.0);
    }
}
