package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Strategy 9: Conformal logZ-residual subtree bound oracle (S8++).
 *
 * S9 keeps the MARK* search policy intact and uses the subtree GNN only as a
 * calibrated bound tightener. It intentionally does not route, prune, or
 * reprioritize internal nodes. This makes it the experimental class for the
 * "certified neural bound oracle" direction.
 *
 * The subtree model output is interpreted as either:
 *   - deltaF:        ΔF = F_CCD - F_emat in kcal/mol (legacy S8/S9 models)
 *   - logZResidual: r = log Z_CCD - log Z_ref (new paper-facing target)
 *
 * The two are equivalent by r = -ΔF/kT. Internally S9 converts everything to
 * r so bound tightening is expressed directly in partition-function units.
 *
 *   Budget allocation:
 *     Sorts by |r| × Z-gap / Z_ref instead of gap alone, so budget
 *     goes to nodes where GNN correction actually shrinks the bound.
 *
 * Same GNN models as S8, no retraining needed.
 */
public class MARKStarBoundGNNS9 extends MARKStarBoundFastQueues {

    private static final double KT = 0.5922;  // kcal/mol at 298K

    private enum SubtreeOutputMode {
        DELTA_F,
        LOG_Z_RESIDUAL
    }

    // --- Leaf GNN (per-conf, same as S7) ---
    private GNNConfEnergyCalculator gnnBatchCalc;

    // --- Subtree GNN (internal nodes) ---
    private GNNSubtreeEnergyCalculator subtreeGNN;

    // --- GNN pool parameters ---
    private int gpuBatchSize = 1000;
    private int budgetMax = 100;
    private double cpQ = 0.06;         // CP quantile for leaf nodes
    private double cpAlpha = 0.001;
    private double cpDelta = 0.10;

    // --- Subtree GNN parameters ---
    private int subtreeBatchSize = Integer.getInteger("osprey.gnn.subtreeBatchSize", 500);
    private double subtreeCpQ = Double.parseDouble(System.getProperty("osprey.gnn.subtreeCpQ", "0.10"));
    private int subtreeBudgetMax = Integer.getInteger("osprey.gnn.subtreeBudgetMax", 50);
    private final SubtreeOutputMode subtreeOutputMode = parseSubtreeOutputMode();
    private final Map<Integer, Double> subtreeCpQByFree = parseSubtreeCpQByFree();
    private final Map<Integer, Integer> subtreeBudgetByFree = parseSubtreeBudgetByFree();

    // --- Leaf GNN pool state ---
    private final List<MARKStarNode> gnnPool = new ArrayList<>();
    private final Set<MARKStarNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private int budgetUsed = 0;

    // --- Subtree GNN pool state ---
    private final List<MARKStarNode> subtreePool = new ArrayList<>();
    private final Set<MARKStarNode> subtreeBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private final Map<Integer, Integer> subtreeBudgetUsedByFree = new TreeMap<>();
    private int subtreeBudgetUsed = 0;

    // --- Statistics ---
    private int s9LeafGNNBounded = 0;
    private int s9SubtreeGNNBounded = 0;
    private int s9CCDFromGNN = 0;
    private int s9CCDFromOriginal = 0;
    private double s9LeafGNNTimeMs = 0;
    private double s9SubtreeGNNTimeMs = 0;
    private int s9LeafOnnxCalls = 0;
    private int s9SubtreeOnnxCalls = 0;

    @Override
    protected void attachExtraStats(edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Result result) {
        result.setStat("s9LeafGNNBounded", s9LeafGNNBounded);
        result.setStat("s9SubtreeGNNBounded", s9SubtreeGNNBounded);
        result.setStat("s9CCDFromGNN", s9CCDFromGNN);
        result.setStat("s9CCDFromOriginal", s9CCDFromOriginal);
        result.setStat("s9LeafOnnxCalls", s9LeafOnnxCalls);
        result.setStat("s9SubtreeOnnxCalls", s9SubtreeOnnxCalls);
        result.setStat("s9LeafGNNTimeMs", (long) s9LeafGNNTimeMs);
        result.setStat("s9SubtreeGNNTimeMs", (long) s9SubtreeGNNTimeMs);
    }
    private final AtomicLong leafConfsEvaluated = new AtomicLong(0);
    private final AtomicLong subtreeNodesEvaluated = new AtomicLong(0);

    private static int s9RoundCounter = 0;

    public MARKStarBoundGNNS9(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat,
                               ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        if (!subtreeCpQByFree.isEmpty()) {
            System.out.println("[Strategy9] Subtree CP buckets by #free: " + subtreeCpQByFree
                    + " (" + subtreePredictionUnits() + ")");
        }
        if (!subtreeBudgetByFree.isEmpty()) {
            System.out.println("[Strategy9] Subtree active-bound budget by #free: " + subtreeBudgetByFree
                    + ", global cap=" + subtreeBudgetMax);
        }
    }

    // --- Setters ---

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
        this.gnnBatchCalc = calc;
    }

    public void setSubtreeGNN(GNNSubtreeEnergyCalculator calc) {
        this.subtreeGNN = calc;
        System.out.println("[Strategy9] Subtree GNN loaded; outputMode=" + subtreeOutputMode
                + ", internalTarget=logZResidual");
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = size;
        System.out.println("[Strategy9] leaf gpuBatchSize=" + size);
    }

    public void setSubtreeBatchSize(int size) {
        this.subtreeBatchSize = size;
        System.out.println("[Strategy9] subtree batchSize=" + size);
    }

    public void setCPParams(double alpha, double delta, double q) {
        this.cpAlpha = alpha;
        this.cpDelta = delta;
        this.cpQ = q;
        this.budgetMax = (int) Math.floor(delta / alpha);
        System.out.println("[Strategy9] Leaf CP: alpha=" + alpha + ", delta=" + delta
                + ", q=" + String.format("%.6f", q) + ", budget=" + budgetMax);
    }

    public void setSubtreeCPParams(double q, int budget) {
        this.subtreeCpQ = q;
        this.subtreeBudgetMax = budget;
        System.out.println("[Strategy9] Subtree CP: q=" + String.format("%.6f", q)
                + " (" + subtreePredictionUnits() + "), budget=" + budget);
    }

    // ========================================================================
    // Lifecycle
    // ========================================================================

    @Override
    public void compute(int maxNumConfs) {
        super.compute(maxNumConfs);

        // Flush remaining pools
        if (gnnBatchCalc != null && !gnnPool.isEmpty()) {
            long t = System.nanoTime();
            flushLeafGNNPool();
            s9LeafGNNTimeMs += (System.nanoTime() - t) / 1e6;
            updateBound();
        }
        if (subtreeGNN != null && !subtreePool.isEmpty()) {
            long t = System.nanoTime();
            flushSubtreePool();
            s9SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;
            updateBound();
        }

        System.out.println("[Strategy9] Done. eps=" + String.format("%.6f", epsilonBound)
                + ", leafGNNBounded=" + s9LeafGNNBounded
                + ", subtreeBounded=" + s9SubtreeGNNBounded
                + ", CCD(fromGNN)=" + s9CCDFromGNN
                + ", CCD(fromOrig)=" + s9CCDFromOriginal
                + ", leafBudget=" + budgetUsed + "/" + budgetMax
                + ", subtreeBudget=" + subtreeBudgetSummary()
                + ", leafGNNTime=" + String.format("%.1fms", s9LeafGNNTimeMs)
                + ", subtreeGNNTime=" + String.format("%.1fms", s9SubtreeGNNTimeMs)
                + ", leafOnnx=" + s9LeafOnnxCalls
                + ", subtreeOnnx=" + s9SubtreeOnnxCalls);

        // Reset per-pfunc state
        gnnPool.clear();
        gnnBoundedNodes.clear();
        subtreePool.clear();
        subtreeBoundedNodes.clear();
        budgetUsed = 0;
        subtreeBudgetUsed = 0;
        subtreeBudgetUsedByFree.clear();
        s9LeafGNNBounded = 0;
        s9SubtreeGNNBounded = 0;
        s9CCDFromGNN = 0;
        s9CCDFromOriginal = 0;
        s9LeafGNNTimeMs = 0;
        s9SubtreeGNNTimeMs = 0;
        s9LeafOnnxCalls = 0;
        s9SubtreeOnnxCalls = 0;
    }

    // ========================================================================
    // Core: conformal neural bound tightening
    // ========================================================================

    /** Count free positions in a partial assignment. */
    private static int countFree(int[] assignments) {
        int n = 0;
        for (int a : assignments) {
            if (a < 0) n++;
        }
        return n;
    }

    private static SubtreeOutputMode parseSubtreeOutputMode() {
        String raw = System.getProperty("osprey.gnn.subtreeOutput", "deltaF").trim();
        if (raw.equalsIgnoreCase("deltaF") || raw.equalsIgnoreCase("delta_f")) {
            return SubtreeOutputMode.DELTA_F;
        }
        if (raw.equalsIgnoreCase("logZResidual")
                || raw.equalsIgnoreCase("log_z_residual")
                || raw.equalsIgnoreCase("logz")
                || raw.equalsIgnoreCase("logZ")) {
            return SubtreeOutputMode.LOG_Z_RESIDUAL;
        }
        throw new IllegalArgumentException("Unknown osprey.gnn.subtreeOutput=" + raw
                + ". Use deltaF or logZResidual.");
    }

    private static Map<Integer, Double> parseSubtreeCpQByFree() {
        String raw = System.getProperty("osprey.gnn.subtreeCpQByFree", "").trim();
        Map<Integer, Double> out = new TreeMap<>();
        if (raw.isEmpty()) {
            return out;
        }
        for (String entry : raw.split(",")) {
            String[] parts = entry.trim().split(":");
            if (parts.length != 2) {
                throw new IllegalArgumentException("Bad osprey.gnn.subtreeCpQByFree entry: " + entry
                        + ". Expected e.g. 0:0.02,1:0.04,2:0.06");
            }
            out.put(Integer.parseInt(parts[0].trim()), Double.parseDouble(parts[1].trim()));
        }
        return Collections.unmodifiableMap(out);
    }

    private static Map<Integer, Integer> parseSubtreeBudgetByFree() {
        String raw = System.getProperty("osprey.gnn.subtreeBudgetByFree", "").trim();
        Map<Integer, Integer> out = new TreeMap<>();
        if (raw.isEmpty()) {
            return out;
        }
        for (String entry : raw.split(",")) {
            String[] parts = entry.trim().split(":");
            if (parts.length != 2) {
                throw new IllegalArgumentException("Bad osprey.gnn.subtreeBudgetByFree entry: " + entry
                        + ". Expected e.g. 0:10,1:20,2:20");
            }
            int nfree = Integer.parseInt(parts[0].trim());
            int budget = Integer.parseInt(parts[1].trim());
            if (budget < 0) {
                throw new IllegalArgumentException("Negative subtree budget for #free=" + nfree);
            }
            out.put(nfree, budget);
        }
        return Collections.unmodifiableMap(out);
    }

    private String subtreePredictionUnits() {
        return subtreeOutputMode == SubtreeOutputMode.DELTA_F ? "kcal/mol DeltaF" : "logZ residual";
    }

    private double toLogZResidual(double modelOutput) {
        if (subtreeOutputMode == SubtreeOutputMode.DELTA_F) {
            return -modelOutput / KT;
        }
        return modelOutput;
    }

    private double toLogZQ(double q) {
        if (subtreeOutputMode == SubtreeOutputMode.DELTA_F) {
            return q / KT;
        }
        return q;
    }

    private double subtreeLogZQ(int nfree) {
        Double bucketQ = subtreeCpQByFree.get(nfree);
        return toLogZQ(bucketQ != null ? bucketQ : subtreeCpQ);
    }

    private BigDecimal scaledZ(BigDecimal z, double logScale) {
        if (Double.isNaN(logScale)) return null;
        if (logScale > 700.0) return null;  // not a safe/tightening update in double space
        if (logScale < -745.0) return BigDecimal.ZERO;
        return z.multiply(BigDecimal.valueOf(Math.exp(logScale)));
    }

    private int subtreeBucketLimit(int nfree) {
        if (subtreeBudgetByFree.isEmpty()) {
            return subtreeBudgetMax;
        }
        return subtreeBudgetByFree.getOrDefault(nfree, 0);
    }

    private int subtreeBucketUsed(int nfree) {
        return subtreeBudgetUsedByFree.getOrDefault(nfree, 0);
    }

    private boolean hasAnySubtreeBudgetSlot() {
        if (subtreeBudgetUsed >= subtreeBudgetMax) {
            return false;
        }
        if (subtreeBudgetByFree.isEmpty()) {
            return true;
        }
        for (Map.Entry<Integer, Integer> entry : subtreeBudgetByFree.entrySet()) {
            if (subtreeBucketUsed(entry.getKey()) < entry.getValue()) {
                return true;
            }
        }
        return false;
    }

    private boolean hasSubtreeBudgetSlot(int nfree) {
        if (subtreeBudgetUsed >= subtreeBudgetMax) {
            return false;
        }
        if (subtreeBudgetByFree.isEmpty()) {
            return true;
        }
        return subtreeBucketUsed(nfree) < subtreeBucketLimit(nfree);
    }

    private void consumeSubtreeBudgetSlot(int nfree) {
        subtreeBudgetUsed++;
        if (!subtreeBudgetByFree.isEmpty()) {
            subtreeBudgetUsedByFree.put(nfree, subtreeBucketUsed(nfree) + 1);
        }
    }

    private void releaseSubtreeBudgetSlot(int nfree) {
        if (subtreeBudgetUsed > 0) {
            subtreeBudgetUsed--;
        }
        if (!subtreeBudgetByFree.isEmpty()) {
            int used = subtreeBucketUsed(nfree);
            if (used <= 1) {
                subtreeBudgetUsedByFree.remove(nfree);
            } else {
                subtreeBudgetUsedByFree.put(nfree, used - 1);
            }
        }
    }

    private String subtreeBudgetSummary() {
        if (subtreeBudgetByFree.isEmpty()) {
            return subtreeBudgetUsed + "/" + subtreeBudgetMax;
        }
        List<String> parts = new ArrayList<>();
        for (Map.Entry<Integer, Integer> entry : subtreeBudgetByFree.entrySet()) {
            int nfree = entry.getKey();
            parts.add(nfree + ":" + subtreeBucketUsed(nfree) + "/" + entry.getValue());
        }
        return "total=" + subtreeBudgetUsed + "/" + subtreeBudgetMax + ", byFree={" + String.join(",", parts) + "}";
    }

    @Override
    protected void tightenBoundInPhases() {
        s9RoundCounter++;

        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};

        Stopwatch loopWatch = new Stopwatch().start();
        Stopwatch phaseWatch = new Stopwatch();
        double epsilonBeforePhase = epsilonBound;

        // Step 1: populateQueues — identical to baseline
        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        updateBound();
        internalZ = ZSums[0];
        leafZ = ZSums[1];

        // Step 2a: Feed leaf nodes into leaf GNN pool (same as S7)
        if (gnnBatchCalc != null) {
            for (MARKStarNode leaf : leafNodes) {
                if (!gnnBoundedNodes.contains(leaf)) {
                    gnnPool.add(leaf);
                }
            }
            List<MARKStarNode> scannedBack = new ArrayList<>();
            while (!leafQueue.isEmpty() && gnnPool.size() < gpuBatchSize) {
                MARKStarNode extra = leafQueue.poll();
                gnnPool.add(extra);
                scannedBack.add(extra);
            }
            leafQueue.addAll(scannedBack);

            if (gnnPool.size() >= gpuBatchSize) {
                long t = System.nanoTime();
                flushLeafGNNPool();
                s9LeafGNNTimeMs += (System.nanoTime() - t) / 1e6;
                updateBound();
                if (epsilonBound <= targetEpsilon) {
                    queue.addAll(leafNodes);
                    queue.addAll(internalNodes);
                    return;
                }
            }
        }

        // Step 2b: Feed internal nodes into subtree GNN pool
        if (subtreeGNN != null) {
            for (MARKStarNode internal : internalNodes) {
                if (!subtreeBoundedNodes.contains(internal) && !internal.getConfSearchNode().isLeaf()) {
                    subtreePool.add(internal);
                }
            }

            if (subtreePool.size() >= subtreeBatchSize) {
                long t = System.nanoTime();
                flushSubtreePool();
                s9SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;

                updateBound();
                if (epsilonBound <= targetEpsilon) {
                    queue.addAll(leafNodes);
                    queue.addAll(internalNodes);
                    return;
                }
            }
        }

        // Step 3: Two-bucket decision — IDENTICAL to baseline MARK*
        boolean isLeafRound = MathTools.isLessThan(internalZ, leafZ);
        if (isLeafRound) {
            phaseWatch.start();

            leafNodes.sort((a, b) -> {
                BigDecimal gapA = a.getUpperBound().subtract(a.getLowerBound());
                BigDecimal gapB = b.getUpperBound().subtract(b.getLowerBound());
                return gapB.compareTo(gapA);
            });

            int ccdCount = Math.min(parallelism.numThreads, leafNodes.size());
            for (int i = 0; i < ccdCount; i++) {
                MARKStarNode leafNode = leafNodes.get(i);
                processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
                if (gnnBoundedNodes.remove(leafNode)) {
                    budgetUsed--;
                    s9CCDFromGNN++;
                } else {
                    s9CCDFromOriginal++;
                }
            }
            loopTasks.waitForFinish();

            for (int i = ccdCount; i < leafNodes.size(); i++) {
                queue.add(leafNodes.get(i));
            }
            queue.addAll(internalNodes);

            phaseWatch.stop();
            totalLeafTime += phaseWatch.getTimeS();
            totalLeafRounds++;

            if (s9RoundCounter % 50 == 0) {
                System.out.println(String.format(
                    "[S9 r%d] LEAF eps=%.6f ccd=%d leafPool=%d subtreePool=%d leafBudget=%d/%d subtreeBudget=%s",
                    s9RoundCounter, epsilonBound, ccdCount,
                    gnnPool.size(), subtreePool.size(),
                    budgetUsed, budgetMax, subtreeBudgetSummary()));
            }
        } else {
            phaseWatch.start();
            for (MARKStarNode internalNode : internalNodes) {
                // If this internal node was subtree-GNN-bounded, release budget when expanded
                if (subtreeBoundedNodes.remove(internalNode)) {
                    releaseSubtreeBudgetSlot(countFree(internalNode.getConfSearchNode().assignments));
                }

                if (!MathTools.isGreaterThan(internalNode.getLowerBound(), BigDecimal.ONE) &&
                        MathTools.isGreaterThan(
                                MathTools.bigDivide(internalNode.getUpperBound(), rootNode.getUpperBound(),
                                        PartitionFunction.decimalPrecision),
                                new BigDecimal(1 - targetEpsilon))) {
                    loopTasks.submit(() -> {
                        boundLowestBoundConfUnderNode(internalNode, newNodes);
                        return null;
                    }, (ignored) -> {});
                } else {
                    processPartialConfNode(newNodes, internalNode, internalNode.getConfSearchNode());
                }
                internalNode.markUpdated();
            }
            loopTasks.waitForFinish();
            numInternalNodesProcessed += internalNodes.size();

            queue.addAll(leafNodes);

            phaseWatch.stop();
            totalInternalTime += phaseWatch.getTimeS();
            totalInternalRounds++;
        }

        if (epsilonBound <= targetEpsilon) {
            double epsilonReduction = Math.max(0, epsilonBeforePhase - epsilonBound);
            if (isLeafRound) totalLeafEpsilonReduction += epsilonReduction;
            else totalInternalEpsilonReduction += epsilonReduction;
            return;
        }

        loopCleanup(newNodes, loopWatch, Math.max(internalNodes.size(), leafNodes.size()));
        double epsilonReduction = Math.max(0, epsilonBeforePhase - epsilonBound);
        if (isLeafRound) totalLeafEpsilonReduction += epsilonReduction;
        else totalInternalEpsilonReduction += epsilonReduction;
    }

    // ========================================================================
    // Leaf GNN Pool (same as S7)
    // ========================================================================

    private void flushLeafGNNPool() {
        if (gnnPool.isEmpty()) return;

        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : gnnPool) {
            if (!gnnBoundedNodes.contains(node)) {
                candidates.add(node);
            }
        }
        gnnPool.clear();
        if (candidates.isEmpty()) return;

        int[][] confs = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            confs[i] = candidates.get(i).getConfSearchNode().assignments;
        }
        double[] gnnEnergies = gnnBatchCalc.calcEnergies(confs);
        leafConfsEvaluated.addAndGet(confs.length);
        s9LeafOnnxCalls++;

        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> Double.compare(gnnEnergies[a], gnnEnergies[b]));

        int bounded = 0;
        int slotsAvailable = budgetMax - budgetUsed;

        for (int rank = 0; rank < indices.length && bounded < slotsAvailable; rank++) {
            int idx = indices[rank];
            MARKStarNode leafNode = candidates.get(idx);
            MARKStarNode.Node node = leafNode.getConfSearchNode();
            double eGNN = gnnEnergies[idx];

            double oldLower = node.getConfLowerBound();
            double oldUpper = node.getConfUpperBound();
            double newLower = Math.max(eGNN - cpQ, oldLower);
            double newUpper = Math.min(eGNN + cpQ, oldUpper);

            if (newLower > newUpper) continue;

            leafNode.setBoundsFromConfLowerAndUpper(newLower, newUpper);
            node.gscore = newLower;
            leafNode.markUpdated();

            gnnBoundedNodes.add(leafNode);
            budgetUsed++;
            s9LeafGNNBounded++;
            bounded++;
        }
    }

    // ========================================================================
    // Subtree GNN Pool (NEW in S8)
    // ========================================================================

    private void flushSubtreePool() {
        if (subtreePool.isEmpty()) return;

        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : subtreePool) {
            if (!subtreeBoundedNodes.contains(node) && !node.getConfSearchNode().isLeaf()) {
                candidates.add(node);
            }
        }
        subtreePool.clear();
        if (candidates.isEmpty()) return;

        // Build assignments array for subtree GNN
        int[][] assignments = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            assignments[i] = candidates.get(i).getConfSearchNode().assignments.clone();
        }

        // Predict subtree-model outputs and convert to logZ residuals:
        //   legacy ΔF models:        r = -ΔF/kT
        //   logZ-residual models:    r = raw output
        double[] rawOutputs = subtreeGNN.predictSubtreeModelOutputs(assignments);
        double[] logZResiduals = new double[rawOutputs.length];
        for (int i = 0; i < rawOutputs.length; i++) {
            logZResiduals[i] = toLogZResidual(rawOutputs[i]);
        }
        subtreeNodesEvaluated.addAndGet(assignments.length);
        s9SubtreeOnnxCalls++;

        int bounded = 0;

        // Sort by |logZ residual| × gap / Z_ref: reward nodes where GNN correction matters
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> {
            MARKStarNode na = candidates.get(a);
            MARKStarNode nb = candidates.get(b);
            double rA = Math.abs(logZResiduals[a]);
            double rB = Math.abs(logZResiduals[b]);
            BigDecimal gapA = na.getUpperBound().subtract(na.getLowerBound());
            BigDecimal gapB = nb.getUpperBound().subtract(nb.getLowerBound());
            double zA = Math.max(na.getUpperBound().doubleValue(), 1e-30);
            double zB = Math.max(nb.getUpperBound().doubleValue(), 1e-30);
            double scoreA = rA * gapA.doubleValue() / zA;
            double scoreB = rB * gapB.doubleValue() / zB;
            return Double.compare(scoreB, scoreA);
        });

        for (int rank = 0; rank < indices.length && hasAnySubtreeBudgetSlot(); rank++) {
            int idx = indices[rank];
            MARKStarNode markNode = candidates.get(idx);
            double logZResidual = logZResiduals[idx];
            int nfree = countFree(markNode.getConfSearchNode().assignments);
            if (!hasSubtreeBudgetSlot(nfree)) continue;
            double logZQ = subtreeLogZQ(nfree);

            BigDecimal oldZUpper = markNode.getUpperBound();
            BigDecimal oldZLower = markNode.getLowerBound();

            // Z_ref: use the current emat-based upper bound as the reference Z
            // (subtreeUpperBound is numConfs * bc(confLowerBound), which is the emat Z upper bound)
            BigDecimal zRef = oldZUpper;

            // Apply correction with CP bound in logZ space:
            //   r_true in [r_hat - q, r_hat + q]
            //   Z_true in Z_ref * [exp(r_hat - q), exp(r_hat + q)]
            BigDecimal newZUpper = scaledZ(zRef, logZResidual + logZQ);
            BigDecimal newZLower = scaledZ(zRef, logZResidual - logZQ);
            if (newZUpper == null || newZLower == null) continue;

            // Only tighten: newUpper must be < oldUpper or newLower must be > oldLower
            boolean tightenedUpper = newZUpper.compareTo(oldZUpper) < 0;
            boolean tightenedLower = newZLower.compareTo(oldZLower) > 0;
            if (!tightenedUpper && !tightenedLower) continue;

            // Sanity: new bounds must not cross
            if (newZLower.compareTo(newZUpper) > 0) continue;

            markNode.tightenSubtreeBounds(
                    tightenedLower ? newZLower : oldZLower,
                    tightenedUpper ? newZUpper : oldZUpper);

            subtreeBoundedNodes.add(markNode);
            consumeSubtreeBudgetSlot(nfree);
            s9SubtreeGNNBounded++;
            bounded++;
        }

        if (s9RoundCounter % 50 == 0 || bounded > 0) {
            System.out.println(String.format(
                "[S9 r%d] Subtree flush: candidates=%d bounded=%d budget=%s (by |logZres|×gap/Z)",
                s9RoundCounter, candidates.size(), bounded, subtreeBudgetSummary()));
        }
    }
}
