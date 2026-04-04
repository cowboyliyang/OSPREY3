package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Strategy 7: Decoupled GNN Pool
 *
 * CCD path is 100% identical to baseline MARK*.
 * GNN runs independently as a side channel:
 *   - Every round, all leaf nodes encountered go into a GNN accumulator pool.
 *   - When pool reaches GPU batch capacity, fire ONE big ONNX inference.
 *   - From results, pick the budgetMax lowest-energy (not yet CCD'd) nodes,
 *     apply CP bounds, update epsilon.
 *   - The CCD decision (leafZ vs internalZ) is NOT affected by GNN pool.
 *
 * This maximizes batch size to amortize ONNX call overhead.
 */
public class MARKStarBoundGNNS7 extends MARKStarBoundFastQueues {

    private GNNConfEnergyCalculator gnnBatchCalc;

    // --- GNN pool parameters ---
    private int gpuBatchSize = 1000;   // fire ONNX when pool reaches this size
    private int budgetMax = 100;       // max simultaneous GNN-bounded nodes
    private double cpQ = 0.06;         // CP quantile (kcal/mol)
    private double cpAlpha = 0.001;
    private double cpDelta = 0.10;

    // --- GNN pool state (reset per pfunc) ---
    private final List<MARKStarNode> gnnPool = new ArrayList<>();
    private final Set<MARKStarNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private final Map<MARKStarNode, Double> gnnPredictions = new IdentityHashMap<>();
    private int budgetUsed = 0;

    // --- Statistics (reset per pfunc) ---
    private int s7GNNBounded = 0;
    private int s7CCDFromGNN = 0;
    private int s7CCDFromOriginal = 0;
    private double s7GNNTimeMs = 0;
    private int s7OnnxCalls = 0;
    private final AtomicLong gnnConfsEvaluated = new AtomicLong(0);

    private static int s7RoundCounter = 0;

    public MARKStarBoundGNNS7(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat,
                               ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
    }

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
        this.gnnBatchCalc = calc;
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = size;
        System.out.println("[Strategy7] gpuBatchSize=" + size);
    }

    public void setCPParams(double alpha, double delta, double q) {
        this.cpAlpha = alpha;
        this.cpDelta = delta;
        this.cpQ = q;
        this.budgetMax = (int) Math.floor(delta / alpha);
        System.out.println("[Strategy7] CP params: alpha=" + alpha + ", delta=" + delta
                + ", q=" + String.format("%.6f", q) + " kcal/mol, budget=" + budgetMax);
    }

    // ========================================================================
    // Lifecycle
    // ========================================================================

    @Override
    public void compute(int maxNumConfs) {
        super.compute(maxNumConfs);

        // Flush any remaining pool at end of pfunc
        if (gnnBatchCalc != null && !gnnPool.isEmpty()) {
            long gnnStart = System.nanoTime();
            flushGNNPool();
            s7GNNTimeMs += (System.nanoTime() - gnnStart) / 1e6;
            updateBound();
        }

        System.out.println("[Strategy7] Done. eps=" + String.format("%.6f", epsilonBound)
                + ", GNNbounded=" + s7GNNBounded
                + ", CCD(fromGNN)=" + s7CCDFromGNN
                + ", CCD(fromOrig)=" + s7CCDFromOriginal
                + ", budgetUsed=" + budgetUsed + "/" + budgetMax
                + ", gnnTime=" + String.format("%.1fms", s7GNNTimeMs)
                + ", onnxCalls=" + s7OnnxCalls
                + ", poolRemain=" + gnnPool.size());

        // Reset per-pfunc state
        gnnPool.clear();
        gnnBoundedNodes.clear();
        gnnPredictions.clear();
        budgetUsed = 0;
        s7GNNBounded = 0;
        s7CCDFromGNN = 0;
        s7CCDFromOriginal = 0;
        s7GNNTimeMs = 0;
        s7OnnxCalls = 0;
    }

    // ========================================================================
    // Core: baseline CCD + independent GNN pool
    // ========================================================================

    @Override
    protected void tightenBoundInPhases() {
        s7RoundCounter++;

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

        // Step 2: Feed leaf nodes into GNN pool (non-invasive)
        if (gnnBatchCalc != null) {
            for (MARKStarNode leaf : leafNodes) {
                if (!gnnBoundedNodes.contains(leaf)) {
                    gnnPool.add(leaf);
                }
            }
            // Also drain leafQueue into pool
            List<MARKStarNode> scannedBack = new ArrayList<>();
            while (!leafQueue.isEmpty() && gnnPool.size() < gpuBatchSize) {
                MARKStarNode extra = leafQueue.poll();
                gnnPool.add(extra);
                scannedBack.add(extra);
            }
            // Put scanned nodes back so baseline CCD can find them later
            leafQueue.addAll(scannedBack);

            // Fire GNN batch when pool is full
            if (gnnPool.size() >= gpuBatchSize) {
                long gnnStart = System.nanoTime();
                flushGNNPool();
                s7GNNTimeMs += (System.nanoTime() - gnnStart) / 1e6;

                updateBound();
                if (epsilonBound <= targetEpsilon) {
                    queue.addAll(leafNodes);
                    queue.addAll(internalNodes);
                    System.out.println(String.format("[S7 r%d] GNN converged eps=%.6f", s7RoundCounter, epsilonBound));
                    return;
                }
            }
        }

        // Step 3: Two-bucket decision — IDENTICAL to baseline MARK*
        boolean isLeafRound = MathTools.isLessThan(internalZ, leafZ);
        if (isLeafRound) {
            // === Leaf (CCD) round ===
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
                    s7CCDFromGNN++;
                } else {
                    s7CCDFromOriginal++;
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

            if (s7RoundCounter % 50 == 0) {
                System.out.println(String.format("[TIMING round %d] LEAF eps=%.6f ccdCount=%d pool=%d budget=%d/%d",
                    s7RoundCounter, epsilonBound, ccdCount, gnnPool.size(), budgetUsed, budgetMax));
            }
        } else {
            // === Internal expansion round ===
            phaseWatch.start();
            for (MARKStarNode internalNode : internalNodes) {
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
    // GNN Pool: one big batch inference + selective CP bounding
    // ========================================================================

    private void flushGNNPool() {
        if (gnnPool.isEmpty()) return;

        // Remove already-bounded nodes
        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : gnnPool) {
            if (!gnnBoundedNodes.contains(node)) {
                candidates.add(node);
            }
        }
        gnnPool.clear();

        if (candidates.isEmpty()) return;

        // ONE big ONNX batch call
        int[][] confs = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            confs[i] = candidates.get(i).getConfSearchNode().assignments;
        }
        double[] gnnEnergies = gnnBatchCalc.calcEnergies(confs);
        gnnConfsEvaluated.addAndGet(confs.length);
        s7OnnxCalls++;

        // Sort by energy ascending (lowest energy = most important for Z)
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> Double.compare(gnnEnergies[a], gnnEnergies[b]));

        // Apply CP bounds to top budgetMax lowest-energy nodes
        int bounded = 0, skipped = 0;
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

            if (newLower > newUpper) {
                skipped++;
                continue;
            }

            leafNode.setBoundsFromConfLowerAndUpper(newLower, newUpper);
            node.gscore = newLower;
            leafNode.markUpdated();

            gnnBoundedNodes.add(leafNode);
            gnnPredictions.put(leafNode, eGNN);
            budgetUsed++;
            s7GNNBounded++;
            bounded++;
        }

        System.out.println(String.format("[S7 r%d] GNN flush: candidates=%d bounded=%d skipped=%d budget=%d/%d",
            s7RoundCounter, candidates.size(), bounded, skipped, budgetUsed, budgetMax));
    }
}
