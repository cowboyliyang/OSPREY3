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
 * Strategy 8: Decoupled GNN Pool with Subtree GNN for Internal Nodes
 *
 * Extends Strategy 7 by also applying GNN predictions to internal (partially assigned) nodes:
 *   - Leaf nodes: same as S7, uses per-conf GNN (GNNConfEnergyCalculator)
 *   - Internal nodes: uses subtree GNN (GNNSubtreeEnergyCalculator) to predict
 *     free energy correction ΔF = F_CCD - F_emat for partial assignments, tightening
 *     bounds before the expand-or-minimize decision.
 *
 * The subtree GNN takes (fixed_rcs, free_mask) and predicts ΔF such that
 * logZ_CCD = logZ_emat - ΔF/kT (exact identity, not Jensen approximation).
 */
public class MARKStarBoundGNNS8 extends MARKStarBoundFastQueues {

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
    private int subtreeBatchSize = 500;    // fire subtree ONNX batch at this size
    private double subtreeCpQ = 0.10;      // wider CP band for subtree (less accurate)
    private int subtreeBudgetMax = 50;

    // --- Leaf GNN pool state ---
    private final List<MARKStarNode> gnnPool = new ArrayList<>();
    private final Set<MARKStarNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private int budgetUsed = 0;

    // --- Subtree GNN pool state ---
    private final List<MARKStarNode> subtreePool = new ArrayList<>();
    private final Set<MARKStarNode> subtreeBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private int subtreeBudgetUsed = 0;

    // --- Statistics ---
    private int s8LeafGNNBounded = 0;
    private int s8SubtreeGNNBounded = 0;
    private int s8CCDFromGNN = 0;
    private int s8CCDFromOriginal = 0;
    private double s8LeafGNNTimeMs = 0;
    private double s8SubtreeGNNTimeMs = 0;
    private int s8LeafOnnxCalls = 0;
    private int s8SubtreeOnnxCalls = 0;
    private final AtomicLong leafConfsEvaluated = new AtomicLong(0);
    private final AtomicLong subtreeNodesEvaluated = new AtomicLong(0);

    private static int s8RoundCounter = 0;

    public MARKStarBoundGNNS8(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat,
                               ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
    }

    // --- Setters ---

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
        this.gnnBatchCalc = calc;
    }

    public void setSubtreeGNN(GNNSubtreeEnergyCalculator calc) {
        this.subtreeGNN = calc;
        System.out.println("[Strategy8] Subtree GNN loaded");
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = size;
        System.out.println("[Strategy8] leaf gpuBatchSize=" + size);
    }

    public void setSubtreeBatchSize(int size) {
        this.subtreeBatchSize = size;
        System.out.println("[Strategy8] subtree batchSize=" + size);
    }

    public void setCPParams(double alpha, double delta, double q) {
        this.cpAlpha = alpha;
        this.cpDelta = delta;
        this.cpQ = q;
        this.budgetMax = (int) Math.floor(delta / alpha);
        System.out.println("[Strategy8] Leaf CP: alpha=" + alpha + ", delta=" + delta
                + ", q=" + String.format("%.6f", q) + ", budget=" + budgetMax);
    }

    public void setSubtreeCPParams(double q, int budget) {
        this.subtreeCpQ = q;
        this.subtreeBudgetMax = budget;
        System.out.println("[Strategy8] Subtree CP: q=" + String.format("%.6f", q) + ", budget=" + budget);
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
            s8LeafGNNTimeMs += (System.nanoTime() - t) / 1e6;
            updateBound();
        }
        if (subtreeGNN != null && !subtreePool.isEmpty()) {
            long t = System.nanoTime();
            flushSubtreePool();
            s8SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;
            updateBound();
        }

        System.out.println("[Strategy8] Done. eps=" + String.format("%.6f", epsilonBound)
                + ", leafGNNBounded=" + s8LeafGNNBounded
                + ", subtreeGNNBounded=" + s8SubtreeGNNBounded
                + ", CCD(fromGNN)=" + s8CCDFromGNN
                + ", CCD(fromOrig)=" + s8CCDFromOriginal
                + ", leafBudget=" + budgetUsed + "/" + budgetMax
                + ", subtreeBudget=" + subtreeBudgetUsed + "/" + subtreeBudgetMax
                + ", leafGNNTime=" + String.format("%.1fms", s8LeafGNNTimeMs)
                + ", subtreeGNNTime=" + String.format("%.1fms", s8SubtreeGNNTimeMs)
                + ", leafOnnx=" + s8LeafOnnxCalls
                + ", subtreeOnnx=" + s8SubtreeOnnxCalls);

        // Reset per-pfunc state
        gnnPool.clear();
        gnnBoundedNodes.clear();
        subtreePool.clear();
        subtreeBoundedNodes.clear();
        budgetUsed = 0;
        subtreeBudgetUsed = 0;
        s8LeafGNNBounded = 0;
        s8SubtreeGNNBounded = 0;
        s8CCDFromGNN = 0;
        s8CCDFromOriginal = 0;
        s8LeafGNNTimeMs = 0;
        s8SubtreeGNNTimeMs = 0;
        s8LeafOnnxCalls = 0;
        s8SubtreeOnnxCalls = 0;
    }

    // ========================================================================
    // Core: baseline CCD + leaf GNN pool + subtree GNN pool
    // ========================================================================

    @Override
    protected void tightenBoundInPhases() {
        s8RoundCounter++;

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
                s8LeafGNNTimeMs += (System.nanoTime() - t) / 1e6;
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
                s8SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;
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
                    s8CCDFromGNN++;
                } else {
                    s8CCDFromOriginal++;
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

            if (s8RoundCounter % 50 == 0) {
                System.out.println(String.format(
                    "[S8 r%d] LEAF eps=%.6f ccd=%d leafPool=%d subtreePool=%d leafBudget=%d/%d subtreeBudget=%d/%d",
                    s8RoundCounter, epsilonBound, ccdCount,
                    gnnPool.size(), subtreePool.size(),
                    budgetUsed, budgetMax, subtreeBudgetUsed, subtreeBudgetMax));
            }
        } else {
            phaseWatch.start();
            for (MARKStarNode internalNode : internalNodes) {
                // If this internal node was subtree-GNN-bounded, release budget when expanded
                if (subtreeBoundedNodes.remove(internalNode)) {
                    subtreeBudgetUsed--;
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
        s8LeafOnnxCalls++;

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
            s8LeafGNNBounded++;
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

        // Predict free energy corrections ΔF
        double[] deltaFs = subtreeGNN.predictSubtreeResiduals(assignments);
        subtreeNodesEvaluated.addAndGet(assignments.length);
        s8SubtreeOnnxCalls++;

        // GNN predicts ΔF = F_CCD - F_emat = -kT * log(Z_CCD/Z_emat). Identity, not approximation.
        // CP bound: |ΔF_pred - ΔF_true| ≤ subtreeCpQ, so:
        //   Z_lower = Z_emat * exp(-(ΔF + cpQ)/kT)    (ΔF_true ≤ ΔF + cpQ → Z smallest)
        //   Z_upper = Z_emat * exp(-(ΔF - cpQ)/kT)    (ΔF_true ≥ ΔF - cpQ → Z largest)
        // We use the current subtreeUpperBound as Z_emat (emat-based upper bound on Z).
        final double kT = 0.5922;  // kcal/mol at 298K
        int bounded = 0;
        int slotsAvailable = subtreeBudgetMax - subtreeBudgetUsed;

        // Sort by upper-lower gap descending (tighten the loosest bounds first)
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> {
            BigDecimal gapA = candidates.get(a).getUpperBound().subtract(candidates.get(a).getLowerBound());
            BigDecimal gapB = candidates.get(b).getUpperBound().subtract(candidates.get(b).getLowerBound());
            return gapB.compareTo(gapA);  // largest gap first
        });

        for (int rank = 0; rank < indices.length && bounded < slotsAvailable; rank++) {
            int idx = indices[rank];
            MARKStarNode markNode = candidates.get(idx);
            double deltaF = deltaFs[idx];

            BigDecimal oldZUpper = markNode.getUpperBound();
            BigDecimal oldZLower = markNode.getLowerBound();

            // Z_emat: use the current emat-based upper bound as the reference Z
            // (subtreeUpperBound is numConfs * bc(confLowerBound), which is the emat Z upper bound)
            BigDecimal zEmat = oldZUpper;

            // Apply correction with CP bound
            double scaleUpper = Math.exp(-(deltaF - subtreeCpQ) / kT);  // ΔF_true ≥ ΔF-cpQ → Z ≤ this
            double scaleLower = Math.exp(-(deltaF + subtreeCpQ) / kT);  // ΔF_true ≤ ΔF+cpQ → Z ≥ this

            BigDecimal newZUpper = zEmat.multiply(new BigDecimal(scaleUpper));
            BigDecimal newZLower = zEmat.multiply(new BigDecimal(scaleLower));

            // Only tighten: newUpper must be < oldUpper, newLower must be > oldLower
            boolean tightenedUpper = newZUpper.compareTo(oldZUpper) < 0;
            boolean tightenedLower = newZLower.compareTo(oldZLower) > 0;
            if (!tightenedUpper && !tightenedLower) continue;

            // Sanity: new bounds must not cross
            if (newZLower.compareTo(newZUpper) > 0) continue;

            markNode.tightenSubtreeBounds(
                    tightenedLower ? newZLower : oldZLower,
                    tightenedUpper ? newZUpper : oldZUpper);

            subtreeBoundedNodes.add(markNode);
            subtreeBudgetUsed++;
            s8SubtreeGNNBounded++;
            bounded++;
        }

        if (s8RoundCounter % 50 == 0 || bounded > 0) {
            System.out.println(String.format(
                "[S8 r%d] Subtree flush: candidates=%d bounded=%d budget=%d/%d",
                s8RoundCounter, candidates.size(), bounded, subtreeBudgetUsed, subtreeBudgetMax));
        }
    }
}
