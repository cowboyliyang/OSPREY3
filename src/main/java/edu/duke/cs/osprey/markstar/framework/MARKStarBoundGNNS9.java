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
 * Strategy 9: Subtree GNN as Search Router (extends S8)
 *
 * Unlike S8 which uses subtree GNN only as a passive bound tightener, S9 uses
 * the ΔF prediction to route internal node processing decisions:
 *
 *   Routing rule (per internal node, after subtree GNN inference):
 *     ΔF - cpQ > 0  → SKIP:      subtree contributes negligibly to Z.
 *                                  Apply GNN tight bound, skip expand entirely.
 *     ΔF + cpQ < 0  → CRITICAL:   CCD more stable than emat estimate.
 *                                  Prioritize for expansion.
 *     otherwise     → UNCERTAIN:  GNN not decisive. Fall back to original MARK* logic.
 *
 *   Depth control:
 *     nfree >= maxFreeForGNN  →  skip GNN, direct expand (GNN unreliable for deep subtrees)
 *
 *   Budget allocation:
 *     Sorts by |ΔF| × Z-gap / Z_emat instead of gap alone, so budget
 *     goes to nodes where GNN correction actually shrinks the bound.
 *
 * Same GNN models as S8, no retraining needed.
 */
public class MARKStarBoundGNNS9 extends MARKStarBoundFastQueues {

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
    private int subtreeBatchSize = 500;
    private double subtreeCpQ = 0.10;
    private int subtreeBudgetMax = 50;

    // --- Routing: depth control ---
    private int maxFreeForGNN = 2;  // only route when nfree <= this

    // --- Leaf GNN pool state ---
    private final List<MARKStarNode> gnnPool = new ArrayList<>();
    private final Set<MARKStarNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private int budgetUsed = 0;

    // --- Subtree GNN pool state ---
    private final List<MARKStarNode> subtreePool = new ArrayList<>();
    private final Set<MARKStarNode> subtreeBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
    private final Map<MARKStarNode, Double> subtreeDeltaFs = new IdentityHashMap<>();
    private int subtreeBudgetUsed = 0;

    // --- Routing classifications ---
    private final Set<MARKStarNode> routedSkip = Collections.newSetFromMap(new IdentityHashMap<>());
    private final Set<MARKStarNode> routedCritical = Collections.newSetFromMap(new IdentityHashMap<>());

    // --- Statistics ---
    private int s9LeafGNNBounded = 0;
    private int s9SubtreeGNNBounded = 0;
    private int s9SubtreeSkipped = 0;
    private int s9SubtreeCritical = 0;
    private int s9SubtreeUncertain = 0;
    private int s9SubtreeDeepSkip = 0;
    private int s9CCDFromGNN = 0;
    private int s9CCDFromOriginal = 0;
    private double s9LeafGNNTimeMs = 0;
    private double s9SubtreeGNNTimeMs = 0;
    private int s9LeafOnnxCalls = 0;
    private int s9SubtreeOnnxCalls = 0;
    private final AtomicLong leafConfsEvaluated = new AtomicLong(0);
    private final AtomicLong subtreeNodesEvaluated = new AtomicLong(0);

    private static int s9RoundCounter = 0;

    public MARKStarBoundGNNS9(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
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
        System.out.println("[Strategy9] Subtree GNN loaded");
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
        System.out.println("[Strategy9] Subtree CP: q=" + String.format("%.6f", q) + ", budget=" + budget);
    }

    public void setMaxFreeForGNN(int maxFree) {
        this.maxFreeForGNN = maxFree;
        System.out.println("[Strategy9] Depth control: maxFreeForGNN=" + maxFree);
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
                + ", skipped=" + s9SubtreeSkipped
                + ", critical=" + s9SubtreeCritical
                + ", uncertain=" + s9SubtreeUncertain
                + ", deepSkip=" + s9SubtreeDeepSkip
                + ", CCD(fromGNN)=" + s9CCDFromGNN
                + ", CCD(fromOrig)=" + s9CCDFromOriginal
                + ", leafBudget=" + budgetUsed + "/" + budgetMax
                + ", subtreeBudget=" + subtreeBudgetUsed + "/" + subtreeBudgetMax
                + ", leafGNNTime=" + String.format("%.1fms", s9LeafGNNTimeMs)
                + ", subtreeGNNTime=" + String.format("%.1fms", s9SubtreeGNNTimeMs)
                + ", leafOnnx=" + s9LeafOnnxCalls
                + ", subtreeOnnx=" + s9SubtreeOnnxCalls);

        // Reset per-pfunc state
        gnnPool.clear();
        gnnBoundedNodes.clear();
        subtreePool.clear();
        subtreeBoundedNodes.clear();
        subtreeDeltaFs.clear();
        routedSkip.clear();
        routedCritical.clear();
        budgetUsed = 0;
        subtreeBudgetUsed = 0;
        s9LeafGNNBounded = 0;
        s9SubtreeGNNBounded = 0;
        s9SubtreeSkipped = 0;
        s9SubtreeCritical = 0;
        s9SubtreeUncertain = 0;
        s9SubtreeDeepSkip = 0;
        s9CCDFromGNN = 0;
        s9CCDFromOriginal = 0;
        s9LeafGNNTimeMs = 0;
        s9SubtreeGNNTimeMs = 0;
        s9LeafOnnxCalls = 0;
        s9SubtreeOnnxCalls = 0;
    }

    // ========================================================================
    // Core: routing-aware bound tightening
    // ========================================================================

    /** Count free positions in a partial assignment. */
    private static int countFree(int[] assignments) {
        int n = 0;
        for (int a : assignments) {
            if (a < 0) n++;
        }
        return n;
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

        // Step 2b: Subtree GNN inference + routing classification
        if (subtreeGNN != null) {
            // Depth filter: only route shallow subtrees
            for (MARKStarNode internal : internalNodes) {
                if (!subtreeBoundedNodes.contains(internal) && !internal.getConfSearchNode().isLeaf()) {
                    int nfree = countFree(internal.getConfSearchNode().assignments);
                    if (nfree <= maxFreeForGNN) {
                        subtreePool.add(internal);
                    } else {
                        s9SubtreeDeepSkip++;
                    }
                }
            }

            if (subtreePool.size() >= subtreeBatchSize) {
                long t = System.nanoTime();
                flushSubtreePool();
                s9SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;

                // --- Route classification ---
                final double kT = 0.5922;
                for (MARKStarNode internal : internalNodes) {
                    if (subtreeBoundedNodes.contains(internal)) continue;
                    Double deltaF = subtreeDeltaFs.get(internal);
                    if (deltaF == null) continue;

                    double pessimistic = deltaF - subtreeCpQ;
                    double optimistic  = deltaF + subtreeCpQ;

                    if (pessimistic > 0) {
                        // Worst-case: CCD less favorable than emat → skip expand
                        routedSkip.add(internal);
                        s9SubtreeSkipped++;
                        // Apply GNN tight bound immediately
                        BigDecimal oldZUpper = internal.getUpperBound();
                        BigDecimal oldZLower = internal.getLowerBound();
                        double scaleUpper = Math.exp(-pessimistic / kT);
                        double scaleLower = Math.exp(-(deltaF + subtreeCpQ) / kT);
                        BigDecimal newZUpper = oldZUpper.multiply(new BigDecimal(scaleUpper));
                        BigDecimal newZLower = oldZUpper.multiply(new BigDecimal(scaleLower));
                        if (newZLower.compareTo(newZUpper) <= 0) {
                            internal.tightenSubtreeBounds(
                                    newZLower.compareTo(oldZLower) > 0 ? newZLower : oldZLower,
                                    newZUpper.compareTo(oldZUpper) < 0 ? newZUpper : oldZUpper);
                        }
                        internal.markUpdated();
                    } else if (optimistic < 0) {
                        // Best-case: CCD more favorable → expand first
                        routedCritical.add(internal);
                        s9SubtreeCritical++;
                    } else {
                        s9SubtreeUncertain++;
                    }
                }

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
                    "[S9 r%d] LEAF eps=%.6f ccd=%d leafPool=%d subtreePool=%d leafBudget=%d/%d subtreeBudget=%d/%d skip=%d crit=%d uncertain=%d deep=%d",
                    s9RoundCounter, epsilonBound, ccdCount,
                    gnnPool.size(), subtreePool.size(),
                    budgetUsed, budgetMax, subtreeBudgetUsed, subtreeBudgetMax,
                    s9SubtreeSkipped, s9SubtreeCritical, s9SubtreeUncertain, s9SubtreeDeepSkip));
            }
        } else {
            // === Internal round with routing (CRITICAL first, SKIP bypassed) ===
            phaseWatch.start();

            internalNodes.sort((a, b) -> {
                boolean aCrit = routedCritical.contains(a);
                boolean bCrit = routedCritical.contains(b);
                if (aCrit != bCrit) return aCrit ? -1 : 1;
                BigDecimal gapA = a.getUpperBound().subtract(a.getLowerBound());
                BigDecimal gapB = b.getUpperBound().subtract(b.getLowerBound());
                return gapB.compareTo(gapA);
            });

            for (MARKStarNode internalNode : internalNodes) {
                if (subtreeBoundedNodes.remove(internalNode)) {
                    subtreeBudgetUsed--;
                }

                // SKIP: GNN guarantees negligible — no expand
                if (routedSkip.contains(internalNode)) {
                    routedSkip.remove(internalNode);
                    continue;
                }

                routedCritical.remove(internalNode);

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

        // Predict free energy corrections ΔF
        double[] deltaFs = subtreeGNN.predictSubtreeResiduals(assignments);
        subtreeNodesEvaluated.addAndGet(assignments.length);
        s9SubtreeOnnxCalls++;

        // Store ΔF for routing classification
        for (int i = 0; i < candidates.size(); i++) {
            subtreeDeltaFs.put(candidates.get(i), deltaFs[i]);
        }

        final double kT = 0.5922;
        int bounded = 0;
        int slotsAvailable = subtreeBudgetMax - subtreeBudgetUsed;

        // Sort by |ΔF| × gap / Z_emat: reward nodes where GNN correction matters
        Integer[] indices = new Integer[candidates.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> {
            MARKStarNode na = candidates.get(a);
            MARKStarNode nb = candidates.get(b);
            double dfA = Math.abs(deltaFs[a]);
            double dfB = Math.abs(deltaFs[b]);
            BigDecimal gapA = na.getUpperBound().subtract(na.getLowerBound());
            BigDecimal gapB = nb.getUpperBound().subtract(nb.getLowerBound());
            double zA = Math.max(na.getUpperBound().doubleValue(), 1e-30);
            double zB = Math.max(nb.getUpperBound().doubleValue(), 1e-30);
            double scoreA = dfA * gapA.doubleValue() / zA;
            double scoreB = dfB * gapB.doubleValue() / zB;
            return Double.compare(scoreB, scoreA);
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
            s9SubtreeGNNBounded++;
            bounded++;
        }

        if (s9RoundCounter % 50 == 0 || bounded > 0) {
            System.out.println(String.format(
                "[S9 r%d] Subtree flush: candidates=%d bounded=%d budget=%d/%d (by |ΔF|×gap/Z)",
                s9RoundCounter, candidates.size(), bounded, subtreeBudgetUsed, subtreeBudgetMax));
        }
    }
}
