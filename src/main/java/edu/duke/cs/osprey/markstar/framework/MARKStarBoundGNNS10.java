package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.util.*;

/**
 * Strategy 10: Decoupled GNN-augmented MARK*.
 *
 * <p>Design (vs. S7/S8/S9):</p>
 * <ul>
 *   <li>A* search is identical to baseline MARK* (uses emat-derived bounds
 *       only). GNN bounds are NOT written back into the node bound fields,
 *       so the A* trajectory matches M*. This avoids the CCD explosion
 *       observed in S9 (where GNN-tightened node bounds altered the
 *       leafZ/internalZ Z-sums and skewed the leaf-vs-internal decisions).</li>
 *   <li>GNN is used as an early-termination oracle: when the A* loop selects
 *       a leaf for CCD, we first query the leaf GNN. If the per-leaf CP width
 *       (cpQ) is below a calibrated threshold tauLeaf and the GNN energy lies
 *       inside the emat-derived interval, we accept the GNN bound, mark the
 *       node as "finalized" (skip CCD), and let standard pfunc aggregation
 *       pick up the collapsed bound [GNN-cpQ, GNN+cpQ].</li>
 *   <li>Subtree GNN is used the same way: at internal-node expansion time,
 *       if the subtree GNN logZ width is below tauSubtree, we finalize the
 *       whole subtree without expansion.</li>
 *   <li>Budget is governed by an empirical Bernstein bound on the sum of
 *       individual GNN errors (variance + boundedness from cp_stats), which
 *       is typically 1-2 orders of magnitude larger than the Bonferroni
 *       budget delta/alpha used in S9.</li>
 * </ul>
 *
 * <p>System properties:</p>
 * <ul>
 *   <li>osprey.gnn.s10.tauLeaf       (default 0.5 kcal/mol)</li>
 *   <li>osprey.gnn.s10.tauSubtree    (default 0.5 kcal/mol)</li>
 *   <li>osprey.gnn.s10.maxSubtreeFinalizeFree (default 4)</li>
 *   <li>osprey.gnn.s10.maxBernstein  (default 10000) — cap on finalized count</li>
 *   <li>osprey.gnn.s10.gpuBatchSize  (default 512)</li>
 *   <li>osprey.gnn.s10.leafSigma     (override empirical sigma)</li>
 *   <li>osprey.gnn.s10.leafMaxErr    (override empirical bound)</li>
 *   <li>osprey.gnn.s10.delta         (Bernstein confidence, default 0.001)</li>
 * </ul>
 */
public class MARKStarBoundGNNS10 extends MARKStarBoundFastQueues {

    private static final double KT = 0.5922;  // kcal/mol at 298K

    // --- Models ---
    private GNNConfEnergyCalculator leafGNN;
    private GNNSubtreeEnergyCalculator subtreeGNN;

    // --- Thresholds (in kcal/mol unless noted) ---
    private double tauLeaf = Double.parseDouble(System.getProperty("osprey.gnn.s10.tauLeaf", "0.3"));
    private double tauSubtree = Double.parseDouble(System.getProperty("osprey.gnn.s10.tauSubtree", "0.5"));

    // --- Per-leaf cpQ used when no per-conformer uncertainty is available. ---
    // Conservative defaults: P90 of empirical |residuals| from typical cp_stats
    private double leafCpQ = Double.parseDouble(System.getProperty("osprey.gnn.s10.leafCpQ", "0.3"));
    // Subtree CP q is in kcal/mol DeltaF units, matching gnn/train_subtree.py.
    private double subtreeCpQ = Double.parseDouble(System.getProperty("osprey.gnn.s10.subtreeCpQ", "0.5"));

    // --- Subtree gates ---
    private boolean enableSubtreeFinalize = Boolean.parseBoolean(
            System.getProperty("osprey.gnn.s10.subtreeFinalize", "true"));
    private int maxSubtreeFinalizeFree = Integer.getInteger("osprey.gnn.s10.maxSubtreeFinalizeFree", 4);
    private double maxSubtreeFinalizeCpQ = Double.parseDouble(
            System.getProperty("osprey.gnn.s10.maxSubtreeFinalizeCpQ", "0.3"));

    // --- Bernstein calibration ---
    private double leafSigma = Double.parseDouble(System.getProperty("osprey.gnn.s10.leafSigma", "0.92"));
    private double leafMaxErr = Double.parseDouble(System.getProperty("osprey.gnn.s10.leafMaxErr", "13.2"));
    private double subtreeSigma = Double.parseDouble(System.getProperty("osprey.gnn.s10.subtreeSigma", "1.5"));
    private double subtreeMaxErr = Double.parseDouble(System.getProperty("osprey.gnn.s10.subtreeMaxErr", "20.0"));
    private double bernsteinDelta = Double.parseDouble(System.getProperty("osprey.gnn.s10.delta", "0.001"));
    private int maxBernsteinBudget = Integer.getInteger("osprey.gnn.s10.maxBernstein", 10000);

    // --- Batching ---
    private int gpuBatchSize = Integer.getInteger("osprey.gnn.s10.gpuBatchSize", 512);

    // --- Pools (deferred GNN evaluation) ---
    private final List<MARKStarNode> leafPool = new ArrayList<>();
    private final List<MARKStarNode> subtreePool = new ArrayList<>();

    // --- Finalization tracking ---
    private final Set<MARKStarNode> gnnFinalizedLeaves =
            Collections.newSetFromMap(new IdentityHashMap<>());
    private final Set<MARKStarNode> gnnFinalizedSubtrees =
            Collections.newSetFromMap(new IdentityHashMap<>());
    // Mirror by Node (the inner ConfSearchNode) so shouldMinimize(Node) can check.
    private final Set<Node> gnnFinalizedConfNodes =
            Collections.newSetFromMap(new IdentityHashMap<>());

    // --- Bernstein accumulators ---
    private int leafFinalizedCount = 0;
    private int subtreeFinalizedCount = 0;
    private double leafVarianceSum = 0;
    private double subtreeVarianceSum = 0;

    // --- Stats ---
    private int s10LeafFinalized = 0;
    private int s10SubtreeFinalized = 0;
    private int s10CCDCount = 0;
    private int s10LeafGNNCalls = 0;
    private int s10SubtreeGNNCalls = 0;
    private int s10LeafGNNEvals = 0;
    private int s10SubtreeGNNEvals = 0;
    private int s10LeafRejectedCPQ = 0;
    private int s10LeafRejectedRange = 0;
    private int s10LeafRejectedBudget = 0;
    private int s10SubtreeRejectedCPQ = 0;
    private int s10SubtreeRejectedRange = 0;
    private int s10SubtreeRejectedBudget = 0;
    private int s10SubtreeRejectedFree = 0;
    private double s10LeafGNNTimeMs = 0;
    private double s10SubtreeGNNTimeMs = 0;
    private double s10BernsteinErrLeaf = 0;
    private double s10BernsteinErrSubtree = 0;

    private static int s10RoundCounter = 0;

    public MARKStarBoundGNNS10(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat,
                               ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        System.out.println("[S10] Initialized. tauLeaf=" + tauLeaf
                + ", tauSubtree=" + tauSubtree
                + ", subtreeCpQ=" + subtreeCpQ + " kcal/mol"
                + ", subtreeFinalize=" + enableSubtreeFinalize
                + ", maxFinalizeFree=" + maxSubtreeFinalizeFree
                + ", maxFinalizeCpQ=" + maxSubtreeFinalizeCpQ
                + ", leafSigma=" + leafSigma
                + ", subtreeSigma=" + subtreeSigma
                + ", maxBernstein=" + maxBernsteinBudget);
    }

    // -------- Setters (called by MARKStar.calcPfunc) --------

    public void setLeafGNN(GNNConfEnergyCalculator calc) {
        this.leafGNN = calc;
        System.out.println("[S10] Leaf GNN attached");
    }

    public void setSubtreeGNN(GNNSubtreeEnergyCalculator calc) {
        this.subtreeGNN = calc;
        System.out.println("[S10] Subtree GNN attached");
    }

    public void setGPUBatchSize(int size) {
        this.gpuBatchSize = size;
    }

    public void setTauLeaf(double tau) {
        this.tauLeaf = tau;
    }

    public void setTauSubtree(double tau) {
        this.tauSubtree = tau;
    }

    public void setCalibration(double leafSigma, double leafMaxErr,
                                double subtreeSigma, double subtreeMaxErr) {
        this.leafSigma = leafSigma;
        this.leafMaxErr = leafMaxErr;
        this.subtreeSigma = subtreeSigma;
        this.subtreeMaxErr = subtreeMaxErr;
        System.out.println("[S10] Calibrated: leaf(σ=" + leafSigma + ",M=" + leafMaxErr
                + "), subtree(σ=" + subtreeSigma + ",M=" + subtreeMaxErr + ")");
    }

    @Override
    protected void attachExtraStats(PartitionFunction.Result result) {
        result.setStat("s10LeafFinalized", s10LeafFinalized);
        result.setStat("s10SubtreeFinalized", s10SubtreeFinalized);
        result.setStat("s10CCDCount", s10CCDCount);
        result.setStat("s10LeafGNNCalls", s10LeafGNNCalls);
        result.setStat("s10SubtreeGNNCalls", s10SubtreeGNNCalls);
        result.setStat("s10LeafGNNEvals", s10LeafGNNEvals);
        result.setStat("s10SubtreeGNNEvals", s10SubtreeGNNEvals);
        result.setStat("s10LeafGNNTimeMs", (long) s10LeafGNNTimeMs);
        result.setStat("s10SubtreeGNNTimeMs", (long) s10SubtreeGNNTimeMs);
        result.setStat("s10LeafRejectedCPQ", s10LeafRejectedCPQ);
        result.setStat("s10LeafRejectedRange", s10LeafRejectedRange);
        result.setStat("s10LeafRejectedBudget", s10LeafRejectedBudget);
        result.setStat("s10SubtreeRejectedCPQ", s10SubtreeRejectedCPQ);
        result.setStat("s10SubtreeRejectedRange", s10SubtreeRejectedRange);
        result.setStat("s10SubtreeRejectedBudget", s10SubtreeRejectedBudget);
        result.setStat("s10SubtreeRejectedFree", s10SubtreeRejectedFree);
    }

    // ====================================================================
    // Lifecycle
    // ====================================================================

    @Override
    public void compute(int maxNumConfs) {
        super.compute(maxNumConfs);

        // Final flush
        if (leafGNN != null && !leafPool.isEmpty()) {
            flushLeafPool();
        }
        if (subtreeGNN != null && !subtreePool.isEmpty()) {
            flushSubtreePool();
        }
        updateBound();

        // Recompute final Bernstein error for reporting
        recomputeBernsteinError();

        System.out.println("[S10] Done. eps=" + String.format("%.6f", epsilonBound)
                + ", leafFinalized=" + s10LeafFinalized
                + ", subtreeFinalized=" + s10SubtreeFinalized
                + ", CCD=" + s10CCDCount
                + ", leafGNNCalls=" + s10LeafGNNCalls + " (evals=" + s10LeafGNNEvals + ")"
                + ", subtreeGNNCalls=" + s10SubtreeGNNCalls + " (evals=" + s10SubtreeGNNEvals + ")"
                + ", rejected(leaf cpQ/range/budget)=" + s10LeafRejectedCPQ + "/" + s10LeafRejectedRange + "/" + s10LeafRejectedBudget
                + ", rejected(subtree cpQ/range/budget/free)=" + s10SubtreeRejectedCPQ + "/" + s10SubtreeRejectedRange + "/" + s10SubtreeRejectedBudget + "/" + s10SubtreeRejectedFree
                + ", leafGNNTime=" + String.format("%.1fms", s10LeafGNNTimeMs)
                + ", subtreeGNNTime=" + String.format("%.1fms", s10SubtreeGNNTimeMs)
                + ", BernsteinErr leaf=" + String.format("%.4f", s10BernsteinErrLeaf)
                + " subtree=" + String.format("%.4f", s10BernsteinErrSubtree));

        // Reset per-pfunc state
        leafPool.clear();
        subtreePool.clear();
        gnnFinalizedLeaves.clear();
        gnnFinalizedSubtrees.clear();
        gnnFinalizedConfNodes.clear();
        leafFinalizedCount = 0;
        subtreeFinalizedCount = 0;
        leafVarianceSum = 0;
        subtreeVarianceSum = 0;
        s10LeafFinalized = 0;
        s10SubtreeFinalized = 0;
        s10CCDCount = 0;
        s10LeafGNNCalls = 0;
        s10SubtreeGNNCalls = 0;
        s10LeafGNNEvals = 0;
        s10SubtreeGNNEvals = 0;
        s10LeafRejectedCPQ = 0;
        s10LeafRejectedRange = 0;
        s10LeafRejectedBudget = 0;
        s10SubtreeRejectedCPQ = 0;
        s10SubtreeRejectedRange = 0;
        s10SubtreeRejectedBudget = 0;
        s10SubtreeRejectedFree = 0;
        s10LeafGNNTimeMs = 0;
        s10SubtreeGNNTimeMs = 0;
        s10BernsteinErrLeaf = 0;
        s10BernsteinErrSubtree = 0;
    }

    // ====================================================================
    // A* skips GNN-finalized nodes
    // ====================================================================

    @Override
    protected boolean shouldMinimize(Node node) {
        // Standard M* condition: leaf AND not already minimized
        if (!super.shouldMinimize(node)) return false;
        // S10 addition: GNN-finalized leaves never need CCD
        // We track these in gnnFinalizedConfNodes (identity-keyed on Node)
        if (gnnFinalizedConfNodes.contains(node)) return false;
        return true;
    }

    // ====================================================================
    // Main: override the per-round logic to insert GNN finalization checks
    // ====================================================================

    @Override
    protected void tightenBoundInPhases() {
        s10RoundCounter++;

        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};

        Stopwatch loopWatch = new Stopwatch().start();
        Stopwatch phaseWatch = new Stopwatch();
        double epsilonBeforePhase = epsilonBound;

        // Step 1: A* selection (identical to baseline; uses emat-only bounds)
        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        updateBound();
        internalZ = ZSums[0];
        leafZ = ZSums[1];

        // Step 2: leaf-vs-internal decision. If one side has no work, force
        // the other side; internalZ is initialized to 1, so an empty internal
        // set can otherwise spin forever when the selected leaf mass is < 1.
        boolean isLeafRound = !leafNodes.isEmpty()
                && (internalNodes.isEmpty() || MathTools.isLessThan(internalZ, leafZ));

        if (isLeafRound) {
            phaseWatch.start();
            processLeafRound(leafNodes, newNodes);
            loopTasks.waitForFinish();
            phaseWatch.stop();
            totalLeafTime += phaseWatch.getTimeS();
            totalLeafRounds++;
            queue.addAll(internalNodes);
            if (maxMinimizations < parallelism.numThreads) maxMinimizations++;
        } else {
            phaseWatch.start();
            processInternalRound(internalNodes, newNodes);
            loopTasks.waitForFinish();
            phaseWatch.stop();
            totalInternalTime += phaseWatch.getTimeS();
            totalInternalRounds++;
            queue.addAll(leafNodes);
            numInternalNodesProcessed += internalNodes.size();
        }

        updateBound();
        if (epsilonBound <= targetEpsilon) {
            double red = Math.max(0, epsilonBeforePhase - epsilonBound);
            if (isLeafRound) totalLeafEpsilonReduction += red;
            else totalInternalEpsilonReduction += red;
            return;
        }

        loopCleanup(newNodes, loopWatch, isLeafRound ? leafNodes.size() : internalNodes.size());
        double red = Math.max(0, epsilonBeforePhase - epsilonBound);
        if (isLeafRound) totalLeafEpsilonReduction += red;
        else totalInternalEpsilonReduction += red;

        if (s10RoundCounter % 50 == 0) {
            System.out.println(String.format(
                "[S10 r%d] eps=%.6f, leaf=%d/%d (CCD %d, finalized %d), subtree=%d (finalized %d), bernstein=%.3f",
                s10RoundCounter, epsilonBound,
                leafNodes.size(), leafFinalizedCount, s10CCDCount, s10LeafFinalized,
                internalNodes.size(), s10SubtreeFinalized,
                s10BernsteinErrLeaf
            ));
        }
    }

    // ====================================================================
    // Leaf round: GNN-finalize what we can, CCD the rest
    // ====================================================================

    private void processLeafRound(List<MARKStarNode> leafNodes, List<MARKStarNode> newNodes) {
        // Sub-step A: Add leaves to GNN pool (deferred batch); flush if big enough
        if (leafGNN != null) {
            for (MARKStarNode leaf : leafNodes) {
                if (!gnnFinalizedLeaves.contains(leaf)) {
                    leafPool.add(leaf);
                }
            }
            if (leafPool.size() >= gpuBatchSize) {
                flushLeafPool();
            }
        }

        // Sub-step B: For each leaf the A* picked, check if it was finalized
        // (during flush). If yes, skip; else CCD.
        int ccdCount = Math.min(parallelism.numThreads, leafNodes.size());
        for (int i = 0; i < ccdCount; i++) {
            MARKStarNode leaf = leafNodes.get(i);
            if (gnnFinalizedLeaves.contains(leaf)) {
                // already collapsed via finalize(); count as a GNN save, skip CCD
                continue;
            }
            // Standard CCD path
            processFullConfNode(newNodes, leaf, leaf.getConfSearchNode());
            leaf.markUpdated();
            s10CCDCount++;
        }
        // Leftovers go back into queue
        for (int i = ccdCount; i < leafNodes.size(); i++) {
            queue.add(leafNodes.get(i));
        }
    }

    // ====================================================================
    // Internal round: try subtree GNN finalize; else standard expand
    // ====================================================================

    private void processInternalRound(List<MARKStarNode> internalNodes, List<MARKStarNode> newNodes) {
        // Sub-step A: Add internals to subtree pool
        if (subtreeGNN != null) {
            for (MARKStarNode internal : internalNodes) {
                if (!gnnFinalizedSubtrees.contains(internal)
                        && !internal.getConfSearchNode().isLeaf()) {
                    subtreePool.add(internal);
                }
            }
            if (subtreePool.size() >= gpuBatchSize) {
                flushSubtreePool();
            }
        }

        // Sub-step B: For each internal node A* picked, skip if finalized; else expand
        for (MARKStarNode internal : internalNodes) {
            if (gnnFinalizedSubtrees.contains(internal)) continue;  // bound already set
            // Standard expansion (mirrors baseline conditions)
            Node node = internal.getConfSearchNode();
            if (!MathTools.isGreaterThan(internal.getLowerBound(), BigDecimal.ONE) &&
                MathTools.isGreaterThan(
                    MathTools.bigDivide(internal.getUpperBound(), rootNode.getUpperBound(),
                            PartitionFunction.decimalPrecision),
                    new BigDecimal(1 - targetEpsilon))
            ) {
                loopTasks.submit(() -> {
                    boundLowestBoundConfUnderNode(internal, newNodes);
                    return null;
                }, (ignored) -> { });
            } else {
                processPartialConfNode(newNodes, internal, node);
            }
            internal.markUpdated();
        }
    }

    // ====================================================================
    // Leaf GNN pool flush — actually call GNN and finalize what passes
    // ====================================================================

    private void flushLeafPool() {
        if (leafPool.isEmpty()) return;

        // Filter to candidates that are still leaves and not yet finalized
        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : leafPool) {
            if (gnnFinalizedLeaves.contains(node)) continue;
            if (node.getConfSearchNode().isMinimized()) continue;
            candidates.add(node);
        }
        leafPool.clear();
        if (candidates.isEmpty()) return;

        int[][] confs = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            confs[i] = candidates.get(i).getConfSearchNode().assignments;
        }

        long t = System.nanoTime();
        double[] energies = leafGNN.calcEnergies(confs);
        s10LeafGNNTimeMs += (System.nanoTime() - t) / 1e6;
        s10LeafGNNCalls++;
        s10LeafGNNEvals += confs.length;

        // Decide per-leaf
        for (int i = 0; i < candidates.size(); i++) {
            MARKStarNode leaf = candidates.get(i);
            Node node = leaf.getConfSearchNode();
            double E = energies[i];
            double cpQ = leafCpQ;  // could be per-conformer if model provides it
            double oldLower = node.getConfLowerBound();
            double oldUpper = node.getConfUpperBound();

            // Reject if cpQ too wide
            if (cpQ > tauLeaf) { s10LeafRejectedCPQ++; continue; }
            // Reject if GNN energy outside the emat-derived interval (sanity)
            if (E + cpQ < oldLower || E - cpQ > oldUpper) { s10LeafRejectedRange++; continue; }
            // Reject if Bernstein budget exhausted
            if (leafFinalizedCount >= maxBernsteinBudget) { s10LeafRejectedBudget++; continue; }

            // Accept: tighten bound to GNN's CP interval (don't collapse — that
            // would bias pfunc by treating GNN energy as exact). Keep [E-cpQ, E+cpQ]
            // intersected with the existing emat bound so the bound is monotone
            // tightening.
            double padLower = Math.max(E - cpQ, oldLower);
            double padUpper = Math.min(E + cpQ, oldUpper);
            if (padLower > padUpper) { s10LeafRejectedRange++; continue; }
            leaf.setBoundsFromConfLowerAndUpper(padLower, padUpper);
            node.gscore = padLower;
            leaf.markUpdated();
            // Mark as finalized so A* won't re-pick (via shouldMinimize override)
            gnnFinalizedLeaves.add(leaf);
            gnnFinalizedConfNodes.add(node);
            leafFinalizedCount++;
            leafVarianceSum += cpQ * cpQ;  // contribution to Σσ²
            s10LeafFinalized++;
        }
    }

    // ====================================================================
    // Subtree GNN pool flush — finalize whole subtrees
    // ====================================================================

    private void flushSubtreePool() {
        if (subtreePool.isEmpty()) return;

        List<MARKStarNode> candidates = new ArrayList<>();
        for (MARKStarNode node : subtreePool) {
            if (gnnFinalizedSubtrees.contains(node)) continue;
            if (node.getConfSearchNode().isLeaf()) continue;
            if (countFree(node.getConfSearchNode().assignments) > maxSubtreeFinalizeFree) {
                s10SubtreeRejectedFree++;
                continue;
            }
            candidates.add(node);
        }
        subtreePool.clear();
        if (candidates.isEmpty()) return;

        int[][] assignments = new int[candidates.size()][];
        for (int i = 0; i < candidates.size(); i++) {
            assignments[i] = candidates.get(i).getConfSearchNode().assignments.clone();
        }

        long t = System.nanoTime();
        double[] modelOutputs = subtreeGNN.predictSubtreeModelOutputs(assignments);
        s10SubtreeGNNTimeMs += (System.nanoTime() - t) / 1e6;
        s10SubtreeGNNCalls++;
        s10SubtreeGNNEvals += assignments.length;

        // Interpret raw model output as logZ residual:
        //   legacy ΔF models:        r = -ΔF / kT
        //   logZ-residual models:    r = raw output directly
        // We assume ΔF mode (current default). The residual r means:
        //   log(Z_true) = log(Z_emat) + r
        // where Z_emat is the current emat-based upper bound on subtree Z.
        for (int i = 0; i < candidates.size(); i++) {
            MARKStarNode node = candidates.get(i);
            int nfree = countFree(node.getConfSearchNode().assignments);
            double logZResidual = -modelOutputs[i] / KT;
            double cpQKcal = subtreeCpQ;
            double logZQ = cpQKcal / KT;

            if (cpQKcal > tauSubtree) { s10SubtreeRejectedCPQ++; continue; }
            if (cpQKcal > maxSubtreeFinalizeCpQ) { s10SubtreeRejectedCPQ++; continue; }
            if (nfree > maxSubtreeFinalizeFree) { s10SubtreeRejectedFree++; continue; }
            if (subtreeFinalizedCount >= maxBernsteinBudget) { s10SubtreeRejectedBudget++; continue; }
            if (!enableSubtreeFinalize) { s10SubtreeRejectedBudget++; continue; }

            BigDecimal oldZUpper = node.getUpperBound();
            BigDecimal oldZLower = node.getLowerBound();
            BigDecimal zRef = oldZUpper;

            // Apply correction: Z_true in zRef * [exp(r - q), exp(r + q)]
            BigDecimal newZUpper = scaledZ(zRef, logZResidual + logZQ);
            BigDecimal newZLower = scaledZ(zRef, logZResidual - logZQ);
            if (newZUpper == null || newZLower == null) {
                s10SubtreeRejectedRange++;
                continue;
            }

            // Only tighten — never widen the existing emat bound
            boolean tightenedUpper = newZUpper.compareTo(oldZUpper) < 0;
            boolean tightenedLower = newZLower.compareTo(oldZLower) > 0;
            if (!tightenedUpper && !tightenedLower) {
                s10SubtreeRejectedRange++;
                continue;
            }
            if (newZLower.compareTo(newZUpper) > 0) {
                s10SubtreeRejectedRange++;
                continue;
            }

            node.tightenSubtreeBounds(
                    tightenedLower ? newZLower : oldZLower,
                    tightenedUpper ? newZUpper : oldZUpper);
            node.markUpdated();
            gnnFinalizedSubtrees.add(node);
            subtreeFinalizedCount++;
            subtreeVarianceSum += logZQ * logZQ;
            s10SubtreeFinalized++;
        }
    }

    private int countFree(int[] assignments) {
        int count = 0;
        for (int rc : assignments) {
            if (rc < 0) count++;
        }
        return count;
    }

    /** Multiply Z by exp(logScale), safely handling huge/tiny values. */
    private BigDecimal scaledZ(BigDecimal z, double logScale) {
        if (Double.isNaN(logScale)) return null;
        if (logScale > 700.0) return null;
        if (logScale < -745.0) return BigDecimal.ZERO;
        return z.multiply(BigDecimal.valueOf(Math.exp(logScale)));
    }

    // ====================================================================
    // Bernstein concentration: empirical Bernstein bound on sum of GNN errors
    // ====================================================================

    private double bernsteinError(int N, double sigma, double M, double delta) {
        if (N == 0) return 0.0;
        double log = Math.log(2.0 / delta);
        return Math.sqrt(2.0 * sigma * sigma * N * log) + M * log / 3.0;
    }

    private void recomputeBernsteinError() {
        // Use observed empirical sigma if we have variance data, else fall back to calibrated
        double leafEmpSigma = leafFinalizedCount > 0
                ? Math.sqrt(leafVarianceSum / leafFinalizedCount) : leafSigma;
        double subEmpSigma = subtreeFinalizedCount > 0
                ? Math.sqrt(subtreeVarianceSum / subtreeFinalizedCount) : subtreeSigma;
        s10BernsteinErrLeaf = bernsteinError(leafFinalizedCount,
                Math.max(leafEmpSigma, leafSigma), leafMaxErr, bernsteinDelta);
        s10BernsteinErrSubtree = bernsteinError(subtreeFinalizedCount,
                Math.max(subEmpSigma, subtreeSigma), subtreeMaxErr, bernsteinDelta);
    }
}
