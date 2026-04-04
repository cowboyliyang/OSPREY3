package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.IdentityHashMap;

/**
 * GNN batch integration strategies for MARK*.
 *
 * Strategy 4 (GNN_SCAN_CCD):
 *   A*-ordered enumeration → batch GNN → convergence check on Z →
 *   selective CCD on top Boltzmann contributors. Bypasses tree search entirely.
 *
 * Strategy 5 (HYBRID_SUBTREE):
 *   Within the MARK* tree loop, dynamically enumerate small subtrees,
 *   batch GNN predict, collapse resolved subtrees. Falls back to normal
 *   expansion for large subtrees.
 *
 * Strategy 6 (GNN_CP_POOL):
 *   Preserves MARK*'s tree expansion. Leaf nodes produced by expansion are
 *   GNN-predicted with Conformal Prediction bounds and placed in a "GNN pool".
 *   Two-bucket decision each round: internalZ vs gnnPoolZ.
 *   CCD is only performed on pool nodes when pool gap dominates.
 *   Small subtrees (numConfs ≤ remaining budget) are enumerated into pool.
 *   A probability budget (δ/α) limits simultaneous GNN-bounded nodes;
 *   CCD releases budget slots. Pool-full → leaves go back to main queue.
 */
public class MARKStarBoundGNNBatch extends MARKStarBoundFastQueues {

    public enum GNNStrategy {
        /** Original single-conf GNN (current behavior) */
        SINGLE,
        /** Strategy 4: A* scan + batch GNN + convergence + selective CCD */
        GNN_SCAN_CCD,
        /** Strategy 5: Hybrid subtree enumeration + batch GNN within tree loop */
        HYBRID_SUBTREE,
        /** Strategy 6: GNN+CP pool integrated with MARK* tree expansion */
        GNN_CP_POOL
    }

    private GNNStrategy strategy = GNNStrategy.SINGLE;
    private GNNConfEnergyCalculator gnnBatchCalc;

    // Statistics
    private final AtomicLong gnnConfsEvaluated = new AtomicLong(0);
    private final AtomicLong ccdConfsEvaluated = new AtomicLong(0);
    private final AtomicLong subtreesCollapsed = new AtomicLong(0);
    private boolean strategy4Done = false; // Strategy 4 runs once in compute()

    // Strategy 4 parameters
    private static final int S4_BATCH_SIZE = 50000;       // GNN batch size per A* round
    private static final int S4_CCD_BATCH = 50;           // CCD this many confs per round before checking epsilon

    // Strategy 5 parameters
    private static final int SUBTREE_MAX_SIZE = 50000;
    private static final int BATCH_TARGET_SIZE = 100000;
    private static final double SIGNIFICANCE_THRESHOLD = 1e-10;

    private static final double kT = 1.9872036e-3 * 298.15; // kcal/mol at 298K

    // Strategy 6: GNN batch parameters
    private static int GNN_MINI_BATCH = 50;               // fire ONNX inference every N good leaves (must be ≤ budgetMax)
    private static final double GNN_EMAT_MIN_CUTOFF = -20.0;   // skip leaves with ematMinimized.confE > this
    private static final double GNN_EMAT_RIGID_CUTOFF = 0.0;  // skip leaves with rigidEmat.confE > this

    // ========================================================================
    // Strategy 6 (GNN_CP_POOL) parameters and state
    // ========================================================================

    // Conformal Prediction parameters — calibrated from validation set
    // Val set CP quantiles (distribution-aligned, finite-sample corrected):
    //   Protein: α=0.01 → q=0.0198, α=0.001 → q=0.0549, max=0.531
    //   Complex: α=0.01 → q=0.0232, α=0.001 → q=0.0582, max=0.384
    // Use α=0.001 (99.9% coverage) as default for safety, take max across states.
    // Budget = floor(δ/α): with δ=0.05, α=0.001 → budget=50
    private double cpAlpha = 0.001;    // per-prediction miscoverage rate (from val set)
    private double cpDelta = 0.10;     // total pfunc failure probability
    private double cpQ = 0.06;         // CP quantile from val set (α=0.001, conservative)
                                       // max(protein=0.0549, complex=0.0582) rounded up
    private int budgetMax = 100;       // floor(δ/α) = floor(0.1/0.001), recomputed in setCPParams
    private int budgetUsed = 0;        // current number of GNN-bounded nodes in pool

    // Track which nodes have been GNN-bounded (for budget accounting)
    private final Set<MARKStarNode> gnnBoundedNodes = Collections.newSetFromMap(new IdentityHashMap<>());

    // Pending leaves: accumulate across rounds, fire every GNN_MINI_BATCH
    private final List<MARKStarNode> pendingGNNLeaves = new ArrayList<>();
    private int s6BudgetExhaustedRound = -1; // round when budget was first exhausted

    // Statistics
    private int s6GNNBounded = 0;
    // CCD verification of GNN accuracy during search
    private static final int CCD_VERIFY_LIMIT = 200;
    private int ccdVerifyCount = 0;
    private final List<Double> ccdVerifyErrors = new ArrayList<>();
    private int s6CCDFromGNN = 0;     // CCD'd nodes that were GNN-bounded (freed budget)
    private int s6CCDFromOriginal = 0; // CCD'd nodes with original bounds
    private double s6GNNTimeMs = 0;   // total GNN inference time in ms
    private int s6OnnxCalls = 0;

    // GNN accuracy tracking: record GNN prediction for nodes before CCD
    private final Map<MARKStarNode, Double> gnnPredictions = new IdentityHashMap<>();
    private final List<double[]> gnnVsCcdErrors = new ArrayList<>(); // [gnnPred, ccdEnergy, error]

    // Optional LUTE comparator for accuracy benchmarking
    private LUTEConfEnergyCalculator luteEcalc = null;
    private final List<Double> luteVerifyErrors = new ArrayList<>();

    public void setLUTECalculator(LUTEConfEnergyCalculator luteEcalc) {
        this.luteEcalc = luteEcalc;
        System.out.println("[LUTE_VERIFY] LUTE comparator enabled");
    }

    public MARKStarBoundGNNBatch(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                  EnergyMatrix minimizingEmat,
                                  ConfEnergyCalculator minimizingConfEcalc,
                                  RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
    }

    public void setGNNStrategy(GNNStrategy s) { this.strategy = s; }

    public void setGNNMiniBatch(int size) {
        GNN_MINI_BATCH = size;
        System.out.println("[Strategy6] GNN_MINI_BATCH=" + size);
    }

    public void setGNNBatchCalculator(GNNConfEnergyCalculator calc) {
        this.gnnBatchCalc = calc;
        // For Strategy 6, do NOT set gnnCalc on parent — we handle GNN ourselves in pool logic.
        // For other strategies, keep existing behavior.
        if (strategy != GNNStrategy.GNN_CP_POOL) {
            setGNNCalculator(calc);
        }
    }

    /** Set Conformal Prediction parameters for Strategy 6. */
    public void setCPParams(double alpha, double delta, double q) {
        this.cpAlpha = alpha;
        this.cpDelta = delta;
        this.cpQ = q;
        this.budgetMax = (int) Math.floor(delta / alpha);
        System.out.println("[Strategy6] CP params: α=" + alpha + ", δ=" + delta
                + ", q=" + String.format("%.6f", q) + " kcal/mol, budget=" + budgetMax);
    }

    public long getGNNConfsEvaluated() { return gnnConfsEvaluated.get(); }
    public long getCCDConfsEvaluated() { return ccdConfsEvaluated.get(); }
    public long getSubtreesCollapsed() { return subtreesCollapsed.get(); }

    @Override
    public void compute(int maxNumConfs) {
        if (strategy == GNNStrategy.GNN_CP_POOL && gnnBatchCalc != null) {
            // Strategy 6 uses normal tree loop but prints stats at end
            super.compute(maxNumConfs);
            System.out.println("[Strategy6] Done. eps=" + String.format("%.6f", epsilonBound)
                    + ", GNNbounded=" + s6GNNBounded
                    + ", CCD(fromGNN)=" + s6CCDFromGNN
                    + ", CCD(fromOrig)=" + s6CCDFromOriginal
                    + ", budgetUsed=" + budgetUsed + "/" + budgetMax
                    + ", budgetExhaustedRound=" + s6BudgetExhaustedRound
                    + ", gnnTime=" + String.format("%.1fms", s6GNNTimeMs)
                    + ", onnxCalls=" + s6OnnxCalls);
            // Print GNN accuracy stats for this pfunc
            if (!gnnVsCcdErrors.isEmpty()) {
                double sumAbsErr = 0, sumErr = 0, maxAbsErr = 0;
                for (double[] e : gnnVsCcdErrors) {
                    double absErr = Math.abs(e[2]);
                    sumAbsErr += absErr;
                    sumErr += e[2];
                    maxAbsErr = Math.max(maxAbsErr, absErr);
                }
                int n = gnnVsCcdErrors.size();
                // Sort errors for percentile reporting
                double[] absErrors = new double[n];
                for (int i = 0; i < n; i++) absErrors[i] = Math.abs(gnnVsCcdErrors.get(i)[2]);
                Arrays.sort(absErrors);
                System.out.println(String.format("[GNN accuracy] n=%d, MAE=%.8f, bias=%.8f, maxErr=%.8f, p50=%.8f, p90=%.8f, p99=%.8f kcal/mol",
                        n, sumAbsErr / n, sumErr / n, maxAbsErr,
                        absErrors[n / 2], absErrors[(int)(n * 0.9)], absErrors[Math.min(n - 1, (int)(n * 0.99))]));
                // Print first 3 examples
                for (int i = 0; i < Math.min(3, n); i++) {
                    double[] e = gnnVsCcdErrors.get(i);
                    System.out.println(String.format("  [example] GNN=%.8f, CCD=%.8f, err=%.8f", e[0], e[1], e[2]));
                }
            }
            // Reset per-pfunc state for next sequence
            gnnBoundedNodes.clear();
            gnnPredictions.clear();
            gnnVsCcdErrors.clear();
            pendingGNNLeaves.clear();
            s6BudgetExhaustedRound = -1;
            budgetUsed = 0;
            s6GNNBounded = 0;
            s6CCDFromGNN = 0;
            s6CCDFromOriginal = 0;
            s6GNNTimeMs = 0;
            s6OnnxCalls = 0;
            return;
        }
        if (strategy == GNNStrategy.GNN_SCAN_CCD && gnnBatchCalc != null && !strategy4Done) {
            computeWithStrategy4();
            strategy4Done = true;
            values.qstar = rootNode.getLowerBound();
            values.pstar = rootNode.getUpperBound();
            values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
            return;
        }
        super.compute(maxNumConfs);
    }

    @Override
    protected void tightenBoundInPhases() {
        if (strategy == GNNStrategy.GNN_CP_POOL) {
            tightenWithStrategy6();
        } else if (strategy == GNNStrategy.HYBRID_SUBTREE) {
            tightenWithStrategy5();
        } else {
            super.tightenBoundInPhases();
        }
    }

    // ========================================================================
    // Strategy 6: GNN+CP bounds within MARK*'s original two-bucket flow
    //
    // Minimal change from original:
    //   1. After populateQueues, GNN-bound leaf nodes (within budget) to tighten
    //      their Z gap. This makes leafZ smaller → more internal expansion rounds.
    //   2. leafZ vs internalZ decision is identical to original.
    //   3. When leaf round triggers, sort by Z gap and CCD numThreads at once.
    //   4. CCD of GNN-bounded node frees its budget slot.
    // ========================================================================

    private static int s6RoundCounter = 0;

    private void tightenWithStrategy6() {
        s6RoundCounter++;
        long t0 = System.nanoTime();
        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};

        Stopwatch loopWatch = new Stopwatch().start();
        Stopwatch phaseWatch = new Stopwatch();
        double epsilonBeforePhase = epsilonBound;

        // Step 1: populateQueues — pulls maxMinimizations leaves + 1000 internals
        System.out.println(String.format("[S6 r%d] START queueSize=%d eps=%.6f", s6RoundCounter, queue.size(), epsilonBound));
        System.out.flush();
        long t1 = System.nanoTime();
        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        long t1b = System.nanoTime();
        System.out.println(String.format("[S6 r%d] populateQ done: %.1fms internals=%d leaves=%d queueRemain=%d",
            s6RoundCounter, (t1b-t1)/1e6, internalNodes.size(), leafNodes.size(), queue.size()));
        System.out.flush();
        updateBound();
        long t1c = System.nanoTime();
        internalZ = ZSums[0];
        leafZ = ZSums[1];

        // Step 2: No longer pull bulk leaves from leafQueue.
        //         Only the leaves from populateQueues (top Z-gap) go to GNN.
        long t2 = System.nanoTime();

        // Step 3: Filter by emat energy, accumulate good leaves, fire every GNN_MINI_BATCH
        System.out.println(String.format("[S6 r%d] step2 done: %.1fms", s6RoundCounter, (System.nanoTime()-t2)/1e6));
        System.out.flush();
        long t3 = System.nanoTime();
        if (gnnBatchCalc != null && budgetUsed < budgetMax) {
            // Collect promising unbounded leaves (filter by emat energy)
            int filtered = 0, skippedEnergy = 0;
            for (MARKStarNode leaf : leafNodes) {
                if (gnnBoundedNodes.contains(leaf)) continue;
                if (budgetUsed + pendingGNNLeaves.size() >= budgetMax) break;
                if (passesGNNEmatFilter(leaf.getConfSearchNode().assignments)) {
                    pendingGNNLeaves.add(leaf);
                    filtered++;
                } else {
                    skippedEnergy++;
                }
            }
            // Also scan leafQueue for promising leaves (stop after scanning a limited number)
            List<MARKStarNode> scannedBack = new ArrayList<>();
            int scanned = 0;
            while (!leafQueue.isEmpty() && scanned < 200
                    && budgetUsed + pendingGNNLeaves.size() < budgetMax) {
                MARKStarNode extra = leafQueue.poll();
                scanned++;
                if (passesGNNEmatFilter(extra.getConfSearchNode().assignments)) {
                    pendingGNNLeaves.add(extra);
                    leafNodes.add(extra);
                    leafZ = leafZ.add(extra.getUpperBound().subtract(extra.getLowerBound()));
                    filtered++;
                } else {
                    scannedBack.add(extra);
                    skippedEnergy++;
                }
            }
            leafQueue.addAll(scannedBack);

            if (filtered > 0 || skippedEnergy > 0) {
                System.out.println(String.format("[S6 r%d] gnnFilter: accepted=%d skippedHigh=%d pending=%d budget=%d/%d",
                    s6RoundCounter, filtered, skippedEnergy, pendingGNNLeaves.size(), budgetUsed, budgetMax));
            }

            // Fire GNN batch when we have enough good leaves
            if (pendingGNNLeaves.size() >= GNN_MINI_BATCH) {
                long gnnStart = System.nanoTime();
                flushPendingGNN();
                s6GNNTimeMs += (System.nanoTime() - gnnStart) / 1e6;

                // Recompute leafZ after GNN tightening
                leafZ = BigDecimal.ZERO;
                for (MARKStarNode leaf : leafNodes) {
                    leafZ = leafZ.add(leaf.getUpperBound().subtract(leaf.getLowerBound()));
                }

                // Check if GNN bounds alone converge epsilon
                updateBound();
                if (epsilonBound <= targetEpsilon) {
                    queue.addAll(leafNodes);
                    queue.addAll(internalNodes);
                    return;
                }
            }
        }
        long t4 = System.nanoTime();
        System.out.println(String.format("[S6 r%d] step3 done: %.1fms eps=%.6f", s6RoundCounter, (t4-t3)/1e6, epsilonBound));
        System.out.flush();
        // Step 4: Two-bucket decision — identical to original MARK*
        boolean isLeafRound = MathTools.isLessThan(internalZ, leafZ);
        if (isLeafRound) {
            // === Leaf (CCD) round ===
            phaseWatch.start();

            // Sort by Z gap descending: CCD the highest-gap leaves first
            long tSort0 = System.nanoTime();
            leafNodes.sort((a, b) -> {
                BigDecimal gapA = a.getUpperBound().subtract(a.getLowerBound());
                BigDecimal gapB = b.getUpperBound().subtract(b.getLowerBound());
                return gapB.compareTo(gapA);
            });
            long tSort1 = System.nanoTime();

            int ccdCount = Math.min(parallelism.numThreads, leafNodes.size());
            long tProc0 = System.nanoTime();
            for (int i = 0; i < ccdCount; i++) {
                MARKStarNode leafNode = leafNodes.get(i);
                processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
                if (gnnBoundedNodes.remove(leafNode)) {
                    budgetUsed--;
                    s6CCDFromGNN++;
                } else {
                    s6CCDFromOriginal++;
                }
            }
            long tProc1 = System.nanoTime();
            loopTasks.waitForFinish();
            long tWait1 = System.nanoTime();

            // Record GNN vs CCD accuracy
            for (int i = 0; i < ccdCount; i++) {
                MARKStarNode leafNode = leafNodes.get(i);
                Double gnnPred = gnnPredictions.remove(leafNode);
                if (gnnPred != null) {
                    double ccdEnergy = leafNode.getConfSearchNode().getConfLowerBound();
                    double error = gnnPred - ccdEnergy;
                    gnnVsCcdErrors.add(new double[]{gnnPred, ccdEnergy, error});
                }
            }

            for (int i = ccdCount; i < leafNodes.size(); i++) {
                queue.add(leafNodes.get(i));
            }
            queue.addAll(internalNodes);

            phaseWatch.stop();
            totalLeafTime += phaseWatch.getTimeS();
            totalLeafRounds++;

            // TIMING: leaf round detail
            System.out.println(String.format("[TIMING round %d] LEAF eps=%.6f ccdCount=%d leafNodes=%d internals=%d | " +
                "sort=%.1fms proc=%.1fms wait=%.1fms phase=%.3fs",
                s6RoundCounter, epsilonBound, ccdCount, leafNodes.size(), internalNodes.size(),
                (tSort1-tSort0)/1e6, (tProc1-tProc0)/1e6, (tWait1-tProc1)/1e6, phaseWatch.getTimeS()));
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

        long t5 = System.nanoTime();

        if (epsilonBound <= targetEpsilon) {
            double epsilonReduction = Math.max(0, epsilonBeforePhase - epsilonBound);
            if (isLeafRound) {
                totalLeafEpsilonReduction += epsilonReduction;
            } else {
                totalInternalEpsilonReduction += epsilonReduction;
            }
            // TIMING: converged
            System.out.println(String.format("[TIMING round %d] CONVERGED total=%.3fs populateQ=%.1fms updateBound=%.1fms step2=%.1fms step3=%.1fms step4=%.3fs",
                s6RoundCounter, (t5-t0)/1e9, (t1b-t1)/1e6, (t1c-t1b)/1e6, (t3-t2)/1e6, (t4-t3)/1e6, (t5-t4)/1e9));
            return;
        }

        long t6a = System.nanoTime();
        loopCleanup(newNodes, loopWatch, Math.max(internalNodes.size(), leafNodes.size()));
        long t6b = System.nanoTime();
        double epsilonReduction = Math.max(0, epsilonBeforePhase - epsilonBound);
        if (isLeafRound) {
            totalLeafEpsilonReduction += epsilonReduction;
        } else {
            totalInternalEpsilonReduction += epsilonReduction;
        }

        // TIMING: full round summary
        System.out.println(String.format("[TIMING round %d] %s eps=%.6f total=%.3fs populateQ=%.1fms updateBound=%.1fms step2=%.1fms step3=%.1fms step4=%.3fs cleanup=%.3fs queueSize=%d newNodes=%d",
            s6RoundCounter, isLeafRound ? "LEAF" : "INTERNAL", epsilonBound,
            (t6b-t0)/1e9, (t1b-t1)/1e6, (t1c-t1b)/1e6, (t3-t2)/1e6, (t4-t3)/1e6, (t5-t4)/1e9, (t6b-t6a)/1e9,
            queue.size(), newNodes.size()));
    }

    /**
     * Flush all pending GNN leaves in ONE batch ONNX call.
     * Apply CP bounds to each, mark as GNN-bounded.
     */

    /** Pre-filter: returns true if this conf is suitable for GNN (low energy region). */
    private boolean passesGNNEmatFilter(int[] assignments) {
        double ematE = minimizingEmat.confE(assignments);
        if (ematE >= GNN_EMAT_MIN_CUTOFF) return false;
        double rigidE = rigidEmat.confE(assignments);
        if (rigidE >= GNN_EMAT_RIGID_CUTOFF) return false;
        return true;
    }

    private void flushPendingGNN() {
        if (pendingGNNLeaves.isEmpty()) return;

        // Deduplicate: skip any that got GNN-bounded already (e.g. from a prior flush)
        List<MARKStarNode> toBound = new ArrayList<>();
        for (MARKStarNode leaf : pendingGNNLeaves) {
            if (!gnnBoundedNodes.contains(leaf) && budgetUsed < budgetMax) {
                toBound.add(leaf);
            }
        }
        pendingGNNLeaves.clear();

        if (toBound.isEmpty()) return;

        // ONE batch ONNX call for all accumulated leaves
        int[][] confs = new int[toBound.size()][];
        for (int i = 0; i < toBound.size(); i++) {
            confs[i] = toBound.get(i).getConfSearchNode().assignments;
        }
        double[] gnnEnergies = gnnBatchCalc.calcEnergies(confs);
        gnnConfsEvaluated.addAndGet(confs.length);
        s6OnnxCalls++;

        int bounded = 0, skipped = 0;
        double sumOldGap = 0, sumNewGap = 0;
        double maxOldGap = 0;
        double exLower = 0, exUpper = 0, exGNN = 0, exNewLower = 0, exNewUpper = 0;
        // Track min/max GNN energy and energy distribution
        double minGNNEnergy = Double.MAX_VALUE, maxGNNEnergy = -Double.MAX_VALUE;
        double minOldLower = Double.MAX_VALUE;
        int countNeg = 0, countLt100 = 0, countLt1000 = 0, countGt1000 = 0;
        for (int i = 0; i < toBound.size(); i++) {
            MARKStarNode leafNode = toBound.get(i);
            MARKStarNode.Node node = leafNode.getConfSearchNode();
            double eGNN = gnnEnergies[i];

            double oldLower = node.getConfLowerBound();
            double oldUpper = node.getConfUpperBound();
            double newLower = Math.max(eGNN - cpQ, oldLower);
            double newUpper = Math.min(eGNN + cpQ, oldUpper);

            if (newLower > newUpper) {
                skipped++;
                double ematE = minimizingEmat.confE(leafNode.getConfSearchNode().assignments);
                double residual = eGNN - ematE;
                if (skipped <= 5) {
                    System.out.println(String.format("[GNN_SKIP] conf=%s oldLower=%.4f oldUpper=%.4f eGNN=%.4f emat=%.4f residual=%.4f",
                        java.util.Arrays.toString(leafNode.getConfSearchNode().assignments),
                        oldLower, oldUpper, eGNN, ematE, residual));
                }
                continue;
            }

            double oldGap = oldUpper - oldLower;
            double newGap = newUpper - newLower;
            sumOldGap += oldGap;
            sumNewGap += newGap;
            if (oldGap > maxOldGap) {
                maxOldGap = oldGap;
                exLower = oldLower; exUpper = oldUpper; exGNN = eGNN;
                exNewLower = newLower; exNewUpper = newUpper;
            }
            if (eGNN < minGNNEnergy) minGNNEnergy = eGNN;
            if (eGNN > maxGNNEnergy) maxGNNEnergy = eGNN;
            if (oldLower < minOldLower) minOldLower = oldLower;
            if (eGNN < 0) countNeg++;
            else if (eGNN < 100) countLt100++;
            else if (eGNN < 1000) countLt1000++;
            else countGt1000++;

            leafNode.setBoundsFromConfLowerAndUpper(newLower, newUpper);
            node.gscore = newLower;
            leafNode.markUpdated();

            gnnBoundedNodes.add(leafNode);
            gnnPredictions.put(leafNode, eGNN);
            budgetUsed++;
            s6GNNBounded++;
            bounded++;
        }
        System.out.println(String.format("[GNN_BATCH] bounded=%d skipped=%d avgOldGap=%.4f avgNewGap=%.4f maxOldGap=%.4f",
            bounded, skipped,
            bounded > 0 ? sumOldGap/bounded : 0,
            bounded > 0 ? sumNewGap/bounded : 0,
            maxOldGap));
        System.out.println(String.format("[GNN_BATCH] worst example: old=[%.4f, %.4f] gap=%.4f → GNN=%.4f → new=[%.4f, %.4f] gap=%.4f",
            exLower, exUpper, exUpper-exLower, exGNN, exNewLower, exNewUpper, exNewUpper-exNewLower));
        System.out.println(String.format("[GNN_BATCH] energy dist: minGNN=%.4f maxGNN=%.4f minOldLower=%.4f | neg=%d <100=%d <1000=%d >=1000=%d",
            minGNNEnergy, maxGNNEnergy, minOldLower, countNeg, countLt100, countLt1000, countGt1000));

        // CCD verification: sample up to 20 bounded confs, run CCD, compare with GNN
        if (bounded > 0 && ccdVerifyCount < CCD_VERIFY_LIMIT) {
            int toVerify = Math.min(20, toBound.size());
            for (int i = 0; i < toVerify && ccdVerifyCount < CCD_VERIFY_LIMIT; i++) {
                MARKStarNode leafNode = toBound.get(i);
                if (!gnnBoundedNodes.contains(leafNode)) continue; // was skipped
                try {
                    int[] assignments = leafNode.getConfSearchNode().assignments;
                    RCTuple tuple = new RCTuple(assignments);
                    double ccdEnergy = minimizingEcalc.calcEnergy(tuple).energy;
                    double gnnE = gnnEnergies[i];
                    double gnnErr = gnnE - ccdEnergy;
                    ccdVerifyErrors.add(Math.abs(gnnErr));

                    // Also compute LUTE prediction if available
                    String luteStr = "";
                    if (luteEcalc != null) {
                        try {
                            double luteE = luteEcalc.calcEnergy(assignments);
                            double luteErr = luteE - ccdEnergy;
                            luteVerifyErrors.add(Math.abs(luteErr));
                            luteStr = String.format(" LUTE=%.6f luteErr=%.6f", luteE, luteErr);
                        } catch (Exception luteEx) {
                            luteStr = " LUTE=PRUNED";
                        }
                    }

                    if (ccdVerifyCount < 10) {
                        System.out.println(String.format("[CCD_VERIFY] conf=%s GNN=%.6f CCD=%.6f err=%.6f%s",
                            java.util.Arrays.toString(assignments), gnnE, ccdEnergy, gnnErr, luteStr));
                    }
                    ccdVerifyCount++;
                } catch (Exception ex) {
                    // skip on error
                }
            }
            if (ccdVerifyCount > 0 && ccdVerifyCount % 50 == 0 || ccdVerifyCount == CCD_VERIFY_LIMIT) {
                double[] sorted = ccdVerifyErrors.stream().mapToDouble(Double::doubleValue).sorted().toArray();
                int n = sorted.length;
                System.out.println(String.format("[CCD_VERIFY_SUMMARY] n=%d GNN: MAE=%.6f P50=%.6f P90=%.6f P99=%.6f max=%.6f",
                    n, ccdVerifyErrors.stream().mapToDouble(Double::doubleValue).average().orElse(0),
                    sorted[n/2], sorted[(int)(n*0.9)], sorted[Math.min(n-1,(int)(n*0.99))], sorted[n-1]));

                if (!luteVerifyErrors.isEmpty()) {
                    double[] luteSorted = luteVerifyErrors.stream().mapToDouble(Double::doubleValue).sorted().toArray();
                    int ln = luteSorted.length;
                    System.out.println(String.format("[CCD_VERIFY_SUMMARY] n=%d LUTE: MAE=%.6f P50=%.6f P90=%.6f P99=%.6f max=%.6f",
                        ln, luteVerifyErrors.stream().mapToDouble(Double::doubleValue).average().orElse(0),
                        luteSorted[ln/2], luteSorted[(int)(ln*0.9)], luteSorted[Math.min(ln-1,(int)(ln*0.99))], luteSorted[ln-1]));
                }
            }
        }

        System.out.flush();
        if (budgetUsed >= budgetMax && s6BudgetExhaustedRound < 0) {
            s6BudgetExhaustedRound = s6RoundCounter;
            System.out.println(String.format("[GNN_BATCH] Budget exhausted at round %d (used=%d/%d)",
                s6RoundCounter, budgetUsed, budgetMax));
        }
    }

    /**
     * After the batch ONNX call, apply GNN bounds to any new leaves
     * that weren't in the original batch (from subsequent expansion).
     * These need a new ONNX call, but should be rare.
     */
    private void applyGNNBoundsFromCache(List<MARKStarNode> leafNodes) {
        List<MARKStarNode> newLeaves = new ArrayList<>();
        for (MARKStarNode leaf : leafNodes) {
            if (!gnnBoundedNodes.contains(leaf) && budgetUsed < budgetMax) {
                newLeaves.add(leaf);
            }
        }
        if (newLeaves.isEmpty()) return;

        // New leaves from expansion after the main batch — one more ONNX call
        int[][] confs = new int[newLeaves.size()][];
        for (int i = 0; i < newLeaves.size(); i++) {
            confs[i] = newLeaves.get(i).getConfSearchNode().assignments;
        }
        double[] gnnEnergies = gnnBatchCalc.calcEnergies(confs);
        gnnConfsEvaluated.addAndGet(confs.length);
        s6OnnxCalls++;

        for (int i = 0; i < newLeaves.size(); i++) {
            MARKStarNode leafNode = newLeaves.get(i);
            MARKStarNode.Node node = leafNode.getConfSearchNode();
            double eGNN = gnnEnergies[i];

            double oldLower = node.getConfLowerBound();
            double oldUpper = node.getConfUpperBound();
            double newLower = Math.max(eGNN - cpQ, oldLower);
            double newUpper = Math.min(eGNN + cpQ, oldUpper);

            if (newLower > newUpper) continue;

            leafNode.setBoundsFromConfLowerAndUpper(newLower, newUpper);
            node.gscore = newLower;
            leafNode.markUpdated();

            gnnBoundedNodes.add(leafNode);
            gnnPredictions.put(leafNode, eGNN);
            budgetUsed++;
            s6GNNBounded++;
        }
    }

    // ========================================================================
    // Strategy 4: A* scan → batch GNN → convergence → selective CCD
    // ========================================================================

    private void computeWithStrategy4() {
        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        Stopwatch timer = new Stopwatch().start();

        // Step 1: Build A* tree for enumeration (ordered by emat lower bound)
        ConfAStarTree astar = new ConfAStarTree.Builder(minimizingEmat, RCs)
                .setTraditional()
                .build();

        // Step 2: Enumerate ALL confs via A*, batch GNN predict
        // For ~100 sequence space this is tractable; for larger spaces
        // we'd use the convergence check to stop early
        BigDecimal zGNN = BigDecimal.ZERO;
        List<double[]> confData = new ArrayList<>(); // [gnnEnergy, ematEnergy] per conf
        List<int[]> allConfs = new ArrayList<>();
        boolean exhausted = false;

        while (!exhausted) {
            List<int[]> batch = new ArrayList<>(S4_BATCH_SIZE);
            double lastEmatBound = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < S4_BATCH_SIZE; i++) {
                ConfSearch.ScoredConf sc = astar.nextConf();
                if (sc == null) { exhausted = true; break; }
                batch.add(sc.getAssignments());
                lastEmatBound = sc.getScore();
            }
            if (batch.isEmpty()) break;

            int[][] batchArray = batch.toArray(new int[0][]);
            double[] gnnEnergies = gnnBatchCalc.calcEnergies(batchArray);
            gnnConfsEvaluated.addAndGet(batchArray.length);

            for (int i = 0; i < batchArray.length; i++) {
                zGNN = zGNN.add(bc.calc(gnnEnergies[i]));
                confData.add(new double[]{gnnEnergies[i], minimizingEmat.confE(batchArray[i])});
                allConfs.add(batchArray[i]);
            }

            // Early stop for large spaces: if remaining Z contribution is negligible
            if (!exhausted && lastEmatBound > Double.NEGATIVE_INFINITY) {
                BigDecimal maxBoltzPerConf = bc.calc(lastEmatBound);
                BigDecimal numRemaining = new BigDecimal(rootNode.getConfSearchNode().numConfs)
                        .subtract(BigDecimal.valueOf(allConfs.size()));
                BigDecimal remainingZ = numRemaining.multiply(maxBoltzPerConf);
                if (MathTools.isGreaterThan(zGNN, BigDecimal.ZERO)) {
                    double ratio = MathTools.bigDivide(remainingZ, zGNN,
                            PartitionFunction.decimalPrecision).doubleValue();
                    if (ratio < targetEpsilon * 0.01) { // well below target
                        System.out.println("[Strategy4] GNN scan converged: " + allConfs.size()
                                + " confs, remaining ratio=" + String.format("%.2e", ratio));
                        break;
                    }
                }
            }
        }

        System.out.println("[Strategy4] GNN scanned " + allConfs.size() + " confs"
                + (exhausted ? " (exhausted)" : " (converged)")
                + ", GNN Z=" + String.format("%.6e", zGNN));

        // Step 3: Sort by GNN Boltzmann weight (low energy first = high weight first)
        Integer[] sortedIdx = new Integer[allConfs.size()];
        for (int i = 0; i < sortedIdx.length; i++) sortedIdx[i] = i;
        Arrays.sort(sortedIdx, (a, b) -> Double.compare(confData.get(a)[0], confData.get(b)[0]));

        // Step 4: Epsilon-driven CCD
        // Z_upper = Σ_enumerated exp(-E_emat/kT) + Z_remaining_bound
        //   As we CCD conf i, replace its emat contribution with CCD contribution.
        // Z_lower = Σ_CCD'd exp(-E_CCD/kT)
        // epsilon = (Z_upper - Z_lower) / Z_upper
        //
        // We CCD confs in GNN-priority order until epsilon <= target.

        // Pre-compute Z_upper from per-conf emat Boltzmann weights
        BigDecimal[] ematBoltzArr = new BigDecimal[allConfs.size()];
        BigDecimal zUpperEnum = BigDecimal.ZERO; // Z_upper for enumerated confs
        for (int i = 0; i < allConfs.size(); i++) {
            ematBoltzArr[i] = bc.calc(confData.get(i)[1]); // emat energy → Boltzmann
            zUpperEnum = zUpperEnum.add(ematBoltzArr[i]);
        }
        // Remaining Z bound for unenumerated confs (0 if exhausted)
        BigDecimal zUpperRemaining = BigDecimal.ZERO;
        if (!exhausted) {
            // Conservative: rootNode upper - enumerated upper
            zUpperRemaining = rootNode.getUpperBound().subtract(zUpperEnum);
            if (MathTools.isLessThan(zUpperRemaining, BigDecimal.ZERO))
                zUpperRemaining = BigDecimal.ZERO;
        }

        BigDecimal zUpper = zUpperEnum.add(zUpperRemaining);
        BigDecimal zLower = BigDecimal.ZERO;
        int ccdCount = 0;

        System.out.println("[Strategy4] Initial Z_upper=" + String.format("%.6e", zUpper)
                + ", Z_upper_enum=" + String.format("%.6e", zUpperEnum)
                + ", Z_remaining=" + String.format("%.6e", zUpperRemaining));

        for (int rank = 0; rank < sortedIdx.length; rank++) {
            // Check epsilon
            if (MathTools.isGreaterThan(zUpper, BigDecimal.ZERO)) {
                double eps = MathTools.bigDivide(zUpper.subtract(zLower), zUpper,
                        PartitionFunction.decimalPrecision).doubleValue();
                if (eps <= targetEpsilon) {
                    System.out.println("[Strategy4] Epsilon converged: " + String.format("%.6f", eps)
                            + " after " + ccdCount + " CCD minimizations");
                    epsilonBound = eps;
                    rootNode.setSubtreeBounds(zLower, zUpper);
                    break;
                }
                if (ccdCount > 0 && ccdCount % 100 == 0) {
                    System.out.println("[Strategy4] CCD progress: " + ccdCount
                            + " confs, eps=" + String.format("%.6f", eps));
                }
            }

            // CCD this conf (GNN-priority order)
            int idx = sortedIdx[rank];
            int[] conf = allConfs.get(idx);
            edu.duke.cs.osprey.confspace.RCTuple tuple = new edu.duke.cs.osprey.confspace.RCTuple(conf);

            double ccdEnergy;
            try {
                ccdEnergy = minimizingEcalc.calcEnergy(tuple).energy;
            } catch (Exception e) {
                ccdEnergy = confData.get(idx)[0]; // fallback to GNN
            }
            ccdCount++;

            BigDecimal ccdBoltz = bc.calc(ccdEnergy);
            // Z_lower grows: add CCD Boltzmann weight
            zLower = zLower.add(ccdBoltz);
            // Z_upper tightens: replace emat Boltzmann with CCD Boltzmann for this conf
            // ccdBoltz <= ematBoltz (since ccdEnergy >= ematEnergy)
            BigDecimal ematBoltz = ematBoltzArr[idx];
            if (MathTools.isGreaterThan(ematBoltz, ccdBoltz)) {
                zUpper = zUpper.subtract(ematBoltz).add(ccdBoltz);
            }

            confData.get(idx)[0] = ccdEnergy;
        }

        ccdConfsEvaluated.addAndGet(ccdCount);

        // Final epsilon if loop exhausted without converging
        BigDecimal finalUpper = zUpper;
        BigDecimal finalLower = zLower;
        if (epsilonBound > targetEpsilon) {
            if (MathTools.isGreaterThan(finalUpper, BigDecimal.ZERO)) {
                epsilonBound = MathTools.bigDivide(finalUpper.subtract(finalLower), finalUpper,
                        PartitionFunction.decimalPrecision).doubleValue();
            }
            rootNode.setSubtreeBounds(finalLower, finalUpper);
        }

        timer.stop();
        System.out.println("[Strategy4] Done. eps=" + String.format("%.6f", epsilonBound)
                + ", time=" + timer.getTime(2)
                + ", GNN=" + gnnConfsEvaluated.get()
                + ", CCD=" + ccdConfsEvaluated.get()
                + " (" + String.format("%.1f%%", 100.0 * ccdCount / Math.max(1, allConfs.size())) + " of scanned)");

        if (epsilonBound <= targetEpsilon) {
            setStatus(Status.Estimated);
        }
    }

    // ========================================================================
    // Strategy 5: Hybrid Subtree Enumeration + Batch GNN
    // ========================================================================

    private void tightenWithStrategy5() {
        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};

        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        updateBound();
        internalZ = ZSums[0];
        leafZ = ZSums[1];

        if (gnnBatchCalc == null) {
            super.tightenBoundInPhases();
            return;
        }

        // Collect work: leaves + small subtrees → batch, large subtrees → expand
        List<int[]> batchConfs = new ArrayList<>();
        List<BatchConfOwner> batchOwners = new ArrayList<>();
        List<MARKStarNode> deferredExpand = new ArrayList<>();

        // All leaves go into batch
        for (MARKStarNode leaf : leafNodes) {
            batchConfs.add(leaf.getConfSearchNode().assignments);
            batchOwners.add(new BatchConfOwner(leaf, null));
        }

        // Check internal nodes: enumerate small subtrees, defer large ones
        for (MARKStarNode internal : internalNodes) {
            MARKStarNode.Node node = internal.getConfSearchNode();
            long subtreeSize = computeSubtreeSize(node);

            BigDecimal zFrac = BigDecimal.ZERO;
            if (MathTools.isGreaterThan(rootNode.getUpperBound(), BigDecimal.ZERO)) {
                zFrac = MathTools.bigDivide(internal.getUpperBound(),
                        rootNode.getUpperBound(), PartitionFunction.decimalPrecision);
            }

            boolean shouldEnumerate = subtreeSize > 0
                    && subtreeSize <= SUBTREE_MAX_SIZE
                    && zFrac.doubleValue() > SIGNIFICANCE_THRESHOLD
                    && batchConfs.size() + subtreeSize <= BATCH_TARGET_SIZE;

            if (shouldEnumerate) {
                List<int[]> subtreeConfs = enumerateSubtree(node);
                for (int[] conf : subtreeConfs) {
                    batchConfs.add(conf);
                    batchOwners.add(new BatchConfOwner(null, internal));
                }
            } else {
                deferredExpand.add(internal);
            }
        }

        // Batch GNN predict
        if (!batchConfs.isEmpty()) {
            int[][] confArray = batchConfs.toArray(new int[0][]);
            double[] energies = gnnBatchCalc.calcEnergies(confArray);
            gnnConfsEvaluated.addAndGet(confArray.length);

            BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

            // Group subtree confs
            Map<MARKStarNode, List<Double>> subtreeEnergies = new LinkedHashMap<>();

            for (int i = 0; i < energies.length; i++) {
                BatchConfOwner owner = batchOwners.get(i);

                if (owner.leafNode != null) {
                    // Single leaf: update bounds
                    MARKStarNode leafNode = owner.leafNode;
                    MARKStarNode.Node node = leafNode.getConfSearchNode();
                    double energy = energies[i];
                    double oldUpper = node.getConfUpperBound();
                    double oldLower = node.getConfLowerBound();
                    double newUpper = Math.min(energy, oldUpper);
                    double newLower = Math.min(energy, oldUpper);
                    leafNode.setBoundsFromConfLowerAndUpper(newLower, newUpper);
                    node.gscore = newLower;
                    leafNode.markUpdated();
                    synchronized (this) { recordReduction(oldLower, oldUpper, energy); }
                    newNodes.add(leafNode);
                } else {
                    subtreeEnergies.computeIfAbsent(owner.subtreeRoot, k -> new ArrayList<>())
                            .add(energies[i]);
                }
            }

            // Collapse subtrees: compute Z_subtree, set as exact bound
            for (Map.Entry<MARKStarNode, List<Double>> entry : subtreeEnergies.entrySet()) {
                MARKStarNode subtreeRoot = entry.getKey();
                List<Double> eList = entry.getValue();

                // Compute Z_subtree via log-sum-exp
                double minE = Double.MAX_VALUE;
                for (double e : eList) if (e < minE) minE = e;

                BigDecimal zSubtree = BigDecimal.ZERO;
                for (double e : eList) {
                    double boltz = Math.exp(-(e - minE) / kT);
                    zSubtree = zSubtree.add(new BigDecimal(boltz));
                }
                zSubtree = zSubtree.multiply(bc.calc(minE));

                // Collapse subtree: set lower = upper = Z_subtree
                // This removes all epsilon contribution from this subtree
                subtreeRoot.setExactSubtreeZ(zSubtree);

                // Propagate: mark ancestors as updated so rootNode recalculates
                subtreeRoot.markUpdated();
                subtreeRoot.computeEpsilonErrorBounds();
                newNodes.add(subtreeRoot);
                subtreesCollapsed.incrementAndGet();
            }
        }

        // Process deferred large subtrees via normal tree expansion
        for (MARKStarNode internalNode : deferredExpand) {
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
        numInternalNodesProcessed += deferredExpand.size();

        if (epsilonBound <= targetEpsilon) return;
        loopCleanup(newNodes, new Stopwatch().start(),
                leafNodes.size() + internalNodes.size());
    }

    // ========================================================================
    // Helpers
    // ========================================================================

    private static class BatchConfOwner {
        final MARKStarNode leafNode;
        final MARKStarNode subtreeRoot;

        BatchConfOwner(MARKStarNode leaf, MARKStarNode subtree) {
            this.leafNode = leaf;
            this.subtreeRoot = subtree;
        }
    }

    private long computeSubtreeSize(MARKStarNode.Node node) {
        long size = 1;
        for (int pos = 0; pos < node.assignments.length; pos++) {
            if (node.assignments[pos] == -1) {
                size *= RCs.getNum(pos);
                if (size > SUBTREE_MAX_SIZE * 10) return size;
            }
        }
        return size;
    }

    private List<int[]> enumerateSubtree(MARKStarNode.Node node) {
        int numPos = node.assignments.length;
        int[] base = Arrays.copyOf(node.assignments, numPos);

        List<Integer> unassigned = new ArrayList<>();
        for (int pos = 0; pos < numPos; pos++) {
            if (base[pos] == -1) unassigned.add(pos);
        }
        if (unassigned.isEmpty()) {
            return Collections.singletonList(base);
        }

        List<int[]> result = new ArrayList<>();
        int[] rcIdx = new int[unassigned.size()];
        int[] rcCnt = new int[unassigned.size()];
        for (int i = 0; i < unassigned.size(); i++) {
            rcCnt[i] = RCs.getNum(unassigned.get(i));
        }

        while (true) {
            int[] conf = Arrays.copyOf(base, numPos);
            for (int i = 0; i < unassigned.size(); i++) {
                conf[unassigned.get(i)] = RCs.get(unassigned.get(i))[rcIdx[i]];
            }
            result.add(conf);

            int carry = unassigned.size() - 1;
            while (carry >= 0) {
                rcIdx[carry]++;
                if (rcIdx[carry] < rcCnt[carry]) break;
                rcIdx[carry] = 0;
                carry--;
            }
            if (carry < 0) break;
        }
        return result;
    }
}
