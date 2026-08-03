package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.bench.ConfSpaces2RL0;
import edu.duke.cs.osprey.branchdp.BranchDPScorer;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.lute.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.util.*;

/**
 * Phase 6: BranchMARK* + BranchBBK* Correctness and Performance Tests
 *
 * Compares two algorithm combinations:
 *
 * 1. MARK* + BBK* (original):
 *    - BBK* uses 1000-conf sampling for MultiSequenceNode scoring
 *    - MARK* (MARKStarBoundFastQueues) for partition function computation
 *
 * 2. BranchMARK* + BranchBBK* (improved):
 *    - BranchBBK* uses BranchDPScorer (branch decomposition DP) for MultiSequenceNode scoring
 *    - BranchMARK* (one-shot DP + Phase 2 minimization) for partition function computation
 *
 * Key correctness criteria:
 * 1. Both algorithms should find the same top-N sequences (same ranking)
 * 2. K* scores should agree within expected tolerance
 * 3. BranchBBK* should be faster due to exact DP scoring replacing 1000-conf sampling
 */
public class TestBranchMARKStar {

    public static final int NUM_CPUs = Integer.getInteger("osprey.branchdp.numCpus", 4);

    // ========== Test: MARK* + BBK* vs BranchMARK* + BranchBBK* ==========



    /**
     * BranchMARK* standalone test (no BBK*, just partition function).
     * Compares standard MARK* vs BranchMARK* on a single wild-type sequence.
     */
    @Test
    public void testBranchMARKStarOnly() {
        int numFlexible = 7;
        double epsilon = 0.68;

        long t0 = System.currentTimeMillis();
        MARKStarResult standardResult = runMARKStarOnly(numFlexible, epsilon, false, "Standard-MARK*");
        long t1 = System.currentTimeMillis();
        MARKStarResult branchResult = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*");
        long t2 = System.currentTimeMillis();

        long stdTime = t1 - t0;
        long branchTime = t2 - t1;
        System.out.println("========== TIMING COMPARISON ==========");
        System.out.println("Standard-MARK* runtime: " + stdTime + " ms");
        System.out.println("BranchMARK*    runtime: " + branchTime + " ms");
        if (stdTime > 0) {
            System.out.println("Speedup (std/branch):   " + String.format("%.2f", (double) stdTime / branchTime) + "x");
        }
        System.out.println("========================================");
    }

    /**
     * BranchMARK* sparse-pfunc standalone run: no original MARK* comparison.
     */
    @Test
    public void testBranchMARKStarSparseOnly() {
        int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
        double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));

        String oldEnergyMode = System.getProperty("branchdp.energyMode");
        System.setProperty("branchdp.energyMode", "sparse");
        try {
            System.out.println("========== BranchMARK* SPARSE PFUNC ONLY ==========");
            long t0 = System.currentTimeMillis();
            MARKStarResult branchResult = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*-SparseOnly");
            long t1 = System.currentTimeMillis();

            System.out.println("BranchMARK* sparse runtime: " + (t1 - t0) + " ms");
            if (!branchResult.scores.isEmpty()) {
                KStarScore score = branchResult.scores.get(0).score;
                System.out.println("BranchMARK* sparse K*: [" + String.format("%.6f", score.lowerBoundLog10())
                        + ", " + String.format("%.6f", score.upperBoundLog10()) + "]");
            }
            System.out.println("==================================================");
        } finally {
            if (oldEnergyMode == null) {
                System.clearProperty("branchdp.energyMode");
            } else {
                System.setProperty("branchdp.energyMode", oldEnergyMode);
            }
        }
    }

    /**
     * Fallback test: dense graph should fall back to standard algorithms.
     */
    @Test
    public void testFallbackOnDenseGraph() {
        int numFlexible = 2;
        double epsilon = 0.68;

        MARKStarResult result = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*-Fallback");
    }


    /**
     * Compare FLAT_SUM BranchMARK* vs original MARK* on the same problem.
     * Both should produce valid K* scores; compare speed and accuracy.
     */
    @Test
    public void testFlatSumVsOriginalMARKStar() {
        int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
        double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));

        System.out.println("========== Original MARK* ==========");
        long t0 = System.currentTimeMillis();
        MARKStarResult markstarResult = runMARKStarOnly(numFlexible, epsilon, false, "Original-MARK*");
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* FLAT_SUM ==========");
        long t2 = System.currentTimeMillis();
        MARKStarResult flatSumResult = runMARKStarOnly(numFlexible, epsilon, true,
                "BranchMARK*");
        long t3 = System.currentTimeMillis();

        long markstarTime = t1 - t0;
        long flatSumTime = t3 - t2;

        System.out.println("\n========== COMPARISON ==========");
        System.out.println("Original MARK*        runtime: " + markstarTime + " ms");
        System.out.println("BranchMARK* FLAT_SUM  runtime: " + flatSumTime + " ms");
        if (markstarTime > 0) {
            System.out.println("Speedup (MARK*/FlatSum): " + String.format("%.2f", (double) markstarTime / flatSumTime) + "x");
        }

        if (!markstarResult.scores.isEmpty() && !flatSumResult.scores.isEmpty()) {
            KStarScore msScore = markstarResult.scores.get(0).score;
            KStarScore fsScore = flatSumResult.scores.get(0).score;
            System.out.println("Original MARK*       K*: [" + String.format("%.6f", msScore.lowerBoundLog10())
                    + ", " + String.format("%.6f", msScore.upperBoundLog10()) + "]");
            System.out.println("BranchMARK* FLAT_SUM K*: [" + String.format("%.6f", fsScore.lowerBoundLog10())
                    + ", " + String.format("%.6f", fsScore.upperBoundLog10()) + "]");
        }
        System.out.println("================================");
    }

    /**
     * Compare the legacy BranchMARK* edge selector (largest lambda-state count)
     * against certificate-gap lookahead on the same wild-type pfunc benchmark.
     */
    @Test
    public void testBranchEdgeSelectionComparison() {
        int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
        double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));

        System.out.println("========== BranchMARK* edge selection: lambdaStates ==========");
        long t0 = System.currentTimeMillis();
        MARKStarResult legacyResult = runBranchWithEdgeSelection(numFlexible, epsilon,
                "lambdaStates", "BranchMARK*-lambdaStates");
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* edge selection: contraction ==========");
        long t2 = System.currentTimeMillis();
        MARKStarResult contractionResult = runBranchWithEdgeSelection(numFlexible, epsilon,
                "contraction", "BranchMARK*-contraction");
        long t3 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* edge selection: contractionPerState ==========");
        long t4 = System.currentTimeMillis();
        MARKStarResult contractionPerStateResult = runBranchWithEdgeSelection(numFlexible, epsilon,
                "contractionPerState", "BranchMARK*-contractionPerState");
        long t5 = System.currentTimeMillis();

        long legacyTime = t1 - t0;
        long contractionTime = t3 - t2;
        long contractionPerStateTime = t5 - t4;

        System.out.println("\n========== EDGE SELECTION COMPARISON ==========");
        System.out.println("BranchMARK* lambdaStates        runtime: " + legacyTime + " ms");
        System.out.println("BranchMARK* contraction         runtime: " + contractionTime + " ms");
        System.out.println("BranchMARK* contractionPerState runtime: " + contractionPerStateTime + " ms");
        if (contractionTime > 0) {
            System.out.println("Speedup (lambdaStates/contraction): "
                    + String.format("%.2f", (double) legacyTime / contractionTime) + "x");
        }
        if (contractionPerStateTime > 0) {
            System.out.println("Speedup (lambdaStates/contractionPerState): "
                    + String.format("%.2f", (double) legacyTime / contractionPerStateTime) + "x");
        }

        if (!legacyResult.scores.isEmpty()
                && !contractionResult.scores.isEmpty()
                && !contractionPerStateResult.scores.isEmpty()) {
            KStarScore legacyScore = legacyResult.scores.get(0).score;
            KStarScore contractionScore = contractionResult.scores.get(0).score;
            KStarScore contractionPerStateScore = contractionPerStateResult.scores.get(0).score;
            System.out.println("BranchMARK* lambdaStates        K*: ["
                    + String.format("%.6f", legacyScore.lowerBoundLog10())
                    + ", " + String.format("%.6f", legacyScore.upperBoundLog10()) + "]");
            System.out.println("BranchMARK* contraction         K*: ["
                    + String.format("%.6f", contractionScore.lowerBoundLog10())
                    + ", " + String.format("%.6f", contractionScore.upperBoundLog10()) + "]");
            System.out.println("BranchMARK* contractionPerState K*: ["
                    + String.format("%.6f", contractionPerStateScore.lowerBoundLog10())
                    + ", " + String.format("%.6f", contractionPerStateScore.upperBoundLog10()) + "]");
        }
        System.out.println("===============================================");
    }

    private MARKStarResult runBranchWithEdgeSelection(int numFlexible, double epsilon,
                                                       String edgeSelection, String label) {
        String oldEdgeSelection = System.getProperty("branchdp.edgeSelection");
        System.setProperty("branchdp.edgeSelection", edgeSelection);
        try {
            return runMARKStarOnly(numFlexible, epsilon, true, label);
        } finally {
            if (oldEdgeSelection == null) {
                System.clearProperty("branchdp.edgeSelection");
            } else {
                System.setProperty("branchdp.edgeSelection", oldEdgeSelection);
            }
        }
    }


    // ========== MARK*-only runner (no BBK*, single sequence) ==========

    public static class MARKStarResult {
        public String label;
        public PartitionFunction.Status status;
        public List<MARKStar.ScoredSequence> scores;
    }

    public static MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition, String label) {
        return runMARKStarOnly(numFlexible, epsilon, useBranchDecomposition, false, label);
    }

    public static MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition,
                                            boolean useGridDP, String label) {
        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildWildTypeConfSpace(numFlexible);

        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
            return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
        };

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setShowPfuncProgress(false)
                .setParallelism(parallelism)
                .setUseBranchDecomposition(useBranchDecomposition)
                .setUseGridDP(useGridDP);
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        markstar.precalcEmats();
        List<MARKStar.ScoredSequence> scores = markstar.run();

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        MARKStarResult result = new MARKStarResult();
        result.label = label;
        result.scores = scores;
        if (!scores.isEmpty()) {
            KStarScore firstScore = scores.get(0).score;
            result.status = firstScore.complex.status;
        } else {
            result.status = PartitionFunction.Status.Estimating;
        }

        return result;
    }



    /** Direct entry point for SLURM — bypasses gradle/JUnit entirely. */
    public static void main(String[] args) throws Exception {
        String mode = args.length > 0 ? args[0] : "benchmark";
        switch (mode) {
            case "export4pos":
                new edu.duke.cs.osprey.markstar.bench.GNNDataExport().exportGNNDataAllMutable();
                break;
            case "exportHighrot8pos":
                new edu.duke.cs.osprey.markstar.bench.GNNDataExport().exportGNNDataHighrot8pos();
                break;
            case "patchRigidHighrot8pos":
                new edu.duke.cs.osprey.markstar.bench.GNNDataExport().patchRigidEmatHighrot8pos();
                break;
            case "egfr_5pos":
                new edu.duke.cs.osprey.markstar.bench.EGFRBench().exportGNNDataEGFR5pos();
                break;
            case "scaling_bench":
                new edu.duke.cs.osprey.markstar.bench.ScalingBench().benchmarkScaling8pos();
                break;
            case "scaling_bench_12pos":
                new edu.duke.cs.osprey.markstar.bench.ScalingBench().benchmarkScaling12pos();
                break;
            case "scaling_bench_16pos":
                new edu.duke.cs.osprey.markstar.bench.ScalingBench().benchmarkScaling16pos();
                break;
            case "scaling_bench_20pos":
                new edu.duke.cs.osprey.markstar.bench.ScalingBench().benchmarkScaling20pos();
                break;
            case "scaling_bench_24pos":
                new edu.duke.cs.osprey.markstar.bench.ScalingBench().benchmarkScaling24pos();
                break;
            case "scaling_pac":
                // PACK* paper Table 2 ("scaling with n"): MARK* vs PACK* on
                // ConfSpaces2RL0.buildWildTypeConfSpace(numFlexible), n in
                // {8,10,12,16,20}. See PackStarScalingBench for property list.
                new edu.duke.cs.osprey.markstar.bench.PackStarScalingBench().benchmarkScalingN();
                break;
            case "flatsum": {
                int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
                double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));
                System.out.println("========== Original MARK* ==========");
                long t0 = System.currentTimeMillis();
                MARKStarResult ms = runMARKStarOnly(numFlexible, epsilon, false, "Original-MARK*");
                long t1 = System.currentTimeMillis();
                System.out.println("\n========== BranchMARK* FLAT_SUM ==========");
                long t2 = System.currentTimeMillis();
                MARKStarResult br = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*");
                long t3 = System.currentTimeMillis();
                long msTime = t1 - t0, brTime = t3 - t2;
                System.out.println("\n========== COMPARISON ==========");
                System.out.println("Original MARK*        runtime: " + msTime + " ms");
                System.out.println("BranchMARK* FLAT_SUM  runtime: " + brTime + " ms");
                if (msTime > 0 && brTime > 0) {
                    System.out.println("Speedup (MARK*/FlatSum): "
                            + String.format("%.2f", (double) msTime / brTime) + "x");
                }
                if (!ms.scores.isEmpty() && !br.scores.isEmpty()) {
                    KStarScore msScore = ms.scores.get(0).score;
                    KStarScore brScore = br.scores.get(0).score;
                    System.out.println("Original MARK*       K*: ["
                            + String.format("%.6f", msScore.lowerBoundLog10()) + ", "
                            + String.format("%.6f", msScore.upperBoundLog10()) + "]");
                    System.out.println("BranchMARK* FLAT_SUM K*: ["
                            + String.format("%.6f", brScore.lowerBoundLog10()) + ", "
                            + String.format("%.6f", brScore.upperBoundLog10()) + "]");
                }
                System.out.println("================================");
                break;
            }
            case "markstar": {
                int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
                double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));
                System.out.println("========== Original MARK* ONLY (branchmarkstar disabled) ==========");
                System.out.println("numFlexible=" + numFlexible + ", epsilon=" + epsilon);
                System.out.println("leafProfile=" + System.getProperty("markstar.leafProfile", "(default false)")
                        + ", leafProfile.allConfs=" + System.getProperty("markstar.leafProfile.allConfs", "(default false)"));
                long t0 = System.currentTimeMillis();
                MARKStarResult ms = runMARKStarOnly(numFlexible, epsilon, false, "Original-MARK*");
                long t1 = System.currentTimeMillis();
                System.out.println("Original MARK* runtime: " + (t1 - t0) + " ms");
                if (!ms.scores.isEmpty()) {
                    KStarScore score = ms.scores.get(0).score;
                    System.out.println("Original MARK* K*: ["
                            + String.format("%.6f", score.lowerBoundLog10()) + ", "
                            + String.format("%.6f", score.upperBoundLog10()) + "]");
                }
                System.out.println("===================================================================");
                break;
            }
            case "branchsparse": {
                int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
                double epsilon = Double.parseDouble(System.getProperty("branchdp.test.epsilon", "0.68"));
                System.setProperty("branchdp.energyMode", "sparse");
                System.out.println("========== BranchMARK* SPARSE PFUNC ONLY ==========");
                System.out.println("numFlexible=" + numFlexible + ", epsilon=" + epsilon);
                System.out.println("residualBudget=" + System.getProperty("branchdp.cutoff.residualBudget", "(default)"));
                long t0 = System.currentTimeMillis();
                MARKStarResult br = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*-SparseOnly");
                long t1 = System.currentTimeMillis();
                System.out.println("BranchMARK* sparse runtime: " + (t1 - t0) + " ms");
                if (!br.scores.isEmpty()) {
                    KStarScore score = br.scores.get(0).score;
                    System.out.println("BranchMARK* sparse K*: ["
                            + String.format("%.6f", score.lowerBoundLog10()) + ", "
                            + String.format("%.6f", score.upperBoundLog10()) + "]");
                }
                System.out.println("==================================================");
                break;
            }
            case "benchmark":
            default:
                new edu.duke.cs.osprey.markstar.bench.GNNStrategyBench().benchmarkGNNStrategies();
                break;
        }
    }

}
