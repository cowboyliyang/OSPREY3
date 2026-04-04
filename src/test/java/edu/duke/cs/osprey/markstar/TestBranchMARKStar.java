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
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound.ComputeMode;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNBatch;
import edu.duke.cs.osprey.markstar.framework.branch.BranchDPScorer;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
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
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.io.File;
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

    private static final int NUM_CPUs = 4;

    // ========== Test: MARK* + BBK* vs BranchMARK* + BranchBBK* ==========

    /**
     * Compare original MARK* + BBK* vs BranchMARK* + BranchBBK* on 2RL0.
     * All 8 design positions (4 mutable + 4 flexible) on the protein strand;
     * ligand strand has no flexibility.
     */
    @Test
    public void testBBKStarComparison_2RL0() {
        TestKStar.ConfSpaces confSpaces = build2RL0ProteinOnly();
        double epsilon = 0.99;
        int numSequences = 5;

        System.out.println("========== MARK* + BBK* (Original) ==========");
        long t0 = System.currentTimeMillis();
        BBKStarResult originalResult = runBBKStar(confSpaces, numSequences, epsilon, false);
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* + BranchBBK* (Improved) ==========");
        long t2 = System.currentTimeMillis();
        BBKStarResult branchResult = runBBKStar(confSpaces, numSequences, epsilon, true);
        long t3 = System.currentTimeMillis();

        // Print timing comparison
        long originalTime = t1 - t0;
        long branchTime = t3 - t2;
        System.out.println("\n========== TIMING COMPARISON ==========");
        System.out.println("MARK* + BBK*          runtime: " + originalTime + " ms");
        System.out.println("BranchMARK* + BranchBBK* runtime: " + branchTime + " ms");
        if (originalTime > 0) {
            System.out.println("Speedup (original/branch):  " + String.format("%.2f", (double) originalTime / branchTime) + "x");
        }

        // Print sequence comparison
        System.out.println("\n========== SEQUENCE COMPARISON ==========");
        printSequenceComparison(originalResult, branchResult);
        System.out.println("==========================================");
    }

    /**
     * Smaller test: use a confspace with fewer mutations for faster testing.
     * Protein: 1 mutation, Ligand: 1 mutation.
     */
    @Test
    public void testBBKStarComparison_Small() {
        TestKStar.ConfSpaces confSpaces = buildSmallMutantConfSpace();
        double epsilon = 0.99;
        int numSequences = 3;

        System.out.println("========== MARK* + BBK* (Original) ==========");
        long t0 = System.currentTimeMillis();
        BBKStarResult originalResult = runBBKStar(confSpaces, numSequences, epsilon, false);
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* + BranchBBK* (Improved) ==========");
        long t2 = System.currentTimeMillis();
        BBKStarResult branchResult = runBBKStar(confSpaces, numSequences, epsilon, true);
        long t3 = System.currentTimeMillis();

        long originalTime = t1 - t0;
        long branchTime = t3 - t2;
        System.out.println("\n========== TIMING COMPARISON ==========");
        System.out.println("MARK* + BBK*          runtime: " + originalTime + " ms");
        System.out.println("BranchMARK* + BranchBBK* runtime: " + branchTime + " ms");
        if (originalTime > 0) {
            System.out.println("Speedup (original/branch):  " + String.format("%.2f", (double) originalTime / branchTime) + "x");
        }

        System.out.println("\n========== SEQUENCE COMPARISON ==========");
        printSequenceComparison(originalResult, branchResult);
        System.out.println("==========================================");
    }

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
     * Fallback test: dense graph should fall back to standard algorithms.
     */
    @Test
    public void testFallbackOnDenseGraph() {
        int numFlexible = 2;
        double epsilon = 0.68;

        MARKStarResult result = runMARKStarOnly(numFlexible, epsilon, true, "BranchMARK*-Fallback");
    }

    /**
     * Phase 7: Compare original MARK* with vs without Grid DP upper bound.
     * Grid DP provides cheap upper bounds (~10ms cached) to skip CCD minimization.
     */
    @Test
    public void testGridDPPhase7() {
        int numFlexible = 10;
        double epsilon = 0.68;

        System.out.println("========== MARK* (Original, no Grid DP) ==========");
        long t0 = System.currentTimeMillis();
        MARKStarResult baseResult = runMARKStarOnly(numFlexible, epsilon, false, null, false, "MARK*-NoGridDP");
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== MARK* + Grid DP (Phase 7) ==========");
        long t2 = System.currentTimeMillis();
        MARKStarResult gridDPResult = runMARKStarOnly(numFlexible, epsilon, false, null, true, "MARK*-GridDP");
        long t3 = System.currentTimeMillis();

        long baseTime = t1 - t0;
        long gridDPTime = t3 - t2;
        System.out.println("\n========== PHASE 7 COMPARISON ==========");
        System.out.println("MARK* (no GridDP)  runtime: " + baseTime + " ms");
        System.out.println("MARK* + GridDP     runtime: " + gridDPTime + " ms");
        if (baseTime > 0) {
            System.out.println("Speedup: " + String.format("%.2f", (double) baseTime / gridDPTime) + "x");
        }

        if (!baseResult.scores.isEmpty() && !gridDPResult.scores.isEmpty()) {
            KStarScore bs = baseResult.scores.get(0).score;
            KStarScore gs = gridDPResult.scores.get(0).score;
            System.out.println("NoGridDP  K*: [" + String.format("%.6f", bs.lowerBoundLog10())
                    + ", " + String.format("%.6f", bs.upperBoundLog10()) + "]");
            System.out.println("GridDP    K*: [" + String.format("%.6f", gs.lowerBoundLog10())
                    + ", " + String.format("%.6f", gs.upperBoundLog10()) + "]");
        }
        System.out.println("=========================================");
    }

    /**
     * Compare FLAT_SUM BranchMARK* vs original MARK* on the same problem.
     * Both should produce valid K* scores; compare speed and accuracy.
     */
    @Test
    public void testFlatSumVsOriginalMARKStar() {
        int numFlexible = 10;
        double epsilon = 0.68;

        System.out.println("========== Original MARK* ==========");
        long t0 = System.currentTimeMillis();
        MARKStarResult markstarResult = runMARKStarOnly(numFlexible, epsilon, false, "Original-MARK*");
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== BranchMARK* FLAT_SUM ==========");
        long t2 = System.currentTimeMillis();
        MARKStarResult flatSumResult = runMARKStarOnly(numFlexible, epsilon, true,
                ComputeMode.FLAT_SUM, "BranchMARK*-FlatSum");
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

    // ========== BBK* runner ==========

    private static class BBKStarResult {
        String label;
        List<KStar.ScoredSequence> sequences;
    }

    /**
     * Run BBK* with either original (MARK* + 1000-conf scoring) or branch
     * (BranchMARK* + BranchDPScorer) algorithms.
     *
     * @param useBranch  false = MARK* + BBK* (original)
     *                   true  = BranchMARK* + BranchBBK* (improved)
     */
    private BBKStarResult runBBKStar(TestKStar.ConfSpaces confSpaces,
                                      int numSequences, double epsilon,
                                      boolean useBranch) {

        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            String label = useBranch ? "BranchMARK* + BranchBBK*" : "MARK* + BBK*";

            KStarScoreWriter.Formatter formatter = (KStarScoreWriter.ScoreInfo info) ->
                    String.format("[%s] seq=%s  K*=[%.6f, %.6f]  protein=%s  ligand=%s  complex=%s",
                            label,
                            info.sequence.toString(Sequence.Renderer.ResType),
                            info.kstarScore.lowerBoundLog10(),
                            info.kstarScore.upperBoundLog10(),
                            info.kstarScore.protein.toString(),
                            info.kstarScore.ligand.toString(),
                            info.kstarScore.complex.toString()
                    );

            // Configure BBK* settings (same for both variants)
            KStar.Settings kstarSettings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(1)
                    .addScoreConsoleWriter(formatter)
                    .build();
            BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
                    .setNumBestSequences(numSequences)
                    .setNumConfsPerBatch(8)
                    .build();

            BBKStar bbkstar = new BBKStar(
                    confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                    kstarSettings, bbkstarSettings);

            // Configure each confspace
            for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {
                SimpleConfSpace confSpace = (SimpleConfSpace) info.confSpace;

                // Minimized conf energy calculator
                info.confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                                .build()
                                .calcReferenceEnergies()
                        ).build();

                // Minimized energy matrix + A* factory
                EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(info.confEcalcMinimized)
                        .build()
                        .calcEnergyMatrix();
                info.confSearchFactoryMinimized = (rcs) ->
                        new ConfAStarTree.Builder(ematMinimized, rcs)
                                .setTraditional()
                                .build();

                // Rigid energy matrix + A* factory (required by BBK*)
                EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                        .setIsMinimizing(false)
                        .build();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(info.confEcalcMinimized, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build()
                        .calcEnergyMatrix();
                info.confSearchFactoryRigid = (rcs) ->
                        new ConfAStarTree.Builder(ematRigid, rcs)
                                .setTraditional()
                                .build();

                // Partition function factory: MARK* (original) or BranchMARK* (improved)
                if (useBranch) {
                    // BranchMARK* as pfunc
                    info.pfuncFactory = rcs -> {
                        BranchMARKStarBound pfunc = new BranchMARKStarBound(
                                confSpace,
                                ematRigid,
                                ematMinimized,
                                info.confEcalcMinimized,
                                rcs,
                                ecalcMinimized.parallelism
                        );
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        pfunc.init(epsilon);
                        return pfunc;
                    };

                    // BranchBBK*: set BranchDPScorer for DP-based MultiSequenceNode scoring
                    info.branchDPScorer = new BranchDPScorer(
                            confSpace, ematRigid, ematMinimized, 8.0, 0.1);
                } else {
                    // Original MARK* as pfunc
                    info.pfuncFactory = rcs -> {
                        MARKStarBound pfunc = new MARKStarBoundFastQueues(
                                confSpace,
                                ematRigid,
                                ematMinimized,
                                info.confEcalcMinimized,
                                rcs,
                                ecalcMinimized.parallelism
                        );
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        pfunc.init(epsilon);
                        return pfunc;
                    };
                    // No BranchDPScorer → original 1000-conf scoring
                }

                // Disable ConfDB for test
                info.confDBFile = null;
            }

            // Run BBK*
            BBKStarResult result = new BBKStarResult();
            result.label = label;
            result.sequences = bbkstar.run(ecalcMinimized.tasks);

            System.out.println("[" + label + "] Found " + result.sequences.size() + " sequences.");
            return result;
        }
    }

    // ========== MARK*-only runner (no BBK*, single sequence) ==========

    private static class MARKStarResult {
        String label;
        PartitionFunction.Status status;
        List<MARKStar.ScoredSequence> scores;
    }

    private MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition, String label) {
        return runMARKStarOnly(numFlexible, epsilon, useBranchDecomposition, null, false, label);
    }

    private MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition,
                                            ComputeMode computeMode, String label) {
        return runMARKStarOnly(numFlexible, epsilon, useBranchDecomposition, computeMode, false, label);
    }

    private MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition,
                                            ComputeMode computeMode, boolean useGridDP, String label) {
        TestKStar.ConfSpaces confSpaces = buildWildTypeConfSpace(numFlexible);

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
        if (computeMode != null) {
            settingsBuilder.setComputeMode(computeMode);
        }
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

    // ========== ConfSpace builders ==========

    /**
     * Build 2RL0 confspace with all design positions on the protein strand (chain A).
     * 4 mutable positions (A156, A172, A192, A193) + 4 flexible-only positions (A157, A160, A170, A190).
     * Ligand strand (chain G) has no flexibility.
     */
    private TestKStar.ConfSpaces build2RL0ProteinOnly() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        // Protein strand (chain A): 4 mutable + 4 flexible positions
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();
        // 4 mutable positions (with mutations)
        protein.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();
        // 4 flexible-only positions (wild-type rotamers only)
        protein.flexibility.get("A157").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A160").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A170").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A190").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Ligand strand (chain G): no flexibility
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build a confspace with the standard 2RL0 protein-ligand system
     * but with fewer mutations for smaller test. Protein: 1 design position,
     * Ligand: 1 design position.
     */
    private TestKStar.ConfSpaces buildSmallMutantConfSpace() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        // Protein strand (chain G): 1 design position with mutation
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Ligand strand (chain A): 1 design position with mutation
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build wild-type-only confspace with flexible residues (no mutations, no BBK*).
     * For testing BranchMARK* partition function computation in isolation.
     */
    private TestKStar.ConfSpaces buildWildTypeConfSpace(int numFlexible) {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        String[] proteinResidues = {"A156", "A157", "A158", "A159", "A160", "A161", "A162", "A163", "A164", "A165",
                "A166", "A167", "A168", "A169", "A170", "A171", "A172", "A173", "A174", "A175"};
        for (int i = 0; i < Math.min(numFlexible, proteinResidues.length); i++) {
            protein.flexibility.get(proteinResidues[i]).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    // ========== GNN Data Export ==========

    /**
     * Export conformation energy data for GNN training.
     * Samples conformations, computes CCD energies in parallel, and writes CSV files.
     */
    @Test
    public void exportGNNData() throws Exception {
        int numFlexible = 8;
        int numSamples = 5000;
        File outputDir = new File("gnn_data/2RL0_flex" + numFlexible);

        TestKStar.ConfSpaces confSpaces = buildWildTypeConfSpace(numFlexible);
        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            SimpleConfSpace confSpace = confSpaces.complex;

            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                            .build()
                            .calcReferenceEnergies()
                    ).build();

            EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                    .build()
                    .calcEnergyMatrix();

            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false)
                    .build();
            ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
            EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                    .build()
                    .calcEnergyMatrix();

            RCs rcs = new RCs(confSpace);
            InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                    confSpace, ematRigid, ematMinimized, rcs, 8.0, 0.1);

            GNNDataExporter exporter = new GNNDataExporter(
                    confEcalc, ematMinimized, ig, rcs);
            exporter.export(numSamples, outputDir);
        }
    }

    // ========== GNN + MARK* Integration ==========

    /**
     * Build conf space with 2 mutable + 11 flexible = 13 positions, rigid ligand.
     * Mutable: A156 (LYS → ALA, VAL, ILE = 4 aa), A193 (ASP → SER, ASN = 3 aa) → 12 sequences.
     * Flexible: A157-A167.
     */
    private TestKStar.ConfSpaces buildMutableConfSpace() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // 2 mutable positions (4×3 = 12 sequences)
        protein.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "ILE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();
        // 11 flexible-only positions
        protein.flexibility.get("A157").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A159").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A160").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A161").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A163").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A165").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A166").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A167").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Ligand strand (chain G): no flexibility
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build confspace with 4 mutable positions, high-rotamer AA = 625 sequences.
     * Same 4 positions as buildAllMutableConfSpace (matching GNN model with 4-pos input).
     * Uses ARG, LYS, MET, GLU, TRP — each has 20-80 rotamers → large conformation space.
     * A156: 5 AA, A164: 5 AA, A172: 5 AA, A192: 5 AA → 5^4 = 625 sequences
     */
    private TestKStar.ConfSpaces buildHighRotamerConfSpace() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // High-rotamer AAs: ARG(81), LYS(81), MET(27), GLU(27), TRP(9) + wildtype
        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();  // wt=PHE
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();  // wt=LYS
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "LEU").addWildTypeRotamers().setContinuous();  // wt=ILE

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build confspace with 4 mutable positions, limited AA = 108 sequences.
     * Same 4 positions as buildAllMutableConfSpace (matching GNN model with 4-pos input).
     * A156: 4 AA, A164: 3 AA, A172: 3 AA, A192: 3 AA → 4×3×3×3 = 108
     */
    private TestKStar.ConfSpaces buildMediumMutableConfSpace() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // 4 positions (same as all20 conf space, matching GNN model)
        protein.flexibility.get("A156").setLibraryRotamers("ALA", "VAL", "ILE", "PHE").addWildTypeRotamers().setContinuous();  // wt=PHE, 4 AA
        protein.flexibility.get("A164").setLibraryRotamers("ALA", "TYR", "PHE").addWildTypeRotamers().setContinuous();          // wt=TYR, 3 AA
        protein.flexibility.get("A172").setLibraryRotamers("LYS", "ARG", "ASP").addWildTypeRotamers().setContinuous();          // wt=LYS, 3 AA
        protein.flexibility.get("A192").setLibraryRotamers("ILE", "VAL", "LEU").addWildTypeRotamers().setContinuous();          // wt=ILE, 3 AA

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build confspace with ALL 20 AA at every position, no flexible-only.
     * 4 mutable positions × 20 AA each = 160,000 sequences.
     * ~300 RCs/pos → ~300^4 ≈ 8 billion conformations.
     */
    private TestKStar.ConfSpaces buildAllMutableConfSpace() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // All 20 standard amino acids (Strand.WildType ensures native AA is included)
        String[] others = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                           "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                           "THR", "TRP", "TYR", "VAL"};

        // 4 positions, all 20 AA each → 20^4 = 160,000 sequences
        protein.flexibility.get("A156").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=PHE
        protein.flexibility.get("A164").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A172").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=LYS
        protein.flexibility.get("A192").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=ILE

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Build confspace with ALL 20 AA at 5 positions.
     * 5 mutable positions × 20 AA each = 20^5 = 3,200,000 sequences.
     * Positions: A156, A162, A164, A172, A192 (spread across binding interface).
     */
    private TestKStar.ConfSpaces buildAllMutableConfSpace5pos() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        String[] others = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                           "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                           "THR", "TRP", "TYR", "VAL"};

        protein.flexibility.get("A156").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=PHE
        protein.flexibility.get("A162").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=LEU
        protein.flexibility.get("A164").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A172").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=LYS
        protein.flexibility.get("A192").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();  // wt=ILE

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Export GNN training data: 5 positions × 20 AA = 3.2M sequences.
     * System properties: osprey.gnn.numSamples, osprey.gnn.outputDir
     */
    @Test
    public void exportGNNDataAllMutable5pos() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_all20_5pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = buildAllMutableConfSpace5pos();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data (5pos × 20AA) ==========");
            SimpleConfSpace proteinCS = confSpaces.protein;
            System.out.println("Protein positions: " + proteinCS.positions.size());
            for (SimpleConfSpace.Position pos : proteinCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator proteinConfEcalc = new ConfEnergyCalculator.Builder(proteinCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(proteinCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix proteinEmatMin = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalc)
                    .build().calcEnergyMatrix();
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();
            ConfEnergyCalculator proteinConfEcalcRigid = new ConfEnergyCalculator(proteinConfEcalc, ecalcRigid);
            EnergyMatrix proteinEmatRigid = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs proteinRCs = new RCs(proteinCS);
            InteractionGraph proteinIG = InteractionGraph.buildWithDualCutoff(
                    proteinCS, proteinEmatRigid, proteinEmatMin, proteinRCs, 8.0, 0.1);
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data (5pos × 20AA) ==========");
            SimpleConfSpace complexCS = confSpaces.complex;
            System.out.println("Complex positions: " + complexCS.positions.size());
            for (SimpleConfSpace.Position pos : complexCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator complexConfEcalc = new ConfEnergyCalculator.Builder(complexCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(complexCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix complexEmatMin = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalc)
                    .build().calcEnergyMatrix();
            ConfEnergyCalculator complexConfEcalcRigid = new ConfEnergyCalculator(complexConfEcalc, ecalcRigid);
            EnergyMatrix complexEmatRigid = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs complexRCs = new RCs(complexCS);
            InteractionGraph complexIG = InteractionGraph.buildWithDualCutoff(
                    complexCS, complexEmatRigid, complexEmatMin, complexRCs, 8.0, 0.1);
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexIG, complexRCs)
                    .export(numSamples, new File(baseDir + "/complex"));

            System.out.println("\n========== Export complete ==========");
        }
    }

    /**
     * Build confspace with ALL 20 AA at every position, no flexible-only.
     * 11 mutable positions × 20 AA each = 20^11 sequences.
     * Positions: A156-A164, A172, A192.
     */
    private TestKStar.ConfSpaces buildAllMutableConfSpace11pos() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // All 20 standard amino acids (Strand.WildType ensures native AA is included)
        String[] others = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                           "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                           "THR", "TRP", "TYR", "VAL"};

        // 11 positions, all 20 AA each → 20^11 sequences
        protein.flexibility.get("A156").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A157").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A159").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A160").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A161").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A163").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers(others).addWildTypeRotamers().setContinuous();

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Export GNN training data: 4 positions × 20 AA = 160k sequences.
     * System properties: osprey.gnn.numSamples, osprey.gnn.outputDir
     */
    @Test
    public void exportGNNDataAllMutable() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_all20_4pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = buildAllMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data (4pos × 20AA) ==========");
            SimpleConfSpace proteinCS = confSpaces.protein;
            System.out.println("Protein positions: " + proteinCS.positions.size());
            for (SimpleConfSpace.Position pos : proteinCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator proteinConfEcalc = new ConfEnergyCalculator.Builder(proteinCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(proteinCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix proteinEmatMin = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalc)
                    .build().calcEnergyMatrix();
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();
            ConfEnergyCalculator proteinConfEcalcRigid = new ConfEnergyCalculator(proteinConfEcalc, ecalcRigid);
            EnergyMatrix proteinEmatRigid = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs proteinRCs = new RCs(proteinCS);
            InteractionGraph proteinIG = InteractionGraph.buildWithDualCutoff(
                    proteinCS, proteinEmatRigid, proteinEmatMin, proteinRCs, 8.0, 0.1);
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data (4pos × 20AA) ==========");
            SimpleConfSpace complexCS = confSpaces.complex;
            System.out.println("Complex positions: " + complexCS.positions.size());
            for (SimpleConfSpace.Position pos : complexCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator complexConfEcalc = new ConfEnergyCalculator.Builder(complexCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(complexCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix complexEmatMin = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalc)
                    .build().calcEnergyMatrix();
            ConfEnergyCalculator complexConfEcalcRigid = new ConfEnergyCalculator(complexConfEcalc, ecalcRigid);
            EnergyMatrix complexEmatRigid = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs complexRCs = new RCs(complexCS);
            InteractionGraph complexIG = InteractionGraph.buildWithDualCutoff(
                    complexCS, complexEmatRigid, complexEmatMin, complexRCs, 8.0, 0.1);
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexIG, complexRCs)
                    .export(numSamples, new File(baseDir + "/complex"));

            System.out.println("\n========== Export complete ==========");
        }
    }

    /**
     * Patch existing confs.csv with E_rigid column using rigid emat.
     * No CCD needed — just emat table lookups.
     * Usage: ./gradlew test --tests "...TestBranchMARKStar.patchRigidEmat" -Dosprey.gnn.dataDir=gnn_data/2RL0_all20_4pos_merged
     */
    @Test
    public void patchRigidEmat() throws Exception {
        String baseDir = System.getProperty("osprey.gnn.dataDir", "gnn_data/2RL0_all20_4pos_merged");

        TestKStar.ConfSpaces confSpaces = buildAllMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(4);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();

            for (String strand : new String[]{"protein", "complex"}) {
                SimpleConfSpace cs = strand.equals("protein") ? confSpaces.protein : confSpaces.complex;
                File confsFile = new File(baseDir + "/" + strand + "/confs.csv");
                if (!confsFile.exists()) {
                    System.out.println("SKIP: " + confsFile + " does not exist");
                    continue;
                }

                System.out.println("\n=== Patching " + strand + " ===");

                // Compute rigid emat
                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMinimized)
                                .build().calcReferenceEnergies())
                        .build();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build().calcEnergyMatrix();
                int numPos = cs.positions.size();

                // Stream-process CSV: read line by line, add E_rigid, write to tmp
                File tmpFile = new File(confsFile.getPath() + ".tmp");
                int count = 0;
                try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(confsFile), 1 << 20);
                     java.io.PrintWriter pw = new java.io.PrintWriter(
                             new java.io.BufferedWriter(new java.io.FileWriter(tmpFile), 1 << 20))) {

                    String header = br.readLine();
                    if (header.contains("E_rigid")) {
                        System.out.println("Already has E_rigid column, skipping.");
                        tmpFile.delete();
                        continue;
                    }
                    // Insert E_rigid after E_emat: ...,E_CCD,E_emat,residual → ...,E_CCD,E_emat,E_rigid,residual
                    pw.println(header.replace("E_emat,residual", "E_emat,E_rigid,residual"));

                    String line;
                    while ((line = br.readLine()) != null) {
                        line = line.trim();
                        if (line.isEmpty()) continue;
                        String[] parts = line.split(",");
                        int[] conf = new int[numPos];
                        for (int p = 0; p < numPos; p++) {
                            conf[p] = Integer.parseInt(parts[p]);
                        }
                        double eRigid = ematRigid.confE(conf);

                        StringBuilder sb = new StringBuilder();
                        for (int p = 0; p < numPos; p++) {
                            sb.append(parts[p]).append(',');
                        }
                        sb.append(parts[numPos]).append(',');      // E_CCD
                        sb.append(parts[numPos + 1]).append(',');  // E_emat
                        sb.append(String.format("%.6f", eRigid)).append(','); // E_rigid
                        sb.append(parts[numPos + 2]);              // residual
                        pw.println(sb);
                        count++;

                        if (count % 500000 == 0) {
                            System.out.println("  processed " + count);
                        }
                    }
                    System.out.println("  total: " + count + " conformations patched");
                }

                // Atomic rename
                java.nio.file.Files.move(tmpFile.toPath(), confsFile.toPath(),
                        java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                System.out.println("  saved: " + confsFile);
            }
        }
        System.out.println("\nDone patching E_rigid.");
    }

    /**
     * Export GNN training data: 11 positions × 20 AA.
     * Positions: A156-A164, A172, A192.
     * System properties: osprey.gnn.numSamples, osprey.gnn.outputDir
     */
    @Test
    public void exportGNNDataAllMutable11pos() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 300000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_all20_11pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = buildAllMutableConfSpace11pos();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data (11pos × 20AA) ==========");
            SimpleConfSpace proteinCS = confSpaces.protein;
            System.out.println("Protein positions: " + proteinCS.positions.size());
            for (SimpleConfSpace.Position pos : proteinCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator proteinConfEcalc = new ConfEnergyCalculator.Builder(proteinCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(proteinCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix proteinEmatMin = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalc)
                    .build().calcEnergyMatrix();
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();
            ConfEnergyCalculator proteinConfEcalcRigid = new ConfEnergyCalculator(proteinConfEcalc, ecalcRigid);
            EnergyMatrix proteinEmatRigid = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs proteinRCs = new RCs(proteinCS);
            InteractionGraph proteinIG = InteractionGraph.buildWithDualCutoff(
                    proteinCS, proteinEmatRigid, proteinEmatMin, proteinRCs, 8.0, 0.1);
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data (11pos × 20AA) ==========");
            SimpleConfSpace complexCS = confSpaces.complex;
            System.out.println("Complex positions: " + complexCS.positions.size());
            for (SimpleConfSpace.Position pos : complexCS.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }
            ConfEnergyCalculator complexConfEcalc = new ConfEnergyCalculator.Builder(complexCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(complexCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix complexEmatMin = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalc)
                    .build().calcEnergyMatrix();
            ConfEnergyCalculator complexConfEcalcRigid = new ConfEnergyCalculator(complexConfEcalc, ecalcRigid);
            EnergyMatrix complexEmatRigid = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs complexRCs = new RCs(complexCS);
            InteractionGraph complexIG = InteractionGraph.buildWithDualCutoff(
                    complexCS, complexEmatRigid, complexEmatMin, complexRCs, 8.0, 0.1);
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexIG, complexRCs)
                    .export(numSamples, new File(baseDir + "/complex"));

            System.out.println("\n========== Export complete ==========");
        }
    }

    /**
     * Export GNN training data for protein and complex conf spaces separately.
     * 13 positions (2 mutable + 11 flexible), rigid ligand, 12 sequences.
     * System properties: osprey.gnn.numSamples, osprey.gnn.outputDir (base dir)
     */
    @Test
    public void exportGNNDataForMARKStar() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 50000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_markstar_8pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = buildMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data ==========");
            SimpleConfSpace proteinCS = confSpaces.protein;
            ConfEnergyCalculator proteinConfEcalc = new ConfEnergyCalculator.Builder(proteinCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(proteinCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix proteinEmatMin = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalc)
                    .build().calcEnergyMatrix();
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();
            ConfEnergyCalculator proteinConfEcalcRigid = new ConfEnergyCalculator(proteinConfEcalc, ecalcRigid);
            EnergyMatrix proteinEmatRigid = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs proteinRCs = new RCs(proteinCS);
            InteractionGraph proteinIG = InteractionGraph.buildWithDualCutoff(
                    proteinCS, proteinEmatRigid, proteinEmatMin, proteinRCs, 8.0, 0.1);
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));
            System.out.println("Protein positions: " + proteinCS.positions.size());

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data ==========");
            SimpleConfSpace complexCS = confSpaces.complex;
            ConfEnergyCalculator complexConfEcalc = new ConfEnergyCalculator.Builder(complexCS, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(complexCS, ecalcMinimized)
                            .build().calcReferenceEnergies())
                    .build();
            EnergyMatrix complexEmatMin = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalc)
                    .build().calcEnergyMatrix();
            ConfEnergyCalculator complexConfEcalcRigid = new ConfEnergyCalculator(complexConfEcalc, ecalcRigid);
            EnergyMatrix complexEmatRigid = new SimplerEnergyMatrixCalculator.Builder(complexConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs complexRCs = new RCs(complexCS);
            InteractionGraph complexIG = InteractionGraph.buildWithDualCutoff(
                    complexCS, complexEmatRigid, complexEmatMin, complexRCs, 8.0, 0.1);
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexIG, complexRCs)
                    .export(numSamples, new File(baseDir + "/complex"));
            System.out.println("Complex positions: " + complexCS.positions.size());

            System.out.println("\n========== Export complete ==========");
            System.out.println("Protein data: " + baseDir + "/protein/");
            System.out.println("Complex data: " + baseDir + "/complex/");
        }
    }

    /**
     * Benchmark MARK* with CCD vs MARK* with GNN.
     * 8 positions (2 mutable + 6 flexible), 12 sequences, rigid ligand.
     * Compares K* scores (= Z_complex / (Z_protein × Z_ligand)).
     *
     * System properties:
     *   osprey.gnn.eval.proteinModelPath - protein GNN ONNX model
     *   osprey.gnn.eval.complexModelPath - complex GNN ONNX model
     */
    @Test
    public void benchmarkMARKStarCCDvsGNN() throws Exception {
        double epsilon = 0.68;

        String proteinModelPath = System.getProperty("osprey.gnn.eval.proteinModelPath");
        String complexModelPath = System.getProperty("osprey.gnn.eval.complexModelPath");
        if (proteinModelPath == null || complexModelPath == null) {
            throw new IllegalArgumentException("Must set -Dosprey.gnn.eval.proteinModelPath and -Dosprey.gnn.eval.complexModelPath");
        }

        // === Run 1: MARK* with GNN ===
        System.out.println("\n==============================================");
        System.out.println("  MARK* + GNN, 13 pos, 12 sequences");
        System.out.println("==============================================");
        long t0 = System.currentTimeMillis();
        MARKStarResult gnnResult = runMARKStarMutable(epsilon,
                proteinModelPath, complexModelPath, "MARK*+GNN");
        long t1 = System.currentTimeMillis();
        long gnnTime = t1 - t0;

        // === Run 2: MARK* with CCD (baseline) ===
        System.out.println("\n==============================================");
        System.out.println("  MARK* + CCD (baseline), 13 pos, 12 sequences");
        System.out.println("==============================================");
        long t2 = System.currentTimeMillis();
        MARKStarResult ccdResult = runMARKStarMutable(epsilon, null, null, "MARK*+CCD");
        long t3 = System.currentTimeMillis();
        long ccdTime = t3 - t2;

        // === Print comparison ===
        System.out.println("\n==============================================");
        System.out.println("  RESULTS COMPARISON");
        System.out.println("==============================================");
        System.out.println("Epsilon: " + epsilon);
        System.out.println("Positions: 13 (2 mutable + 11 flexible)");
        System.out.println("Sequences: 12 (4 × 3)");
        System.out.println();

        System.out.println("MARK*+GNN runtime:  " + gnnTime + " ms");
        System.out.println("MARK*+CCD runtime:  " + ccdTime + " ms");
        if (gnnTime > 0) {
            System.out.println("Speedup (CCD/GNN):  " + String.format("%.2f", (double)ccdTime / gnnTime) + "x");
        }

        System.out.println();
        printKStarComparison("MARK*+GNN", gnnResult, "MARK*+CCD", ccdResult);
        System.out.println("==============================================");
    }

    /**
     * Run MARK* with the mutable conf space (13 positions, 12 sequences).
     * If proteinModelPath and complexModelPath are non-null, use GNN; otherwise use CCD.
     */
    private MARKStarResult runMARKStarMutable(double epsilon,
                                               String proteinModelPath, String complexModelPath,
                                               String label) {
        TestKStar.ConfSpaces confSpaces = buildMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) ->
                new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                .build().calcReferenceEnergies())
                        .build();

        MARKStar.Settings settings = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(2)
                .setShowPfuncProgress(false)
                .setParallelism(parallelism)
                .build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        markstar.precalcEmats();

        // Optionally load GNN models
        GNNConfEnergyCalculator proteinGNN = null;
        GNNConfEnergyCalculator complexGNN = null;
        if (proteinModelPath != null && complexModelPath != null) {
            proteinGNN = new GNNConfEnergyCalculator(
                    new File(proteinModelPath),
                    markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size());
            markstar.protein.gnnCalc = proteinGNN;
            System.out.println("Protein GNN loaded: " + proteinModelPath
                    + " (" + confSpaces.protein.positions.size() + " positions)");

            complexGNN = new GNNConfEnergyCalculator(
                    new File(complexModelPath),
                    markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size());
            markstar.complex.gnnCalc = complexGNN;
            System.out.println("Complex GNN loaded: " + complexModelPath
                    + " (" + confSpaces.complex.positions.size() + " positions)");

            System.out.println("Ligand: no flexibility, skipping GNN");
        } else {
            System.out.println("Running with CCD (no GNN)");
        }

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

        try {
            if (proteinGNN != null) proteinGNN.close();
            if (complexGNN != null) complexGNN.close();
        } catch (Exception e) {
            System.err.println("Warning: error closing GNN calculators: " + e.getMessage());
        }

        return result;
    }

    private void printKStarComparison(String label1, MARKStarResult r1,
                                       String label2, MARKStarResult r2) {
        int n = Math.max(r1.scores.size(), r2.scores.size());
        for (int i = 0; i < n; i++) {
            System.out.println("--- Sequence " + (i + 1) + " ---");
            if (i < r1.scores.size()) {
                MARKStar.ScoredSequence s = r1.scores.get(i);
                System.out.println("  [" + label1 + "]  K*(log10): " + s.score);
            }
            if (i < r2.scores.size()) {
                MARKStar.ScoredSequence s = r2.scores.get(i);
                System.out.println("  [" + label2 + "]  K*(log10): " + s.score);
            }
            if (i < r1.scores.size() && i < r2.scores.size()) {
                double diff = Math.abs(
                        r1.scores.get(i).score.lowerBoundLog10() - r2.scores.get(i).score.lowerBoundLog10());
                System.out.println("  |diff(log10)| = " + String.format("%.6f", diff));
            }
        }
    }

    // ========== GNN Batch Strategy Benchmark ==========

    /**
     * Benchmark three approaches on 4pos × 20AA = 160K sequences:
     *   1. MARK* + single-conf GNN (current)
     *   2. Strategy 4: GNN scan + selective CCD
     *   3. Strategy 5: Hybrid subtree enumeration + batch GNN
     *
     * System properties:
     *   osprey.gnn.eval.proteinModelPath - protein ONNX model
     *   osprey.gnn.eval.complexModelPath - complex ONNX model
     */
    /**
     * Parameterized benchmark: run a single strategy per invocation.
     * System property osprey.gnn.benchStrategy selects which to run:
     *   ccd, single, strategy4, strategy5 (default: ccd)
     */
    @Test
    public void benchmarkGNNStrategies() throws Exception {
        double epsilon = 0.683;
        int numSeqs = Integer.getInteger("osprey.gnn.numSeqs", 20);
        String which = System.getProperty("osprey.gnn.benchStrategy", "ccd");

        String proteinModelPath = System.getProperty("osprey.gnn.eval.proteinModelPath",
                "gnn_data/2RL0_all20_4pos/protein/model/gnn_model.onnx");
        String complexModelPath = System.getProperty("osprey.gnn.eval.complexModelPath",
                "gnn_data/2RL0_all20_4pos/complex/model/gnn_model.onnx");

        System.out.println("==============================================");
        System.out.println("  GNN Batch Strategy Benchmark");
        System.out.println("  strategy=" + which + ", epsilon=" + epsilon + ", numSeqs=" + numSeqs);
        System.out.println("==============================================");

        String label;
        MARKStarBoundGNNBatch.GNNStrategy strategy;
        String pModel, cModel;

        switch (which) {
            case "ccd":
                label = "MARK*+CCD";
                strategy = null;
                pModel = null;
                cModel = null;
                break;
            case "single":
                label = "MARK*+GNN-single";
                strategy = null;
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "strategy4":
                label = "Strategy4:scan+CCD";
                strategy = MARKStarBoundGNNBatch.GNNStrategy.GNN_SCAN_CCD;
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "strategy5":
                label = "Strategy5:subtree";
                strategy = MARKStarBoundGNNBatch.GNNStrategy.HYBRID_SUBTREE;
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "strategy6":
                label = "Strategy6:CP-pool";
                strategy = MARKStarBoundGNNBatch.GNNStrategy.GNN_CP_POOL;
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "strategy7":
                label = "Strategy7:decoupled";
                strategy = null; // S7 uses its own class, not GNNStrategy enum
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "branch_ccd":
                label = "BranchMARK*+CCD";
                strategy = null;
                pModel = null;
                cModel = null;
                break;
            case "branch_gnn":
                label = "BranchMARK*+GNN";
                strategy = null; // uses BranchMARKStarBound's built-in GNN pool
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            default:
                throw new IllegalArgumentException("Unknown strategy: " + which
                        + ". Use one of: ccd, single, strategy4, strategy5, strategy6, strategy7, branch_ccd, branch_gnn");
        }

        long t0 = System.currentTimeMillis();
        GNNBenchResult result = runGNNBenchmark(epsilon, numSeqs,
                pModel, cModel, strategy, label);
        long elapsed = System.currentTimeMillis() - t0;

        // Print results
        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + label);
        System.out.println("==============================================");
        System.out.println(String.format("%-25s %10s %12s %12s %10s",
                "Method", "Time(s)", "GNN confs", "CCD confs", "#Seqs"));
        System.out.println("----------------------------------------------------------------------");
        printBenchRow(label, elapsed, result);

        // K* scores: filter out unconverged (Infinity bounds), sort by lower bound descending
        List<MARKStar.ScoredSequence> converged = new ArrayList<>();
        int unconvergedCount = 0;
        for (MARKStar.ScoredSequence s : result.scores) {
            double lb = s.score.lowerBoundLog10();
            double ub = s.score.upperBoundLog10();
            if (Double.isInfinite(lb) || Double.isInfinite(ub)) {
                unconvergedCount++;
            } else {
                converged.add(s);
            }
        }
        converged.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));

        System.out.println("\n=== K* Scores (top 30 converged by lower bound, descending) ===");
        System.out.println("  Total: " + result.scores.size() + " sequences, " + converged.size() + " converged, " + unconvergedCount + " unconverged");
        for (int i = 0; i < Math.min(30, converged.size()); i++) {
            MARKStar.ScoredSequence s = converged.get(i);
            System.out.println(String.format("  Rank %3d: %s  K*=[%.4f, %.4f]",
                    i + 1,
                    s.sequence.toString(Sequence.Renderer.ResType),
                    s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
        }
    }

    private void printBenchRow(String label, long timeMs, GNNBenchResult r) {
        System.out.println(String.format("%-25s %10.1f %12d %12d %10d",
                label, timeMs / 1000.0, r.gnnConfs, r.ccdConfs, r.scores.size()));
    }

    private static class GNNBenchResult {
        List<MARKStar.ScoredSequence> scores;
        long gnnConfs;
        long ccdConfs;
    }

    private GNNBenchResult runGNNBenchmark(double epsilon, int numSeqs,
                                            String proteinModelPath, String complexModelPath,
                                            MARKStarBoundGNNBatch.GNNStrategy strategy,
                                            String label) {
        // Select conf space: "medium" (108 seqs) or "highrot" (625 seqs, large conf space)
        String confSpaceType = System.getProperty("osprey.gnn.confSpace", "medium");
        TestKStar.ConfSpaces confSpaces;
        String ematDefault;
        switch (confSpaceType) {
            case "highrot":
                confSpaces = buildHighRotamerConfSpace();
                ematDefault = "emat_cache/highrot_4pos.*.dat";
                break;
            case "medium":
            default:
                confSpaces = buildMediumMutableConfSpace();
                ematDefault = "emat_cache/med10_2pos.*.dat";
                break;
        }
        int cpus = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) ->
                new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                .build().calcReferenceEnergies())
                        .build();

        // Use emat cache to avoid recomputing across runs
        String ematCachePattern = System.getProperty("osprey.gnn.ematCache", ematDefault);

        String which = System.getProperty("osprey.gnn.benchStrategy", "ccd");
        boolean useBranch = which.startsWith("branch_");

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(4)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .setEnergyMatrixCachePattern(ematCachePattern);
        if (useBranch) {
            settingsBuilder.setUseBranchDecomposition(true);
        }
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        long ematStart = System.currentTimeMillis();
        markstar.precalcEmats();
        long ematTime = System.currentTimeMillis() - ematStart;
        System.out.println("  [" + label + "] Emat time: " + ematTime + " ms (cached: " + ematCachePattern + ")");

        // Load GNN models
        GNNConfEnergyCalculator proteinGNN = null;
        GNNConfEnergyCalculator complexGNN = null;
        if (proteinModelPath != null && complexModelPath != null) {
            proteinGNN = new GNNConfEnergyCalculator(
                    new File(proteinModelPath),
                    markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size());
            markstar.protein.gnnCalc = proteinGNN;

            complexGNN = new GNNConfEnergyCalculator(
                    new File(complexModelPath),
                    markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size());
            markstar.complex.gnnCalc = complexGNN;

            // Build RC mapping if inference confspace differs from training confspace
            String gnnTrainingCS = System.getProperty("osprey.gnn.trainingConfSpace", "");
            if (!gnnTrainingCS.isEmpty() && !gnnTrainingCS.equals(confSpaceType)) {
                System.out.println("  Building RC mapping: inference=" + confSpaceType + " → training=" + gnnTrainingCS);
                TestKStar.ConfSpaces trainingConfSpaces;
                switch (gnnTrainingCS) {
                    case "all20": trainingConfSpaces = buildAllMutableConfSpace(); break;
                    case "highrot": trainingConfSpaces = buildHighRotamerConfSpace(); break;
                    case "medium": trainingConfSpaces = buildMediumMutableConfSpace(); break;
                    default: throw new IllegalArgumentException("Unknown trainingConfSpace: " + gnnTrainingCS);
                }
                int[][] proteinMap = GNNConfEnergyCalculator.buildRCMapping(confSpaces.protein, trainingConfSpaces.protein);
                int[][] complexMap = GNNConfEnergyCalculator.buildRCMapping(confSpaces.complex, trainingConfSpaces.complex);
                proteinGNN.setRCMapping(proteinMap);
                complexGNN.setRCMapping(complexMap);
            }

            // Set strategy if specified
            if (which.equals("strategy7") || which.equals("branch_gnn")) {
                markstar.protein.useStrategy7 = true;
                markstar.complex.useStrategy7 = true;
                int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
                markstar.protein.s7GPUBatchSize = gpuBatch;
                markstar.complex.s7GPUBatchSize = gpuBatch;
                System.out.println("  " + (useBranch ? "BranchMARK*+GNN" : "Strategy7") + ": gpuBatchSize=" + gpuBatch);
            } else if (strategy != null) {
                markstar.protein.gnnStrategy = strategy;
                markstar.complex.gnnStrategy = strategy;

                // GNN mini batch size (for Strategy 6)
                int miniBatch = Integer.getInteger("osprey.gnn.miniBatch", -1);
                if (miniBatch > 0) {
                    markstar.protein.gnnMiniBatch = miniBatch;
                    markstar.complex.gnnMiniBatch = miniBatch;
                    System.out.println("  GNN mini batch size: " + miniBatch);
                }
            }

            System.out.println("  GNN loaded: protein=" + confSpaces.protein.positions.size()
                    + " pos, complex=" + confSpaces.complex.positions.size() + " pos");
            if (strategy != null) {
                System.out.println("  Strategy: " + strategy);
            }
        } else {
            System.out.println("  Running with CCD (no GNN)");
        }

        // Load or train LUTE for accuracy comparison (optional)
        boolean enableLute = Boolean.parseBoolean(System.getProperty("osprey.lute.enable", "false"));
        boolean luteNoDEE = Boolean.parseBoolean(System.getProperty("osprey.lute.noDEE", "false"));
        if (enableLute) {
            double luteMaxRMSE = Double.parseDouble(System.getProperty("osprey.lute.maxRMSE", "0.1"));
            String luteDir = luteNoDEE ? "lute_bench_nodee" : "lute_bench";
            System.out.println("  Loading/training LUTE for accuracy comparison (maxRMSE=" + luteMaxRMSE
                    + ", noDEE=" + luteNoDEE + ", dir=" + luteDir + ")...");
            new File(luteDir).mkdirs();

            for (String spaceName : new String[]{"protein", "complex"}) {
                SimpleConfSpace cs = spaceName.equals("protein") ? confSpaces.protein : confSpaces.complex;
                MARKStar.ConfSpaceInfo info = spaceName.equals("protein") ? markstar.protein : markstar.complex;

                File luteFile = new File(String.format("%s/LUTE.%s.dat", luteDir, spaceName));
                LUTEConfEnergyCalculator luteEcalc;

                if (luteFile.exists()) {
                    System.out.println("    Loading cached LUTE " + spaceName + " from " + luteFile);
                    LUTEState luteState = LUTEIO.read(luteFile);
                    luteEcalc = new LUTEConfEnergyCalculator(cs, luteState);
                } else {
                    System.out.println("    Training LUTE " + spaceName + " (noDEE=" + luteNoDEE + ")...");
                    PruningMatrix pmat;
                    if (luteNoDEE) {
                        // No DEE: empty pruning matrix (nothing pruned)
                        pmat = new PruningMatrix(cs);
                    } else {
                        pmat = new SimpleDEE.Runner()
                                .setSinglesThreshold(100.0)
                                .setPairsThreshold(100.0)
                                .setGoldsteinDiffThreshold(10.0)
                                .setShowProgress(true)
                                .setCacheFile(new File(String.format("%s/LUTE.%s.pmat.dat", luteDir, spaceName)))
                                .setParallelism(parallelism)
                                .run(cs, info.minimizingEmat);
                    }

                    File confDBFile = new File(String.format("%s/LUTE.%s.conf.db", luteDir, spaceName));
                    try (ConfDB confdb = new ConfDB(cs, confDBFile)) {
                        ConfDB.ConfTable confTable = confdb.new ConfTable("lute");
                        LUTE lute = new LUTE(cs);
                        ConfSampler sampler = new RandomizedDFSConfSampler(cs, pmat, 12345);
                        lute.sampleTuplesAndFit(info.minimizingConfEcalc, info.minimizingEmat, pmat,
                                confTable, sampler, LUTE.Fitter.OLSCG, 1.5, luteMaxRMSE);
                        lute.reportConfSpaceSize(pmat);
                        lute.save(luteFile);
                        System.out.println("    LUTE " + spaceName + " saved to " + luteFile);
                        luteEcalc = new LUTEConfEnergyCalculator(cs, new LUTEState(lute.getTrainingSystem()));
                    }
                }
                info.luteCalc = luteEcalc;
            }
            System.out.println("  LUTE comparison enabled for both protein and complex");
        }

        List<MARKStar.ScoredSequence> scores = markstar.run();

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        GNNBenchResult result = new GNNBenchResult();
        result.scores = scores;
        result.gnnConfs = 0;
        result.ccdConfs = 0;

        // TODO: collect stats from pfunc instances
        // For now, approximate from the pfunc results
        System.out.println("  [" + label + "] Computed " + scores.size() + " sequences");

        try {
            if (proteinGNN != null) proteinGNN.close();
            if (complexGNN != null) complexGNN.close();
        } catch (Exception e) {
            System.err.println("Warning: error closing GNN: " + e.getMessage());
        }

        return result;
    }

    // ========== Comparison output ==========

    private void printSequenceComparison(BBKStarResult original, BBKStarResult branch) {
        int n = Math.max(original.sequences.size(), branch.sequences.size());
        for (int i = 0; i < n; i++) {
            System.out.println("--- Rank " + (i + 1) + " ---");
            if (i < original.sequences.size()) {
                KStar.ScoredSequence s = original.sequences.get(i);
                System.out.println("  [MARK*+BBK*]          " + s.sequence.toString(Sequence.Renderer.ResType)
                        + "  K*=[" + String.format("%.6f", s.score.lowerBoundLog10())
                        + ", " + String.format("%.6f", s.score.upperBoundLog10()) + "]");
            } else {
                System.out.println("  [MARK*+BBK*]          (no result)");
            }
            if (i < branch.sequences.size()) {
                KStar.ScoredSequence s = branch.sequences.get(i);
                System.out.println("  [BranchMARK*+BranchBBK*] " + s.sequence.toString(Sequence.Renderer.ResType)
                        + "  K*=[" + String.format("%.6f", s.score.lowerBoundLog10())
                        + ", " + String.format("%.6f", s.score.upperBoundLog10()) + "]");
            } else {
                System.out.println("  [BranchMARK*+BranchBBK*] (no result)");
            }
        }
    }

    // ========== LUTE vs CCD Accuracy Benchmark ==========

    /**
     * Train LUTE on the same 4pos × 20AA conf space, then sample conformations
     * and compare LUTE predicted energy vs actual CCD minimized energy.
     *
     * System properties:
     *   osprey.lute.numCPUs      - number of CPUs (default: 4)
     *   osprey.lute.numTestConfs - number of conformations to test (default: 500)
     *   osprey.lute.maxRMSE     - LUTE training target RMSE (default: 0.1)
     */
    @Test
    public void benchmarkLUTEvsCCD() throws Exception {
        int numCPUs = Integer.getInteger("osprey.lute.numCPUs", NUM_CPUs);
        int numTestConfs = Integer.getInteger("osprey.lute.numTestConfs", 500);
        double maxRMSE = Double.parseDouble(System.getProperty("osprey.lute.maxRMSE", "0.1"));

        System.out.println("==============================================");
        System.out.println("  LUTE vs CCD Accuracy Benchmark");
        System.out.println("  4pos × 20AA, numCPUs=" + numCPUs + ", numTestConfs=" + numTestConfs);
        System.out.println("==============================================");

        TestKStar.ConfSpaces confSpaces = buildAllMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        for (String spaceName : new String[]{"protein", "complex"}) {
            SimpleConfSpace confSpace = spaceName.equals("protein") ? confSpaces.protein : confSpaces.complex;

            System.out.println("\n========== " + spaceName.toUpperCase() + " ==========");
            System.out.println("Positions: " + confSpace.positions.size());
            for (SimpleConfSpace.Position pos : confSpace.positions) {
                System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
            }

            try (EnergyCalculator ecalcMin = new EnergyCalculator.Builder(confSpace, confSpaces.ffparams)
                    .setParallelism(parallelism)
                    .build()) {

                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalcMin)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMin)
                                .build().calcReferenceEnergies())
                        .build();

                // Compute energy matrix
                String ematFile = String.format("lute_bench/LUTE.%s.emat.dat", spaceName);
                new File("lute_bench").mkdirs();
                EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                        .setCacheFile(new File(ematFile))
                        .build()
                        .calcEnergyMatrix();

                // DEE pruning (important for good LUTE fits)
                String pmatFile = String.format("lute_bench/LUTE.%s.pmat.dat", spaceName);
                PruningMatrix pmat = new SimpleDEE.Runner()
                        .setSinglesThreshold(100.0)
                        .setPairsThreshold(100.0)
                        .setGoldsteinDiffThreshold(10.0)
                        .setShowProgress(true)
                        .setCacheFile(new File(pmatFile))
                        .setParallelism(parallelism)
                        .run(confSpace, emat);

                // Train LUTE (or load from cache)
                File luteFile = new File(String.format("lute_bench/LUTE.%s.dat", spaceName));
                LUTEConfEnergyCalculator luteEcalc;
                long trainTime;

                if (luteFile.exists()) {
                    System.out.println("\nLoading cached LUTE model from " + luteFile);
                    long t0 = System.currentTimeMillis();
                    LUTEState luteState = LUTEIO.read(luteFile);
                    luteEcalc = new LUTEConfEnergyCalculator(confSpace, luteState);
                    trainTime = System.currentTimeMillis() - t0;
                    System.out.println("LUTE loaded in " + trainTime + " ms");
                } else {
                    System.out.println("\nTraining LUTE (maxRMSE=" + maxRMSE + ") ...");
                    long t0 = System.currentTimeMillis();

                    File confDBFile = new File(String.format("lute_bench/LUTE.%s.conf.db", spaceName));
                    try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
                        ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

                        LUTE lute = new LUTE(confSpace);
                        ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, 12345);
                        lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler,
                                LUTE.Fitter.OLSCG, 1.5, maxRMSE);
                        lute.reportConfSpaceSize(pmat);

                        // Save for reuse
                        lute.save(luteFile);
                        System.out.println("LUTE model saved to " + luteFile);

                        luteEcalc = new LUTEConfEnergyCalculator(confSpace, new LUTEState(lute.getTrainingSystem()));
                    }
                    trainTime = System.currentTimeMillis() - t0;
                    System.out.println("LUTE training done in " + trainTime + " ms");
                }

                // Sample conformations via A* (low-energy first) and compare
                System.out.println("\nSampling " + numTestConfs + " conformations via A* and comparing...");
                ConfAStarTree astar = new ConfAStarTree.Builder(emat, pmat)
                        .setTraditional()
                        .build();

                double sumSqErr = 0;
                double sumAbsErr = 0;
                double maxAbsErr = 0;
                int count = 0;

                System.out.println(String.format("%-6s %40s %12s %12s %12s %12s",
                        "Rank", "Conf", "CCD", "LUTE", "Diff", "|Diff|"));
                System.out.println("------------------------------------------" +
                        "------------------------------------------------------------");

                for (int i = 0; i < numTestConfs; i++) {
                    ConfSearch.ScoredConf scoredConf = astar.nextConf();
                    if (scoredConf == null) break;

                    int[] assignments = scoredConf.getAssignments();

                    // Real CCD minimized energy
                    double ccdEnergy = confEcalc.calcEnergy(new RCTuple(assignments)).energy;

                    // LUTE predicted energy
                    double luteEnergy = luteEcalc.calcEnergy(assignments);

                    double diff = luteEnergy - ccdEnergy;
                    double absDiff = Math.abs(diff);

                    sumSqErr += diff * diff;
                    sumAbsErr += absDiff;
                    if (absDiff > maxAbsErr) maxAbsErr = absDiff;
                    count++;

                    if (i < 20 || i % 50 == 0 || i == numTestConfs - 1) {
                        System.out.println(String.format("%-6d %40s %12.4f %12.4f %12.4f %12.4f",
                                i + 1, Arrays.toString(assignments),
                                ccdEnergy, luteEnergy, diff, absDiff));
                    }
                }

                double rmse = Math.sqrt(sumSqErr / count);
                double mae = sumAbsErr / count;

                System.out.println("\n========== " + spaceName.toUpperCase() + " RESULTS ==========");
                System.out.println("  Conformations tested: " + count);
                System.out.println(String.format("  RMSE:     %.6f kcal/mol", rmse));
                System.out.println(String.format("  MAE:      %.6f kcal/mol", mae));
                System.out.println(String.format("  Max|err|: %.6f kcal/mol", maxAbsErr));
                System.out.println("  LUTE training time: " + trainTime + " ms");
            }
        }
        System.out.println("\n========== LUTE vs CCD benchmark complete ==========");
    }

    /** Direct entry point for SLURM — bypasses gradle/JUnit entirely. */
    public static void main(String[] args) throws Exception {
        String mode = args.length > 0 ? args[0] : "benchmark";
        switch (mode) {
            case "export4pos":
                new TestBranchMARKStar().exportGNNDataAllMutable();
                break;
            case "luteVsCCD":
                new TestBranchMARKStar().benchmarkLUTEvsCCD();
                break;
            case "benchmark":
            default:
                new TestBranchMARKStar().benchmarkGNNStrategies();
                break;
        }
    }
}
