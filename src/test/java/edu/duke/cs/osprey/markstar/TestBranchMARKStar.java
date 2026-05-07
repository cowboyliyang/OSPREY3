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
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
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
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import org.junit.jupiter.api.Test;

import java.io.File;
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
        MARKStarResult baseResult = runMARKStarOnly(numFlexible, epsilon, false, false, "MARK*-NoGridDP");
        long t1 = System.currentTimeMillis();

        System.out.println("\n========== MARK* + Grid DP (Phase 7) ==========");
        long t2 = System.currentTimeMillis();
        MARKStarResult gridDPResult = runMARKStarOnly(numFlexible, epsilon, false, true, "MARK*-GridDP");
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
        return runMARKStarOnly(numFlexible, epsilon, useBranchDecomposition, false, label);
    }

    private MARKStarResult runMARKStarOnly(int numFlexible, double epsilon,
                                            boolean useBranchDecomposition,
                                            boolean useGridDP, String label) {
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
     * Build confspace with 8 mutable positions, high-rotamer AAs.
     * Same AA palette as buildHighRotamerConfSpace (ARG, LYS, MET, GLU + wildtype)
     * but expanded to 8 positions for a more challenging benchmark.
     *
     * 8 positions × 5 AA each = 5^8 = 390,625 sequences.
     * Each position ~170-220 RCs → total conformations ~200^8 ≈ 2.56 × 10^16.
     *
     * Positions: A156(PHE), A158(HIS), A162(THR), A164(TYR),
     *            A168(GLU), A172(LYS), A174(TYR), A192(ILE)
     */
    private TestKStar.ConfSpaces buildHighRotamerConfSpace8pos() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // High-rotamer AAs: ARG, LYS, MET, GLU + wildtype at each position
        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();  // wt=PHE
        protein.flexibility.get("A158").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "HIS").addWildTypeRotamers().setContinuous();  // wt=HIS
        protein.flexibility.get("A162").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "THR").addWildTypeRotamers().setContinuous();  // wt=THR
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A168").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ASP").addWildTypeRotamers().setContinuous();  // wt=GLU (in list), +ASP for 5 AA
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();  // wt=LYS (in list), +TRP for 5 AA
        protein.flexibility.get("A174").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ILE").addWildTypeRotamers().setContinuous();  // wt=ILE

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
     * Build confspace with 8 mutable highrot positions (same as buildHighRotamerConfSpace8pos)
     * + 8 flexible-only WT positions = 16 total positions.
     *
     * Designed to study BranchMARK* vs MARK* scaling: sequence count stays the same as
     * the 8-pos benchmark (34 seqs at maxMut=1), but the per-sequence conformation space is
     * larger (16 conformational positions instead of 8), so this isolates the
     * effect of conformational complexity on the branch-decomposition speedup.
     *
     * Mutable (highrot, 5 AA each): A156, A158, A162, A164, A168, A172, A174, A192
     * Flexible (WT only, 4 pos):    A157, A159, A165, A167
     */
    private TestKStar.ConfSpaces buildHighRotamerConfSpace12pos() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // 8 mutable highrot positions (identical to buildHighRotamerConfSpace8pos)
        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ILE").addWildTypeRotamers().setContinuous();

        // 4 flexible-only WT positions (added to enlarge conformational space)
        String[] flexResidues = {"A157", "A159", "A165", "A167"};
        for (String res : flexResidues) {
            protein.flexibility.get(res).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

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
     * Build confspace with 4 mutable + 4 flexible = 8 positions.
     * Same 8 positions as buildHighRotamerConfSpace8pos, but only 4 are mutable (highrot),
     * the other 4 are flexible-only (wildtype rotamers).
     *
     * Mutable (highrot):  A156(5 AA), A164(5 AA), A172(5 AA), A192(5 AA) → 5^4 = 625 sequences
     * Flexible (wt only): A158(wt=HIS), A162(wt=THR), A168(wt=GLU), A174(wt=TYR)
     *
     * GNN models trained on full 8-pos highrot can be reused with RC mapping:
     * mutable positions map normally; flexible positions have only wt rotamers
     * which are a subset of the training confspace's rotamers.
     */
    private TestKStar.ConfSpaces buildHybridConfSpace8pos() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // 4 mutable positions (same highrot AAs as 8pos full)
        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();  // wt=PHE
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();  // wt=TYR
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();  // wt=LYS
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ILE").addWildTypeRotamers().setContinuous();  // wt=ILE

        // 4 flexible-only positions (wt rotamers only)
        protein.flexibility.get("A158").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();  // wt=HIE
        protein.flexibility.get("A162").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();  // wt=THR
        protein.flexibility.get("A168").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();  // wt=GLU
        protein.flexibility.get("A174").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();  // wt=TYR

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
     * Build the TRAINING confspace for GNN RC mapping.
     * This is the full 8-pos highrot confspace that the GNN was trained on.
     * Used only for building RC mapping between hybrid (4mut+4flex) and training confspace.
     */
    private SimpleConfSpace buildTrainingConfSpaceProtein8pos() {
        ForcefieldParams ffparams = new ForcefieldParams();
        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ILE").addWildTypeRotamers().setContinuous();

        return new SimpleConfSpace.Builder().addStrand(protein).build();
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
     * Export GNN training data for 8-position high-rotamer conf space.
     * 8 positions × 5 AA each = 390,625 sequences.
     * Exports both protein and complex data.
     *
     * System properties:
     *   osprey.gnn.numSamples - number of conformations to sample (default: 300000)
     *   osprey.gnn.outputDir  - output directory (default: gnn_data/2RL0_highrot_8pos)
     *   osprey.gnn.numCPUs    - number of CPUs (default: 4)
     */
    @Test
    public void exportGNNDataHighrot8pos() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 300000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_highrot_8pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        System.out.println("==============================================");
        System.out.println("  GNN Data Export: 8-pos highrot");
        System.out.println("  numSamples=" + numSamples + ", outputDir=" + baseDir);
        System.out.println("  numCPUs=" + numCPUs);
        System.out.println("==============================================");

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data (8pos highrot) ==========");
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
            ConfEnergyCalculator proteinConfEcalcRigid = new ConfEnergyCalculator(proteinConfEcalc, ecalcRigid);
            EnergyMatrix proteinEmatRigid = new SimplerEnergyMatrixCalculator.Builder(proteinConfEcalcRigid)
                    .build().calcEnergyMatrix();
            RCs proteinRCs = new RCs(proteinCS);
            InteractionGraph proteinIG = InteractionGraph.buildWithDualCutoff(
                    proteinCS, proteinEmatRigid, proteinEmatMin, proteinRCs, 8.0, 0.1);
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinEmatRigid, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data (8pos highrot) ==========");
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
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexEmatRigid, complexIG, complexRCs)
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
     * Patch existing highrot 8pos confs.csv: replace NaN E_rigid with real values from rigid emat.
     * No CCD needed — just emat table lookups.
     * Usage: java -cp "$CP" --add-modules jdk.incubator.foreign -Xmx64g \
     *   edu.duke.cs.osprey.markstar.TestBranchMARKStar patchRigidHighrot8pos
     */
    public void patchRigidEmatHighrot8pos() throws Exception {
        String baseDir = System.getProperty("osprey.gnn.dataDir", "gnn_data/2RL0_highrot_8pos");

        TestKStar.ConfSpaces confSpaces = buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(
                Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs));

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

                System.out.println("\n=== Patching " + strand + " (highrot 8pos) ===");

                // Compute rigid emat
                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMinimized)
                                .build().calcReferenceEnergies())
                        .build();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build().calcEnergyMatrix();
                int numPos = cs.positions.size();

                // Stream-process CSV: replace NaN E_rigid with real values
                File tmpFile = new File(confsFile.getPath() + ".tmp");
                int count = 0;
                try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(confsFile), 1 << 20);
                     java.io.PrintWriter pw = new java.io.PrintWriter(
                             new java.io.BufferedWriter(new java.io.FileWriter(tmpFile), 1 << 20))) {

                    String header = br.readLine();
                    // Header already has E_rigid column: rc_0,...,E_CCD,E_emat,E_rigid,residual
                    pw.println(header);

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

                        // parts: rc_0..rc_7, E_CCD, E_emat, E_rigid(NaN), residual
                        StringBuilder sb = new StringBuilder();
                        for (int p = 0; p < numPos; p++) {
                            sb.append(parts[p]).append(',');
                        }
                        sb.append(parts[numPos]).append(',');      // E_CCD
                        sb.append(parts[numPos + 1]).append(',');  // E_emat
                        sb.append(String.format("%.6f", eRigid)).append(','); // E_rigid (replace NaN)
                        sb.append(parts[numPos + 3]);              // residual (skip old E_rigid at numPos+2)
                        pw.println(sb);
                        count++;

                        if (count % 100000 == 0) {
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
        System.out.println("\nDone patching E_rigid (highrot 8pos).");
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
            case "strategy8":
                label = "Strategy8:S7+subtreeGNN";
                strategy = null; // S8 uses its own class
                pModel = proteinModelPath;
                cModel = complexModelPath;
                break;
            case "strategy9":
                label = "Strategy9:subtreeRouter";
                strategy = null; // S9 uses its own class
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
                        + ". Use one of: ccd, single, strategy4, strategy5, strategy6, strategy7, strategy8, strategy9, branch_ccd, branch_gnn");
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
            if (which.equals("strategy8")) {
                markstar.protein.useStrategy8 = true;
                markstar.complex.useStrategy8 = true;
                int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
                markstar.protein.s7GPUBatchSize = gpuBatch;
                markstar.complex.s7GPUBatchSize = gpuBatch;

                // Load subtree GNN models
                String proteinSubtreeModelPath = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/protein/model_subtree/subtree_model.onnx");
                String complexSubtreeModelPath = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/complex/model_subtree/subtree_model.onnx");
                GNNSubtreeEnergyCalculator proteinSubtreeGNN = new GNNSubtreeEnergyCalculator(
                        new java.io.File(proteinSubtreeModelPath), markstar.protein.minimizingEmat,
                        confSpaces.protein.positions.size());
                GNNSubtreeEnergyCalculator complexSubtreeGNN = new GNNSubtreeEnergyCalculator(
                        new java.io.File(complexSubtreeModelPath), markstar.complex.minimizingEmat,
                        confSpaces.complex.positions.size());

                // Apply RC mapping to subtree GNN if needed
                if (!gnnTrainingCS.isEmpty() && !gnnTrainingCS.equals(confSpaceType)) {
                    TestKStar.ConfSpaces trainCS;
                    switch (gnnTrainingCS) {
                        case "all20": trainCS = buildAllMutableConfSpace(); break;
                        case "highrot": trainCS = buildHighRotamerConfSpace(); break;
                        case "medium": trainCS = buildMediumMutableConfSpace(); break;
                        default: throw new IllegalArgumentException("Unknown trainingConfSpace: " + gnnTrainingCS);
                    }
                    proteinSubtreeGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.protein, trainCS.protein));
                    complexSubtreeGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.complex, trainCS.complex));
                }

                markstar.protein.subtreeGnnCalc = proteinSubtreeGNN;
                markstar.complex.subtreeGnnCalc = complexSubtreeGNN;
                System.out.println("  Strategy8: gpuBatchSize=" + gpuBatch
                        + ", subtree models loaded");
            } else if (which.equals("strategy9")) {
                markstar.protein.useStrategy9 = true;
                markstar.complex.useStrategy9 = true;
                int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
                markstar.protein.s7GPUBatchSize = gpuBatch;
                markstar.complex.s7GPUBatchSize = gpuBatch;

                // Load subtree GNN models (same as S8)
                String proteinSubtreeModelPath = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/protein/model_subtree/subtree_model.onnx");
                String complexSubtreeModelPath = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/complex/model_subtree/subtree_model.onnx");
                GNNSubtreeEnergyCalculator protSubGNN = new GNNSubtreeEnergyCalculator(
                        new java.io.File(proteinSubtreeModelPath), markstar.protein.minimizingEmat,
                        confSpaces.protein.positions.size());
                GNNSubtreeEnergyCalculator compSubGNN = new GNNSubtreeEnergyCalculator(
                        new java.io.File(complexSubtreeModelPath), markstar.complex.minimizingEmat,
                        confSpaces.complex.positions.size());

                if (!gnnTrainingCS.isEmpty() && !gnnTrainingCS.equals(confSpaceType)) {
                    TestKStar.ConfSpaces trainCS;
                    switch (gnnTrainingCS) {
                        case "all20": trainCS = buildAllMutableConfSpace(); break;
                        case "highrot": trainCS = buildHighRotamerConfSpace(); break;
                        case "medium": trainCS = buildMediumMutableConfSpace(); break;
                        default: throw new IllegalArgumentException("Unknown trainingConfSpace: " + gnnTrainingCS);
                    }
                    protSubGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.protein, trainCS.protein));
                    compSubGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.complex, trainCS.complex));
                }

                markstar.protein.subtreeGnnCalc = protSubGNN;
                markstar.complex.subtreeGnnCalc = compSubGNN;
                System.out.println("  Strategy9: gpuBatchSize=" + gpuBatch
                        + ", subtree models loaded (router mode)");
            } else if (which.equals("strategy7") || which.equals("branch_gnn")) {
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

    // ========== Hybrid 8-pos Benchmark: 4 mutable + 4 flexible ==========

    /**
     * Benchmark on 8-pos highrot confspace (same confspace GNN was trained on).
     * No RC mapping needed — confspace matches training exactly.
     *
     * Methods: markstar_ccd, markstar_s8, branch_ccd, branch_s8
     *
     * System properties:
     *   osprey.hybrid.method   — which method (default: markstar_ccd)
     *   osprey.hybrid.epsilon  — partition function epsilon (default: 0.683)
     *   osprey.hybrid.numCPUs  — number of CPUs (default: 16)
     *   osprey.hybrid.maxMut   — max simultaneous mutations (default: 1)
     *   osprey.gnn.eval.proteinModelPath         — protein leaf GNN model (8pos highrot)
     *   osprey.gnn.eval.complexModelPath         — complex leaf GNN model (8pos highrot)
     *   osprey.gnn.eval.proteinSubtreeModelPath  — protein subtree GNN model
     *   osprey.gnn.eval.complexSubtreeModelPath  — complex subtree GNN model
     *   osprey.gnn.gpuBatchSize                  — GNN batch size (default: 1000)
     */
    public void benchmarkHybrid8pos() throws Exception {
        String method = System.getProperty("osprey.hybrid.method", "markstar_ccd");
        double epsilon = Double.parseDouble(System.getProperty("osprey.hybrid.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.hybrid.numCPUs", NUM_CPUs);
        int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
        int maxMut = Integer.getInteger("osprey.hybrid.maxMut", 1);

        boolean useBranch = method.startsWith("branch_");
        boolean useGNN = method.endsWith("_s8");

        System.out.println("==============================================");
        System.out.println("  8-pos Highrot Benchmark (full confspace)");
        System.out.println("  method=" + method + ", maxMut=" + maxMut);
        System.out.println("  epsilon=" + epsilon + ", cpus=" + cpus);
        System.out.println("==============================================");

        // GNN model paths
        String proteinModelPath = null, complexModelPath = null;
        String proteinSubtreeModelPath = null, complexSubtreeModelPath = null;
        if (useGNN) {
            proteinModelPath = System.getProperty("osprey.gnn.eval.proteinModelPath");
            complexModelPath = System.getProperty("osprey.gnn.eval.complexModelPath");
            proteinSubtreeModelPath = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath");
            complexSubtreeModelPath = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath");
            if (proteinModelPath == null || complexModelPath == null
                    || proteinSubtreeModelPath == null || complexSubtreeModelPath == null) {
                throw new IllegalArgumentException(
                    "S8 methods require all 4 model paths via system properties");
            }
        }

        // Use full 8-pos highrot confspace (same as GNN training — no RC mapping needed)
        TestKStar.ConfSpaces confSpaces = buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        System.out.println("\nConf space:");
        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.protein.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.complex.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }

        long t0 = System.currentTimeMillis();

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

        // MARKStar settings
        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(maxMut)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism);
        if (useBranch) {
            settingsBuilder.setUseBranchDecomposition(true);
        }
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        markstar.precalcEmats();

        // Load GNN models — no RC mapping needed (same confspace as training)
        GNNConfEnergyCalculator proteinGNN = null, complexGNN = null;
        if (useGNN) {
            // Protein GNN (leaf + subtree) — spread across GPUs to avoid OOM
            proteinGNN = new GNNConfEnergyCalculator(
                    new File(proteinModelPath),
                    markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size(), 0);
            markstar.protein.gnnCalc = proteinGNN;

            GNNSubtreeEnergyCalculator proteinSubtreeGNN = new GNNSubtreeEnergyCalculator(
                    new File(proteinSubtreeModelPath),
                    markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size(), 1);
            markstar.protein.subtreeGnnCalc = proteinSubtreeGNN;

            // Complex GNN (leaf + subtree)
            GNNConfEnergyCalculator complexGNN2 = new GNNConfEnergyCalculator(
                    new File(complexModelPath),
                    markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size(), 2);
            markstar.complex.gnnCalc = complexGNN2;

            GNNSubtreeEnergyCalculator complexSubtreeGNN = new GNNSubtreeEnergyCalculator(
                    new File(complexSubtreeModelPath),
                    markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size(), 3);
            markstar.complex.subtreeGnnCalc = complexSubtreeGNN;

            // Strategy 8 settings
            markstar.protein.useStrategy8 = true;
            markstar.complex.useStrategy8 = true;
            markstar.protein.s7GPUBatchSize = gpuBatch;
            markstar.complex.s7GPUBatchSize = gpuBatch;

            System.out.println("  GNN loaded: protein=" + confSpaces.protein.positions.size()
                    + " pos, complex=" + confSpaces.complex.positions.size() + " pos");
            System.out.println("  Strategy8: gpuBatchSize=" + gpuBatch + ", no RC mapping (same confspace)");
        } else {
            System.out.println("  Running with CCD (no GNN)");
        }

        List<MARKStar.ScoredSequence> sequences = markstar.run();
        long elapsed = System.currentTimeMillis() - t0;

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + method);
        System.out.println("==============================================");
        System.out.println("Total runtime: " + String.format("%.1f", elapsed / 1000.0) + " s");
        System.out.println("Sequences found: " + sequences.size());

        sequences.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
        System.out.println(String.format("%-5s %-50s %12s %12s", "Rank", "Sequence", "K*(lb)", "K*(ub)"));
        System.out.println("-----------------------------------------------------------------------");
        for (int i = 0; i < Math.min(30, sequences.size()); i++) {
            MARKStar.ScoredSequence s = sequences.get(i);
            System.out.println(String.format("%-5d %-50s %12.4f %12.4f",
                    i + 1,
                    s.sequence.toString(Sequence.Renderer.ResType),
                    s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
        }
        System.out.println("==============================================");

        // Close GNN resources
        try {
            if (proteinGNN != null) proteinGNN.close();
        } catch (Exception e) {
            System.err.println("Warning: error closing GNN: " + e.getMessage());
        }
    }

    // ========== 12-pos Scaling Benchmark (CCD only): BranchMARK* vs MARK* ==========

    // ========== Scaling ConfSpace builders for RECOMB submission ==========

    /**
     * Build a parameterized scaling confspace: 8 mutable highrot positions
     * (identical to buildHighRotamerConfSpace8pos) + a configurable list of
     * flexible-only WT positions. Used to construct the n in {16,20,24} scaling
     * series for the BranchMARK* paper.
     *
     * @param flexResidues list of residue numbers (e.g., "A157") to add as
     *                     flexible-only (WT rotamers, continuous CCD).
     */
    private TestKStar.ConfSpaces buildScalingConfSpace(String[] flexResidues) {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();

        // 8 mutable highrot positions (identical to buildHighRotamerConfSpace8pos)
        protein.flexibility.get("A156").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG", "LYS", "MET", "GLU", "ILE").addWildTypeRotamers().setContinuous();

        // Flexible-only WT positions (no mutations)
        for (String res : flexResidues) {
            protein.flexibility.get(res).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /** 16-pos scaling confspace: 8 mutable + 8 flexible WT. */
    private TestKStar.ConfSpaces buildScalingConfSpace16pos() {
        return buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170"});
    }

    /** 20-pos scaling confspace: 8 mutable + 12 flexible WT. */
    private TestKStar.ConfSpaces buildScalingConfSpace20pos() {
        return buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170",
                "A171", "A173", "A175", "A176"});
    }

    /** 24-pos scaling confspace: 8 mutable + 16 flexible WT. */
    private TestKStar.ConfSpaces buildScalingConfSpace24pos() {
        return buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170",
                "A171", "A173", "A175", "A176", "A177", "A178", "A180", "A184"});
    }

    // ========== Generic scaling benchmark runner (any n) ==========

    /**
     * Generic scaling benchmark: runs MARK* or BranchMARK* on a given confspace
     * and reports total runtime, sequences scored, and K* table. Driven by
     * system properties identical to benchmarkHybrid12pos().
     *
     * @param confSpaces  protein/ligand/complex confspaces (any size)
     * @param numPos      position count, for log labeling
     */
    private void runScalingBenchmark(TestKStar.ConfSpaces confSpaces, int numPos) throws Exception {
        String method = System.getProperty("osprey.hybrid.method", "markstar_ccd");
        double epsilon = Double.parseDouble(System.getProperty("osprey.hybrid.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.hybrid.numCPUs", NUM_CPUs);
        int maxMut = Integer.getInteger("osprey.hybrid.maxMut", 1);

        boolean useBranch = method.startsWith("branch_");
        if (method.endsWith("_s8")) {
            throw new IllegalArgumentException(
                    "Scaling benchmark is CCD-only (no GNN model). Use markstar_ccd or branch_ccd.");
        }

        System.out.println("==============================================");
        System.out.println("  " + numPos + "-pos Scaling Benchmark (CCD only)");
        System.out.println("  method=" + method + ", maxMut=" + maxMut);
        System.out.println("  epsilon=" + epsilon + ", cpus=" + cpus);
        System.out.println("==============================================");

        Parallelism parallelism = Parallelism.makeCpu(cpus);

        System.out.println("\nConf space:");
        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.protein.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());

        long t0 = System.currentTimeMillis();

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

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(maxMut)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism);
        if (useBranch) {
            settingsBuilder.setUseBranchDecomposition(true);
        }
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        markstar.precalcEmats();
        List<MARKStar.ScoredSequence> sequences = markstar.run();
        long elapsed = System.currentTimeMillis() - t0;

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + method + "  (" + numPos + " positions)");
        System.out.println("==============================================");
        System.out.println("Total runtime: " + String.format("%.1f", elapsed / 1000.0) + " s");
        System.out.println("Sequences found: " + sequences.size());

        sequences.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
        System.out.println(String.format("%-5s %-50s %12s %12s", "Rank", "Sequence", "K*(lb)", "K*(ub)"));
        System.out.println("-----------------------------------------------------------------------");
        for (int i = 0; i < Math.min(30, sequences.size()); i++) {
            MARKStar.ScoredSequence s = sequences.get(i);
            System.out.println(String.format("%-5d %-50s %12.4f %12.4f",
                    i + 1, s.sequence.toString(Sequence.Renderer.ResType),
                    s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
        }
        System.out.println("==============================================");
    }

    /** 16-position scaling benchmark for the RECOMB submission. */
    public void benchmarkHybrid16pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace16pos(), 16);
    }

    /** 20-position scaling benchmark for the RECOMB submission. */
    public void benchmarkHybrid20pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace20pos(), 20);
    }

    /** 24-position scaling benchmark for the RECOMB submission. */
    public void benchmarkHybrid24pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace24pos(), 24);
    }

    /**
     * Mirror of benchmarkHybrid8pos(), but on the 12-pos confspace
     * (8 mutable highrot + 16 flexible WT) and without any GNN — pure CCD.
     *
     * Sequence count is identical to the 8-pos benchmark (34 seqs at maxMut=1), so the
     * scaling axis isolated by this benchmark is conformational complexity (8 → 12 positions).
     *
     * Methods (osprey.hybrid.method):
     *   markstar_ccd  — MARK* + CCD baseline (12 positions per expansion)
     *   branch_ccd    — BranchMARK* + CCD (edge-level expansion via branch decomposition)
     *
     * Other system properties match benchmarkHybrid8pos exactly:
     *   osprey.hybrid.epsilon (default 0.683), osprey.hybrid.numCPUs (default 16),
     *   osprey.hybrid.maxMut (default 1).
     */
    public void benchmarkHybrid12pos() throws Exception {
        String method = System.getProperty("osprey.hybrid.method", "markstar_ccd");
        double epsilon = Double.parseDouble(System.getProperty("osprey.hybrid.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.hybrid.numCPUs", NUM_CPUs);
        int maxMut = Integer.getInteger("osprey.hybrid.maxMut", 1);

        boolean useBranch = method.startsWith("branch_");
        if (method.endsWith("_s8")) {
            throw new IllegalArgumentException(
                "12-pos benchmark is CCD-only (no GNN model trained on 12-pos confspace). " +
                "Use markstar_ccd or branch_ccd.");
        }

        System.out.println("==============================================");
        System.out.println("  12-pos Scaling Benchmark (CCD only)");
        System.out.println("  method=" + method + ", maxMut=" + maxMut);
        System.out.println("  epsilon=" + epsilon + ", cpus=" + cpus);
        System.out.println("==============================================");

        TestKStar.ConfSpaces confSpaces = buildHighRotamerConfSpace12pos();
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        System.out.println("\nConf space:");
        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.protein.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.complex.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }

        long t0 = System.currentTimeMillis();

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

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(maxMut)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism);
        if (useBranch) {
            settingsBuilder.setUseBranchDecomposition(true);
        }
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        markstar.precalcEmats();

        System.out.println("  Running with CCD (no GNN, 12-pos scaling test)");

        List<MARKStar.ScoredSequence> sequences = markstar.run();
        long elapsed = System.currentTimeMillis() - t0;

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + method + "  (12 positions)");
        System.out.println("==============================================");
        System.out.println("Total runtime: " + String.format("%.1f", elapsed / 1000.0) + " s");
        System.out.println("Sequences found: " + sequences.size());

        sequences.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
        System.out.println(String.format("%-5s %-50s %12s %12s", "Rank", "Sequence", "K*(lb)", "K*(ub)"));
        System.out.println("-----------------------------------------------------------------------");
        for (int i = 0; i < Math.min(30, sequences.size()); i++) {
            MARKStar.ScoredSequence s = sequences.get(i);
            System.out.println(String.format("%-5d %-50s %12.4f %12.4f",
                    i + 1,
                    s.sequence.toString(Sequence.Renderer.ResType),
                    s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
        }
        System.out.println("==============================================");
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
            case "exportHighrot8pos":
                new TestBranchMARKStar().exportGNNDataHighrot8pos();
                break;
            case "patchRigidHighrot8pos":
                new TestBranchMARKStar().patchRigidEmatHighrot8pos();
                break;
            case "bbk_gnn_bench":
                new TestBranchMARKStar().benchmarkBBKStarGNN8pos();
                break;
            case "hybrid_bench":
                new TestBranchMARKStar().benchmarkHybrid8pos();
                break;
            case "hybrid_bench_12pos":
                new TestBranchMARKStar().benchmarkHybrid12pos();
                break;
            case "hybrid_bench_16pos":
                new TestBranchMARKStar().benchmarkHybrid16pos();
                break;
            case "hybrid_bench_20pos":
                new TestBranchMARKStar().benchmarkHybrid20pos();
                break;
            case "hybrid_bench_24pos":
                new TestBranchMARKStar().benchmarkHybrid24pos();
                break;
            case "benchmark":
            default:
                new TestBranchMARKStar().benchmarkGNNStrategies();
                break;
        }
    }

    // ========== BBK* + GNN Benchmark: 8-pos highrot ==========

    /**
     * Benchmark 4 methods on 8-position high-rotamer conf space (390,625 sequences):
     *   1. MARK* + CCD       (baseline)
     *   2. MARK* + GNN S8    (leaf + subtree GNN)
     *   3. BranchMARK* + CCD
     *   4. BranchMARK* + GNN S8
     *
     * All methods use BBK* for sequence selection/pruning.
     *
     * Parameterized via system property osprey.bbkgnn.method:
     *   markstar_ccd, markstar_s8, branch_ccd, branch_s8
     *
     * System properties:
     *   osprey.bbkgnn.method             - which method to run (default: markstar_ccd)
     *   osprey.bbkgnn.numSeqs            - number of best sequences to find (default: 10)
     *   osprey.bbkgnn.epsilon            - partition function epsilon (default: 0.683)
     *   osprey.bbkgnn.numCPUs            - number of CPUs (default: 16)
     *   osprey.bbkgnn.maxMut             - max simultaneous mutations (default: 8)
     *   osprey.gnn.eval.proteinModelPath - protein GNN ONNX model
     *   osprey.gnn.eval.complexModelPath - complex GNN ONNX model
     *   osprey.gnn.eval.proteinSubtreeModelPath - protein subtree GNN ONNX model
     *   osprey.gnn.eval.complexSubtreeModelPath - complex subtree GNN ONNX model
     *   osprey.gnn.gpuBatchSize          - GNN batch size (default: 1000)
     */
    @Test
    public void benchmarkBBKStarGNN8pos() throws Exception {
        String method = System.getProperty("osprey.bbkgnn.method", "markstar_ccd");
        int numSeqs = Integer.getInteger("osprey.bbkgnn.numSeqs", 10);
        double epsilon = Double.parseDouble(System.getProperty("osprey.bbkgnn.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.bbkgnn.numCPUs", NUM_CPUs);
        int maxMut = Integer.getInteger("osprey.bbkgnn.maxMut", 8);

        System.out.println("==============================================");
        System.out.println("  BBK* + GNN Benchmark: 8-pos highrot");
        System.out.println("  method=" + method + ", numSeqs=" + numSeqs);
        System.out.println("  epsilon=" + epsilon + ", cpus=" + cpus + ", maxMut=" + maxMut);
        System.out.println("==============================================");

        boolean useBranch = method.startsWith("branch_");
        boolean useGNN = method.endsWith("_s8");

        // GNN model paths (only needed for S8 methods)
        String proteinModelPath = null, complexModelPath = null;
        String proteinSubtreeModelPath = null, complexSubtreeModelPath = null;
        int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);

        if (useGNN) {
            proteinModelPath = System.getProperty("osprey.gnn.eval.proteinModelPath");
            complexModelPath = System.getProperty("osprey.gnn.eval.complexModelPath");
            proteinSubtreeModelPath = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath");
            complexSubtreeModelPath = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath");
            if (proteinModelPath == null || complexModelPath == null
                    || proteinSubtreeModelPath == null || complexSubtreeModelPath == null) {
                throw new IllegalArgumentException(
                    "S8 methods require all 4 model paths: "
                    + "osprey.gnn.eval.{protein,complex}ModelPath and "
                    + "osprey.gnn.eval.{protein,complex}SubtreeModelPath");
            }
        }

        // Build conf space
        TestKStar.ConfSpaces confSpaces = buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        System.out.println("\nConf space:");
        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.protein.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.complex.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }

        long t0 = System.currentTimeMillis();

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            String label = method;

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

            KStar.Settings kstarSettings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(maxMut)
                    .addScoreConsoleWriter(formatter)
                    .build();
            BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
                    .setNumBestSequences(numSeqs)
                    .setNumConfsPerBatch(8)
                    .build();

            BBKStar bbkstar = new BBKStar(
                    confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                    kstarSettings, bbkstarSettings);

            // Load GNN models (shared across all pfunc instances)
            GNNConfEnergyCalculator proteinGNN = null, complexGNN = null;
            GNNSubtreeEnergyCalculator proteinSubtreeGNN = null, complexSubtreeGNN = null;

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

                // Rigid energy matrix + A* factory
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

                // Determine if this is protein or complex (for GNN model selection)
                boolean isProtein = (confSpace == confSpaces.protein);
                boolean isComplex = (confSpace == confSpaces.complex);

                // Load GNN models for this conf space
                GNNConfEnergyCalculator gnn = null;
                GNNSubtreeEnergyCalculator subtreeGnn = null;
                if (useGNN && (isProtein || isComplex)) {
                    String modelPath = isProtein ? proteinModelPath : complexModelPath;
                    String subtreePath = isProtein ? proteinSubtreeModelPath : complexSubtreeModelPath;
                    int leafGpu = isProtein ? 0 : 1;
                    int subtreeGpu = isProtein ? 2 : 3;

                    gnn = new GNNConfEnergyCalculator(
                            new File(modelPath), ematMinimized, confSpace.positions.size(), leafGpu);
                    subtreeGnn = new GNNSubtreeEnergyCalculator(
                            new File(subtreePath), ematMinimized, confSpace.positions.size(), subtreeGpu);

                    if (isProtein) { proteinGNN = gnn; proteinSubtreeGNN = subtreeGnn; }
                    else { complexGNN = gnn; complexSubtreeGNN = subtreeGnn; }

                    System.out.println("  GNN loaded for " + (isProtein ? "protein" : "complex")
                            + ": " + confSpace.positions.size() + " pos, model=" + modelPath);
                }

                // Capture for lambda
                final GNNConfEnergyCalculator finalGnn = gnn;
                final GNNSubtreeEnergyCalculator finalSubtreeGnn = subtreeGnn;

                // Partition function factory
                if (useBranch && useGNN) {
                    // BranchMARK* + S8
                    info.pfuncFactory = rcs -> {
                        BranchMARKStarBound pfunc = new BranchMARKStarBound(
                                confSpace, ematRigid, ematMinimized,
                                info.confEcalcMinimized, rcs, ecalcMinimized.parallelism);
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        if (finalGnn != null) {
                            pfunc.setGNNBatchCalculator(finalGnn);
                            pfunc.setGPUBatchSize(gpuBatch);
                            pfunc.setCPParams(0.001, 0.10, 0.06);
                        }
                        if (finalSubtreeGnn != null) {
                            pfunc.setSubtreeGNN(finalSubtreeGnn);
                        }
                        pfunc.init(epsilon);
                        return pfunc;
                    };
                    // BranchBBK* scoring
                    info.branchDPScorer = new BranchDPScorer(
                            confSpace, ematRigid, ematMinimized, 8.0, 0.1);
                } else if (useBranch) {
                    // BranchMARK* + CCD
                    info.pfuncFactory = rcs -> {
                        BranchMARKStarBound pfunc = new BranchMARKStarBound(
                                confSpace, ematRigid, ematMinimized,
                                info.confEcalcMinimized, rcs, ecalcMinimized.parallelism);
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        pfunc.init(epsilon);
                        return pfunc;
                    };
                    info.branchDPScorer = new BranchDPScorer(
                            confSpace, ematRigid, ematMinimized, 8.0, 0.1);
                } else if (useGNN) {
                    // MARK* + S8
                    info.pfuncFactory = rcs -> {
                        edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNS8 pfunc =
                                new edu.duke.cs.osprey.markstar.framework.MARKStarBoundGNNS8(
                                        confSpace, ematRigid, ematMinimized,
                                        info.confEcalcMinimized, rcs, ecalcMinimized.parallelism);
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        if (finalGnn != null) {
                            pfunc.setGNNBatchCalculator(finalGnn);
                            pfunc.setGPUBatchSize(gpuBatch);
                            pfunc.setCPParams(0.001, 0.10, 0.06);
                        }
                        if (finalSubtreeGnn != null) {
                            pfunc.setSubtreeGNN(finalSubtreeGnn);
                        }
                        pfunc.init(epsilon);
                        return pfunc;
                    };
                } else {
                    // MARK* + CCD
                    info.pfuncFactory = rcs -> {
                        MARKStarBoundFastQueues pfunc = new MARKStarBoundFastQueues(
                                confSpace, ematRigid, ematMinimized,
                                info.confEcalcMinimized, rcs, ecalcMinimized.parallelism);
                        pfunc.setCorrections(new UpdatingEnergyMatrix(confSpace, ematMinimized, info.confEcalcMinimized));
                        pfunc.init(epsilon);
                        return pfunc;
                    };
                }

                info.confDBFile = null;
            }

            // Run BBK*
            List<KStar.ScoredSequence> sequences = bbkstar.run(ecalcMinimized.tasks);
            long elapsed = System.currentTimeMillis() - t0;

            // Print results
            System.out.println("\n==============================================");
            System.out.println("  RESULTS: " + method);
            System.out.println("==============================================");
            System.out.println("Total runtime: " + String.format("%.1f", elapsed / 1000.0) + " s");
            System.out.println("Sequences found: " + sequences.size());
            System.out.println();

            // Sort by K* score
            sequences.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
            System.out.println(String.format("%-5s %-40s %12s %12s", "Rank", "Sequence", "K*(lb)", "K*(ub)"));
            System.out.println("-----------------------------------------------------------------------");
            for (int i = 0; i < Math.min(30, sequences.size()); i++) {
                KStar.ScoredSequence s = sequences.get(i);
                System.out.println(String.format("%-5d %-40s %12.4f %12.4f",
                        i + 1,
                        s.sequence.toString(Sequence.Renderer.ResType),
                        s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
            }
            System.out.println("==============================================");

            // Close GNN resources
            try {
                if (proteinGNN != null) proteinGNN.close();
                if (complexGNN != null) complexGNN.close();
            } catch (Exception e) {
                System.err.println("Warning: error closing GNN: " + e.getMessage());
            }
        }
    }

    // ========================================================================
    // EGFR 5pos all20 GNN data export
    // ========================================================================

    /**
     * Build confspace for EGFR with a small-molecule ligand.
     * Loads GAFF parameters (prepi + frcmod) for the ligand from prepped files.
     *
     * System properties:
     *   osprey.gnn.egfr.proteinPdb    — protein-only PDB
     *   osprey.gnn.egfr.ligandPdb     — ligand-only PDB (HETATM)
     *   osprey.gnn.egfr.prepFile      — GAFF .prepi file
     *   osprey.gnn.egfr.frcmodFile    — GAFF .frcmod file
     *   osprey.gnn.egfr.tcFile        — template coords (.tc)
     *   osprey.gnn.egfr.rotFile       — rotatable dihedrals (.rot)
     *   osprey.gnn.egfr.ligandResname — 3-letter ligand resname (e.g., AQ4)
     *   osprey.gnn.egfr.hotspotResNums — comma-separated PDB residue numbers
     *   osprey.gnn.egfr.chain         — chain ID (default: A)
     *   osprey.gnn.egfr.outputDir     — output base dir (default: gnn_data/egfr_5pos)
     *   osprey.gnn.numSamples         — max samples (default: 200000)
     */
    @Test
    public void exportGNNDataEGFR5pos() throws Exception {
        String proteinPdb = System.getProperty("osprey.gnn.egfr.proteinPdb",
                "gnn_data/egfr_erlotinib/protein.pdb");
        String ligandPdb = System.getProperty("osprey.gnn.egfr.ligandPdb",
                "gnn_data/egfr_erlotinib/ligand.h.pdb");
        String prepiPath = System.getProperty("osprey.gnn.egfr.prepFile",
                "gnn_data/egfr_erlotinib/ligand.prepi");
        String frcmodPath = System.getProperty("osprey.gnn.egfr.frcmodFile",
                "gnn_data/egfr_erlotinib/ligand.frcmod");
        String tcPath = System.getProperty("osprey.gnn.egfr.tcFile",
                "gnn_data/egfr_erlotinib/ligand.tc");
        String rotPath = System.getProperty("osprey.gnn.egfr.rotFile",
                "gnn_data/egfr_erlotinib/ligand.rot");
        String ligandResname = System.getProperty("osprey.gnn.egfr.ligandResname", "AQ4");
        String hotspotStr = System.getProperty("osprey.gnn.egfr.hotspotResNums", "694,766,768,772,773");
        String chain = System.getProperty("osprey.gnn.egfr.chain", "A");
        String baseDir = System.getProperty("osprey.gnn.egfr.outputDir", "gnn_data/egfr_5pos");
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", NUM_CPUs);

        int[] hotspotResNums = Arrays.stream(hotspotStr.split(","))
                .mapToInt(s -> Integer.parseInt(s.trim())).toArray();

        System.out.println("=== EGFR 5pos GNN Data Export ===");
        System.out.println("Protein PDB: " + proteinPdb);
        System.out.println("Ligand PDB: " + ligandPdb);
        System.out.println("Hotspot residues: " + hotspotStr);
        System.out.println("Ligand resname: " + ligandResname);
        System.out.println("Output: " + baseDir);

        File proteinPdbFile = new File(proteinPdb);
        File ligandPdbFile = new File(ligandPdb);
        File prepiFile = new File(prepiPath);
        File frcmodFile = new File(frcmodPath);
        File tcFile = new File(tcPath);
        File rotFile = new File(rotPath);

        TestKStar.ConfSpaces confSpaces = buildConfSpaceEGFR(
                proteinPdbFile, ligandPdbFile,
                prepiFile, frcmodFile, tcFile, rotFile,
                ligandResname, chain, hotspotResNums);

        System.out.println("\nProtein positions: " + confSpaces.protein.positions.size());
        for (var pos : confSpaces.protein.positions) {
            System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }
        System.out.println("Complex positions: " + confSpaces.complex.positions.size());
        for (var pos : confSpaces.complex.positions) {
            System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }

        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            // --- Export protein data ---
            System.out.println("\n========== Exporting PROTEIN data (EGFR 5pos) ==========");
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
            new GNNDataExporter(proteinConfEcalc, proteinEmatMin, proteinEmatRigid, proteinIG, proteinRCs)
                    .export(numSamples, new File(baseDir + "/protein"));

            // --- Export complex data ---
            System.out.println("\n========== Exporting COMPLEX data (EGFR 5pos) ==========");
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
            new GNNDataExporter(complexConfEcalc, complexEmatMin, complexEmatRigid, complexIG, complexRCs)
                    .export(numSamples, new File(baseDir + "/complex"));

            System.out.println("\n========== EGFR 5pos GNN data export complete ==========");
        }
    }

    /**
     * Build a TestKStar.ConfSpaces for an EGFR system with small-molecule ligand.
     * The ligand gets GAFF parameters from prepi+frcmod files.
     * Specified hotspot residues are made flexible+mutable to all 20 AA.
     */
    private TestKStar.ConfSpaces buildConfSpaceEGFR(
            File proteinPdb, File ligandPdb,
            File prepiFile, File frcmodFile, File tcFile, File rotFile,
            String ligandResname, String chain, int[] hotspotResNums
    ) throws Exception {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

        // Forcefield: AMBER ff14SB + GAFF frcmod for ligand
        ForcefieldParams.Forcefield ff = ForcefieldParams.Forcefield.AMBER;
        ForcefieldFileParser parser;
        if (frcmodFile.exists()) {
            parser = new ForcefieldFileParser(
                    ForcefieldParams.class.getResourceAsStream(ff.paramsPath),
                    frcmodFile.toPath()
            );
        } else {
            parser = new ForcefieldFileParser(
                    ForcefieldParams.class.getResourceAsStream(ff.paramsPath)
            );
        }
        confSpaces.ffparams = new ForcefieldParams(ff, parser);

        // Load GAFF template files (prepi + tc + rot)
        String prepiContent = Files.readString(prepiFile.toPath());
        String tcContent = tcFile.exists() ? Files.readString(tcFile.toPath()) : "";
        String rotContent = rotFile.exists() ? Files.readString(rotFile.toPath()) : "";

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addTemplates(prepiContent)
                .addTemplateCoords(tcContent)
                .addRotamers(rotContent)
                .build();

        // Load PDBs
        Molecule proteinMol = PDBIO.readFile(proteinPdb);
        Molecule ligandMol = PDBIO.readFile(ligandPdb);

        // Find residue range for protein strand
        int firstRes = Integer.MAX_VALUE, lastRes = Integer.MIN_VALUE;
        for (Residue res : proteinMol.residues) {
            String resNum = Residues.normalizeResNum(res.getPDBResNumber());
            // Strip leading chain ID (e.g. "A672" -> "672")
            int i = 0;
            while (i < resNum.length() && !Character.isDigit(resNum.charAt(i))) i++;
            int rn = Integer.parseInt(resNum.substring(i).replaceAll("[^0-9-].*$", ""));
            if (rn < firstRes) firstRes = rn;
            if (rn > lastRes) lastRes = rn;
        }
        System.out.println("Protein residues: " + chain + firstRes + " to " + chain + lastRes);

        // Build protein strand
        Strand protein = new Strand.Builder(proteinMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(chain + firstRes, chain + lastRes)
                .build();

        // Make hotspot residues flexible + mutable to all 20 AA
        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
                          "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
        for (int resNum : hotspotResNums) {
            String id = chain + resNum;
            Strand.ResidueFlex rf = protein.flexibility.get(id);
            if (rf != null) {
                rf.setLibraryRotamers(all20)
                        .addWildTypeRotamers()
                        .setContinuous();
                System.out.println("  Flexible: " + id + " (wt=" + rf.wildType + ")");
            } else {
                System.out.println("  WARNING: " + id + " not found in protein flexibility");
            }
        }

        // Build ligand strand (fixed — no flexibility)
        int ligFirstRes = Integer.MAX_VALUE, ligLastRes = Integer.MIN_VALUE;
        String ligChain = chain;
        for (Residue res : ligandMol.residues) {
            try {
                String rnStr = Residues.normalizeResNum(res.getPDBResNumber());
                int i = 0;
                while (i < rnStr.length() && !Character.isDigit(rnStr.charAt(i))) i++;
                int rn = Integer.parseInt(rnStr.substring(i).replaceAll("[^0-9-].*$", ""));
                if (rn < ligFirstRes) { ligFirstRes = rn; ligChain = String.valueOf(res.getChainId()); }
                if (rn > ligLastRes) ligLastRes = rn;
            } catch (NumberFormatException e) {
                // skip non-numeric residue IDs (e.g., waters)
            }
        }
        Strand ligand = new Strand.Builder(ligandMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(ligChain + ligFirstRes, ligChain + ligLastRes)
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }
}
