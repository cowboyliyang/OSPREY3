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
import edu.duke.cs.osprey.markstar.framework.branch.BranchDPScorer;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.List;

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
}
