package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.List;

/**
 * Scaling benchmark series for the BranchMARK* paper.
 *
 * Compares MARK* vs BranchMARK* (both ± GNN-S8) as conformational complexity grows
 * from 8 → 12 → 16 → 20 → 24 positions while sequence count stays bounded by
 * {@code osprey.scaling.maxMut} (default 1, i.e. 34 sequences for 8 mutables).
 *
 * Entry points (callable from {@link TestBranchMARKStar} main dispatch):
 *   - {@link #benchmarkScaling8pos()}  — full 8-pos highrot, supports CCD or S8
 *   - {@link #benchmarkScaling12pos()} — 8 mut + 4 flex WT, CCD-only
 *   - {@link #benchmarkScaling16pos()} / 20 / 24 — 8 mut + N-8 flex WT, CCD-only
 *
 * System properties:
 *   osprey.scaling.method  — markstar_ccd | branch_ccd | markstar_s8 | branch_s8 (8pos only)
 *   osprey.scaling.epsilon — partition function epsilon (default 0.683)
 *   osprey.scaling.numCPUs — number of CPUs (default {@link TestBranchMARKStar#NUM_CPUs})
 *   osprey.scaling.maxMut  — max simultaneous mutations (default 1)
 *   For S8: osprey.gnn.eval.{protein,complex}{,Subtree}ModelPath + osprey.gnn.gpuBatchSize
 */
public class ScalingBench {

    // ========== Scaling 8-pos benchmark (full GNN-capable) ==========

    public void benchmarkScaling8pos() throws Exception {
        String method = System.getProperty("osprey.scaling.method", "markstar_ccd");
        double epsilon = Double.parseDouble(System.getProperty("osprey.scaling.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.scaling.numCPUs", TestBranchMARKStar.NUM_CPUs);
        int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
        int maxMut = Integer.getInteger("osprey.scaling.maxMut", 1);

        boolean useBranch = method.startsWith("branch_");
        boolean useGNN = method.endsWith("_s8");

        System.out.println("==============================================");
        System.out.println("  8-pos Highrot Benchmark (full confspace)");
        System.out.println("  method=" + method + ", maxMut=" + maxMut);
        System.out.println("  epsilon=" + epsilon + ", cpus=" + cpus);
        System.out.println("==============================================");

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

        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        printConfSpaceSummary(confSpaces);
        long t0 = System.currentTimeMillis();

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) ->
                new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                .build().calcReferenceEnergies())
                        .build();

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon).setMaxSimultaneousMutations(maxMut)
                .setShowPfuncProgress(true).setParallelism(parallelism);
        if (useBranch) settingsBuilder.setUseBranchDecomposition(true);
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
        markstar.precalcEmats();

        GNNConfEnergyCalculator proteinGNN = null;
        if (useGNN) {
            proteinGNN = new GNNConfEnergyCalculator(
                    new File(proteinModelPath), markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size(), 0);
            markstar.protein.gnnCalc = proteinGNN;
            markstar.protein.subtreeGnnCalc = new GNNSubtreeEnergyCalculator(
                    new File(proteinSubtreeModelPath), markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size(), 1);

            markstar.complex.gnnCalc = new GNNConfEnergyCalculator(
                    new File(complexModelPath), markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size(), 2);
            markstar.complex.subtreeGnnCalc = new GNNSubtreeEnergyCalculator(
                    new File(complexSubtreeModelPath), markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size(), 3);

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

        printResults(method, 8, elapsed, sequences);

        try { if (proteinGNN != null) proteinGNN.close(); }
        catch (Exception e) { System.err.println("Warning: error closing GNN: " + e.getMessage()); }
    }

    // ========== 12-pos benchmark (CCD-only, custom confspace) ==========

    public void benchmarkScaling12pos() throws Exception {
        String method = System.getProperty("osprey.scaling.method", "markstar_ccd");
        if (method.endsWith("_s8")) {
            throw new IllegalArgumentException(
                "12-pos benchmark is CCD-only (no GNN model trained on 12-pos confspace).");
        }
        runScalingBenchmark(ConfSpaces2RL0.buildHighRotamerConfSpace12pos(), 12);
    }

    // ========== 16 / 20 / 24-pos scaling benchmarks (CCD-only) ==========

    public void benchmarkScaling16pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170"}), 16);
    }

    public void benchmarkScaling20pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170",
                "A171", "A173", "A175", "A176"}), 20);
    }

    public void benchmarkScaling24pos() throws Exception {
        runScalingBenchmark(buildScalingConfSpace(new String[]{
                "A157", "A159", "A160", "A161", "A165", "A166", "A167", "A170",
                "A171", "A173", "A175", "A176", "A177", "A178", "A180", "A184"}), 24);
    }

    // ========== Shared helpers ==========

    /**
     * Build a parameterized scaling confspace: 8 mutable highrot positions
     * (identical to {@link ConfSpaces2RL0#buildHighRotamerConfSpace8pos()})
     * plus a configurable list of flexible-only WT positions.
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
        protein.flexibility.get("A156").setLibraryRotamers("ARG","LYS","MET","GLU","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG","LYS","MET","GLU","HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG","LYS","MET","GLU","THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG","LYS","MET","GLU","ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG","LYS","MET","GLU","TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG","LYS","MET","GLU","ILE").addWildTypeRotamers().setContinuous();

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

    private void runScalingBenchmark(TestKStar.ConfSpaces confSpaces, int numPos) throws Exception {
        String method = System.getProperty("osprey.scaling.method", "markstar_ccd");
        double epsilon = Double.parseDouble(System.getProperty("osprey.scaling.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.scaling.numCPUs", TestBranchMARKStar.NUM_CPUs);
        int maxMut = Integer.getInteger("osprey.scaling.maxMut", 1);

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
        printConfSpaceSummary(confSpaces);

        long t0 = System.currentTimeMillis();

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) ->
                new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                .build().calcReferenceEnergies())
                        .build();

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon).setMaxSimultaneousMutations(maxMut)
                .setShowPfuncProgress(true).setParallelism(parallelism);
        if (useBranch) settingsBuilder.setUseBranchDecomposition(true);
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
        markstar.precalcEmats();
        List<MARKStar.ScoredSequence> sequences = markstar.run();
        long elapsed = System.currentTimeMillis() - t0;

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        printResults(method, numPos, elapsed, sequences);
    }

    private static void printConfSpaceSummary(TestKStar.ConfSpaces confSpaces) {
        System.out.println("\nConf space:");
        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        for (SimpleConfSpace.Position pos : confSpaces.protein.positions) {
            System.out.println("    " + pos.resNum + ": " + pos.resTypes.size() + " AA, "
                    + pos.resConfs.size() + " RCs");
        }
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());
    }

    private static void printResults(String method, int numPos, long elapsedMs,
                                     List<MARKStar.ScoredSequence> sequences) {
        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + method + "  (" + numPos + " positions)");
        System.out.println("==============================================");
        System.out.println("Total runtime: " + String.format("%.1f", elapsedMs / 1000.0) + " s");
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
}
