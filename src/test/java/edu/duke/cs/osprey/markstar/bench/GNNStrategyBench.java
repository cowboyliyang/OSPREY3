package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.lute.LUTE;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEIO;
import edu.duke.cs.osprey.lute.LUTEState;
import edu.duke.cs.osprey.lute.ConfSampler;
import edu.duke.cs.osprey.lute.RandomizedDFSConfSampler;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Parameterized GNN-strategy benchmark on the 2RL0 medium / highrot conf spaces.
 *
 * Runs one strategy per invocation; selected by {@code osprey.gnn.benchStrategy}:
 *   ccd, strategy7, strategy8, strategy9, branch_ccd, branch_gnn
 *
 * Driven from {@link edu.duke.cs.osprey.markstar.RunBenchmark} (which is invoked
 * by {@code bench_s8.sh}).
 */
public class GNNStrategyBench {

    private static class GNNBenchResult {
        List<MARKStar.ScoredSequence> scores;
        long gnnConfs;
        long ccdConfs;
    }

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
        String pModel, cModel;

        switch (which) {
            case "ccd":
                label = "MARK*+CCD"; pModel = null; cModel = null; break;
            case "strategy7":
                label = "Strategy7:decoupled"; pModel = proteinModelPath; cModel = complexModelPath; break;
            case "strategy8":
                label = "Strategy8:S7+subtreeGNN"; pModel = proteinModelPath; cModel = complexModelPath; break;
            case "strategy9":
                label = "Strategy9:logZBoundOracle"; pModel = proteinModelPath; cModel = complexModelPath; break;
            case "branch_ccd":
                label = "BranchMARK*+CCD"; pModel = null; cModel = null; break;
            case "branch_gnn":
                label = "BranchMARK*+GNN"; pModel = proteinModelPath; cModel = complexModelPath; break;
            default:
                throw new IllegalArgumentException("Unknown strategy: " + which
                        + ". Use one of: ccd, strategy7, strategy8, strategy9, branch_ccd, branch_gnn");
        }

        long t0 = System.currentTimeMillis();
        GNNBenchResult result = runGNNBenchmark(epsilon, numSeqs, pModel, cModel, label);
        long elapsed = System.currentTimeMillis() - t0;

        System.out.println("\n==============================================");
        System.out.println("  RESULTS: " + label);
        System.out.println("==============================================");
        System.out.println(String.format("%-25s %10s %12s %12s %10s",
                "Method", "Time(s)", "GNN confs", "CCD confs", "#Seqs"));
        System.out.println("----------------------------------------------------------------------");
        printBenchRow(label, elapsed, result);

        List<MARKStar.ScoredSequence> converged = new ArrayList<>();
        int unconvergedCount = 0;
        for (MARKStar.ScoredSequence s : result.scores) {
            double lb = s.score.lowerBoundLog10();
            double ub = s.score.upperBoundLog10();
            if (Double.isInfinite(lb) || Double.isInfinite(ub)) unconvergedCount++;
            else converged.add(s);
        }
        converged.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));

        System.out.println("\n=== K* Scores (top 30 converged by lower bound, descending) ===");
        System.out.println("  Total: " + result.scores.size() + " sequences, " + converged.size()
                + " converged, " + unconvergedCount + " unconverged");
        for (int i = 0; i < Math.min(30, converged.size()); i++) {
            MARKStar.ScoredSequence s = converged.get(i);
            System.out.println(String.format("  Rank %3d: %s  K*=[%.4f, %.4f]",
                    i + 1, s.sequence.toString(Sequence.Renderer.ResType),
                    s.score.lowerBoundLog10(), s.score.upperBoundLog10()));
        }
    }

    private void printBenchRow(String label, long timeMs, GNNBenchResult r) {
        System.out.println(String.format("%-25s %10.1f %12d %12d %10d",
                label, timeMs / 1000.0, r.gnnConfs, r.ccdConfs, r.scores.size()));
    }

    private GNNBenchResult runGNNBenchmark(double epsilon, int numSeqs,
                                           String proteinModelPath, String complexModelPath,
                                           String label) {
        String confSpaceType = System.getProperty("osprey.gnn.confSpace", "medium");
        TestKStar.ConfSpaces confSpaces;
        String ematDefault;
        switch (confSpaceType) {
            case "highrot":
                confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace();
                ematDefault = "emat_cache/highrot_4pos.*.dat";
                break;
            case "medium":
            default:
                confSpaces = ConfSpaces2RL0.buildMediumMutableConfSpace();
                ematDefault = "emat_cache/med10_2pos.*.dat";
                break;
        }
        int cpus = Integer.getInteger("osprey.gnn.numCPUs", TestBranchMARKStar.NUM_CPUs);
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

        String ematCachePattern = System.getProperty("osprey.gnn.ematCache", ematDefault);

        String which = System.getProperty("osprey.gnn.benchStrategy", "ccd");
        boolean useBranch = which.startsWith("branch_");

        MARKStar.Settings.Builder settingsBuilder = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(4)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .setEnergyMatrixCachePattern(ematCachePattern);
        if (useBranch) settingsBuilder.setUseBranchDecomposition(true);
        MARKStar.Settings settings = settingsBuilder.build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        long ematStart = System.currentTimeMillis();
        markstar.precalcEmats();
        long ematTime = System.currentTimeMillis() - ematStart;
        System.out.println("  [" + label + "] Emat time: " + ematTime + " ms (cached: " + ematCachePattern + ")");

        GNNConfEnergyCalculator proteinGNN = null;
        GNNConfEnergyCalculator complexGNN = null;
        if (proteinModelPath != null && complexModelPath != null) {
            proteinGNN = new GNNConfEnergyCalculator(
                    new File(proteinModelPath), markstar.protein.minimizingEmat,
                    confSpaces.protein.positions.size());
            markstar.protein.gnnCalc = proteinGNN;

            complexGNN = new GNNConfEnergyCalculator(
                    new File(complexModelPath), markstar.complex.minimizingEmat,
                    confSpaces.complex.positions.size());
            markstar.complex.gnnCalc = complexGNN;

            String gnnTrainingCS = System.getProperty("osprey.gnn.trainingConfSpace", "");
            if (!gnnTrainingCS.isEmpty() && !gnnTrainingCS.equals(confSpaceType)) {
                System.out.println("  Building RC mapping: inference=" + confSpaceType + " → training=" + gnnTrainingCS);
                TestKStar.ConfSpaces trainingConfSpaces = trainingConfSpacesFor(gnnTrainingCS);
                int[][] proteinMap = GNNConfEnergyCalculator.buildRCMapping(confSpaces.protein, trainingConfSpaces.protein);
                int[][] complexMap = GNNConfEnergyCalculator.buildRCMapping(confSpaces.complex, trainingConfSpaces.complex);
                proteinGNN.setRCMapping(proteinMap);
                complexGNN.setRCMapping(complexMap);
            }

            if (which.equals("strategy8") || which.equals("strategy9")) {
                boolean s8 = which.equals("strategy8");
                if (s8) {
                    markstar.protein.useStrategy8 = true;
                    markstar.complex.useStrategy8 = true;
                } else {
                    markstar.protein.useStrategy9 = true;
                    markstar.complex.useStrategy9 = true;
                }
                int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
                markstar.protein.s7GPUBatchSize = gpuBatch;
                markstar.complex.s7GPUBatchSize = gpuBatch;

                String pSubPath = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/protein/model_subtree/subtree_model.onnx");
                String cSubPath = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath",
                        "gnn_data/2RL0_all20_4pos_merged/complex/model_subtree/subtree_model.onnx");
                GNNSubtreeEnergyCalculator pSubGNN = new GNNSubtreeEnergyCalculator(
                        new File(pSubPath), markstar.protein.minimizingEmat,
                        confSpaces.protein.positions.size());
                GNNSubtreeEnergyCalculator cSubGNN = new GNNSubtreeEnergyCalculator(
                        new File(cSubPath), markstar.complex.minimizingEmat,
                        confSpaces.complex.positions.size());

                if (!gnnTrainingCS.isEmpty() && !gnnTrainingCS.equals(confSpaceType)) {
                    TestKStar.ConfSpaces trainCS = trainingConfSpacesFor(gnnTrainingCS);
                    pSubGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.protein, trainCS.protein));
                    cSubGNN.setRCMapping(
                            GNNConfEnergyCalculator.buildRCMapping(confSpaces.complex, trainCS.complex));
                }

                markstar.protein.subtreeGnnCalc = pSubGNN;
                markstar.complex.subtreeGnnCalc = cSubGNN;
                System.out.println("  " + (s8 ? "Strategy8" : "Strategy9")
                        + ": gpuBatchSize=" + gpuBatch + ", subtree models loaded");
            } else if (which.equals("strategy7") || which.equals("branch_gnn")) {
                markstar.protein.useStrategy7 = true;
                markstar.complex.useStrategy7 = true;
                int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);
                markstar.protein.s7GPUBatchSize = gpuBatch;
                markstar.complex.s7GPUBatchSize = gpuBatch;
                System.out.println("  " + (useBranch ? "BranchMARK*+GNN" : "Strategy7")
                        + ": gpuBatchSize=" + gpuBatch);
            }

            System.out.println("  GNN loaded: protein=" + confSpaces.protein.positions.size()
                    + " pos, complex=" + confSpaces.complex.positions.size() + " pos");
        } else {
            System.out.println("  Running with CCD (no GNN)");
        }

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
        System.out.println("  [" + label + "] Computed " + scores.size() + " sequences");

        try {
            if (proteinGNN != null) proteinGNN.close();
            if (complexGNN != null) complexGNN.close();
        } catch (Exception e) {
            System.err.println("Warning: error closing GNN: " + e.getMessage());
        }
        return result;
    }

    private static TestKStar.ConfSpaces trainingConfSpacesFor(String name) {
        switch (name) {
            case "all20": return ConfSpaces2RL0.buildAllMutableConfSpace();
            case "highrot": return ConfSpaces2RL0.buildHighRotamerConfSpace();
            case "medium": return ConfSpaces2RL0.buildMediumMutableConfSpace();
            default: throw new IllegalArgumentException("Unknown trainingConfSpace: " + name);
        }
    }
}
