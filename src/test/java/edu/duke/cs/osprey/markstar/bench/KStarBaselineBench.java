package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * K* baseline benchmark for comparison with MARK* and GNN S9.
 *
 * Uses the same 2RL0 conf spaces from {@link ConfSpaces2RL0} so results are
 * directly comparable with {@link GNNStrategyBench}.
 *
 * System properties:
 *   osprey.kstar.confSpace  — medium | highrot | highrot8 | all20
 *   osprey.kstar.epsilon    — approximation ratio (default 0.683)
 *   osprey.kstar.numSeqs    — max sequences to enumerate (default 20)
 *   osprey.kstar.numCPUs    — parallelism (default from TestBranchMARKStar)
 *   osprey.kstar.ematCache  — emat cache pattern (default auto by confSpace)
 *   osprey.kstar.outputCSV  — CSV output path (default stdout only)
 */
public class KStarBaselineBench {

    @Test
    public void benchmarkKStar() {
        double epsilon = Double.parseDouble(
                System.getProperty("osprey.kstar.epsilon", "0.683"));
        int numSeqs = Integer.getInteger("osprey.kstar.numSeqs", 20);
        int cpus = Integer.getInteger("osprey.kstar.numCPUs",
                TestBranchMARKStar.NUM_CPUs);
        String confSpaceType = System.getProperty("osprey.kstar.confSpace", "medium");

        System.out.println("==============================================");
        System.out.println("  K* Baseline Benchmark");
        System.out.println("  confSpace=" + confSpaceType + ", epsilon=" + epsilon
                + ", numSeqs=" + numSeqs + ", cpus=" + cpus);
        System.out.println("==============================================");

        TestKStar.ConfSpaces confSpaces;
        String ematDefault;
        switch (confSpaceType) {
            case "highrot":
                confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace();
                ematDefault = "emat_cache/highrot_4pos.*.dat";
                break;
            case "highrot8":
                confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace8pos();
                ematDefault = "emat_cache/highrot_8pos.*.dat";
                break;
            case "all20":
                confSpaces = ConfSpaces2RL0.buildAllMutableConfSpace();
                ematDefault = "emat_cache/all20_4pos.*.dat";
                break;
            case "medium":
            default:
                confSpaces = ConfSpaces2RL0.buildMediumMutableConfSpace();
                ematDefault = "emat_cache/med10_2pos.*.dat";
                break;
        }
        String ematCachePattern = System.getProperty("osprey.kstar.ematCache", ematDefault);

        Parallelism parallelism = Parallelism.makeCpu(cpus);

        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            KStar.Settings settings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(4)
                    .setShowPfuncProgress(true)
                    .build();

            KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand,
                    confSpaces.complex, settings);

            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                SimpleConfSpace cs = (SimpleConfSpace) info.confSpace;

                info.confEcalc = new ConfEnergyCalculator.Builder(cs, ecalc)
                        .setReferenceEnergies(
                                new SimplerEnergyMatrixCalculator.Builder(cs, ecalc)
                                        .build().calcReferenceEnergies())
                        .build();

                EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
                        .setCacheFile(new File(String.format(
                                ematCachePattern.replace("*", "%s"),
                                info.type.name().toLowerCase())))
                        .build()
                        .calcEnergyMatrix();

                info.pfuncFactory = (rcs) -> {
                    GradientDescentPfunc pfunc = new GradientDescentPfunc(
                            info.confEcalc,
                            new ConfAStarTree.Builder(emat, rcs)
                                    .setTraditional().build(),
                            new ConfAStarTree.Builder(emat, rcs)
                                    .setTraditional().build(),
                            rcs.getNumConformations()
                    );
                    return pfunc;
                };
                info.confDBFile = null;
            }

            long t0 = System.currentTimeMillis();
            List<KStar.ScoredSequence> scores = kstar.run(ecalc.tasks);
            long elapsed = System.currentTimeMillis() - t0;

            // Print results
            System.out.println("\n==============================================");
            System.out.println("  K* RESULTS");
            System.out.println("==============================================");
            System.out.println(String.format("  Total time: %.1f s", elapsed / 1000.0));
            System.out.println(String.format("  Sequences: %d", scores.size()));

            List<KStar.ScoredSequence> converged = new ArrayList<>();
            int unconverged = 0;
            for (KStar.ScoredSequence s : scores) {
                if (s.score.lowerBoundLog10() == Double.NEGATIVE_INFINITY
                        || s.score.upperBoundLog10() == Double.POSITIVE_INFINITY) {
                    unconverged++;
                } else {
                    converged.add(s);
                }
            }
            converged.sort((a, b) -> Double.compare(
                    b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));

            System.out.println("  Converged: " + converged.size()
                    + ", Unconverged: " + unconverged);
            System.out.println("\n=== K* Scores (top 30 by lower bound) ===");
            for (int i = 0; i < Math.min(30, converged.size()); i++) {
                KStar.ScoredSequence s = converged.get(i);
                System.out.println(String.format(
                        "  Rank %3d: %s  K*=[%.4f, %.4f]",
                        i + 1,
                        s.sequence.toString(Sequence.Renderer.ResType),
                        s.score.lowerBoundLog10(),
                        s.score.upperBoundLog10()));
            }

            // Write CSV if requested
            String csvPath = System.getProperty("osprey.kstar.outputCSV");
            if (csvPath != null) {
                try (PrintWriter pw = new PrintWriter(new FileWriter(csvPath))) {
                    pw.println("rank,sequence,lb_log10,ub_log10,elapsed_s");
                    for (int i = 0; i < converged.size(); i++) {
                        KStar.ScoredSequence s = converged.get(i);
                        pw.println(String.format("%d,%s,%.6f,%.6f,%.1f",
                                i + 1,
                                s.sequence.toString(Sequence.Renderer.ResType),
                                s.score.lowerBoundLog10(),
                                s.score.upperBoundLog10(),
                                elapsed / 1000.0));
                    }
                } catch (IOException e) {
                    System.err.println("Error writing CSV: " + e.getMessage());
                }
                System.out.println("CSV written to " + csvPath);
            }
        }
    }
}
