package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.packstar.PackStarPartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * PACK* paper Table 2 ("scaling with n") driver.
 *
 * Runs MARK* (deterministic, branch-decomposition DP) and PACK* (PAC two-stage
 * estimator) on the SAME {@link ConfSpaces2RL0#buildWildTypeConfSpace(int)}
 * system already used by {@link TestBranchMARKStar}'s numFlexible harness (the
 * "10-position benchmark" referenced in the PACK* paper's Introduction), so the
 * two methods' epsilon/CCD-count scaling with n is directly comparable at
 * n in {8,10,12,16,20} (SI Theorem 3 / main paper Table 2).
 *
 * This class intentionally duplicates (rather than refactors) the KStar +
 * PackStarPartitionFunction wiring already proven in
 * {@link GenericPDBBench#runPackStar} so the real 38-system benchmark path is
 * left untouched.
 *
 * Entry point: {@link TestBranchMARKStar} main dispatch, mode "scaling_pac".
 *
 * System properties:
 *   branchdp.test.numFlexible   — n (default 10)
 *   osprey.scalingpac.method    — markstar | packstar | both (default both)
 *   osprey.scalingpac.epsilon   — target epsilon (default 0.683)
 *   osprey.scalingpac.numCPUs   — CPU count (default TestBranchMARKStar.NUM_CPUs)
 *   osprey.scalingpac.outputCsv — output CSV path (appended; header written once)
 *   plus every packstar.pac.* / packstar.dp.* property (samples, trainSamples,
 *   pilotSamples, confidence, residualBound, clip, iterate, etaEnabled, dumpDir,
 *   randomSeed, dp.gpu, ...) — these are read directly by PackStarEstimator /
 *   PackStarPartitionFunction from system properties, exactly as in the real
 *   38-system benchmark, so no extra plumbing is needed here.
 */
public class PackStarScalingBench {

    private static final String HEADER =
        "n,method,rank,sequence,target_eps,score_log10,lb_log10,ub_log10," +
        "prot_lb_log10,prot_ub_log10,prot_status,prot_eps,prot_nconf,prot_nscored," +
        "lig_lb_log10,lig_ub_log10,lig_status,lig_eps,lig_nconf,lig_nscored," +
        "comp_lb_log10,comp_ub_log10,comp_status,comp_eps,comp_nconf,comp_nscored," +
        "total_time_s,n_s,eta_enabled,seed";

    public void benchmarkScalingN() throws Exception {
        int numFlexible = Integer.getInteger("branchdp.test.numFlexible", 10);
        double epsilon = Double.parseDouble(System.getProperty("osprey.scalingpac.epsilon", "0.683"));
        String method = System.getProperty("osprey.scalingpac.method", "both");
        String outputCsv = System.getProperty("osprey.scalingpac.outputCsv",
                "/usr/xtmp/lz280/bench_comparison/results/scaling_n.csv");
        int cpus = Integer.getInteger("osprey.scalingpac.numCPUs", TestBranchMARKStar.NUM_CPUs);

        System.out.println("==============================================");
        System.out.println("  PACK* paper Table 2: scaling with n");
        System.out.println("  numFlexible=" + numFlexible + ", epsilon=" + epsilon
                + ", method=" + method + ", cpus=" + cpus);
        System.out.println("  outputCsv=" + outputCsv);
        System.out.println("==============================================");

        File outFile = new File(outputCsv);
        if (outFile.getParentFile() != null) {
            outFile.getParentFile().mkdirs();
        }

        if (method.equals("markstar") || method.equals("both")) {
            runMarkstarLeg(numFlexible, epsilon, outputCsv);
        }
        if (method.equals("packstar") || method.equals("both")) {
            runPackstarLeg(numFlexible, epsilon, cpus, outputCsv);
        }
    }

    private void runMarkstarLeg(int numFlexible, double epsilon, String outputCsv) throws Exception {
        long t0 = System.currentTimeMillis();
        TestBranchMARKStar.MARKStarResult result =
                TestBranchMARKStar.runMARKStarOnly(numFlexible, epsilon, true, "MARKStar-scalingN");
        double elapsedS = (System.currentTimeMillis() - t0) / 1000.0;
        writeMarkstarRows(numFlexible, result.scores, epsilon, elapsedS, outputCsv);
        System.out.println("[scaling_n] MARK* n=" + numFlexible + " done in "
                + String.format("%.1f", elapsedS) + " s, " + result.scores.size() + " sequence(s) scored");
    }

    private void runPackstarLeg(int numFlexible, double epsilon, int cpus, String outputCsv) throws Exception {
        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildWildTypeConfSpace(numFlexible);
        Parallelism parallelism = Parallelism.makeCpu(cpus);

        boolean reduceMinimizations = Boolean.parseBoolean(
                System.getProperty("osprey.packstar.reduceMinimizations", "true"));
        boolean correctionTighteningEnabled = Boolean.parseBoolean(
                System.getProperty("osprey.packstar.correctionTightening", "true"));

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        try {
            KStar.Settings settings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(1)
                    .setShowPfuncProgress(true)
                    .build();
            KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);

            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                SimpleConfSpace cs = (SimpleConfSpace) info.confSpace;
                String stateName = info.type.name();

                ConfEnergyCalculator minimizingConfEcalc = new ConfEnergyCalculator.Builder(cs, minimizingEcalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, minimizingEcalc)
                                .build().calcReferenceEnergies())
                        .build();
                ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(cs, rigidEcalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, rigidEcalc)
                                .build().calcReferenceEnergies())
                        .build();
                info.confEcalc = minimizingConfEcalc;

                EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                        .build().calcEnergyMatrix();
                EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
                        .build().calcEnergyMatrix();
                UpdatingEnergyMatrix corrections = new UpdatingEnergyMatrix(cs, minimizingEmat, minimizingConfEcalc);

                info.pfuncFactory = (rcs) -> {
                    PackStarPartitionFunction pfunc = new PackStarPartitionFunction(
                            cs, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism, stateName);
                    pfunc.setCorrections(corrections);
                    pfunc.setReduceMinimizations(reduceMinimizations);
                    pfunc.setCorrectionTighteningEnabled(correctionTighteningEnabled);
                    return pfunc;
                };
                info.confDBFile = null;
            }

            long t0 = System.currentTimeMillis();
            List<KStar.ScoredSequence> scores = kstar.run(minimizingEcalc.tasks);
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
            double elapsedS = (System.currentTimeMillis() - t0) / 1000.0;
            writePackstarRows(numFlexible, scores, epsilon, elapsedS, outputCsv);
            System.out.println("[scaling_n] PACK* n=" + numFlexible + " done in "
                    + String.format("%.1f", elapsedS) + " s, " + scores.size() + " sequence(s) scored");
        } finally {
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
        }
    }

    private void writeMarkstarRows(int n, List<MARKStar.ScoredSequence> scores,
                                    double targetEps, double totalTimeS, String outputCsv) {
        boolean writeHeader = !new File(outputCsv).exists();
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputCsv, true))) {
            if (writeHeader) pw.println(HEADER);
            scores.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
            for (int i = 0; i < scores.size(); i++) {
                MARKStar.ScoredSequence s = scores.get(i);
                pw.println(formatRow(n, "markstar", i + 1, s.sequence, s.score, targetEps, totalTimeS));
            }
        } catch (IOException e) {
            System.err.println("[scaling_n] CSV write failed: " + e.getMessage());
        }
    }

    private void writePackstarRows(int n, List<KStar.ScoredSequence> scores,
                                    double targetEps, double totalTimeS, String outputCsv) {
        boolean writeHeader = !new File(outputCsv).exists();
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputCsv, true))) {
            if (writeHeader) pw.println(HEADER);
            scores.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
            for (int i = 0; i < scores.size(); i++) {
                KStar.ScoredSequence s = scores.get(i);
                pw.println(formatRow(n, "packstar", i + 1, s.sequence, s.score, targetEps, totalTimeS));
            }
        } catch (IOException e) {
            System.err.println("[scaling_n] CSV write failed: " + e.getMessage());
        }
    }

    /**
     * Mirrors GenericPDBBench#formatScoreRow, prefixed with the scaling variable n and
     * suffixed with the sweep identifiers (n_s, eta_enabled, seed) so that Table 2
     * (scaling with n), Fig 1 (sample-size convergence), Table 3 (no-eta ablation) and
     * the PAC coverage runs can all append to a shared or per-sweep CSV and still be
     * disambiguated by column instead of by filename/job bookkeeping alone.
     */
    private static String formatRow(int n, String method, int rank, Sequence sequence,
                                      KStarScore score, double targetEps, double totalTimeS) {
        Double scoreLog = score.scoreLog10();
        String scoreStr = (scoreLog == null || scoreLog.isNaN()) ? "" : String.format("%.6f", scoreLog);
        Double lbLog = score.lowerBoundLog10();
        Double ubLog = score.upperBoundLog10();
        String lbStr = (lbLog == null || lbLog.isNaN()) ? "" : String.format("%.6f", lbLog);
        String ubStr = (ubLog == null || ubLog.isNaN()) ? "" : String.format("%.6f", ubLog);
        String nSamples = System.getProperty("packstar.pac.samples", "");
        String etaEnabled = System.getProperty("packstar.pac.etaEnabled", "true");
        String seed = System.getProperty("packstar.pac.randomSeed", "");
        return String.format("%d,%s,%d,%s,%.6f,%s,%s,%s,%s,%s,%s,%.1f,%s,%s,%s",
                n, method, rank,
                sequence.toString(Sequence.Renderer.ResType),
                targetEps,
                scoreStr, lbStr, ubStr,
                formatPfunc(score.protein),
                formatPfunc(score.ligand),
                formatPfunc(score.complex),
                totalTimeS,
                nSamples, etaEnabled, seed);
    }

    /** Mirrors GenericPDBBench#formatPfunc: 6 columns "lb,ub,status,eps,nconf,nscored". */
    private static String formatPfunc(edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Result r) {
        if (r == null) return ",,N/A,,0,0";
        Double lb = KStarScore.scoreToLog10(r.values.calcLowerBound());
        Double ub = KStarScore.scoreToLog10(r.values.calcUpperBound());
        String lbStr = (lb == null || lb.isNaN()) ? "" : String.format("%.6f", lb);
        String ubStr = (ub == null || ub.isNaN()) ? "" : String.format("%.6f", ub);
        String epsStr = "";
        try {
            double eps = r.values.getEffectiveEpsilon();
            epsStr = Double.isNaN(eps) || Double.isInfinite(eps) ? "" : String.format("%.6f", eps);
        } catch (RuntimeException e) {
            // MARK*/PACK* bounds can use MagicBigDecimal infinities; keep the row
            // writable even when epsilon is not numerically meaningful.
        }
        return String.format("%s,%s,%s,%s,%d,%d",
                lbStr, ubStr, r.status.name(), epsStr, r.numConfs,
                r.getStat("numConfsScored"));
    }
}
