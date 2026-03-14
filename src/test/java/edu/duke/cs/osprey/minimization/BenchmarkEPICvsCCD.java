package edu.duke.cs.osprey.minimization;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.NewEPICMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Benchmark: EPIC vs CCD on 3 residues (no SAPE).
 * Compares precomputation cost, per-conformation minimization time, and energy quality.
 */
public class BenchmarkEPICvsCCD {

    @Test
    public void benchmarkEPICvsCCD_3res() {
        // 1. Build conf space: 7 flexible residues (2 mutable + 5 flexible)
        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        System.out.println("=== EPIC vs CCD Benchmark (7 residues, no SAPE) ===");
        System.out.println("Conf space: " + confSpace.positions.size() + " positions");
        for (int i = 0; i < confSpace.positions.size(); i++) {
            SimpleConfSpace.Position pos = confSpace.positions.get(i);
            System.out.println("  pos " + i + " (" + pos.resNum + "): " + pos.resConfs.size() + " RCs");
        }

        Parallelism parallelism = Parallelism.makeCpu(1);

        // 2. Energy calculators and matrices
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, rigidEcalc);

        System.out.println("\nComputing energy matrices...");
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build().calcEnergyMatrix();
        System.out.println("Energy matrices computed.");

        // 3. Compute EPIC matrix (no SAPE)
        System.out.println("\nComputing EPIC matrix (useSAPE=false)...");
        EPICSettings epicSettings = EPICSettings.defaultEPIC();
        try {
            java.lang.reflect.Field f = EPICSettings.class.getDeclaredField("useSAPE");
            f.setAccessible(true);
            f.setBoolean(epicSettings, false);
        } catch (Exception e) {
            throw new RuntimeException("Failed to disable SAPE", e);
        }
        PruningMatrix pruneMat = new PruningMatrix(confSpace);

        long epicPrecompStart = System.nanoTime();
        NewEPICMatrixCalculator epicCalc = new NewEPICMatrixCalculator(
                confSpace, confEcalcMinimized, pruneMat, epicSettings);
        epicCalc.calcPEM();
        NewEPICMatrix epicMat = epicCalc.getEPICMatrix();
        long epicPrecompNs = System.nanoTime() - epicPrecompStart;
        System.out.println(String.format("EPIC precomputation: %.1f ms", epicPrecompNs / 1e6));

        // 4. Get conformations via A*
        int numConfs = 20;
        RCs rcs = new RCs(confSpace);
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        List<ScoredConf> scoredConfs = new ArrayList<>();
        for (int i = 0; i < numConfs; i++) {
            ScoredConf sc = astar.nextConf();
            if (sc == null) break;
            scoredConfs.add(sc);
        }
        int n = scoredConfs.size();
        System.out.println("\nGot " + n + " conformations.\n");

        // JIT warm-up
        if (!scoredConfs.isEmpty()) {
            int[] warmConf = scoredConfs.get(0).getAssignments();
            confEcalcMinimized.calcEnergy(new RCTuple(warmConf));
            epicMat.minimizeEnergy(warmConf);
        }

        // 5. Benchmark
        double[] rigidEnergies = new double[n];
        double[] ccdEnergies = new double[n];
        double[] epicEnergies = new double[n];
        long[] ccdTimesNs = new long[n];
        long[] epicTimesNs = new long[n];

        // Rigid energies
        for (int i = 0; i < n; i++) {
            rigidEnergies[i] = confEcalcRigid.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
        }

        // CCD (default maxIter=30)
        long ccdTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            long t0 = System.nanoTime();
            ccdEnergies[i] = confEcalcMinimized.calcEnergy(new RCTuple(scoredConfs.get(i).getAssignments())).energy;
            ccdTimesNs[i] = System.nanoTime() - t0;
        }
        long ccdTotalNs = System.nanoTime() - ccdTotalStart;

        // EPIC minimization
        long epicTotalStart = System.nanoTime();
        for (int i = 0; i < n; i++) {
            long t0 = System.nanoTime();
            epicEnergies[i] = epicMat.minimizeEnergy(scoredConfs.get(i).getAssignments());
            epicTimesNs[i] = System.nanoTime() - t0;
        }
        long epicTotalNs = System.nanoTime() - epicTotalStart;

        // 6. Print results
        System.out.println(String.format("%-5s| %11s | %11s | %8s | %8s | %8s | %8s",
                "Conf", "CCD_E", "EPIC_E", "CCD_drop", "CCD_ms", "EPIC_ms", "Gap"));
        System.out.println("-".repeat(80));

        double sumGap = 0, sumAbsGap = 0;
        double sumCCDMs = 0, sumEPICMs = 0;
        for (int i = 0; i < n; i++) {
            double ccdDrop = rigidEnergies[i] - ccdEnergies[i];
            double gap = epicEnergies[i] - ccdEnergies[i];
            double ccdMs = ccdTimesNs[i] / 1e6;
            double epicMs = epicTimesNs[i] / 1e6;

            sumGap += gap;
            sumAbsGap += Math.abs(gap);
            sumCCDMs += ccdMs;
            sumEPICMs += epicMs;

            System.out.println(String.format("%-5d| %11.4f | %11.4f | %8.4f | %8.2f | %8.2f | %+8.4f",
                    i, ccdEnergies[i], epicEnergies[i], ccdDrop, ccdMs, epicMs, gap));
        }
        System.out.println("-".repeat(80));

        System.out.println(String.format("\nSummary (%d confs):", n));
        System.out.println(String.format("  EPIC precomputation:   %10.1f ms", epicPrecompNs / 1e6));
        System.out.println(String.format("  CCD total minimize:    %10.1f ms  (avg %.2f ms/conf)", ccdTotalNs / 1e6, sumCCDMs / n));
        System.out.println(String.format("  EPIC total minimize:   %10.1f ms  (avg %.2f ms/conf)", epicTotalNs / 1e6, sumEPICMs / n));
        System.out.println(String.format("  EPIC speedup vs CCD:   %.2fx (minimize only)", (double) ccdTotalNs / epicTotalNs));
        System.out.println(String.format("  Avg gap (EPIC-CCD):    %+.4f kcal/mol", sumGap / n));
        System.out.println(String.format("  Avg |gap|:             %.4f kcal/mol", sumAbsGap / n));

        // Break-even: how many confs before EPIC precomp + EPIC min < CCD min?
        double avgCCDMs = sumCCDMs / n;
        double avgEPICMs = sumEPICMs / n;
        if (avgCCDMs > avgEPICMs) {
            double breakEven = (epicPrecompNs / 1e6) / (avgCCDMs - avgEPICMs);
            System.out.println(String.format("  Break-even point:      %.0f confs (EPIC precomp amortized)", breakEven));
        } else {
            System.out.println("  Break-even point:      N/A (EPIC not faster per conf)");
        }

        minimizingEcalc.close();
        rigidEcalc.close();
    }
}
