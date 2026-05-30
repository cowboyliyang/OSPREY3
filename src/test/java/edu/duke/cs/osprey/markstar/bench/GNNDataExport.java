package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.io.File;

/**
 * GNN training-data exporters for the 2RL0 system.
 *
 * Three entry points (all callable from {@link TestBranchMARKStar} main dispatch):
 *   - {@link #exportGNNDataAllMutable()}    — 4 pos × 20 AA, all-mutable
 *   - {@link #exportGNNDataHighrot8pos()}   — 8 pos × 5 AA, highrot
 *   - {@link #patchRigidEmatHighrot8pos()}  — backfill E_rigid column on cached CSV
 */
public class GNNDataExport {

    /** Export GNN training data: 4 positions × 20 AA = 160k sequences. */
    @Test
    public void exportGNNDataAllMutable() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_all20_4pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", TestBranchMARKStar.NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildAllMutableConfSpace();
        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            System.out.println("\n========== Exporting PROTEIN data (4pos × 20AA) ==========");
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

            System.out.println("\n========== Exporting COMPLEX data (4pos × 20AA) ==========");
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

            System.out.println("\n========== Export complete ==========");
        }
    }

    /** Export GNN training data: 8 positions × 5 AA = 390,625 sequences. */
    @Test
    public void exportGNNDataHighrot8pos() throws Exception {
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 300000);
        String baseDir = System.getProperty("osprey.gnn.outputDir", "gnn_data/2RL0_highrot_8pos");
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", TestBranchMARKStar.NUM_CPUs);

        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace8pos();
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

            for (String strand : new String[]{"protein", "complex"}) {
                System.out.println("\n========== Exporting " + strand.toUpperCase()
                        + " data (8pos highrot) ==========");
                SimpleConfSpace cs = strand.equals("protein") ? confSpaces.protein : confSpaces.complex;
                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMinimized)
                                .build().calcReferenceEnergies())
                        .build();
                EnergyMatrix ematMin = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                        .build().calcEnergyMatrix();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build().calcEnergyMatrix();
                RCs rcs = new RCs(cs);
                InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                        cs, ematRigid, ematMin, rcs, 8.0, 0.1);
                new GNNDataExporter(confEcalc, ematMin, ematRigid, ig, rcs)
                        .export(numSamples, new File(baseDir + "/" + strand));
            }
            System.out.println("\n========== Export complete ==========");
        }
    }

    /**
     * Backfill E_rigid column in cached confs.csv for the 8-pos highrot dataset.
     * Reads from osprey.gnn.dataDir (default gnn_data/2RL0_highrot_8pos), patches
     * NaN E_rigid values using the rigid emat, and atomically rewrites the CSV.
     */
    public void patchRigidEmatHighrot8pos() throws Exception {
        String baseDir = System.getProperty("osprey.gnn.dataDir", "gnn_data/2RL0_highrot_8pos");

        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildHighRotamerConfSpace8pos();
        Parallelism parallelism = Parallelism.makeCpu(
                Integer.getInteger("osprey.gnn.numCPUs", TestBranchMARKStar.NUM_CPUs));

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

                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMinimized)
                                .build().calcReferenceEnergies())
                        .build();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build().calcEnergyMatrix();
                int numPos = cs.positions.size();

                File tmpFile = new File(confsFile.getPath() + ".tmp");
                int count = 0;
                try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(confsFile), 1 << 20);
                     java.io.PrintWriter pw = new java.io.PrintWriter(
                             new java.io.BufferedWriter(new java.io.FileWriter(tmpFile), 1 << 20))) {
                    String header = br.readLine();
                    pw.println(header);
                    String line;
                    while ((line = br.readLine()) != null) {
                        line = line.trim();
                        if (line.isEmpty()) continue;
                        String[] parts = line.split(",");
                        int[] conf = new int[numPos];
                        for (int p = 0; p < numPos; p++) conf[p] = Integer.parseInt(parts[p]);
                        double eRigid = ematRigid.confE(conf);

                        StringBuilder sb = new StringBuilder();
                        for (int p = 0; p < numPos; p++) sb.append(parts[p]).append(',');
                        sb.append(parts[numPos]).append(',');      // E_CCD
                        sb.append(parts[numPos + 1]).append(',');  // E_emat
                        sb.append(String.format("%.6f", eRigid)).append(',');
                        sb.append(parts[numPos + 3]);              // residual
                        pw.println(sb);
                        count++;
                        if (count % 100000 == 0) System.out.println("  processed " + count);
                    }
                    System.out.println("  total: " + count + " conformations patched");
                }
                java.nio.file.Files.move(tmpFile.toPath(), confsFile.toPath(),
                        java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                System.out.println("  saved: " + confsFile);
            }
        }
        System.out.println("\nDone patching E_rigid (highrot 8pos).");
    }
}
