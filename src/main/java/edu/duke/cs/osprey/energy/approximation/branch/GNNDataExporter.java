package edu.duke.cs.osprey.energy.approximation.branch;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Progress;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Exports conformation energy data for GNN training (v2).
 *
 * Outputs:
 *   1. confs.csv           — rc assignments + E_CCD + E_emat + residual
 *   2. graph.csv           — interaction graph edge list
 *   3. meta.csv            — per-position RC counts
 *   4. rc_features.csv     — per-RC: AA type index, chi angles
 *   5. pairwise_energies.csv — emat pairwise energies for all RC pairs
 *   6. ca_distances.csv    — Cα-Cα distances between position pairs
 */
public class GNNDataExporter {

    // Must match Python side (train.py AA_TYPES dict)
    private static final Map<String, Integer> AA_TYPE_MAP = new LinkedHashMap<>();
    static {
        AA_TYPE_MAP.put("ALA", 0);  AA_TYPE_MAP.put("ARG", 1);
        AA_TYPE_MAP.put("ASN", 2);  AA_TYPE_MAP.put("ASP", 3);
        AA_TYPE_MAP.put("CYS", 4);  AA_TYPE_MAP.put("GLN", 5);
        AA_TYPE_MAP.put("GLU", 6);  AA_TYPE_MAP.put("GLY", 7);
        AA_TYPE_MAP.put("HIS", 8);  AA_TYPE_MAP.put("HIE", 8);
        AA_TYPE_MAP.put("HID", 8);  AA_TYPE_MAP.put("HIP", 8);
        AA_TYPE_MAP.put("ILE", 9);
        AA_TYPE_MAP.put("LEU", 10); AA_TYPE_MAP.put("LYS", 11);
        AA_TYPE_MAP.put("MET", 12); AA_TYPE_MAP.put("PHE", 13);
        AA_TYPE_MAP.put("PRO", 14); AA_TYPE_MAP.put("SER", 15);
        AA_TYPE_MAP.put("THR", 16); AA_TYPE_MAP.put("TRP", 17);
        AA_TYPE_MAP.put("TYR", 18); AA_TYPE_MAP.put("VAL", 19);
    }
    private static final int MAX_CHI = 4;

    private final SimpleConfSpace confSpace;
    private final ConfEnergyCalculator confEcalc;
    private final EnergyMatrix ematMinimized;
    private final InteractionGraph interactionGraph;
    private final RCs rcs;

    public GNNDataExporter(
            ConfEnergyCalculator confEcalc,
            EnergyMatrix ematMinimized,
            InteractionGraph interactionGraph,
            RCs rcs
    ) {
        this.confSpace = confEcalc.confSpace;
        this.confEcalc = confEcalc;
        this.ematMinimized = ematMinimized;
        this.interactionGraph = interactionGraph;
        this.rcs = rcs;
    }

    public void export(int numSamples, File outputDir) throws IOException {
        outputDir.mkdirs();

        // 1. Sample conformations via A* (low energy first)
        System.out.println("Enumerating " + numSamples + " lowest-energy conformations via A*...");
        List<int[]> confs = sampleByAStar(numSamples);
        System.out.println("Enumerated " + confs.size() + " conformations.");

        // 2. Compute CCD energies in parallel
        double[] eCCD = new double[confs.size()];
        double[] eEmat = new double[confs.size()];
        AtomicInteger failures = new AtomicInteger(0);

        Progress progress = new Progress(confs.size());
        progress.setReportMemory(true);
        System.out.println("Computing CCD energies for " + confs.size() + " conformations...");

        for (int i = 0; i < confs.size(); i++) {
            int[] conf = confs.get(i);
            eEmat[i] = ematMinimized.confE(conf);

            final int idx = i;
            confEcalc.calcEnergyAsync(new RCTuple(conf), (EnergyCalculator.EnergiedParametricMolecule epm) -> {
                if (epm != null && Double.isFinite(epm.energy)) {
                    eCCD[idx] = epm.energy;
                } else {
                    eCCD[idx] = Double.NaN;
                    failures.incrementAndGet();
                }
                progress.incrementProgress();
            });
        }
        confEcalc.tasks.waitForFinish();

        int numFailures = failures.get();
        if (numFailures > 0) {
            System.out.println("WARNING: " + numFailures + " CCD minimizations failed (NaN).");
        }

        // 3. Write all output files
        int written = writeConfs(new File(outputDir, "confs.csv"), confs, eCCD, eEmat);
        writeGraph(new File(outputDir, "graph.csv"));
        writeMeta(new File(outputDir, "meta.csv"));
        writeRCFeatures(new File(outputDir, "rc_features.csv"));
        writePairwiseEnergies(new File(outputDir, "pairwise_energies.csv"));
        writeCaDistances(new File(outputDir, "ca_distances.csv"));

        System.out.println("Export complete: " + outputDir.getAbsolutePath());
        System.out.println("  confs.csv:             " + written + " valid samples");
        System.out.println("  graph.csv:             " + interactionGraph.getNumEdges() + " edges");
        System.out.println("  meta.csv:              " + confSpace.positions.size() + " positions");
        System.out.println("  rc_features.csv:       RC features exported");
        System.out.println("  pairwise_energies.csv: pairwise energies exported");
        System.out.println("  ca_distances.csv:      Cα distances exported");
    }

    private List<int[]> sampleByAStar(int numSamples) {
        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional()
                .build();

        List<int[]> confs = new ArrayList<>(numSamples);
        for (int i = 0; i < numSamples; i++) {
            ConfSearch.ScoredConf sc = astar.nextConf();
            if (sc == null) {
                System.out.println("A* exhausted after " + i + " conformations.");
                break;
            }
            confs.add(sc.getAssignments());
        }
        return confs;
    }

    private int writeConfs(File file, List<int[]> confs, double[] eCCD, double[] eEmat) throws IOException {
        int n = confSpace.positions.size();
        int written = 0;
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
            StringBuilder header = new StringBuilder();
            for (int p = 0; p < n; p++) {
                header.append("rc_").append(p).append(',');
            }
            header.append("E_CCD,E_emat,residual");
            pw.println(header);

            for (int i = 0; i < confs.size(); i++) {
                if (Double.isNaN(eCCD[i])) continue;
                int[] conf = confs.get(i);
                StringBuilder sb = new StringBuilder();
                for (int p = 0; p < n; p++) {
                    sb.append(conf[p]).append(',');
                }
                sb.append(String.format("%.6f", eCCD[i])).append(',');
                sb.append(String.format("%.6f", eEmat[i])).append(',');
                sb.append(String.format("%.6f", eCCD[i] - eEmat[i]));
                pw.println(sb);
                written++;
            }
        }
        return written;
    }

    private void writeGraph(File file) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
            pw.println("src,dst");
            for (int[] edge : interactionGraph.getEdgeList()) {
                pw.println(edge[0] + "," + edge[1]);
                pw.println(edge[1] + "," + edge[0]);
            }
        }
    }

    private void writeMeta(File file) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
            pw.println("pos,num_rcs,res_type");
            for (SimpleConfSpace.Position pos : confSpace.positions) {
                pw.println(pos.index + "," + pos.resConfs.size() + "," + pos.resNum);
            }
        }
    }

    private void writeRCFeatures(File file) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
            pw.println("pos,rc,aa_type_idx,chi1,chi2,chi3,chi4");

            for (SimpleConfSpace.Position pos : confSpace.positions) {
                for (int r = 0; r < pos.resConfs.size(); r++) {
                    SimpleConfSpace.ResidueConf rc = pos.resConfs.get(r);
                    String aaName = rc.template.name;
                    Integer aaIdx = AA_TYPE_MAP.get(aaName);
                    if (aaIdx == null) {
                        // Unknown AA type, use -1
                        System.err.println("WARNING: Unknown AA type " + aaName + ", using index -1");
                        aaIdx = -1;
                    }

                    double[] chis = new double[MAX_CHI];
                    Arrays.fill(chis, 0.0);
                    if (rc.rotamerIndex != null && rc.template.numDihedrals > 0) {
                        try {
                            double[] libChis = rc.template.getRotamericDihedrals(rc.rotamerIndex);
                            for (int d = 0; d < Math.min(libChis.length, MAX_CHI); d++) {
                                chis[d] = libChis[d];
                            }
                        } catch (Exception e) {
                            // keep zeros
                        }
                    }

                    pw.println(pos.index + "," + r + "," + aaIdx + ","
                            + String.format("%.4f", chis[0]) + ","
                            + String.format("%.4f", chis[1]) + ","
                            + String.format("%.4f", chis[2]) + ","
                            + String.format("%.4f", chis[3]));
                }
            }
        }
    }

    private void writePairwiseEnergies(File file) throws IOException {
        int numPos = confSpace.positions.size();
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file), 1 << 20))) {
            pw.println("pos1,rc1,pos2,rc2,E_pair_min");

            for (int p1 = 0; p1 < numPos; p1++) {
                int nrc1 = confSpace.positions.get(p1).resConfs.size();
                for (int p2 = p1 + 1; p2 < numPos; p2++) {
                    int nrc2 = confSpace.positions.get(p2).resConfs.size();
                    for (int r1 = 0; r1 < nrc1; r1++) {
                        for (int r2 = 0; r2 < nrc2; r2++) {
                            double e = ematMinimized.getPairwise(p1, r1, p2, r2);
                            pw.println(p1 + "," + r1 + "," + p2 + "," + r2 + ","
                                    + String.format("%.6f", e));
                        }
                    }
                }
            }
        }
    }

    private void writeCaDistances(File file) throws IOException {
        int numPos = confSpace.positions.size();

        // Build a molecule with first RC at each position to get Cα coords
        int[] defaultConf = new int[numPos];
        Molecule mol = confSpace.makeMolecule(new RCTuple(defaultConf)).mol;

        // Collect Cα coordinates
        double[][] caCoords = new double[numPos][];
        for (int p = 0; p < numPos; p++) {
            String resNum = confSpace.positions.get(p).resNum;
            Residue res = mol.getResByPDBResNumber(resNum);
            Atom ca = res.getAtomByName("CA");
            if (ca != null) {
                caCoords[p] = ca.getCoords();
            } else {
                System.err.println("WARNING: No CA atom for position " + p + " (" + resNum + ")");
                caCoords[p] = new double[]{0, 0, 0};
            }
        }

        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
            pw.println("pos1,pos2,distance");
            for (int p1 = 0; p1 < numPos; p1++) {
                for (int p2 = p1 + 1; p2 < numPos; p2++) {
                    double dx = caCoords[p1][0] - caCoords[p2][0];
                    double dy = caCoords[p1][1] - caCoords[p2][1];
                    double dz = caCoords[p1][2] - caCoords[p2][2];
                    double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    pw.println(p1 + "," + p2 + "," + String.format("%.4f", dist));
                }
            }
        }
    }
}
