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
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Progress;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Exports conformation energy data for GNN training.
 *
 * Uses A* enumeration to sample conformations in emat energy order,
 * ensuring only low-energy (partition-function-relevant) conformations are included.
 *
 * CCD energies are computed in parallel via calcEnergyAsync.
 *
 * Outputs:
 *   1. confs.csv   — rc assignments + E_CCD + E_emat + residual
 *   2. graph.csv   — interaction graph edge list
 *   3. meta.csv    — per-position RC counts
 */
public class GNNDataExporter {

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
        if (!confs.isEmpty()) {
            double firstEmat = ematMinimized.confE(confs.get(0));
            double lastEmat = ematMinimized.confE(confs.get(confs.size() - 1));
            System.out.println("  emat range: [" + String.format("%.4f", firstEmat)
                    + ", " + String.format("%.4f", lastEmat) + "]");
        }

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

        // 3. Write outputs
        File confsFile = new File(outputDir, "confs.csv");
        int written = writeConfs(confsFile, confs, eCCD, eEmat);

        File graphFile = new File(outputDir, "graph.csv");
        writeGraph(graphFile);

        File metaFile = new File(outputDir, "meta.csv");
        writeMeta(metaFile);

        System.out.println("Export complete: " + outputDir.getAbsolutePath());
        System.out.println("  confs.csv:  " + written + " valid samples (of " + confs.size() + " enumerated)");
        System.out.println("  graph.csv:  " + interactionGraph.getNumEdges() + " edges");
        System.out.println("  meta.csv:   " + confSpace.positions.size() + " positions");
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
                if (Double.isNaN(eCCD[i])) {
                    continue;
                }
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
}
