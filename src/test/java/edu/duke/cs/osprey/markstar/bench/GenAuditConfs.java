package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.TestKStar;

import java.io.PrintWriter;
import java.util.Random;

/** Emit N random valid confs for a design/state as an audit-leaf CSV (for CCD throughput probing). */
public class GenAuditConfs {
    public static void main(String[] args) throws Exception {
        String pdb = System.getProperty("osprey.bench.pdbPath");
        String did = System.getProperty("osprey.bench.designId", "unknown");
        String mut = System.getProperty("osprey.bench.mutable", "");
        String flex = System.getProperty("osprey.bench.flexible", "");
        String state = System.getProperty("osprey.audit.state", "complex");
        String out = System.getProperty("osprey.audit.output");
        int n = Integer.getInteger("osprey.gen.n", 600);
        long seed = Long.getLong("osprey.gen.seed", 42L);

        TestKStar.ConfSpaces cs = AuditLeafCCD.buildConfSpaces(pdb, mut, flex);
        SimpleConfSpace conf = AuditLeafCCD.selectConfSpace(cs, state);
        int npos = conf.positions.size();
        int[] sizes = new int[npos];
        for (int i = 0; i < npos; i++) sizes[i] = conf.positions.get(i).resConfs.size();
        System.out.println("design=" + did + " state=" + state + " positions=" + npos);
        StringBuilder sz = new StringBuilder("rcs per pos:");
        long space = 1;
        for (int i = 0; i < npos; i++) { sz.append(' ').append(sizes[i]); space *= sizes[i]; }
        System.out.println(sz + "  (conf space=" + space + ")");

        Random rng = new Random(seed);
        try (PrintWriter pw = new PrintWriter(out)) {
            pw.println("run_id,design_id,pdb,state,sequence_index,sequence,conf_id,assignments,gnn_energy_kcal");
            for (int k = 0; k < n; k++) {
                StringBuilder a = new StringBuilder();
                for (int i = 0; i < npos; i++) {
                    if (i > 0) a.append(';');
                    a.append(rng.nextInt(sizes[i]));
                }
                pw.println("gen," + did + "," + did + "," + state + ",0,WT,c" + k + "," + a + ",0.0");
            }
        }
        System.out.println("wrote " + n + " confs -> " + out);
    }
}
