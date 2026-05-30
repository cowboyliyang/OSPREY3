package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Post hoc CCD audit for S11 leaf-GNN replacement rows.
 */
public class AuditLeafCCD {

    private static final double KT = 0.5922;

    public static void main(String[] args) throws Exception {

        String input = requireProperty("osprey.audit.input");
        String output = requireProperty("osprey.audit.output");
        String state = System.getProperty("osprey.audit.state", "complex").trim();
        String designId = System.getProperty("osprey.bench.designId", "unknown");
        String pdbPath = requireProperty("osprey.bench.pdbPath");
        String mutableStr = System.getProperty("osprey.bench.mutable", "");
        String flexibleStr = System.getProperty("osprey.bench.flexible", "");
        int cpus = Integer.getInteger("osprey.bench.numCPUs", 8);
        int shardIndex = Integer.getInteger("osprey.audit.shardIndex", 0);
        int numShards = Math.max(1, Integer.getInteger("osprey.audit.numShards", 1));

        System.out.println("==============================================");
        System.out.println("  S11 Leaf CCD Audit");
        System.out.println("  Design: " + designId);
        System.out.println("  State: " + state);
        System.out.println("  Input: " + input);
        System.out.println("  Output: " + output);
        System.out.println("  Shard: " + shardIndex + " / " + numShards);
        System.out.println("==============================================");

        TestKStar.ConfSpaces confSpaces = buildConfSpaces(pdbPath, mutableStr, flexibleStr);
        SimpleConfSpace confSpace = selectConfSpace(confSpaces, state);
        if (confSpace.positions.isEmpty()) {
            throw new IllegalArgumentException("Selected state has no positions: " + state);
        }

        Parallelism parallelism = Parallelism.makeCpu(cpus);
        long startedMs = System.currentTimeMillis();
        int seen = 0;
        int audited = 0;
        int failed = 0;

        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, confSpaces.ffparams)
                .setParallelism(parallelism).build()) {

            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                            .build().calcReferenceEnergies())
                    .build();

            File outFile = new File(output);
            File parent = outFile.getParentFile();
            if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
                throw new IOException("Failed to create output directory: " + parent.getAbsolutePath());
            }

            try (PrintWriter pw = new PrintWriter(new FileWriter(outFile))) {
                pw.println("run_id,design_id,pdb,state,sequence_index,sequence,conf_id,assignments,"
                        + "gnn_energy_kcal,ccd_energy_kcal,delta_kcal,log_zhat,log_ztrue,"
                        + "audit_time_ms,status,error,audit_partition,audit_inclusion_prob,audit_ht_weight");

                List<File> inputFiles = collectInputFiles(new File(input));
                for (File csv : inputFiles) {
                    try (BufferedReader br = new BufferedReader(new FileReader(csv))) {
                        String header = br.readLine();
                        if (header == null) continue;
                        List<String> columns = parseCsvLine(header);
                        Map<String, Integer> index = indexColumns(columns);

                        String line;
                        while ((line = br.readLine()) != null) {
                            if (line.trim().isEmpty()) continue;
                            Map<String, String> row = rowMap(index, parseCsvLine(line));

                            if (!matches(row.get("design_id"), designId)) continue;
                            if (!matchesState(row.get("state"), state)) continue;

                            int rowIndex = seen++;
                            if (Math.floorMod(rowIndex, numShards) != shardIndex) continue;

                            long t0 = System.nanoTime();
                            try {
                                int[] assignments = parseAssignments(row.get("assignments"));
                                double gnnEnergy = parseDouble(row.get("gnn_energy_kcal"));
                                ConfSearch.ScoredConf scored = new ConfSearch.ScoredConf(assignments, gnnEnergy);
                                double ccdEnergy = confEcalc.calcEnergy(scored).getEnergy();
                                double elapsedMs = (System.nanoTime() - t0) / 1e6;
                                audited++;
                                pw.println(formatOutput(row, gnnEnergy, ccdEnergy, elapsedMs, "ok", ""));
                            } catch (Exception e) {
                                double elapsedMs = (System.nanoTime() - t0) / 1e6;
                                failed++;
                                pw.println(formatOutput(row, Double.NaN, Double.NaN, elapsedMs,
                                        "error", e.getClass().getSimpleName() + ": " + e.getMessage()));
                            }
                        }
                    }
                }
            }
        }

        double elapsedS = (System.currentTimeMillis() - startedMs) / 1000.0;
        System.out.println("Audit done: seen=" + seen
                + ", audited=" + audited
                + ", failed=" + failed
                + ", elapsed=" + String.format(Locale.US, "%.1fs", elapsedS));
    }

    private static TestKStar.ConfSpaces buildConfSpaces(String pdbPath, String mutableStr,
                                                        String flexibleStr) {

        String[] mutableResidues = mutableStr.isEmpty() ? new String[0] : mutableStr.split(";");
        String[] flexibleResidues = flexibleStr.isEmpty() ? new String[0] : flexibleStr.split(";");

        ForcefieldParams ffparams = new ForcefieldParams();
        Molecule mol = PDBIO.readFile(pdbPath);
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY",
                "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
                "THR","TRP","TYR","VAL"};

        Map<String, String[]> chainRanges = getChainRanges(mol);
        Set<String> proteinChainSet = new LinkedHashSet<>();
        for (String res : mutableResidues) {
            res = res.trim();
            if (!res.isEmpty() && res.length() >= 2) {
                proteinChainSet.add(String.valueOf(res.charAt(0)));
            }
        }
        Set<String> ligandChainSet = new LinkedHashSet<>();
        for (String ch : chainRanges.keySet()) {
            if (!proteinChainSet.contains(ch) && !ch.trim().isEmpty()) {
                ligandChainSet.add(ch);
            }
        }

        String pFirst = null;
        String pLast = null;
        for (String ch : proteinChainSet) {
            String[] range = chainRanges.get(ch);
            if (range != null) {
                if (pFirst == null) pFirst = range[0];
                pLast = range[1];
            }
        }
        if (pFirst == null || pLast == null) {
            throw new IllegalArgumentException("Could not infer protein strand from mutable residues: " + mutableStr);
        }

        Strand proteinStrand = new Strand.Builder(mol).setTemplateLibrary(templateLib)
                .setResidues(pFirst, pLast).build();

        for (String res : mutableResidues) {
            res = res.trim();
            if (res.isEmpty()) continue;
            proteinStrand.flexibility.get(res).setLibraryRotamers(all20)
                    .addWildTypeRotamers().setContinuous();
        }
        for (String res : flexibleResidues) {
            res = res.trim();
            if (res.isEmpty()) continue;
            if (!proteinChainSet.contains(String.valueOf(res.charAt(0)))) continue;
            proteinStrand.flexibility.get(res).setLibraryRotamers(Strand.WildType)
                    .addWildTypeRotamers().setContinuous();
        }

        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = ffparams;

        if (!ligandChainSet.isEmpty()) {
            List<Strand> ligandStrands = new ArrayList<>();
            for (String ch : ligandChainSet) {
                String[] range = chainRanges.get(ch);
                if (range == null) continue;
                Strand ligStrand = new Strand.Builder(mol).setTemplateLibrary(templateLib)
                        .setResidues(range[0], range[1]).build();
                for (String res : flexibleResidues) {
                    res = res.trim();
                    if (res.isEmpty()) continue;
                    if (!ch.equals(String.valueOf(res.charAt(0)))) continue;
                    ligStrand.flexibility.get(res).setLibraryRotamers(Strand.WildType)
                            .addWildTypeRotamers().setContinuous();
                }
                ligandStrands.add(ligStrand);
            }
            confSpaces.protein = new SimpleConfSpace.Builder().addStrand(proteinStrand).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().addStrands(ligandStrands).build();
            SimpleConfSpace.Builder complexBuilder = new SimpleConfSpace.Builder().addStrand(proteinStrand);
            for (Strand ligandStrand : ligandStrands) complexBuilder.addStrand(ligandStrand);
            confSpaces.complex = complexBuilder.build();
        } else {
            confSpaces.protein = new SimpleConfSpace.Builder().addStrand(proteinStrand).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().build();
            confSpaces.complex = confSpaces.protein;
        }

        return confSpaces;
    }

    private static SimpleConfSpace selectConfSpace(TestKStar.ConfSpaces confSpaces, String state) {
        String normalized = normalizeState(state);
        if (normalized.equals("protein")) return confSpaces.protein;
        if (normalized.equals("ligand")) return confSpaces.ligand;
        if (normalized.equals("complex")) return confSpaces.complex;
        throw new IllegalArgumentException("Unknown state: " + state);
    }

    private static String formatOutput(Map<String, String> row, double gnnEnergyOverride,
                                       double ccdEnergy, double auditTimeMs,
                                       String status, String error) {
        double gnnEnergy = Double.isNaN(gnnEnergyOverride)
                ? parseDouble(row.get("gnn_energy_kcal"))
                : gnnEnergyOverride;
        double delta = Double.isNaN(ccdEnergy) || Double.isNaN(gnnEnergy)
                ? Double.NaN : ccdEnergy - gnnEnergy;
        double logZtrue = Double.isNaN(ccdEnergy) ? Double.NaN : -ccdEnergy / KT;
        return csv(row.get("run_id")) + ","
                + csv(row.get("design_id")) + ","
                + csv(row.get("pdb")) + ","
                + csv(row.get("state")) + ","
                + csv(row.get("sequence_index")) + ","
                + csv(row.get("sequence")) + ","
                + csv(row.get("conf_id")) + ","
                + csv(row.get("assignments")) + ","
                + formatDouble(gnnEnergy) + ","
                + formatDouble(ccdEnergy) + ","
                + formatDouble(delta) + ","
                + csv(row.get("log_zhat")) + ","
                + formatDouble(logZtrue) + ","
                + formatDouble(auditTimeMs) + ","
                + csv(status) + ","
                + csv(error) + ","
                + csv(row.get("audit_partition")) + ","
                + csv(row.get("audit_inclusion_prob")) + ","
                + csv(row.get("audit_ht_weight"));
    }

    private static List<File> collectInputFiles(File input) {
        List<File> files = new ArrayList<>();
        collectInputFiles(input, files);
        files.sort((a, b) -> a.getAbsolutePath().compareTo(b.getAbsolutePath()));
        return files;
    }

    private static void collectInputFiles(File input, List<File> files) {
        if (!input.exists()) {
            throw new IllegalArgumentException("Audit input does not exist: " + input.getAbsolutePath());
        }
        if (input.isFile()) {
            if (input.getName().endsWith(".csv") && !input.getName().endsWith(".summary.csv")) {
                files.add(input);
            }
            return;
        }
        File[] children = input.listFiles();
        if (children == null) return;
        for (File child : children) collectInputFiles(child, files);
    }

    private static Map<String, Integer> indexColumns(List<String> columns) {
        Map<String, Integer> index = new HashMap<>();
        for (int i = 0; i < columns.size(); i++) index.put(columns.get(i), i);
        return index;
    }

    private static Map<String, String> rowMap(Map<String, Integer> index, List<String> values) {
        Map<String, String> row = new HashMap<>();
        for (Map.Entry<String, Integer> entry : index.entrySet()) {
            int i = entry.getValue();
            row.put(entry.getKey(), i < values.size() ? values.get(i) : "");
        }
        return row;
    }

    private static List<String> parseCsvLine(String line) {
        List<String> out = new ArrayList<>();
        StringBuilder cell = new StringBuilder();
        boolean quoted = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (quoted) {
                if (c == '"') {
                    if (i + 1 < line.length() && line.charAt(i + 1) == '"') {
                        cell.append('"');
                        i++;
                    } else {
                        quoted = false;
                    }
                } else {
                    cell.append(c);
                }
            } else {
                if (c == '"') {
                    quoted = true;
                } else if (c == ',') {
                    out.add(cell.toString());
                    cell.setLength(0);
                } else {
                    cell.append(c);
                }
            }
        }
        out.add(cell.toString());
        return out;
    }

    private static int[] parseAssignments(String value) {
        if (value == null || value.trim().isEmpty()) {
            throw new IllegalArgumentException("empty assignments");
        }
        String[] parts = value.split(";");
        int[] assignments = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            assignments[i] = Integer.parseInt(parts[i].trim());
        }
        return assignments;
    }

    private static boolean matches(String rowValue, String expected) {
        return expected == null || expected.equals("unknown") || expected.trim().isEmpty()
                || expected.equals(rowValue);
    }

    private static boolean matchesState(String rowValue, String expected) {
        return normalizeState(rowValue).equals(normalizeState(expected));
    }

    private static String normalizeState(String state) {
        return state == null ? "" : state.trim().toLowerCase(Locale.US);
    }

    private static String requireProperty(String key) {
        String value = System.getProperty(key);
        if (value == null || value.trim().isEmpty()) {
            throw new IllegalArgumentException("Missing required property: " + key);
        }
        return value;
    }

    private static double parseDouble(String value) {
        if (value == null || value.trim().isEmpty()) return Double.NaN;
        return Double.parseDouble(value.trim());
    }

    private static String csv(String value) {
        if (value == null) return "";
        boolean quote = value.indexOf(',') >= 0 || value.indexOf('"') >= 0
                || value.indexOf('\n') >= 0 || value.indexOf('\r') >= 0;
        if (!quote) return value;
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    private static String formatDouble(double value) {
        if (Double.isNaN(value)) return "NaN";
        if (value == Double.POSITIVE_INFINITY) return "Infinity";
        if (value == Double.NEGATIVE_INFINITY) return "-Infinity";
        return String.format(Locale.US, "%.12g", value);
    }

    private static Map<String, String[]> getChainRanges(Molecule mol) {
        Map<String, String[]> ranges = new LinkedHashMap<>();
        Map<String, String> firstRes = new LinkedHashMap<>();
        Map<String, String> lastRes = new LinkedHashMap<>();
        for (edu.duke.cs.osprey.structure.Residue res : mol.residues) {
            String chain = String.valueOf(res.getChainId());
            String fullId = res.getPDBResNumber().trim();
            if (!firstRes.containsKey(chain)) firstRes.put(chain, fullId);
            lastRes.put(chain, fullId);
        }
        for (String ch : firstRes.keySet()) {
            ranges.put(ch, new String[]{firstRes.get(ch), lastRes.get(ch)});
        }
        return ranges;
    }
}
