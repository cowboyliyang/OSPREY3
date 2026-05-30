package edu.duke.cs.osprey.markstar.framework;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Buffers S11 leaf-GNN replacements and writes the subset selected for post hoc CCD audit.
 */
public class AuditLeafLogger implements AutoCloseable {

    public enum AuditMode {
        OFF(0),
        FULL(1),
        TOP_MASS(2),
        TOP_MASS_SAMPLE(3);

        private final int code;

        AuditMode(int code) {
            this.code = code;
        }

        public int code() {
            return code;
        }

        public static AuditMode parse(String value) {
            if (value == null) return FULL;
            String normalized = value.trim().toLowerCase();
            if (normalized.isEmpty()) return FULL;
            if (normalized.equals("off") || normalized.equals("none") || normalized.equals("false")) {
                return OFF;
            }
            if (normalized.equals("full") || normalized.equals("all")) {
                return FULL;
            }
            if (normalized.equals("top") || normalized.equals("top_mass") || normalized.equals("mass")) {
                return TOP_MASS;
            }
            if (normalized.equals("sample") || normalized.equals("tail_sample")
                    || normalized.equals("top_mass_sample") || normalized.equals("weighted")) {
                return TOP_MASS_SAMPLE;
            }
            throw new IllegalArgumentException("Unknown S11 audit mode: " + value
                    + " (expected off, full, top_mass, or top_mass_sample)");
        }
    }

    private static final double KT = 0.5922;

    private final String runId;
    private final String designId;
    private final String pdb;
    private final String state;
    private final int sequenceIndex;
    private final String sequence;
    private final AuditMode mode;
    private final double topMassFraction;
    private final int tailSampleCount;
    private final long randomSeed;
    private final File outputFile;
    private final File summaryFile;
    private final List<Record> records = new ArrayList<>();

    private boolean closed = false;
    private long eventCounter = 0;
    private int selectedCount = 0;
    private int topSelectedCount = 0;
    private int tailSampleSelectedCount = 0;
    private double totalLogZhat = Double.NEGATIVE_INFINITY;
    private double selectedLogZhat = Double.NEGATIVE_INFINITY;

    public AuditLeafLogger(String designId, String state, int sequenceIndex, String sequence) {
        this.runId = System.getProperty("osprey.gnn.s11.runId",
                System.getProperty("osprey.bench.runId", Long.toString(System.currentTimeMillis())));
        this.designId = emptyToUnknown(designId);
        this.pdb = pdbName(System.getProperty("osprey.bench.pdbPath", "unknown"));
        this.state = emptyToUnknown(state);
        this.sequenceIndex = sequenceIndex;
        this.sequence = emptyToUnknown(sequence);
        this.mode = AuditMode.parse(System.getProperty("osprey.gnn.s11.auditMode", "full"));
        this.topMassFraction = clamp01(Double.parseDouble(
                System.getProperty("osprey.gnn.s11.auditTopMass", "0.999")));
        this.tailSampleCount = Math.max(0, Integer.getInteger("osprey.gnn.s11.auditTailSamples", 1000));
        this.randomSeed = Long.getLong("osprey.gnn.s11.auditSeed", 20260529L);

        String auditDirProp = System.getProperty("osprey.gnn.s11.auditDir");
        String outputDir = System.getProperty("osprey.bench.outputDir", ".");
        File auditDir = auditDirProp == null || auditDirProp.trim().isEmpty()
                ? new File(outputDir, "audit_leaves")
                : new File(auditDirProp);
        File stateDir = new File(new File(auditDir, sanitize(this.designId)), sanitize(this.state));
        this.outputFile = new File(stateDir, String.format("seq_%05d.csv", sequenceIndex));
        this.summaryFile = new File(stateDir, String.format("seq_%05d.summary.csv", sequenceIndex));

        if (mode == AuditMode.OFF) {
            System.out.println("[S11_AUDIT] disabled for " + this.designId + "/" + this.state
                    + " seq=" + sequenceIndex);
        } else {
            System.out.println("[S11_AUDIT] mode=" + mode.name().toLowerCase()
                    + " output=" + outputFile.getAbsolutePath()
                    + " topMass=" + topMassFraction
                    + " tailSamples=" + tailSampleCount);
        }
    }

    public boolean isEnabled() {
        return mode != AuditMode.OFF;
    }

    public int getModeCode() {
        return mode.code();
    }

    public synchronized int getRecordedCount() {
        return records.size();
    }

    public synchronized int getSelectedCount() {
        return selectedCount;
    }

    public synchronized int getTopSelectedCount() {
        return topSelectedCount;
    }

    public synchronized int getTailSampleSelectedCount() {
        return tailSampleSelectedCount;
    }

    public synchronized void record(int[] assignments, int nodeLevel,
                                    double oldLower, double oldUpper,
                                    double gnnEnergy,
                                    double newLower, double newUpper,
                                    boolean rejectedAboveUpper) {
        if (mode == AuditMode.OFF) return;
        if (closed) {
            throw new IllegalStateException("S11 audit logger is already closed");
        }
        int[] copy = Arrays.copyOf(assignments, assignments.length);
        String assignmentsText = assignmentsToString(copy);
        Record rec = new Record();
        rec.eventIndex = ++eventCounter;
        rec.assignments = copy;
        rec.assignmentsText = assignmentsText;
        rec.confId = state + "_" + sequenceIndex + "_" + Integer.toHexString(Arrays.hashCode(copy));
        rec.oldLower = oldLower;
        rec.oldUpper = oldUpper;
        rec.gnnEnergy = gnnEnergy;
        rec.newLower = newLower;
        rec.newUpper = newUpper;
        rec.rejectedAboveUpper = rejectedAboveUpper;
        rec.logZhat = -newUpper / KT;
        rec.nodeLevel = nodeLevel;
        rec.timestampMs = System.currentTimeMillis();
        records.add(rec);
    }

    @Override
    public synchronized void close() {
        if (closed) return;
        closed = true;
        if (mode == AuditMode.OFF) return;

        List<Record> selected = selectRecords();
        selectedCount = selected.size();
        selectedLogZhat = logSumExp(selected);

        File parent = outputFile.getParentFile();
        if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
            throw new RuntimeException("Failed to create S11 audit directory: " + parent.getAbsolutePath());
        }

        try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
            pw.println("run_id,design_id,pdb,state,sequence_index,sequence,conf_id,assignments,"
                    + "old_lower_kcal,old_upper_kcal,gnn_energy_kcal,new_lower_kcal,new_upper_kcal,"
                    + "rejected_above_upper,log_zhat,zhat_rank_hint,event_index,node_level,timestamp_ms,"
                    + "audit_mode,audit_partition,audit_inclusion_prob,audit_ht_weight");
            for (Record rec : selected) {
                pw.println(formatRecord(rec));
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to write S11 audit leaves: " + outputFile.getAbsolutePath(), e);
        }

        try (PrintWriter pw = new PrintWriter(new FileWriter(summaryFile))) {
            pw.println("run_id,design_id,pdb,state,sequence_index,sequence,audit_mode,total_records,"
                    + "selected_records,top_selected,tail_sample_selected,top_mass_fraction,tail_sample_count,"
                    + "log_zhat_total,log_zhat_selected,audit_file");
            pw.println(csv(runId) + "," + csv(designId) + "," + csv(pdb) + "," + csv(state) + ","
                    + sequenceIndex + "," + csv(sequence) + "," + csv(mode.name().toLowerCase()) + ","
                    + records.size() + "," + selectedCount + "," + topSelectedCount + ","
                    + tailSampleSelectedCount + "," + formatDouble(topMassFraction) + ","
                    + tailSampleCount + "," + formatDouble(totalLogZhat) + ","
                    + formatDouble(selectedLogZhat) + "," + csv(outputFile.getAbsolutePath()));
        } catch (IOException e) {
            throw new RuntimeException("Failed to write S11 audit summary: " + summaryFile.getAbsolutePath(), e);
        }

        System.out.println("[S11_AUDIT] wrote " + selectedCount + "/" + records.size()
                + " records to " + outputFile.getAbsolutePath());
    }

    private List<Record> selectRecords() {
        List<Record> sorted = new ArrayList<>(records);
        sorted.sort(Comparator.comparingDouble((Record r) -> r.logZhat).reversed());
        totalLogZhat = logSumExp(sorted);
        for (int i = 0; i < sorted.size(); i++) {
            sorted.get(i).rank = i + 1;
        }

        if (mode == AuditMode.FULL) {
            for (Record rec : sorted) select(rec, "full", 1.0, 1.0);
            return sorted;
        }

        Set<Record> selected = new LinkedHashSet<>();
        List<Record> tail = selectTopMass(sorted, selected);
        if (mode == AuditMode.TOP_MASS_SAMPLE) {
            sampleTail(tail, selected);
        }

        return new ArrayList<>(selected);
    }

    private List<Record> selectTopMass(List<Record> sorted, Set<Record> selected) {
        List<Record> tail = new ArrayList<>();
        if (sorted.isEmpty()) return tail;
        double maxLog = sorted.get(0).logZhat;
        double total = 0.0;
        for (Record rec : sorted) total += Math.exp(rec.logZhat - maxLog);
        double target = topMassFraction * total;
        double cumulative = 0.0;
        boolean selectingTop = target > 0.0;
        for (Record rec : sorted) {
            if (selectingTop && cumulative < target) {
                select(rec, "top_mass", 1.0, 1.0);
                selected.add(rec);
                topSelectedCount++;
                cumulative += Math.exp(rec.logZhat - maxLog);
            } else {
                tail.add(rec);
            }
        }
        return tail;
    }

    private void sampleTail(List<Record> tail, Set<Record> selected) {
        if (tail.isEmpty() || tailSampleCount <= 0) return;
        double logTailZ = logSumExp(tail);
        Random rng = new Random(randomSeed ^ outputFile.getAbsolutePath().hashCode());
        for (Record rec : tail) {
            double expectedProb = tailSampleCount * Math.exp(rec.logZhat - logTailZ);
            double inclusionProb = Math.min(1.0, expectedProb);
            if (inclusionProb >= 1.0 || rng.nextDouble() < inclusionProb) {
                select(rec, "tail_sample", inclusionProb, 1.0 / inclusionProb);
                selected.add(rec);
                tailSampleSelectedCount++;
            }
        }
    }

    private void select(Record rec, String partition, double inclusionProb, double htWeight) {
        rec.auditPartition = partition;
        rec.auditInclusionProb = inclusionProb;
        rec.auditHtWeight = htWeight;
    }

    private String formatRecord(Record rec) {
        return csv(runId) + ","
                + csv(designId) + ","
                + csv(pdb) + ","
                + csv(state) + ","
                + sequenceIndex + ","
                + csv(sequence) + ","
                + csv(rec.confId) + ","
                + csv(rec.assignmentsText) + ","
                + formatDouble(rec.oldLower) + ","
                + formatDouble(rec.oldUpper) + ","
                + formatDouble(rec.gnnEnergy) + ","
                + formatDouble(rec.newLower) + ","
                + formatDouble(rec.newUpper) + ","
                + rec.rejectedAboveUpper + ","
                + formatDouble(rec.logZhat) + ","
                + rec.rank + ","
                + rec.eventIndex + ","
                + rec.nodeLevel + ","
                + rec.timestampMs + ","
                + csv(mode.name().toLowerCase()) + ","
                + csv(rec.auditPartition) + ","
                + formatDouble(rec.auditInclusionProb) + ","
                + formatDouble(rec.auditHtWeight);
    }

    private static double logSumExp(List<Record> records) {
        if (records.isEmpty()) return Double.NEGATIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (Record rec : records) max = Math.max(max, rec.logZhat);
        if (!Double.isFinite(max)) return max;
        double sum = 0.0;
        for (Record rec : records) sum += Math.exp(rec.logZhat - max);
        return max + Math.log(sum);
    }

    private static String assignmentsToString(int[] assignments) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < assignments.length; i++) {
            if (i > 0) sb.append(';');
            sb.append(assignments[i]);
        }
        return sb.toString();
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
        return String.format(java.util.Locale.US, "%.12g", value);
    }

    private static String sanitize(String value) {
        if (value == null || value.trim().isEmpty()) return "unknown";
        return value.replaceAll("[^A-Za-z0-9._-]", "_");
    }

    private static String emptyToUnknown(String value) {
        return value == null || value.trim().isEmpty() ? "unknown" : value;
    }

    private static String pdbName(String path) {
        if (path == null || path.trim().isEmpty()) return "unknown";
        return new File(path).getName();
    }

    private static double clamp01(double value) {
        if (Double.isNaN(value)) return 0.999;
        return Math.max(0.0, Math.min(1.0, value));
    }

    private static class Record {
        long eventIndex;
        int[] assignments;
        String assignmentsText;
        String confId;
        double oldLower;
        double oldUpper;
        double gnnEnergy;
        double newLower;
        double newUpper;
        boolean rejectedAboveUpper;
        double logZhat;
        int nodeLevel;
        long timestampMs;
        int rank;
        String auditPartition = "unselected";
        double auditInclusionProb = 0.0;
        double auditHtWeight = 0.0;
    }
}
