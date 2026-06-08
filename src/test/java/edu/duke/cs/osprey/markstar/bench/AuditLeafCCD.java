package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
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
    static final int DEFAULT_GPU_STREAMS_PER_GPU = 64;
    private static final int DEFAULT_CPU_BATCH_FACTOR = 8;
    private static final int DEFAULT_GPU_BATCH_FACTOR = 16;

    public static void main(String[] args) throws Exception {

        String input = requireProperty("osprey.audit.input");
        String output = requireProperty("osprey.audit.output");
        String state = System.getProperty("osprey.audit.state", "complex").trim();
        String designId = System.getProperty("osprey.bench.designId", "unknown");
        String pdbPath = requireProperty("osprey.bench.pdbPath");
        String mutableStr = System.getProperty("osprey.bench.mutable", "");
        String flexibleStr = System.getProperty("osprey.bench.flexible", "");
        String device = normalizeAuditDevice(System.getProperty("osprey.audit.device", "cpu"));
        int cpus = Integer.getInteger("osprey.bench.numCPUs", 8);
        int numGpus = Math.max(1, Integer.getInteger("osprey.audit.numGpus", 1));
        int streamsPerGpu = Math.max(1, Integer.getInteger(
                "osprey.audit.streamsPerGpu", DEFAULT_GPU_STREAMS_PER_GPU));
        int shardIndex = Integer.getInteger("osprey.audit.shardIndex", 0);
        int numShards = Math.max(1, Integer.getInteger("osprey.audit.numShards", 1));
        EnergyCalculator.Type energyType = auditEnergyType(device);
        Parallelism parallelism = auditParallelism(device, cpus, numGpus, streamsPerGpu);
        int batchSize = Math.max(1, Integer.getInteger("osprey.audit.batchSize",
                defaultBatchSize(device, parallelism)));
        int warmupConfs = Math.max(0, Integer.getInteger("osprey.audit.warmupConfs",
                device.equals("gpu") ? parallelism.getParallelism() : 0));

        System.out.println("==============================================");
        System.out.println("  S11 Leaf CCD Audit");
        System.out.println("  Design: " + designId);
        System.out.println("  State: " + state);
        System.out.println("  Device: " + device + "  energyType=" + energyType
                + "  parallelism=" + parallelism.getParallelism()
                + (device.equals("gpu") ? "  gpus=" + numGpus + "  streams/gpu=" + streamsPerGpu : ""));
        System.out.println("  Input: " + input);
        System.out.println("  Output: " + output);
        System.out.println("  Shard: " + shardIndex + " / " + numShards);
        System.out.println("  Batch size: " + batchSize);
        System.out.println("  Warmup confs: " + warmupConfs);
        System.out.println("==============================================");

        AuditTiming timing = new AuditTiming();
        long confSpaceStartNs = System.nanoTime();
        TestKStar.ConfSpaces confSpaces = buildConfSpaces(pdbPath, mutableStr, flexibleStr);
        SimpleConfSpace confSpace = selectConfSpace(confSpaces, state);
        if (confSpace.positions.isEmpty()) {
            throw new IllegalArgumentException("Selected state has no positions: " + state);
        }
        timing.confSpaceNs = System.nanoTime() - confSpaceStartNs;

        long startedNs = System.nanoTime();
        int seen = 0;
        int audited = 0;
        int failed = 0;

        long ecalcBuildStartNs = System.nanoTime();
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setType(energyType)
                .build()) {
            timing.ecalcBuildNs = System.nanoTime() - ecalcBuildStartNs;

            long referenceStartNs = System.nanoTime();
            SimpleReferenceEnergies referenceEnergies = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                    .build().calcReferenceEnergies();
            timing.referenceNs = System.nanoTime() - referenceStartNs;
            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                    .setReferenceEnergies(referenceEnergies)
                    .build();
            System.out.println("  Actual energy type: " + ecalc.type);

            List<File> inputFiles = collectInputFiles(new File(input));
            if (warmupConfs > 0) {
                long warmupStartNs = System.nanoTime();
                BatchStats warmupStats = runWarmup(confEcalc, inputFiles, designId, state,
                        shardIndex, numShards, warmupConfs, timing);
                timing.warmupTotalNs = System.nanoTime() - warmupStartNs;
                System.out.println("  Warmup done: ok=" + warmupStats.ok
                        + ", failed=" + warmupStats.error
                        + ", elapsed=" + seconds(timing.warmupTotalNs));
            }

            File outFile = new File(output);
            File parent = outFile.getParentFile();
            if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
                throw new IOException("Failed to create output directory: " + parent.getAbsolutePath());
            }

            try (PrintWriter pw = new PrintWriter(new FileWriter(outFile))) {
                pw.println("run_id,design_id,pdb,state,sequence_index,sequence,conf_id,assignments,"
                        + "gnn_energy_kcal,ccd_energy_kcal,delta_kcal,log_zhat,log_ztrue,"
                        + "audit_time_ms,status,error,audit_partition,audit_inclusion_prob,audit_ht_weight");

                List<AuditWork> batch = new ArrayList<>(batchSize);
                for (File csv : inputFiles) {
                    try (BufferedReader br = new BufferedReader(new FileReader(csv))) {
                        String header = br.readLine();
                        if (header == null) continue;
                        List<String> columns = parseCsvLine(header);
                        Map<String, Integer> index = indexColumns(columns);

                        String line;
                        while ((line = br.readLine()) != null) {
                            if (line.trim().isEmpty()) continue;
                            long rowStartNs = System.nanoTime();
                            Map<String, String> row = rowMap(index, parseCsvLine(line));

                            if (!matches(row.get("design_id"), designId)) {
                                timing.inputNs += System.nanoTime() - rowStartNs;
                                continue;
                            }
                            if (!matchesState(row.get("state"), state)) {
                                timing.inputNs += System.nanoTime() - rowStartNs;
                                continue;
                            }

                            int rowIndex = seen++;
                            if (Math.floorMod(rowIndex, numShards) != shardIndex) {
                                timing.inputNs += System.nanoTime() - rowStartNs;
                                continue;
                            }

                            try {
                                int[] assignments = parseAssignments(row.get("assignments"));
                                double gnnEnergy = parseDouble(row.get("gnn_energy_kcal"));
                                AuditWork work = new AuditWork(row,
                                        new ConfSearch.ScoredConf(assignments, gnnEnergy), gnnEnergy);
                                timing.inputNs += System.nanoTime() - rowStartNs;
                                batch.add(work);
                                if (batch.size() >= batchSize) {
                                    BatchStats stats = flushBatch(confEcalc, pw, batch, timing);
                                    audited += stats.ok;
                                    failed += stats.error;
                                    batch.clear();
                                }
                            } catch (Exception e) {
                                double elapsedMs = (System.nanoTime() - rowStartNs) / 1e6;
                                timing.inputNs += System.nanoTime() - rowStartNs;
                                failed++;
                                long writeStartNs = System.nanoTime();
                                pw.println(formatOutput(row, Double.NaN, Double.NaN, elapsedMs,
                                        "error", e.getClass().getSimpleName() + ": " + e.getMessage()));
                                timing.writeNs += System.nanoTime() - writeStartNs;
                            }
                        }
                    }
                }
                BatchStats stats = flushBatch(confEcalc, pw, batch, timing);
                audited += stats.ok;
                failed += stats.error;
            }
        }

        timing.totalAfterConfSpaceNs = System.nanoTime() - startedNs;
        System.out.println("Audit done: seen=" + seen
                + ", audited=" + audited
                + ", failed=" + failed
                + ", elapsed=" + seconds(timing.totalAfterConfSpaceNs));
        printTimingSummary(timing);
    }

    static String normalizeAuditDevice(String device) {
        String normalized = device == null ? "cpu" : device.trim().toLowerCase(Locale.US);
        if (normalized.equals("cpu") || normalized.equals("gpu")) return normalized;
        throw new IllegalArgumentException("Unknown audit device: " + device + " (expected cpu or gpu)");
    }

    static Parallelism auditParallelism(String device, int cpus, int numGpus, int streamsPerGpu) {
        device = normalizeAuditDevice(device);
        if (device.equals("gpu")) {
            return Parallelism.make(Math.max(1, cpus), Math.max(1, numGpus), Math.max(1, streamsPerGpu));
        }
        return Parallelism.makeCpu(Math.max(1, cpus));
    }

    static int defaultBatchSize(String device, Parallelism parallelism) {
        device = normalizeAuditDevice(device);
        int factor = device.equals("gpu") ? DEFAULT_GPU_BATCH_FACTOR : DEFAULT_CPU_BATCH_FACTOR;
        return Math.max(1, parallelism.getParallelism() * factor);
    }

    static EnergyCalculator.Type auditEnergyType(String device) {
        device = normalizeAuditDevice(device);
        String explicit = System.getProperty("osprey.audit.energyType", "").trim();
        EnergyCalculator.Type type = explicit.isEmpty()
                ? (device.equals("gpu") ? EnergyCalculator.Type.ResidueCudaCCD : EnergyCalculator.Type.Cpu)
                : parseEnergyType(explicit);
        if (device.equals("gpu") && !isGpuEnergyType(type)) {
            throw new IllegalArgumentException("osprey.audit.device=gpu requires a CUDA energy type, got " + type);
        }
        if (!type.isSupported()) {
            throw new IllegalStateException("Requested audit energy type is not supported: " + type);
        }
        return type;
    }

    private static EnergyCalculator.Type parseEnergyType(String value) {
        for (EnergyCalculator.Type type : EnergyCalculator.Type.values()) {
            if (type.name().equalsIgnoreCase(value)) return type;
        }
        throw new IllegalArgumentException("Unknown audit energy type: " + value);
    }

    private static boolean isGpuEnergyType(EnergyCalculator.Type type) {
        return type == EnergyCalculator.Type.Cuda
                || type == EnergyCalculator.Type.ResidueCuda
                || type == EnergyCalculator.Type.CudaCCD
                || type == EnergyCalculator.Type.ResidueCudaCCD;
    }

    private static BatchStats flushBatch(ConfEnergyCalculator confEcalc, PrintWriter pw,
                                         List<AuditWork> batch, AuditTiming timing) {
        if (batch.isEmpty()) return new BatchStats(0, 0);

        long submitStartNs = System.nanoTime();
        for (AuditWork work : batch) {
            final AuditWork submitted = work;
            confEcalc.tasks.submit(() -> runAuditWork(confEcalc, submitted.scored), (result) -> {
                submitted.ccdEnergy = result.ccdEnergy;
                submitted.auditNs = result.auditNs;
                submitted.status = result.status;
                submitted.error = result.error;
            });
        }
        timing.submitNs += System.nanoTime() - submitStartNs;

        long waitStartNs = System.nanoTime();
        confEcalc.tasks.waitForFinish();
        timing.waitNs += System.nanoTime() - waitStartNs;

        int ok = 0;
        int error = 0;
        long writeStartNs = System.nanoTime();
        for (AuditWork work : batch) {
            if ("ok".equals(work.status)) ok++;
            else error++;
            if (work.auditNs >= 0) {
                timing.taskNs += work.auditNs;
                timing.maxTaskNs = Math.max(timing.maxTaskNs, work.auditNs);
                timing.taskCount++;
            }
            pw.println(formatOutput(work.row, work.gnnEnergy, work.ccdEnergy, work.auditTimeMs(),
                    work.status, work.error));
        }
        timing.writeNs += System.nanoTime() - writeStartNs;
        timing.batches++;
        timing.submitted += batch.size();
        return new BatchStats(ok, error);
    }

    private static BatchStats runWarmup(ConfEnergyCalculator confEcalc, List<File> inputFiles,
                                        String designId, String state, int shardIndex,
                                        int numShards, int maxConfs, AuditTiming timing) throws IOException {

        long loadStartNs = System.nanoTime();
        List<ConfSearch.ScoredConf> confs = loadWarmupConfs(inputFiles, designId, state,
                shardIndex, numShards, maxConfs);
        timing.warmupLoadNs += System.nanoTime() - loadStartNs;
        if (confs.isEmpty()) return new BatchStats(0, 0);

        AuditResult[] results = new AuditResult[confs.size()];
        long submitStartNs = System.nanoTime();
        for (int i = 0; i < confs.size(); i++) {
            final int idx = i;
            confEcalc.tasks.submit(() -> runAuditWork(confEcalc, confs.get(idx)),
                    (result) -> results[idx] = result);
        }
        timing.warmupSubmitNs += System.nanoTime() - submitStartNs;

        long waitStartNs = System.nanoTime();
        confEcalc.tasks.waitForFinish();
        timing.warmupWaitNs += System.nanoTime() - waitStartNs;

        int ok = 0;
        int error = 0;
        for (AuditResult result : results) {
            if (result != null && "ok".equals(result.status)) ok++;
            else error++;
            if (result != null) {
                timing.warmupTaskNs += result.auditNs;
                timing.warmupMaxTaskNs = Math.max(timing.warmupMaxTaskNs, result.auditNs);
                timing.warmupTaskCount++;
            }
        }
        return new BatchStats(ok, error);
    }

    private static List<ConfSearch.ScoredConf> loadWarmupConfs(List<File> inputFiles,
                                                                String designId,
                                                                String state,
                                                                int shardIndex,
                                                                int numShards,
                                                                int maxConfs) throws IOException {

        List<ConfSearch.ScoredConf> confs = new ArrayList<>(maxConfs);
        int seen = 0;
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

                    try {
                        int[] assignments = parseAssignments(row.get("assignments"));
                        double gnnEnergy = parseDouble(row.get("gnn_energy_kcal"));
                        confs.add(new ConfSearch.ScoredConf(assignments, gnnEnergy));
                    } catch (Exception e) {
                        // Ignore malformed warmup rows; the real audit path will report them in the CSV.
                    }

                    if (confs.size() >= maxConfs) return confs;
                }
            }
        }
        return confs;
    }

    private static AuditResult runAuditWork(ConfEnergyCalculator confEcalc, ConfSearch.ScoredConf scored) {
        long t0 = System.nanoTime();
        try {
            double ccdEnergy = confEcalc.calcEnergy(scored).getEnergy();
            return new AuditResult(ccdEnergy, System.nanoTime() - t0, "ok", "");
        } catch (Throwable t) {
            return new AuditResult(Double.NaN, System.nanoTime() - t0,
                    "error", t.getClass().getSimpleName() + ": " + t.getMessage());
        }
    }

    private static class AuditResult {
        final double ccdEnergy;
        final long auditNs;
        final String status;
        final String error;

        AuditResult(double ccdEnergy, long auditNs, String status, String error) {
            this.ccdEnergy = ccdEnergy;
            this.auditNs = auditNs;
            this.status = status;
            this.error = error;
        }
    }

    private static class BatchStats {
        final int ok;
        final int error;

        BatchStats(int ok, int error) {
            this.ok = ok;
            this.error = error;
        }
    }

    private static class AuditWork {
        final Map<String, String> row;
        final ConfSearch.ScoredConf scored;
        final double gnnEnergy;
        double ccdEnergy = Double.NaN;
        long auditNs = -1;
        String status = "error";
        String error = "audit task did not finish";

        AuditWork(Map<String, String> row, ConfSearch.ScoredConf scored, double gnnEnergy) {
            this.row = row;
            this.scored = scored;
            this.gnnEnergy = gnnEnergy;
        }

        double auditTimeMs() {
            return auditNs < 0 ? Double.NaN : auditNs / 1e6;
        }
    }

    private static class AuditTiming {
        long confSpaceNs;
        long totalAfterConfSpaceNs;
        long ecalcBuildNs;
        long referenceNs;
        long warmupTotalNs;
        long warmupLoadNs;
        long warmupSubmitNs;
        long warmupWaitNs;
        long warmupTaskNs;
        long warmupMaxTaskNs;
        int warmupTaskCount;
        long inputNs;
        long submitNs;
        long waitNs;
        long writeNs;
        long taskNs;
        long maxTaskNs;
        int taskCount;
        int batches;
        int submitted;
    }

    private static void printTimingSummary(AuditTiming timing) {
        System.out.println("Timing summary:");
        System.out.println("  confspace_setup_s=" + seconds(timing.confSpaceNs)
                + " total_after_confspace_s=" + seconds(timing.totalAfterConfSpaceNs)
                + " ecalc_build_s=" + seconds(timing.ecalcBuildNs)
                + " reference_energies_s=" + seconds(timing.referenceNs));
        if (timing.warmupTotalNs > 0) {
            System.out.println("  warmup_total_s=" + seconds(timing.warmupTotalNs)
                    + " warmup_load_s=" + seconds(timing.warmupLoadNs)
                    + " warmup_submit_s=" + seconds(timing.warmupSubmitNs)
                    + " warmup_wait_s=" + seconds(timing.warmupWaitNs)
                    + " warmup_task_sum_s=" + seconds(timing.warmupTaskNs)
                    + " warmup_task_avg_ms=" + ms(avgNs(timing.warmupTaskNs, timing.warmupTaskCount))
                    + " warmup_task_max_ms=" + ms(timing.warmupMaxTaskNs));
        }
        System.out.println("  audit_batches=" + timing.batches
                + " audit_submitted=" + timing.submitted
                + " input_parse_s=" + seconds(timing.inputNs)
                + " submit_s=" + seconds(timing.submitNs)
                + " wait_s=" + seconds(timing.waitNs)
                + " output_write_s=" + seconds(timing.writeNs)
                + " task_sum_s=" + seconds(timing.taskNs)
                + " task_avg_ms=" + ms(avgNs(timing.taskNs, timing.taskCount))
                + " task_max_ms=" + ms(timing.maxTaskNs));
    }

    private static long avgNs(long totalNs, int count) {
        return count <= 0 ? 0 : totalNs / count;
    }

    private static String seconds(long ns) {
        return String.format(Locale.US, "%.3f", ns / 1e9);
    }

    private static String ms(long ns) {
        return String.format(Locale.US, "%.3f", ns / 1e6);
    }

    static TestKStar.ConfSpaces buildConfSpaces(String pdbPath, String mutableStr,
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
            Strand.ResidueFlex resFlex = proteinStrand.flexibility.get(res);
            if (resFlex == null) {
                System.out.println("confspace.warning missing mutable residue after template matching: " + res);
                continue;
            }
            resFlex.setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
        }
        for (String res : flexibleResidues) {
            res = res.trim();
            if (res.isEmpty()) continue;
            if (!proteinChainSet.contains(String.valueOf(res.charAt(0)))) continue;
            Strand.ResidueFlex resFlex = proteinStrand.flexibility.get(res);
            if (resFlex == null) {
                System.out.println("confspace.warning missing protein flexible residue after template matching: " + res);
                continue;
            }
            resFlex.setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
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
                    Strand.ResidueFlex resFlex = ligStrand.flexibility.get(res);
                    if (resFlex == null) {
                        System.out.println("confspace.warning missing ligand flexible residue after template matching: " + res);
                        continue;
                    }
                    resFlex.setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
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

    static SimpleConfSpace selectConfSpace(TestKStar.ConfSpaces confSpaces, String state) {
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
                ? parseDoubleOrNaN(row.get("gnn_energy_kcal"))
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

    static List<File> collectInputFiles(File input) {
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

    static Map<String, Integer> indexColumns(List<String> columns) {
        Map<String, Integer> index = new HashMap<>();
        for (int i = 0; i < columns.size(); i++) index.put(columns.get(i), i);
        return index;
    }

    static Map<String, String> rowMap(Map<String, Integer> index, List<String> values) {
        Map<String, String> row = new HashMap<>();
        for (Map.Entry<String, Integer> entry : index.entrySet()) {
            int i = entry.getValue();
            row.put(entry.getKey(), i < values.size() ? values.get(i) : "");
        }
        return row;
    }

    static List<String> parseCsvLine(String line) {
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

    static int[] parseAssignments(String value) {
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

    private static double parseDoubleOrNaN(String value) {
        try {
            return parseDouble(value);
        } catch (RuntimeException e) {
            return Double.NaN;
        }
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
