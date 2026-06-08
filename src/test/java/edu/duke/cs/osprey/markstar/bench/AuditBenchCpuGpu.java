package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * CPU-vs-GPU CCD throughput micro-benchmark for the S11 leaf audit.
 *
 * Re-minimizes the confs recorded in an audit-leaf CSV in a parallel batch
 * (same calcEnergyAsync + tasks.waitForFinish path MARK* uses), and reports
 * the ONE-TIME setup overhead (confspace build + reference energies + ecalc /
 * GPU-context init -- which a large case amortizes) SEPARATELY from the
 * steady-state per-conf throughput. Warmup (first GPU pass = PTX JIT) is
 * excluded by running several timed passes on the same ecalc and reading off
 * the median of passes >= 1.
 *
 * Properties:
 *   osprey.audit.device        cpu | gpu          (default cpu)
 *   osprey.audit.input         CSV file or dir of audit-leaf CSVs
 *   osprey.audit.state         Protein|Ligand|Complex (default complex)
 *   osprey.bench.designId      design id filter   (default unknown = no filter)
 *   osprey.bench.pdbPath       prepped PDB
 *   osprey.bench.mutable       ';'-separated mutable residues
 *   osprey.bench.flexible      ';'-separated flexible residues
 *   osprey.bench.numCPUs       CPU threads        (default 8)
 *   osprey.audit.numGpus       GPUs               (default 1)
 *   osprey.audit.streamsPerGpu GPU streams        (default 64)
 *   osprey.audit.batchSize     submitted confs before waitForFinish
 *                              (default 8x CPU parallelism, 16x GPU parallelism)
 *   osprey.audit.maxConfs      cap N              (default all)
 *   osprey.audit.repeats       timed passes       (default 5)
 */
public class AuditBenchCpuGpu {

    public static void main(String[] args) throws Exception {
        String input    = req("osprey.audit.input");
        String state    = System.getProperty("osprey.audit.state", "complex").trim();
        String designId = System.getProperty("osprey.bench.designId", "unknown");
        String pdbPath  = req("osprey.bench.pdbPath");
        String mutable  = System.getProperty("osprey.bench.mutable", "");
        String flexible = System.getProperty("osprey.bench.flexible", "");
        String device   = System.getProperty("osprey.audit.device", "cpu").trim().toLowerCase(Locale.US);
        int cpus     = Integer.getInteger("osprey.bench.numCPUs", 8);
        int numGpus  = Integer.getInteger("osprey.audit.numGpus", 1);
        int streams  = Integer.getInteger("osprey.audit.streamsPerGpu",
                AuditLeafCCD.DEFAULT_GPU_STREAMS_PER_GPU);
        int maxConfs = Integer.getInteger("osprey.audit.maxConfs", 1_000_000);
        int repeats  = Integer.getInteger("osprey.audit.repeats", 5);
        EnergyCalculator.Type energyType = AuditLeafCCD.auditEnergyType(device);
        Parallelism par = AuditLeafCCD.auditParallelism(device, cpus, numGpus, streams);
        int batchSize = Math.max(1, Integer.getInteger("osprey.audit.batchSize",
                AuditLeafCCD.defaultBatchSize(device, par)));

        System.out.println("================ AuditBenchCpuGpu ================");
        System.out.println("  device=" + device + " cpus=" + cpus
                + (device.equals("gpu") ? "  gpus=" + numGpus + "  streams/gpu=" + streams : "")
                + "  energyType=" + energyType);
        System.out.println("  parallelism=" + par.getParallelism() + "  batchSize=" + batchSize);
        System.out.println("  design=" + designId + "  state=" + state);
        System.out.println("  input=" + input);

        TestKStar.ConfSpaces confSpaces = AuditLeafCCD.buildConfSpaces(pdbPath, mutable, flexible);
        SimpleConfSpace confSpace = AuditLeafCCD.selectConfSpace(confSpaces, state);
        if (confSpace.positions.isEmpty()) {
            throw new IllegalArgumentException("Selected state has no positions: " + state);
        }

        List<ConfSearch.ScoredConf> confs = loadConfs(input, designId, state, maxConfs);
        int n = confs.size();
        System.out.println("  confs=" + n + "  flexible positions=" + confSpace.positions.size());
        if (n == 0) throw new IllegalArgumentException("No confs matched design/state");

        long tSetup0 = System.nanoTime();
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, confSpaces.ffparams)
                .setParallelism(par)
                .setType(energyType)
                .build()) {

            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                            .build().calcReferenceEnergies())
                    .build();
            double setupMs = (System.nanoTime() - tSetup0) / 1e6;

            System.out.println("  ecalc.type=" + ecalc.type);
            System.out.printf(Locale.US,
                    "  SETUP_MS=%.1f   (one-time: confspace + reference-energy emat + ecalc/GPU init; amortizes in big cases)%n",
                    setupMs);

            double[] energies = new double[n];
            double[] passMs = new double[repeats];

            for (int r = 0; r < repeats; r++) {
                Arrays.fill(energies, 0.0);
                long t0 = System.nanoTime();
                for (int start = 0; start < n; start += batchSize) {
                    int end = Math.min(n, start + batchSize);
                    for (int i = start; i < end; i++) {
                        final int idx = i;
                        confEcalc.calcEnergyAsync(confs.get(i), (econf) -> energies[idx] = econf.getEnergy());
                    }
                    confEcalc.tasks.waitForFinish();
                }
                passMs[r] = (System.nanoTime() - t0) / 1e6;

                double checksum = 0;
                for (double e : energies) checksum += e;
                System.out.printf(Locale.US,
                        "  PASS %d  phase_ms=%.1f  per_conf_ms=%.4f  throughput=%.1f conf/s  checksum=%.4f%s%n",
                        r, passMs[r], passMs[r] / n, n / (passMs[r] / 1000.0), checksum,
                        (r == 0 ? "   <- includes warmup/JIT" : ""));
            }

            // steady state = median of passes >= 1 (warmup excluded); fall back to pass 0 if repeats==1
            double steadyMs;
            if (repeats >= 2) {
                double[] tail = Arrays.copyOfRange(passMs, 1, repeats);
                Arrays.sort(tail);
                steadyMs = tail[tail.length / 2];
            } else {
                steadyMs = passMs[0];
            }
            System.out.println("  -------- SUMMARY (overhead excluded) --------");
            System.out.printf(Locale.US,
                    "  RESULT device=%s streams=%d batchSize=%d n=%d steady_phase_ms=%.1f steady_per_conf_ms=%.4f steady_throughput=%.1f conf/s setup_ms=%.1f%n",
                    device, streams, batchSize, n, steadyMs, steadyMs / n, n / (steadyMs / 1000.0), setupMs);
        }
    }

    private static List<ConfSearch.ScoredConf> loadConfs(String input, String designId, String state, int maxConfs)
            throws Exception {
        List<ConfSearch.ScoredConf> confs = new ArrayList<>();
        for (File csv : AuditLeafCCD.collectInputFiles(new File(input))) {
            try (BufferedReader br = new BufferedReader(new FileReader(csv))) {
                String header = br.readLine();
                if (header == null) continue;
                Map<String, Integer> index = AuditLeafCCD.indexColumns(AuditLeafCCD.parseCsvLine(header));
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.trim().isEmpty()) continue;
                    Map<String, String> row = AuditLeafCCD.rowMap(index, AuditLeafCCD.parseCsvLine(line));
                    if (!matches(row.get("design_id"), designId)) continue;
                    if (!norm(row.get("state")).equals(norm(state))) continue;
                    int[] a = AuditLeafCCD.parseAssignments(row.get("assignments"));
                    double g = parseD(row.get("gnn_energy_kcal"));
                    confs.add(new ConfSearch.ScoredConf(a, g));
                    if (confs.size() >= maxConfs) return confs;
                }
            }
        }
        return confs;
    }

    private static boolean matches(String v, String expected) {
        return expected == null || expected.equals("unknown") || expected.trim().isEmpty() || expected.equals(v);
    }

    private static String norm(String s) { return s == null ? "" : s.trim().toLowerCase(Locale.US); }

    private static double parseD(String s) {
        return (s == null || s.trim().isEmpty()) ? Double.NaN : Double.parseDouble(s.trim());
    }

    private static String req(String k) {
        String v = System.getProperty(k);
        if (v == null || v.trim().isEmpty()) throw new IllegalArgumentException("Missing required property: " + k);
        return v;
    }
}
