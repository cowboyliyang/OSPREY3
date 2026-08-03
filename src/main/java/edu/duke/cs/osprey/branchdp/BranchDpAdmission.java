package edu.duke.cs.osprey.branchdp;

import edu.duke.cs.osprey.astar.conf.RCs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Allocation-free production-admission data for branch DP.
 *
 * <p>The backend produces one {@link Prediction} after it has exhaustively
 * previewed all requested roots, but before it materializes enumeration arrays
 * or DP tables. PACK* can therefore sum predictions for every unique K* state,
 * selectively ask expensive contributors for a stronger decomposition, and
 * reject the complete case without creating TiB-scale storage.</p>
 */
public final class BranchDpAdmission {

    public static final String GPU_WORK_PER_SECOND_PER_GPU_PROPERTY =
            "branchdp.admission.gpuWorkPerSecondPerGpu";
    public static final String GPU_COUNT_PROPERTY =
            "branchdp.admission.gpuCount";
    public static final String OOC_BYTES_PER_SECOND_PROPERTY =
            "branchdp.admission.oocBytesPerSecond";
    public static final String SAFETY_FACTOR_PROPERTY =
            "branchdp.admission.safetyFactor";
    public static final String SOFT_STATE_HOURS_PROPERTY =
            "branchdp.admission.softStateHours";
    public static final String DP_SWEEPS_PROPERTY =
            "branchdp.admission.dpSweeps";
    public static final String CASE_SLA_HOURS_PROPERTY =
            "branchdp.admission.caseSlaHours";
    public static final String FINAL_MAX_STATES_PROPERTY =
            "branchdp.admission.finalMaxStates";
    public static final String FINAL_EXACT_MAX_MILLIS_PROPERTY =
            "branchdp.admission.finalExactMaxMillis";

    private static final Map<String, ExactPolicy> EXACT_POLICIES =
            new ConcurrentHashMap<>();

    private BranchDpAdmission() {
    }

    /** Hardware conversion used for both root comparison and SLA accounting. */
    public static final class Hardware {

        public final double gpuWorkPerSecondPerGpu;
        public final int gpuCount;
        public final double oocBytesPerSecond;
        public final double safetyFactor;

        public Hardware(double gpuWorkPerSecondPerGpu, int gpuCount,
                        double oocBytesPerSecond) {
            this(gpuWorkPerSecondPerGpu, gpuCount, oocBytesPerSecond, 1.0);
        }

        public Hardware(double gpuWorkPerSecondPerGpu, int gpuCount,
                        double oocBytesPerSecond, double safetyFactor) {
            this.gpuWorkPerSecondPerGpu = positiveOrZero(gpuWorkPerSecondPerGpu);
            this.gpuCount = Math.max(0, gpuCount);
            this.oocBytesPerSecond = positiveOrZero(oocBytesPerSecond);
            this.safetyFactor = finiteAtLeastOne(safetyFactor);
        }

        public static Hardware fromBackendConfig() {
            return new Hardware(
                    BranchDpConfig.getBackendDouble(
                            GPU_WORK_PER_SECOND_PER_GPU_PROPERTY, 0.0,
                            BranchDpConfig.getBackendLogPrefix()),
                    BranchDpConfig.getBackendInteger(
                            GPU_COUNT_PROPERTY, 0,
                            BranchDpConfig.getBackendLogPrefix()),
                    BranchDpConfig.getBackendDouble(
                            OOC_BYTES_PER_SECOND_PROPERTY, 0.0,
                            BranchDpConfig.getBackendLogPrefix()),
                    BranchDpConfig.getBackendDouble(
                            SAFETY_FACTOR_PROPERTY, 1.0,
                            BranchDpConfig.getBackendLogPrefix()));
        }

        public boolean hasComputeRate() {
            return gpuWorkPerSecondPerGpu > 0.0 && gpuCount > 0;
        }

        public double computeSeconds(BigInteger gpuWork) {
            if (gpuWork == null || gpuWork.signum() <= 0) return 0.0;
            if (!hasComputeRate()) return Double.POSITIVE_INFINITY;
            return safeMultiply(safeRatio(gpuWork,
                    gpuWorkPerSecondPerGpu * gpuCount), safetyFactor);
        }

        public double oocSeconds(BigInteger trafficBytes, boolean trafficAvailable) {
            if (!trafficAvailable) return Double.POSITIVE_INFINITY;
            if (trafficBytes == null || trafficBytes.signum() <= 0) return 0.0;
            if (oocBytesPerSecond <= 0.0) return Double.POSITIVE_INFINITY;
            return safeMultiply(safeRatio(trafficBytes, oocBytesPerSecond),
                    safetyFactor);
        }

        public double totalSeconds(BigInteger gpuWork, BigInteger trafficBytes,
                                   boolean trafficAvailable,
                                   int gpuUnsupportedEdges) {
            if (gpuUnsupportedEdges > 0) return Double.POSITIVE_INFINITY;
            return safeAdd(computeSeconds(gpuWork),
                    oocSeconds(trafficBytes, trafficAvailable));
        }

        private static double positiveOrZero(double value) {
            return Double.isFinite(value) && value > 0.0 ? value : 0.0;
        }

        private static double finiteAtLeastOne(double value) {
            return Double.isFinite(value) && value >= 1.0 ? value : 1.0;
        }
    }

    /** One unique state/sequence estimate, before large DP storage exists. */
    public static final class Prediction {

        public final String stateName;
        public final String stateKey;
        public final int branchwidth;
        public final int rootSplitEdge;
        public final BigInteger gpuWork;
        public final BigInteger oocTrafficBytes;
        public final boolean oocTrafficAvailable;
        public final int gpuUnsupportedEdges;
        public final int dpSweeps;
        public final double gpuWorkPerSecondPerGpu;
        public final int gpuCount;
        public final double oocBytesPerSecond;
        public final double safetyFactor;
        public final double computeSeconds;
        public final double oocSeconds;
        public final double predictedSeconds;
        public final int firstStageBranchwidth;
        public final BigInteger firstStageGpuWork;
        public final double firstStagePredictedSeconds;
        public final boolean adaptiveAttempted;
        public final boolean adaptiveAccepted;

        public Prediction(String stateName, String stateKey,
                          int branchwidth, int rootSplitEdge,
                          BigInteger gpuWork, BigInteger oocTrafficBytes,
                          boolean oocTrafficAvailable, int gpuUnsupportedEdges,
                          Hardware hardware, int dpSweeps,
                          int firstStageBranchwidth,
                          BigInteger firstStageGpuWork,
                          double firstStagePredictedSeconds,
                          boolean adaptiveAttempted,
                          boolean adaptiveAccepted) {
            this.stateName = stateName == null ? "unknown" : stateName;
            this.stateKey = stateKey;
            this.branchwidth = branchwidth;
            this.rootSplitEdge = rootSplitEdge;
            this.gpuWork = gpuWork == null ? BigInteger.ZERO : gpuWork;
            this.oocTrafficBytes = oocTrafficBytes == null
                    ? BigInteger.ZERO : oocTrafficBytes;
            this.oocTrafficAvailable = oocTrafficAvailable;
            this.gpuUnsupportedEdges = Math.max(0, gpuUnsupportedEdges);
            this.dpSweeps = Math.max(1, dpSweeps);
            this.gpuWorkPerSecondPerGpu = hardware.gpuWorkPerSecondPerGpu;
            this.gpuCount = hardware.gpuCount;
            this.oocBytesPerSecond = hardware.oocBytesPerSecond;
            this.safetyFactor = hardware.safetyFactor;
            BigInteger totalGpuWork = totalGpuWork();
            BigInteger totalOocTraffic = totalOocTrafficBytes();
            this.computeSeconds = this.gpuUnsupportedEdges == 0
                    ? hardware.computeSeconds(totalGpuWork)
                    : Double.POSITIVE_INFINITY;
            this.oocSeconds = this.gpuUnsupportedEdges == 0
                    ? hardware.oocSeconds(totalOocTraffic,
                    this.oocTrafficAvailable)
                    : Double.POSITIVE_INFINITY;
            this.predictedSeconds = safeAdd(computeSeconds, oocSeconds);
            this.firstStageBranchwidth = firstStageBranchwidth;
            this.firstStageGpuWork = firstStageGpuWork == null
                    ? this.gpuWork : firstStageGpuWork;
            this.firstStagePredictedSeconds = firstStagePredictedSeconds;
            this.adaptiveAttempted = adaptiveAttempted;
            this.adaptiveAccepted = adaptiveAccepted;
        }

        public double predictedHours() {
            return predictedSeconds / 3600.0;
        }

        public BigInteger totalGpuWork() {
            return gpuWork.multiply(BigInteger.valueOf(dpSweeps));
        }

        public BigInteger totalOocTrafficBytes() {
            return oocTrafficBytes.multiply(BigInteger.valueOf(dpSweeps));
        }

        public double computeHours() {
            return computeSeconds / 3600.0;
        }

        public double oocHours() {
            return oocSeconds / 3600.0;
        }

        public int acceptedBranchwidthDrop() {
            return adaptiveAccepted
                    ? Math.max(0, firstStageBranchwidth - branchwidth)
                    : 0;
        }

        public boolean hasFinitePrediction() {
            return Double.isFinite(predictedSeconds);
        }

        public String format() {
            return String.format(Locale.ROOT,
                    "state=%s bw=%d root=%d gpuWork=%s oocBytes=%s "
                            + "dpSweeps=%d totalGpuWork=%s totalOocBytes=%s "
                            + "gpuRatePerGpu=%.3f gpuCount=%d oocRate=%.3f safetyFactor=%.3f "
                            + "computeHours=%.4f oocHours=%.4f predictedHours=%.4f "
                            + "adaptiveAttempted=%s adaptiveAccepted=%s",
                    stateName, branchwidth, rootSplitEdge, gpuWork,
                    oocTrafficAvailable ? oocTrafficBytes.toString() : "unavailable",
                    dpSweeps, totalGpuWork(),
                    oocTrafficAvailable ? totalOocTrafficBytes().toString() : "unavailable",
                    gpuWorkPerSecondPerGpu, gpuCount, oocBytesPerSecond,
                    safetyFactor,
                    computeHours(), oocHours(), predictedHours(),
                    adaptiveAttempted, adaptiveAccepted);
        }
    }

    /** Immutable whole-case sum. Each prediction must represent one unique state. */
    public static final class CaseSummary {

        public final List<Prediction> predictions;
        public final double totalSeconds;
        public final double slaHours;

        public CaseSummary(List<Prediction> predictions, double slaHours) {
            this.predictions = Collections.unmodifiableList(
                    new ArrayList<>(predictions));
            double total = 0.0;
            for (Prediction prediction : predictions) {
                total = safeAdd(total, prediction.predictedSeconds);
            }
            this.totalSeconds = total;
            this.slaHours = slaHours;
        }

        public double totalHours() {
            return totalSeconds / 3600.0;
        }

        public boolean withinSla() {
            return slaHours <= 0.0
                    || (Double.isFinite(totalSeconds)
                    && totalHours() <= slaHours);
        }

        public List<Prediction> contributorsDescending() {
            List<Prediction> sorted = new ArrayList<>(predictions);
            sorted.sort(Comparator.comparingDouble(
                    (Prediction p) -> p.predictedSeconds).reversed());
            return sorted;
        }
    }

    /** Per-state bounded exact retry retained from preflight into the real run. */
    public static final class ExactPolicy {

        public final int minDrop;
        public final int maxDrop;
        public final long maxMillis;
        public final double maxPredictedSeconds;

        public ExactPolicy(int minDrop, int maxDrop, long maxMillis) {
            this(minDrop, maxDrop, maxMillis, Double.POSITIVE_INFINITY);
        }

        public ExactPolicy(int minDrop, int maxDrop, long maxMillis,
                           double maxPredictedSeconds) {
            this.minDrop = Math.max(1, minDrop);
            this.maxDrop = Math.max(this.minDrop, maxDrop);
            this.maxMillis = Math.max(1L, maxMillis);
            this.maxPredictedSeconds = maxPredictedSeconds > 0.0
                    ? maxPredictedSeconds : Double.POSITIVE_INFINITY;
        }

        public ExactPolicy withPredictionCeiling(double seconds) {
            return new ExactPolicy(minDrop, maxDrop, maxMillis, seconds);
        }
    }

    /** Exact RC identity; avoids applying a final-pass decision to another mutant. */
    public static String stateKey(String stateName, RCs rcs) {
        StringBuilder key = new StringBuilder(
                stateName == null ? "unknown" : stateName);
        key.append('|').append(rcs.getNumPos());
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            key.append('|').append(pos).append(':');
            int[] allowed = rcs.get(pos);
            for (int rc : allowed) {
                key.append(rc).append(',');
            }
        }
        return key.toString();
    }

    public static ExactPolicy getExactPolicy(String stateKey) {
        return stateKey == null ? null : EXACT_POLICIES.get(stateKey);
    }

    public static void putExactPolicy(String stateKey, ExactPolicy policy) {
        if (stateKey == null) throw new IllegalArgumentException("stateKey is required");
        if (policy == null) {
            EXACT_POLICIES.remove(stateKey);
        } else {
            EXACT_POLICIES.put(stateKey, policy);
        }
    }

    public static void clearExactPolicies() {
        EXACT_POLICIES.clear();
    }

    public static int exactPolicyCount() {
        return EXACT_POLICIES.size();
    }

    /**
     * Serialize the retained exact policies to a durable file so a standalone
     * preflight search can hand its winning branchwidth/root plan to a later
     * production JVM.  Format (tab-separated, stateKey last so it may contain
     * '|', ':', ','):  minDrop maxDrop maxMillis maxPredictedSeconds stateKey
     */
    public static int writeExactPolicies(File file) {
        if (file == null) {
            throw new IllegalArgumentException("policy file is required");
        }
        File parent = file.getAbsoluteFile().getParentFile();
        if (parent != null) {
            parent.mkdirs();
        }
        int written = 0;
        try (BufferedWriter out = new BufferedWriter(new FileWriter(file))) {
            out.write("# packstar exact-policy dump v1: "
                    + "minDrop\tmaxDrop\tmaxMillis\tmaxPredictedSeconds\tstateKey");
            out.newLine();
            List<Map.Entry<String, ExactPolicy>> entries =
                    new ArrayList<>(EXACT_POLICIES.entrySet());
            entries.sort(Map.Entry.comparingByKey());
            for (Map.Entry<String, ExactPolicy> entry : entries) {
                ExactPolicy p = entry.getValue();
                out.write(p.minDrop + "\t" + p.maxDrop + "\t" + p.maxMillis
                        + "\t" + Double.toString(p.maxPredictedSeconds)
                        + "\t" + entry.getKey());
                out.newLine();
                written++;
            }
        } catch (IOException ex) {
            throw new UncheckedIOException(
                    "failed to write exact-policy file " + file, ex);
        }
        return written;
    }

    /**
     * Load exact policies previously written by {@link #writeExactPolicies}.
     * Existing entries with the same state key are overwritten.
     */
    public static int loadExactPolicies(File file) {
        if (file == null) {
            throw new IllegalArgumentException("policy file is required");
        }
        int loaded = 0;
        try (BufferedReader in = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }
                String[] f = line.split("\t", 5);
                if (f.length < 5) {
                    throw new IllegalStateException(
                            "malformed exact-policy line: " + line);
                }
                EXACT_POLICIES.put(f[4], new ExactPolicy(
                        Integer.parseInt(f[0].trim()),
                        Integer.parseInt(f[1].trim()),
                        Long.parseLong(f[2].trim()),
                        Double.parseDouble(f[3].trim())));
                loaded++;
            }
        } catch (IOException ex) {
            throw new UncheckedIOException(
                    "failed to read exact-policy file " + file, ex);
        }
        return loaded;
    }

    /**
     * Require the formal allocation-free preview to reproduce an exact plan
     * retained by whole-case preflight before any DP table is materialized.
     */
    public static void enforceRetainedPredictionCeiling(
            Prediction prediction, String logPrefix) {
        if (prediction == null) return;
        ExactPolicy policy = getExactPolicy(prediction.stateKey);
        if (policy == null || !Double.isFinite(policy.maxPredictedSeconds)) {
            return;
        }
        double tolerance = Math.max(1.0,
                policy.maxPredictedSeconds * 1.0e-9);
        if (!Double.isFinite(prediction.predictedSeconds)
                || prediction.predictedSeconds
                > policy.maxPredictedSeconds + tolerance) {
            String prefix = logPrefix == null ? "Branch-DP:" : logPrefix;
            throw new IllegalStateException(prefix
                    + " retained whole-case decomposition could not be reproduced for state="
                    + prediction.stateName + ": predictedHours="
                    + String.format(Locale.ROOT, "%.4f",
                    prediction.predictedHours())
                    + " exceeds preflight ceilingHours="
                    + String.format(Locale.ROOT, "%.4f",
                    policy.maxPredictedSeconds / 3600.0)
                    + ". This check ran before DP-table materialization.");
        }
    }

    private static double safeRatio(BigInteger numerator, double denominator) {
        if (!(denominator > 0.0) || !Double.isFinite(denominator)) {
            return Double.POSITIVE_INFINITY;
        }
        double value = numerator.doubleValue() / denominator;
        return Double.isFinite(value) && value >= 0.0
                ? value : Double.POSITIVE_INFINITY;
    }

    private static double safeAdd(double a, double b) {
        if (!Double.isFinite(a) || !Double.isFinite(b)) {
            return Double.POSITIVE_INFINITY;
        }
        double sum = a + b;
        return Double.isFinite(sum) && sum >= 0.0
                ? sum : Double.POSITIVE_INFINITY;
    }

    private static double safeMultiply(double value, double multiplier) {
        if (!Double.isFinite(value) || !Double.isFinite(multiplier)
                || value < 0.0 || multiplier < 0.0) {
            return Double.POSITIVE_INFINITY;
        }
        double product = value * multiplier;
        return Double.isFinite(product) && product >= 0.0
                ? product : Double.POSITIVE_INFINITY;
    }
}
