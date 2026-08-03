package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarSequenceSharding;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Base64;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Whole-case PACK* production admission before any real DP table is built. */
public final class PackStarCasePreflight {

    private static final String LOG_PREFIX = "PACK* preflight:";
    private static final int DEFAULT_FINAL_MAX_STATES = 3;
    private static final long DEFAULT_FINAL_EXACT_MAX_MILLIS = 300_000L;
    private static final int DEFAULT_MAX_ROUNDS = 1;
    private static final int DEFAULT_PREVIEW_THREADS = 1;
    private static final String DISTRIBUTED_PROPERTY =
            "packstar.admission.distributed";
    private static final String SHARD_INDEX_PROPERTY =
            "packstar.admission.shardIndex";
    private static final String SHARD_COUNT_PROPERTY =
            "packstar.admission.shardCount";
    private static final String SHARD_SLA_HOURS_PROPERTY =
            "packstar.admission.shardSlaHours";
    private static final String FORCE_FINAL_PASS_PROPERTY =
            "packstar.admission.forceFinalPass";
    private static final String SHARD_RESULT_OUT_PROPERTY =
            "packstar.admission.shardResultOut";
    private static final String LOCKED_POLICY_ONLY_PROPERTY =
            "packstar.admission.lockedPolicyOnly";

    private PackStarCasePreflight() {
    }

    /**
     * Run preflight when caseSlaHours is positive; otherwise leave normal K*
     * construction unchanged and return {@code null}.
     */
    public static BranchDpAdmission.CaseSummary runIfConfigured(KStar kstar) {
        double slaHours = Math.max(0.0, PackStarConfig.getDouble(
                BranchDpAdmission.CASE_SLA_HOURS_PROPERTY, 0.0, LOG_PREFIX));
        if (slaHours <= 0.0) {
            return null;
        }

        requireAdaptiveDecomposition();
        int defaultShardCount = Math.max(1,
                environmentInteger("SLURM_NTASKS",
                        environmentInteger("SLURM_NPROCS", 1)));
        int shardCount = Math.max(1, PackStarConfig.getInteger(
                SHARD_COUNT_PROPERTY, defaultShardCount, LOG_PREFIX));
        int shardIndex = PackStarConfig.getInteger(
                SHARD_INDEX_PROPERTY,
                environmentInteger("SLURM_PROCID", 0), LOG_PREFIX);
        boolean distributed = shardCount > 1 || PackStarConfig.getBoolean(
                DISTRIBUTED_PROPERTY, false);
        if (!distributed) {
            shardCount = 1;
            shardIndex = 0;
        }
        if (shardIndex < 0 || shardIndex >= shardCount) {
            throw new IllegalArgumentException(LOG_PREFIX
                    + " invalid shard index=" + shardIndex
                    + " for shardCount=" + shardCount);
        }

        BranchDpAdmission.clearExactPolicies();
        String policyIn = PackStarConfig.getProperty(
                "packstar.admission.policyIn", "").trim();
        String policyOut = PackStarConfig.getProperty(
                "packstar.admission.policyOut", "").trim();
        boolean lockLoadedPolicies = false;
        if (!policyIn.isEmpty()) {
            int loaded = BranchDpAdmission.loadExactPolicies(new File(policyIn));
            lockLoadedPolicies = true;
            System.out.println(LOG_PREFIX + " loaded " + loaded
                    + " exact policies from " + policyIn
                    + " (lock mode: apply saved plan, no re-search)");
        }
        boolean lockedPolicyOnly = PackStarConfig.getBoolean(
                LOCKED_POLICY_ONLY_PROPERTY, false);
        if (lockedPolicyOnly) {
            if (!lockLoadedPolicies) {
                throw new IllegalStateException(LOG_PREFIX
                        + " lockedPolicyOnly requires packstar.admission.policyIn");
            }
            System.out.println(LOG_PREFIX
                    + " locked-policy-only mode: exact policies are loaded; "
                    + "preview and admission are skipped because the parent "
                    + "Slurm step already admitted this frozen plan.");
            return new BranchDpAdmission.CaseSummary(List.of(), slaHours);
        }
        List<Request> allRequests;
        List<Request> requests;
        List<Request> decisionRequests;
        int totalBundles;
        int assignedBundles;
        int replicatedBundleOrdinal;
        if (distributed) {
            // A distributed rank owns complete global-sequence bundles.  This
            // preserves K*'s protein -> ligand -> complex dependency and keeps
            // any short-circuit decision local to the node.  The distributed
            // path intentionally does not use the single-node filtered-state
            // de-duplication: a rank must be able to execute its sequences
            // independently, without a cross-rank pfunc cache.
            List<SequenceBundle> allBundles = collectSequenceBundles(kstar);
            boolean hasWildType = kstar.complex.confSpace.seqSpace()
                    .containsWildTypeSequence();
            List<SequenceBundle> selectedBundles = selectShard(
                    allBundles, hasWildType, shardIndex, shardCount);
            allRequests = flatten(allBundles);
            requests = flatten(selectedBundles);
            totalBundles = allBundles.size();
            assignedBundles = selectedBundles.size();
            replicatedBundleOrdinal = KStarSequenceSharding
                    .replicatedBundleOrdinal(hasWildType);
        } else {
            allRequests = collectUniqueRequests(kstar);
            if (allRequests.isEmpty()) {
                throw new IllegalStateException(LOG_PREFIX
                        + " no unique K* states were generated; refusing an empty production admission.");
            }
            requests = allRequests;
            totalBundles = 0;
            assignedBundles = 0;
            replicatedBundleOrdinal =
                    KStarSequenceSharding.NO_REPLICATED_BUNDLE;
        }
        // A sequence bundle is the scheduling unit, but a rank still avoids
        // previewing the same filtered state more than once.  This mirrors the
        // normal K* pfunc cache while preserving one output row per
        // (global-sequence, protein/ligand/complex) request for reduction.
        decisionRequests = distributed
                ? deduplicateRequests(requests) : requests;

        if (distributed) {
            System.out.println(LOG_PREFIX + " distributed mode: shard="
                    + shardIndex + "/" + shardCount
                    + " assignedSequences=" + assignedBundles
                    + " totalSequences=" + totalBundles
                    + " assignedStates=" + requests.size()
                    + " uniqueAssignedStates=" + decisionRequests.size()
                    + " totalStates=" + allRequests.size()
                    + " (sequence-bundle scale-out; each DP stays node-local)");
        }

        double decisionSlaHours = distributed
                ? distributedShardSlaHours()
                : slaHours;
        if (requests.isEmpty()) {
            BranchDpAdmission.CaseSummary empty =
                    new BranchDpAdmission.CaseSummary(List.of(), decisionSlaHours);
            writeShardResultIfConfigured(empty, requests, List.of(),
                    shardIndex, shardCount,
                    totalBundles, assignedBundles, allRequests.size(),
                    replicatedBundleOrdinal,
                    slaHours, decisionSlaHours, distributed);
            return empty;
        }

        System.out.println(LOG_PREFIX + " phase=initial uniqueStates="
                + decisionRequests.size() + ", assignedStates="
                + requests.size() + ", caseSlaHours="
                + String.format(Locale.ROOT, "%.4f", decisionSlaHours)
                + ", allocationPhase=before-materialization");

        int defaultFinalMaxStates = distributed
                ? Math.max(1, decisionRequests.size()) : DEFAULT_FINAL_MAX_STATES;
        int maxStates = Math.max(0, PackStarConfig.getInteger(
                BranchDpAdmission.FINAL_MAX_STATES_PROPERTY,
                defaultFinalMaxStates, LOG_PREFIX));
        long maxMillis = Math.max(1L, PackStarConfig.getLong(
                BranchDpAdmission.FINAL_EXACT_MAX_MILLIS_PROPERTY,
                DEFAULT_FINAL_EXACT_MAX_MILLIS, LOG_PREFIX));
        int maxRounds = Math.max(1, PackStarConfig.getInteger(
                "packstar.admission.maxRounds",
                distributed ? 4 : DEFAULT_MAX_ROUNDS, LOG_PREFIX));
        int defaultPreviewThreads = distributed
                ? Math.max(1, environmentInteger("SLURM_CPUS_PER_TASK",
                        Runtime.getRuntime().availableProcessors()))
                : DEFAULT_PREVIEW_THREADS;
        int previewThreads = Math.max(1, PackStarConfig.getInteger(
                "packstar.admission.previewThreads",
                defaultPreviewThreads, LOG_PREFIX));
        boolean forceFinalPass = distributed && PackStarConfig.getBoolean(
                FORCE_FINAL_PASS_PROPERTY, true);
        System.out.println(LOG_PREFIX + " search budget: finalMaxStates="
                + maxStates + ", finalExactMaxMillis=" + maxMillis
                + ", maxRounds=" + maxRounds
                + ", previewThreads=" + previewThreads
                + ", forceFinalPass=" + forceFinalPass);
        PackStarAdmissionDecision.Settings settings =
                new PackStarAdmissionDecision.Settings(
                        decisionSlaHours, maxStates, maxMillis, maxRounds,
                        previewThreads, lockLoadedPolicies, forceFinalPass,
                        BranchDecomposition.configuredExactImproveOptions());
        List<PackStarAdmissionDecision.State> states = new ArrayList<>();
        for (Request request : decisionRequests) {
            states.add(new PackStarAdmissionDecision.State(
                    request.sequence.toString(), () -> preview(request)));
        }

        BranchDpAdmission.CaseSummary summary;
        try (DryRunScope ignored = new DryRunScope()) {
            summary = PackStarAdmissionDecision.decide(states, settings)
                    .finalSummary;
        } finally {
            if (!policyOut.isEmpty()) {
                try {
                    int wrote = BranchDpAdmission.writeExactPolicies(
                            new File(policyOut));
                    System.out.println(LOG_PREFIX + " wrote " + wrote
                            + " exact policies to " + policyOut);
                } catch (RuntimeException ex) {
                    System.err.println(LOG_PREFIX
                            + " failed to write exact policies to " + policyOut
                            + ": " + ex.getMessage());
                }
            }
        }
        List<BranchDpAdmission.Prediction> outputPredictions = distributed
                ? expandPredictions(decisionRequests, requests, summary.predictions)
                : summary.predictions;
        writeShardResultIfConfigured(summary, requests, outputPredictions,
                shardIndex, shardCount,
                totalBundles, assignedBundles, allRequests.size(),
                replicatedBundleOrdinal,
                slaHours, decisionSlaHours, distributed);
        return summary;
    }

    private static int environmentInteger(String name, int defaultValue) {
        String value = System.getenv(name);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value.trim());
        } catch (NumberFormatException ignored) {
            return defaultValue;
        }
    }

    private static double distributedShardSlaHours() {
        double value = PackStarConfig.getDouble(SHARD_SLA_HOURS_PROPERTY,
                Double.MAX_VALUE, LOG_PREFIX);
        if (!(value > 0.0) || !Double.isFinite(value)) {
            return Double.MAX_VALUE;
        }
        return value;
    }

    private static List<SequenceBundle> selectShard(
            List<SequenceBundle> bundles, boolean hasWildType,
            int shardIndex, int shardCount) {
        if (shardCount <= 1) {
            return bundles;
        }
        List<SequenceBundle> selected = new ArrayList<>();
        for (SequenceBundle bundle : bundles) {
            if (KStarSequenceSharding.isAssignedToShard(
                    bundle.ordinal, hasWildType, shardIndex, shardCount)) {
                selected.add(bundle);
            }
        }
        return selected;
    }

    private static List<Request> flatten(List<SequenceBundle> bundles) {
        List<Request> requests = new ArrayList<>();
        for (SequenceBundle bundle : bundles) {
            requests.addAll(bundle.requests);
        }
        return requests;
    }

    private static List<Request> deduplicateRequests(List<Request> requests) {
        Map<KStar.ConfSpaceInfo, LinkedHashMap<Sequence, Request>> byState =
                new IdentityHashMap<>();
        for (Request request : requests) {
            byState.computeIfAbsent(request.info, ignored -> new LinkedHashMap<>())
                    .putIfAbsent(request.sequence, request);
        }
        List<Request> unique = new ArrayList<>();
        for (LinkedHashMap<Sequence, Request> states : byState.values()) {
            unique.addAll(states.values());
        }
        return unique;
    }

    private static List<BranchDpAdmission.Prediction> expandPredictions(
            List<Request> decisionRequests, List<Request> outputRequests,
            List<BranchDpAdmission.Prediction> predictions) {
        if (decisionRequests.size() != predictions.size()) {
            throw new IllegalStateException(LOG_PREFIX
                    + " decision prediction size mismatch: requests="
                    + decisionRequests.size() + ", predictions="
                    + predictions.size());
        }
        Map<KStar.ConfSpaceInfo, LinkedHashMap<Sequence,
                BranchDpAdmission.Prediction>> byState = new IdentityHashMap<>();
        for (int i = 0; i < decisionRequests.size(); i++) {
            Request request = decisionRequests.get(i);
            byState.computeIfAbsent(request.info,
                    ignored -> new LinkedHashMap<>())
                    .put(request.sequence, predictions.get(i));
        }
        List<BranchDpAdmission.Prediction> expanded = new ArrayList<>();
        for (Request request : outputRequests) {
            BranchDpAdmission.Prediction prediction = byState
                    .getOrDefault(request.info, new LinkedHashMap<>())
                    .get(request.sequence);
            if (prediction == null) {
                throw new IllegalStateException(LOG_PREFIX
                        + " missing prediction for sequence-bundle request="
                        + request.sequence);
            }
            expanded.add(prediction);
        }
        return expanded;
    }

    private static void writeShardResultIfConfigured(
            BranchDpAdmission.CaseSummary summary, List<Request> requests,
            List<BranchDpAdmission.Prediction> outputPredictions,
            int shardIndex, int shardCount, int totalBundles,
            int assignedBundles, int totalStates,
            int replicatedBundleOrdinal,
            double globalSlaHours, double localSlaHours, boolean distributed) {
        if (!distributed) {
            return;
        }
        String resultPath = PackStarConfig.getProperty(
                SHARD_RESULT_OUT_PROPERTY, "").trim();
        if (resultPath.isEmpty()) {
            System.out.println(LOG_PREFIX
                    + " distributed mode has no shardResultOut; "
                    + "predictions will not be reducible across nodes.");
            return;
        }
        if (outputPredictions.size() != requests.size()) {
            throw new IllegalStateException(LOG_PREFIX
                    + " shard result size mismatch: predictions="
                    + outputPredictions.size() + ", requests="
                    + requests.size());
        }
        File file = new File(resultPath);
        File parent = file.getAbsoluteFile().getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()
                && !parent.isDirectory()) {
            throw new IllegalStateException(LOG_PREFIX
                    + " cannot create shard-result directory " + parent);
        }
        try (BufferedWriter out = new BufferedWriter(new FileWriter(file))) {
            out.write("# packstar preflight shard v3");
            out.newLine();
            out.write("# shardIndex=" + shardIndex);
            out.newLine();
            out.write("# shardCount=" + shardCount);
            out.newLine();
            out.write("# totalStates=" + totalStates);
            out.newLine();
            out.write("# totalBundles=" + totalBundles);
            out.newLine();
            int statesPerBundle = totalBundles == 0
                    ? 0 : totalStates / totalBundles;
            out.write("# statesPerBundle=" + statesPerBundle);
            out.newLine();
            out.write("# assignedBundles=" + assignedBundles);
            out.newLine();
            out.write("# replicatedBundleOrdinal="
                    + replicatedBundleOrdinal);
            out.newLine();
            out.write("# globalSlaHours=" + Double.toString(globalSlaHours));
            out.newLine();
            out.write("# localSlaHours=" + Double.toString(localSlaHours));
            out.newLine();
            out.write("bundleOrdinal\tglobalSequenceB64\tordinal"
                    + "\tstateNameB64\tsequenceB64\tstateKeyB64"
                    + "\tbranchwidth\trootSplitEdge\tpredictedSeconds"
                    + "\tgpuWork\toocTrafficBytes\toocTrafficAvailable"
                    + "\tdpSweeps\tadaptiveAttempted\tadaptiveAccepted");
            out.newLine();
            for (int i = 0; i < outputPredictions.size(); i++) {
                BranchDpAdmission.Prediction prediction =
                        outputPredictions.get(i);
                Request request = requests.get(i);
                out.write(Integer.toString(request.bundleOrdinal));
                out.write('\t');
                out.write(encode(request.globalSequence));
                out.write('\t');
                out.write(Integer.toString(request.ordinal));
                out.write('\t');
                out.write(encode(prediction.stateName));
                out.write('\t');
                out.write(encode(request.sequence.toString()));
                out.write('\t');
                out.write(encode(prediction.stateKey));
                out.write('\t');
                out.write(Integer.toString(prediction.branchwidth));
                out.write('\t');
                out.write(Integer.toString(prediction.rootSplitEdge));
                out.write('\t');
                out.write(Double.toString(prediction.predictedSeconds));
                out.write('\t');
                out.write(prediction.gpuWork.toString());
                out.write('\t');
                out.write(prediction.oocTrafficBytes.toString());
                out.write('\t');
                out.write(Boolean.toString(prediction.oocTrafficAvailable));
                out.write('\t');
                out.write(Integer.toString(prediction.dpSweeps));
                out.write('\t');
                out.write(Boolean.toString(prediction.adaptiveAttempted));
                out.write('\t');
                out.write(Boolean.toString(prediction.adaptiveAccepted));
                out.newLine();
            }
        } catch (IOException ex) {
            throw new RuntimeException(LOG_PREFIX
                    + " failed to write shard result " + file, ex);
        }
        System.out.println(LOG_PREFIX + " wrote shard result " + file
                + " bundles=" + assignedBundles
                + " states=" + summary.predictions.size()
                + " localPredictedCaseHours=" + String.format(Locale.ROOT,
                "%.4f", summary.totalHours()));
    }

    private static String encode(String value) {
        return Base64.getEncoder().encodeToString(
                (value == null ? "" : value).getBytes(
                        java.nio.charset.StandardCharsets.UTF_8));
    }

    private static void requireAdaptiveDecomposition() {
        String strategy = PackStarConfig.getProperty(
                "branchdp.decomp.strategy", "hicks");
        if (!strategy.trim().equalsIgnoreCase("adaptive")
                && !strategy.trim().equalsIgnoreCase("production")
                && !strategy.trim().equalsIgnoreCase("adaptive_exact")
                && !strategy.trim().equalsIgnoreCase("adaptiveexact")) {
            throw new IllegalStateException(LOG_PREFIX
                    + " case SLA requires -Dpackstar.decomp.strategy=adaptive"
                    + " (or the branchdp alias); configured strategy=" + strategy);
        }
    }

    private static List<Request> collectUniqueRequests(KStar kstar) {
        List<Sequence> globalSequences = collectGlobalSequences(kstar);

        Map<KStar.ConfSpaceInfo, LinkedHashMap<Sequence, Request>> byState =
                new IdentityHashMap<>();
        for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
            byState.put(info, new LinkedHashMap<>());
        }
        for (Sequence global : globalSequences) {
            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                Sequence filtered = global.filter(info.confSpace.seqSpace());
                byState.get(info).computeIfAbsent(filtered,
                        seq -> new Request(info, seq, 0, global.toString()));
            }
        }

        List<Request> requests = new ArrayList<>();
        int ordinal = 0;
        for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
            for (Request request : byState.get(info).values()) {
                request.ordinal = ordinal++;
                requests.add(request);
            }
        }
        return requests;
    }

    private static List<SequenceBundle> collectSequenceBundles(KStar kstar) {
        List<Sequence> globalSequences = collectGlobalSequences(kstar);
        List<SequenceBundle> bundles = new ArrayList<>();
        int requestOrdinal = 0;
        for (int bundleOrdinal = 0; bundleOrdinal < globalSequences.size();
             bundleOrdinal++) {
            Sequence global = globalSequences.get(bundleOrdinal);
            List<Request> requests = new ArrayList<>();
            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                Sequence filtered = global.filter(info.confSpace.seqSpace());
                requests.add(new Request(info, filtered, bundleOrdinal,
                        global.toString(), requestOrdinal++));
            }
            bundles.add(new SequenceBundle(bundleOrdinal, global.toString(),
                    requests));
        }
        return bundles;
    }

    private static List<Sequence> collectGlobalSequences(KStar kstar) {
        List<Sequence> globalSequences = new ArrayList<>();
        if (kstar.complex.confSpace.seqSpace().containsWildTypeSequence()) {
            globalSequences.add(kstar.complex.confSpace.seqSpace()
                    .makeWildTypeSequence());
        }
        globalSequences.addAll(kstar.complex.confSpace.seqSpace().getMutants(
                kstar.settings.maxSimultaneousMutations, true));
        return globalSequences;
    }

    private static BranchDpAdmission.Prediction preview(Request request) {
        RCs rcs = request.sequence.makeRCs(request.info.confSpace);
        PartitionFunction pfunc = request.info.pfuncFactory.make(rcs);
        try {
            if (!(pfunc instanceof PackStarPartitionFunction)) {
                throw new IllegalStateException(LOG_PREFIX
                        + " expected PackStarPartitionFunction for state="
                        + request.info.type + ", got "
                        + pfunc.getClass().getName());
            }
            BranchDpAdmission.Prediction prediction =
                    ((PackStarPartitionFunction)pfunc)
                            .getAdmissionPrediction();
            if (prediction == null) {
                throw new IllegalStateException(LOG_PREFIX
                        + " backend produced no admission prediction for state="
                        + request.info.type);
            }
            return prediction;
        } finally {
            close(pfunc);
        }
    }

    private static void close(PartitionFunction pfunc) {
        if (!(pfunc instanceof AutoCloseable)) return;
        try {
            ((AutoCloseable)pfunc).close();
        } catch (RuntimeException ex) {
            throw ex;
        } catch (Exception ex) {
            throw new RuntimeException("failed to close preflight pfunc", ex);
        }
    }

    private static final class Request {

        final KStar.ConfSpaceInfo info;
        final Sequence sequence;
        final int bundleOrdinal;
        final String globalSequence;
        int ordinal;

        Request(KStar.ConfSpaceInfo info, Sequence sequence,
                int bundleOrdinal, String globalSequence) {
            this(info, sequence, bundleOrdinal, globalSequence, -1);
        }

        Request(KStar.ConfSpaceInfo info, Sequence sequence,
                int bundleOrdinal, String globalSequence, int ordinal) {
            this.info = info;
            this.sequence = sequence;
            this.bundleOrdinal = bundleOrdinal;
            this.globalSequence = globalSequence;
            this.ordinal = ordinal;
        }
    }

    private static final class SequenceBundle {

        final int ordinal;
        final String globalSequence;
        final List<Request> requests;

        SequenceBundle(int ordinal, String globalSequence,
                       List<Request> requests) {
            this.ordinal = ordinal;
            this.globalSequence = globalSequence;
            this.requests = List.copyOf(requests);
        }
    }

    private static final class DryRunScope implements AutoCloseable {

        private static final String BRANCH_KEY = "branchdp.dp.dryRun";
        private static final String PACK_KEY = "packstar.dp.dryRun";
        private final String oldBranch = System.getProperty(BRANCH_KEY);
        private final String oldPack = System.getProperty(PACK_KEY);

        DryRunScope() {
            System.setProperty(BRANCH_KEY, "true");
            System.setProperty(PACK_KEY, "true");
        }

        @Override
        public void close() {
            restore(BRANCH_KEY, oldBranch);
            restore(PACK_KEY, oldPack);
        }

        private static void restore(String key, String value) {
            if (value == null) {
                System.clearProperty(key);
            } else {
                System.setProperty(key, value);
            }
        }
    }
}
