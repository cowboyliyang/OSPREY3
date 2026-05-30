package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar.MARKStarResult;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

/**
 * Focused parameter sweep for BranchMARK*'s pending-edge lookahead.
 *
 * Iterates over (edge-selection rule × max-pending-edges) on a small wild-type
 * confspace and records runtime + K* bounds for each config. Reports the
 * fastest setting per edge-selection rule and the smallest pendingEdges that
 * comes within 5% of the best runtime.
 *
 * Deprecated parallel paths (internal-batch, enumeration) are pinned off.
 *
 * Driven by {@code slurm/scripts/run_parallelization_ablation.slurm}.
 *
 * System properties:
 *   osprey.branchmarkstar.numFlexible        / branchmarkstar.test.numFlexible  (default 10)
 *   osprey.branchmarkstar.epsilon            / branchmarkstar.test.epsilon      (default 0.68)
 *   osprey.branchmarkstar.pendingEdgeSweep   / branchmarkstar.pendingEdgeSweep  (default "2,3,4,5,6,8,12")
 *   osprey.branchmarkstar.edgeSelections     / branchmarkstar.edgeSelections    (default "contraction,contractionPerState")
 */
public class PendingEdgeSweep {

    private static class SweepRecord {
        final String edgeSelection;
        final int maxPendingEdges;
        final long runtimeMs;
        final String kstarBounds;

        SweepRecord(String edgeSelection, int maxPendingEdges, long runtimeMs, String kstarBounds) {
            this.edgeSelection = edgeSelection;
            this.maxPendingEdges = maxPendingEdges;
            this.runtimeMs = runtimeMs;
            this.kstarBounds = kstarBounds;
        }
    }

    @Test
    public void testPendingEdgeLookaheadSweep() {
        int numFlexible = Integer.getInteger("osprey.branchmarkstar.numFlexible",
                Integer.getInteger("branchmarkstar.test.numFlexible", 10));
        double epsilon = Double.parseDouble(System.getProperty("osprey.branchmarkstar.epsilon",
                System.getProperty("branchmarkstar.test.epsilon", "0.68")));
        int[] pendingEdgeValues = parseIntListProperty("osprey.branchmarkstar.pendingEdgeSweep",
                "branchmarkstar.pendingEdgeSweep", "2,3,4,5,6,8,12");
        String[] edgeSelections = parseStringListProperty("osprey.branchmarkstar.edgeSelections",
                "branchmarkstar.edgeSelections", "contraction,contractionPerState");

        System.out.println("\n========== PENDING EDGE LOOKAHEAD SWEEP ==========");
        System.out.println("Deprecated modes fixed off: parallelInternal=false, parallelEnumeration=false");
        System.out.println("parallelLookahead=true"
                + "   numFlexible=" + numFlexible
                + "   epsilon=" + epsilon
                + "   NUM_CPUs=" + TestBranchMARKStar.NUM_CPUs);
        System.out.println("pendingEdgeValues=" + Arrays.toString(pendingEdgeValues));
        System.out.println("edgeSelections=" + Arrays.toString(edgeSelections));

        List<SweepRecord> records = new ArrayList<>();

        for (String edgeSelection : edgeSelections) {
            for (int maxPendingEdges : pendingEdgeValues) {
                String label = "BranchMARK*-" + edgeSelection + "-pending" + maxPendingEdges;
                System.out.println("\n========== pending-edge sweep: edgeSelection="
                        + edgeSelection + " maxPendingEdges=" + maxPendingEdges + " ==========");

                long t0 = System.currentTimeMillis();
                MARKStarResult res = runBranchWithPendingEdgeConfig(numFlexible, epsilon,
                        edgeSelection, maxPendingEdges, label);
                long t1 = System.currentTimeMillis();

                String kstarBounds = "n/a";
                if (!res.scores.isEmpty()) {
                    KStarScore score = res.scores.get(0).score;
                    kstarBounds = String.format("[%.6f, %.6f]",
                            score.lowerBoundLog10(), score.upperBoundLog10());
                }

                long runtimeMs = t1 - t0;
                records.add(new SweepRecord(edgeSelection, maxPendingEdges, runtimeMs, kstarBounds));
                System.out.println(String.format(
                        "pendingEdgeSweep selection=%s pendingEdges=%d runtime=%d ms K*=%s",
                        edgeSelection, maxPendingEdges, runtimeMs, kstarBounds));
            }
        }

        System.out.println("\n========== PENDING EDGE LOOKAHEAD SWEEP SUMMARY ==========");
        for (SweepRecord record : records) {
            System.out.println(String.format("%-22s pendingEdges=%2d runtime=%7d ms   K*=%s",
                    record.edgeSelection, record.maxPendingEdges,
                    record.runtimeMs, record.kstarBounds));
        }

        for (String edgeSelection : edgeSelections) {
            long bestRuntime = Long.MAX_VALUE;
            int bestPendingEdges = -1;
            for (SweepRecord record : records) {
                if (record.edgeSelection.equals(edgeSelection) && record.runtimeMs < bestRuntime) {
                    bestRuntime = record.runtimeMs;
                    bestPendingEdges = record.maxPendingEdges;
                }
            }
            if (bestPendingEdges < 0) continue;

            long nearOptimalLimit = Math.round(bestRuntime * 1.05);
            int nearOptimalPendingEdges = bestPendingEdges;
            for (SweepRecord record : records) {
                if (record.edgeSelection.equals(edgeSelection)
                        && record.runtimeMs <= nearOptimalLimit
                        && record.maxPendingEdges < nearOptimalPendingEdges) {
                    nearOptimalPendingEdges = record.maxPendingEdges;
                }
            }

            System.out.println(String.format(
                    "pendingEdgeSweep best selection=%s pendingEdges=%d runtime=%d ms nearOptimalWithin5pct=%d",
                    edgeSelection, bestPendingEdges, bestRuntime, nearOptimalPendingEdges));
        }
        System.out.println("===========================================================");
    }

    private MARKStarResult runBranchWithPendingEdgeConfig(int numFlexible, double epsilon,
                                                          String edgeSelection,
                                                          int maxPendingEdges,
                                                          String label) {
        Map<String, String> backup = new HashMap<>();
        String[] keys = {
                "branchmarkstar.edgeSelection",
                "branchmarkstar.parallel.internal",
                "branchmarkstar.edgeSelection.parallelLookahead",
                "branchmarkstar.parallel.enumeration",
                "branchmarkstar.edgeSelection.maxPendingEdges",
        };
        for (String k : keys) backup.put(k, System.getProperty(k));

        System.setProperty("branchmarkstar.edgeSelection", edgeSelection);
        System.setProperty("branchmarkstar.parallel.internal", "false");
        System.setProperty("branchmarkstar.edgeSelection.parallelLookahead", "true");
        System.setProperty("branchmarkstar.parallel.enumeration", "false");
        System.setProperty("branchmarkstar.edgeSelection.maxPendingEdges",
                Integer.toString(maxPendingEdges));

        try {
            return TestBranchMARKStar.runMARKStarOnly(numFlexible, epsilon, true, label);
        } finally {
            for (Map.Entry<String, String> e : backup.entrySet()) {
                if (e.getValue() == null) System.clearProperty(e.getKey());
                else System.setProperty(e.getKey(), e.getValue());
            }
        }
    }

    private static int[] parseIntListProperty(String primaryKey, String fallbackKey,
                                              String defaultValue) {
        String value = System.getProperty(primaryKey,
                System.getProperty(fallbackKey, defaultValue));
        LinkedHashSet<Integer> parsed = new LinkedHashSet<>();
        for (String token : value.split("[,\\s]+")) {
            if (token.isEmpty()) continue;
            parsed.add(Integer.parseInt(token));
        }
        return parsed.stream().mapToInt(Integer::intValue).toArray();
    }

    private static String[] parseStringListProperty(String primaryKey, String fallbackKey,
                                                    String defaultValue) {
        String value = System.getProperty(primaryKey,
                System.getProperty(fallbackKey, defaultValue));
        LinkedHashSet<String> parsed = new LinkedHashSet<>();
        for (String token : value.split("[,\\s]+")) {
            String trimmed = token.trim();
            if (!trimmed.isEmpty()) parsed.add(trimmed);
        }
        return parsed.toArray(new String[0]);
    }
}
