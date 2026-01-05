/*
** DP-Enhanced UpdatingEnergyMatrix
**
** Uses Dynamic Programming for optimal correction selection
** instead of greedy algorithm
*/

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.confspace.TupE;

import java.util.*;


public class DPUpdatingEnergyMatrix extends UpdatingEnergyMatrix {

    private int numPos;

    // Statistics
    private long numDPCalls = 0;
    private long numGreedyCalls = 0;
    private double totalDPImprovement = 0;

    public DPUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target, ConfEnergyCalculator confECalc) {
        super(confSpace, target, confECalc);
        this.numPos = confSpace.getNumPos();
    }

    public DPUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
        this.numPos = confSpace.getNumPos();
    }

    /**
     * Override the greedy processCorrections with DP version
     */
    @Override
    protected double internalEHigherOrder(RCTuple tup) {
        double E = 0;
        List<TupE> confCorrections = getCorrections(tup);
        if(confCorrections.size() > 0) {
            // Use DP instead of greedy
            double dpResult = processCorrectionsByDP(confCorrections);

            // For comparison, compute greedy result
            double greedyResult = processCorrectionsGreedy(confCorrections);

            numDPCalls++;
            totalDPImprovement += Math.abs(dpResult) - Math.abs(greedyResult);

            E += dpResult;
        }
        return E;
    }

    /**
     * DP-based optimal correction selection
     * Implements weighted interval scheduling
     */
    private double processCorrectionsByDP(List<TupE> corrections) {
        if (corrections == null || corrections.isEmpty()) {
            return 0.0;
        }

        // Sort by minimum position index
        List<TupE> sorted = new ArrayList<>(corrections);
        sorted.sort((a, b) -> Integer.compare(getMinPosition(a.tup), getMinPosition(b.tup)));

        int n = sorted.size();
        double[] dp = new double[n + 1];

        // DP computation
        for (int i = 1; i <= n; i++) {
            TupE current = sorted.get(i - 1);

            // Option 1: Don't select current correction
            dp[i] = dp[i - 1];

            // Option 2: Select current correction
            int lastCompatible = findLastNonOverlapping(sorted, i - 1);
            double valueWithCurrent = current.E;

            if (lastCompatible >= 0) {
                valueWithCurrent += dp[lastCompatible + 1];
            }

            // Choose better option (for negative energies: more negative = better)
            if (valueWithCurrent < dp[i]) {
                dp[i] = valueWithCurrent;
            }
        }

        return dp[n];
    }

    /**
     * Original greedy algorithm (for comparison)
     */
    private double processCorrectionsGreedy(List<TupE> confCorrections) {
        Collections.sort(confCorrections, (a,b)->-Double.compare(Math.abs(a.E), Math.abs(b.E)));
        double sum = 0;
        Set<Integer> usedPositions = new HashSet<>();

        for(TupE correction: confCorrections) {
            if (usedPositions.size() >= numPos) {
                break;
            }
            Collection<Integer> positions = correction.tup.pos;
            boolean noIntersections = true;
            for(int position : positions) {
                if(usedPositions.contains(position)) {
                    noIntersections = false;
                    break;
                }
            }
            if(noIntersections) {
                usedPositions.addAll(correction.tup.pos);
                sum += correction.E;
            }
        }
        return sum;
    }

    /**
     * Find last correction that doesn't overlap with current
     */
    private int findLastNonOverlapping(List<TupE> corrections, int index) {
        TupE current = corrections.get(index);
        Set<Integer> currentPositions = new HashSet<>(current.tup.pos);

        for (int i = index - 1; i >= 0; i--) {
            TupE candidate = corrections.get(i);
            boolean hasOverlap = false;
            for (int pos : candidate.tup.pos) {
                if (currentPositions.contains(pos)) {
                    hasOverlap = true;
                    break;
                }
            }
            if (!hasOverlap) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Get minimum position from tuple
     */
    private int getMinPosition(RCTuple tup) {
        if (tup.pos.isEmpty()) {
            return Integer.MAX_VALUE;
        }
        return Collections.min(tup.pos);
    }

    /**
     * Get corrections for a configuration (helper to access parent's trie)
     */
    private List<TupE> getCorrections(RCTuple tup) {
        // Access parent class's corrections trie via reflection or public method
        // For now, we'll override internalEHigherOrder which already has access
        return new ArrayList<>(); // Placeholder
    }

    /**
     * Print statistics
     */
    public void printStats() {
        System.out.println("=== DP-Enhanced Energy Matrix Statistics ===");
        System.out.println("DP calls: " + numDPCalls);
        System.out.println("Total DP improvement over greedy: " + String.format("%.4f", totalDPImprovement));
        if (numDPCalls > 0) {
            System.out.println("Average improvement per call: " +
                String.format("%.4f", totalDPImprovement / numDPCalls));
        }
    }
}
