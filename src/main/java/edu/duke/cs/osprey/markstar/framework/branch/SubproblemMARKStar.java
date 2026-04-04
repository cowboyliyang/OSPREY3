/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
**
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;

/**
 * DEPRECATED — Superseded by incremental lambda enumeration in RootedTreeEdge.
 *
 * This class implemented per-edge A* expansion over lambda positions with partial
 * (M+lambda) minimization. It has been replaced by a two-phase approach:
 *   Phase 1: Incremental lambda enumeration in RootedTreeEdge (analogous to A* but
 *            using sorted lambda-states with tail bounds instead of H-score).
 *   Phase 2: Full-conformation minimization in BranchMARKStarBound with Z correction.
 *
 * The new approach avoids partial minimization (which violates minimization indivisibility)
 * and supports lazy evaluation (only enumerate high-contribution lambda-states).
 *
 * Kept for reference. Not used by BranchMARKStarBound.
 *
 * @deprecated Use RootedTreeEdge.enumerateNextLambda() instead.
 */
@Deprecated
public class SubproblemMARKStar {

    private final Set<Integer> lambdaPositions;
    private final Set<Integer> mSetPositions;
    private final EnergyMatrix globalRigidEmat;
    private final EnergyMatrix globalMinEmat;
    private final RCs globalRCs;
    private final int numLambdaPos;

    private final Map<MSetAssignment, SubproblemBound> boundsByMSet = new HashMap<>();

    private static final BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    /**
     * A leaf conformation pending continuous minimization.
     */
    public static class PendingMinimization {
        public final MSetAssignment mSetAssignment;
        public final int[] localConf;     // local RC indices for lambda positions
        public final int[] globalConf;    // full global conformation (all positions)
        public final double rigidEnergy;
        public final double minEnergy;

        PendingMinimization(MSetAssignment mSetAssignment, int[] localConf,
                            int[] globalConf, double rigidEnergy, double minEnergy) {
            this.mSetAssignment = mSetAssignment;
            this.localConf = localConf;
            this.globalConf = globalConf;
            this.rigidEnergy = rigidEnergy;
            this.minEnergy = minEnergy;
        }
    }

    /**
     * Bound state for a single M-set assignment.
     */
    public static class SubproblemBound {
        SubproblemEnergyMatrix rigidEmat;
        SubproblemEnergyMatrix minEmat;
        PriorityQueue<SearchNode> queue;
        BigDecimal lowerBound;
        BigDecimal upperBound;
        int numExpansions;
        List<LeafConformation> leafConformations;
        // New leaves found since last check (for minimization)
        List<LeafConformation> pendingLeaves;

        public SubproblemBound() {
            queue = new PriorityQueue<>();
            lowerBound = BigDecimal.ZERO;
            upperBound = BigDecimal.ZERO;
            numExpansions = 0;
            leafConformations = new ArrayList<>();
            pendingLeaves = new ArrayList<>();
        }
    }

    /**
     * A fully assigned leaf conformation with energy bounds.
     * After minimization, trueEnergy is set and bounds tighten.
     */
    public static class LeafConformation {
        int[] localConf;
        double rigidEnergy;
        double minEnergy;
        double trueEnergy;  // NaN until minimized
        boolean minimized;

        LeafConformation(int[] localConf, double rigidEnergy, double minEnergy) {
            this.localConf = localConf;
            this.rigidEnergy = rigidEnergy;
            this.minEnergy = minEnergy;
            this.trueEnergy = Double.NaN;
            this.minimized = false;
        }
    }

    /**
     * A node in the subproblem A* search tree.
     */
    public static class SearchNode implements Comparable<SearchNode> {
        int[] assignments;   // local position -> RC index (-1 = unassigned)
        int level;           // number of assigned lambda positions
        double confLower;    // pairwise lower bound on energy
        double confUpper;    // pairwise upper bound (rigid) on energy
        BigDecimal boltzLower;  // boltz(confUpper) * numConfs
        BigDecimal boltzUpper;  // boltz(confLower) * numConfs
        BigInteger numConfs;

        SearchNode(int numPos) {
            assignments = new int[numPos];
            Arrays.fill(assignments, -1);
            level = 0;
        }

        SearchNode assign(int localPos, int rc) {
            SearchNode child = new SearchNode(assignments.length);
            System.arraycopy(assignments, 0, child.assignments, 0, assignments.length);
            child.assignments[localPos] = rc;
            child.level = this.level + 1;
            return child;
        }

        boolean isLeaf() {
            return level == assignments.length;
        }

        BigDecimal getError() {
            if (boltzUpper == null || boltzLower == null) return BigDecimal.ZERO;
            return boltzUpper.subtract(boltzLower);
        }

        @Override
        public int compareTo(SearchNode other) {
            return other.getError().compareTo(this.getError());
        }
    }

    public SubproblemMARKStar(Set<Integer> lambdaPositions, Set<Integer> mSetPositions,
                               EnergyMatrix globalRigidEmat, EnergyMatrix globalMinEmat,
                               RCs globalRCs) {
        this.lambdaPositions = lambdaPositions;
        this.mSetPositions = mSetPositions;
        this.globalRigidEmat = globalRigidEmat;
        this.globalMinEmat = globalMinEmat;
        this.globalRCs = globalRCs;
        this.numLambdaPos = lambdaPositions.size();
    }

    /** Initialize search tree for a new M-set assignment. */
    public void initializeForMSet(MSetAssignment mSetAssignment) {
        if (boundsByMSet.containsKey(mSetAssignment)) return;

        SubproblemBound bound = new SubproblemBound();

        bound.rigidEmat = new SubproblemEnergyMatrix(globalRigidEmat, lambdaPositions, mSetAssignment, globalRCs);
        bound.minEmat = new SubproblemEnergyMatrix(globalMinEmat, lambdaPositions, mSetAssignment, globalRCs);

        if (numLambdaPos == 0) {
            bound.lowerBound = BigDecimal.ONE;
            bound.upperBound = BigDecimal.ONE;
            boundsByMSet.put(mSetAssignment, bound);
            return;
        }

        SearchNode root = new SearchNode(numLambdaPos);
        scoreNode(root, bound);
        bound.queue.add(root);

        bound.lowerBound = root.boltzLower;
        bound.upperBound = root.boltzUpper;

        boundsByMSet.put(mSetAssignment, bound);
    }

    /**
     * Tighten bounds for a specific M-set assignment by expanding
     * the highest-error node in the subproblem queue.
     * Returns a PendingMinimization if a new leaf was found, null otherwise.
     */
    public PendingMinimization tightenBound(MSetAssignment mSetAssignment) {
        SubproblemBound bound = boundsByMSet.get(mSetAssignment);
        if (bound == null || bound.queue.isEmpty()) return null;

        SearchNode node = bound.queue.poll();
        if (node.isLeaf()) {
            // Fully assigned — record leaf for minimization
            double eRigid = bound.rigidEmat.confE(node.assignments);
            double eMin = bound.minEmat.confE(node.assignments);
            LeafConformation leaf = new LeafConformation(node.assignments.clone(), eRigid, eMin);
            bound.leafConformations.add(leaf);
            bound.pendingLeaves.add(leaf);
            recomputeBounds(bound);

            // Build global conformation for minimization
            int[] globalConf = bound.minEmat.toGlobalConf(node.assignments);
            return new PendingMinimization(mSetAssignment, node.assignments.clone(),
                    globalConf, eRigid, eMin);
        }

        // Choose next position to expand
        int nextLocalPos = chooseNextPosition(node, bound);

        // Expand all RCs at this position
        for (int rci = 0; rci < bound.minEmat.getNumConfAtPos(nextLocalPos); rci++) {
            SearchNode child = node.assign(nextLocalPos, rci);
            scoreNode(child, bound);
            if (child.boltzUpper.compareTo(BigDecimal.ZERO) > 0) {
                bound.queue.add(child);
            }
        }

        bound.numExpansions++;
        recomputeBounds(bound);
        return null;
    }

    /**
     * Update a leaf conformation's energy after continuous minimization.
     * The true energy replaces the rigid/min bounds for this conformation.
     */
    public void updateLeafEnergy(MSetAssignment mSetAssignment, int[] localConf, double trueEnergy) {
        SubproblemBound bound = boundsByMSet.get(mSetAssignment);
        if (bound == null) return;

        for (LeafConformation leaf : bound.leafConformations) {
            if (Arrays.equals(leaf.localConf, localConf)) {
                leaf.trueEnergy = trueEnergy;
                leaf.minimized = true;
                break;
            }
        }

        recomputeBounds(bound);
    }

    /**
     * Get and clear pending leaves for a given M-set assignment.
     */
    public List<LeafConformation> drainPendingLeaves(MSetAssignment mSetAssignment) {
        SubproblemBound bound = boundsByMSet.get(mSetAssignment);
        if (bound == null) return Collections.emptyList();
        List<LeafConformation> result = new ArrayList<>(bound.pendingLeaves);
        bound.pendingLeaves.clear();
        return result;
    }

    /** Score a search node using pairwise energy bounds. */
    private void scoreNode(SearchNode node, SubproblemBound bound) {
        double gscoreMin = 0.0;
        double gscoreRigid = 0.0;
        for (int i = 0; i < numLambdaPos; i++) {
            if (node.assignments[i] != -1) {
                gscoreMin += bound.minEmat.getOneBody(i, node.assignments[i]);
                gscoreRigid += bound.rigidEmat.getOneBody(i, node.assignments[i]);
                for (int j = i + 1; j < numLambdaPos; j++) {
                    if (node.assignments[j] != -1) {
                        gscoreMin += bound.minEmat.getPairwise(i, node.assignments[i], j, node.assignments[j]);
                        gscoreRigid += bound.rigidEmat.getPairwise(i, node.assignments[i], j, node.assignments[j]);
                    }
                }
            }
        }

        double hscoreMin = 0.0;
        double hscoreRigid = 0.0;
        for (int i = 0; i < numLambdaPos; i++) {
            if (node.assignments[i] == -1) {
                double bestMin = Double.MAX_VALUE;
                double bestRigid = Double.MAX_VALUE;
                for (int rci = 0; rci < bound.minEmat.getNumConfAtPos(i); rci++) {
                    double eMin = bound.minEmat.getOneBody(i, rci);
                    double eRigid = bound.rigidEmat.getOneBody(i, rci);
                    for (int j = 0; j < numLambdaPos; j++) {
                        if (node.assignments[j] != -1) {
                            eMin += bound.minEmat.getPairwise(i, rci, j, node.assignments[j]);
                            eRigid += bound.rigidEmat.getPairwise(i, rci, j, node.assignments[j]);
                        }
                    }
                    bestMin = Math.min(bestMin, eMin);
                    bestRigid = Math.min(bestRigid, eRigid);
                }
                if (bestMin < Double.MAX_VALUE) hscoreMin += bestMin;
                if (bestRigid < Double.MAX_VALUE) hscoreRigid += bestRigid;
            }
        }

        node.confLower = gscoreMin + hscoreMin;
        node.confUpper = gscoreRigid + hscoreRigid;

        if (node.confLower > node.confUpper) {
            double temp = node.confLower;
            node.confLower = node.confUpper;
            node.confUpper = temp;
        }

        node.numConfs = BigInteger.ONE;
        for (int i = 0; i < numLambdaPos; i++) {
            if (node.assignments[i] == -1) {
                node.numConfs = node.numConfs.multiply(BigInteger.valueOf(bound.minEmat.getNumConfAtPos(i)));
            }
        }

        node.boltzLower = bc.calc(node.confUpper).multiply(new BigDecimal(node.numConfs));
        node.boltzUpper = bc.calc(node.confLower).multiply(new BigDecimal(node.numConfs));
    }

    private int chooseNextPosition(SearchNode node, SubproblemBound bound) {
        int bestPos = -1;
        int bestCount = Integer.MAX_VALUE;
        for (int i = 0; i < numLambdaPos; i++) {
            if (node.assignments[i] == -1) {
                int count = bound.minEmat.getNumConfAtPos(i);
                if (count < bestCount) {
                    bestCount = count;
                    bestPos = i;
                }
            }
        }
        return bestPos;
    }

    /** Recompute total bounds by summing over all nodes in queue + leaves. */
    private void recomputeBounds(SubproblemBound bound) {
        BigDecimal lower = BigDecimal.ZERO;
        BigDecimal upper = BigDecimal.ZERO;

        for (SearchNode node : bound.queue) {
            lower = lower.add(node.boltzLower);
            upper = upper.add(node.boltzUpper);
        }

        for (LeafConformation leaf : bound.leafConformations) {
            if (leaf.minimized) {
                // After minimization: both bounds use the true energy
                BigDecimal boltz = bc.calc(leaf.trueEnergy);
                lower = lower.add(boltz);
                upper = upper.add(boltz);
            } else {
                // Before minimization: use rigid/min bounds
                double eHigh = Math.max(leaf.rigidEnergy, leaf.minEnergy);
                double eLow = Math.min(leaf.rigidEnergy, leaf.minEnergy);
                lower = lower.add(bc.calc(eHigh));
                upper = upper.add(bc.calc(eLow));
            }
        }

        bound.lowerBound = lower;
        bound.upperBound = upper;
    }

    // ========== Query methods ==========

    public boolean isInitialized(MSetAssignment mSetAssignment) {
        return boundsByMSet.containsKey(mSetAssignment);
    }

    public SubproblemBound getBound(MSetAssignment mSetAssignment) {
        return boundsByMSet.get(mSetAssignment);
    }

    public BigDecimal getLowerBound(MSetAssignment mSetAssignment) {
        SubproblemBound b = boundsByMSet.get(mSetAssignment);
        return b != null ? b.lowerBound : BigDecimal.ZERO;
    }

    public BigDecimal getUpperBound(MSetAssignment mSetAssignment) {
        SubproblemBound b = boundsByMSet.get(mSetAssignment);
        return b != null ? b.upperBound : BigDecimal.ZERO;
    }

    public double getEpsilon(MSetAssignment mSetAssignment) {
        SubproblemBound b = boundsByMSet.get(mSetAssignment);
        if (b == null) return 1.0;
        if (MathTools.isZero(b.upperBound)) return 0.0;
        BigDecimal gap = b.upperBound.subtract(b.lowerBound);
        return MathTools.bigDivide(gap, b.upperBound, PartitionFunction.decimalPrecision).doubleValue();
    }

    public Set<Integer> getLambdaPositions() { return lambdaPositions; }
    public Set<Integer> getMSetPositions() { return mSetPositions; }
    public int getNumLambdaPositions() { return numLambdaPos; }
}
