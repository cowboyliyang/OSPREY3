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

package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;

/**
 * Search node for decomposition-guided flat sum mode in BranchMARK*.
 *
 * Each node represents a partial assignment through the branch decomposition tree.
 * "Assigned edges" have their lambda-sets fully assigned; "pending edges" have
 * their M-states determined (by ancestor assignments) but lambda-states not yet chosen.
 *
 * Z bounds are maintained like MARK*'s A* tree:
 * - subtreeUpperBound = numCompletions × boltz(confLowerBound)
 * - subtreeLowerBound = numCompletions × boltz(confUpperBound)
 * where confLowerBound/confUpperBound are energy bounds from g+h scoring.
 */
public class DecompSearchNode implements Comparable<DecompSearchNode> {

    /** An assigned edge: edge + chosen lambda-state index */
    public static class EdgeAssignment {
        public final RootedTreeEdge edge;
        public final int mIdx;
        public final int lambdaIdx;

        public EdgeAssignment(RootedTreeEdge edge, int mIdx, int lambdaIdx) {
            this.edge = edge;
            this.mIdx = mIdx;
            this.lambdaIdx = lambdaIdx;
        }
    }

    /** A pending edge: edge + determined M-state index */
    public static class PendingEdge {
        public final RootedTreeEdge edge;
        public final int mIdx;

        public PendingEdge(RootedTreeEdge edge, int mIdx) {
            this.edge = edge;
            this.mIdx = mIdx;
        }
    }

    // Assigned edges (in order of assignment)
    public final List<EdgeAssignment> assignedEdges;

    // Pending edges (M-states determined, lambda not yet assigned)
    public final List<PendingEdge> pendingEdges;

    // Partial conformation: global RC array (-1 for unassigned positions)
    public final int[] partialConf;

    // Energy bounds (lower = best possible, upper = worst possible)
    // confLowerBound: g_minimizing + h_minimizing (lowest possible total energy)
    // confUpperBound: g_rigid + h_negated_rigid (highest possible total energy)
    public double confLowerBound;
    public double confUpperBound;

    // Number of full conformations consistent with current partial assignment
    public BigInteger numCompletions;

    // Z contribution bounds
    // subtreeUpperBound = numCompletions × boltz(confLowerBound)  [low energy → high boltz → Z upper]
    // subtreeLowerBound = numCompletions × boltz(confUpperBound)  [high energy → low boltz → Z lower]
    public BigDecimal subtreeUpperBound;
    public BigDecimal subtreeLowerBound;

    // Error = subtreeUpperBound - subtreeLowerBound
    public BigDecimal errorBound;

    // For leaves (all edges assigned): minimization state
    public boolean minimized = false;
    public double minimizedEnergy = Double.NaN;

    /**
     * Create a root search node (no edges assigned, all positions unassigned).
     */
    public static DecompSearchNode makeRoot(
            RootedTreeEdge rootEdge,
            int numPos,
            RCs rcs,
            AStarScorer gScorerMin, AStarScorer hScorerMin,
            AStarScorer gScorerRigid, AStarScorer hScorerNegRigid,
            BoltzmannCalculator bc) {

        DecompSearchNode node = new DecompSearchNode(numPos);

        // Root edge is the first pending edge (M-state = 0 for root)
        node.pendingEdges.add(new PendingEdge(rootEdge, 0));

        // Compute initial bounds using scorers on empty assignment
        ConfIndex emptyIndex = buildConfIndex(node.partialConf, numPos);
        node.confLowerBound = gScorerMin.calc(emptyIndex, rcs) + hScorerMin.calc(emptyIndex, rcs);
        // MARK* uses: confUpperBound = rigidScore + (-negatedHScore)
        // negatedHScorer gives -(max of original), so we negate it to get max of original
        node.confUpperBound = gScorerRigid.calc(emptyIndex, rcs) - hScorerNegRigid.calc(emptyIndex, rcs);

        // All conformations are completions
        node.numCompletions = computeNumCompletions(node.partialConf, rcs);

        // Compute Z bounds
        node.computeZBounds(bc);

        return node;
    }

    /**
     * Create a child node by assigning one pending edge's lambda-state.
     */
    public static DecompSearchNode makeChild(
            DecompSearchNode parent,
            int pendingEdgeIndex,
            int lambdaIdx,
            RCs rcs,
            AStarScorer gScorerMin, AStarScorer hScorerMin,
            AStarScorer gScorerRigid, AStarScorer hScorerNegRigid,
            BoltzmannCalculator bc) {

        int numPos = parent.partialConf.length;
        DecompSearchNode child = new DecompSearchNode(numPos);

        // Copy parent's assigned edges
        child.assignedEdges.addAll(parent.assignedEdges);

        // Get the pending edge being assigned
        PendingEdge pe = parent.pendingEdges.get(pendingEdgeIndex);
        RootedTreeEdge edge = pe.edge;
        int mIdx = pe.mIdx;

        // Add this edge assignment
        child.assignedEdges.add(new EdgeAssignment(edge, mIdx, lambdaIdx));

        // Copy parent's partial conf and add new assignments
        System.arraycopy(parent.partialConf, 0, child.partialConf, 0, numPos);

        // Assign M-positions (may already be assigned by ancestors, but set them to be safe)
        int[] mPositions = edge.getMPositionsSorted();
        int[] mGlobalRCs = edge.getMGlobalRCs(mIdx);
        for (int i = 0; i < mPositions.length; i++) {
            child.partialConf[mPositions[i]] = mGlobalRCs[i];
        }

        // Assign lambda-positions
        int[] lambdaPositions = edge.getLambdaPositionsSorted();
        int[] lambdaGlobalRCs = edge.getLambdaGlobalRCs(lambdaIdx);
        for (int i = 0; i < lambdaPositions.length; i++) {
            child.partialConf[lambdaPositions[i]] = lambdaGlobalRCs[i];
        }

        // Copy pending edges (excluding the one we just assigned)
        for (int i = 0; i < parent.pendingEdges.size(); i++) {
            if (i != pendingEdgeIndex) {
                child.pendingEdges.add(parent.pendingEdges.get(i));
            }
        }

        // Add F-set children of the assigned edge as new pending edges
        if (edge.getFset() != null && !edge.getFset().isEmpty()) {
            // Need local RCs to compute child M-states
            int[] mRCsLocal = edge.decodeMStatePublic(mIdx);
            int[] lambdaRCsLocal = edge.decodeLambdaStatePublic(lambdaIdx);

            for (RootedTreeEdge fEdge : edge.getFset()) {
                int[] fMRCs = edge.getMstateForFullState(mRCsLocal, lambdaRCsLocal, fEdge);
                int fMIdx = fEdge.computeIndexInA(fMRCs);
                child.pendingEdges.add(new PendingEdge(fEdge, fMIdx));
            }
        }

        // Compute energy bounds using scorers
        ConfIndex childIndex = buildConfIndex(child.partialConf, numPos);
        child.confLowerBound = gScorerMin.calc(childIndex, rcs) + hScorerMin.calc(childIndex, rcs);
        // MARK* uses: confUpperBound = rigidScore + (-negatedHScore)
        child.confUpperBound = gScorerRigid.calc(childIndex, rcs) - hScorerNegRigid.calc(childIndex, rcs);

        // Compute number of completions
        child.numCompletions = computeNumCompletions(child.partialConf, rcs);

        // Compute Z bounds
        child.computeZBounds(bc);

        return child;
    }

    /** Is this a leaf node (all edges assigned = full conformation)? */
    public boolean isLeaf() {
        return pendingEdges.isEmpty();
    }

    /**
     * Pick the pending edge with the largest error contribution.
     * Returns index into pendingEdges list.
     */
    public int pickHighestErrorPendingEdge() {
        if (pendingEdges.size() == 1) return 0;

        // Heuristic: pick the edge with the most lambda-states (largest branching factor)
        // This is a simple heuristic; could be improved with per-edge error estimates.
        int bestIdx = 0;
        int bestSize = pendingEdges.get(0).edge.getTotalLambdaStates();
        for (int i = 1; i < pendingEdges.size(); i++) {
            int size = pendingEdges.get(i).edge.getTotalLambdaStates();
            if (size > bestSize) {
                bestSize = size;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    /** Priority queue ordering: descending by errorBound (largest error first) */
    @Override
    public int compareTo(DecompSearchNode other) {
        return other.errorBound.compareTo(this.errorBound);
    }

    // ========== Private implementation ==========

    private DecompSearchNode(int numPos) {
        this.assignedEdges = new ArrayList<>();
        this.pendingEdges = new ArrayList<>();
        this.partialConf = new int[numPos];
        Arrays.fill(this.partialConf, -1);
        this.numCompletions = BigInteger.ONE;
        this.subtreeUpperBound = BigDecimal.ZERO;
        this.subtreeLowerBound = BigDecimal.ZERO;
        this.errorBound = BigDecimal.ZERO;
    }

    /**
     * Compute Z bounds from energy bounds and numCompletions.
     *
     * confLowerBound → subtreeUpperBound (low energy → high Boltzmann → Z upper)
     * confUpperBound → subtreeLowerBound (high energy → low Boltzmann → Z lower)
     *
     * Special care: energy matrices can contain +Infinity (hard clashes), and the
     * negated rigid h-scorer can produce -Infinity. BoltzmannCalculator.calc() uses
     * ExpFunction.exp() which crashes on such values (BigDecimal.pow overflow).
     * We guard against this here.
     */
    private void computeZBounds(BoltzmannCalculator bc) {
        BigDecimal numComp = new BigDecimal(numCompletions);

        // Z upper bound from confLowerBound
        // confLowerBound = g_min + h_min; should always be finite (MPLP gives finite results)
        // But guard just in case: if not finite, use 0 (node will have 0 error, gets ignored)
        subtreeUpperBound = safeBoltzCalc(bc, confLowerBound).multiply(numComp);

        // Z lower bound from confUpperBound
        // confUpperBound = g_rigid + h_neg_rigid; can be -Infinity when rigid emat has +Inf clashes
        // If -Infinity: boltz(-Inf) = exp(+Inf) = overflow. Return 0 (conservative: no Z lower from this node)
        subtreeLowerBound = safeBoltzCalc(bc, confUpperBound).multiply(numComp);

        // Ensure lower ≤ upper (can happen due to negated h-score math)
        if (subtreeLowerBound.compareTo(subtreeUpperBound) > 0) {
            BigDecimal temp = subtreeLowerBound;
            subtreeLowerBound = subtreeUpperBound;
            subtreeUpperBound = temp;
        }

        // Error = gap between upper and lower
        errorBound = subtreeUpperBound.subtract(subtreeLowerBound);
        if (errorBound.signum() < 0) errorBound = BigDecimal.ZERO;
    }

    /**
     * Safe Boltzmann weight calculation that handles extreme/infinite energies.
     *
     * bc.calc() uses ExpFunction.exp() which crashes when the exponent exceeds
     * BigDecimal.pow's limit of 999999999. This happens when:
     * - energy is -Infinity → exponent = +Infinity → (int)Infinity = MAX_VALUE > 999999999
     * - energy is extremely negative → exponent > 999999999
     *
     * Returns BigDecimal.ZERO for any non-computable case (conservative bound).
     */
    private static BigDecimal safeBoltzCalc(BoltzmannCalculator bc, double energy) {
        if (Double.isNaN(energy) || Double.isInfinite(energy)) {
            // +Infinity energy → Boltzmann weight = 0
            // -Infinity energy → Boltzmann weight would be +Infinity, not representable → 0
            // NaN → no useful bound → 0
            return BigDecimal.ZERO;
        }
        // Guard extreme finite values: -energy/RT must fit in int and be < 999999999
        double RT = BoltzmannCalculator.RClassic * 298.0;
        double exponent = -energy / RT;
        if (exponent > 900000000) {
            return BigDecimal.ZERO;
        }
        return bc.calc(energy);
    }

    /**
     * Build a ConfIndex from a partial conformation array.
     * partialConf[pos] = rc if assigned, -1 if unassigned.
     */
    static ConfIndex buildConfIndex(int[] partialConf, int numPos) {
        ConfIndex index = new ConfIndex(numPos);
        // Collect defined and undefined positions
        // ConfIndex expects definedPos to be sorted (ascending)
        for (int pos = 0; pos < numPos; pos++) {
            if (partialConf[pos] >= 0) {
                index.definedPos[index.numDefined] = pos;
                index.definedRCs[index.numDefined] = partialConf[pos];
                index.numDefined++;
            } else {
                index.undefinedPos[index.numUndefined] = pos;
                index.numUndefined++;
            }
        }
        return index;
    }

    /**
     * Compute the number of full conformations consistent with a partial assignment.
     * = product of numRCs at each unassigned position.
     */
    private static BigInteger computeNumCompletions(int[] partialConf, RCs rcs) {
        BigInteger count = BigInteger.ONE;
        for (int pos = 0; pos < partialConf.length; pos++) {
            if (partialConf[pos] < 0) {
                int numRCs = rcs.getNum(pos);
                if (numRCs > 0) {
                    count = count.multiply(BigInteger.valueOf(numRCs));
                }
            }
        }
        return count;
    }

    @Override
    public String toString() {
        int assignedCount = 0;
        for (int rc : partialConf) {
            if (rc >= 0) assignedCount++;
        }
        return String.format("DecompSearchNode[assigned=%d/%d, pending=%d edges, leaf=%b, error=%.4e]",
                assignedCount, partialConf.length, pendingEdges.size(), isLeaf(),
                errorBound.doubleValue());
    }
}
