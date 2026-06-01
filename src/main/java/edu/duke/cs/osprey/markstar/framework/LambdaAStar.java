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

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;

import java.util.*;

/**
 * A* search over the lambda-state space of a single branch decomposition edge.
 *
 * Given a fixed M-state, this searches over the lambda positions one at a time,
 * producing full lambda-states in order of increasing minimizing energy (g+h).
 *
 * This avoids enumerating all |RCs|^branchwidth lambda-states upfront.
 * Instead, states are produced lazily on demand via next().
 *
 * Energy scoring:
 *   g = sum of one-body + pairwise energies for assigned lambda positions
 *       (lambda-lambda and lambda-M interactions)
 *   h = sum over unassigned lambda positions of min_{rc} (one-body + pairwise with assigned)
 *       (traditional pairwise h-score)
 *
 * Uses the minimizing energy matrix for the priority (lower energy = better Z contribution).
 */
public class LambdaAStar {

    /** A node in the lambda A* search tree */
    private static class Node implements Comparable<Node> {
        /** Local RC assignments for lambda positions (-1 = unassigned) */
        final int[] lambdaLocalRCs;
        /** Number of lambda positions assigned so far */
        final int numAssigned;
        /** g-score: energy contribution from assigned positions */
        final double gScore;
        /** h-score: optimistic estimate for remaining positions */
        final double hScore;

        double fScore() { return gScore + hScore; }

        Node(int numLambdaPos) {
            this.lambdaLocalRCs = new int[numLambdaPos];
            Arrays.fill(this.lambdaLocalRCs, -1);
            this.numAssigned = 0;
            this.gScore = 0.0;
            this.hScore = 0.0;
        }

        Node(int[] lambdaLocalRCs, int numAssigned, double gScore, double hScore) {
            this.lambdaLocalRCs = lambdaLocalRCs;
            this.numAssigned = numAssigned;
            this.gScore = gScore;
            this.hScore = hScore;
        }

        boolean isLeaf() { return numAssigned >= lambdaLocalRCs.length; }

        @Override
        public int compareTo(Node other) {
            return Double.compare(this.fScore(), other.fScore());
        }
    }

    // Edge info
    private final RootedTreeEdge edge;
    private final int[] mLocalRCs;          // fixed M-state (local RCs)
    private final int[] lambdaPositions;    // = edge.getLambdaPositionsSorted()
    private final int[] mPositions;         // = edge.getMPositionsSorted()
    private final RCs rcs;
    private final EnergyMatrix minEmat;

    // Priority queue
    private final PriorityQueue<Node> openList = new PriorityQueue<>();

    // Number of lambda-states produced so far
    private int numProduced = 0;

    // Total possible lambda-states
    private final int totalStates;

    // Best f-score on the frontier (for aggregate bound)
    private double frontierBestF = Double.POSITIVE_INFINITY;

    /**
     * @param edge       the branch decomposition edge
     * @param mIdx       the fixed M-state index
     * @param minEmat    minimizing energy matrix (for g+h scoring)
     */
    public LambdaAStar(RootedTreeEdge edge, long mIdx, EnergyMatrix minEmat) {
        this.edge = edge;
        this.mLocalRCs = edge.decodeMStatePublic(mIdx);
        this.lambdaPositions = edge.getLambdaPositionsSorted();
        this.mPositions = edge.getMPositionsSorted();
        this.rcs = edge.getRCs();
        this.minEmat = minEmat;
        this.totalStates = edge.getTotalLambdaStates();

        // Create root node (no lambda positions assigned)
        Node root = new Node(lambdaPositions.length);
        double h = computeHScore(root);
        Node rootWithH = new Node(root.lambdaLocalRCs, 0, 0.0, h);
        openList.add(rootWithH);
    }

    /** Get the next lowest-energy full lambda-state, or -1 if exhausted. */
    public int next() {
        while (!openList.isEmpty()) {
            Node node = openList.poll();

            if (node.isLeaf()) {
                // Full lambda-state — encode to flat index and return
                numProduced++;
                frontierBestF = openList.isEmpty() ? Double.POSITIVE_INFINITY : openList.peek().fScore();
                return encodeLambdaState(node.lambdaLocalRCs);
            }

            // Expand: assign the next lambda position
            int posIdx = node.numAssigned;  // index into lambdaPositions
            int globalPos = lambdaPositions[posIdx];
            int numRCsAtPos = rcs.getNum(globalPos);

            for (int localRC = 0; localRC < numRCsAtPos; localRC++) {
                int[] childRCs = node.lambdaLocalRCs.clone();
                childRCs[posIdx] = localRC;
                double childG = node.gScore + computeGIncrement(childRCs, posIdx, localRC);
                int childNumAssigned = node.numAssigned + 1;

                Node child = new Node(childRCs, childNumAssigned, childG, 0.0);
                double childH = computeHScore(child);
                child = new Node(childRCs, childNumAssigned, childG, childH);

                openList.add(child);
            }
        }
        return -1; // exhausted
    }

    /** Number of lambda-states produced so far. */
    public int getNumProduced() { return numProduced; }

    /** Number of remaining (unproduced) lambda-states. */
    public int getNumRemaining() { return totalStates - numProduced; }

    /** Best f-score currently on the frontier (lower bound for all remaining states). */
    public double getFrontierBestF() { return frontierBestF; }

    /** Is the search exhausted? */
    public boolean isExhausted() { return openList.isEmpty() && numProduced >= totalStates; }

    // ========== Energy scoring ==========

    /**
     * Incremental g-score when assigning lambdaPositions[posIdx] = localRC.
     * = one-body(pos, rc) + pairwise with all previously assigned lambda positions
     *   + pairwise with all M positions
     */
    private double computeGIncrement(int[] lambdaLocalRCs, int posIdx, int localRC) {
        int globalPos = lambdaPositions[posIdx];
        int globalRC = rcs.get(globalPos, localRC);
        double energy = minEmat.getOneBody(globalPos, globalRC);

        // Pairwise with assigned lambda positions (indices 0..posIdx-1)
        for (int j = 0; j < posIdx; j++) {
            int otherGlobalPos = lambdaPositions[j];
            int otherGlobalRC = rcs.get(otherGlobalPos, lambdaLocalRCs[j]);
            energy += minEmat.getPairwise(globalPos, globalRC, otherGlobalPos, otherGlobalRC);
        }

        // Pairwise with M positions
        for (int j = 0; j < mPositions.length; j++) {
            int mGlobalPos = mPositions[j];
            int mGlobalRC = rcs.get(mGlobalPos, mLocalRCs[j]);
            energy += minEmat.getPairwise(globalPos, globalRC, mGlobalPos, mGlobalRC);
        }

        return energy;
    }

    /**
     * H-score: for each unassigned lambda position, take the minimum over RCs of
     * (one-body + pairwise with all assigned lambda + pairwise with M).
     * This is the traditional pairwise h-score.
     */
    private double computeHScore(Node node) {
        double h = 0.0;
        for (int i = node.numAssigned; i < lambdaPositions.length; i++) {
            int globalPos = lambdaPositions[i];
            int numRCsAtPos = rcs.getNum(globalPos);
            double minE = Double.POSITIVE_INFINITY;

            for (int localRC = 0; localRC < numRCsAtPos; localRC++) {
                int globalRC = rcs.get(globalPos, localRC);
                double e = minEmat.getOneBody(globalPos, globalRC);

                // Pairwise with assigned lambda positions
                for (int j = 0; j < node.numAssigned; j++) {
                    int otherGlobalPos = lambdaPositions[j];
                    int otherGlobalRC = rcs.get(otherGlobalPos, node.lambdaLocalRCs[j]);
                    e += minEmat.getPairwise(globalPos, globalRC, otherGlobalPos, otherGlobalRC);
                }

                // Pairwise with M positions
                for (int k = 0; k < mPositions.length; k++) {
                    int mGlobalPos = mPositions[k];
                    int mGlobalRC = rcs.get(mGlobalPos, mLocalRCs[k]);
                    e += minEmat.getPairwise(globalPos, globalRC, mGlobalPos, mGlobalRC);
                }

                // Pairwise with other unassigned lambda positions (min over their RCs too)
                // NOTE: standard traditional h-score does NOT include unassigned-unassigned pairwise.
                // Including it would require MPLP. We keep it simple here.

                minE = Math.min(minE, e);
            }
            h += minE;
        }
        return h;
    }

    /**
     * Encode local lambda RCs to flat lambda-state index.
     * Inverse of RootedTreeEdge.decodeLambdaState().
     */
    private int encodeLambdaState(int[] localRCs) {
        int index = 0;
        int stride = 1;
        for (int i = lambdaPositions.length - 1; i >= 0; i--) {
            index += localRCs[i] * stride;
            stride *= rcs.getNum(lambdaPositions[i]);
        }
        return index;
    }
}
