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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * DEPRECATED — Superseded by direct emat access in RootedTreeEdge.
 *
 * Was used by SubproblemMARKStar to project global energy matrices to subproblem
 * lambda positions. No longer needed: RootedTreeEdge now computes lambda-only and
 * lambda-M energies directly from the global emat.
 *
 * Kept for reference. Not used by BranchMARKStarBound.
 *
 * @deprecated No longer used.
 */
@Deprecated
public class SubproblemEnergyMatrix {

    private final EnergyMatrix globalEmat;
    private final int[] localToGlobal;   // local position index → global position index
    private final int[] globalToLocal;   // global position index → local position index (-1 if not lambda)
    private final MSetAssignment mSetAssignment;
    private final RCs globalRCs;
    private final int numLocalPos;

    /**
     * @param globalEmat      the global energy matrix
     * @param lambdaPositions the lambda (free) positions for this subproblem (global indices)
     * @param mSetAssignment  the fixed M-set assignment
     * @param globalRCs       the global rotamer conformations
     */
    public SubproblemEnergyMatrix(EnergyMatrix globalEmat, Set<Integer> lambdaPositions,
                                   MSetAssignment mSetAssignment, RCs globalRCs) {
        this.globalEmat = globalEmat;
        this.mSetAssignment = mSetAssignment;
        this.globalRCs = globalRCs;
        this.numLocalPos = lambdaPositions.size();

        // Build index mappings
        List<Integer> sortedLambda = new ArrayList<>(lambdaPositions);
        sortedLambda.sort(Integer::compareTo);
        localToGlobal = new int[numLocalPos];
        for (int i = 0; i < numLocalPos; i++) {
            localToGlobal[i] = sortedLambda.get(i);
        }

        globalToLocal = new int[globalEmat.getNumPos()];
        java.util.Arrays.fill(globalToLocal, -1);
        for (int i = 0; i < numLocalPos; i++) {
            globalToLocal[localToGlobal[i]] = i;
        }
    }

    /** Number of lambda (free) positions in this subproblem. */
    public int getNumPos() {
        return numLocalPos;
    }

    /** Number of RCs at a local position. */
    public int getNumConfAtPos(int localPos) {
        return globalRCs.getNum(localToGlobal[localPos]);
    }

    /** Get the global RC index for (localPos, localRCIndex). */
    public int getGlobalRC(int localPos, int localRCIndex) {
        return globalRCs.get(localToGlobal[localPos], localRCIndex);
    }

    /**
     * Effective one-body energy for a lambda position:
     * E_self(pos, rc) + Σ_{m ∈ M-set} E_pair(pos, rc, m, rc_m)
     */
    public double getOneBody(int localPos, int localRCIndex) {
        int gPos = localToGlobal[localPos];
        int gRC = globalRCs.get(gPos, localRCIndex);

        double energy = globalEmat.getOneBody(gPos, gRC);

        // Add cross-terms with fixed M-set positions
        int[] mPositions = mSetAssignment.getPositions();
        int[] mRCs = mSetAssignment.getRCs();
        for (int i = 0; i < mPositions.length; i++) {
            energy += globalEmat.getPairwise(gPos, gRC, mPositions[i], mRCs[i]);
        }

        return energy;
    }

    /** Pairwise energy between two lambda positions. */
    public double getPairwise(int localPos1, int localRCIndex1, int localPos2, int localRCIndex2) {
        int gPos1 = localToGlobal[localPos1];
        int gRC1 = globalRCs.get(gPos1, localRCIndex1);
        int gPos2 = localToGlobal[localPos2];
        int gRC2 = globalRCs.get(gPos2, localRCIndex2);
        return globalEmat.getPairwise(gPos1, gRC1, gPos2, gRC2);
    }

    /**
     * Total pairwise energy for a full local conformation.
     * localConf[i] = RC index at local position i (as index into RCs.get(globalPos, .)).
     */
    public double confE(int[] localConf) {
        double energy = 0.0;
        for (int i = 0; i < numLocalPos; i++) {
            energy += getOneBody(i, localConf[i]);
            for (int j = i + 1; j < numLocalPos; j++) {
                energy += getPairwise(i, localConf[i], j, localConf[j]);
            }
        }
        return energy;
    }

    /** Convert local position index to global. */
    public int globalPos(int localPos) {
        return localToGlobal[localPos];
    }

    /** Convert global position index to local (-1 if not in subproblem). */
    public int localPos(int globalPos) {
        if (globalPos < 0 || globalPos >= globalToLocal.length) return -1;
        return globalToLocal[globalPos];
    }

    /** Get the M-set assignment this matrix is conditioned on. */
    public MSetAssignment getMSetAssignment() {
        return mSetAssignment;
    }

    /**
     * Convert a local conformation to a global partial assignment.
     * Returns an array of size globalNumPos, with -1 for unassigned positions.
     */
    public int[] toGlobalConf(int[] localConf) {
        int[] globalConf = new int[globalEmat.getNumPos()];
        java.util.Arrays.fill(globalConf, -1);
        // Fill lambda positions
        for (int i = 0; i < numLocalPos; i++) {
            globalConf[localToGlobal[i]] = globalRCs.get(localToGlobal[i], localConf[i]);
        }
        // Fill M-set positions
        int[] mPositions = mSetAssignment.getPositions();
        int[] mRCs = mSetAssignment.getRCs();
        for (int i = 0; i < mPositions.length; i++) {
            globalConf[mPositions[i]] = mRCs[i];
        }
        return globalConf;
    }
}
