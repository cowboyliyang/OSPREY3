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

import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.util.Arrays;

/**
 * DEPRECATED — Superseded by flat M-state indexing in RootedTreeEdge.
 *
 * Was used by SubproblemMARKStar as a key for per-M-state bound tracking.
 * RootedTreeEdge now uses integer mIdx (flat index into logZ arrays) instead.
 *
 * Kept for reference. Not used by BranchMARKStarBound.
 *
 * @deprecated Use RootedTreeEdge.decodeMStatePublic(mIdx) instead.
 */
@Deprecated
public class MSetAssignment implements Comparable<MSetAssignment> {

    private final int[] positions;
    private final int[] rcs;
    private final int hashCode;

    public MSetAssignment(int[] positions, int[] rcs) {
        this.positions = positions.clone();
        this.rcs = rcs.clone();
        this.hashCode = Arrays.hashCode(positions) * 31 + Arrays.hashCode(rcs);
    }

    public int[] getPositions() { return positions; }
    public int[] getRCs() { return rcs; }
    public int size() { return positions.length; }

    /** Get RC for a specific global position, or -1 if not in this assignment. */
    public int getRC(int globalPos) {
        for (int i = 0; i < positions.length; i++) {
            if (positions[i] == globalPos) return rcs[i];
        }
        return -1;
    }

    /**
     * Compute pairwise energy within the M-set positions only.
     * E_M(m) = Σ_i E_self(pos_i, rc_i) + Σ_{i<j} E_pair(pos_i, rc_i, pos_j, rc_j)
     */
    public double computeMSetEnergy(EnergyMatrix emat) {
        double energy = 0.0;
        for (int i = 0; i < positions.length; i++) {
            energy += emat.getOneBody(positions[i], rcs[i]);
            for (int j = i + 1; j < positions.length; j++) {
                energy += emat.getPairwise(positions[i], rcs[i], positions[j], rcs[j]);
            }
        }
        return energy;
    }

    /**
     * Compute pairwise energy between a lambda position/RC and all M-set positions.
     * Σ_{m ∈ M-set} E_pair(lambdaPos, lambdaRC, m, rc_m)
     */
    public double computeCrossEnergy(int lambdaPos, int lambdaRC, EnergyMatrix emat) {
        double energy = 0.0;
        for (int i = 0; i < positions.length; i++) {
            energy += emat.getPairwise(lambdaPos, lambdaRC, positions[i], rcs[i]);
        }
        return energy;
    }

    /** Check if two assignments agree on all shared positions. */
    public boolean isCompatible(MSetAssignment other) {
        int i = 0, j = 0;
        while (i < positions.length && j < other.positions.length) {
            if (positions[i] == other.positions[j]) {
                if (rcs[i] != other.rcs[j]) return false;
                i++; j++;
            } else if (positions[i] < other.positions[j]) {
                i++;
            } else {
                j++;
            }
        }
        return true;
    }

    /** Merge two assignments on disjoint positions. Positions must be disjoint. */
    public MSetAssignment merge(MSetAssignment other) {
        int[] mergedPos = new int[positions.length + other.positions.length];
        int[] mergedRCs = new int[mergedPos.length];
        int i = 0, j = 0, k = 0;
        while (i < positions.length && j < other.positions.length) {
            if (positions[i] < other.positions[j]) {
                mergedPos[k] = positions[i];
                mergedRCs[k] = rcs[i];
                i++;
            } else {
                mergedPos[k] = other.positions[j];
                mergedRCs[k] = other.rcs[j];
                j++;
            }
            k++;
        }
        while (i < positions.length) { mergedPos[k] = positions[i]; mergedRCs[k] = rcs[i]; i++; k++; }
        while (j < other.positions.length) { mergedPos[k] = other.positions[j]; mergedRCs[k] = other.rcs[j]; j++; k++; }
        return new MSetAssignment(mergedPos, mergedRCs);
    }

    @Override
    public int hashCode() { return hashCode; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MSetAssignment)) return false;
        MSetAssignment other = (MSetAssignment) o;
        return Arrays.equals(positions, other.positions) && Arrays.equals(rcs, other.rcs);
    }

    @Override
    public int compareTo(MSetAssignment o) {
        // Compare by positions first, then by RCs
        int cmp = Arrays.compare(positions, o.positions);
        if (cmp != 0) return cmp;
        return Arrays.compare(rcs, o.rcs);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("MSet{");
        for (int i = 0; i < positions.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(positions[i]).append("=").append(rcs[i]);
        }
        sb.append("}");
        return sb.toString();
    }
}
