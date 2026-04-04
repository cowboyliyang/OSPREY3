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

/**
 * Node in the global BranchMARK* scheduling queue.
 * Represents a (lambda-edge, M-state index) pair with its error contribution.
 * The global scheduler expands the node with the largest log-space error first.
 */
public class BranchMARKStarNode implements Comparable<BranchMARKStarNode> {

    private final RootedTreeEdge edge;
    private final int mIdx;
    private double errorBound;  // logZUpper - logZLower for this (edge, mState)

    public BranchMARKStarNode(RootedTreeEdge edge, int mIdx) {
        this.edge = edge;
        this.mIdx = mIdx;
        this.errorBound = edge.getErrorForMState(mIdx);
    }

    public RootedTreeEdge getEdge() { return edge; }
    public int getMIdx() { return mIdx; }
    public double getErrorBound() { return errorBound; }

    public void refreshError() {
        this.errorBound = edge.getErrorForMState(mIdx);
    }

    @Override
    public int compareTo(BranchMARKStarNode other) {
        // Largest error first (descending)
        return Double.compare(other.errorBound, this.errorBound);
    }

    @Override
    public String toString() {
        return String.format("BranchNode{edge=%s, mIdx=%d, err=%.4f}", edge, mIdx, errorBound);
    }
}
