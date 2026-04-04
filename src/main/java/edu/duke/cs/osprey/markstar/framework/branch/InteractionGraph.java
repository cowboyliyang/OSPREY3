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
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.structure.Residue;

import java.util.*;

/**
 * Represents a residue interaction graph built from OSPREY's energy matrices.
 * Positions are vertices; edges connect positions with significant pairwise energy interactions.
 */
public class InteractionGraph {

    private final int numPositions;
    private final boolean[][] adjacency;
    private final List<int[]> edgeList;

    private InteractionGraph(int numPositions, boolean[][] adjacency) {
        this.numPositions = numPositions;
        this.adjacency = adjacency;
        this.edgeList = new ArrayList<>();
        for (int i = 0; i < numPositions; i++) {
            for (int j = i + 1; j < numPositions; j++) {
                if (adjacency[i][j]) {
                    edgeList.add(new int[]{i, j});
                }
            }
        }
    }

    /**
     * Build an interaction graph from energy matrices.
     * An edge (i,j) exists if the maximum absolute pairwise energy across all RC pairs exceeds the threshold.
     *
     * @param rigidEmat       rigid (upper-bound) energy matrix
     * @param minimizingEmat  minimizing (lower-bound) energy matrix
     * @param rcs             rotamer conformations
     * @param threshold       minimum |pairwise energy| to create an edge (e.g. 0.1 kcal/mol)
     */
    public static InteractionGraph buildFromEnergyMatrix(EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                                          RCs rcs, double threshold) {
        int numPos = rcs.getNumPos();
        boolean[][] adj = new boolean[numPos][numPos];

        for (int i = 0; i < numPos; i++) {
            for (int j = i + 1; j < numPos; j++) {
                double maxAbsEnergy = 0.0;
                for (int rci = 0; rci < rcs.getNum(i); rci++) {
                    int rc1 = rcs.get(i, rci);
                    for (int rcj = 0; rcj < rcs.getNum(j); rcj++) {
                        int rc2 = rcs.get(j, rcj);
                        double rigidE = rigidEmat.getPairwise(i, rc1, j, rc2);
                        double minE = minimizingEmat.getPairwise(i, rc1, j, rc2);
                        maxAbsEnergy = Math.max(maxAbsEnergy, Math.abs(rigidE));
                        maxAbsEnergy = Math.max(maxAbsEnergy, Math.abs(minE));
                        maxAbsEnergy = Math.max(maxAbsEnergy, Math.abs(rigidE - minE));
                    }
                }
                if (maxAbsEnergy > threshold) {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }

        return new InteractionGraph(numPos, adj);
    }

    /**
     * Build an interaction graph using BWM-style dual cutoff: distance AND energy.
     * An edge (i,j) exists only if minDist <= distCutoff AND maxAbsEnergy > energyCutoff.
     *
     * Distance is computed from template coordinates (closest atom pair between residues).
     * Energy uses the same max-absolute logic as {@link #buildFromEnergyMatrix}.
     *
     * @param confSpace       the conformation space (used for template residue coordinates)
     * @param rigidEmat       rigid (upper-bound) energy matrix
     * @param minimizingEmat  minimizing (lower-bound) energy matrix
     * @param rcs             rotamer conformations
     * @param distCutoff      max distance (Angstroms) to keep an edge
     * @param energyCutoff    min |pairwise energy| (kcal/mol) to keep an edge
     */
    public static InteractionGraph buildWithDualCutoff(
            SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
            RCs rcs, double distCutoff, double energyCutoff) {

        int numPos = rcs.getNumPos();
        boolean[][] adj = new boolean[numPos][numPos];
        int edgesKept = 0;
        int edgesCut = 0;

        for (int i = 0; i < numPos; i++) {
            for (int j = i + 1; j < numPos; j++) {
                // Distance: use template coordinates (min distance over all atom pairs)
                SimpleConfSpace.Position posI = confSpace.positions.get(i);
                SimpleConfSpace.Position posJ = confSpace.positions.get(j);
                Residue resI = posI.strand.mol.getResByPDBResNumber(posI.resNum);
                Residue resJ = posJ.strand.mol.getResByPDBResNumber(posJ.resNum);
                double minDist = resI.distanceTo(resJ);

                // Energy: max absolute pairwise energy across all RC pairs
                double maxAbsEnergy = 0.0;
                for (int rci = 0; rci < rcs.getNum(i); rci++) {
                    int rc1 = rcs.get(i, rci);
                    for (int rcj = 0; rcj < rcs.getNum(j); rcj++) {
                        int rc2 = rcs.get(j, rcj);
                        double rigidE = rigidEmat.getPairwise(i, rc1, j, rc2);
                        double minE = minimizingEmat.getPairwise(i, rc1, j, rc2);
                        maxAbsEnergy = Math.max(maxAbsEnergy, Math.abs(rigidE));
                        maxAbsEnergy = Math.max(maxAbsEnergy, Math.abs(minE));
                    }
                }

                // Dual cutoff: both conditions must be met to keep the edge
                if (minDist <= distCutoff && maxAbsEnergy > energyCutoff) {
                    adj[i][j] = true;
                    adj[j][i] = true;
                    edgesKept++;
                } else {
                    System.out.println("Cutting (" + i + "," + j + "): distance " +
                            String.format("%.2f", minDist) + ", energy " +
                            String.format("%.4f", maxAbsEnergy));
                    edgesCut++;
                }
            }
        }

        int totalPossible = numPos * (numPos - 1) / 2;
        System.out.println("InteractionGraph dual cutoff: kept " + edgesKept + "/" + totalPossible
                + " edges (cut " + edgesCut + "), distCutoff=" + distCutoff + ", energyCutoff=" + energyCutoff);

        return new InteractionGraph(numPos, adj);
    }

    public boolean hasEdge(int pos1, int pos2) {
        return adjacency[pos1][pos2];
    }

    public Set<Integer> getNeighbors(int pos) {
        Set<Integer> neighbors = new LinkedHashSet<>();
        for (int j = 0; j < numPositions; j++) {
            if (adjacency[pos][j])
                neighbors.add(j);
        }
        return neighbors;
    }

    public int getNumPositions() {
        return numPositions;
    }

    public int getNumEdges() {
        return edgeList.size();
    }

    public List<int[]> getEdgeList() {
        return edgeList;
    }

    /**
     * Build a complete interaction graph where all position pairs are edges.
     * Matches the BWM* approach of using the full interaction graph without filtering.
     */
    public static InteractionGraph buildComplete(int numPositions) {
        boolean[][] adj = new boolean[numPositions][numPositions];
        for (int i = 0; i < numPositions; i++) {
            for (int j = i + 1; j < numPositions; j++) {
                adj[i][j] = true;
                adj[j][i] = true;
            }
        }
        return new InteractionGraph(numPositions, adj);
    }

    /** Returns graph density as fraction of possible edges present. */
    public double getDensity() {
        int maxEdges = numPositions * (numPositions - 1) / 2;
        return maxEdges > 0 ? (double) edgeList.size() / maxEdges : 1.0;
    }
}
