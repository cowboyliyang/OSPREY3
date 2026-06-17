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
**
** Contact Info:
**    Bruce Donald, Duke University, Department of Computer Science
**    e-mail: www.cs.duke.edu/brd/
*/

package edu.duke.cs.osprey.wmb;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * A log-domain graphical-model view of an energy matrix over the unassigned
 * positions of a partial conformation.
 *
 * <p>The partition-function contribution of the subtree rooted at a partial
 * assignment is {@code Z = sum_completion exp(-E(completion)/RT)}.  Writing the
 * pairwise-decomposable energy out, the assigned positions contribute a fixed
 * constant, each assigned-to-unassigned interaction folds into the one-body
 * term of the unassigned position, and the unassigned-to-unassigned
 * interactions remain pairwise.  This class exposes those pieces in the log
 * domain (theta = -E/RT) so a mini-bucket pass can bound {@code log Z} without
 * ever enumerating a completion.</p>
 *
 * <p>The unassigned positions become variables {@code 0..numVars-1}; variable
 * {@code v} ranges over the rotamers available at {@link #position(int)}, and
 * its domain value {@code k} maps to the energy-matrix rotamer index
 * {@link #rotamer(int, int)}.</p>
 */
public class WmbModel {

	private final int numVars;
	private final int[] posOfVar;        // variable -> energy-matrix position
	private final int[][] rcOfVar;       // variable, domain value -> energy-matrix rotamer
	private final int[] domains;         // variable -> domain size
	private final double[][] logUnary;     // variable, domain value -> -E/RT (folded)
	private final double[][][][] logPair;  // a < b -> table[ka][kb] = -E/RT, else null
	private final double logConstant;    // -E(assigned part)/RT
	private final double rt;

	public WmbModel(EnergyMatrix emat, RCs rcs, int[] assignments, double rt) {
		this.rt = rt;

		int numPos = rcs.getNumPos();
		List<Integer> unassigned = new ArrayList<>();
		List<Integer> assigned = new ArrayList<>();
		for (int pos = 0; pos < numPos; pos++) {
			if (assignments[pos] < 0) {
				unassigned.add(pos);
			} else {
				assigned.add(pos);
			}
		}

		this.numVars = unassigned.size();
		this.posOfVar = new int[numVars];
		this.rcOfVar = new int[numVars][];
		this.domains = new int[numVars];
		for (int v = 0; v < numVars; v++) {
			int pos = unassigned.get(v);
			posOfVar[v] = pos;
			int n = rcs.getNum(pos);
			int[] rcMap = new int[n];
			for (int k = 0; k < n; k++) {
				rcMap[k] = rcs.get(pos, k);
			}
			rcOfVar[v] = rcMap;
			domains[v] = n;
		}

		// constant: full energy of the assigned sub-assignment
		// the energy-matrix constant offset is part of every conformation's energy,
		// just like the A* g-score adds emat.getConstTerm(), so it belongs here
		double assignedEnergy = emat.getConstTerm();
		for (int a = 0; a < assigned.size(); a++) {
			int posA = assigned.get(a);
			int rcA = assignments[posA];
			assignedEnergy += emat.getEnergy(posA, rcA);
			for (int b = a + 1; b < assigned.size(); b++) {
				int posB = assigned.get(b);
				int rcB = assignments[posB];
				assignedEnergy += emat.getEnergy(posA, rcA, posB, rcB);
			}
		}
		this.logConstant = -assignedEnergy / rt;

		// one-body, folding in every assigned-to-unassigned interaction
		this.logUnary = new double[numVars][];
		for (int v = 0; v < numVars; v++) {
			int pos = posOfVar[v];
			double[] table = new double[domains[v]];
			for (int k = 0; k < domains[v]; k++) {
				int rc = rcOfVar[v][k];
				double e = emat.getEnergy(pos, rc);
				for (int posS : assigned) {
					e += emat.getEnergy(pos, rc, posS, assignments[posS]);
				}
				table[k] = -e / rt;
			}
			logUnary[v] = table;
		}

		// pairwise among unassigned positions
		this.logPair = new double[numVars][numVars][][];
		for (int a = 0; a < numVars; a++) {
			for (int b = a + 1; b < numVars; b++) {
				int posA = posOfVar[a];
				int posB = posOfVar[b];
				double[][] table = new double[domains[a]][domains[b]];
				for (int ka = 0; ka < domains[a]; ka++) {
					int rcA = rcOfVar[a][ka];
					for (int kb = 0; kb < domains[b]; kb++) {
						int rcB = rcOfVar[b][kb];
						table[ka][kb] = -emat.getEnergy(posA, rcA, posB, rcB) / rt;
					}
				}
				logPair[a][b] = table;
			}
		}
	}

	public int numVars() {
		return numVars;
	}

	public int[] domains() {
		return domains;
	}

	public double logConstant() {
		return logConstant;
	}

	public double rt() {
		return rt;
	}

	public int position(int var) {
		return posOfVar[var];
	}

	public int rotamer(int var, int domainValue) {
		return rcOfVar[var][domainValue];
	}

	public double[] logUnary(int var) {
		return logUnary[var];
	}

	/** Pairwise log-potential for {@code a < b}, or null when out of order. */
	public double[][] logPair(int a, int b) {
		return a < b ? logPair[a][b] : null;
	}

	/** Energy-matrix graphs are dense, so every pair of variables interacts. */
	public List<int[]> edges() {
		List<int[]> edges = new ArrayList<>();
		for (int a = 0; a < numVars; a++) {
			for (int b = a + 1; b < numVars; b++) {
				edges.add(new int[]{a, b});
			}
		}
		return edges;
	}
}
