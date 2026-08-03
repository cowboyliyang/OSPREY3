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
import java.util.Arrays;
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

	/**
	 * One sparse pairwise log-potential for a generic factor graph.
	 *
	 * <p>The table is indexed as {@code logValues[leftValue][rightValue]}.
	 * Variables must be supplied in ascending order so that the same storage
	 * convention is used as {@link #logPair(int, int)}.</p>
	 */
	public static final class PairPotential {
		public final int left;
		public final int right;
		private final double[][] logValues;

		public PairPotential(int left, int right, double[][] logValues) {
			if (left < 0 || right <= left) {
				throw new IllegalArgumentException("pair variables must satisfy 0 <= left < right");
			}
			if (logValues == null) {
				throw new IllegalArgumentException("pair log-potential table is required");
			}
			this.left = left;
			this.right = right;
			this.logValues = copyMatrix(logValues);
		}

		public double[][] logValues() {
			return copyMatrix(logValues);
		}
	}

	private final int numVars;
	private final int[] posOfVar;        // variable -> energy-matrix position
	private final int[][] rcOfVar;       // variable, domain value -> energy-matrix rotamer
	private final int[] domains;         // variable -> domain size
	private final double[][] logUnary;     // variable, domain value -> -E/RT (folded)
	private final double[][][][] logPair;  // a < b -> table[ka][kb] = -E/RT, else null
	private final double logConstant;    // -E(assigned part)/RT
	private final double rt;

	/**
	 * Construct a generic sparse pairwise log-factor model.
	 *
	 * <p>This is the non-OSPREY adapter point used by COHERE-IDP. Variables
	 * map to themselves ({@code position(v) == v}) and domain values map to
	 * themselves ({@code rotamer(v,k) == k}). Potentials are already in the
	 * log domain and use the caller's declared base measure. Missing pair
	 * factors are true graph non-edges, rather than dense zero tables.</p>
	 */
	public static WmbModel fromLogPotentials(
			double logConstant,
			double[][] logUnary,
			List<PairPotential> logPairs
	) {
		return new WmbModel(logConstant, logUnary, logPairs, 1.0);
	}

	private WmbModel(
			double logConstant,
			double[][] sourceLogUnary,
			List<PairPotential> logPairs,
			double rt
		) {
		this(logConstant, sourceLogUnary, logPairs, rt, null, null);
	}

	private WmbModel(
			double logConstant,
			double[][] sourceLogUnary,
			List<PairPotential> logPairs,
			double rt,
			int[] sourcePositions,
			int[][] sourceRotamers
		) {
		if (!Double.isFinite(logConstant)) {
			throw new IllegalArgumentException("logConstant must be finite");
		}
		if (sourceLogUnary == null) {
			throw new IllegalArgumentException("unary log-potentials are required");
		}
		if (logPairs == null) {
			throw new IllegalArgumentException("pair log-potentials are required");
		}
		if (!(rt > 0.0) || Double.isNaN(rt)) {
			throw new IllegalArgumentException("rt must be positive");
		}
		if ((sourcePositions == null) != (sourceRotamers == null)) {
			throw new IllegalArgumentException("position and rotamer maps must be supplied together");
		}
		if (sourcePositions != null
				&& (sourcePositions.length != sourceLogUnary.length
					|| sourceRotamers.length != sourceLogUnary.length)) {
			throw new IllegalArgumentException("position and rotamer maps must match variable count");
		}

		this.numVars = sourceLogUnary.length;
		this.posOfVar = new int[numVars];
		this.rcOfVar = new int[numVars][];
		this.domains = new int[numVars];
		this.logUnary = new double[numVars][];
		for (int var = 0; var < numVars; var++) {
			double[] unary = sourceLogUnary[var];
			if (unary == null || unary.length == 0) {
				throw new IllegalArgumentException("every variable needs a nonempty unary table");
			}
			validateLogValues(unary, "unary log-potential");
			domains[var] = unary.length;
			if (sourcePositions == null) {
				posOfVar[var] = var;
				rcOfVar[var] = new int[unary.length];
				for (int value = 0; value < unary.length; value++) {
					rcOfVar[var][value] = value;
				}
			} else {
				if (sourcePositions[var] < 0
					|| sourceRotamers[var] == null
					|| sourceRotamers[var].length != unary.length) {
					throw new IllegalArgumentException("invalid position or rotamer map");
				}
				posOfVar[var] = sourcePositions[var];
				rcOfVar[var] = sourceRotamers[var].clone();
			}
			this.logUnary[var] = unary.clone();
		}

		this.logPair = new double[numVars][numVars][][];
		for (PairPotential pair : logPairs) {
			if (pair == null) {
				throw new IllegalArgumentException("pair log-potential entries are required");
			}
			if (pair.right >= numVars) {
				throw new IllegalArgumentException("pair variable is out of range");
			}
			if (logPair[pair.left][pair.right] != null) {
				throw new IllegalArgumentException(
					"duplicate pair log-potential for " + pair.left + "," + pair.right
				);
			}
			double[][] table = pair.logValues;
			if (table.length != domains[pair.left]) {
				throw new IllegalArgumentException("pair table has the wrong left domain size");
			}
			double[][] copy = new double[table.length][];
			for (int leftValue = 0; leftValue < table.length; leftValue++) {
				if (table[leftValue] == null || table[leftValue].length != domains[pair.right]) {
					throw new IllegalArgumentException("pair table has the wrong right domain size");
				}
				validateLogValues(table[leftValue], "pair log-potential");
				copy[leftValue] = table[leftValue].clone();
			}
			logPair[pair.left][pair.right] = copy;
		}

		this.logConstant = logConstant;
		this.rt = rt;
	}

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

	/**
	 * Return the same factor graph with every log-potential multiplied by
	 * {@code inverseTemperature}. At zero inverse temperature, all log
	 * potentials are exactly zero and the resulting model is uniform over its
	 * finite domain.
	 *
	 * <p>For an energy-backed model this is equivalent to changing
	 * {@code rt} to {@code rt / inverseTemperature}; the uniform
	 * infinite-temperature endpoint records positive-infinite {@code rt}. It also enables tempered
	 * WMB portfolios for generic sparse models that have no {@link EnergyMatrix}
	 * or {@link RCs} representation.</p>
	 */
	public WmbModel scaledLogPotentials(double inverseTemperature) {
		if (inverseTemperature < 0.0 || !Double.isFinite(inverseTemperature)) {
			throw new IllegalArgumentException("inverseTemperature must be nonnegative and finite");
		}
		double[][] scaledUnary = new double[numVars][];
		for (int var = 0; var < numVars; var++) {
			scaledUnary[var] = scale(logUnary[var], inverseTemperature);
		}
		List<PairPotential> scaledPairs = new ArrayList<>();
		for (int left = 0; left < numVars; left++) {
			for (int right = left + 1; right < numVars; right++) {
				double[][] pair = logPair[left][right];
				if (pair == null) {
					continue;
				}
				double[][] scaledPair = new double[pair.length][];
				for (int leftValue = 0; leftValue < pair.length; leftValue++) {
					scaledPair[leftValue] = scale(pair[leftValue], inverseTemperature);
				}
				scaledPairs.add(new PairPotential(left, right, scaledPair));
			}
		}
		return new WmbModel(
			logConstant * inverseTemperature,
			scaledUnary,
			scaledPairs,
			inverseTemperature == 0.0 ? Double.POSITIVE_INFINITY : rt / inverseTemperature,
			posOfVar,
			rcOfVar
		);
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

	public double logValue(int[] domainValues) {
		if (domainValues == null || domainValues.length != numVars) {
			throw new IllegalArgumentException("domain assignment has wrong length");
		}
		double val = logConstant;
		for (int v = 0; v < numVars; v++) {
			int k = domainValues[v];
			if (k < 0 || k >= domains[v]) {
				throw new IllegalArgumentException("domain value out of range");
			}
			val += logUnary[v][k];
		}
		for (int a = 0; a < numVars; a++) {
			for (int b = a + 1; b < numVars; b++) {
				double[][] pair = logPair[a][b];
				if (pair != null) {
					val += pair[domainValues[a]][domainValues[b]];
				}
			}
		}
		return val;
	}

	/** Return exactly the non-null pairwise interactions in this model. */
	public List<int[]> edges() {
		List<int[]> edges = new ArrayList<>();
		for (int a = 0; a < numVars; a++) {
			for (int b = a + 1; b < numVars; b++) {
				if (logPair[a][b] != null) {
					edges.add(new int[]{a, b});
				}
			}
		}
		return edges;
	}

	private static double[][] copyMatrix(double[][] source) {
		double[][] copy = new double[source.length][];
		for (int row = 0; row < source.length; row++) {
			if (source[row] == null) {
				throw new IllegalArgumentException("log-potential matrix rows are required");
			}
			copy[row] = source[row].clone();
		}
		return copy;
	}

	private static void validateLogValues(double[] values, String label) {
		for (double value : values) {
			if (Double.isNaN(value) || value == Double.POSITIVE_INFINITY) {
				throw new IllegalArgumentException(label + " values must be finite or -infinity");
			}
		}
	}

	private static double[] scale(double[] values, double multiplier) {
		double[] out = Arrays.copyOf(values, values.length);
		for (int i = 0; i < out.length; i++) {
			out[i] *= multiplier;
		}
		return out;
	}
}
