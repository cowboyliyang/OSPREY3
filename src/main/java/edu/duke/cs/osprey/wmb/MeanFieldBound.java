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

/**
 * Mean-field (Gibbs-Bogoliubov) lower bound on a subtree partition function.
 *
 * <p>For any distribution {@code q} over conformations, Jensen's inequality gives
 * {@code log Z = log E_q[e^theta/q] >= E_q[theta] + H(q)} -- an "average
 * log-potential plus entropy" lower bound that holds for <i>every</i> {@code q},
 * so it is valid regardless of how the optimization turns out.  Restricting
 * {@code q} to a fully factorized form {@code q(x) = prod_i q_i(x_i)} makes the
 * bound a sum of one- and two-variable expectations plus per-position entropies,
 * computable in {@code O(edges * q^2)} -- independent of the induced width.</p>
 *
 * <p>Unlike the min-the-rest mini-bucket lower bound, which takes a worst-case
 * configuration per bucket and collapses toward zero on rigid energies with
 * clashes, this bound is an expectation over the mass and stays numerically
 * tame.  Coordinate ascent (the mean-field fixed point) tightens it monotonically
 * toward {@code log Z}; it is exact when the true distribution factorizes.</p>
 *
 * <p>To lower-bound the <i>true</i> partition function this is run on the rigid
 * energy matrix, whose energies over-estimate the minimized ones, so its Z
 * under-estimates the true Z.</p>
 */
public class MeanFieldBound {

	private MeanFieldBound() {}

	public static final int DEFAULT_MAX_SWEEPS = 200;
	public static final double DEFAULT_TOLERANCE = 1e-9;

	/** Natural-log mean-field lower bound on {@code log Z} for the given model. */
	public static double logZLower(WmbModel model) {
		return logZLower(model, DEFAULT_MAX_SWEEPS, DEFAULT_TOLERANCE);
	}

	public static double logZLower(WmbModel model, int maxSweeps, double tolerance) {
		if (maxSweeps < 1) {
			throw new IllegalArgumentException("maxSweeps must be >= 1");
		}
		int numVars = model.numVars();
		if (numVars == 0) {
			return model.logConstant();
		}
		int[] domains = model.domains();

		// factorized beliefs q_i, initialized uniform
		double[][] q = new double[numVars][];
		for (int i = 0; i < numVars; i++) {
			q[i] = new double[domains[i]];
			java.util.Arrays.fill(q[i], 1.0 / domains[i]);
		}

		double prevF = Double.NEGATIVE_INFINITY;
		double freeEnergy = Double.NEGATIVE_INFINITY;
		for (int sweep = 0; sweep < maxSweeps; sweep++) {
			// coordinate ascent: refresh each belief from its neighbors' current beliefs
			for (int i = 0; i < numVars; i++) {
				double[] field = localField(model, q, i);
				q[i] = softmax(field);
			}
			freeEnergy = freeEnergy(model, q);
			if (freeEnergy - prevF <= tolerance && sweep > 0) {
				break;
			}
			prevF = freeEnergy;
		}
		return model.logConstant() + freeEnergy;
	}

	/**
	 * Lower bound on the true subtree partition function (natural log), run on the
	 * rigid energy matrix.  The rigid energies over-estimate the true ones, so this
	 * Z under-estimates the true Z -- a valid lower bound.
	 */
	public static double lowerLogZ(EnergyMatrix rigidEmat, RCs rcs, int[] assignments,
	                               int maxSweeps, double tolerance, double rt) {
		return logZLower(new WmbModel(rigidEmat, rcs, assignments, rt), maxSweeps, tolerance);
	}

	// ---- internals --------------------------------------------------------

	/** field_i(a) = theta_i(a) + sum_{j != i} sum_b q_j(b) theta_ij(a, b). */
	private static double[] localField(WmbModel model, double[][] q, int i) {
		int numVars = model.numVars();
		double[] theta = model.logUnary(i);
		int dimI = theta.length;
		double[] field = theta.clone();
		for (int j = 0; j < numVars; j++) {
			if (j == i) {
				continue;
			}
			boolean iFirst = i < j;
			double[][] pair = iFirst ? model.logPair(i, j) : model.logPair(j, i);
			double[] qj = q[j];
			for (int a = 0; a < dimI; a++) {
				double s = 0.0;
				for (int b = 0; b < qj.length; b++) {
					double qb = qj[b];
					if (qb == 0.0) {
						continue; // 0 * theta = 0, even when theta is -inf
					}
					s += qb * (iFirst ? pair[a][b] : pair[b][a]);
				}
				field[a] += s;
			}
		}
		return field;
	}

	private static double[] softmax(double[] field) {
		double max = Double.NEGATIVE_INFINITY;
		for (double v : field) {
			if (v > max) {
				max = v;
			}
		}
		double[] out = new double[field.length];
		if (max == Double.NEGATIVE_INFINITY) {
			java.util.Arrays.fill(out, 1.0 / field.length);
			return out;
		}
		double sum = 0.0;
		for (int a = 0; a < field.length; a++) {
			double e = Math.exp(field[a] - max);
			out[a] = e;
			sum += e;
		}
		for (int a = 0; a < out.length; a++) {
			out[a] /= sum;
		}
		return out;
	}

	/** F(q) = E_q[theta] + H(q): average log-potential plus factorized entropy. */
	private static double freeEnergy(WmbModel model, double[][] q) {
		int numVars = model.numVars();
		double energy = 0.0;
		double entropy = 0.0;
		for (int i = 0; i < numVars; i++) {
			double[] theta = model.logUnary(i);
			double[] qi = q[i];
			for (int a = 0; a < qi.length; a++) {
				double p = qi[a];
				if (p > 0.0) {
					energy += p * theta[a];
					entropy -= p * Math.log(p);
				}
			}
		}
		for (int i = 0; i < numVars; i++) {
			double[] qi = q[i];
			for (int j = i + 1; j < numVars; j++) {
				double[][] pair = model.logPair(i, j);
				double[] qj = q[j];
				for (int a = 0; a < qi.length; a++) {
					double pa = qi[a];
					if (pa == 0.0) {
						continue;
					}
					for (int b = 0; b < qj.length; b++) {
						double pb = qj[b];
						if (pb > 0.0) {
							energy += pa * pb * pair[a][b];
						}
					}
				}
			}
		}
		return energy + entropy;
	}
}
