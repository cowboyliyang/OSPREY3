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
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
**
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
**
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;

/**
 * L-BFGS-B minimizer with sparse numerical gradients.
 *
 * Uses per-DOF sparse energy evaluation (getValForDOF) for cheap finite-difference
 * gradients, combined with L-BFGS two-loop recursion and projected backtracking
 * line search for bound-constrained optimization.
 *
 * Designed as a drop-in alternative to SimpleCCDMinimizer for use with GridDP warm start.
 */
public class SimpleLBFGSMinimizer implements Minimizer.NeedsCleanup, Minimizer.Reusable {

	public static final int DefaultMaxIterations = 30;
	private static final double ConvergenceThreshold = 0.001;

	// L-BFGS parameters
	private static final int MEMORY_DEPTH = 7;
	private static final double ARMIJO_C1 = 1e-4;
	private static final double BACKTRACK_RHO = 0.5;
	private static final int MAX_LINESEARCH_STEPS = 20;
	private static final double GAMMA_MIN = 1e-6;
	private static final double GAMMA_MAX = 1e6;

	private int maxIterations = DefaultMaxIterations;
	private ObjectiveFunction f;
	private int lastIterations = 0;

	// Bounds
	private double[] lo;
	private double[] hi;
	private boolean[] pinned; // DOFs where lo == hi

	public SimpleLBFGSMinimizer() {
	}

	public SimpleLBFGSMinimizer(ObjectiveFunction f) {
		init(f);
	}

	public void setMaxIterations(int val) {
		this.maxIterations = val;
	}

	@Override
	public void init(ObjectiveFunction f) {
		this.f = f;
		int n = f.getNumDOFs();
		DoubleMatrix1D[] constraints = f.getConstraints();
		lo = new double[n];
		hi = new double[n];
		pinned = new boolean[n];
		for (int d = 0; d < n; d++) {
			lo[d] = constraints[0].get(d);
			hi[d] = constraints[1].get(d);
			pinned[d] = (lo[d] >= hi[d]);
		}
	}

	@Override
	public Minimizer.Result minimizeFromCenter() {
		return minimizeFrom(f.getDOFsCenter());
	}

	@Override
	public Minimizer.Result minimizeFrom(DoubleMatrix1D startx) {
		int n = f.getNumDOFs();

		// Current point
		double[] x = new double[n];
		for (int d = 0; d < n; d++) {
			x[d] = clamp(d, startx.get(d));
		}

		// Evaluate initial function value and gradient
		DoubleMatrix1D xvec = toVector(x);
		double fx = f.getValue(xvec);
		double[] grad = computeGradient(x);

		// L-BFGS history (circular buffer)
		double[][] sHistory = new double[MEMORY_DEPTH][];
		double[][] yHistory = new double[MEMORY_DEPTH][];
		double[] rhoHistory = new double[MEMORY_DEPTH];
		int historyCount = 0;
		int historyStart = 0;

		int actualIterations = 0;

		for (int iter = 0; iter < maxIterations; iter++) {
			actualIterations = iter + 1;

			// Compute search direction via L-BFGS two-loop recursion
			double[] d = lbfgsTwoLoop(grad, sHistory, yHistory, rhoHistory, historyCount, historyStart, n);

			// Project direction: zero out components that would push beyond bounds
			projectDirection(x, d);

			// Check for descent direction
			double dirDeriv = dot(d, grad);
			if (dirDeriv >= 0) {
				// Not a descent direction, fallback to projected steepest descent
				for (int i = 0; i < n; i++) {
					d[i] = -grad[i];
				}
				projectDirection(x, d);
				dirDeriv = dot(d, grad);
				if (dirDeriv >= 0) {
					// Even steepest descent is not descent (at boundary), stop
					break;
				}
			}

			// Backtracking line search with Armijo condition
			double alpha = 1.0;
			double[] xNew = new double[n];
			double fxNew = Double.POSITIVE_INFINITY;
			boolean lsSuccess = false;

			for (int ls = 0; ls < MAX_LINESEARCH_STEPS; ls++) {
				// x_new = project(x + alpha * d)
				for (int i = 0; i < n; i++) {
					xNew[i] = clamp(i, x[i] + alpha * d[i]);
				}
				DoubleMatrix1D xNewVec = toVector(xNew);
				fxNew = f.getValue(xNewVec);

				if (Double.isNaN(fxNew) || Double.isInfinite(fxNew)) {
					alpha *= BACKTRACK_RHO;
					continue;
				}

				// Armijo condition: f(x_new) <= f(x) + c1 * alpha * d'g
				if (fxNew <= fx + ARMIJO_C1 * alpha * dirDeriv) {
					lsSuccess = true;
					break;
				}
				alpha *= BACKTRACK_RHO;
			}

			if (!lsSuccess) {
				// Line search failed, try accepting any improvement
				if (fxNew < fx) {
					lsSuccess = true;
				} else {
					break;
				}
			}

			// Compute gradient at new point
			double[] gradNew = computeGradient(xNew);

			// L-BFGS update: s = x_new - x, y = g_new - g
			double[] s = new double[n];
			double[] y = new double[n];
			for (int i = 0; i < n; i++) {
				s[i] = xNew[i] - x[i];
				y[i] = gradNew[i] - grad[i];
			}

			double sy = dot(s, y);
			if (sy > 1e-12) {
				// Curvature condition satisfied, update history
				int idx;
				if (historyCount < MEMORY_DEPTH) {
					idx = historyCount;
					historyCount++;
				} else {
					idx = historyStart;
					historyStart = (historyStart + 1) % MEMORY_DEPTH;
				}
				sHistory[idx] = s;
				yHistory[idx] = y;
				rhoHistory[idx] = 1.0 / sy;
			}
			// else: skip update (negative curvature)

			double improvement = fx - fxNew;

			// Accept step
			System.arraycopy(xNew, 0, x, 0, n);
			fx = fxNew;
			System.arraycopy(gradNew, 0, grad, 0, n);

			// Convergence check
			if (improvement >= 0 && improvement < ConvergenceThreshold) {
				break;
			}
		}

		// Apply final DOF values to the molecule
		DoubleMatrix1D result = toVector(x);
		f.setDOFs(result);
		lastIterations = actualIterations;

		return new Minimizer.Result(result, fx);
	}

	/**
	 * L-BFGS two-loop recursion to compute search direction d = -H_k * g.
	 */
	private double[] lbfgsTwoLoop(double[] grad, double[][] sHist, double[][] yHist,
								   double[] rhoHist, int count, int start, int n) {
		double[] q = grad.clone();
		double[] alphas = new double[count];

		// First loop: from most recent to oldest
		for (int j = count - 1; j >= 0; j--) {
			int idx = (start + j) % MEMORY_DEPTH;
			double alpha = rhoHist[idx] * dot(sHist[idx], q);
			alphas[j] = alpha;
			for (int i = 0; i < n; i++) {
				q[i] -= alpha * yHist[idx][i];
			}
		}

		// Initial Hessian scaling: gamma = s'y / y'y (most recent pair)
		double gamma = 1.0;
		if (count > 0) {
			int newest = (start + count - 1) % MEMORY_DEPTH;
			double yy = dot(yHist[newest], yHist[newest]);
			if (yy > 0) {
				gamma = dot(sHist[newest], yHist[newest]) / yy;
				gamma = Math.max(GAMMA_MIN, Math.min(GAMMA_MAX, gamma));
			}
		}

		// r = gamma * q
		double[] r = new double[n];
		for (int i = 0; i < n; i++) {
			r[i] = gamma * q[i];
		}

		// Second loop: from oldest to most recent
		for (int j = 0; j < count; j++) {
			int idx = (start + j) % MEMORY_DEPTH;
			double beta = rhoHist[idx] * dot(yHist[idx], r);
			for (int i = 0; i < n; i++) {
				r[i] += (alphas[j] - beta) * sHist[idx][i];
			}
		}

		// d = -r
		for (int i = 0; i < n; i++) {
			r[i] = -r[i];
		}
		return r;
	}

	/**
	 * Project search direction to respect bounds.
	 * Zero out components that would move a DOF beyond its bound.
	 */
	private void projectDirection(double[] x, double[] d) {
		for (int i = 0; i < x.length; i++) {
			if (pinned[i]) {
				d[i] = 0;
			} else if (d[i] < 0 && x[i] <= lo[i]) {
				d[i] = 0; // at lower bound, can't go lower
			} else if (d[i] > 0 && x[i] >= hi[i]) {
				d[i] = 0; // at upper bound, can't go higher
			}
		}
	}

	/**
	 * Compute sparse numerical gradient using central differences.
	 * Uses getValForDOF for per-DOF sparse energy evaluation.
	 */
	private double[] computeGradient(double[] x) {
		int n = x.length;
		double[] grad = new double[n];

		// First set all DOFs to current values
		f.setDOFs(toVector(x));

		for (int d = 0; d < n; d++) {
			if (pinned[d]) {
				grad[d] = 0;
				continue;
			}

			double range = hi[d] - lo[d];
			double h = Math.max(1e-5, range * 1e-4);

			double xd = x[d];
			double xdPlus = Math.min(hi[d], xd + h);
			double xdMinus = Math.max(lo[d], xd - h);
			double actualH = xdPlus - xdMinus;

			if (actualH < 1e-12) {
				grad[d] = 0;
				continue;
			}

			double fPlus = f.getValForDOF(d, xdPlus);
			double fMinus = f.getValForDOF(d, xdMinus);

			// Restore DOF to original value
			f.setDOF(d, xd);

			if (Double.isNaN(fPlus) || Double.isInfinite(fPlus)
					|| Double.isNaN(fMinus) || Double.isInfinite(fMinus)) {
				grad[d] = 0;
			} else {
				grad[d] = (fPlus - fMinus) / actualH;
			}
		}

		return grad;
	}

	/** Return the number of iterations used in the last minimizeFrom() call. */
	public int getLastIterations() {
		return lastIterations;
	}

	@Override
	public void clean() {
		// No resources to clean
	}

	private double clamp(int d, double val) {
		return Math.max(lo[d], Math.min(hi[d], val));
	}

	private DoubleMatrix1D toVector(double[] x) {
		DoubleMatrix1D v = DoubleFactory1D.dense.make(x.length);
		for (int i = 0; i < x.length; i++) {
			v.set(i, x[i]);
		}
		return v;
	}

	private static double dot(double[] a, double[] b) {
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i] * b[i];
		}
		return sum;
	}
}
