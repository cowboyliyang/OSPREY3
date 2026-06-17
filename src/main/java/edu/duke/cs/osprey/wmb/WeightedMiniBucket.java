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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Weighted mini-bucket elimination: a deterministic, tunable {@code [Z-, Z+]}
 * bound on an energy-matrix subtree partition function.
 *
 * <p>Exact bucket elimination eliminates the variables in a fixed order; the
 * bucket of a variable collects every factor mentioning it, the factors are
 * multiplied, and the variable is summed out.  Its cost is exponential in the
 * induced width.  Mini-bucket elimination caps that cost with a single knob,
 * {@code iBound}: a bucket whose combined scope would exceed {@code iBound + 1}
 * variables is split into mini-buckets, each holding at most that many
 * variables, so every message table has at most {@code q^(iBound+1)} cells
 * regardless of the induced width or the number of completions.</p>
 *
 * <ul>
 *   <li><b>Upper bound:</b> each mini-bucket eliminates the variable with a
 *   weighted (Hoelder) power sum, the weights {@code 1/p} summing to one.
 *   Hoelder's inequality makes the product of the per-mini-bucket power sums an
 *   upper bound on the exact bucket message, so {@code Z+ >= Z}.</li>
 *   <li><b>Lower bound:</b> one mini-bucket eliminates the variable by an exact
 *   sum and the rest by {@code min}.  Replacing a duplicated factor by its
 *   minimum over the eliminated variable can only shrink the message, so
 *   {@code Z- <= Z}.</li>
 * </ul>
 *
 * <p>When {@code iBound >= induced width} no bucket is split and both variants
 * reduce to exact bucket elimination, so the bracket collapses onto the exact
 * {@code log Z}.  No completion is ever enumerated and no energy is ever
 * minimized along this path.</p>
 */
public class WeightedMiniBucket {

	private WeightedMiniBucket() {}

	// ---- log-domain factor ------------------------------------------------

	/** A log-domain factor: a row-major table over a strictly-ascending scope. */
	static final class LogFactor {
		final int[] scope;
		final int[] dims;
		final double[] table;

		LogFactor(int[] scope, int[] dims, double[] table) {
			this.scope = scope;
			this.dims = dims;
			this.table = table;
		}

		boolean isConstant() {
			return scope.length == 0;
		}

		double constant() {
			return table[0];
		}
	}

	private static int[] strides(int[] dims) {
		int[] s = new int[dims.length];
		int acc = 1;
		for (int i = dims.length - 1; i >= 0; i--) {
			s[i] = acc;
			acc *= dims[i];
		}
		return s;
	}

	private static int product(int[] dims) {
		int p = 1;
		for (int d : dims) {
			p *= d;
		}
		return p;
	}

	/** Sorted union of the factor scopes. */
	private static int[] unionScope(List<LogFactor> factors) {
		TreeMap<Integer, Boolean> seen = new TreeMap<>();
		for (LogFactor f : factors) {
			for (int v : f.scope) {
				seen.put(v, Boolean.TRUE);
			}
		}
		int[] out = new int[seen.size()];
		int i = 0;
		for (int v : seen.keySet()) {
			out[i++] = v;
		}
		return out;
	}

	/** Multiply factors == add log tables over the union of their scopes. */
	static LogFactor combine(List<LogFactor> factors, int[] domainOfVar) {
		int[] scope = unionScope(factors);
		if (scope.length == 0) {
			double total = 0.0;
			for (LogFactor f : factors) {
				total += f.constant();
			}
			return new LogFactor(new int[0], new int[0], new double[]{total});
		}
		int[] dims = new int[scope.length];
		for (int i = 0; i < scope.length; i++) {
			dims[i] = domainOfVar[scope[i]];
		}
		int[] str = strides(dims);
		int size = product(dims);

		// For each factor, map its scope onto union axes with the factor's own strides.
		int n = factors.size();
		int[][] axisStride = new int[n][scope.length];
		for (int fi = 0; fi < n; fi++) {
			LogFactor f = factors.get(fi);
			int[] fstr = strides(f.dims);
			for (int j = 0; j < f.scope.length; j++) {
				int var = f.scope[j];
				int uAxis = indexOf(scope, var);
				axisStride[fi][uAxis] = fstr[j];
			}
		}

		double[] table = new double[size];
		int[] coords = new int[scope.length];
		for (int idx = 0; idx < size; idx++) {
			double sum = 0.0;
			for (int fi = 0; fi < n; fi++) {
				int fidx = 0;
				int[] as = axisStride[fi];
				for (int a = 0; a < scope.length; a++) {
					fidx += coords[a] * as[a];
				}
				sum += factors.get(fi).table[fidx];
			}
			table[idx] = sum;
			increment(coords, dims);
		}
		return new LogFactor(scope, dims, table);
	}

	/**
	 * Remove {@code var} from {@code factor}.  {@code weight == 1} is an exact
	 * sum, a weight in {@code (0, 1)} is a Hoelder power sum, and a null weight
	 * is a {@code min} (the lower-bound dual).
	 */
	static LogFactor powerEliminate(LogFactor factor, int var, Double weight) {
		int axis = indexOf(factor.scope, var);
		int dimVar = factor.dims[axis];

		int[] outScope = new int[factor.scope.length - 1];
		int[] outDims = new int[factor.scope.length - 1];
		for (int i = 0, j = 0; i < factor.scope.length; i++) {
			if (i != axis) {
				outScope[j] = factor.scope[i];
				outDims[j] = factor.dims[i];
				j++;
			}
		}

		int[] inStr = strides(factor.dims);
		int axisStr = inStr[axis];
		int outSize = product(outDims);
		double[] out = new double[outSize];

		int[] coords = new int[outScope.length];
		double[] scratch = new double[dimVar];
		for (int idx = 0; idx < outSize; idx++) {
			// base flat index in the input table with var = 0
			int base = 0;
			for (int a = 0, ia = 0; a < factor.scope.length; a++) {
				if (a == axis) {
					continue;
				}
				base += coords[ia] * inStr[a];
				ia++;
			}
			for (int k = 0; k < dimVar; k++) {
				scratch[k] = factor.table[base + k * axisStr];
			}
			out[idx] = reduce(scratch, weight);
			increment(coords, outDims);
		}
		return new LogFactor(outScope, outDims, out);
	}

	private static double reduce(double[] values, Double weight) {
		if (weight == null) {
			double m = Double.POSITIVE_INFINITY;
			for (double v : values) {
				if (v < m) {
					m = v;
				}
			}
			return m;
		}
		if (weight == 1.0) {
			return logSumExp(values, 1.0);
		}
		return weight * logSumExp(values, weight);
	}

	/** log( sum_k exp(values[k] / weight) ) * (handled stably). */
	private static double logSumExp(double[] values, double weight) {
		double max = Double.NEGATIVE_INFINITY;
		for (double v : values) {
			double scaled = v / weight;
			if (scaled > max) {
				max = scaled;
			}
		}
		if (max == Double.NEGATIVE_INFINITY) {
			return Double.NEGATIVE_INFINITY;
		}
		double sum = 0.0;
		for (double v : values) {
			sum += Math.exp(v / weight - max);
		}
		return max + Math.log(sum);
	}

	private static void increment(int[] coords, int[] dims) {
		for (int a = coords.length - 1; a >= 0; a--) {
			coords[a]++;
			if (coords[a] < dims[a]) {
				return;
			}
			coords[a] = 0;
		}
	}

	private static int indexOf(int[] arr, int value) {
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] == value) {
				return i;
			}
		}
		return -1;
	}

	// ---- mini-bucket partitioning ----------------------------------------

	/** Greedy scope-based partition: each mini-bucket scope has <= iBound+1 vars. */
	static List<List<LogFactor>> partition(List<LogFactor> factors, int iBound) {
		List<LogFactor> ordered = new ArrayList<>(factors);
		ordered.sort(Comparator.comparingInt((LogFactor f) -> -f.scope.length));
		List<int[]> scopes = new ArrayList<>();   // mini-bucket scopes (sorted)
		List<List<LogFactor>> buckets = new ArrayList<>();
		for (LogFactor f : ordered) {
			boolean placed = false;
			for (int i = 0; i < buckets.size(); i++) {
				int[] merged = mergeSorted(scopes.get(i), f.scope);
				if (merged.length <= iBound + 1) {
					scopes.set(i, merged);
					buckets.get(i).add(f);
					placed = true;
					break;
				}
			}
			if (!placed) {
				List<LogFactor> b = new ArrayList<>();
				b.add(f);
				buckets.add(b);
				scopes.add(f.scope.clone());
			}
		}
		return buckets;
	}

	private static int[] mergeSorted(int[] a, int[] b) {
		TreeMap<Integer, Boolean> seen = new TreeMap<>();
		for (int v : a) {
			seen.put(v, Boolean.TRUE);
		}
		for (int v : b) {
			seen.put(v, Boolean.TRUE);
		}
		int[] out = new int[seen.size()];
		int i = 0;
		for (int v : seen.keySet()) {
			out[i++] = v;
		}
		return out;
	}

	// ---- elimination driver ----------------------------------------------

	private static final class Stats {
		int maxVars = 0;
		long maxCells = 0;
	}

	private static double run(WmbModel model, int iBound, int[] order, boolean upper, Stats stats) {
		int numVars = model.numVars();
		int[] domainOfVar = model.domains();
		int[] pos = new int[numVars];
		for (int k = 0; k < order.length; k++) {
			pos[order[k]] = k;
		}

		Map<Integer, List<LogFactor>> buckets = new TreeMap<>();
		for (int v = 0; v < numVars; v++) {
			buckets.put(v, new ArrayList<>());
		}
		double logConst = model.logConstant();

		for (LogFactor f : initialFactors(model)) {
			placeFactor(f, buckets, pos);
		}

		for (int v : order) {
			List<LogFactor> fs = buckets.get(v);
			if (fs == null || fs.isEmpty()) {
				continue;
			}
			List<List<LogFactor>> parts = partition(fs, iBound);
			int p = parts.size();
			for (int k = 0; k < p; k++) {
				List<LogFactor> part = parts.get(k);
				Double weight;
				if (p == 1) {
					weight = 1.0;
				} else if (upper) {
					weight = 1.0 / p;
				} else {
					weight = (k == 0) ? Double.valueOf(1.0) : null;
				}
				LogFactor combined = combine(part, domainOfVar);
				stats.maxVars = Math.max(stats.maxVars, combined.scope.length);
				stats.maxCells = Math.max(stats.maxCells, combined.table.length);
				LogFactor msg = powerEliminate(combined, v, weight);
				if (msg.isConstant()) {
					logConst += msg.constant();
				} else {
					placeFactor(msg, buckets, pos);
				}
			}
		}
		return logConst;
	}

	private static void placeFactor(LogFactor f, Map<Integer, List<LogFactor>> buckets, int[] pos) {
		if (f.scope.length == 0) {
			return;
		}
		int earliest = f.scope[0];
		for (int v : f.scope) {
			if (pos[v] < pos[earliest]) {
				earliest = v;
			}
		}
		buckets.get(earliest).add(f);
	}

	private static List<LogFactor> initialFactors(WmbModel model) {
		List<LogFactor> factors = new ArrayList<>();
		int numVars = model.numVars();
		int[] domainOfVar = model.domains();
		for (int v = 0; v < numVars; v++) {
			factors.add(new LogFactor(new int[]{v}, new int[]{domainOfVar[v]}, model.logUnary(v).clone()));
		}
		for (int a = 0; a < numVars; a++) {
			for (int b = a + 1; b < numVars; b++) {
				double[][] pair = model.logPair(a, b);
				if (pair == null) {
					continue;
				}
				int da = domainOfVar[a];
				int db = domainOfVar[b];
				double[] flat = new double[da * db];
				for (int ka = 0; ka < da; ka++) {
					System.arraycopy(pair[ka], 0, flat, ka * db, db);
				}
				factors.add(new LogFactor(new int[]{a, b}, new int[]{da, db}, flat));
			}
		}
		return factors;
	}

	// ---- induced width / elimination order --------------------------------

	static int inducedWidth(WmbModel model, int[] order) {
		int numVars = model.numVars();
		List<java.util.Set<Integer>> adj = new ArrayList<>();
		for (int v = 0; v < numVars; v++) {
			adj.add(new java.util.TreeSet<>());
		}
		for (int[] e : model.edges()) {
			adj.get(e[0]).add(e[1]);
			adj.get(e[1]).add(e[0]);
		}
		boolean[] removed = new boolean[numVars];
		int width = 0;
		for (int v : order) {
			java.util.Set<Integer> nb = new java.util.TreeSet<>();
			for (int u : adj.get(v)) {
				if (!removed[u]) {
					nb.add(u);
				}
			}
			width = Math.max(width, nb.size());
			for (int a : nb) {
				for (int b : nb) {
					if (a != b) {
						adj.get(a).add(b);
					}
				}
			}
			removed[v] = true;
		}
		return width;
	}

	static int[] naturalOrder(int numVars) {
		int[] order = new int[numVars];
		for (int i = 0; i < numVars; i++) {
			order[i] = i;
		}
		return order;
	}

	// ---- public API -------------------------------------------------------

	/** Deterministic {@code [Z-, Z+]} bracket (natural log) for one model at {@code iBound}. */
	public static MiniBucketBound boundsForModel(WmbModel model, int iBound) {
		return boundsForModel(model, iBound, naturalOrder(model.numVars()));
	}

	public static MiniBucketBound boundsForModel(WmbModel model, int iBound, int[] order) {
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		int numVars = model.numVars();
		if (numVars == 0) {
			double c = model.logConstant();
			return new MiniBucketBound(c, c, iBound, 0, 0, 0);
		}
		int width = inducedWidth(model, order);
		int iEff = Math.min(iBound, Math.max(width, 1));

		Stats su = new Stats();
		double upper = run(model, iEff, order, true, su);
		Stats sl = new Stats();
		double lower = run(model, iEff, order, false, sl);

		return new MiniBucketBound(lower, upper, iBound, width,
				Math.max(su.maxVars, sl.maxVars),
				Math.max(su.maxCells, sl.maxCells));
	}

	/**
	 * Deterministic bracket on the true subtree partition function: the upper
	 * bound runs on the minimized energy matrix (whose energies under-estimate
	 * the true ones, so its Z over-estimates), and the lower bound runs on the
	 * rigid energy matrix (whose energies over-estimate, so its Z
	 * under-estimates).
	 */
	public static MiniBucketBound bounds(edu.duke.cs.osprey.ematrix.EnergyMatrix minimizedEmat,
	                                     edu.duke.cs.osprey.ematrix.EnergyMatrix rigidEmat,
	                                     edu.duke.cs.osprey.astar.conf.RCs rcs,
	                                     int[] assignments, int iBound, double rt) {
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		WmbModel minModel = new WmbModel(minimizedEmat, rcs, assignments, rt);
		WmbModel rigidModel = new WmbModel(rigidEmat, rcs, assignments, rt);

		int numVars = minModel.numVars();
		if (numVars == 0) {
			double up = minModel.logConstant();
			double lo = rigidModel.logConstant();
			return new MiniBucketBound(lo, up, iBound, 0, 0, 0);
		}
		int[] order = naturalOrder(numVars);
		int width = inducedWidth(minModel, order);
		int iEff = Math.min(iBound, Math.max(width, 1));

		Stats su = new Stats();
		double upper = run(minModel, iEff, order, true, su);
		Stats sl = new Stats();
		double lower = run(rigidModel, iEff, order, false, sl);

		return new MiniBucketBound(lower, upper, iBound, width,
				Math.max(su.maxVars, sl.maxVars),
				Math.max(su.maxCells, sl.maxCells));
	}

	/**
	 * Just the upper bound on the subtree partition function (natural log), run on
	 * the minimized energy matrix.  This is the well-behaved direction: the
	 * minimized energies under-estimate the true ones, so this Z over-estimates,
	 * giving a valid upper bound.  The rigid-matrix lower bound is intentionally
	 * not computed here -- the rigid extremes collapse it toward zero, so it never
	 * beats the accumulated-exact lower bound the caller already holds.
	 */
	public static double upperLogZ(edu.duke.cs.osprey.ematrix.EnergyMatrix minimizedEmat,
	                               edu.duke.cs.osprey.astar.conf.RCs rcs,
	                               int[] assignments, int iBound, double rt) {
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		WmbModel model = new WmbModel(minimizedEmat, rcs, assignments, rt);
		if (model.numVars() == 0) {
			return model.logConstant();
		}
		int[] order = naturalOrder(model.numVars());
		int iEff = Math.min(iBound, Math.max(inducedWidth(model, order), 1));
		return run(model, iEff, order, true, new Stats());
	}
}
