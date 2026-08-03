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
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
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

	private static long productLong(int[] dims) {
		long p = 1L;
		for (int d : dims) {
			if (d < 0) {
				throw new IllegalArgumentException("negative domain size");
			}
			if (d != 0 && p > Long.MAX_VALUE / d) {
				return Long.MAX_VALUE;
			}
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
		long sizeLong = productLong(dims);
		if (sizeLong > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("WMB factor table too large: " + sizeLong + " cells");
		}
		int size = (int) sizeLong;
		int[] str = strides(dims);

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
		long outSizeLong = productLong(outDims);
		if (outSizeLong > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("WMB factor table too large: " + outSizeLong + " cells");
		}
		int outSize = (int) outSizeLong;
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

	static LogFactor maxEliminate(LogFactor factor, int var) {
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
		long outSizeLong = productLong(outDims);
		if (outSizeLong > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("WMB factor table too large: " + outSizeLong + " cells");
		}
		int outSize = (int) outSizeLong;
		double[] out = new double[outSize];

		int[] coords = new int[outScope.length];
		for (int idx = 0; idx < outSize; idx++) {
			int base = 0;
			for (int a = 0, ia = 0; a < factor.scope.length; a++) {
				if (a == axis) {
					continue;
				}
				base += coords[ia] * inStr[a];
				ia++;
			}
			double max = Double.NEGATIVE_INFINITY;
			for (int k = 0; k < dimVar; k++) {
				max = Math.max(max, factor.table[base + k * axisStr]);
			}
			out[idx] = max;
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

	public static final long UNLIMITED_TABLE_CELLS = Long.MAX_VALUE;

	private static int effectiveIBound(WmbModel model, int requestedIBound, int[] order, long maxTableCells) {
		if (model == null) {
			throw new IllegalArgumentException("model is required");
		}
		if (requestedIBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		if (maxTableCells < 1) {
			throw new IllegalArgumentException("max table cells must be >= 1");
		}
		int width = model.numVars() == 0 ? 0 : inducedWidth(model, order);
		int exactIBound = Math.max(width, 1);
		int iEff = Math.min(requestedIBound, exactIBound);
		if (maxTableCells == UNLIMITED_TABLE_CELLS) {
			return iEff;
		}

		while (iEff > 1 && maxDomainProduct(model, iEff + 1) > maxTableCells) {
			iEff--;
		}
		if (maxDomainProduct(model, iEff + 1) > maxTableCells) {
			throw new IllegalArgumentException("WMB iBound=" + requestedIBound
				+ " cannot satisfy maxTableCells=" + maxTableCells
				+ "; the largest pair table needs " + maxDomainProduct(model, 2) + " cells");
		}
		return iEff;
	}

	private static long maxDomainProduct(WmbModel model, int maxVars) {
		int[] domains = model.domains().clone();
		Arrays.sort(domains);
		long product = 1L;
		int used = 0;
		for (int i = domains.length - 1; i >= 0 && used < maxVars; i--, used++) {
			int d = domains[i];
			if (d != 0 && product > Long.MAX_VALUE / d) {
				return Long.MAX_VALUE;
			}
			product *= d;
		}
		return product;
	}

	private static double run(WmbModel model, int iBound, int[] order, boolean upper, Stats stats) {
		return run(model, iBound, order, upper, stats, null);
	}

	private static double run(WmbModel model, int iBound, int[] order, boolean upper, Stats stats,
	                          List<LogFactor>[] conditionalBuckets) {
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
			if (conditionalBuckets != null) {
				conditionalBuckets[v] = new ArrayList<>(fs);
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

	private static double maxUpper(List<LogFactor> inputFactors, int[] domains,
	                               int[] order, int iBound) {
		Stats stats = new Stats();
		int numVars = domains.length;
		int[] pos = new int[numVars];
		for (int k = 0; k < order.length; k++) {
			pos[order[k]] = k;
		}

		Map<Integer, List<LogFactor>> buckets = new TreeMap<>();
		for (int v = 0; v < numVars; v++) {
			buckets.put(v, new ArrayList<>());
		}
		double constant = 0.0;
		for (LogFactor factor : inputFactors) {
			if (factor.isConstant()) {
				constant += factor.constant();
			} else {
				placeFactor(factor, buckets, pos);
			}
		}

		for (int v : order) {
			List<LogFactor> fs = buckets.get(v);
			if (fs == null || fs.isEmpty()) {
				continue;
			}
			for (List<LogFactor> part : partition(fs, iBound)) {
				LogFactor combined = combine(part, domains);
				stats.maxVars = Math.max(stats.maxVars, combined.scope.length);
				stats.maxCells = Math.max(stats.maxCells, combined.table.length);
				LogFactor msg = maxEliminate(combined, v);
				if (msg.isConstant()) {
					constant += msg.constant();
				} else {
					placeFactor(msg, buckets, pos);
				}
			}
		}
		return constant;
	}

	private static double value(LogFactor factor, int[] assignment) {
		if (factor.scope.length == 0) {
			return factor.constant();
		}
		int idx = 0;
		int stride = 1;
		for (int a = factor.scope.length - 1; a >= 0; a--) {
			int val = assignment[factor.scope[a]];
			if (val < 0) {
				throw new IllegalArgumentException("factor scope contains an unassigned variable");
			}
			idx += val * stride;
			stride *= factor.dims[a];
		}
		return factor.table[idx];
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

	public static final class Sample {
		public final int[] domainValues;
		public final double logQ;

		private Sample(int[] domainValues, double logQ) {
			this.domainValues = domainValues;
			this.logQ = logQ;
		}
	}

	public static final class Proposal {
		private final WmbModel model;
		private final int[] order;
		private final int[] domains;
		private final int iBound;
		private final List<LogFactor>[] conditionalBuckets;

		private Proposal(WmbModel model, int[] order, int[] domains,
		                 int iBound, List<LogFactor>[] conditionalBuckets) {
			this.model = model;
			this.order = order;
			this.domains = domains;
			this.iBound = iBound;
			this.conditionalBuckets = conditionalBuckets;
		}

		public WmbModel model() {
			return model;
		}

		public static final class LogWeightCap {
			public final boolean exact;
			public final String reason;
			public final long assignments;
			public final double logWeightUpper;
			public final double fallbackLogWeightUpper;
			public final double miniBucketLogWeightUpper;
			public final String miniBucketReason;

			private LogWeightCap(boolean exact, String reason, long assignments,
			                     double logWeightUpper, double fallbackLogWeightUpper,
			                     double miniBucketLogWeightUpper, String miniBucketReason) {
				this.exact = exact;
				this.reason = reason;
				this.assignments = assignments;
				this.logWeightUpper = logWeightUpper;
				this.fallbackLogWeightUpper = fallbackLogWeightUpper;
				this.miniBucketLogWeightUpper = miniBucketLogWeightUpper;
				this.miniBucketReason = miniBucketReason;
			}
		}

		private static final class MiniBucketCap {
			final double logWeightUpper;
			final String reason;

			MiniBucketCap(double logWeightUpper, String reason) {
				this.logWeightUpper = logWeightUpper;
				this.reason = reason;
			}
		}

		public Sample sample(Random rng) {
			if (rng == null) {
				throw new IllegalArgumentException("random generator is required");
			}
			int[] assignment = new int[model.numVars()];
			Arrays.fill(assignment, -1);
			double logQ = 0.0;
			for (int oi = order.length - 1; oi >= 0; oi--) {
				int var = order[oi];
				double[] scores = conditionalScores(var, assignment);
				double logNorm = logSumExp(scores, 1.0);
				int chosen;
				if (Double.isFinite(logNorm)) {
					chosen = draw(scores, logNorm, rng);
					logQ += scores[chosen] - logNorm;
				} else {
					chosen = rng.nextInt(domains[var]);
					logQ -= Math.log(domains[var]);
				}
				assignment[var] = chosen;
			}
			return new Sample(assignment, logQ);
		}

		public double logProbability(int[] domainValues) {
			if (domainValues == null || domainValues.length != model.numVars()) {
				throw new IllegalArgumentException("domain assignment has wrong length");
			}
			int[] partial = new int[domainValues.length];
			Arrays.fill(partial, -1);
			double logQ = 0.0;
			for (int oi = order.length - 1; oi >= 0; oi--) {
				int var = order[oi];
				double[] scores = conditionalScores(var, partial);
				double logNorm = logSumExp(scores, 1.0);
				int val = domainValues[var];
				if (val < 0 || val >= domains[var]) {
					throw new IllegalArgumentException("domain value out of range");
				}
				if (Double.isFinite(logNorm)) {
					logQ += scores[val] - logNorm;
				} else {
					logQ -= Math.log(domains[var]);
				}
				partial[var] = val;
			}
			return logQ;
		}

		/**
		 * Conservative lower bound on {@code log q(c)} for every full assignment
		 * sampled by this proposal.
		 */
		public double logProbabilityLowerBound() {
			double logQLower = 0.0;
			for (int oi = order.length - 1; oi >= 0; oi--) {
				int var = order[oi];
				List<LogFactor> factors = conditionalBuckets[var];
				if (factors == null || factors.isEmpty()) {
					logQLower -= Math.log(domains[var]);
					continue;
				}
				double scoreLower = 0.0;
				double scoreUpper = 0.0;
				for (LogFactor factor : factors) {
					double min = minFiniteValue(factor);
					double max = maxFiniteValue(factor);
					if (!Double.isFinite(min) || !Double.isFinite(max)) {
						continue;
					}
					scoreLower += min;
					scoreUpper += max;
				}
				logQLower += scoreLower - scoreUpper - Math.log(domains[var]);
			}
			return logQLower;
		}

		/**
		 * Upper-bound the local importance weight
		 * {@code exp(theta_min(c)) / q(c)} for this model and proposal.
		 *
		 * <p>The fallback bound separates {@code theta_min(c)} and {@code q(c)}:
		 * {@code log Z_upper - min_c log q(c)}.  When the local support is small
		 * enough to enumerate, the exact max over
		 * {@code theta_min(c) - log q(c)} keeps low-probability steric/sentinel
		 * assignments paired with their own low Boltzmann weight.</p>
		 */
		public LogWeightCap logWeightUpperBound(double fallbackLogZUpper, long maxAssignments) {
			double fallback = fallbackLogZUpper - logProbabilityLowerBound();
			MiniBucketCap miniBucket = logWeightMiniBucketUpperBound();
			long assignments = countAssignments();
			if (assignments < 0) {
				return fallbackOrMiniBucket("support size overflow", -1L, miniBucket, fallback);
			}
			if (maxAssignments >= 0 && assignments > maxAssignments) {
				return fallbackOrMiniBucket("support too large for exact local cap",
					assignments, miniBucket, fallback);
			}
			if (assignments == 0L) {
				return fallbackOrMiniBucket("empty proposal support", assignments, miniBucket, fallback);
			}

			int[] values = new int[domains.length];
			double exact = Double.NEGATIVE_INFINITY;
			for (long i = 0; i < assignments; i++) {
				double theta = model.logValue(values);
				double logQ = logProbability(values);
				if (Double.isFinite(theta) && Double.isFinite(logQ)) {
					exact = Math.max(exact, theta - logQ);
				} else if (theta != Double.NEGATIVE_INFINITY && logQ == Double.NEGATIVE_INFINITY) {
					exact = Double.POSITIVE_INFINITY;
					break;
				}
				incrementAssignment(values);
			}
			if (!Double.isFinite(exact)) {
				return fallbackOrMiniBucket("non-finite exact local cap", assignments, miniBucket, fallback);
			}
			return new LogWeightCap(true, "exact local support cap", assignments,
				Math.min(exact, fallback), fallback, miniBucket.logWeightUpper, miniBucket.reason);
		}

		private LogWeightCap fallbackOrMiniBucket(String exactReason, long assignments,
		                                         MiniBucketCap miniBucket, double fallback) {
			if (Double.isFinite(miniBucket.logWeightUpper)) {
				return new LogWeightCap(false, "mini-bucket local cap after " + exactReason,
					assignments, Math.min(miniBucket.logWeightUpper, fallback), fallback,
					miniBucket.logWeightUpper, miniBucket.reason);
			}
			return new LogWeightCap(false, exactReason + "; mini-bucket local cap unavailable: "
				+ miniBucket.reason, assignments, fallback, fallback,
				miniBucket.logWeightUpper, miniBucket.reason);
		}

		private MiniBucketCap logWeightMiniBucketUpperBound() {
			try {
				List<LogFactor> factors = new ArrayList<>();
				factors.addAll(initialFactors(model));
				for (int var = 0; var < domains.length; var++) {
					List<LogFactor> bucket = conditionalBuckets[var];
					if (bucket == null || bucket.isEmpty()) {
						factors.add(constantFactor(Math.log(domains[var])));
						continue;
					}
					for (List<LogFactor> part : partition(bucket, iBound)) {
						factors.add(conditionalPenaltyFactor(part, var));
					}
				}
				double cap = maxUpper(factors, domains, order, iBound);
				return Double.isFinite(cap)
					? new MiniBucketCap(cap, "mini-bucket local cap")
					: new MiniBucketCap(cap, "non-finite mini-bucket objective");
			} catch (RuntimeException ex) {
				return new MiniBucketCap(Double.POSITIVE_INFINITY,
					ex.getClass().getSimpleName() + ": " + ex.getMessage());
			}
		}

		private static LogFactor constantFactor(double value) {
			return new LogFactor(new int[0], new int[0], new double[] { value });
		}

		private LogFactor conditionalPenaltyFactor(List<LogFactor> factors, int var) {
			LogFactor combined = combine(factors, domains);
			int axis = indexOf(combined.scope, var);
			if (axis < 0) {
				throw new IllegalArgumentException("conditional bucket is missing eliminated variable");
			}
			int dimVar = combined.dims[axis];
			int[] str = strides(combined.dims);
			int axisStr = str[axis];
			double[] out = new double[combined.table.length];
			int[] contextDims = new int[combined.dims.length - 1];
			for (int i = 0, j = 0; i < combined.dims.length; i++) {
				if (i != axis) {
					contextDims[j++] = combined.dims[i];
				}
			}

			int[] context = new int[contextDims.length];
			long contexts = productLong(contextDims);
			for (long c = 0; c < contexts; c++) {
				int base = 0;
				for (int a = 0, ia = 0; a < combined.dims.length; a++) {
					if (a == axis) {
						continue;
					}
					base += context[ia] * str[a];
					ia++;
				}

				double[] scores = new double[dimVar];
				for (int k = 0; k < dimVar; k++) {
					scores[k] = combined.table[base + k * axisStr];
				}
				double norm = logSumExp(scores, 1.0);
				for (int k = 0; k < dimVar; k++) {
					int idx = base + k * axisStr;
					if (Double.isFinite(norm)) {
						out[idx] = Double.isFinite(scores[k])
							? norm - scores[k]
							: Double.NEGATIVE_INFINITY;
					} else {
						out[idx] = Math.log(dimVar);
					}
				}
				increment(context, contextDims);
			}
			return new LogFactor(combined.scope.clone(), combined.dims.clone(), out);
		}

		private long countAssignments() {
			long n = 1L;
			for (int domain : domains) {
				if (domain <= 0) {
					return 0L;
				}
				if (n > Long.MAX_VALUE / domain) {
					return -1L;
				}
				n *= domain;
			}
			return n;
		}

		private void incrementAssignment(int[] values) {
			for (int v = values.length - 1; v >= 0; v--) {
				values[v]++;
				if (values[v] < domains[v]) {
					return;
				}
				values[v] = 0;
			}
		}

		private double[] conditionalScores(int var, int[] assignment) {
			double[] scores = new double[domains[var]];
			List<LogFactor> factors = conditionalBuckets[var];
			for (int k = 0; k < domains[var]; k++) {
				assignment[var] = k;
				double score = 0.0;
				if (factors != null) {
					for (LogFactor factor : factors) {
						score += value(factor, assignment);
					}
				}
				scores[k] = score;
			}
			assignment[var] = -1;
			return scores;
		}

		private static double minValue(LogFactor factor) {
			double min = Double.POSITIVE_INFINITY;
			for (double val : factor.table) {
				min = Math.min(min, val);
			}
			return min;
		}

		private static double maxValue(LogFactor factor) {
			double max = Double.NEGATIVE_INFINITY;
			for (double val : factor.table) {
				max = Math.max(max, val);
			}
			return max;
		}

		private static double minFiniteValue(LogFactor factor) {
			double min = Double.POSITIVE_INFINITY;
			for (double val : factor.table) {
				if (Double.isFinite(val)) {
					min = Math.min(min, val);
				}
			}
			return min;
		}

		private static double maxFiniteValue(LogFactor factor) {
			double max = Double.NEGATIVE_INFINITY;
			for (double val : factor.table) {
				if (Double.isFinite(val)) {
					max = Math.max(max, val);
				}
			}
			return max;
		}

		private static int draw(double[] scores, double logNorm, Random rng) {
			double u = rng.nextDouble();
			double cdf = 0.0;
			for (int k = 0; k < scores.length; k++) {
				cdf += Math.exp(scores[k] - logNorm);
				if (u <= cdf || k == scores.length - 1) {
					return k;
				}
			}
			return scores.length - 1;
		}
	}

	public static Proposal proposalForModel(WmbModel model, int iBound) {
		return proposalForModel(model, iBound, naturalOrder(model.numVars()));
	}

	public static Proposal proposalForModel(WmbModel model, int iBound, long maxTableCells) {
		return proposalForModel(model, iBound, naturalOrder(model.numVars()), maxTableCells);
	}

	@SuppressWarnings("unchecked")
	public static Proposal proposalForModel(WmbModel model, int iBound, int[] order) {
		return proposalForModel(model, iBound, order, UNLIMITED_TABLE_CELLS);
	}

	@SuppressWarnings("unchecked")
	public static Proposal proposalForModel(WmbModel model, int iBound, int[] order, long maxTableCells) {
		if (model == null) {
			throw new IllegalArgumentException("model is required");
		}
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		int iEff = effectiveIBound(model, iBound, order, maxTableCells);
		List<LogFactor>[] conditionalBuckets = (List<LogFactor>[]) new List[model.numVars()];
		run(model, iEff, order, true, new Stats(), conditionalBuckets);
		return new Proposal(model, order.clone(), model.domains().clone(), iEff, conditionalBuckets);
	}

	/** Deterministic {@code [Z-, Z+]} bracket (natural log) for one model at {@code iBound}. */
	public static MiniBucketBound boundsForModel(WmbModel model, int iBound) {
		return boundsForModel(model, iBound, naturalOrder(model.numVars()));
	}

	public static MiniBucketBound boundsForModel(WmbModel model, int iBound, long maxTableCells) {
		return boundsForModel(model, iBound, naturalOrder(model.numVars()), maxTableCells);
	}

	public static MiniBucketBound boundsForModel(WmbModel model, int iBound, int[] order) {
		return boundsForModel(model, iBound, order, UNLIMITED_TABLE_CELLS);
	}

	public static MiniBucketBound boundsForModel(WmbModel model, int iBound, int[] order, long maxTableCells) {
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		int numVars = model.numVars();
		if (numVars == 0) {
			double c = model.logConstant();
			return new MiniBucketBound(c, c, iBound, 0, 0, 0);
		}
		int width = inducedWidth(model, order);
		int iEff = effectiveIBound(model, iBound, order, maxTableCells);

		Stats su = new Stats();
		double upper = run(model, iEff, order, true, su);
		Stats sl = new Stats();
		double lower = run(model, iEff, order, false, sl);

		return new MiniBucketBound(lower, upper, iEff, width,
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
		return bounds(minimizedEmat, rigidEmat, rcs, assignments, iBound, UNLIMITED_TABLE_CELLS, rt);
	}

	public static MiniBucketBound bounds(edu.duke.cs.osprey.ematrix.EnergyMatrix minimizedEmat,
	                                     edu.duke.cs.osprey.ematrix.EnergyMatrix rigidEmat,
	                                     edu.duke.cs.osprey.astar.conf.RCs rcs,
	                                     int[] assignments, int iBound, long maxTableCells, double rt) {
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
		int iEff = effectiveIBound(minModel, iBound, order, maxTableCells);

		Stats su = new Stats();
		double upper = run(minModel, iEff, order, true, su);
		Stats sl = new Stats();
		double lower = run(rigidModel, iEff, order, false, sl);

		return new MiniBucketBound(lower, upper, iEff, width,
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
		return upperLogZ(minimizedEmat, rcs, assignments, iBound, UNLIMITED_TABLE_CELLS, rt);
	}

	public static double upperLogZ(edu.duke.cs.osprey.ematrix.EnergyMatrix minimizedEmat,
	                               edu.duke.cs.osprey.astar.conf.RCs rcs,
	                               int[] assignments, int iBound, long maxTableCells, double rt) {
		if (iBound < 1) {
			throw new IllegalArgumentException("iBound must be >= 1");
		}
		WmbModel model = new WmbModel(minimizedEmat, rcs, assignments, rt);
		if (model.numVars() == 0) {
			return model.logConstant();
		}
		int[] order = naturalOrder(model.numVars());
		int iEff = effectiveIBound(model, iBound, order, maxTableCells);
		return run(model, iEff, order, true, new Stats());
	}
}
