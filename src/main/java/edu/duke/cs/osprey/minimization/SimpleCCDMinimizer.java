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

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.Factory;

public class SimpleCCDMinimizer implements Minimizer.NeedsCleanup, Minimizer.Reusable {

	public static final int DefaultMaxIterations = 30; // same as CCDMinimizer
	private int maxIterations = DefaultMaxIterations;
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer

	// DOF value logging configuration
	public static boolean ENABLE_DOF_VALUE_LOGGING = false;

	// ThreadLocal to pass RC tuple information from EnergyCalculator to SimpleCCDMinimizer
	public static ThreadLocal<String> currentRCTuple = new ThreadLocal<>();

	// ThreadLocal to pass unmatched DOF indices (for warm start priority optimization)
	public static ThreadLocal<int[]> unmatchedDOFIndices = new ThreadLocal<>();

	// Lock for thread-safe logging (prevent interleaved output from multiple threads)
	private static final Object LOG_LOCK = new Object();

	private Factory<LineSearcher,Void> lineSearcherFactory;
	private ObjectiveFunction f;
	private List<LineSearcher> lineSearchers;
	private int lastIterations = 0;

	public SimpleCCDMinimizer() {
		this((context) -> new SurfingLineSearcher());
	}
	
	public SimpleCCDMinimizer(ObjectiveFunction f) {
		this();
		init(f);
	}
	
	public SimpleCCDMinimizer(Factory<LineSearcher,Void> lineSearcherFactory) {
		this.lineSearcherFactory = lineSearcherFactory;

		lineSearchers = new ArrayList<>();
	}

	public void setMaxIterations(int val) {
		this.maxIterations = val;
	}

	@Override
	public void init(ObjectiveFunction f) {
		
		this.f = f;
		
		// build the dofs
		lineSearchers.clear();
		
		for (int d=0; d<f.getNumDOFs(); d++) {
			
			ObjectiveFunction.OneDof fd = new ObjectiveFunction.OneDof(f, d);

			if (fd.getXMin() < fd.getXMax()) {
				LineSearcher lineSearcher = lineSearcherFactory.make(null);
				lineSearcher.init(fd);
				lineSearchers.add(lineSearcher);
			} else {
				lineSearchers.add(null);
			}
		}
	}
	
	@Override
	public Minimizer.Result minimizeFromCenter() {
		return minimizeFrom(f.getDOFsCenter());
	}

	@Override
	public Minimizer.Result minimizeFrom(DoubleMatrix1D startx) {

		// Save initial DOF values and RC tuple for logging
		// ALWAYS capture these, regardless of ENABLE_DOF_VALUE_LOGGING flag
		// (the flag may change during execution due to multi-threading)
		String rcTuple = currentRCTuple.get();
		String dofBeforeStr = (rcTuple != null) ? formatDOFs(startx) : null;
		int numDOFs = f.getNumDOFs();

		// Get unmatched DOF indices for warm start priority optimization
		int[] unmatchedDOFs = unmatchedDOFIndices.get();

		// Build priority DOF order: unmatched DOFs first, then the rest
		// This reorders but does NOT skip any DOFs
		int n = f.getNumDOFs();
		int[] dofOrder = buildDOFOrder(n, unmatchedDOFs);

		// Boost step size for unmatched DOFs so they explore more aggressively
		// Unmatched DOFs start from voxel center and need to adapt to cached DOFs
		if (unmatchedDOFs != null && unmatchedDOFs.length > 0) {
			for (int dofIdx : unmatchedDOFs) {
				if (dofIdx >= 0 && dofIdx < n) {
					LineSearcher ls = lineSearchers.get(dofIdx);
					if (ls instanceof SurfingLineSearcher) {
						((SurfingLineSearcher) ls).setStepScale(8.0);
					}
				}
			}
		}

		DoubleMatrix1D herex = startx.copy();
		DoubleMatrix1D nextx = startx.copy();

		// ccd is pretty simple actually
		// just do a line search along each dimension until we stop improving
		// we deal with cycles by just capping the number of iterations

		// get the current objective function value
		double herefx = f.getValue(herex);

		int actualIterations = 0; // Track actual number of iterations
		for (int iter=0; iter<maxIterations; iter++) {
			actualIterations = iter + 1;

			// update all the dofs using line search
			// dofOrder puts unmatched DOFs first so they adapt before matched DOFs
			for (int di=0; di<n; di++) {
				int d = dofOrder[di];
				LineSearcher lineSearcher = lineSearchers.get(d);
				if (lineSearcher != null) {
					double xd = nextx.get(d);
					xd = lineSearcher.search(xd);
					nextx.set(d, xd);
				}
			}

			// how much did we improve?
			double nextfx = f.getValue(nextx);
			double improvement = herefx - nextfx;

			if (improvement > 0) {

				// take the step
				herex.assign(nextx);
				herefx = nextfx;

				if (improvement < ConvergenceThreshold) {
					break;
				}

			} else {
				break;
			}
		}

		// update the protein conf, one last time
		f.setDOFs(herex);
		lastIterations = actualIterations;

		// Clean up ThreadLocal to prevent memory leaks
		if (unmatchedDOFs != null) {
			unmatchedDOFIndices.remove();
		}

		// Log complete minimization record atomically (thread-safe)
		// Only log when ENABLE_DOF_VALUE_LOGGING is explicitly turned on
		if (ENABLE_DOF_VALUE_LOGGING && rcTuple != null && dofBeforeStr != null) {
			synchronized (LOG_LOCK) {
				String rcInfo = " RC=" + rcTuple;
				System.out.println("\n============ SimpleCCDMinimizer: DOF LOGGING ============");
				System.out.println("[MINIMIZATION-START] NumDOFs=" + numDOFs + rcInfo);
				System.out.println("[DOF-BEFORE] DOFs=" + dofBeforeStr + rcInfo);
				System.out.println("[DOF-AFTER] DOFs=" + formatDOFs(herex) + " E=" + String.format("%.4f", herefx) + " Iterations=" + actualIterations + rcInfo);
				System.out.println("==========================================================");
				System.out.flush();
			}
		}

		return new Minimizer.Result(herex, herefx);
	}
	
	/**
	 * Build DOF traversal order: unmatched DOFs first, then remaining DOFs.
	 * All DOFs are included exactly once - only the order changes.
	 */
	private static int[] buildDOFOrder(int n, int[] unmatchedDOFs) {
		if (unmatchedDOFs == null || unmatchedDOFs.length == 0) {
			// No priority info: default sequential order
			int[] order = new int[n];
			for (int i = 0; i < n; i++) order[i] = i;
			return order;
		}

		// Build set of unmatched indices for fast lookup
		boolean[] isUnmatched = new boolean[n];
		int validCount = 0;
		for (int dofIdx : unmatchedDOFs) {
			if (dofIdx >= 0 && dofIdx < n) {
				isUnmatched[dofIdx] = true;
				validCount++;
			}
		}

		int[] order = new int[n];
		int idx = 0;

		// First: unmatched DOFs (these need to adapt to cached values)
		for (int dofIdx : unmatchedDOFs) {
			if (dofIdx >= 0 && dofIdx < n) {
				order[idx++] = dofIdx;
			}
		}

		// Then: all other DOFs in original order
		for (int d = 0; d < n; d++) {
			if (!isUnmatched[d]) {
				order[idx++] = d;
			}
		}

		return order;
	}

	/** Return the number of iterations used in the last minimizeFrom() call. */
	public int getLastIterations() {
		return lastIterations;
	}

	@Override
	public void clean() {
		for (LineSearcher lineSearcher : lineSearchers) {
			if (lineSearcher instanceof LineSearcher.NeedsCleanup) {
				((LineSearcher.NeedsCleanup)lineSearcher).cleanup();
			}
		}
	}

	/**
	 * Format DOF values for logging
	 */
	private static String formatDOFs(DoubleMatrix1D dofs) {
		if (dofs == null) {
			return "null";
		}
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < dofs.size(); i++) {
			if (i > 0) sb.append(", ");
			sb.append(String.format("%.4f", dofs.get(i)));
		}
		sb.append("]");
		return sb.toString();
	}
}
