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

import java.math.BigDecimal;
import java.math.MathContext;

import edu.duke.cs.osprey.tools.ExpFunction;

/**
 * A deterministic bracket {@code [Z-, Z+]} on a subtree partition function,
 * reported in the natural-log domain.  The bracket is valid for every
 * {@code iBound}; it collapses onto the exact log Z once {@code iBound} reaches
 * the induced width.
 */
public class MiniBucketBound {

	public final double logZLower;
	public final double logZUpper;
	public final int iBound;
	public final int inducedWidth;
	public final int maxMiniBucketVars;
	public final long maxTableCells;

	public MiniBucketBound(double logZLower, double logZUpper, int iBound,
	                       int inducedWidth, int maxMiniBucketVars, long maxTableCells) {
		this.logZLower = logZLower;
		this.logZUpper = logZUpper;
		this.iBound = iBound;
		this.inducedWidth = inducedWidth;
		this.maxMiniBucketVars = maxMiniBucketVars;
		this.maxTableCells = maxTableCells;
	}

	public double logGap() {
		return logZUpper - logZLower;
	}

	public boolean isExact() {
		return iBound >= inducedWidth;
	}

	public BigDecimal zUpper(ExpFunction ef) {
		return ef.exp(logZUpper);
	}

	public BigDecimal zLower(ExpFunction ef) {
		return ef.exp(logZLower);
	}

	@Override
	public String toString() {
		return String.format("MiniBucketBound[logZ in [%.6f, %.6f], gap=%.6f, iBound=%d, width=%d]",
				logZLower, logZUpper, logGap(), iBound, inducedWidth);
	}

	static BigDecimal expLog(double logValue, MathContext context) {
		return new ExpFunction(context).exp(logValue);
	}
}
