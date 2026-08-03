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

package edu.duke.cs.osprey.packstar;

import java.util.Arrays;

/**
 * Read-only sample record emitted by the PACK* estimator.
 */
public class PackStarSampleTrace {

	public enum Stage {
		TRAIN,
		PILOT,
		ESTIMATION
	}

	public final Stage stage;
	public final int sampleIndex;
	public final double eTrue;
	public final double eMin;
	public final double eProposal;
	public final double xi;
	public final double logWeight;
	public final double logZCorrected;
	public final double logZMinDet;
	public final double clipLogCap;
	public final boolean clipped;

	private final int[] conf;

	public PackStarSampleTrace(Stage stage,
	                           int sampleIndex,
	                           int[] conf,
	                           double eTrue,
	                           double eMin,
	                           double eProposal,
	                           double xi,
	                           double logWeight,
	                           double logZCorrected,
	                           double logZMinDet,
	                           double clipLogCap,
	                           boolean clipped) {
		if (stage == null) {
			throw new IllegalArgumentException("stage is required");
		}
		if (sampleIndex < 0) {
			throw new IllegalArgumentException("sample index must be nonnegative");
		}
		if (conf == null) {
			throw new IllegalArgumentException("conformation is required");
		}
		this.stage = stage;
		this.sampleIndex = sampleIndex;
		this.conf = Arrays.copyOf(conf, conf.length);
		this.eTrue = eTrue;
		this.eMin = eMin;
		this.eProposal = eProposal;
		this.xi = xi;
		this.logWeight = logWeight;
		this.logZCorrected = logZCorrected;
		this.logZMinDet = logZMinDet;
		this.clipLogCap = clipLogCap;
		this.clipped = clipped;
	}

	public int[] getConf() {
		return Arrays.copyOf(conf, conf.length);
	}
}
