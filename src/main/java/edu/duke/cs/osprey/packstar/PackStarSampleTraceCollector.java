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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Simple in-memory trace sink for callers that want PACK* samples after a run.
 */
public class PackStarSampleTraceCollector implements PackStarSampleListener {

	private final List<PackStarSampleTrace> samples =
		Collections.synchronizedList(new ArrayList<>());

	@Override
	public void onSample(PackStarSampleTrace sample) {
		if (sample == null) {
			throw new IllegalArgumentException("sample is required");
		}
		samples.add(sample);
	}

	public List<PackStarSampleTrace> snapshot() {
		synchronized (samples) {
			return new ArrayList<>(samples);
		}
	}

	public void clear() {
		samples.clear();
	}
}
