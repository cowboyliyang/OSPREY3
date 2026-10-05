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

/**
 * Receives one trace per draw when registered, including duplicates and cache
 * hits. Calls are serial, in draw order, on the thread running compute(), after
 * each CCD batch and before its certificate checks. Listener exceptions
 * propagate to the caller. Paths that do not sample emit no traces.
 */
public interface PackStarSampleListener {

	void onSample(PackStarSampleTrace sample);
}
