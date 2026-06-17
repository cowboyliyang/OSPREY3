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

package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDpBackend;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;

/**
 * BranchMARK*: public BranchMARK* entry point.
 *
 * <p>The shared branch-DP implementation lives in {@link BranchDpBackend} so
 * other branch-DP backends can reuse it without inheriting BranchMARK* policy.</p>
 */
public class BranchMARKStarBound extends BranchDpBackend {

    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism) {
        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism, null);
    }

    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                               EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                               RCs rcs, Parallelism parallelism, String stateNameOverride) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc,
                rcs, parallelism, stateNameOverride);
    }
}
