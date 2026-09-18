/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
*/

package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.structure.Molecule;

/**
 * A predicate that identifies a functional conformation in a PACK* sample.
 *
 * <p>The predicate is evaluated only on the independent final sample, after
 * CCD minimization and after the proposal has been frozen.  Discrete events
 * can inspect the conformation assignment.  Geometric events can call
 * {@link Sample#getMolecule()} to inspect the minimized coordinates.</p>
 */
@FunctionalInterface
public interface PackStarFunctionalEvent {

    /** Return true when this sample satisfies the functional condition. */
    boolean test(Sample sample);

    /**
     * Data made available to one event evaluation.
     *
     * <p>The molecule returned by {@link #getMolecule()} is owned by the
     * sample evaluator.  Event code should inspect it and should not retain or
     * mutate it after returning from {@link #test(Sample)}.</p>
     */
    final class Sample {

        private final int sampleIndex;
        private final int[] conformation;
        private final double trueEnergy;
        private final double minimizedEnergy;
        private final EnergyCalculator.EnergiedParametricMolecule epmol;
        private Molecule restoredMolecule;

        Sample(int sampleIndex, int[] conformation,
               double trueEnergy, double minimizedEnergy,
               EnergyCalculator.EnergiedParametricMolecule epmol) {
            this.sampleIndex = sampleIndex;
            this.conformation = conformation.clone();
            this.trueEnergy = trueEnergy;
            this.minimizedEnergy = minimizedEnergy;
            this.epmol = epmol;
        }

        /** Index in the logical final sample, including repeated assignments. */
        public int getSampleIndex() {
            return sampleIndex;
        }

        /** Return a copy of the discrete RC assignment. */
        public int[] getConformation() {
            return conformation.clone();
        }

        /** Return one RC assignment by position. */
        public int getRC(int position) {
            return conformation[position];
        }

        /** Return the full minimized energy used by PACK*. */
        public double getTrueEnergy() {
            return trueEnergy;
        }

        /** Return the sparse minimizing-matrix energy used by the proposal. */
        public double getMinimizedEnergy() {
            return minimizedEnergy;
        }

        /**
         * Return the minimized molecule for geometric predicates.
         *
         * <p>For the ordinary Java energy path, the CCD result stores the
         * minimized DOF vector separately from the parametric molecule.  The
         * vector is applied lazily here so existing CCD and sampling paths do
         * not pay a coordinate-copy cost when no geometric event is installed.
         * The returned object is valid for the duration of the predicate.</p>
         */
        public Molecule getMolecule() {
            if (restoredMolecule != null) {
                return restoredMolecule;
            }
            if (epmol == null || epmol.pmol == null || epmol.pmol.mol == null) {
                throw new IllegalStateException("sample has no minimized molecule");
            }
            int count = epmol.pmol.dofs.size();
            if ((epmol.params == null && count != 0)
                    || (epmol.params != null && epmol.params.size() != count)) {
                throw new IllegalStateException("minimized DOF vector does not match molecule");
            }
            // Validate the entire vector before changing any coordinates.
            if (epmol.params != null) {
                for (int i = 0; i < count; i++) {
                    if (!Double.isFinite(epmol.params.get(i))) {
                        throw new IllegalStateException("non-finite minimized DOF at index " + i);
                    }
                }
                for (int i = 0; i < count; i++) {
                    epmol.pmol.dofs.get(i).apply(epmol.params.get(i));
                }
            }
            restoredMolecule = epmol.pmol.mol;
            return restoredMolecule;
        }

        /** Return the underlying energy result for advanced predicates. */
        public EnergyCalculator.EnergiedParametricMolecule
        getEnergiedParametricMolecule() {
            return epmol;
        }

        /** Return the number of stored minimized DOF values, or zero. */
        public int getNumDofs() {
            return epmol == null || epmol.params == null ? 0 : epmol.params.size();
        }

        /** Return one stored minimized DOF value. */
        public double getDofValue(int index) {
            if (epmol == null || epmol.params == null) {
                throw new IllegalStateException("this sample has no minimized DOF vector");
            }
            return epmol.params.get(index);
        }
    }
}
