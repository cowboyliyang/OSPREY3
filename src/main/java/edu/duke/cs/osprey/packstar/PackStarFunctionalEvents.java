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

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

import java.util.Arrays;
import java.util.Objects;

/** Common discrete and geometric PACK* functional-event predicates. */
public final class PackStarFunctionalEvents {

    private PackStarFunctionalEvents() {
    }

    /** Match one complete discrete conformation assignment. */
    public static PackStarFunctionalEvent conformationEquals(int[] expected) {
        Objects.requireNonNull(expected, "expected conformation");
        int[] copy = expected.clone();
        return sample -> Arrays.equals(sample.getConformation(), copy);
    }

    /** Match one RC at one position. */
    public static PackStarFunctionalEvent conformationAt(int position, int rc) {
        if (position < 0) {
            throw new IllegalArgumentException("position must be nonnegative");
        }
        return sample -> position < sample.getConformation().length
                && sample.getRC(position) == rc;
    }

    /** Combine events with a logical AND. */
    public static PackStarFunctionalEvent allOf(PackStarFunctionalEvent... events) {
        Objects.requireNonNull(events, "events");
        PackStarFunctionalEvent[] copy = events.clone();
        for (PackStarFunctionalEvent event : copy) {
            Objects.requireNonNull(event, "event");
        }
        return sample -> {
            for (PackStarFunctionalEvent event : copy) {
                if (!event.test(sample)) {
                    return false;
                }
            }
            return true;
        };
    }

    /** Combine events with a logical OR. */
    public static PackStarFunctionalEvent anyOf(PackStarFunctionalEvent... events) {
        Objects.requireNonNull(events, "events");
        PackStarFunctionalEvent[] copy = events.clone();
        for (PackStarFunctionalEvent event : copy) {
            Objects.requireNonNull(event, "event");
        }
        return sample -> {
            for (PackStarFunctionalEvent event : copy) {
                if (event.test(sample)) {
                    return true;
                }
            }
            return false;
        };
    }

    /** Negate an event. */
    public static PackStarFunctionalEvent negate(PackStarFunctionalEvent event) {
        Objects.requireNonNull(event, "event");
        return sample -> !event.test(sample);
    }

    /**
     * Match an atom-atom distance interval in the minimized molecule.
     * Bounds are inclusive and measured in angstroms.
     * Missing geometry raises an event error, rather than counting as a miss.
     * A sequence-dependent loss of an atom must be handled explicitly by the
     * caller's event definition before invoking this strict geometric predicate.
     */
    public static PackStarFunctionalEvent atomDistanceWithin(
            String residue1, String atom1,
            String residue2, String atom2,
            double minDistance, double maxDistance) {
        Objects.requireNonNull(residue1, "residue1");
        Objects.requireNonNull(atom1, "atom1");
        Objects.requireNonNull(residue2, "residue2");
        Objects.requireNonNull(atom2, "atom2");
        if (!Double.isFinite(minDistance) || !Double.isFinite(maxDistance)
                || minDistance < 0.0
                || minDistance > maxDistance) {
            throw new IllegalArgumentException(
                    "distance bounds must be finite, nonnegative and ordered");
        }

        return sample -> {
            Molecule molecule = sample.getMolecule();
            Residue res1 = molecule.getResByPDBResNumberOrNull(residue1);
            Residue res2 = molecule.getResByPDBResNumberOrNull(residue2);
            if (res1 == null || res2 == null) {
                throw new IllegalStateException("event residue missing: "
                        + (res1 == null ? residue1 : residue2));
            }
            Atom first = res1.getAtomByName(atom1);
            Atom second = res2.getAtomByName(atom2);
            if (first == null || second == null) {
                throw new IllegalStateException("event atom missing: "
                        + (first == null ? residue1 + ":" + atom1 : residue2 + ":" + atom2));
            }
            double[] firstCoords = first.getCoords();
            double[] secondCoords = second.getCoords();
            if (firstCoords == null || secondCoords == null
                    || firstCoords.length != 3 || secondCoords.length != 3) {
                throw new IllegalStateException("event coordinates missing or malformed");
            }
            for (int i = 0; i < 3; i++) {
                if (!Double.isFinite(firstCoords[i]) || !Double.isFinite(secondCoords[i])) {
                    throw new IllegalStateException("event coordinates are non-finite");
                }
            }
            double dx = firstCoords[0] - secondCoords[0];
            double dy = firstCoords[1] - secondCoords[1];
            double dz = firstCoords[2] - secondCoords[2];
            double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (!Double.isFinite(distance)) {
                throw new IllegalStateException("event distance is non-finite");
            }
            return distance >= minDistance
                    && distance <= maxDistance;
        };
    }
}
