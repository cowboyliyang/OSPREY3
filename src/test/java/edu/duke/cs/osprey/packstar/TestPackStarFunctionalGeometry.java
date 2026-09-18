package edu.duke.cs.osprey.packstar;

import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.markstar.bench.ConfSpaces2RL0;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/** Pre-freeze software checks, not functional-design or coverage evidence. */
public class TestPackStarFunctionalGeometry {

    private static PackStarFunctionalEvent.Sample sample(
            EnergyCalculator.EnergiedParametricMolecule energy) {
        return new PackStarFunctionalEvent.Sample(0, new int[]{0}, 0, 0, energy);
    }

    private static class CountingDof extends DegreeOfFreedom {
        int calls;
        public void apply(double value) { calls++; }
        public DOFBlock getBlock() { return null; }
        public String getName() { return "count"; }
    }

    @Test
    public void missingGeometryIsAnEventErrorIncludingUnderNegation() {
        PackStarFunctionalEvent distance = PackStarFunctionalEvents.atomDistanceWithin(
                "A1", "CA", "A2", "CA", 0, 5);
        for (PackStarFunctionalEvent event : List.of(distance,
                PackStarFunctionalEvents.negate(distance))) {
            PackStarFunctionalObservableResult result =
                    PackStarFunctionalObservableEvaluator.evaluate(
                            "missing", event, 0.6, new double[]{0}, i -> sample(null));
            assertEquals(PackStarFunctionalObservableResult.Status.UNRESOLVED_EVENT_ERROR,
                    result.getStatus());
            assertTrue(Double.isNaN(result.getProbability()));
            assertTrue(Double.isNaN(result.getRawHitFraction()));
            assertFalse(result.isResolved());
        }
    }

    @Test
    public void invalidDofVectorIsRejectedBeforeAnyCoordinatesChange() {
        CountingDof first = new CountingDof();
        CountingDof second = new CountingDof();
        ParametricMolecule pmol = new ParametricMolecule(
                new Molecule(), List.of(first, second), null);
        for (double[] values : List.of(new double[]{1}, new double[]{1, Double.NaN},
                new double[]{1, Double.POSITIVE_INFINITY})) {
            assertThrows(IllegalStateException.class, () -> sample(
                    new EnergyCalculator.EnergiedParametricMolecule(
                            pmol, null, new DenseDoubleMatrix1D(values), 0)).getMolecule());
        }
        assertThrows(IllegalStateException.class, () -> sample(
                new EnergyCalculator.EnergiedParametricMolecule(pmol, null, 0)).getMolecule());
        assertEquals(0, first.calls);
        assertEquals(0, second.calls);
    }

    @Test
    public void geometryIsRestoredOncePerLogicalSample() {
        CountingDof dof = new CountingDof();
        ParametricMolecule pmol = new ParametricMolecule(new Molecule(), List.of(dof), null);
        EnergyCalculator.EnergiedParametricMolecule energy =
                new EnergyCalculator.EnergiedParametricMolecule(
                        pmol, null, new DenseDoubleMatrix1D(new double[]{10}), 0);
        PackStarFunctionalEvent.Sample first = sample(energy);
        assertSame(first.getMolecule(), first.getMolecule());
        assertEquals(1, dof.calls);
        sample(energy).getMolecule();
        assertEquals(2, dof.calls, "repeated draws must restore their own evaluation");
    }

    @Test
    public void rigidGeometryNeedsNoDofVectorAndMissingResidueIsAnError() {
        PackStarFunctionalEvent.Sample rigid = sample(
                new EnergyCalculator.EnergiedParametricMolecule(
                        new ParametricMolecule(new Molecule()), null, 0));
        assertNotNull(rigid.getMolecule());
        assertThrows(IllegalStateException.class, () ->
                PackStarFunctionalEvents.atomDistanceWithin("A1", "CA", "A2", "CA", 0, 5)
                        .test(rigid));
    }

    @Test
    public void distanceBoundsMustBePhysicalAndFinite() {
        for (double[] bounds : List.of(new double[]{-1, 2}, new double[]{2, 1},
                new double[]{0, Double.POSITIVE_INFINITY}, new double[]{Double.NaN, 2})) {
            assertThrows(IllegalArgumentException.class, () ->
                    PackStarFunctionalEvents.atomDistanceWithin(
                            "A1", "CA", "A2", "CA", bounds[0], bounds[1]));
        }
    }

    @Test
    public void dpEarlyExitDoesNotClaimAnEventEstimate() {
        SimpleConfSpace space = ConfSpaces2RL0.buildWildTypeConfSpace(2).protein;
        Parallelism parallelism = Parallelism.makeCpu(1);
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(space, new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(false)
                .setParallelism(parallelism).build()) {
            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(space, ecalc).build();
            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                    .build().calcEnergyMatrix();
            try (PackStarPartitionFunction bound = new PackStarPartitionFunction(
                    space, emat, emat, confEcalc, new RCs(space), parallelism, "pc-rigid-regression")) {
                bound.setFunctionalEvent("must-not-evaluate", s -> {
                    fail("DP exit must not fabricate a final sample");
                    return true;
                });
                bound.init(0.1);
                bound.compute(100);
                PackStarResult saved = bound.makeResult();
                bound.close();
                PackStarFunctionalObservableResult result = saved.getFunctionalObservableResult();
                assertEquals(PackStarFunctionalObservableResult.Status.NOT_COMPUTED, result.getStatus());
                assertTrue(Double.isNaN(result.getProbability()));
                assertEquals(0, result.getSampleCount());
                assertTrue(result.getDiagnostic().contains("initial DP bounds"), result.getDiagnostic());
            }
        }
    }

    @Test
    public void realCcdRestorationAndEnumeratedProbabilityAgreeWithIndependentReference() {
        // Isolated human protein strand, two WT positions, no sequence design.
        // The distance grid is a numerical fixture, not a biological event.
        SimpleConfSpace space = ConfSpaces2RL0.buildWildTypeConfSpace(2).protein;
        RCs rcs = new RCs(space);
        assertEquals(2, rcs.getNumPos());
        int count = rcs.getNumConformations().intValueExact();
        assertTrue(count > 1 && count <= 256, "bounded exhaustive regression fixture");
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(space, new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true)
                .setParallelism(Parallelism.makeCpu(1)).build()) {
            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(space, ecalc).build();
            List<EnergyCalculator.EnergiedParametricMolecule> energies = new ArrayList<>();
            List<int[]> assignments = new ArrayList<>();
            List<Double> distances = new ArrayList<>();
            double minEnergy = Double.POSITIVE_INFINITY;
            for (int a : rcs.get(0)) for (int b : rcs.get(1)) {
                int[] assignment = {a, b};
                EnergyCalculator.EnergiedParametricMolecule energy =
                        confEcalc.calcEnergy(new RCTuple(assignment));
                assertTrue(Double.isFinite(energy.energy));
                assertNotNull(energy.params);
                assertTrue(energy.params.size() > 0);
                // Independent reference: apply the saved vector directly.
                for (int i = 0; i < energy.params.size(); i++) {
                    energy.pmol.dofs.get(i).apply(energy.params.get(i));
                }
                double[] expected = energy.pmol.mol.getResByPDBResNumberOrNull("A156").coords.clone();
                distances.add(distance(energy.pmol.mol));
                try (EnergyFunction function = ecalc.makeEnergyFunction(energy)) {
                    assertEquals(energy.energy, function.getEnergy(), 1e-6);
                }
                // Make stale coordinates deliberately; the public Sample must restore CCD.
                energy.pmol.dofs.get(0).apply(energy.params.get(0) + 5.0);
                Molecule restored = sample(energy).getMolecule();
                assertArrayEquals(expected, restored.getResByPDBResNumberOrNull("A156").coords, 1e-8);
                try (EnergyFunction function = ecalc.makeEnergyFunction(energy)) {
                    assertEquals(energy.energy, function.getEnergy(), 1e-6);
                }
                energies.add(energy);
                assignments.add(assignment);
                minEnergy = Math.min(minEnergy, energy.energy);
            }
            assertEquals(count, energies.size());
            double rt = 0.592;
            double[] logWeights = new double[count];
            double total = 0;
            for (int i = 0; i < count; i++) {
                logWeights[i] = (minEnergy - energies.get(i).energy) / rt;
                total += Math.exp(logWeights[i]);
            }
            boolean mixed = false;
            for (double cutoff : new double[]{4, 6, 8, 10, 12}) {
                double hit = 0;
                int hits = 0;
                for (int i = 0; i < count; i++) if (distances.get(i) <= cutoff) {
                    hit += Math.exp(logWeights[i]);
                    hits++;
                }
                mixed |= hits > 0 && hits < count;
                PackStarFunctionalEvent event = PackStarFunctionalEvents.atomDistanceWithin(
                        "A156", "CZ", "A157", "CA", 0, cutoff);
                PackStarFunctionalObservableResult result = PackStarFunctionalObservableEvaluator.evaluate(
                        "numerical-distance", event, rt, logWeights, i ->
                                new PackStarFunctionalEvent.Sample(i, assignments.get(i),
                                        energies.get(i).energy, 0, energies.get(i)));
                assertEquals(hit / total, result.getProbability(), 1e-12);
                assertEquals(hits, result.getHitCount());
                // A nonuniform deterministic census with inverse-frequency weights
                // checks repeated-draw accounting; this is not a coverage experiment.
                List<Integer> draws = new ArrayList<>();
                List<Double> repeatedWeights = new ArrayList<>();
                for (int i = 0; i < count; i++) for (int j = 0; j < 1 + i % 3; j++) {
                    draws.add(i);
                    repeatedWeights.add(logWeights[i] - Math.log(1 + i % 3));
                }
                PackStarFunctionalObservableResult repeated = PackStarFunctionalObservableEvaluator.evaluate(
                        "repeated-distance", event, rt,
                        repeatedWeights.stream().mapToDouble(Double::doubleValue).toArray(), i -> {
                            int k = draws.get(i);
                            return new PackStarFunctionalEvent.Sample(i, assignments.get(k),
                                    energies.get(k).energy, 0, energies.get(k));
                        });
                assertEquals(hit / total, repeated.getProbability(), 1e-12);
                assertEquals(draws.size(), repeated.getSampleCount());
                System.out.println("[PC-GEOMETRY-REFERENCE] count=" + count + " cutoff=" + cutoff
                        + " hits=" + hits + " reference=" + hit / total
                        + " observed=" + result.getProbability());
            }
            assertTrue(mixed, "at least one predeclared threshold must split the fixture");
            assertThrows(IllegalStateException.class, () ->
                    PackStarFunctionalEvents.atomDistanceWithin("A156", "MISSING", "A157", "CA", 0, 5)
                            .test(sample(energies.get(0))));
            energies.get(0).pmol.mol.getResByPDBResNumberOrNull("A156").coords[0] = Double.NaN;
            assertThrows(IllegalStateException.class, () ->
                    PackStarFunctionalEvents.atomDistanceWithin("A156", "N", "A157", "CA", 0, 5)
                            .test(sample(energies.get(0))));
        }
    }

    private static double distance(Molecule molecule) {
        double[] a = molecule.getResByPDBResNumberOrNull("A156").getAtomByName("CZ").getCoords();
        double[] b = molecule.getResByPDBResNumberOrNull("A157").getAtomByName("CA").getCoords();
        double squared = 0;
        for (int i = 0; i < 3; i++) squared += (a[i] - b[i]) * (a[i] - b[i]);
        return Math.sqrt(squared);
    }
}
