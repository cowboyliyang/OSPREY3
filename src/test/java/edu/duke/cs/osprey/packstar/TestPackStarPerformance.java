package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.*;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.ByteBoundedCache;
import org.junit.jupiter.api.Test;

import java.lang.reflect.*;
import java.util.*;
import static org.junit.jupiter.api.Assertions.*;

/** Exercise real CCD, ordered sample multiplicities, and proposal-cache invalidation. */
public class TestPackStarPerformance {
    private static Object field(Object value, String name) throws Exception {
        Field field = value.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.get(value);
    }

    private static Object invoke(PackStarEstimator estimator, String name,
                                 Class<?>[] types, Object... args) throws Exception {
        Method method = PackStarEstimator.class.getDeclaredMethod(name, types);
        method.setAccessible(true);
        try { return method.invoke(estimator, args); }
        catch (InvocationTargetException ex) { throw new RuntimeException(ex.getCause()); }
    }

    @SuppressWarnings("unchecked")
    private static List<Object> ccd(PackStarEstimator estimator, List<int[]> samples,
                                     EnergyMatrix proposal, double logZ, String id,
                                     boolean pose) throws Exception {
        return (List<Object>) invoke(estimator, "runParallelCCD",
                new Class<?>[]{List.class, EnergyMatrix.class, double.class, String.class, boolean.class},
                samples, proposal, logZ, id, pose);
    }

    private static final class PropertiesScope implements AutoCloseable {
        private final Properties saved = (Properties) System.getProperties().clone();
        void set(String key, String value) { System.setProperty(key, value); }
        @Override public void close() { System.setProperties(saved); }
    }

    @Test public void targetCachePreservesMultiplicityAndRefreshesProposalProvenance() throws Exception {
        try (var properties = new PropertiesScope()) {
            properties.set("packstar.pac.ccd.cacheBytes", "1048576");
            properties.set("packstar.pac.ccd.compactResults", "true");
            properties.set("packstar.pac.ccd.submissionBatchSize", "2");
            properties.set("packstar.pac.ccd.streaming", "true");
            Strand strand = new Strand.Builder(PDBIO.readFile("src/test/resources/1CC8.ss.pdb")).build();
            for (String residue : List.of("A16", "A19")) {
                strand.flexibility.get(residue).setLibraryRotamers(Strand.WildType)
                        .addWildTypeRotamers().setContinuous();
            }
            var space = new SimpleConfSpace.Builder().addStrand(strand).build();
            var parallelism = Parallelism.makeCpu(2);
            try (var energy = new EnergyCalculator.Builder(space, new ForcefieldParams())
                    .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true)
                    .setParallelism(parallelism).build()) {
                var calculator = new ConfEnergyCalculator.Builder(space, energy).build();
                var base = new SimplerEnergyMatrixCalculator.Builder(calculator).build().calcEnergyMatrix();
                var graph = InteractionGraph.buildFromEdges(2, List.of(new int[]{0, 1}));
                var estimator = new PackStarEstimator(null, null, base, base, graph,
                        calculator, new RCs(space), space);
                List<int[]> samples = List.of(new int[]{0, 0}, new int[]{1, 0}, new int[]{0, 0});
                List<Object> first = ccd(estimator, samples, base, 1.25, "old", false);
                assertEquals(2, estimator.getTotalCCDCalls());
                assertSame(first.get(0), first.get(2));
                for (Object sample : first) {
                    assertNull(field(sample, "epmol"));
                    assertNotNull(field(sample, "features"));
                }
                var proposal = new EnergyMatrix(base);
                proposal.setOneBody(0, 0, proposal.getOneBody(0, 0) + 0.75);
                List<Object> second = ccd(estimator, samples, proposal, 9.0, "new", false);
                assertEquals(2, estimator.getTotalCCDCalls(), "cross-batch hits must avoid CCD");
                assertEquals(samples.size(), second.size());
                assertNotSame(first.get(0), second.get(0), "provenance belongs to each batch");
                for (int i = 0; i < samples.size(); i++) {
                    assertEquals(field(first.get(i), "eTrue"), field(second.get(i), "eTrue"));
                    assertEquals("old", field(first.get(i), "sourceProposalId"));
                    assertEquals("new", field(second.get(i), "sourceProposalId"));
                    assertEquals(9.0, field(second.get(i), "sourceProposalLogZ"));
                    assertSame(field(first.get(i), "features"), field(second.get(i), "features"));
                }
                assertEquals((double) field(first.get(0), "sourceProposalEnergy") + 0.75,
                        (double) field(second.get(0), "sourceProposalEnergy"), 1e-12);
                List<Object> geometry = ccd(estimator, samples, null, Double.NaN, null, true);
                assertEquals(4, estimator.getTotalCCDCalls(), "geometric events need actual poses");
                var pose = (EnergyCalculator.EnergiedParametricMolecule) field(geometry.get(0), "epmol");
                assertNotNull(new PackStarFunctionalEvent.Sample(0, samples.get(0),
                        (double) field(geometry.get(0), "eTrue"),
                        (double) field(geometry.get(0), "eMin"), pose).getMolecule());
                invoke(estimator, "shutdownSamplingResources", new Class<?>[0]);
                assertEquals(0, ((ByteBoundedCache<?, ?>) field(estimator, "ccdEnergyCache")).bytes());
            }
        }
    }

    @Test public void cdfReusePreservesSeededDrawsAndInvalidatesInPlaceRefits() throws Exception {
        for (int threads : new int[]{1, 2}) {
            try (var properties = new PropertiesScope()) {
                properties.set("packstar.pac.sampling.threads", Integer.toString(threads));
                properties.set("packstar.pac.sampling.largeLambdaThreshold", "1");
                properties.set("packstar.pac.sampling.gpu", "false");
                properties.set("packstar.pac.sampling.cdfCacheBytes", "1048576");
                properties.set("branchdp.dp.gpu", "false");
                var rcs = new RCs(new int[][]{{0, 1}, {0, 1}, {0, 1}});
                var base = new EnergyMatrix(3, new int[]{2, 2, 2}, 0);
                var graph = InteractionGraph.buildFromEdges(3, List.of(new int[]{0, 1}, new int[]{1, 2}));
                var root = PackStarTripleDecompositionCosts.rootProposalGraph(graph, rcs, null, true).selected.root;
                RootedTreeEdge.postOrderInitIncremental(root, base, base, graph, 1.0);
                RootedTreeEdge.postOrderComputeFullDP(root);
                var estimator = new PackStarEstimator(root, root.getLeftChild().getChildOfEdge(),
                        base, base, graph, null, rcs, null);
                try {
                    List<int[]> first = draws(estimator, 123L);
                    List<int[]> second = draws(estimator, 123L);
                    equalDraws(first, second);
                    var cache = (ByteBoundedCache<?, ?>) field(estimator, "conditionalCdfCache");
                    assertTrue(cache.hits() > 0);
                    base.setOneBody(0, 0, 1000.0);
                    invoke(estimator, "recomputeDP", new Class<?>[]{EnergyMatrix.class}, base);
                    assertEquals(0, cache.bytes());
                    List<int[]> changed = draws(estimator, 123L);
                    assertTrue(changed.stream().allMatch(conf -> conf[0] == 1));
                    invoke(estimator, "invalidateSamplingCaches", new Class<?>[0]);
                    equalDraws(changed, draws(estimator, 123L));
                } finally {
                    invoke(estimator, "shutdownSamplingResources", new Class<?>[0]);
                    RootedTreeEdge.postOrderReleaseLargeMemory(root);
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static List<int[]> draws(PackStarEstimator estimator, long seed) throws Exception {
        return (List<int[]>) invoke(estimator, "sampleConformationsFromDP",
                new Class<?>[]{int.class, Random.class}, 128, new Random(seed));
    }

    private static void equalDraws(List<int[]> first, List<int[]> second) {
        assertEquals(first.size(), second.size());
        for (int i = 0; i < first.size(); i++) assertArrayEquals(first.get(i), second.get(i));
    }
}
