package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

public class TestPackStarSampleAccounting {
    private Path output;

    private static PackStarResult result(long records, long callbacks) {
        var source = PartitionFunction.Result.makeAborted()
                .setStat(PackStarResult.CCD_SAMPLE_RECORDS_STAT, records)
                .setStat(PackStarResult.SAMPLE_CALLBACKS_STAT, callbacks);
        return new PackStarResult(source, PackStarFunctionalObservableResult.notConfigured());
    }

    @Test public void coverageRejectsMissingPartialAndExtraTraces() {
        assertThrows(IllegalStateException.class,
                () -> PackStarSmallCoverage.checkSampleEnergies(result(180, 0), 0, 0));
        assertThrows(IllegalStateException.class,
                () -> PackStarSmallCoverage.checkSampleEnergies(result(180, 179), 179, 0));
        assertThrows(IllegalStateException.class,
                () -> PackStarSmallCoverage.checkSampleEnergies(result(180, 181), 181, 0));
        var missing = new PackStarResult(PartitionFunction.Result.makeAborted(),
                PackStarFunctionalObservableResult.notConfigured());
        assertThrows(IllegalStateException.class,
                () -> PackStarSmallCoverage.checkSampleEnergies(missing, 0, 0));
    }

    @Test public void coverageDistinguishesUncheckedFromPassedAndRejectsBadEnergies() {
        assertEquals("NOT_SAMPLED", PackStarSmallCoverage.checkSampleEnergies(result(0, 0), 0, 0));
        assertEquals("PASSED", PackStarSmallCoverage.checkSampleEnergies(result(180, 180), 180, 1e-8));
        for (double error : new double[]{1e-5, Double.NaN, Double.POSITIVE_INFINITY}) {
            assertThrows(IllegalStateException.class,
                    () -> PackStarSmallCoverage.checkSampleEnergies(result(180, 180), 180, error));
        }
    }

    @Test public void realCcdAccountingSurvivesCacheHitsAbortAndNoSampling() throws Exception {
        // Retain diagnostic artifacts with the Slurm run instead of deleting them
        // during JUnit teardown (which also races NFS directory visibility).
        output = Files.createTempDirectory(Path.of(System.getProperty("java.io.tmpdir")),
                "packstar-accounting-");
        System.out.println("ACCOUNTING_ARTIFACTS=" + output);
        Properties saved = (Properties) System.getProperties().clone();
        try {
            System.setProperty("packstar.cutoff.strategy", "COMPLETE");
            System.setProperty("packstar.dp.gpu", "false");
            System.setProperty("packstar.dp.parallel.threads", "2");
            System.setProperty("packstar.pac.sampling.threads", "2");
            System.setProperty("packstar.pac.sampling.gpu", "false");
            System.setProperty("packstar.pac.sampling.progress", "false");
            System.setProperty("packstar.pac.trainSamples", "200");
            System.setProperty("packstar.pac.maxEstSamples", "4000");
            System.setProperty("packstar.pac.frequencySeverity.maxRefits", "2");
            System.setProperty("packstar.pac.frequencySeverity.tripleEta", "false");
            System.setProperty("packstar.pac.randomSeed", "42");
            System.setProperty("packstar.pac.ccd.cacheBytes", "67108864");
            Strand strand = new Strand.Builder(PDBIO.readFile("src/test/resources/1CC8.ss.pdb")).build();
            for (String residue : List.of("A16", "A19")) {
                strand.flexibility.get(residue).setLibraryRotamers(Strand.WildType)
                        .addWildTypeRotamers().setContinuous();
            }
            var space = new SimpleConfSpace.Builder().addStrand(strand).build();
            var parallelism = Parallelism.makeCpu(2);
            try (var energy = new EnergyCalculator.Builder(space, new ForcefieldParams())
                    .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true)
                    .setParallelism(parallelism).build();
                 var rigid = new EnergyCalculator.SharedBuilder(energy).setIsMinimizing(false).build()) {
                var calculator = new ConfEnergyCalculator.Builder(space, energy).build();
                var min = new SimplerEnergyMatrixCalculator.Builder(calculator).build().calcEnergyMatrix();
                // Avoid an unrelated ~1e-14 lower-bound cancellation at an identical pose.
                for (int pos = 0; pos < min.getNumPos(); pos++) {
                    for (int rc = 0; rc < min.getNumConfAtPos(pos); rc++) {
                        min.setOneBody(pos, rc, min.getOneBody(pos, rc) - 1e-8);
                    }
                }
                var rig = new SimplerEnergyMatrixCalculator.Builder(new ConfEnergyCalculator(calculator, rigid))
                        .build().calcEnergyMatrix();
                PackStarResult untraced = run("untraced", space, calculator, parallelism, min, rig, 0.683, null);
                assertEquals(PartitionFunction.Status.Estimated, untraced.status);
                assertTrue(untraced.numConfs > 0);
                assertEquals(0, untraced.getStat(PackStarResult.SAMPLE_CALLBACKS_STAT));

                List<PackStarSampleTrace> cachedTraces = new ArrayList<>();
                PackStarResult cached = run("cached", space, calculator, parallelism, min, rig, 0.683, cachedTraces);
                sameEstimate(untraced, cached);
                assertEquals(untraced.numConfs, cached.numConfs);
                Set<PackStarSampleTrace.Stage> stages = EnumSet.noneOf(PackStarSampleTrace.Stage.class);
                Map<String, Double> independentEnergies = new HashMap<>();
                double maxError = 0;
                double rt = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
                for (int i = 0; i < cachedTraces.size(); i++) {
                    PackStarSampleTrace trace = cachedTraces.get(i);
                    assertEquals(i, trace.sampleIndex);
                    stages.add(trace.stage);
                    double expected = independentEnergies.computeIfAbsent(Arrays.toString(trace.getConf()),
                            key -> calculator.calcEnergy(new RCTuple(trace.getConf())).energy);
                    maxError = Math.max(maxError, Math.abs(expected - trace.eTrue));
                    assertEquals(trace.eTrue - trace.eProposal, trace.xi, 1e-12);
                    assertEquals(-trace.xi / rt, trace.logWeight, 1e-12);
                    assertTrue(Double.isFinite(trace.logZCorrected));
                    assertEquals(trace.logWeight > trace.clipLogCap, trace.clipped);
                    if (trace.stage == PackStarSampleTrace.Stage.TRAIN) assertTrue(Double.isNaN(trace.clipLogCap));
                    else assertTrue(Double.isFinite(trace.clipLogCap));
                }
                assertEquals(EnumSet.of(PackStarSampleTrace.Stage.TRAIN, PackStarSampleTrace.Stage.DISCOVERY,
                        PackStarSampleTrace.Stage.VALIDATION, PackStarSampleTrace.Stage.ESTIMATION), stages);
                assertEquals("PASSED", PackStarSmallCoverage.checkSampleEnergies(cached, cachedTraces.size(), maxError));
                assertTrue(cachedTraces.size() > independentEnergies.size(), "duplicate draws must remain visible");

                System.setProperty("packstar.pac.ccd.cacheBytes", "0");
                List<PackStarSampleTrace> uncachedTraces = new ArrayList<>();
                PackStarResult uncached = run("uncached", space, calculator, parallelism, min, rig, 0.683, uncachedTraces);
                sameEstimate(cached, uncached);
                assertTrue(cached.numConfs < uncached.numConfs, "cache hits reduce actual CCD, not callbacks");
                assertEquals(cachedTraces.size(), uncachedTraces.size());
                for (int i = 0; i < cachedTraces.size(); i++) {
                    var a = cachedTraces.get(i);
                    var b = uncachedTraces.get(i);
                    assertArrayEquals(a.getConf(), b.getConf());
                    assertEquals(a.stage, b.stage);
                    assertEquals(a.eTrue, b.eTrue);
                    assertEquals(a.eProposal, b.eProposal);
                    assertEquals(a.logZCorrected, b.logZCorrected);
                }

                System.setProperty("packstar.pac.maxEstSamples", "16");
                System.setProperty("packstar.pac.frequencySeverity.maxRefits", "0");
                List<PackStarSampleTrace> abortedTraces = new ArrayList<>();
                PackStarResult aborted = run("aborted", space, calculator, parallelism, min, rig, 1e-8, abortedTraces);
                assertEquals(PartitionFunction.Status.Aborted, aborted.status);
                assertTrue(aborted.numConfs > 0);
                assertTrue(abortedTraces.size() > 0);

                System.setProperty("packstar.pac.frequencySeverity.outputDir", output.resolve("listener-error").toString());
                var listenerError = new IllegalStateException("deliberate energy-check failure");
                try (var contexts = energy.tasks.contextGroup();
                     var pf = new PackStarPartitionFunction(space, rig, min, calculator,
                             new RCs(space), parallelism, "performance-2")) {
                    pf.setSampleListener(trace -> { throw listenerError; });
                    pf.init(1e-8);
                    pf.putTaskContexts(contexts);
                    assertSame(listenerError, assertThrows(IllegalStateException.class,
                            () -> pf.compute(Integer.MAX_VALUE)));
                    assertTrue(pf.getNumConfsEvaluated() > 0);
                    PackStarResult interrupted = pf.makeResult();
                    assertEquals(pf.getNumConfsEvaluated(), interrupted.numConfs);
                    assertEquals(200, interrupted.getStat(PackStarResult.CCD_SAMPLE_RECORDS_STAT));
                    assertEquals(0, interrupted.getStat(PackStarResult.SAMPLE_CALLBACKS_STAT));
                    assertThrows(IllegalStateException.class,
                            () -> PackStarSmallCoverage.checkSampleEnergies(interrupted, 0, 0));
                }

                List<PackStarSampleTrace> skippedTraces = new ArrayList<>();
                PackStarResult skipped = run("skipped", space, calculator, parallelism, min, min, 0.683, skippedTraces);
                assertEquals(PartitionFunction.Status.Estimated, skipped.status);
                assertEquals(0, skipped.numConfs);
                assertTrue(skippedTraces.isEmpty());
                assertEquals("NOT_SAMPLED", PackStarSmallCoverage.checkSampleEnergies(skipped, 0, 0));
            }
        } finally {
            System.setProperties(saved);
        }
    }

    private PackStarResult run(String name, SimpleConfSpace space, ConfEnergyCalculator calculator,
                               Parallelism parallelism, EnergyMatrix min, EnergyMatrix rig,
                               double epsilon, List<PackStarSampleTrace> traces) throws Exception {
        System.setProperty("packstar.pac.frequencySeverity.outputDir", output.resolve(name).toString());
        Thread computeThread = Thread.currentThread();
        try (var contexts = calculator.ecalc.tasks.contextGroup();
             var pf = new PackStarPartitionFunction(space, rig, min, calculator,
                     new RCs(space), parallelism, "performance-2")) {
            if (traces != null) pf.setSampleListener(trace -> {
                assertSame(computeThread, Thread.currentThread());
                assertTrue(pf.getNumConfsEvaluated() > 0);
                traces.add(trace);
                // A consumer cannot corrupt the estimator's assignment array.
                trace.getConf()[0] = -999;
            });
            pf.init(epsilon);
            pf.putTaskContexts(contexts);
            pf.compute(Integer.MAX_VALUE);
            PackStarResult result = pf.makeResult();
            assertEquals(pf.getNumConfsEvaluated(), result.numConfs);
            assertEquals(result.numConfs, result.getStat(PackStarResult.FULL_CCD_STAT));
            if (traces != null) {
                assertEquals(traces.size(), result.getStat(PackStarResult.CCD_SAMPLE_RECORDS_STAT));
                assertEquals(traces.size(), result.getStat(PackStarResult.SAMPLE_CALLBACKS_STAT));
            }
            System.out.println("ACCOUNTING " + name + " status=" + result.status
                    + " CCD=" + result.numConfs
                    + " records=" + result.getStat(PackStarResult.CCD_SAMPLE_RECORDS_STAT)
                    + " callbacks=" + result.getStat(PackStarResult.SAMPLE_CALLBACKS_STAT));
            return result;
        }
    }

    private static void sameEstimate(PackStarResult a, PackStarResult b) {
        assertEquals(a.status, b.status);
        assertEquals(a.values.calcLowerBound(), b.values.calcLowerBound());
        assertEquals(a.values.calcUpperBound(), b.values.calcUpperBound());
    }
}
