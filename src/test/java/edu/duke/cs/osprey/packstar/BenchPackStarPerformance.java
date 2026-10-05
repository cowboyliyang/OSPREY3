package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import java.lang.management.ManagementFactory;
import java.nio.file.*;
import java.util.*;

/** Paired benchmark using only APIs also present at the pre-optimization baseline. */
public final class BenchPackStarPerformance {
    public static void main(String[] args) throws Exception {
        Path out = Path.of(args[0]);
        String mode = args[1];
        int positions = Integer.parseInt(args[2]);
        long seed = Long.parseLong(args[3]);
        boolean optimized = mode.equals("optimized");
        Files.createDirectories(out);
        System.setProperty("packstar.pac.ccd.cacheBytes", optimized ? "67108864" : "0");
        System.setProperty("packstar.pac.sampling.cdfCacheBytes", optimized ? "67108864" : "0");
        System.setProperty("packstar.pac.sampling.gpu.preparedCacheBytes", optimized ? "67108864" : "0");
        System.setProperty("packstar.pac.ccd.compactResults", Boolean.toString(optimized));
        System.setProperty("packstar.pac.ccd.streaming", Boolean.toString(optimized));
        System.setProperty("packstar.pac.ccd.submissionBatchSize", "32");
        System.setProperty("packstar.cutoff.strategy", "COMPLETE");
        System.setProperty("packstar.dp.gpu", "false");
        System.setProperty("packstar.dp.parallel.threads", "8");
        System.setProperty("packstar.pac.sampling.threads", "8");
        System.setProperty("packstar.pac.sampling.progress", "false");
        System.setProperty("packstar.pac.trainSamples", "200");
        System.setProperty("packstar.pac.maxEstSamples", "4000");
        System.setProperty("packstar.pac.frequencySeverity.maxRefits", "2");
        System.setProperty("packstar.pac.frequencySeverity.tripleEta", "false");
        System.setProperty("packstar.pac.randomSeed", Long.toString(seed));

        Strand strand = new Strand.Builder(PDBIO.readFile("src/test/resources/1CC8.ss.pdb")).build();
        String[] residues = {"A16", "A19", "A23", "A25"};
        for (int i = 0; i < positions; i++) {
            strand.flexibility.get(residues[i]).setLibraryRotamers(Strand.WildType)
                    .addWildTypeRotamers().setContinuous();
        }
        var space = new SimpleConfSpace.Builder().addStrand(strand).build();
        var parallelism = Parallelism.makeCpu(8);
        try (var energy = new EnergyCalculator.Builder(space, new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true).setParallelism(parallelism).build();
             var rigid = new EnergyCalculator.SharedBuilder(energy).setIsMinimizing(false).build()) {
            var calculator = new ConfEnergyCalculator.Builder(space, energy).build();
            var min = new SimplerEnergyMatrixCalculator.Builder(calculator).build().calcEnergyMatrix();
            // Benchmark input only: round the lower envelope conservatively.
            // The unmodified baseline otherwise aborts on ~1e-14 cancellation
            // at an exactly matching CCD/EMAT pose. Both versions receive the
            // identical lowered matrix; the target CCD energies are unchanged.
            for (int pos = 0; pos < min.getNumPos(); pos++) {
                for (int rc = 0; rc < min.getNumConfAtPos(pos); rc++) {
                    min.setOneBody(pos, rc, min.getOneBody(pos, rc) - 1e-8);
                }
            }
            var rig = new SimplerEnergyMatrixCalculator.Builder(new ConfEnergyCalculator(calculator, rigid))
                    .build().calcEnergyMatrix();
            // Warm both code versions before the measured invocation. Each run owns fresh caches.
            for (int run = 0; run < 2; run++) {
                Path directory = out.resolve(run == 0 ? "warmup" : "measured");
                System.setProperty("packstar.pac.frequencySeverity.outputDir", directory.toString());
                System.gc(); // outside timing, discard the previous run's transient heap
                for (var pool : ManagementFactory.getMemoryPoolMXBeans()) pool.resetPeakUsage();
                long start = System.nanoTime();
                PackStarResult result;
                try (var contexts = energy.tasks.contextGroup();
                     var pf = new PackStarPartitionFunction(space, rig, min, calculator,
                             new RCs(space), parallelism, "performance-" + positions)) {
                    pf.init(0.683);
                    pf.putTaskContexts(contexts);
                    pf.compute(Integer.MAX_VALUE);
                    result = pf.makeResult();
                }
                long nanos = System.nanoTime() - start;
                long peakHeap = ManagementFactory.getMemoryPoolMXBeans().stream()
                        .filter(pool -> pool.getType() == java.lang.management.MemoryType.HEAP)
                        .mapToLong(pool -> pool.getPeakUsage().getUsed()).sum();
                if (run == 1) {
                    String header = "mode\tpositions\tseed\tnanos\tpeakHeapBytes\tstatus\tlower\tupper\tepsilon\n";
                    String epsilon;
                    try { epsilon = Double.toString(result.values.getEffectiveEpsilon()); }
                    catch (RuntimeException ex) { epsilon = "unavailable"; }
                    String row = mode + "\t" + positions + "\t" + seed + "\t" + nanos + "\t" + peakHeap
                            + "\t" + result.status + "\t" + result.values.calcLowerBound()
                            + "\t" + result.values.calcUpperBound() + "\t" + epsilon + "\n";
                    Files.writeString(out.resolve("result.tsv"), header + row);
                    System.out.print("BENCH_RESULT " + row);
                }
            }
        }
    }
}
