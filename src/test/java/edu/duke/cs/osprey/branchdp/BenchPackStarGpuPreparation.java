package edu.duke.cs.osprey.branchdp;

import java.lang.reflect.*;
import java.nio.file.*;
import java.util.*;

/** Same seeded requests on the old and new classes; no CUDA kernel changes. */
public final class BenchPackStarGpuPreparation {
    public static void main(String[] args) throws Exception {
        Path out = Path.of(args[0]);
        boolean optimized = args[1].equals("optimized");
        Files.createDirectories(out);
        System.setProperty("branchdp.dp.foldChildren", "true");
        System.setProperty("branchdp.pac.sampling.gpu.persistentContext", "true");
        System.setProperty("branchdp.pac.sampling.gpu.residentChildTables", "true");
        Method buildParent = TestGpuFullDP.class.getDeclaredMethod("buildParent",
                int[].class, int[].class, int[].class, int[][].class);
        buildParent.setAccessible(true);
        RootedTreeEdge edge = (RootedTreeEdge) buildParent.invoke(new TestGpuFullDP(),
                new int[]{256, 128, 32, 4}, new int[]{0, 1}, new int[]{2, 3}, new int[][]{{0, 1, 2}});
        if (!edge.canUseGpuSampling()) throw new AssertionError("GPU sampling structural gate failed");
        Method prepare;
        Object cache = null;
        if (optimized) {
            Class<?> cacheClass = Class.forName(RootedTreeEdge.class.getName() + "$SamplingRequestCache");
            cache = cacheClass.getConstructor(long.class).newInstance(64L << 20);
            prepare = RootedTreeEdge.class.getDeclaredMethod("prepareGpuSamplingRequest",
                    long[].class, long.class, boolean.class, cacheClass);
        } else {
            prepare = RootedTreeEdge.class.getDeclaredMethod("buildGpuSamplingRequest",
                    long[].class, long.class, boolean.class);
        }
        prepare.setAccessible(true);
        long[] assignments = new long[512];
        for (int i = 0; i < assignments.length; i++) assignments[i] = i % 256;
        StringBuilder timings = new StringBuilder("mode\tmethod\trepeat\tprepareNanos\tsampleNanos\n");
        for (String method : List.of("gumbel", "cdfDedup")) {
            SamplingGpuPhase1.resetForTesting();
            if (cache != null) cache.getClass().getMethod("clear").invoke(cache);
            System.setProperty("branchdp.pac.sampling.gpu.method", method);
            for (int repeat = -2; repeat < 10; repeat++) {
                long start = System.nanoTime();
                SamplingGpuPhase1.Request request = (SamplingGpuPhase1.Request) (optimized
                        ? prepare.invoke(edge, assignments, 123L + repeat, false, cache)
                        : prepare.invoke(edge, assignments, 123L + repeat, false));
                long prepared = System.nanoTime();
                request.multiGpu = false;
                request.maxGpus = 1;
                int[] samples = SamplingGpuPhase1.sample(request);
                long finished = System.nanoTime();
                if (samples == null) throw new AssertionError("GPU sampling unexpectedly fell back");
                if (repeat >= 0) {
                    timings.append(args[1]).append('\t').append(method).append('\t').append(repeat)
                            .append('\t').append(prepared - start).append('\t').append(finished - prepared).append('\n');
                    Files.writeString(out.resolve(method + "-" + repeat + ".txt"), Arrays.toString(samples));
                }
            }
        }
        if (cache != null) cache.getClass().getMethod("clear").invoke(cache);
        SamplingGpuPhase1.resetForTesting();
        for (RootedTreeEdge child : edge.getFset()) child.releaseLargeMemory();
        edge.releaseLargeMemory();
        Files.writeString(out.resolve("timings.tsv"), timings);
        System.out.print(timings);
    }
}
