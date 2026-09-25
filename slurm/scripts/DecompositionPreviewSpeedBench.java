package edu.duke.cs.osprey.branchdp;

import java.util.*;

/** Compare exact ordered trees against the frozen implementation in another JVM. */
public final class DecompositionPreviewSpeedBench {
    public static void main(String[] args) {
        for (int trial = 0; trial < 9; trial++) {
            int n = 6 + trial % 3 * 3;
            Random random = new Random(9371 + trial);
            List<int[]> edges = new ArrayList<>();
            for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++)
                if (q == p + 1 || (p == 0 && q == n - 1) || random.nextDouble() < 0.3)
                    edges.add(new int[]{p, q});
            int[] counts = new int[n];
            for (int p = 0; p < n; p++) counts[p] = trial >= 6 ? 50000 + p : 2 + random.nextInt(30);
            InteractionGraph graph = InteractionGraph.buildFromEdges(n, edges);
            long elapsed = 0;
            String signature = null;
            for (int repeat = 0; repeat < 4; repeat++) {
                BranchDecomposition d = new BranchDecomposition(graph,
                        BranchDecomposition.Strategy.WEIGHTED_HICKS, counts, null, false);
                long start = System.nanoTime();
                d.compute();
                long time = System.nanoTime() - start;
                if (repeat > 0) elapsed += time;
                StringBuilder tree = new StringBuilder().append(d.getBranchwidth());
                BranchTree t = d.getTree();
                for (int e = 0; e < t.getNumEdges(); e++) {
                    BranchEdge edge = t.getEdge(e);
                    tree.append(';').append(edge.getn1().getIndex()).append(',')
                            .append(edge.getn2().getIndex()).append(':').append(edge.getM());
                }
                String current = tree.toString();
                if (signature != null && !signature.equals(current)) throw new AssertionError("non-deterministic tree");
                signature = current;
            }
            System.out.printf(Locale.ROOT, "RESULT\t%d\t%.9f\t%s%n", trial, elapsed / 1e9, signature);
        }
    }
}
