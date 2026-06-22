package edu.duke.cs.osprey.branchdp;

import java.util.Arrays;

/**
 * Exact vertex-elimination treewidth for small interaction graphs.
 *
 * <p>This is intended as a diagnostic for branch-decomposition experiments.
 * It does not build a branch decomposition, but it gives a rigorous lower bound
 * on branchwidth through {@code tw + 1 <= floor(3*bw/2)}.</p>
 */
public final class ExactTreewidth {

    private ExactTreewidth() {
    }

    public static final class Result {
        public final int treewidth;
        public final int branchwidthLowerBound;
        public final int branchwidthUpperBound;
        public final int numPositions;
        public final int numEdges;
        public final long elapsedNanos;

        private Result(int treewidth, int branchwidthLowerBound, int branchwidthUpperBound,
                       int numPositions, int numEdges, long elapsedNanos) {
            this.treewidth = treewidth;
            this.branchwidthLowerBound = branchwidthLowerBound;
            this.branchwidthUpperBound = branchwidthUpperBound;
            this.numPositions = numPositions;
            this.numEdges = numEdges;
            this.elapsedNanos = elapsedNanos;
        }

        public double elapsedMillis() {
            return elapsedNanos / 1_000_000.0;
        }
    }

    public static Result compute(InteractionGraph graph) {
        int n = graph.getNumPositions();
        if (n > 30) {
            throw new IllegalArgumentException("Exact treewidth diagnostic supports at most 30 positions, got " + n);
        }

        long started = System.nanoTime();
        int treewidth = computeTreewidth(graph);
        int bwLower = graph.getNumEdges() == 0 ? 0 : minBranchwidthFromTreewidth(treewidth);
        int bwUpper = graph.getNumEdges() == 0 ? 0 : treewidth + 1;
        long elapsed = System.nanoTime() - started;

        return new Result(treewidth, bwLower, bwUpper,
                graph.getNumPositions(), graph.getNumEdges(), elapsed);
    }

    private static int computeTreewidth(InteractionGraph graph) {
        int n = graph.getNumPositions();
        if (n <= 1) {
            return 0;
        }

        int[] adjacency = buildAdjacencyMasks(graph);
        int numMasks = 1 << n;
        int fullMask = numMasks - 1;
        int[] dp = new int[numMasks];
        Arrays.fill(dp, n);
        dp[0] = 0;

        for (int eliminatedMask = 0; eliminatedMask < numMasks; eliminatedMask++) {
            int prefixWidth = dp[eliminatedMask];
            if (prefixWidth >= n) {
                continue;
            }
            int remaining = fullMask ^ eliminatedMask;
            while (remaining != 0) {
                int vertexBit = remaining & -remaining;
                int vertex = Integer.numberOfTrailingZeros(vertexBit);
                remaining ^= vertexBit;

                int width = reachableRemainingNeighborCount(adjacency, eliminatedMask, vertex);
                int nextMask = eliminatedMask | vertexBit;
                int nextWidth = Math.max(prefixWidth, width);
                if (nextWidth < dp[nextMask]) {
                    dp[nextMask] = nextWidth;
                }
            }
        }

        return dp[fullMask];
    }

    private static int[] buildAdjacencyMasks(InteractionGraph graph) {
        int n = graph.getNumPositions();
        int[] adjacency = new int[n];
        for (int[] edge : graph.getEdgeList()) {
            int a = edge[0];
            int b = edge[1];
            adjacency[a] |= 1 << b;
            adjacency[b] |= 1 << a;
        }
        return adjacency;
    }

    /**
     * Degree of {@code vertex} in the filled graph after eliminating
     * {@code eliminatedMask}. A remaining vertex is adjacent if it is reachable
     * by a path whose internal vertices are all already eliminated.
     */
    private static int reachableRemainingNeighborCount(int[] adjacency, int eliminatedMask, int vertex) {
        int vertexBit = 1 << vertex;
        int fullMask = (1 << adjacency.length) - 1;
        int remainingMask = fullMask ^ eliminatedMask;
        int visited = vertexBit;
        int queue = adjacency[vertex];
        int reachableRemaining = 0;

        while (queue != 0) {
            int bit = queue & -queue;
            queue ^= bit;
            if ((visited & bit) != 0) {
                continue;
            }
            visited |= bit;

            if ((eliminatedMask & bit) != 0) {
                int eliminatedVertex = Integer.numberOfTrailingZeros(bit);
                queue |= adjacency[eliminatedVertex] & ~visited;
            } else if ((remainingMask & bit) != 0) {
                reachableRemaining |= bit;
            }
        }

        reachableRemaining &= ~vertexBit;
        return Integer.bitCount(reachableRemaining);
    }

    private static int minBranchwidthFromTreewidth(int treewidth) {
        if (treewidth <= 1) {
            return 1;
        }
        int twPlusOne = treewidth + 1;
        return (2 * twPlusOne + 2) / 3;
    }
}
