package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDecomposition;
import edu.duke.cs.osprey.branchdp.BranchDpBackend;
import edu.duke.cs.osprey.branchdp.BranchDpConfig;
import edu.duke.cs.osprey.branchdp.DPTableTooLargeException;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.branchdp.RootedTreeEdge;
import edu.duke.cs.osprey.branchdp.RootedTreeNode;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Structural previews only: no energy evaluation, enumeration arrays or DP tables. */
final class PackStarTripleDecompositionCosts {
    static final String PREFIX = "packstar.pac.frequencySeverity.tripleEta";
    static final String STRATEGY_PROPERTY = PREFIX + "SelectionStrategy";

    static boolean configuredEnabled() {
        String strategy = configuredStrategy();
        if (strategy.equals("budget-forward") || strategy.equals("decomposition-cost")) return true;
        if (strategy.equals("fill-edge")) return false;
        throw new IllegalArgumentException("Unknown " + STRATEGY_PROPERTY + ": " + strategy);
    }

    static String configuredStrategy() {
        return PackStarConfig.getProperty(STRATEGY_PROPERTY, "budget-forward").trim();
    }

    /** A deterministic structural tie-break, never a wall-time estimate. */
    static int compare(Cost a, Cost b) {
        int c = a.work.compareTo(b.work);
        if (c == 0) c = Long.compare(a.hostBytes, b.hostBytes);
        if (c == 0) c = Long.compare(a.fileBytes, b.fileBytes);
        if (c == 0) c = Long.compare(a.totalTableBytes, b.totalTableBytes);
        if (c == 0) c = Long.compare(a.maxTableBytes, b.maxTableBytes);
        if (c == 0) c = Long.compare(a.maxMStates, b.maxMStates);
        return c;
    }

    @FunctionalInterface
    interface Previewer {
        Cost preview(Set<Long> fillEdges);
    }

    static final class Cost {
        final int branchwidth;
        final int rootSplitIndex;
        final long maxMStates, totalMStates, maxTableBytes, totalTableBytes;
        final long hostBytes, fileBytes;
        final BigInteger work;

        Cost(int branchwidth, long maxMStates, long totalMStates,
             long maxTableBytes, long totalTableBytes, long hostBytes, long fileBytes,
             BigInteger work) {
            this(branchwidth, maxMStates, totalMStates, maxTableBytes, totalTableBytes,
                    hostBytes, fileBytes, work, -1);
        }

        Cost(int branchwidth, long maxMStates, long totalMStates,
             long maxTableBytes, long totalTableBytes, long hostBytes, long fileBytes,
             BigInteger work, int rootSplitIndex) {
            this.branchwidth = branchwidth;
            this.rootSplitIndex = rootSplitIndex;
            this.maxMStates = maxMStates;
            this.totalMStates = totalMStates;
            this.maxTableBytes = maxTableBytes;
            this.totalTableBytes = totalTableBytes;
            this.hostBytes = hostBytes;
            this.fileBytes = fileBytes;
            this.work = work;
        }

        /** Dimensionless growth; a zero-cost candidate never gets infinite priority. */
        double growthFrom(Cost before) {
            double growth = Math.max(0, work.doubleValue() / Math.max(1, before.work.doubleValue()) - 1);
            growth = Math.max(growth, growth(maxMStates, before.maxMStates));
            growth = Math.max(growth, growth(totalMStates, before.totalMStates));
            growth = Math.max(growth, growth(maxTableBytes, before.maxTableBytes));
            growth = Math.max(growth, growth(totalTableBytes, before.totalTableBytes));
            growth = Math.max(growth, growth(hostBytes, before.hostBytes));
            // Crossing into mmap is normalized by the previous total storage,
            // not one byte when the previous file footprint was zero.
            return Math.max(growth, Math.max(0, (double) fileBytes - before.fileBytes)
                    / Math.max(1, before.totalTableBytes));
        }

        private static double growth(long after, long before) {
            return Math.max(0, (double) after - before) / Math.max(1, before);
        }

        @Override public String toString() {
            return "bw=" + branchwidth + ";rootSplit=" + rootSplitIndex + ";maxMStates=" + maxMStates
                    + ";totalMStates=" + totalMStates + ";maxTableBytes=" + maxTableBytes
                    + ";totalTableBytes=" + totalTableBytes + ";hostBytes=" + hostBytes
                    + ";fileBytes=" + fileBytes + ";work=" + work;
        }
    }

    static final class Limits {
        final long maxMStates, maxTableBytes, totalTableBytes, hostBytes, fileBytes;
        final BigInteger maxWork;

        Limits(long maxMStates, long maxTableBytes, long totalTableBytes,
               long hostBytes, long fileBytes) {
            this(maxMStates, maxTableBytes, totalTableBytes, hostBytes, fileBytes,
                    BigInteger.valueOf(Long.MAX_VALUE));
        }

        Limits(long maxMStates, long maxTableBytes, long totalTableBytes,
               long hostBytes, long fileBytes, BigInteger maxWork) {
            if (maxWork == null || maxWork.signum() < 1)
                throw new IllegalArgumentException("triple MaxWork must be positive");
            this.maxWork = maxWork;
            if (maxMStates < 1 || maxTableBytes < 1 || totalTableBytes < 1
                    || hostBytes < 1 || fileBytes < 0) {
                throw new IllegalArgumentException("Invalid triple decomposition resource limits");
            }
            this.maxMStates = maxMStates;
            this.maxTableBytes = maxTableBytes;
            this.totalTableBytes = totalTableBytes;
            this.hostBytes = hostBytes;
            this.fileBytes = fileBytes;
        }

        static Limits configured() {
            long defaultHost = PackStarConfig.getBytes("packstar.rootSplit.hostBudgetBytes",
                    0, "PACK*");
            if (defaultHost <= 0) defaultHost = Math.max(1, (long) (Runtime.getRuntime().maxMemory() * 0.7));
            return new Limits(
                    PackStarConfig.getLong(PREFIX + "MaxMStates", Long.MAX_VALUE, "PACK*"),
                    PackStarConfig.getBytes(PREFIX + "MaxTableBytes", Long.MAX_VALUE, "PACK*"),
                    PackStarConfig.getBytes(PREFIX + "MaxTotalTableBytes", Long.MAX_VALUE, "PACK*"),
                    PackStarConfig.getBytes(PREFIX + "MaxHostBytes", defaultHost, "PACK*"),
                    PackStarConfig.getBytes(PREFIX + "MaxFileBytes", Long.MAX_VALUE, "PACK*"),
                    new BigInteger(PackStarConfig.getProperty(PREFIX + "MaxWork",
                            Long.toString(Long.MAX_VALUE)).trim()));
        }

        String rejection(Cost cost) {
            if (cost.maxMStates == Long.MAX_VALUE || cost.totalMStates == Long.MAX_VALUE
                    || cost.totalTableBytes == Long.MAX_VALUE || cost.hostBytes == Long.MAX_VALUE
                    || cost.fileBytes == Long.MAX_VALUE) return "size-overflow";
            if (cost.work.compareTo(maxWork) > 0) return "work";
            if (cost.maxMStates > maxMStates) return "max-M-states";
            if (cost.maxTableBytes > maxTableBytes) return "max-table-bytes";
            if (cost.totalTableBytes > totalTableBytes) return "total-table-bytes";
            if (cost.hostBytes > hostBytes) return "host-bytes";
            if (cost.fileBytes > fileBytes) return "file-bytes";
            return null;
        }

        @Override public String toString() {
            return "maxMStates=" + maxMStates + ";maxTableBytes=" + maxTableBytes
                    + ";totalTableBytes=" + totalTableBytes + ";hostBytes=" + hostBytes
                    + ";fileBytes=" + fileBytes + ";maxWork=" + maxWork;
        }
    }

    static final class ProposalRoot {
        final BranchDpBackend.RootingCandidate selected;
        final int branchwidth;
        final long rootSelectionNanos;
        final boolean materialized;

        ProposalRoot(BranchDpBackend.RootingCandidate selected, int branchwidth, long nanos,
                boolean materialized) {
            this.selected = selected;
            this.branchwidth = branchwidth;
            this.rootSelectionNanos = nanos;
            this.materialized = materialized;
        }

        Cost cost() {
            Cost c = summarize(selected.root, !materialized);
            return new Cost(branchwidth, c.maxMStates, c.totalMStates, c.maxTableBytes,
                    c.totalTableBytes, c.hostBytes, c.fileBytes, c.work, selected.splitEdgeIndex);
        }
    }

    /** The same graph, decomposition and configured root policy for preview and execution. */
    static ProposalRoot rootProposalGraph(InteractionGraph graph, RCs rcs,
            SimpleConfSpace confSpace, boolean materialize) {
        try (BranchDpConfig.Scope ignored = BranchDpConfig.enterPackStarAliasScope()) {
            int[] counts = new int[rcs.getNumPos()];
            for (int p = 0; p < counts.length; p++) counts[p] = rcs.getNum(p);
            BranchDecomposition decomposition = new BranchDecomposition(graph,
                    BranchDecomposition.Strategy.WEIGHTED_HICKS, counts, null, false);
            decomposition.compute();
            long started = System.nanoTime();
            BranchDpBackend.RootingCandidate selected = BranchDpBackend.selectConfiguredRoot(
                    decomposition, graph, rcs, confSpace, materialize);
            if (selected == null) throw new IllegalStateException("Proposal graph has no branch tree");
            return new ProposalRoot(selected, decomposition.getBranchwidth(),
                    System.nanoTime() - started, materialize);
        }
    }

    /** One cache per immutable pfunc graph/RC domain, shared by all folds and refits. */
    static final class Cache implements Previewer {
        private final RCs rcs;
        private final InteractionGraph graph;
        private final RootedTreeNode baseRoot;
        private final SimpleConfSpace confSpace;
        private final int capacity;
        private final Map<List<Long>, Cost> costs = new LinkedHashMap<>(16, 0.75f, true);
        private Cost baseCost;
        long requests, hits, builds, nanos;
        long rootSelectionNanos;

        Cache(RCs rcs, InteractionGraph graph, RootedTreeNode baseRoot, int capacity) {
            this(rcs, graph, baseRoot, capacity, null);
        }

        Cache(RCs rcs, InteractionGraph graph, RootedTreeNode baseRoot, int capacity,
                SimpleConfSpace confSpace) {
            if (capacity < 1) throw new IllegalArgumentException("preview cache capacity must be positive");
            this.rcs = rcs;
            this.graph = graph;
            this.baseRoot = baseRoot;
            this.confSpace = confSpace;
            this.capacity = capacity;
        }

        @Override public Cost preview(Set<Long> fillEdges) {
            requests++;
            List<Long> key = new ArrayList<>(fillEdges);
            Collections.sort(key);
            if (key.isEmpty() && baseCost != null) { hits++; return baseCost; }
            Cost cached = costs.get(key);
            if (cached != null) { hits++; return cached; }
            long started = System.nanoTime();
            Cost result;
            // Match PACK* aliases even when a caller is outside a backend scope.
            try (BranchDpConfig.Scope ignored = BranchDpConfig.enterPackStarAliasScope()) {
                if (key.isEmpty() && baseRoot != null) {
                    result = summarize(baseRoot, false);
                } else {
                    List<int[]> edges = new ArrayList<>(graph.getEdgeList());
                    for (long packed : key) edges.add(new int[]{(int) (packed >>> 32), (int) packed});
                    InteractionGraph desired = InteractionGraph.buildFromEdges(graph.getNumPositions(), edges);
                    RootedTreeNode root = null;
                    try {
                        ProposalRoot proposal = rootProposalGraph(desired, rcs, confSpace, false);
                        root = proposal.selected.root;
                        rootSelectionNanos += proposal.rootSelectionNanos;
                        result = proposal.cost();
                    } catch (DPTableTooLargeException | BranchDpBackend.RootSelectionException ex) {
                        // Capacity failures are an optional-candidate rejection;
                        // cache them too. Programming/configuration errors propagate.
                        result = new Cost(Integer.MAX_VALUE, Long.MAX_VALUE,
                                Long.MAX_VALUE, Long.MAX_VALUE, Long.MAX_VALUE,
                                Long.MAX_VALUE, Long.MAX_VALUE, BigInteger.valueOf(Long.MAX_VALUE));
                        System.out.println("[PACK*-triple-m2-preview] rejected fill=" + key
                                + ", reason=" + ex.getMessage());
                    } finally {
                        if (root != null) RootedTreeEdge.postOrderReleaseLargeMemory(root);
                    }
                }
                builds++;
            } finally {
                nanos += System.nanoTime() - started;
            }
            if (key.isEmpty()) baseCost = result;
            else {
                costs.put(List.copyOf(key), result);
                if (costs.size() > capacity) costs.remove(costs.keySet().iterator().next());
            }
            return result;
        }

        String audit() {
            return "requests=" + requests + ", hits=" + hits + ", builds=" + builds
                    + ", cachedGraphs=" + costs.size() + ", previewMs=" + nanos / 1e6
                    + ", rootSelectionMs=" + rootSelectionNanos / 1e6;
        }
    }

    static Cost summarize(RootedTreeNode root, boolean requireNoTables) {
        List<RootedTreeEdge> edges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(root, edges);
        long maxM = 0, totalM = 0, maxTable = 0, totalTable = 0, host = 0, file = 0;
        int width = maximumWidth(root);
        BigInteger work = BigInteger.ZERO;
        for (RootedTreeEdge edge : edges) {
            if (requireNoTables && edge.hasDPTable())
                throw new IllegalStateException("triple preview allocated a DP table");
            long m = edge.getMStateCount(), lambda = edge.getTotalLambdaStates();
            if (m < 0 || lambda < 0) throw new IllegalStateException("negative preview state count");
            long table = BranchDpBackend.estimateDPTableBytes(m);
            boolean mapped = RootedTreeEdge.shouldUseFileBackedDPTable(m);
            host = add(host, mapped ? BranchDpBackend.estimateDPAuxHostBytes(m, lambda)
                    : BranchDpBackend.estimateDPHostBytes(m, lambda));
            if (mapped) file = add(file, table);
            maxM = Math.max(maxM, m); totalM = add(totalM, m);
            maxTable = Math.max(maxTable, table); totalTable = add(totalTable, table);
            work = work.add(BigInteger.valueOf(Math.max(1, m)).multiply(BigInteger.valueOf(Math.max(1, lambda))));
        }
        return new Cost(width, maxM, totalM, maxTable, totalTable, host, file, work);
    }

    private static long add(long a, long b) {
        return a > Long.MAX_VALUE - b ? Long.MAX_VALUE : a + b;
    }

    private static int maximumWidth(RootedTreeNode node) {
        if (node == null) return 0;
        int width = node.getChildOfEdge() == null ? 0 : node.getChildOfEdge().getM().size();
        return Math.max(width, Math.max(maximumWidth(node.getLeftChild()), maximumWidth(node.getRightChild())));
    }
}
