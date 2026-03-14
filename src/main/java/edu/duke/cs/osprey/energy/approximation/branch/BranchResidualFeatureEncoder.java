package edu.duke.cs.osprey.energy.approximation.branch;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeNode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Hash-based branch/edge local feature encoder for full conformations.
 *
 * Features are designed to be additive over local edge factors:
 * psi(conf) = sum_e phi_e(conf)
 */
public class BranchResidualFeatureEncoder {

    public static class SparseFeatures {
        public final int[] indices;
        public final double[] values;

        public SparseFeatures(int[] indices, double[] values) {
            this.indices = indices;
            this.values = values;
        }

        public int size() {
            return indices.length;
        }
    }

    private final SimpleConfSpace confSpace;
    private final int hashDim;
    private final List<int[]> edges;
    private final int[] posDegrees;
    private final int[] separatorCounts;
    private final long hashSeed;

    public BranchResidualFeatureEncoder(
            SimpleConfSpace confSpace,
            InteractionGraph interactionGraph,
            RootedTreeEdge rootEdge,
            int hashDim
    ) {
        if (hashDim <= 0) {
            throw new IllegalArgumentException("hashDim must be > 0");
        }
        this.confSpace = confSpace;
        this.hashDim = hashDim;
        this.hashSeed = 0x9e3779b97f4a7c15L ^ ((long) hashDim * 0x100000001b3L);

        this.edges = new ArrayList<>();
        if (interactionGraph != null) {
            for (int[] edge : interactionGraph.getEdgeList()) {
                this.edges.add(new int[] { edge[0], edge[1] });
            }
        }

        this.posDegrees = new int[confSpace.positions.size()];
        for (int[] edge : edges) {
            posDegrees[edge[0]]++;
            posDegrees[edge[1]]++;
        }

        this.separatorCounts = collectSeparatorCounts(confSpace.positions.size(), rootEdge);
    }

    public int getHashDim() {
        return hashDim;
    }

    public int getNumEdges() {
        return edges.size();
    }

    public SparseFeatures encode(int[] conf) {
        if (conf == null || conf.length != confSpace.positions.size()) {
            throw new IllegalArgumentException(String.format(
                    "expected full conf assignment of length %d, got %d",
                    confSpace.positions.size(),
                    conf == null ? -1 : conf.length
            ));
        }

        Map<Integer, Double> accum = new HashMap<>(Math.max(16, conf.length + edges.size() * 2));

        add(accum, hashFeature(1, conf.length), 1.0); // bias

        for (int pos = 0; pos < conf.length; pos++) {
            int rc = conf[pos];
            int sepBin = bucket(separatorCounts[pos], 7);
            int degBin = bucket(posDegrees[pos], 7);
            add(accum, hashFeature(2, pos, rc), 1.0);
            add(accum, hashFeature(3, pos, rc, sepBin, degBin), 1.0);
        }

        for (int ei = 0; ei < edges.size(); ei++) {
            int[] edge = edges.get(ei);
            int p = edge[0];
            int q = edge[1];
            int rp = conf[p];
            int rq = conf[q];

            int loPos = Math.min(p, q);
            int hiPos = Math.max(p, q);
            int loRc = (loPos == p) ? rp : rq;
            int hiRc = (loPos == p) ? rq : rp;

            int sepBin = bucket(separatorCounts[p] + separatorCounts[q], 7);
            int degBin = bucket(posDegrees[p] + posDegrees[q], 7);

            add(accum, hashFeature(4, loPos, hiPos, loRc, hiRc), 1.0);
            add(accum, hashFeature(5, loPos, hiPos, sepBin, degBin), 0.25);
            add(accum, hashFeature(6, ei, loRc, hiRc, sepBin, degBin), 0.5);
        }

        List<Map.Entry<Integer, Double>> entries = new ArrayList<>(accum.entrySet());
        entries.sort(Comparator.comparingInt(Map.Entry::getKey));

        int[] indices = new int[entries.size()];
        double[] values = new double[entries.size()];
        for (int i = 0; i < entries.size(); i++) {
            indices[i] = entries.get(i).getKey();
            values[i] = entries.get(i).getValue();
        }

        return new SparseFeatures(indices, values);
    }

    private static void add(Map<Integer, Double> accum, int index, double value) {
        accum.merge(index, value, Double::sum);
    }

    private int hashFeature(int type, int... parts) {
        long h = hashSeed;
        h = mix64(h ^ (type * 0x9e3779b97f4a7c15L));
        for (int part : parts) {
            h = mix64(h ^ ((long) part * 0xbf58476d1ce4e5b9L));
        }
        return (int) Long.remainderUnsigned(h, hashDim);
    }

    private static long mix64(long x) {
        x ^= (x >>> 33);
        x *= 0xff51afd7ed558ccdL;
        x ^= (x >>> 33);
        x *= 0xc4ceb9fe1a85ec53L;
        x ^= (x >>> 33);
        return x;
    }

    private static int bucket(int v, int maxBucket) {
        if (v <= 0) {
            return 0;
        }
        return Math.min(maxBucket, v);
    }

    private static int[] collectSeparatorCounts(int numPositions, RootedTreeEdge rootEdge) {
        int[] counts = new int[numPositions];
        if (rootEdge == null) {
            return counts;
        }

        Set<RootedTreeEdge> seen = Collections.newSetFromMap(new java.util.IdentityHashMap<>());
        collectSeparatorCountsRec(rootEdge, counts, seen);
        return counts;
    }

    private static void collectSeparatorCountsRec(RootedTreeEdge edge, int[] counts, Set<RootedTreeEdge> seen) {
        if (edge == null || !seen.add(edge)) {
            return;
        }

        if (edge.getLambda() != null) {
            for (Integer pos : edge.getLambda()) {
                if (pos != null && pos >= 0 && pos < counts.length) {
                    counts[pos]++;
                }
            }
        }

        RootedTreeNode child = edge.getChild();
        if (child == null) {
            return;
        }
        if (child.getLeftChild() != null) {
            collectSeparatorCountsRec(child.getLeftChild().getChildOfEdge(), counts, seen);
        }
        if (child.getRightChild() != null) {
            collectSeparatorCountsRec(child.getRightChild().getChildOfEdge(), counts, seen);
        }
    }
}
