/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
**
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

import java.util.*;

/**
 * Per-confspace DP evaluator for branch decomposition Z bounds.
 *
 * Built once per confspace (protein / ligand / complex). The interaction graph and
 * branch tree structure depend only on the confspace geometry, not on the sequence.
 * Different sequences produce different RCs (filtering rotamers by AA type), but
 * the tree topology is the same.
 *
 * Usage:
 *   BranchDPScorer scorer = new BranchDPScorer(confSpace, rigidEmat, minEmat, 8.0, 0.1);
 *   double[] bounds = scorer.computeZBounds(rcs);  // [logZLower, logZUpper]
 *
 * Cross-sequence caching (改动 3):
 *   Edges whose M-set and λ-set contain NO design positions produce the same DP table
 *   for all sequences → cached permanently after first computation.
 *   Edges that DO touch design positions are cached by AA-type key at those positions.
 */
public class BranchDPScorer {

    private final SimpleConfSpace confSpace;
    private final EnergyMatrix rigidEmat;
    private final EnergyMatrix minimizingEmat;
    private final InteractionGraph interactionGraph;
    private final BranchDecomposition branchDecomposition;
    private final double RT;

    // Branch tree structure (topology, built once)
    private final RootedTreeNode rootedRoot;
    private final RootedTreeEdge rootedRootEdge;
    private final List<RootedTreeEdge> allLambdaEdges;

    // Design positions: positions where different sequences have different AA types
    private final Set<Integer> designPositions;

    // Cross-sequence DP cache: edgeIndex → (cacheKey → [logZLower[], logZUpper[]])
    // cacheKey = AA types at design positions within this edge's M+λ set
    private final Map<Integer, Map<String, double[][]>> dpCache;

    // Edge indexing for cache lookup
    private final Map<RootedTreeEdge, Integer> edgeIndices;

    // Whether this scorer is usable (branch decomposition succeeded)
    private final boolean usable;

    // Cutoff parameters (stored for logging)
    private final double distCutoff;
    private final double energyCutoff;

    private static final double TESS_FALLBACK_THRESHOLD = 0.5;

    /**
     * Build a BranchDPScorer for a given confspace.
     * Constructs the interaction graph and branch decomposition once.
     * Returns a scorer that can compute Z bounds for any sequence's RCs.
     *
     * @param confSpace      the conformation space
     * @param rigidEmat      rigid energy matrix (for Z lower bound via energy upper bound)
     * @param minimizingEmat minimized energy matrix (for Z upper bound via energy lower bound)
     * @param distCutoff     distance cutoff for interaction graph edges (Angstroms)
     * @param energyCutoff   energy cutoff for interaction graph edges (kcal/mol)
     */
    public BranchDPScorer(SimpleConfSpace confSpace,
                          EnergyMatrix rigidEmat,
                          EnergyMatrix minimizingEmat,
                          double distCutoff, double energyCutoff) {
        this.confSpace = confSpace;
        this.rigidEmat = rigidEmat;
        this.minimizingEmat = minimizingEmat;
        this.RT = BoltzmannCalculator.RClassic * 298.0;
        this.distCutoff = distCutoff;
        this.energyCutoff = energyCutoff;
        this.dpCache = new HashMap<>();
        this.edgeIndices = new HashMap<>();

        // Identify design positions (positions with mutations)
        designPositions = new HashSet<>();
        for (int i = 0; i < confSpace.positions.size(); i++) {
            if (confSpace.positions.get(i).resTypes.size() > 1) {
                designPositions.add(i);
            }
        }

        // Build interaction graph with full RCs (all rotamers, unfiltered)
        RCs fullRCs = new RCs(confSpace);
        interactionGraph = InteractionGraph.buildWithDualCutoff(
                confSpace, rigidEmat, minimizingEmat, fullRCs,
                distCutoff, energyCutoff);

        // Small graph guard: branch decomposition requires enough edges to form a proper binary tree.
        // With ≤ 3 interaction edges the tree may have degree-2 internal nodes, which breaks rooting.
        int numIGEdges = interactionGraph.getNumEdges();
        if (numIGEdges <= 3) {
            System.out.println("BranchDPScorer: Interaction graph too small (" + numIGEdges
                    + " edges), scorer not usable. Falling back to standard scoring.");
            this.branchDecomposition = null;
            this.rootedRoot = null;
            this.rootedRootEdge = null;
            this.allLambdaEdges = Collections.emptyList();
            this.usable = false;
            return;
        }

        // Compute branch decomposition
        branchDecomposition = new BranchDecomposition(interactionGraph);
        branchDecomposition.compute();

        // Root the tree with full RCs (topology only — actual DP uses per-sequence RCs)
        RootedTreeNode tempRoot = branchDecomposition.rootBranchTree(fullRCs);
        if (tempRoot == null) {
            System.out.println("BranchDPScorer: Empty tree, scorer not usable.");
            this.rootedRoot = null;
            this.rootedRootEdge = null;
            this.allLambdaEdges = Collections.emptyList();
            this.usable = false;
            return;
        }

        // Compute L/lambda/F-sets (topology depends on tree structure, not on RCs values)
        RootedTreeEdge.postOrderCompLlambda(tempRoot);
        RootedTreeEdge tempRootEdge = tempRoot.getLeftChild().getChildOfEdge();
        tempRootEdge.compactTree();

        // TESS check
        double logTESS = tempRootEdge.computeLogTESS();
        double logNaive = 0.0;
        for (int pos = 0; pos < fullRCs.getNumPos(); pos++) {
            if (fullRCs.getNum(pos) > 0) {
                logNaive += Math.log(fullRCs.getNum(pos));
            }
        }
        double tessRatio = Math.exp(logTESS - logNaive);

        if (tessRatio > TESS_FALLBACK_THRESHOLD) {
            System.out.println("BranchDPScorer: TESS ratio " + String.format("%.4f", tessRatio)
                    + " > threshold, scorer not usable.");
            this.rootedRoot = null;
            this.rootedRootEdge = null;
            this.allLambdaEdges = Collections.emptyList();
            this.usable = false;
            return;
        }

        this.rootedRoot = tempRoot;
        this.rootedRootEdge = tempRootEdge;

        // Collect lambda edges and index them
        List<RootedTreeEdge> edges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(rootedRoot, edges);
        this.allLambdaEdges = edges;
        for (int i = 0; i < edges.size(); i++) {
            edgeIndices.put(edges.get(i), i);
        }

        this.usable = true;

        System.out.println("BranchDPScorer: Built. branchwidth=" + branchDecomposition.getBranchwidth()
                + ", " + allLambdaEdges.size() + " lambda-edges"
                + ", " + designPositions.size() + " design positions"
                + ", TESS ratio=" + String.format("%.4f", tessRatio));
    }

    /** Whether this scorer can be used (branch decomposition succeeded and TESS check passed). */
    public boolean isUsable() {
        return usable;
    }

    /**
     * Compute Z bounds for a given set of RCs (sequence-specific).
     *
     * Rebuilds the rooted tree with the given RCs (which may have fewer rotamers
     * per position than the full confspace), then runs one-shot bottom-up DP.
     *
     * Cross-sequence caching: edges whose positions don't overlap with design positions
     * are cached permanently. Others are cached by AA-type key.
     *
     * @param rcs sequence-specific rotamer conformations
     * @return [logZLower, logZUpper] for the root edge (mIdx=0)
     */
    public double[] computeZBounds(RCs rcs) {
        if (!usable) {
            return new double[]{ Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        }

        // Rebuild tree with these RCs (new RCs → new M-state sizes, new energies)
        RootedTreeNode seqRoot = branchDecomposition.rootBranchTree(rcs);
        if (seqRoot == null) {
            return new double[]{ Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        }

        RootedTreeEdge.postOrderCompLlambda(seqRoot);
        RootedTreeEdge seqRootEdge = seqRoot.getLeftChild().getChildOfEdge();
        seqRootEdge.compactTree();

        // Initialize enumeration data (pre-compute energies, allocate arrays)
        RootedTreeEdge.postOrderInitIncremental(seqRoot, rigidEmat, minimizingEmat,
                interactionGraph, RT);

        // Collect lambda edges for this tree
        List<RootedTreeEdge> seqEdges = new ArrayList<>();
        RootedTreeEdge.collectLambdaEdges(seqRoot, seqEdges);

        // One-shot bottom-up DP with caching
        postOrderComputeWithCache(seqRoot, rcs);

        // Return root edge bounds
        double[] rootLower = seqRootEdge.getLogZLower();
        double[] rootUpper = seqRootEdge.getLogZUpper();
        if (rootLower != null && rootLower.length > 0) {
            return new double[]{ rootLower[0], rootUpper[0] };
        }
        return new double[]{ Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
    }

    /**
     * Post-order DP with cross-sequence caching.
     * For each lambda-edge, compute a cache key based on AA types at design positions.
     * If cache hit, copy cached logZ values. Otherwise, compute and cache.
     */
    private void postOrderComputeWithCache(RootedTreeNode node, RCs rcs) {
        if (node == null) return;
        postOrderComputeWithCache(node.getLeftChild(), rcs);
        postOrderComputeWithCache(node.getRightChild(), rcs);

        if (node.getChildOfEdge() != null && node.getChildOfEdge().getIsLambdaEdge()) {
            RootedTreeEdge edge = node.getChildOfEdge();

            // Compute cache key for this edge
            String cacheKey = computeCacheKey(edge, rcs);
            int edgeIdx = getEdgeTopologyIndex(edge);

            // Check cache
            Map<String, double[][]> edgeCache = dpCache.get(edgeIdx);
            if (edgeCache != null) {
                double[][] cached = edgeCache.get(cacheKey);
                if (cached != null && cached[0].length == edge.getMArraySize()) {
                    // Cache hit: copy values
                    double[] logZLower = edge.getLogZLower();
                    double[] logZUpper = edge.getLogZUpper();
                    System.arraycopy(cached[0], 0, logZLower, 0, logZLower.length);
                    System.arraycopy(cached[1], 0, logZUpper, 0, logZUpper.length);
                    return;
                }
            }

            // Cache miss: compute
            edge.computeFullDP();

            // Store in cache
            if (edgeCache == null) {
                edgeCache = new HashMap<>();
                dpCache.put(edgeIdx, edgeCache);
            }
            double[] logZLower = edge.getLogZLower();
            double[] logZUpper = edge.getLogZUpper();
            double[][] toCache = new double[2][];
            toCache[0] = Arrays.copyOf(logZLower, logZLower.length);
            toCache[1] = Arrays.copyOf(logZUpper, logZUpper.length);
            edgeCache.put(cacheKey, toCache);
        }
    }

    /**
     * Compute cache key for an edge based on AA types at design positions
     * within the edge's M-set and λ-set.
     *
     * If no design positions overlap with M ∪ λ, returns empty string (all sequences same).
     * Otherwise, returns "pos=AAtype,pos=AAtype,..." for design positions in sorted order.
     */
    private String computeCacheKey(RootedTreeEdge edge, RCs rcs) {
        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        StringBuilder key = new StringBuilder();
        // Check M positions
        for (int pos : mPositions) {
            if (designPositions.contains(pos)) {
                key.append(pos).append('=').append(getAATypeKey(rcs, pos)).append(',');
            }
        }
        // Check lambda positions
        for (int pos : lambdaPositions) {
            if (designPositions.contains(pos)) {
                key.append(pos).append('=').append(getAATypeKey(rcs, pos)).append(',');
            }
        }
        return key.toString();
    }

    /**
     * Get a string key representing the AA type at a position for the given RCs.
     * Uses the number of RCs and their indices as a proxy for AA type identity.
     */
    private String getAATypeKey(RCs rcs, int pos) {
        int numRCs = rcs.getNum(pos);
        if (numRCs == 0) return "X";
        // Use first RC's global index as AA type identifier
        // (all RCs at a design position for a given sequence have the same AA type)
        return String.valueOf(rcs.get(pos, 0));
    }

    /**
     * Get a topology index for the edge based on its position in the post-order traversal.
     * Since we rebuild the tree for each sequence, we use M+λ position sets as the identity.
     */
    private int getEdgeTopologyIndex(RootedTreeEdge edge) {
        // Use hash of M positions + lambda positions as topology index
        int hash = Arrays.hashCode(edge.getMPositionsSorted()) * 31
                 + Arrays.hashCode(edge.getLambdaPositionsSorted());
        return hash;
    }

    /** Get the number of cached edge entries (for diagnostics). */
    public int getCacheSize() {
        int total = 0;
        for (Map<String, double[][]> edgeCache : dpCache.values()) {
            total += edgeCache.size();
        }
        return total;
    }

    /** Clear the cache. */
    public void clearCache() {
        dpCache.clear();
    }
}
