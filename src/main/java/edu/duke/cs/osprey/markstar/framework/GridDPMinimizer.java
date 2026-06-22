package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.branchdp.RootedTreeEdge;
import edu.duke.cs.osprey.branchdp.RootedTreeNode;
import edu.duke.cs.osprey.structure.AtomConnectivity;

import org.apache.commons.math3.linear.*;

import java.util.*;
import java.util.concurrent.*;

/**
 * Grid DP Minimizer: approximates CCD minimization using discrete grid + branch decomposition DP.
 *
 * Given a fixed RC tuple, discretizes each position's continuous DOFs into g grid points,
 * precomputes one-body (intra + shell) and pairwise energy tables, then uses the branch
 * decomposition tree to find the optimal grid combination via bottom-up DP.
 *
 * Result is a valid upper bound on the true minimized energy (true minimum is between grid points).
 */
public class GridDPMinimizer {

    private final SimpleConfSpace confSpace;
    private final InteractionGraph interactionGraph;
    private final RootedTreeEdge rootEdge;
    private final int gridSize;
    private final ForcefieldParams ffparams;
    private final AtomConnectivity atomConnectivity;
    private final ResPairCache resPairCache;
    private final SimpleReferenceEnergies eref;  // null = no reference energy correction
    private final int numThreads;

    // Energy table cache (lazy, across conformations)
    private final boolean useCache;
    private HashMap<Long, double[]> oneBodyCache;      // key: (pos, rc) → double[gridStates]
    private HashMap<Long, double[][]> pairCache;        // key: (pos_i, rc_i, pos_j, rc_j) → double[gs_i][gs_j]
    // Cache stats (per minimize call)
    private int oneBodyHits, oneBodyMisses;
    private int pairHits, pairMisses;
    // Subtree DP cache: key = RC assignments at positions in subtree
    private HashMap<RootedTreeEdge, int[]> edgeSubtreePositions;  // precomputed: edge → positions affecting its DP
    private HashMap<RootedTreeEdge, HashMap<Long, double[]>> dpCache;
    private int dpCacheHits, dpCacheMisses;

    // Per-minimization state
    private int numPositions;
    private ParametricMolecule pmol;
    private int[] confAssignments;  // RC index at each position for current minimization
    private int[] bestGridStatePerPosition;  // best grid state index for each position (after DP)

    // Position -> DOF mapping
    private int[][] positionDOFIndices;   // positionDOFIndices[pos] = {global dof index, ...}
    private int[] numDOFsPerPosition;
    private int[] gridStatesPerPosition;  // g^numDOFs per position

    // Grid point values: gridValues[globalDofIdx][gridPointIdx]
    private double[][] gridValues;

    // Precomputed energy tables
    private double[][] oneBodyEnergy;     // [pos][gridState]
    // Pairwise: stored as flat map from edge key to 2D table
    private HashMap<Long, double[][]> pairEnergyTables;

    // DP tables per edge
    private HashMap<RootedTreeEdge, double[]> dpTables;

    // Root M-state indices sorted by energy (ascending), for top-N access
    private int[] sortedRootMStates;
    private double bestTotalEnergy;

    // Top-N enumeration mode: true = priority queue (all edges), false = brute-force (root edge only)
    private boolean useTopNPriorityQueue = true;

    // Lazy top-N enumeration state (BWM*-style on-demand pop)
    private boolean enumInitialized = false;
    private int enumPopCount = 0;
    private Map<RootedTreeEdge, int[]> enumSortedLambdas;
    private Map<RootedTreeEdge, double[]> enumLambdaContribs;
    private Map<RootedTreeEdge, Double> enumBestEdgeEnergies;
    private PriorityQueue<Object[]> enumPQ;

    // Timing
    public long precomputeTimeNs;
    public long dpTimeNs;

    public GridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                           RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                           SimpleReferenceEnergies eref, boolean useCache, int numThreads) {
        this.confSpace = confSpace;
        this.interactionGraph = interactionGraph;
        this.rootEdge = rootEdge;
        this.gridSize = gridSize;
        this.ffparams = ffparams;
        this.eref = eref;
        this.useCache = useCache;
        this.numThreads = Math.max(1, numThreads);
        this.numPositions = confSpace.positions.size();

        if (useCache) {
            this.oneBodyCache = new HashMap<>();
            this.pairCache = new HashMap<>();
            this.dpCache = new HashMap<>();
            this.edgeSubtreePositions = new HashMap<>();
            precomputeEdgeSubtreePositions(rootEdge);
        }

        this.atomConnectivity = new AtomConnectivity.Builder().build();
        this.resPairCache = new ResPairCache(ffparams, atomConnectivity);
    }

    /** Constructor with default numThreads=1 */
    public GridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                           RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                           SimpleReferenceEnergies eref, boolean useCache) {
        this(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, useCache, 1);
    }

    /** Constructor with default useCache=false (backward compatible) */
    public GridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                           RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams,
                           SimpleReferenceEnergies eref) {
        this(confSpace, interactionGraph, rootEdge, gridSize, ffparams, eref, false, 1);
    }

    /** Convenience constructor without reference energies */
    public GridDPMinimizer(SimpleConfSpace confSpace, InteractionGraph interactionGraph,
                           RootedTreeEdge rootEdge, int gridSize, ForcefieldParams ffparams) {
        this(confSpace, interactionGraph, rootEdge, gridSize, ffparams, null, false, 1);
    }

    public static class Result {
        public final double energy;
        public final long precomputeTimeNs;
        public final long dpTimeNs;
        public final int[] bestGridAssignment;  // grid state index per position
        public final int cacheOneBodyHits, cacheOneBodyMisses;
        public final int cachePairHits, cachePairMisses;
        public final int cacheDPHits, cacheDPMisses;
        public final double[] bestDOFValues;  // optimal DOF values from backtracking

        public Result(double energy, long precomputeNs, long dpNs, int[] bestAssignment,
                      int oneBodyHits, int oneBodyMisses, int pairHits, int pairMisses,
                      int dpHits, int dpMisses, double[] bestDOFValues) {
            this.energy = energy;
            this.precomputeTimeNs = precomputeNs;
            this.dpTimeNs = dpNs;
            this.bestGridAssignment = bestAssignment;
            this.cacheOneBodyHits = oneBodyHits;
            this.cacheOneBodyMisses = oneBodyMisses;
            this.cachePairHits = pairHits;
            this.cachePairMisses = pairMisses;
            this.cacheDPHits = dpHits;
            this.cacheDPMisses = dpMisses;
            this.bestDOFValues = bestDOFValues;
        }

        /** Total energy cache hits (one-body + pair) */
        public int totalEnergyCacheHits() {
            return cacheOneBodyHits + cachePairHits;
        }

        /** Total energy cache lookups */
        public int totalEnergyCacheLookups() {
            return cacheOneBodyHits + cacheOneBodyMisses + cachePairHits + cachePairMisses;
        }
    }

    /**
     * Minimize energy for a fixed RC tuple using Grid DP.
     */
    public Result minimize(int[] conf) {
        this.confAssignments = conf;
        // Reset lazy enumeration state
        enumInitialized = false;
        enumPopCount = 0;
        RCTuple rcTuple = new RCTuple(conf);
        pmol = confSpace.makeMolecule(rcTuple);

        // Reset per-call cache stats
        oneBodyHits = 0; oneBodyMisses = 0;
        pairHits = 0; pairMisses = 0;
        dpCacheHits = 0; dpCacheMisses = 0;

        buildPositionDOFMap();
        generateGridPoints();

        long t0 = System.nanoTime();
        precomputeEnergies();
        precomputeTimeNs = System.nanoTime() - t0;

        long t1 = System.nanoTime();
        double energy = runDP();
        dpTimeNs = System.nanoTime() - t1;

        // Convert best grid states to DOF values
        double[] bestDOFs = getBestDOFValues();

        return new Result(energy, precomputeTimeNs, dpTimeNs, bestGridStatePerPosition,
                oneBodyHits, oneBodyMisses, pairHits, pairMisses,
                dpCacheHits, dpCacheMisses, bestDOFs);
    }

    /** Get the optimal DOF values from the best grid assignment */
    public double[] getBestDOFValues() {
        if (bestGridStatePerPosition == null || pmol == null) return null;
        double[] dofValues = new double[pmol.dofs.size()];
        for (int pos = 0; pos < numPositions; pos++) {
            int[] dofIndices = positionDOFIndices[pos];
            if (dofIndices == null) continue;
            int gridState = bestGridStatePerPosition[pos];
            // Decode grid state to per-DOF grid indices
            int remaining = gridState;
            for (int d = dofIndices.length - 1; d >= 0; d--) {
                int gridIdx = remaining % gridSize;
                remaining /= gridSize;
                dofValues[dofIndices[d]] = gridValues[dofIndices[d]][gridIdx];
            }
        }
        return dofValues;
    }

    public void setUseTopNPriorityQueue(boolean use) { this.useTopNPriorityQueue = use; }

    /**
     * Get DOF values for the top-N grid assignments.
     * PQ mode: uses lazy pop (init + N pops). Root-only mode: brute-force enumeration.
     */
    public List<double[]> getTopNDOFValues(int n) {
        if (useTopNPriorityQueue) {
            List<double[]> results = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                double[] dofs = popNextDOFValues();
                if (dofs == null) break;
                results.add(dofs);
            }
            return results;
        } else {
            return getTopNDOFValuesRootOnly(n);
        }
    }

    // ========== Lazy top-N enumeration (BWM*-style on-demand pop) ==========

    /**
     * Initialize top-N enumeration. Pre-sorts lambda states at all edges
     * and seeds the priority queue. Must be called after minimize().
     */
    public void initTopNEnumeration() {
        if (enumInitialized) return;
        if (bestGridStatePerPosition == null || pmol == null || dpTables == null) {
            enumInitialized = true;
            return;
        }

        enumSortedLambdas = new HashMap<>();
        enumLambdaContribs = new HashMap<>();
        enumBestEdgeEnergies = new HashMap<>();
        precomputeSortedLambdasForAllEdges(rootEdge, enumSortedLambdas, enumLambdaContribs, enumBestEdgeEnergies);

        enumPQ = new PriorityQueue<>(Comparator.comparingDouble(a -> (Double) a[0]));
        // Seed: for each edge, push rank=1 (2nd-best lambda; rank=0 is the best)
        for (Map.Entry<RootedTreeEdge, int[]> entry : enumSortedLambdas.entrySet()) {
            RootedTreeEdge edge = entry.getKey();
            int[] sorted = entry.getValue();
            if (sorted.length > 1) {
                double bestE = enumBestEdgeEnergies.get(edge);
                double newE = enumLambdaContribs.get(edge)[sorted[1]];
                enumPQ.add(new Object[]{bestTotalEnergy + (newE - bestE), edge, 1});
            }
        }
        enumPopCount = 0;
        enumInitialized = true;
    }

    /**
     * Pop the next-best DOF values from the enumeration.
     * First call returns the best solution (rank 0). Subsequent calls pop from PQ.
     * Returns null when no more solutions are available.
     */
    public double[] popNextDOFValues() {
        if (!enumInitialized) initTopNEnumeration();
        if (bestGridStatePerPosition == null || pmol == null) return null;

        if (enumPopCount == 0) {
            enumPopCount++;
            return getBestDOFValues();
        }

        if (enumPQ == null || enumPQ.isEmpty()) return null;

        Object[] top = enumPQ.poll();
        RootedTreeEdge edge = (RootedTreeEdge) top[1];
        int rank = (Integer) top[2];

        int[] sorted = enumSortedLambdas.get(edge);
        int lambdaState = sorted[rank];
        double[] result = buildAssignmentForEdgeLambda(edge, lambdaState);

        // Push next rank at same edge
        if (rank + 1 < sorted.length) {
            double bestE = enumBestEdgeEnergies.get(edge);
            double nextE = enumLambdaContribs.get(edge)[sorted[rank + 1]];
            enumPQ.add(new Object[]{bestTotalEnergy + (nextE - bestE), edge, rank + 1});
        }

        enumPopCount++;
        return result;
    }

    /** Pre-sort lambda states at each edge for the M-state from the best solution. */
    private void precomputeSortedLambdasForAllEdges(
            RootedTreeEdge edge,
            Map<RootedTreeEdge, int[]> sortedLambdas,
            Map<RootedTreeEdge, double[]> lambdaContribsMap,
            Map<RootedTreeEdge, Double> bestEdgeEnergies) {
        if (edge == null || !edge.getIsLambdaEdge()) return;

        double[] dpTable = dpTables.get(edge);
        if (dpTable == null) return;

        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        // Find M-state from best assignment
        int[] mGridIndices = new int[mPositions.length];
        for (int i = 0; i < mPositions.length; i++) {
            mGridIndices[i] = bestGridStatePerPosition[mPositions[i]];
        }
        int mState = encodeGridState(mPositions, mGridIndices);
        bestEdgeEnergies.put(edge, dpTable[mState]);

        int lambdaStates = computeGridMStates(lambdaPositions);

        RootedTreeNode lc = edge.getCompactLeftChild();
        RootedTreeNode rc = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (lc != null) ? lc.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rc != null) ? rc.getChildOfEdge() : null;

        // Compute energy contribution for each lambda state
        double[] contribs = new double[lambdaStates];
        for (int l = 0; l < lambdaStates; l++) {
            int[] lgIdx = decodeGridState(lambdaPositions, l);
            double e = 0.0;
            for (int i = 0; i < lambdaPositions.length; i++)
                e += oneBodyEnergy[lambdaPositions[i]][lgIdx[i]];
            for (int a = 0; a < lambdaPositions.length; a++)
                for (int b = a + 1; b < lambdaPositions.length; b++)
                    if (interactionGraph.hasEdge(lambdaPositions[a], lambdaPositions[b]))
                        e += getPairEnergy(lambdaPositions[a], lambdaPositions[b], lgIdx[a], lgIdx[b]);
            for (int a = 0; a < lambdaPositions.length; a++)
                for (int b = 0; b < mPositions.length; b++)
                    if (interactionGraph.hasEdge(lambdaPositions[a], mPositions[b]))
                        e += getPairEnergy(lambdaPositions[a], mPositions[b], lgIdx[a], mGridIndices[b]);
            if (leftEdge != null && leftEdge.getIsLambdaEdge()) {
                int childM = projectToChildMState(edge, mGridIndices, lgIdx, leftEdge);
                double[] childDP = dpTables.get(leftEdge);
                if (childDP != null) e += childDP[childM];
            }
            if (rightEdge != null && rightEdge.getIsLambdaEdge()) {
                int childM = projectToChildMState(edge, mGridIndices, lgIdx, rightEdge);
                double[] childDP = dpTables.get(rightEdge);
                if (childDP != null) e += childDP[childM];
            }
            contribs[l] = e;
        }
        lambdaContribsMap.put(edge, contribs);

        // Sort lambda indices by energy
        Integer[] idx = new Integer[lambdaStates];
        for (int i = 0; i < lambdaStates; i++) idx[i] = i;
        Arrays.sort(idx, Comparator.comparingDouble(i -> contribs[i]));
        int[] sorted = new int[lambdaStates];
        for (int i = 0; i < lambdaStates; i++) sorted[i] = idx[i];
        sortedLambdas.put(edge, sorted);

        // Recurse
        if (leftEdge != null) precomputeSortedLambdasForAllEdges(leftEdge, sortedLambdas, lambdaContribsMap, bestEdgeEnergies);
        if (rightEdge != null) precomputeSortedLambdasForAllEdges(rightEdge, sortedLambdas, lambdaContribsMap, bestEdgeEnergies);
    }

    /** Build a complete DOF assignment by overriding lambda at one edge and re-backtracking descendants. */
    private double[] buildAssignmentForEdgeLambda(RootedTreeEdge edge, int lambdaState) {
        int[] gridAssignment = bestGridStatePerPosition.clone();
        int[] lambdaPos = edge.getLambdaPositionsSorted();
        int[] lgIdx = decodeGridState(lambdaPos, lambdaState);
        for (int i = 0; i < lambdaPos.length; i++)
            gridAssignment[lambdaPos[i]] = lgIdx[i];
        // Re-backtrack descendants (child M-states changed)
        RootedTreeNode lc = edge.getCompactLeftChild();
        RootedTreeNode rc = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (lc != null) ? lc.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rc != null) ? rc.getChildOfEdge() : null;
        if (leftEdge != null) backtrackToAssignment(leftEdge, gridAssignment);
        if (rightEdge != null) backtrackToAssignment(rightEdge, gridAssignment);
        return gridAssignmentToDOFs(gridAssignment);
    }

    // ========== Root-edge-only top-N (brute-force, original) ==========

    /**
     * Brute-force top-N: re-enumerates all lambda states at the root edge only.
     */
    private List<double[]> getTopNDOFValuesRootOnly(int n) {
        if (sortedRootMStates == null || pmol == null || dpTables == null)
            return Collections.emptyList();

        // Best root M-state (typically rootMStates=1 when root M-positions have no DOFs)
        int bestRootMState = sortedRootMStates[0];
        int[] rootM = rootEdge.getMPositionsSorted();
        int[] rootMGridIndices = decodeGridState(rootM, bestRootMState);

        // Root M energy (constant across all lambda states)
        double rootMEnergy = 0.0;
        for (int i = 0; i < rootM.length; i++) {
            rootMEnergy += oneBodyEnergy[rootM[i]][rootMGridIndices[i]];
        }
        for (int a = 0; a < rootM.length; a++) {
            for (int b = a + 1; b < rootM.length; b++) {
                if (interactionGraph.hasEdge(rootM[a], rootM[b])) {
                    rootMEnergy += getPairEnergy(rootM[a], rootM[b],
                            rootMGridIndices[a], rootMGridIndices[b]);
                }
            }
        }

        // Re-enumerate ALL lambda states at the root edge
        int[] rootLambda = rootEdge.getLambdaPositionsSorted();
        int lambdaStates = computeGridMStates(rootLambda);

        // Get child edges
        RootedTreeNode leftChild = rootEdge.getCompactLeftChild();
        RootedTreeNode rightChild = rootEdge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;

        // Compute total energy for each lambda state
        int count = Math.min(n, lambdaStates);
        Integer[] lambdaRanking = new Integer[lambdaStates];
        double[] lambdaEnergies = new double[lambdaStates];

        for (int lambdaState = 0; lambdaState < lambdaStates; lambdaState++) {
            lambdaRanking[lambdaState] = lambdaState;
            int[] lambdaGridIndices = decodeGridState(rootLambda, lambdaState);

            double e = rootMEnergy;

            // Lambda one-body
            for (int i = 0; i < rootLambda.length; i++) {
                e += oneBodyEnergy[rootLambda[i]][lambdaGridIndices[i]];
            }
            // Lambda-Lambda pairwise
            for (int a = 0; a < rootLambda.length; a++) {
                for (int b = a + 1; b < rootLambda.length; b++) {
                    if (interactionGraph.hasEdge(rootLambda[a], rootLambda[b])) {
                        e += getPairEnergy(rootLambda[a], rootLambda[b],
                                lambdaGridIndices[a], lambdaGridIndices[b]);
                    }
                }
            }
            // Lambda-M pairwise
            for (int a = 0; a < rootLambda.length; a++) {
                for (int b = 0; b < rootM.length; b++) {
                    if (interactionGraph.hasEdge(rootLambda[a], rootM[b])) {
                        e += getPairEnergy(rootLambda[a], rootM[b],
                                lambdaGridIndices[a], rootMGridIndices[b]);
                    }
                }
            }
            // Children DP tables
            if (leftEdge != null && leftEdge.getIsLambdaEdge()) {
                int leftMState = projectToChildMState(rootEdge, rootMGridIndices, lambdaGridIndices, leftEdge);
                double[] leftDP = dpTables.get(leftEdge);
                if (leftDP != null) e += leftDP[leftMState];
            }
            if (rightEdge != null && rightEdge.getIsLambdaEdge()) {
                int rightMState = projectToChildMState(rootEdge, rootMGridIndices, lambdaGridIndices, rightEdge);
                double[] rightDP = dpTables.get(rightEdge);
                if (rightDP != null) e += rightDP[rightMState];
            }

            lambdaEnergies[lambdaState] = e;
        }

        // Sort lambda states by energy (ascending)
        Arrays.sort(lambdaRanking, Comparator.comparingDouble(idx -> lambdaEnergies[idx]));

        // Build top-N DOF value arrays
        List<double[]> results = new ArrayList<>(count);
        for (int k = 0; k < count; k++) {
            int lambdaState = lambdaRanking[k];
            int[] gridAssignment = new int[numPositions];

            // Set root M-positions
            for (int i = 0; i < rootM.length; i++) {
                gridAssignment[rootM[i]] = rootMGridIndices[i];
            }
            // Set root lambda-positions
            int[] lambdaGridIndices = decodeGridState(rootLambda, lambdaState);
            for (int i = 0; i < rootLambda.length; i++) {
                gridAssignment[rootLambda[i]] = lambdaGridIndices[i];
            }
            // Backtrack children to fill remaining positions
            if (leftEdge != null) backtrackToAssignment(leftEdge, gridAssignment);
            if (rightEdge != null) backtrackToAssignment(rightEdge, gridAssignment);

            results.add(gridAssignmentToDOFs(gridAssignment));
        }

        return results;
    }

    /**
     * Backtrack through the DP tree to fill in optimal lambda states for a given grid assignment.
     * Uses the full numPositions-sized array for correct indexing.
     */
    private void backtrackToAssignment(RootedTreeEdge edge, int[] gridAssignment) {
        if (edge == null || !edge.getIsLambdaEdge()) return;

        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        // Encode M-state from the assignment
        int[] mGridIndices = new int[mPositions.length];
        for (int i = 0; i < mPositions.length; i++) {
            mGridIndices[i] = gridAssignment[mPositions[i]];
        }
        int mState = encodeGridState(mPositions, mGridIndices);

        // Get best lambda state for this M-state
        int[] bestLambdaTable = dpBestLambdaState.get(edge);
        int bestLambda = (bestLambdaTable != null) ? bestLambdaTable[mState] : 0;
        int[] lambdaGridIndices = decodeGridState(lambdaPositions, bestLambda);

        // Record lambda positions
        for (int i = 0; i < lambdaPositions.length; i++) {
            gridAssignment[lambdaPositions[i]] = lambdaGridIndices[i];
        }

        // Recurse on children
        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;

        if (leftEdge != null) backtrackToAssignment(leftEdge, gridAssignment);
        if (rightEdge != null) backtrackToAssignment(rightEdge, gridAssignment);
    }

    /** Convert a grid assignment (per-position grid state indices) to DOF values. */
    private double[] gridAssignmentToDOFs(int[] gridAssignment) {
        double[] dofValues = new double[pmol.dofs.size()];
        for (int pos = 0; pos < numPositions; pos++) {
            int[] dofIndices = positionDOFIndices[pos];
            if (dofIndices == null || dofIndices.length == 0) continue;
            int gridState = gridAssignment[pos];
            int remaining = gridState;
            for (int d = dofIndices.length - 1; d >= 0; d--) {
                int gridIdx = remaining % gridSize;
                remaining /= gridSize;
                dofValues[dofIndices[d]] = gridValues[dofIndices[d]][gridIdx];
            }
        }
        return dofValues;
    }

    // ========== Position -> DOF mapping ==========

    private void buildPositionDOFMap() {
        positionDOFIndices = new int[numPositions][];
        numDOFsPerPosition = new int[numPositions];
        gridStatesPerPosition = new int[numPositions];

        // Map each DOF to its position by matching residue PDB number
        Map<String, List<Integer>> resNumToDofIndices = new LinkedHashMap<>();
        for (int d = 0; d < pmol.dofs.size(); d++) {
            DegreeOfFreedom dof = pmol.dofs.get(d);
            String resNum = null;
            if (dof.getResidue() != null) {
                resNum = dof.getResidue().getPDBResNumber();
            }
            if (resNum != null) {
                resNumToDofIndices.computeIfAbsent(resNum, k -> new ArrayList<>()).add(d);
            }
        }

        for (int pos = 0; pos < numPositions; pos++) {
            String posResNum = confSpace.positions.get(pos).resNum;
            List<Integer> dofIndices = resNumToDofIndices.getOrDefault(posResNum, Collections.emptyList());
            positionDOFIndices[pos] = dofIndices.stream().mapToInt(Integer::intValue).toArray();
            numDOFsPerPosition[pos] = positionDOFIndices[pos].length;
            gridStatesPerPosition[pos] = intPow(gridSize, numDOFsPerPosition[pos]);
        }
    }

    // ========== Grid generation ==========

    private void generateGridPoints() {
        gridValues = new double[pmol.dofs.size()][gridSize];
        for (int d = 0; d < pmol.dofs.size(); d++) {
            double lo = pmol.dofBounds.getMin(d);
            double hi = pmol.dofBounds.getMax(d);
            if (gridSize == 1) {
                gridValues[d][0] = (lo + hi) / 2.0;
            } else {
                for (int g = 0; g < gridSize; g++) {
                    gridValues[d][g] = lo + g * (hi - lo) / (gridSize - 1);
                }
            }
        }
    }

    // ========== Energy precomputation ==========

    private void precomputeEnergies() {
        // TODO: parallel precompute disabled — overhead dominates at small g with cache
        // To re-enable: if (numThreads > 1) { precomputeEnergiesParallel(); } else { ... }
        precomputeOneBody();
        precomputePairwise();
    }

    /**
     * Parallel precomputation of one-body and pairwise energy tables.
     * Uses ThreadLocal so each thread creates only ONE ParametricMolecule + ResPairCache,
     * reused across all tasks assigned to that thread.
     */
    private void precomputeEnergiesParallel() {
        RCTuple rcTuple = new RCTuple(confAssignments);
        ThreadLocal<ParametricMolecule> threadPmol = ThreadLocal.withInitial(() -> confSpace.makeMolecule(rcTuple));
        ThreadLocal<ResPairCache> threadCache = ThreadLocal.withInitial(() -> new ResPairCache(ffparams, atomConnectivity));

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        try {
            precomputeOneBodyParallel(executor, threadPmol, threadCache);
            precomputePairwiseParallel(executor, threadPmol, threadCache);
        } finally {
            executor.shutdown();
        }
    }

    /**
     * Hook for one-body energy evaluation during grid precompute.
     * Subclasses can override this to replace forcefield calls with surrogate models.
     */
    protected double evalOneBodyEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                       int pos, int rc, int gridState, double[] posDofValues,
                                       double erefOffset) {
        return efunc.getEnergy() + erefOffset;
    }

    /**
     * Hook for pairwise energy evaluation during grid precompute.
     * Subclasses can override this to replace forcefield calls with surrogate models.
     */
    protected double evalPairEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                    int pos1, int rc1, int gridState1, double[] pos1DofValues,
                                    int pos2, int rc2, int gridState2, double[] pos2DofValues) {
        return efunc.getEnergy();
    }

    /**
     * Hook for final one-body energy evaluation in continuous relaxation.
     * Default keeps current forcefield behavior.
     */
    protected double evalFinalOneBodyEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                            int pos, int rc, double[] posDofValues, double erefOffset) {
        return efunc.getEnergy() + erefOffset;
    }

    /**
     * Hook for final pairwise energy evaluation in continuous relaxation.
     * Default keeps current forcefield behavior.
     */
    protected double evalFinalPairEnergy(ParametricMolecule mol, ResidueForcefieldEnergy efunc,
                                         int pos1, int rc1, double[] pos1DofValues,
                                         int pos2, int rc2, double[] pos2DofValues) {
        return efunc.getEnergy();
    }

    private void precomputeOneBodyParallel(ExecutorService executor,
                                            ThreadLocal<ParametricMolecule> threadPmol,
                                            ThreadLocal<ResPairCache> threadCache) {
        oneBodyEnergy = new double[numPositions][];

        // Identify which positions need computation (not cached)
        List<Integer> toCompute = new ArrayList<>();
        for (int pos = 0; pos < numPositions; pos++) {
            if (useCache) {
                long key = oneBodyKey(pos, confAssignments[pos]);
                double[] cached = oneBodyCache.get(key);
                if (cached != null) {
                    oneBodyEnergy[pos] = cached;
                    oneBodyHits++;
                    continue;
                }
                oneBodyMisses++;
            }
            toCompute.add(pos);
        }

        if (toCompute.isEmpty()) return;

        // Submit parallel tasks — each task reuses thread-local molecule
        Map<Integer, Future<double[]>> futures = new LinkedHashMap<>();

        for (int pos : toCompute) {
            final int p = pos;
            final int numStates = gridStatesPerPosition[p];
            final int[] dofIndices = positionDOFIndices[p];
            double erefOffset = 0;
            if (eref != null) {
                erefOffset = eref.getOffset(confSpace, p, confAssignments[p]);
            }
            final double erefOff = erefOffset;

            futures.put(p, executor.submit(() -> {
                ParametricMolecule mol = threadPmol.get();
                ResPairCache cache = threadCache.get();

                ResidueInteractions inters = ResInterGen.of(confSpace)
                        .addIntra(p)
                        .addShell(p)
                        .make();
                ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(cache, inters, mol.mol);

                double[] energies = new double[numStates];
                for (int gs = 0; gs < numStates; gs++) {
                    setPositionGridStateOnMol(mol, dofIndices, gs);
                    double[] posDofValues = decodeGridStateToDofValues(dofIndices, gs);
                    energies[gs] = evalOneBodyEnergy(mol, efunc, p, confAssignments[p], gs, posDofValues, erefOff);
                }
                return energies;
            }));
        }

        // Collect results on main thread
        for (Map.Entry<Integer, Future<double[]>> entry : futures.entrySet()) {
            int pos = entry.getKey();
            try {
                oneBodyEnergy[pos] = entry.getValue().get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel one-body precompute failed for pos " + pos, e);
            }
            if (useCache) {
                long key = oneBodyKey(pos, confAssignments[pos]);
                oneBodyCache.put(key, oneBodyEnergy[pos]);
            }
        }
    }

    private void precomputePairwiseParallel(ExecutorService executor,
                                             ThreadLocal<ParametricMolecule> threadPmol,
                                             ThreadLocal<ResPairCache> threadCache) {
        pairEnergyTables = new HashMap<>();

        // Collect edges and check cache
        List<int[]> edgesToCompute = new ArrayList<>();
        for (int i = 0; i < numPositions; i++) {
            for (int j = i + 1; j < numPositions; j++) {
                if (!interactionGraph.hasEdge(i, j)) continue;

                if (useCache) {
                    long key = pairKey(i, confAssignments[i], j, confAssignments[j]);
                    double[][] cached = pairCache.get(key);
                    if (cached != null) {
                        pairEnergyTables.put(edgeKey(i, j), cached);
                        pairHits++;
                        continue;
                    }
                    pairMisses++;
                }
                edgesToCompute.add(new int[]{i, j});
            }
        }

        if (edgesToCompute.isEmpty()) return;

        // Submit parallel tasks
        Map<long[], Future<double[][]>> futures = new LinkedHashMap<>();

        for (int[] edge : edgesToCompute) {
            final int ei = edge[0], ej = edge[1];
            final int statesI = gridStatesPerPosition[ei];
            final int statesJ = gridStatesPerPosition[ej];
            final int[] dofIndicesI = positionDOFIndices[ei];
            final int[] dofIndicesJ = positionDOFIndices[ej];

            long[] edgeInfo = new long[]{ei, ej};
            futures.put(edgeInfo, executor.submit(() -> {
                ParametricMolecule mol = threadPmol.get();
                ResPairCache cache = threadCache.get();

                ResidueInteractions inters = ResInterGen.of(confSpace)
                        .addInter(ei, ej)
                        .make();
                ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(cache, inters, mol.mol);

                double[][] table = new double[statesI][statesJ];
                for (int gi = 0; gi < statesI; gi++) {
                    setPositionGridStateOnMol(mol, dofIndicesI, gi);
                    for (int gj = 0; gj < statesJ; gj++) {
                        setPositionGridStateOnMol(mol, dofIndicesJ, gj);
                        double[] dofsI = decodeGridStateToDofValues(dofIndicesI, gi);
                        double[] dofsJ = decodeGridStateToDofValues(dofIndicesJ, gj);
                        table[gi][gj] = evalPairEnergy(
                                mol, efunc,
                                ei, confAssignments[ei], gi, dofsI,
                                ej, confAssignments[ej], gj, dofsJ
                        );
                    }
                }
                return table;
            }));
        }

        // Collect results
        for (Map.Entry<long[], Future<double[][]>> entry : futures.entrySet()) {
            int ei = (int) entry.getKey()[0];
            int ej = (int) entry.getKey()[1];
            try {
                double[][] table = entry.getValue().get();
                pairEnergyTables.put(edgeKey(ei, ej), table);
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel pairwise precompute failed for edge (" + ei + "," + ej + ")", e);
            }
            if (useCache) {
                long key = pairKey(ei, confAssignments[ei], ej, confAssignments[ej]);
                pairCache.put(key, pairEnergyTables.get(edgeKey(ei, ej)));
            }
        }
    }

    /** Apply a grid state to a specific molecule (thread-safe: only touches the given pmol). */
    private void setPositionGridStateOnMol(ParametricMolecule mol, int[] dofIndices, int gridState) {
        int numDOFs = dofIndices.length;
        if (numDOFs == 0) return;
        int remaining = gridState;
        for (int d = numDOFs - 1; d >= 0; d--) {
            int gridPointIdx = remaining % gridSize;
            remaining /= gridSize;
            mol.dofs.get(dofIndices[d]).apply(gridValues[dofIndices[d]][gridPointIdx]);
        }
    }

    private static long oneBodyKey(int pos, int rc) {
        return ((long) pos << 32) | (rc & 0xFFFFFFFFL);
    }

    private static long pairKey(int pos_i, int rc_i, int pos_j, int rc_j) {
        return ((long) pos_i << 48) | ((long) rc_i << 32) | ((long) pos_j << 16) | (rc_j & 0xFFFFL);
    }

    private void precomputeOneBody() {
        oneBodyEnergy = new double[numPositions][];
        for (int pos = 0; pos < numPositions; pos++) {
            int numStates = gridStatesPerPosition[pos];

            // Try cache lookup
            if (useCache) {
                long key = oneBodyKey(pos, confAssignments[pos]);
                double[] cached = oneBodyCache.get(key);
                if (cached != null) {
                    oneBodyEnergy[pos] = cached;
                    oneBodyHits++;
                    continue;
                }
                oneBodyMisses++;
            }

            oneBodyEnergy[pos] = new double[numStates];

            // Reference energy offset (negated, so adding it subtracts the reference)
            double erefOffset = 0;
            if (eref != null) {
                erefOffset = eref.getOffset(confSpace, pos, confAssignments[pos]);
            }

            // Create energy function for this position (intra + shell)
            ResidueInteractions inters = ResInterGen.of(confSpace)
                    .addIntra(pos)
                    .addShell(pos)
                    .make();
            ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, pmol.mol);

            for (int gs = 0; gs < numStates; gs++) {
                setPositionGridState(pos, gs);
                double[] posDofValues = decodeGridStateToDofValues(positionDOFIndices[pos], gs);
                oneBodyEnergy[pos][gs] = evalOneBodyEnergy(
                        pmol, efunc, pos, confAssignments[pos], gs, posDofValues, erefOffset
                );
            }

            // Store in cache
            if (useCache) {
                long key = oneBodyKey(pos, confAssignments[pos]);
                oneBodyCache.put(key, oneBodyEnergy[pos]);
            }
        }
    }

    private void precomputePairwise() {
        pairEnergyTables = new HashMap<>();
        for (int i = 0; i < numPositions; i++) {
            for (int j = i + 1; j < numPositions; j++) {
                if (!interactionGraph.hasEdge(i, j)) continue;

                // Try cache lookup
                if (useCache) {
                    long key = pairKey(i, confAssignments[i], j, confAssignments[j]);
                    double[][] cached = pairCache.get(key);
                    if (cached != null) {
                        pairEnergyTables.put(edgeKey(i, j), cached);
                        pairHits++;
                        continue;
                    }
                    pairMisses++;
                }

                int statesI = gridStatesPerPosition[i];
                int statesJ = gridStatesPerPosition[j];
                double[][] table = new double[statesI][statesJ];

                ResidueInteractions inters = ResInterGen.of(confSpace)
                        .addInter(i, j)
                        .make();
                ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, pmol.mol);

                for (int gi = 0; gi < statesI; gi++) {
                    setPositionGridState(i, gi);
                    for (int gj = 0; gj < statesJ; gj++) {
                        setPositionGridState(j, gj);
                        double[] dofsI = decodeGridStateToDofValues(positionDOFIndices[i], gi);
                        double[] dofsJ = decodeGridStateToDofValues(positionDOFIndices[j], gj);
                        table[gi][gj] = evalPairEnergy(
                                pmol, efunc,
                                i, confAssignments[i], gi, dofsI,
                                j, confAssignments[j], gj, dofsJ
                        );
                    }
                }

                pairEnergyTables.put(edgeKey(i, j), table);

                // Store in cache
                if (useCache) {
                    long key = pairKey(i, confAssignments[i], j, confAssignments[j]);
                    pairCache.put(key, table);
                }
            }
        }
    }

    // ========== DOF setting ==========

    private void setPositionGridState(int pos, int gridState) {
        int[] dofIndices = positionDOFIndices[pos];
        int numDOFs = dofIndices.length;
        if (numDOFs == 0) return;

        // Decode gridState into per-DOF grid point indices (mixed-radix)
        int remaining = gridState;
        for (int d = numDOFs - 1; d >= 0; d--) {
            int gridPointIdx = remaining % gridSize;
            remaining /= gridSize;
            pmol.dofs.get(dofIndices[d]).apply(gridValues[dofIndices[d]][gridPointIdx]);
        }
    }

    /**
     * Decode a position-local grid state into DOF values.
     * Returned values are in the same order as dofIndices.
     */
    private double[] decodeGridStateToDofValues(int[] dofIndices, int gridState) {
        int numDOFs = dofIndices.length;
        double[] values = new double[numDOFs];
        if (numDOFs == 0) return values;

        int remaining = gridState;
        for (int d = numDOFs - 1; d >= 0; d--) {
            int gridPointIdx = remaining % gridSize;
            remaining /= gridSize;
            values[d] = gridValues[dofIndices[d]][gridPointIdx];
        }
        return values;
    }

    private double[] extractPositionDofValuesFromAllDofs(int pos, double[] allDofs) {
        int[] dofIndices = positionDOFIndices[pos];
        double[] values = new double[dofIndices.length];
        for (int d = 0; d < dofIndices.length; d++) {
            values[d] = allDofs[dofIndices[d]];
        }
        return values;
    }

    // ========== Subtree DP cache support ==========

    /** Precompute which positions affect each edge's DP table (M ∪ λ ∪ descendants). */
    private void precomputeEdgeSubtreePositions(RootedTreeEdge edge) {
        if (edge == null || !edge.getIsLambdaEdge()) return;

        TreeSet<Integer> positions = new TreeSet<>();
        collectSubtreePositions(edge, positions);
        edgeSubtreePositions.put(edge, positions.stream().mapToInt(Integer::intValue).toArray());

        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;
        if (leftEdge != null) precomputeEdgeSubtreePositions(leftEdge);
        if (rightEdge != null) precomputeEdgeSubtreePositions(rightEdge);
    }

    private void collectSubtreePositions(RootedTreeEdge edge, Set<Integer> positions) {
        if (edge == null || !edge.getIsLambdaEdge()) return;
        for (int pos : edge.getMPositionsSorted()) positions.add(pos);
        for (int pos : edge.getLambdaPositionsSorted()) positions.add(pos);

        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;
        if (leftEdge != null) collectSubtreePositions(leftEdge, positions);
        if (rightEdge != null) collectSubtreePositions(rightEdge, positions);
    }

    /** Polynomial hash of RC assignments at the given positions. */
    private long subtreeRCKey(int[] subtreePos) {
        long key = 0;
        for (int pos : subtreePos) {
            key = key * 1000003L + confAssignments[pos];
        }
        return key;
    }

    // ========== DP on branch decomposition tree ==========

    private double runDP() {
        dpTables = new HashMap<>();
        bestGridStatePerPosition = new int[numPositions];

        // Bottom-up DP on compact tree — also track best lambda states per edge
        dpBestLambdaState = new HashMap<>();
        dpPostOrder(rootEdge);

        // At root: enumerate root M-states, add M one-body and M-M pairwise
        int[] rootM = rootEdge.getMPositionsSorted();
        int rootMStates = computeGridMStates(rootM);
        double[] rootDP = dpTables.get(rootEdge);

        // Compute energy for each root M-state
        double[] rootMEnergies = new double[rootMStates];
        for (int mState = 0; mState < rootMStates; mState++) {
            double e = (rootDP != null) ? rootDP[mState] : 0.0;

            // Decode M-state to per-position grid indices
            int[] mGridIndices = decodeGridState(rootM, mState);

            // Add root M one-body
            for (int i = 0; i < rootM.length; i++) {
                e += oneBodyEnergy[rootM[i]][mGridIndices[i]];
            }

            // Add root M-M pairwise
            for (int a = 0; a < rootM.length; a++) {
                for (int b = a + 1; b < rootM.length; b++) {
                    if (interactionGraph.hasEdge(rootM[a], rootM[b])) {
                        e += getPairEnergy(rootM[a], rootM[b], mGridIndices[a], mGridIndices[b]);
                    }
                }
            }

            rootMEnergies[mState] = e;
        }

        // Sort root M-states by energy (ascending) for top-N access
        Integer[] sortedIdx = new Integer[rootMStates];
        for (int i = 0; i < rootMStates; i++) sortedIdx[i] = i;
        Arrays.sort(sortedIdx, Comparator.comparingDouble(idx -> rootMEnergies[idx]));
        this.sortedRootMStates = new int[rootMStates];
        for (int i = 0; i < rootMStates; i++) this.sortedRootMStates[i] = sortedIdx[i];

        // Best solution
        int bestRootMState = sortedRootMStates[0];
        double bestEnergy = rootMEnergies[bestRootMState];
        this.bestTotalEnergy = bestEnergy;

        // Backtrack: record best grid states for root M positions
        int[] bestRootMGrid = decodeGridState(rootM, bestRootMState);
        for (int i = 0; i < rootM.length; i++) {
            bestGridStatePerPosition[rootM[i]] = bestRootMGrid[i];
        }

        // Backtrack through tree to recover lambda states
        backtrackBestStates(rootEdge, bestRootMGrid);

        return bestEnergy;
    }

    // Per-edge: for each M-state, which lambda-state was best
    private HashMap<RootedTreeEdge, int[]> dpBestLambdaState;

    private void backtrackBestStates(RootedTreeEdge edge, int[] parentMGridIndices) {
        if (edge == null || !edge.getIsLambdaEdge()) return;

        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();

        // Find this edge's M-state from parent's perspective
        int mState = encodeGridState(mPositions, projectPositionGridIndices(mPositions, parentMGridIndices));

        // Get best lambda state for this M-state
        int[] bestLambdaTable = dpBestLambdaState.get(edge);
        int bestLambda = (bestLambdaTable != null) ? bestLambdaTable[mState] : 0;
        int[] lambdaGridIndices = decodeGridState(lambdaPositions, bestLambda);

        // Record lambda positions' best grid states
        for (int i = 0; i < lambdaPositions.length; i++) {
            bestGridStatePerPosition[lambdaPositions[i]] = lambdaGridIndices[i];
        }

        // Build combined M + lambda grid indices for children
        int[] combinedGrid = new int[numPositions];
        for (int i = 0; i < mPositions.length; i++) {
            combinedGrid[mPositions[i]] = projectPositionGridIndices(mPositions, parentMGridIndices)[i];
        }
        for (int i = 0; i < lambdaPositions.length; i++) {
            combinedGrid[lambdaPositions[i]] = lambdaGridIndices[i];
        }

        // Recurse on children
        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;

        if (leftEdge != null) backtrackBestStates(leftEdge, combinedGrid);
        if (rightEdge != null) backtrackBestStates(rightEdge, combinedGrid);
    }

    /** Extract grid indices for specific positions from a full-size array */
    /** Extract grid indices for specific positions from a full-size (numPositions) array */
    private int[] projectPositionGridIndices(int[] positions, int[] fullGridIndices) {
        int[] result = new int[positions.length];
        for (int i = 0; i < positions.length; i++) {
            result[i] = fullGridIndices[positions[i]];
        }
        return result;
    }

    private void dpPostOrder(RootedTreeEdge edge) {
        if (edge == null) return;

        // If not a lambda edge, skip (no new positions, DP table is all zeros)
        if (!edge.getIsLambdaEdge()) return;

        // Check DP cache: if the RC assignments at all subtree positions match, reuse the table
        if (useCache) {
            int[] subtreePos = edgeSubtreePositions.get(edge);
            long rcKey = subtreeRCKey(subtreePos);
            HashMap<Long, double[]> edgeDPMap = dpCache.get(edge);
            if (edgeDPMap != null) {
                double[] cached = edgeDPMap.get(rcKey);
                if (cached != null) {
                    dpTables.put(edge, cached);
                    dpCacheHits++;
                    return;  // Entire subtree result is cached, skip children
                }
            }
            dpCacheMisses++;
        }

        // Recurse on compact children
        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();

        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;

        if (leftEdge != null) dpPostOrder(leftEdge);
        if (rightEdge != null) dpPostOrder(rightEdge);

        // Compute DP table for this edge
        int[] mPositions = edge.getMPositionsSorted();
        int[] lambdaPositions = edge.getLambdaPositionsSorted();
        int mStates = computeGridMStates(mPositions);
        int lambdaStates = computeGridMStates(lambdaPositions);

        double[] dpTable = new double[mStates];
        int[] bestLambdaForM = new int[mStates];
        Arrays.fill(dpTable, Double.POSITIVE_INFINITY);

        for (int mState = 0; mState < mStates; mState++) {
            int[] mGridIndices = decodeGridState(mPositions, mState);

            for (int lambdaState = 0; lambdaState < lambdaStates; lambdaState++) {
                int[] lambdaGridIndices = decodeGridState(lambdaPositions, lambdaState);

                double e = 0.0;

                // Lambda one-body
                for (int i = 0; i < lambdaPositions.length; i++) {
                    e += oneBodyEnergy[lambdaPositions[i]][lambdaGridIndices[i]];
                }

                // Lambda-Lambda pairwise
                for (int a = 0; a < lambdaPositions.length; a++) {
                    for (int b = a + 1; b < lambdaPositions.length; b++) {
                        if (interactionGraph.hasEdge(lambdaPositions[a], lambdaPositions[b])) {
                            e += getPairEnergy(lambdaPositions[a], lambdaPositions[b],
                                    lambdaGridIndices[a], lambdaGridIndices[b]);
                        }
                    }
                }

                // Lambda-M pairwise
                for (int a = 0; a < lambdaPositions.length; a++) {
                    for (int b = 0; b < mPositions.length; b++) {
                        if (interactionGraph.hasEdge(lambdaPositions[a], mPositions[b])) {
                            e += getPairEnergy(lambdaPositions[a], mPositions[b],
                                    lambdaGridIndices[a], mGridIndices[b]);
                        }
                    }
                }

                // Children DP tables
                if (leftEdge != null && leftEdge.getIsLambdaEdge()) {
                    int leftMState = projectToChildMState(edge, mGridIndices, lambdaGridIndices, leftEdge);
                    double[] leftDP = dpTables.get(leftEdge);
                    if (leftDP != null) {
                        e += leftDP[leftMState];
                    }
                }
                if (rightEdge != null && rightEdge.getIsLambdaEdge()) {
                    int rightMState = projectToChildMState(edge, mGridIndices, lambdaGridIndices, rightEdge);
                    double[] rightDP = dpTables.get(rightEdge);
                    if (rightDP != null) {
                        e += rightDP[rightMState];
                    }
                }

                if (e < dpTable[mState]) {
                    dpTable[mState] = e;
                    bestLambdaForM[mState] = lambdaState;
                }
            }
        }

        dpTables.put(edge, dpTable);
        dpBestLambdaState.put(edge, bestLambdaForM);

        // Store in DP cache
        if (useCache) {
            int[] subtreePos = edgeSubtreePositions.get(edge);
            long rcKey = subtreeRCKey(subtreePos);
            dpCache.computeIfAbsent(edge, k -> new HashMap<>()).put(rcKey, dpTable);
        }
    }

    // ========== State encoding/decoding ==========

    private int computeGridMStates(int[] positions) {
        int total = 1;
        for (int pos : positions) {
            total *= gridStatesPerPosition[pos];
        }
        return total;
    }

    /**
     * Decode a state index into per-position grid state indices.
     * Mixed-radix decoding: positions[0] is the most significant digit.
     */
    private int[] decodeGridState(int[] positions, int stateIndex) {
        int[] gridIndices = new int[positions.length];
        int remaining = stateIndex;
        for (int i = positions.length - 1; i >= 0; i--) {
            int states = gridStatesPerPosition[positions[i]];
            gridIndices[i] = remaining % states;
            remaining /= states;
        }
        return gridIndices;
    }

    private int encodeGridState(int[] positions, int[] gridIndices) {
        int index = 0;
        int stride = 1;
        for (int i = positions.length - 1; i >= 0; i--) {
            index += gridIndices[i] * stride;
            stride *= gridStatesPerPosition[positions[i]];
        }
        return index;
    }

    /**
     * Project parent's (M-state, lambda-state) to a child edge's M-state.
     * Child's M-positions are a subset of parent's M ∪ lambda positions.
     */
    private int projectToChildMState(RootedTreeEdge parentEdge,
                                     int[] parentMGridIndices, int[] parentLambdaGridIndices,
                                     RootedTreeEdge childEdge) {
        int[] parentMPos = parentEdge.getMPositionsSorted();
        int[] parentLambdaPos = parentEdge.getLambdaPositionsSorted();
        int[] childMPos = childEdge.getMPositionsSorted();

        int[] childMGridIndices = new int[childMPos.length];

        for (int i = 0; i < childMPos.length; i++) {
            int targetPos = childMPos[i];
            boolean found = false;

            // Search in parent's M-positions
            for (int j = 0; j < parentMPos.length; j++) {
                if (parentMPos[j] == targetPos) {
                    childMGridIndices[i] = parentMGridIndices[j];
                    found = true;
                    break;
                }
            }
            if (found) continue;

            // Search in parent's lambda-positions
            for (int j = 0; j < parentLambdaPos.length; j++) {
                if (parentLambdaPos[j] == targetPos) {
                    childMGridIndices[i] = parentLambdaGridIndices[j];
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new RuntimeException("Child M-position " + targetPos +
                        " not found in parent M or lambda");
            }
        }

        return encodeGridState(childMPos, childMGridIndices);
    }

    // ========== Pairwise energy lookup ==========

    private double getPairEnergy(int pos1, int pos2, int gridState1, int gridState2) {
        if (pos1 > pos2) {
            return getPairEnergy(pos2, pos1, gridState2, gridState1);
        }
        double[][] table = pairEnergyTables.get(edgeKey(pos1, pos2));
        if (table == null) return 0.0;
        return table[gridState1][gridState2];
    }

    private static long edgeKey(int i, int j) {
        return ((long) Math.min(i, j) << 32) | Math.max(i, j);
    }

    // ========== Utility ==========

    private static int intPow(int base, int exp) {
        int result = 1;
        for (int i = 0; i < exp; i++) result *= base;
        return result;
    }

    // ========== Continuous DP: Quadratic Relaxation ==========

    /** Quadratic energy model: E(x) = x^T A x + g^T x + c */
    static class QuadraticForm {
        final int[] dofIndices;  // sorted global DOF indices this QF covers
        final double[][] A;      // n×n symmetric Hessian
        final double[] g;        // n-vector linear term
        double c;                // constant

        QuadraticForm(int[] dofIndices, double[][] A, double[] g, double c) {
            this.dofIndices = dofIndices;
            this.A = A;
            this.g = g;
            this.c = c;
        }

        int dim() { return dofIndices.length; }
    }

    /** Stores Schur complement data for back-substitution at one edge */
    static class SchurResult {
        final int[] mDofs;       // M DOF indices (global)
        final int[] lambdaDofs;  // lambda DOF indices (global)
        final RealMatrix A_LL_inv;
        final RealMatrix A_LM;
        final RealVector g_L;

        SchurResult(int[] mDofs, int[] lambdaDofs,
                    RealMatrix A_LL_inv, RealMatrix A_LM, RealVector g_L) {
            this.mDofs = mDofs;
            this.lambdaDofs = lambdaDofs;
            this.A_LL_inv = A_LL_inv;
            this.A_LM = A_LM;
            this.g_L = g_L;
        }
    }

    /** Result of continuous DP relaxation */
    public static class ContinuousResult {
        public final double energy;           // true energy after full forcefield eval
        public final double quadraticEnergy;   // energy predicted by quadratic model
        public final double[] dofValues;       // optimal DOF values (clamped to bounds)
        public final long continuousDPTimeNs;  // time for fitting + DP + back-sub
        public final long evalTimeNs;          // time for final full-energy eval

        public ContinuousResult(double energy, double quadraticEnergy, double[] dofValues,
                                long continuousDPTimeNs, long evalTimeNs) {
            this.energy = energy;
            this.quadraticEnergy = quadraticEnergy;
            this.dofValues = dofValues;
            this.continuousDPTimeNs = continuousDPTimeNs;
            this.evalTimeNs = evalTimeNs;
        }
    }

    // ========== Quadratic Fitting ==========

    /** Decode a grid state to actual DOF values for a position */
    private double[] decodeGridStateToDOFValues(int pos, int gridState) {
        int[] dofIdx = positionDOFIndices[pos];
        int nDofs = dofIdx.length;
        double[] vals = new double[nDofs];
        int remaining = gridState;
        for (int d = nDofs - 1; d >= 0; d--) {
            int gridPointIdx = remaining % gridSize;
            remaining /= gridSize;
            vals[d] = gridValues[dofIdx[d]][gridPointIdx];
        }
        return vals;
    }

    /**
     * Fit a quadratic E(x) = x^T A x + g^T x + c to one-body energy data at a position.
     * Uses QR decomposition for least squares.
     */
    private QuadraticForm fitOneBodyQuadratic(int pos) {
        int[] dofIdx = positionDOFIndices[pos];
        int nDofs = dofIdx.length;
        int nStates = gridStatesPerPosition[pos];

        if (nDofs == 0) {
            return new QuadraticForm(new int[0], new double[0][0], new double[0],
                    oneBodyEnergy[pos][0]);
        }

        int nParams = nDofs * (nDofs + 1) / 2 + nDofs + 1;
        double[][] design = new double[nStates][nParams];
        double[] energies = new double[nStates];

        for (int gs = 0; gs < nStates; gs++) {
            double[] x = decodeGridStateToDOFValues(pos, gs);
            int col = 0;
            // Quadratic terms: A_aa * x_a^2, 2*A_ab * x_a*x_b (a<b)
            for (int a = 0; a < nDofs; a++) {
                for (int b = a; b < nDofs; b++) {
                    design[gs][col++] = (a == b) ? x[a] * x[b] : 2.0 * x[a] * x[b];
                }
            }
            // Linear terms
            for (int a = 0; a < nDofs; a++) {
                design[gs][col++] = x[a];
            }
            // Constant
            design[gs][col] = 1.0;
            energies[gs] = oneBodyEnergy[pos][gs];
        }

        RealVector params = solveLS(design, energies);
        return extractQuadraticForm(dofIdx, nDofs, params);
    }

    /**
     * Fit a quadratic to pairwise energy data between two positions.
     * The resulting QF is over the combined DOFs of both positions.
     */
    private QuadraticForm fitPairwiseQuadratic(int posI, int posJ) {
        // Normalize so posI < posJ (pair tables are stored as [statesMin][statesMax])
        if (posI > posJ) return fitPairwiseQuadratic(posJ, posI);

        int[] dofIdxI = positionDOFIndices[posI];
        int[] dofIdxJ = positionDOFIndices[posJ];
        int nI = dofIdxI.length;
        int nJ = dofIdxJ.length;
        int nDofs = nI + nJ;

        // Combined DOF indices
        int[] combinedDofs = new int[nDofs];
        System.arraycopy(dofIdxI, 0, combinedDofs, 0, nI);
        System.arraycopy(dofIdxJ, 0, combinedDofs, nI, nJ);

        int statesI = gridStatesPerPosition[posI];
        int statesJ = gridStatesPerPosition[posJ];

        if (nDofs == 0) {
            // Both positions have no DOFs
            double[][] pairTable = pairEnergyTables.get(edgeKey(posI, posJ));
            double e = (pairTable != null) ? pairTable[0][0] : 0.0;
            return new QuadraticForm(new int[0], new double[0][0], new double[0], e);
        }

        int nDataPoints = statesI * statesJ;
        int nParams = nDofs * (nDofs + 1) / 2 + nDofs + 1;
        double[][] design = new double[nDataPoints][nParams];
        double[] energies = new double[nDataPoints];

        double[][] pairTable = pairEnergyTables.get(edgeKey(posI, posJ));

        int row = 0;
        for (int gi = 0; gi < statesI; gi++) {
            double[] xI = decodeGridStateToDOFValues(posI, gi);
            for (int gj = 0; gj < statesJ; gj++) {
                double[] xJ = decodeGridStateToDOFValues(posJ, gj);

                double[] x = new double[nDofs];
                System.arraycopy(xI, 0, x, 0, nI);
                System.arraycopy(xJ, 0, x, nI, nJ);

                int col = 0;
                for (int a = 0; a < nDofs; a++) {
                    for (int b = a; b < nDofs; b++) {
                        design[row][col++] = (a == b) ? x[a] * x[b] : 2.0 * x[a] * x[b];
                    }
                }
                for (int a = 0; a < nDofs; a++) design[row][col++] = x[a];
                design[row][col] = 1.0;

                energies[row] = (pairTable != null) ? pairTable[gi][gj] : 0.0;
                row++;
            }
        }

        RealVector params = solveLS(design, energies);
        return extractQuadraticForm(combinedDofs, nDofs, params);
    }

    /** Solve least squares via QR decomposition */
    private RealVector solveLS(double[][] design, double[] energies) {
        RealMatrix D = new Array2DRowRealMatrix(design);
        RealVector e = new ArrayRealVector(energies);
        return new QRDecomposition(D).getSolver().solve(e);
    }

    /** Extract A, g, c from the parameter vector */
    private QuadraticForm extractQuadraticForm(int[] dofIndices, int nDofs, RealVector params) {
        double[][] A = new double[nDofs][nDofs];
        double[] gVec = new double[nDofs];
        int col = 0;
        for (int a = 0; a < nDofs; a++) {
            for (int b = a; b < nDofs; b++) {
                A[a][b] = params.getEntry(col);
                A[b][a] = params.getEntry(col);
                col++;
            }
        }
        for (int a = 0; a < nDofs; a++) {
            gVec[a] = params.getEntry(col++);
        }
        double cVal = params.getEntry(col);
        return new QuadraticForm(dofIndices.clone(), A, gVec, cVal);
    }

    // ========== Continuous DP on Branch Tree ==========

    /**
     * Add a sub-QF into a larger combined QF (dense arrays).
     * Maps subQF's DOF indices into combinedDofs ordering and adds A, g, c.
     */
    private static double addSubQF(QuadraticForm subQF, int[] combinedDofs,
                                    double[][] A, double[] g, double cAccum) {
        if (subQF.dim() == 0) {
            return cAccum + subQF.c;
        }
        // Build index mapping: subQF.dofIndices[i] → position in combinedDofs
        int[] mapping = new int[subQF.dim()];
        for (int i = 0; i < subQF.dim(); i++) {
            mapping[i] = indexOfDof(combinedDofs, subQF.dofIndices[i]);
        }
        for (int i = 0; i < subQF.dim(); i++) {
            for (int j = 0; j < subQF.dim(); j++) {
                A[mapping[i]][mapping[j]] += subQF.A[i][j];
            }
            g[mapping[i]] += subQF.g[i];
        }
        return cAccum + subQF.c;
    }

    /** Find index of a DOF in an array */
    private static int indexOfDof(int[] arr, int val) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == val) return i;
        }
        throw new RuntimeException("DOF " + val + " not found in combined DOF list");
    }

    /**
     * Bottom-up continuous DP: fits quadratics and does Schur complement elimination.
     * Returns a reduced QuadraticForm over the edge's M DOFs.
     */
    private QuadraticForm continuousDPPostOrder(RootedTreeEdge edge,
                                                 Map<RootedTreeEdge, SchurResult> schurResults) {
        if (edge == null || !edge.getIsLambdaEdge()) return null;

        // Recurse children
        RootedTreeNode leftChild = edge.getCompactLeftChild();
        RootedTreeNode rightChild = edge.getCompactRightChild();
        RootedTreeEdge leftEdge = (leftChild != null) ? leftChild.getChildOfEdge() : null;
        RootedTreeEdge rightEdge = (rightChild != null) ? rightChild.getChildOfEdge() : null;

        QuadraticForm leftQF = (leftEdge != null) ? continuousDPPostOrder(leftEdge, schurResults) : null;
        QuadraticForm rightQF = (rightEdge != null) ? continuousDPPostOrder(rightEdge, schurResults) : null;

        // Collect DOFs for M and lambda positions
        int[] mPos = edge.getMPositionsSorted();
        int[] lambdaPos = edge.getLambdaPositionsSorted();

        List<Integer> mDofList = new ArrayList<>();
        for (int pos : mPos) for (int d : positionDOFIndices[pos]) mDofList.add(d);
        List<Integer> lambdaDofList = new ArrayList<>();
        for (int pos : lambdaPos) for (int d : positionDOFIndices[pos]) lambdaDofList.add(d);

        int nM = mDofList.size();
        int nL = lambdaDofList.size();
        int nTotal = nM + nL;

        // Combined DOF ordering: [M DOFs..., lambda DOFs...]
        int[] combinedDofs = new int[nTotal];
        for (int i = 0; i < nM; i++) combinedDofs[i] = mDofList.get(i);
        for (int i = 0; i < nL; i++) combinedDofs[nM + i] = lambdaDofList.get(i);

        // Initialize combined QF arrays
        double[][] A = new double[nTotal][nTotal];
        double[] g = new double[nTotal];
        double c = 0;

        // Add lambda one-body quadratics
        for (int pos : lambdaPos) {
            c = addSubQF(fitOneBodyQuadratic(pos), combinedDofs, A, g, c);
        }

        // Add lambda-lambda pairwise quadratics
        for (int a = 0; a < lambdaPos.length; a++) {
            for (int b = a + 1; b < lambdaPos.length; b++) {
                if (interactionGraph.hasEdge(lambdaPos[a], lambdaPos[b])) {
                    c = addSubQF(fitPairwiseQuadratic(lambdaPos[a], lambdaPos[b]),
                            combinedDofs, A, g, c);
                }
            }
        }

        // Add lambda-M pairwise quadratics
        for (int a = 0; a < lambdaPos.length; a++) {
            for (int b = 0; b < mPos.length; b++) {
                if (interactionGraph.hasEdge(lambdaPos[a], mPos[b])) {
                    c = addSubQF(fitPairwiseQuadratic(lambdaPos[a], mPos[b]),
                            combinedDofs, A, g, c);
                }
            }
        }

        // Add children's reduced QFs
        if (leftQF != null) c = addSubQF(leftQF, combinedDofs, A, g, c);
        if (rightQF != null) c = addSubQF(rightQF, combinedDofs, A, g, c);

        // If no lambda DOFs to eliminate, return combined QF
        if (nL == 0) {
            int[] mDofs = new int[nM];
            System.arraycopy(combinedDofs, 0, mDofs, 0, nM);
            return new QuadraticForm(mDofs, A, g, c);
        }

        // Schur complement: eliminate lambda DOFs
        RealMatrix Afull = new Array2DRowRealMatrix(A);

        RealMatrix A_MM = (nM > 0) ? Afull.getSubMatrix(0, nM - 1, 0, nM - 1) : null;
        RealMatrix A_ML = (nM > 0) ? Afull.getSubMatrix(0, nM - 1, nM, nTotal - 1) : null;
        RealMatrix A_LM = (nM > 0) ? Afull.getSubMatrix(nM, nTotal - 1, 0, nM - 1) : null;
        RealMatrix A_LL = Afull.getSubMatrix(nM, nTotal - 1, nM, nTotal - 1);

        // Regularize A_LL
        for (int i = 0; i < nL; i++) {
            A_LL.addToEntry(i, i, 1e-10);
        }

        RealVector g_M = (nM > 0) ? new ArrayRealVector(Arrays.copyOfRange(g, 0, nM)) : null;
        RealVector g_L = new ArrayRealVector(Arrays.copyOfRange(g, nM, nTotal));

        // Invert A_LL
        RealMatrix A_LL_inv;
        try {
            A_LL_inv = new LUDecomposition(A_LL).getSolver().getInverse();
        } catch (SingularMatrixException e2) {
            // Fallback: stronger regularization
            for (int i = 0; i < nL; i++) A_LL.addToEntry(i, i, 1e-6);
            A_LL_inv = new LUDecomposition(A_LL).getSolver().getInverse();
        }

        // Store SchurResult for back-substitution
        int[] mDofs = new int[nM];
        System.arraycopy(combinedDofs, 0, mDofs, 0, nM);
        int[] lambdaDofs = new int[nL];
        System.arraycopy(combinedDofs, nM, lambdaDofs, 0, nL);
        schurResults.put(edge, new SchurResult(mDofs, lambdaDofs, A_LL_inv,
                A_LM, g_L));

        if (nM == 0) {
            // All DOFs eliminated — return scalar
            double c_star = c - g_L.dotProduct(A_LL_inv.operate(g_L)) / 4.0;
            return new QuadraticForm(new int[0], new double[0][0], new double[0], c_star);
        }

        // A* = A_MM - A_ML * A_LL_inv * A_LM
        RealMatrix schur = A_MM.subtract(A_ML.multiply(A_LL_inv).multiply(A_LM));
        // g* = g_M - A_ML * A_LL_inv * g_L
        RealVector g_star = g_M.subtract(A_ML.operate(A_LL_inv.operate(g_L)));
        // c* = c - g_L^T * A_LL_inv * g_L / 4
        double c_star = c - g_L.dotProduct(A_LL_inv.operate(g_L)) / 4.0;

        return new QuadraticForm(mDofs, schur.getData(), g_star.toArray(), c_star);
    }

    /**
     * Top-down back-substitution: given M DOF values, recover optimal lambda DOFs.
     */
    private void backSubstituteContinuous(RootedTreeEdge edge, double[] allDOFs,
                                           Map<RootedTreeEdge, SchurResult> schurResults) {
        if (edge == null || !edge.getIsLambdaEdge()) return;

        SchurResult sr = schurResults.get(edge);
        if (sr != null && sr.lambdaDofs.length > 0) {
            // x_lambda* = -A_LL_inv * (A_LM * x_M + g_L / 2)
            RealVector x_M;
            if (sr.mDofs.length > 0) {
                double[] mVals = new double[sr.mDofs.length];
                for (int i = 0; i < sr.mDofs.length; i++) {
                    mVals[i] = allDOFs[sr.mDofs[i]];
                }
                x_M = new ArrayRealVector(mVals);
            } else {
                x_M = new ArrayRealVector(0);
            }

            RealVector rhs = sr.g_L.mapMultiply(0.5);
            if (sr.mDofs.length > 0) {
                rhs = rhs.add(sr.A_LM.operate(x_M));
            }
            RealVector x_L = sr.A_LL_inv.operate(rhs).mapMultiply(-1.0);

            for (int i = 0; i < sr.lambdaDofs.length; i++) {
                allDOFs[sr.lambdaDofs[i]] = x_L.getEntry(i);
            }
        }

        // Recurse children
        RootedTreeNode lc = edge.getCompactLeftChild();
        RootedTreeNode rc = edge.getCompactRightChild();
        if (lc != null && lc.getChildOfEdge() != null)
            backSubstituteContinuous(lc.getChildOfEdge(), allDOFs, schurResults);
        if (rc != null && rc.getChildOfEdge() != null)
            backSubstituteContinuous(rc.getChildOfEdge(), allDOFs, schurResults);
    }

    /**
     * Continuous DP relaxation: fits quadratic models to existing grid data,
     * minimizes analytically via Schur complement on the branch tree,
     * then evaluates the true energy at the optimal DOF values.
     *
     * Must be called after minimize() — reuses precomputed grid energy tables.
     * Requires gridSize >= 3 for meaningful quadratic fitting.
     */
    public ContinuousResult continuousRelax() {
        if (oneBodyEnergy == null || pairEnergyTables == null || pmol == null) {
            throw new IllegalStateException("Must call minimize() first");
        }
        if (gridSize < 3) {
            throw new IllegalArgumentException("continuousRelax requires gridSize >= 3, got " + gridSize);
        }

        long t0 = System.nanoTime();

        // Bottom-up continuous DP with Schur complement elimination
        Map<RootedTreeEdge, SchurResult> schurResults = new HashMap<>();
        QuadraticForm rootQF = continuousDPPostOrder(rootEdge, schurResults);

        // At root: add M one-body + M-M pairwise to get the full root QF
        int[] rootM = rootEdge.getMPositionsSorted();
        List<Integer> rootMDofList = new ArrayList<>();
        for (int pos : rootM) for (int d : positionDOFIndices[pos]) rootMDofList.add(d);
        int nRootM = rootMDofList.size();

        // Build the full root QF by adding M one-body and M-M pairwise
        double[][] rootA;
        double[] rootG;
        double rootC;
        int[] rootMDofs;

        if (rootQF != null && rootQF.dim() > 0) {
            rootA = new double[nRootM][nRootM];
            rootG = new double[nRootM];
            rootMDofs = rootQF.dofIndices;
            // Copy rootQF into arrays
            for (int i = 0; i < rootQF.dim(); i++) {
                for (int j = 0; j < rootQF.dim(); j++) {
                    int ii = indexOfDof(rootMDofs, rootQF.dofIndices[i]);
                    int jj = indexOfDof(rootMDofs, rootQF.dofIndices[j]);
                    rootA[ii][jj] += rootQF.A[i][j];
                }
                int ii = indexOfDof(rootMDofs, rootQF.dofIndices[i]);
                rootG[ii] += rootQF.g[i];
            }
            rootC = rootQF.c;
        } else {
            rootMDofs = new int[nRootM];
            for (int i = 0; i < nRootM; i++) rootMDofs[i] = rootMDofList.get(i);
            rootA = new double[nRootM][nRootM];
            rootG = new double[nRootM];
            rootC = (rootQF != null) ? rootQF.c : 0;
        }

        // Add root M one-body quadratics
        for (int pos : rootM) {
            rootC = addSubQF(fitOneBodyQuadratic(pos), rootMDofs, rootA, rootG, rootC);
        }

        // Add root M-M pairwise quadratics
        for (int a = 0; a < rootM.length; a++) {
            for (int b = a + 1; b < rootM.length; b++) {
                if (interactionGraph.hasEdge(rootM[a], rootM[b])) {
                    rootC = addSubQF(fitPairwiseQuadratic(rootM[a], rootM[b]),
                            rootMDofs, rootA, rootG, rootC);
                }
            }
        }

        // Minimize: 2Ax + g = 0  =>  x* = -A^{-1} g / 2
        double[] allDOFs = new double[pmol.dofs.size()];
        double quadEnergy = rootC;

        if (nRootM > 0) {
            RealMatrix Amat = new Array2DRowRealMatrix(rootA);
            RealVector gvec = new ArrayRealVector(rootG);

            // Regularize
            for (int i = 0; i < nRootM; i++) Amat.addToEntry(i, i, 1e-10);

            try {
                RealVector x_star = new LUDecomposition(Amat).getSolver()
                        .solve(gvec.mapMultiply(-0.5));
                for (int i = 0; i < nRootM; i++) {
                    allDOFs[rootMDofs[i]] = x_star.getEntry(i);
                }
                // Quadratic energy at minimum
                quadEnergy = rootC + gvec.dotProduct(x_star)
                        + x_star.dotProduct(Amat.operate(x_star));
            } catch (SingularMatrixException e) {
                // Fallback: use grid DP solution
                double[] gridDOFs = getBestDOFValues();
                if (gridDOFs != null) System.arraycopy(gridDOFs, 0, allDOFs, 0, gridDOFs.length);
            }
        }

        // Back-substitute to recover lambda DOFs
        backSubstituteContinuous(rootEdge, allDOFs, schurResults);

        long continuousDPTimeNs = System.nanoTime() - t0;

        // Clamp to DOF bounds
        for (int d = 0; d < allDOFs.length; d++) {
            double lo = pmol.dofBounds.getMin(d);
            double hi = pmol.dofBounds.getMax(d);
            allDOFs[d] = Math.max(lo, Math.min(hi, allDOFs[d]));
        }

        // Full energy evaluation
        long t1 = System.nanoTime();
        for (int d = 0; d < allDOFs.length; d++) {
            pmol.dofs.get(d).apply(allDOFs[d]);
        }

        double trueEnergy = 0;
        for (int pos = 0; pos < numPositions; pos++) {
            ResidueInteractions inters = ResInterGen.of(confSpace)
                    .addIntra(pos).addShell(pos).make();
            ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, pmol.mol);
            double erefOffset = (eref != null) ? eref.getOffset(confSpace, pos, confAssignments[pos]) : 0;
            double[] posDofValues = extractPositionDofValuesFromAllDofs(pos, allDOFs);
            trueEnergy += evalFinalOneBodyEnergy(
                    pmol, efunc, pos, confAssignments[pos], posDofValues, erefOffset
            );
        }
        for (int i = 0; i < numPositions; i++) {
            for (int j = i + 1; j < numPositions; j++) {
                if (!interactionGraph.hasEdge(i, j)) continue;
                ResidueInteractions inters = ResInterGen.of(confSpace).addInter(i, j).make();
                ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, pmol.mol);
                double[] dofsI = extractPositionDofValuesFromAllDofs(i, allDOFs);
                double[] dofsJ = extractPositionDofValuesFromAllDofs(j, allDOFs);
                trueEnergy += evalFinalPairEnergy(
                        pmol, efunc, i, confAssignments[i], dofsI, j, confAssignments[j], dofsJ
                );
            }
        }
        long evalTimeNs = System.nanoTime() - t1;

        return new ContinuousResult(trueEnergy, quadEnergy, allDOFs, continuousDPTimeNs, evalTimeNs);
    }
}
