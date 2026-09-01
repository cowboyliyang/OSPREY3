package edu.duke.cs.osprey.branchdp;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.AbstractTupleMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** Exact n-body local-factor accounting in the branch-DP M/lambda split. */
public class TestHigherOrderFullDP {

    private static final double RT = 1.7;

    private static void set(Object target, String field, Object value) {
        try {
            Field declared = RootedTreeEdge.class.getDeclaredField(field);
            declared.setAccessible(true);
            declared.set(target, value);
        } catch (ReflectiveOperationException ex) {
            throw new RuntimeException("set " + field, ex);
        }
    }

    private static void setLegacyNullHigherOrderDefault(EnergyMatrix matrix) {
        try {
            Field declared = AbstractTupleMatrix.class.getDeclaredField(
                    "defaultHigherInteraction");
            declared.setAccessible(true);
            declared.set(matrix, null);
        } catch (ReflectiveOperationException ex) {
            throw new RuntimeException("set legacy higher-order default", ex);
        }
    }

    private static RCs makeRCs(int[] cards) {
        int[][] allowed = new int[cards.length][];
        for (int pos = 0; pos < cards.length; pos++) {
            allowed[pos] = new int[cards[pos]];
            for (int rc = 0; rc < cards[pos]; rc++) {
                allowed[pos][rc] = rc;
            }
        }
        return new RCs(allowed);
    }

    private static InteractionGraph fullyConnected(int numPos) {
        boolean[][] adjacency = new boolean[numPos][numPos];
        for (int pos1 = 0; pos1 < numPos; pos1++) {
            for (int pos2 = 0; pos2 < numPos; pos2++) {
                adjacency[pos1][pos2] = pos1 != pos2;
            }
        }
        try {
            Constructor<InteractionGraph> constructor =
                    InteractionGraph.class.getDeclaredConstructor(
                            int.class, boolean[][].class);
            constructor.setAccessible(true);
            return constructor.newInstance(numPos, adjacency);
        } catch (ReflectiveOperationException ex) {
            throw new RuntimeException("build InteractionGraph", ex);
        }
    }

    private static EnergyMatrix makeEnergyMatrix(int[] cards, double seed) {
        EnergyMatrix emat = new EnergyMatrix(cards.length, cards, 0.0);
        for (int pos = 0; pos < cards.length; pos++) {
            for (int rc = 0; rc < cards[pos]; rc++) {
                emat.setOneBody(pos, rc,
                        seed + 0.11 * pos - 0.07 * rc);
            }
        }
        for (int pos1 = 0; pos1 < cards.length; pos1++) {
            for (int pos2 = pos1 + 1; pos2 < cards.length; pos2++) {
                for (int rc1 = 0; rc1 < cards[pos1]; rc1++) {
                    for (int rc2 = 0; rc2 < cards[pos2]; rc2++) {
                        emat.setPairwise(pos1, rc1, pos2, rc2,
                                0.03 * seed + 0.05 * (rc1 + 1)
                                        * (rc2 + 2) - 0.01 * (pos1 + pos2));
                    }
                }
            }
        }
        for (int rc0 = 0; rc0 < cards[0]; rc0++) {
            for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                    emat.setHigherOrder(new RCTuple(
                            0, rc0, 1, rc1, 2, rc2),
                            0.2 * seed + 0.13 * rc0
                                    - 0.09 * rc1 + 0.17 * rc2);
                }
            }
        }
        for (int rc1 = 0; rc1 < cards[1]; rc1++) {
            for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                for (int rc3 = 0; rc3 < cards[3]; rc3++) {
                    emat.setHigherOrder(new RCTuple(
                            1, rc1, 2, rc2, 3, rc3),
                            -0.15 * seed + 0.08 * rc1
                                    + 0.04 * rc2 - 0.06 * rc3);
                }
            }
        }
        return emat;
    }

    private static double logSumExp(double current, double value) {
        if (current == Double.NEGATIVE_INFINITY) return value;
        double maximum = Math.max(current, value);
        return maximum + Math.log(
                Math.exp(current - maximum) + Math.exp(value - maximum));
    }

    private static double expectedConditionalLogZ(
            EnergyMatrix emat, int fixedRc0, int[] cards) {
        double result = Double.NEGATIVE_INFINITY;
        for (int rc1 = 0; rc1 < cards[1]; rc1++) {
            for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                for (int rc3 = 0; rc3 < cards[3]; rc3++) {
                    int[] conf = {fixedRc0, rc1, rc2, rc3};
                    double localEnergy = emat.getInternalEnergy(
                            new RCTuple(conf))
                            - emat.getOneBody(0, fixedRc0);
                    result = logSumExp(result, -localEnergy / RT);
                }
            }
        }
        return result;
    }

    private static double expectedFullLogZ(
            EnergyMatrix emat, int[] cards) {
        double result = Double.NEGATIVE_INFINITY;
        for (int rc0 = 0; rc0 < cards[0]; rc0++) {
            for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                    for (int rc3 = 0; rc3 < cards[3]; rc3++) {
                        int[] conf = {rc0, rc1, rc2, rc3};
                        double energy = emat.getInternalEnergy(
                                new RCTuple(conf));
                        result = logSumExp(result, -energy / RT);
                    }
                }
            }
        }
        return result;
    }

    private static EnergyMatrix makePairwiseOnlyEnergyMatrix(
            int[] cards, double seed) {
        EnergyMatrix emat = new EnergyMatrix(cards.length, cards, 0.0);
        for (int pos = 0; pos < cards.length; pos++) {
            for (int rc = 0; rc < cards[pos]; rc++) {
                emat.setOneBody(pos, rc,
                        seed + 0.11 * pos - 0.07 * rc);
            }
        }
        for (int pos1 = 0; pos1 < cards.length; pos1++) {
            for (int pos2 = pos1 + 1; pos2 < cards.length; pos2++) {
                for (int rc1 = 0; rc1 < cards[pos1]; rc1++) {
                    for (int rc2 = 0; rc2 < cards[pos2]; rc2++) {
                        emat.setPairwise(pos1, rc1, pos2, rc2,
                                0.03 * seed + 0.05 * (rc1 + 1)
                                        * (rc2 + 2) - 0.01 * (pos1 + pos2));
                    }
                }
            }
        }
        return emat;
    }

    private static double expectedSparseFullLogZ(
            EnergyMatrix pairwise, int[] cards) {
        double result = Double.NEGATIVE_INFINITY;
        for (int rc0 = 0; rc0 < cards[0]; rc0++) {
            for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                    for (int rc3 = 0; rc3 < cards[3]; rc3++) {
                        int[] conf = {rc0, rc1, rc2, rc3};
                        double energy = pairwise.getInternalEnergy(
                                new RCTuple(conf));
                        if (rc0 == 0 && rc1 == 1 && rc2 == 2) {
                            energy += 0.73;
                        }
                        if (rc1 == 0 && rc2 == 1 && rc3 == 1) {
                            energy += 0.41;
                        }
                        result = logSumExp(result, -energy / RT);
                    }
                }
            }
        }
        return result;
    }

    @Test
    public void exactTripleFactorsMatchBruteForceAcrossMAndLambda() {
        System.setProperty("branchdp.dp.gpu", "false");
        System.setProperty("branchdp.dp.nativeKernel", "false");
        System.setProperty("branchdp.dp.parallel", "false");
        System.setProperty("branchdp.dp.progress", "false");

        int[] cards = {2, 2, 3, 2};
        RCs rcs = makeRCs(cards);
        EnergyMatrix rigid = makeEnergyMatrix(cards, 1.0);
        EnergyMatrix minimizing = makeEnergyMatrix(cards, 2.0);

        RootedTreeEdge edge = new RootedTreeEdge(
                null, null, new LinkedHashSet<>(), false, rcs);
        set(edge, "isLambdaEdge", Boolean.TRUE);
        set(edge, "mPositionsSorted", new int[]{0});
        set(edge, "lambdaPositionsSorted", new int[]{1, 2, 3});
        set(edge, "mStateCount", Long.valueOf(cards[0]));
        set(edge, "mArraySize", Integer.valueOf(cards[0]));
        set(edge, "totalLambdaStates",
                Integer.valueOf(cards[1] * cards[2] * cards[3]));
        set(edge, "dpTable", new DenseDPTable(cards[0]));
        set(edge, "Fset", new LinkedHashSet<RootedTreeEdge>());
        set(edge, "cachedRigidEmat", rigid);
        set(edge, "cachedMinEmat", minimizing);
        set(edge, "cachedG", fullyConnected(cards.length));
        set(edge, "cachedRT", Double.valueOf(RT));

        edge.computeFullDP();

        for (int rc0 = 0; rc0 < cards[0]; rc0++) {
            assertEquals(expectedConditionalLogZ(rigid, rc0, cards),
                    edge.getLogZLower(rc0), 1.0e-12,
                    "rigid conditional logZ at M rc " + rc0);
            assertEquals(expectedConditionalLogZ(minimizing, rc0, cards),
                    edge.getLogZUpper(rc0), 1.0e-12,
                    "minimizing conditional logZ at M rc " + rc0);
        }
    }

    @Test
    public void completeBranchTreeCountsEveryTripleFactorExactlyOnce() {
        System.setProperty("branchdp.dp.gpu", "false");
        System.setProperty("branchdp.dp.nativeKernel", "false");
        System.setProperty("branchdp.dp.parallel", "false");
        System.setProperty("branchdp.dp.progress", "false");

        int[] cards = {2, 2, 3, 2};
        RCs rcs = makeRCs(cards);
        EnergyMatrix rigid = makeEnergyMatrix(cards, 1.0);
        EnergyMatrix minimizing = makeEnergyMatrix(cards, 2.0);
        InteractionGraph graph = InteractionGraph.buildComplete(cards.length);

        BranchDecomposition decomposition = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.GREEDY_MERGE, cards);
        decomposition.compute();
        RootedTreeNode root = decomposition.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(root);
        RootedTreeEdge.postOrderInitIncremental(
                root, rigid, minimizing, graph, RT);
        RootedTreeEdge.postOrderComputeFullDP(root);

        RootedTreeEdge rootEdge =
                root.getLeftChild().getChildOfEdge();
        assertEquals(expectedFullLogZ(rigid, cards),
                rootEdge.getLogZLower(0), 1.0e-12,
                "rigid complete-tree logZ");
        assertEquals(expectedFullLogZ(minimizing, cards),
                rootEdge.getLogZUpper(0), 1.0e-12,
                "minimizing complete-tree logZ");
    }

    @Test
    public void copiedLegacyMatrixSparseTripleFactorsTreatMissingAssignmentsAsZero() {
        System.setProperty("branchdp.dp.gpu", "false");
        System.setProperty("branchdp.dp.nativeKernel", "false");
        System.setProperty("branchdp.dp.parallel", "false");
        System.setProperty("branchdp.dp.progress", "false");

        int[] cards = {2, 2, 3, 2};
        RCs rcs = makeRCs(cards);
        EnergyMatrix pairwise = makePairwiseOnlyEnergyMatrix(cards, 1.5);
        // Reproduce matrices serialized by the old ConfSpaceIteration path:
        // their HigherTupleFinder default is null rather than energy zero.
        setLegacyNullHigherOrderDefault(pairwise);
        EnergyMatrix sparse = new EnergyMatrix(pairwise);
        sparse.setHigherOrder(new RCTuple(0, 0, 1, 1, 2, 2), 0.73);
        sparse.setHigherOrder(new RCTuple(1, 0, 2, 1, 3, 1), 0.41);

        int[] absent = {1, 1, 2, 0};
        assertEquals(pairwise.getInternalEnergy(new RCTuple(absent)),
                sparse.getInternalEnergy(new RCTuple(absent)), 0.0,
                "an unstored sparse triple assignment must contribute zero");

        InteractionGraph graph = InteractionGraph.buildComplete(cards.length);
        BranchDecomposition decomposition = new BranchDecomposition(
                graph, BranchDecomposition.Strategy.GREEDY_MERGE, cards);
        decomposition.compute();
        RootedTreeNode root = decomposition.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(root);
        RootedTreeEdge.postOrderInitIncremental(
                root, sparse, sparse, graph, RT);
        RootedTreeEdge.postOrderComputeFullDP(root);

        RootedTreeEdge rootEdge =
                root.getLeftChild().getChildOfEdge();
        double expected = expectedSparseFullLogZ(pairwise, cards);
        assertEquals(expected, rootEdge.getLogZLower(0), 1.0e-12,
                "rigid sparse complete-tree logZ");
        assertEquals(expected, rootEdge.getLogZUpper(0), 1.0e-12,
                "minimizing sparse complete-tree logZ");
    }

    @Test
    public void zeroEnergyFillEdgesMakeANonCliqueTripleExact() {
        System.setProperty("branchdp.dp.gpu", "false");
        System.setProperty("branchdp.dp.nativeKernel", "false");
        System.setProperty("branchdp.dp.parallel", "false");
        System.setProperty("branchdp.dp.progress", "false");

        int[] cards = {2, 2, 2, 2};
        RCs rcs = makeRCs(cards);
        List<int[]> baseEdges = List.of(
                new int[]{0, 1}, new int[]{1, 2}, new int[]{2, 3});
        InteractionGraph baseGraph = InteractionGraph.buildFromEdges(
                cards.length, baseEdges);
        List<int[]> proposalEdges = new ArrayList<>(baseEdges);
        proposalEdges.add(new int[]{0, 3});
        proposalEdges.add(new int[]{1, 3});
        InteractionGraph proposalGraph = InteractionGraph.buildFromEdges(
                cards.length, proposalEdges);

        EnergyMatrix emat = new EnergyMatrix(cards.length, cards, 0.0);
        for (int pos = 0; pos < cards.length; pos++) {
            for (int rc = 0; rc < cards[pos]; rc++) {
                emat.setOneBody(pos, rc, 0.08 * pos - 0.03 * rc);
            }
        }
        for (int[] edge : baseEdges) {
            for (int rc1 = 0; rc1 < cards[edge[0]]; rc1++) {
                for (int rc2 = 0; rc2 < cards[edge[1]]; rc2++) {
                    emat.setPairwise(edge[0], rc1, edge[1], rc2,
                            0.04 * (rc1 + 1) * (rc2 + 2));
                }
            }
        }
        // The two proposal-only fill pairs remain exactly zero.
        for (int rc0 = 0; rc0 < cards[0]; rc0++) {
            for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                for (int rc3 = 0; rc3 < cards[3]; rc3++) {
                    emat.setHigherOrder(new RCTuple(
                            0, rc0, 1, rc1, 3, rc3),
                            ((rc0 ^ rc1 ^ rc3) == 0) ? -0.7 : 0.9);
                }
            }
        }

        assertEquals(3, baseGraph.getNumEdges());
        assertEquals(5, proposalGraph.getNumEdges());
        BranchDecomposition decomposition = new BranchDecomposition(
                proposalGraph, BranchDecomposition.Strategy.GREEDY_MERGE,
                cards);
        decomposition.compute();
        RootedTreeNode root = decomposition.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(root, true);
        RootedTreeEdge rootEdge = root.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();
        RootedTreeEdge.postOrderInitIncremental(
                root, emat, emat, proposalGraph, RT);
        RootedTreeEdge.postOrderComputeFullDP(root);

        double expected = expectedFullLogZ(emat, cards);
        assertEquals(expected, rootEdge.getLogZLower(0), 1.0e-12,
                "fill-edge rigid logZ");
        assertEquals(expected, rootEdge.getLogZUpper(0), 1.0e-12,
                "fill-edge minimizing logZ");
    }
}
