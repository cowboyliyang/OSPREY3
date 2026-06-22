package edu.duke.cs.osprey.energy.approximation.branch;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Sparse conf sampler with a mixed budget:
 *   - branch coverage
 *   - low-energy pool
 *   - global random
 */
public class BranchSparseSampler {

    public static class Split {
        public final double branchCoverageFraction;
        public final double lowEnergyFraction;
        public final double randomFraction;

        public Split(double branchCoverageFraction, double lowEnergyFraction, double randomFraction) {
            double sum = branchCoverageFraction + lowEnergyFraction + randomFraction;
            if (sum <= 0.0) {
                throw new IllegalArgumentException("sampling split fractions must sum to > 0");
            }
            this.branchCoverageFraction = branchCoverageFraction / sum;
            this.lowEnergyFraction = lowEnergyFraction / sum;
            this.randomFraction = randomFraction / sum;
        }

        public static Split defaults() {
            return new Split(0.50, 0.30, 0.20);
        }
    }

    private final SimpleConfSpace confSpace;
    private final InteractionGraph interactionGraph;
    private final PruningMatrix pmat;
    private final Random rand;

    public BranchSparseSampler(
            SimpleConfSpace confSpace,
            InteractionGraph interactionGraph,
            PruningMatrix pmat,
            long randomSeed
    ) {
        this.confSpace = confSpace;
        this.interactionGraph = interactionGraph;
        this.pmat = pmat;
        this.rand = new Random(randomSeed);
    }

    public List<int[]> sample(
            int totalSamples,
            Split split,
            EnergyMatrix lowEnergyEmat,
            int lowEnergyPoolLimit
    ) {
        int total = Math.max(0, totalSamples);
        if (total == 0) {
            return new ArrayList<>();
        }

        int branchBudget = (int) Math.round(total * split.branchCoverageFraction);
        int lowBudget = (int) Math.round(total * split.lowEnergyFraction);
        if (branchBudget + lowBudget > total) {
            lowBudget = Math.max(0, total - branchBudget);
        }
        int randomBudget = total - branchBudget - lowBudget;

        Conf.Set unique = new Conf.Set();
        List<int[]> out = new ArrayList<>(total);

        sampleBranchCoverage(branchBudget, unique, out);
        sampleLowEnergy(lowBudget, lowEnergyEmat, lowEnergyPoolLimit, unique, out);
        sampleRandom(randomBudget, unique, out);

        while (out.size() < total) {
            int[] conf = randomConf();
            if (tryAdd(conf, unique, out)) {
                // continue
            }
        }

        return out;
    }

    private void sampleBranchCoverage(int budget, Conf.Set unique, List<int[]> out) {
        if (budget <= 0) {
            return;
        }
        int targetSize = out.size() + budget;

        List<int[]> edges = interactionGraph == null ? new ArrayList<>() : interactionGraph.getEdgeList();
        if (edges.isEmpty()) {
            sampleRandom(budget, unique, out);
            return;
        }

        Set<Long> seenEdgeStates = new HashSet<>();
        int maxAttempts = Math.max(budget * 20, 200);
        int attempts = 0;

        while (out.size() < targetSize && attempts < maxAttempts) {
            attempts++;

            int edgeIndex = rand.nextInt(edges.size());
            int[] edge = edges.get(edgeIndex);
            int p = edge[0];
            int q = edge[1];

            int rp = rand.nextInt(confSpace.positions.get(p).resConfs.size());
            int rq = rand.nextInt(confSpace.positions.get(q).resConfs.size());

            long edgeState = packEdgeState(edgeIndex, rp, rq);
            int[] conf = randomConf();
            conf[p] = rp;
            conf[q] = rq;

            if (seenEdgeStates.contains(edgeState) && rand.nextDouble() < 0.7) {
                continue;
            }
            if (tryAdd(conf, unique, out)) {
                seenEdgeStates.add(edgeState);
            }
        }

        while (out.size() < targetSize) {
            int[] conf = randomConf();
            tryAdd(conf, unique, out);
        }
    }

    private void sampleLowEnergy(
            int budget,
            EnergyMatrix lowEnergyEmat,
            int lowEnergyPoolLimit,
            Conf.Set unique,
            List<int[]> out
    ) {
        if (budget <= 0 || lowEnergyEmat == null) {
            return;
        }
        int targetSize = out.size() + budget;

        RCs rcs = (pmat != null) ? new RCs(pmat) : new RCs(confSpace);
        ConfAStarTree astar = new ConfAStarTree.Builder(lowEnergyEmat, rcs)
                .setTraditional()
                .build();

        int limit = Math.max(budget, lowEnergyPoolLimit);
        int pulled = 0;
        while (out.size() < targetSize && pulled < limit) {
            ConfSearch.ScoredConf sc = astar.nextConf();
            if (sc == null) {
                break;
            }
            pulled++;
            int[] conf = sc.getAssignments().clone();
            tryAdd(conf, unique, out);
        }
    }

    private void sampleRandom(int budget, Conf.Set unique, List<int[]> out) {
        if (budget <= 0) {
            return;
        }

        int startSize = out.size();
        int maxAttempts = Math.max(budget * 20, 200);
        int attempts = 0;

        while (out.size() < startSize + budget && attempts < maxAttempts) {
            attempts++;
            int[] conf = randomConf();
            tryAdd(conf, unique, out);
        }
    }

    private int[] randomConf() {
        int[] conf = new int[confSpace.positions.size()];
        for (SimpleConfSpace.Position pos : confSpace.positions) {
            conf[pos.index] = rand.nextInt(pos.resConfs.size());
        }
        return conf;
    }

    private boolean tryAdd(int[] conf, Conf.Set unique, List<int[]> out) {
        if (pmat != null && pmat.isPruned(new RCTuple(conf))) {
            return false;
        }
        if (unique.add(conf)) {
            out.add(conf);
            return true;
        }
        return false;
    }

    private static long packEdgeState(int edgeIndex, int rc1, int rc2) {
        long out = (edgeIndex & 0xffffffffL);
        out = (out << 16) ^ (rc1 & 0xffffL);
        out = (out << 16) ^ (rc2 & 0xffffL);
        return out;
    }
}
