package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;

/**
 * Exact small-instance audit of the sparse-vs-full partition-function bias
 * introduced by BranchMARK*'s sparse interaction graph.
 *
 * Enumerates every conformation of a small wild-type-only confspace
 * (configurable via {@code branchmarkstar.sparseFull.numFlexible}) and
 * compares two partition functions over the same conformations:
 *   - full   = sum using all pairwise terms in the minimized energy matrix
 *   - sparse = sum using only pairwise terms retained by the dual-cutoff
 *              {@link InteractionGraph} that BranchMARK* actually uses
 *
 * Reports per-strand log10 Z and the resulting K* gap (sparse − full),
 * which is the systematic error introduced purely by edge pruning.
 *
 * Driven by {@code slurm/scripts/run_sparse_full_gap.slurm}.
 */
public class SparseFullCheck {

    private static final int NUM_CPUs = Integer.getInteger("osprey.branchmarkstar.numCpus", 4);

    private static class SparseFullZ {
        String label;
        int positions;
        long confs;
        int sparseEdges;
        int possibleEdges;
        BigDecimal fullZ = BigDecimal.ZERO;
        BigDecimal sparseZ = BigDecimal.ZERO;
        double fullLog10Z;
        double sparseLog10Z;
    }

    @Test
    public void compareSparseFullPartitionFunctionGap() {
        int numFlexible = Integer.getInteger("branchmarkstar.sparseFull.numFlexible", 3);
        long maxConfs = Long.getLong("branchmarkstar.sparseFull.maxConfs", 1_000_000L);

        TestKStar.ConfSpaces confSpaces = ConfSpaces2RL0.buildWildTypeConfSpace(numFlexible);
        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            SparseFullZ protein = compareSparseFullZ("protein", confSpaces.protein, ecalcMinimized, maxConfs);
            SparseFullZ ligand = compareSparseFullZ("ligand", confSpaces.ligand, ecalcMinimized, maxConfs);
            SparseFullZ complex = compareSparseFullZ("complex", confSpaces.complex, ecalcMinimized, maxConfs);

            double fullK = complex.fullLog10Z - protein.fullLog10Z - ligand.fullLog10Z;
            double sparseK = complex.sparseLog10Z - protein.sparseLog10Z - ligand.sparseLog10Z;

            System.out.println("\n========== SPARSE/FULL PARTITION FUNCTION GAP ==========");
            printSparseFullZ(protein);
            printSparseFullZ(ligand);
            printSparseFullZ(complex);
            System.out.println(String.format(
                    "K* log10(full)=%.6f, log10(sparse)=%.6f, sparse-full=%.6f",
                    fullK, sparseK, sparseK - fullK));
            System.out.println("========================================================");
        }
    }

    private SparseFullZ compareSparseFullZ(String label, SimpleConfSpace confSpace,
                                           EnergyCalculator ecalcMinimized, long maxConfs) {
        ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies())
                .build();
        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                .build()
                .calcEnergyMatrix();

        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
        EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build()
                .calcEnergyMatrix();

        RCs rcs = new RCs(confSpace);
        long confCount = countConfs(rcs);
        if (confCount > maxConfs) {
            throw new IllegalArgumentException(label + " has " + confCount
                    + " conformations; raise -Dbranchmarkstar.sparseFull.maxConfs="
                    + confCount + " if you really want exact enumeration.");
        }

        InteractionGraph sparseGraph = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, 8.0, 0.1);

        SparseFullZ result = new SparseFullZ();
        result.label = label;
        result.positions = rcs.getNumPos();
        result.confs = confCount;
        result.sparseEdges = sparseGraph.getNumEdges();
        result.possibleEdges = result.positions * (result.positions - 1) / 2;

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        int[] conf = new int[rcs.getNumPos()];
        enumerateSparseFullZ(0, conf, rcs, ematMinimized, sparseGraph, bc, result);
        result.fullLog10Z = log10(result.fullZ);
        result.sparseLog10Z = log10(result.sparseZ);
        return result;
    }

    private void enumerateSparseFullZ(int pos, int[] conf, RCs rcs, EnergyMatrix emat,
                                      InteractionGraph sparseGraph, BoltzmannCalculator bc,
                                      SparseFullZ result) {
        if (pos == rcs.getNumPos()) {
            double fullEnergy = confEnergy(conf, emat, null);
            double sparseEnergy = confEnergy(conf, emat, sparseGraph);
            if (Double.isFinite(fullEnergy)) {
                result.fullZ = result.fullZ.add(bc.calc(fullEnergy));
            }
            if (Double.isFinite(sparseEnergy)) {
                result.sparseZ = result.sparseZ.add(bc.calc(sparseEnergy));
            }
            return;
        }
        for (int rci = 0; rci < rcs.getNum(pos); rci++) {
            conf[pos] = rcs.get(pos, rci);
            enumerateSparseFullZ(pos + 1, conf, rcs, emat, sparseGraph, bc, result);
        }
    }

    private double confEnergy(int[] conf, EnergyMatrix emat, InteractionGraph sparseGraph) {
        double energy = emat.getConstTerm();
        for (int i = 0; i < conf.length; i++) {
            energy += emat.getOneBody(i, conf[i]);
            for (int j = i + 1; j < conf.length; j++) {
                if (sparseGraph == null || sparseGraph.hasEdge(i, j)) {
                    energy += emat.getPairwise(i, conf[i], j, conf[j]);
                }
            }
        }
        return energy;
    }

    private long countConfs(RCs rcs) {
        long count = 1;
        for (int pos = 0; pos < rcs.getNumPos(); pos++) {
            count = Math.multiplyExact(count, Math.max(1, rcs.getNum(pos)));
        }
        return count;
    }

    private static double log10(BigDecimal z) {
        if (z.signum() <= 0) return Double.NEGATIVE_INFINITY;
        BigExp exp = new BigExp(z);
        return Math.log10(exp.fp) + exp.exp;
    }

    private static void printSparseFullZ(SparseFullZ z) {
        System.out.println(String.format(
                "%-8s pos=%2d confs=%8d sparseEdges=%3d/%3d log10Z(full)=%12.6f log10Z(sparse)=%12.6f sparse-full=%10.6f",
                z.label, z.positions, z.confs, z.sparseEdges, z.possibleEdges,
                z.fullLog10Z, z.sparseLog10Z, z.sparseLog10Z - z.fullLog10Z));
    }
}
