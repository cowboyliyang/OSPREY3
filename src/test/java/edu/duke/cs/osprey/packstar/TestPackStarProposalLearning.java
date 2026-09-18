package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import org.junit.jupiter.api.Test;
import java.util.*;
import static org.junit.jupiter.api.Assertions.*;

class TestPackStarProposalLearning {
    @Test void jointFitCorrectsTheSumRatherThanEachTermIndependently() {
        int[][] features = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};
        double[] residual = {1, 3, 4, 6};
        var fit = PackStarProposalLearning.jointFit(features, residual,
                new double[4], new double[4], 10, 100);
        for (int i = 0; i < 4; i++) {
            double predicted = fit.offset;
            for (int j : features[i]) predicted += fit.delta[j];
            assertEquals(residual[i], predicted, 1e-7);
        }
    }

    @Test void jointFitPreservesMultiplicityAndSourceWeights() {
        var weighted = PackStarProposalLearning.jointFit(
                new int[][]{{}, {}}, new double[]{0, 4}, new double[]{Math.log(3), 0},
                new double[0], 10, 10);
        var repeated = PackStarProposalLearning.jointFit(
                new int[][]{{}, {}, {}, {}}, new double[]{0, 0, 0, 4}, new double[4],
                new double[0], 10, 10);
        assertEquals(1, weighted.offset, 1e-12);
        assertEquals(weighted.offset, repeated.offset, 1e-12);
    }

    @Test void momentFitRespondsToExponentialTailAndNotMeanEnergy() {
        var cell = new PackStarProposalLearning.MomentCell();
        cell.add(new int[]{0, 0}, 0, -1, 1);
        cell.add(new int[]{0, 1}, 0, 1, 1);
        double expected = -0.5 * Math.log((Math.exp(2) + Math.exp(-2)) / 2);
        assertEquals(expected, cell.correction(0, 1, 2, 0, 10), 1e-12);
        assertTrue(expected < -0.6); // the mean residual is zero
        for (int i = 0; i < 100; i++) cell.add(new int[]{0, 0}, 0, -1, 1);
        assertEquals(0, cell.correction(0, 1, 3, 0, 10)); // still only 2 contexts
    }

    @Test void conditionalMomentUsesSourceWeightsAndStableExponents() {
        var weighted = new PackStarProposalLearning.MomentCell();
        weighted.add(new int[]{0}, Math.log(3), -1000, 1);
        weighted.add(new int[]{1}, 0, 1000, 1);
        var repeated = new PackStarProposalLearning.MomentCell();
        for (int i = 0; i < 3; i++) repeated.add(new int[]{0}, 10000, -1000, 1);
        repeated.add(new int[]{1}, 10000, 1000, 1);
        assertEquals(weighted.correction(0, 1, 1, 0, 2000),
                repeated.correction(0, 1, 1, 0, 2000), 1e-9);
        assertEquals(-3, weighted.correction(0, 1, 1, 0, 3));
    }

    @Test void logRhoIncludesChangedNormalizerAndIsGaugeInvariant() {
        var data = new PackStarProposalLearning.Data(new int[][]{{0}, {1}},
                new double[]{-1, 1}, new double[2]);
        double[] h = {-1, 1};
        assertEquals(0, PackStarProposalLearning.logRho(data, h, -1, 1), 1e-12);
        assertTrue(PackStarProposalLearning.logRho(data, new double[2], -1, 1) > 0);
        assertEquals(0, PackStarProposalLearning.logRho(data,
                new double[]{999, 1001}, -1, 1), 1e-10);
        var shiftedSources = new PackStarProposalLearning.Data(data.conf, data.residual,
                new double[]{10000, 10000});
        assertEquals(0, PackStarProposalLearning.logRho(shiftedSources, h, -1, 1), 1e-10);
    }

    @Test void sparseTripleSelectionRepairsParityWithoutDenseCandidateTables() {
        int[][] allowed = new int[6][];
        for (int p = 0; p < 6; p++) allowed[p] = new int[]{0, 1};
        // Unobserved RC has a defined zero fallback in the selected table.
        allowed[0] = new int[]{0, 1, 2};
        RCs rcs = new RCs(allowed);
        EnergyMatrix base = new EnergyMatrix(6, new int[]{3, 2, 2, 2, 2, 2}, 0);
        for (int a = 0; a < 3; a++) for (int b = 0; b < a; b++)
            for (int ra : allowed[a]) for (int rb : allowed[b]) base.setPairwise(a, ra, b, rb, 1.0);
        InteractionGraph graph = InteractionGraph.buildFromEnergyMatrix(base, base, rcs, 0.5);
        int[][] conf = new int[64][6];
        double[] residual = new double[64];
        for (int i = 0; i < 64; i++) {
            for (int p = 0; p < 6; p++) conf[i][p] = (i >>> p) & 1;
            residual[i] = ((conf[i][0] ^ conf[i][1] ^ conf[i][2]) == 0) ? -1 : 1;
        }
        var data = new PackStarProposalLearning.Data(conf, residual, new double[64]);
        var model = PackStarTripleEtaCorrections.emptyFitted(6).fitSelectedSecondMoment(
                rcs, graph, data, new PackStarProposalLearning.Data[]{data, data},
                1, 1, 0, 1, 0, 10, 100);
        assertEquals(1, model.positionTriples().size());
        assertEquals("0,1,2", model.positionTriples().get(0).toString());
        double[] h = new double[64];
        for (int i = 0; i < 64; i++) h[i] = model.scoreResidual(conf[i], (a, ra, b, rb) -> 0);
        assertEquals(0, PackStarProposalLearning.logRho(data, h, -1, 1), 1e-10);
        assertEquals(0, model.scoreResidual(new int[]{2, 0, 0, 0, 0, 0}, (a, ra, b, rb) -> 0));
        EnergyMatrix corrected = new EnergyMatrix(base);
        model.applyResidualTo(corrected, 1, (a, ra, b, rb) -> 0);
        assertTrue(corrected.hasHigherOrderTerms());
        for (int i = 0; i < conf.length; i++)
            assertEquals(h[i], corrected.confE(conf[i]) - base.confE(conf[i]), 1e-10);
        // A full-data training gain must not override contradictory held-out
        // evidence. Each nested fit still sees the original training labels.
        PackStarProposalLearning.Data[] contradictory = new PackStarProposalLearning.Data[2];
        for (int f = 0; f < 2; f++) {
            double[] changed = residual.clone();
            for (int i = 0; i < conf.length; i++)
                if (PackStarProposalLearning.fold(conf[i]) == f) changed[i] = -changed[i];
            contradictory[f] = new PackStarProposalLearning.Data(conf, changed, new double[64]);
        }
        var rejected = PackStarTripleEtaCorrections.emptyFitted(6).fitSelectedSecondMoment(
                rcs, graph, data, contradictory, 1, 1, 0, 1, 0, 10, 100);
        assertTrue(rejected.positionTriples().isEmpty());
        var capped = PackStarTripleEtaCorrections.emptyFitted(6).fitSelectedSecondMoment(
                rcs, graph, data, new PackStarProposalLearning.Data[]{data, data},
                1, 1, 0, 1, 0, 10, 11); // selected table requires 3*2*2=12
        assertTrue(capped.positionTriples().isEmpty());
    }

    @Test void noSignalOrNoBudgetKeepsPairOnly() {
        RCs rcs = new RCs(new int[][]{{0, 1}, {0, 1}, {0, 1}});
        EnergyMatrix base = new EnergyMatrix(3, new int[]{2, 2, 2}, 0);
        var graph = InteractionGraph.buildFromEnergyMatrix(base, base, rcs, 0.5);
        int[][] conf = new int[8][3];
        for (int i = 0; i < 8; i++) for (int p = 0; p < 3; p++) conf[i][p] = (i >>> p) & 1;
        var data = new PackStarProposalLearning.Data(conf, new double[8], new double[8]);
        var empty = PackStarTripleEtaCorrections.emptyFitted(3).fitSelectedSecondMoment(
                rcs, graph, data, new PackStarProposalLearning.Data[]{data, data},
                1, 2, 3, 1, 1, 3, 100);
        assertTrue(empty.positionTriples().isEmpty());
    }
}
