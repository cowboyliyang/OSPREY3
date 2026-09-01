package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** Exact agreement between signed selected triples and EnergyMatrix. */
public class TestPackStarTripleEtaCorrections {

    @Test
    public void assignmentCapFailureIsTypedForPairOnlyFallback() {
        PackStarTripleEtaCorrections.AssignmentCapExceededException error =
                assertThrows(
                        PackStarTripleEtaCorrections.AssignmentCapExceededException.class,
                        () -> PackStarTripleEtaCorrections.ensureAssignmentCapacity(
                                500001L, 500000L));

        assertEquals(500001L, error.assignments);
        assertEquals(500000L, error.maximumAssignments);
        assertTrue(error.getMessage().contains("assignments=500001"));
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

    private static EnergyMatrix makeTriangleBase(int[] cards) {
        EnergyMatrix emat = new EnergyMatrix(cards.length, cards, 0.0);
        int[][] edges = {{0, 1}, {0, 2}, {1, 2}};
        for (int[] edge : edges) {
            for (int rc1 = 0; rc1 < cards[edge[0]]; rc1++) {
                for (int rc2 = 0; rc2 < cards[edge[1]]; rc2++) {
                    emat.setPairwise(edge[0], rc1, edge[1], rc2, 1.0);
                }
            }
        }
        return emat;
    }

    @Test
    public void fittedResidualIsAddedDirectlyWithoutPairAllocation() {
        int[] cards = {2, 2, 2, 2};
        RCs rcs = makeRCs(cards);
        EnergyMatrix base = makeTriangleBase(cards);
        InteractionGraph graph = InteractionGraph.buildFromEnergyMatrix(
                base, base, rcs, 0.5);

        assertEquals(1L,
                PackStarTripleEtaCorrections.countCliquePositionTriples(
                        cards.length, graph));
        assertEquals(8L,
                PackStarTripleEtaCorrections.countCliqueAssignments(
                        rcs, graph));

        List<PackStarTripleEtaCorrections.Entry> entries = new ArrayList<>();
        for (int rc2 = 0; rc2 < cards[2]; rc2++) {
            for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                for (int rc0 = 0; rc0 < cards[0]; rc0++) {
                    PackStarTripleEtaCorrections.Entry entry =
                            new PackStarTripleEtaCorrections.Entry(
                                    entries.size(),
                                    2, rc2, 1, rc1, 0, rc0, 0.0);
                    entry.complete(rc0 + 2.0 * rc1 + 4.0 * rc2 - 1.0);
                    entries.add(entry);
                }
            }
        }
        PackStarTripleEtaCorrections corrections =
                new PackStarTripleEtaCorrections(
                        cards.length, 0.5, 1L, entries, 0.0);

        assertEquals(6L, corrections.positiveFactors);
        assertEquals(2L, corrections.nonpositiveFactors);
        assertEquals(0L, corrections.nonfiniteFactors);

        PackStarTripleEtaCorrections.PairEtaLookup pairEta =
                (pos1, rc1, pos2, rc2) -> 2.0 / 3.0;
        PackStarTripleEtaCorrections.ResidualSummary summary =
                corrections.summarizeResidual(pairEta);
        assertEquals(6L, summary.positiveFactors);
        assertEquals(1L, summary.negativeFactors);
        assertEquals(1L, summary.zeroFactors);
        assertEquals(-1.0, summary.minimumCorrectionKcal, 1.0e-12);
        assertEquals(6.0, summary.maximumCorrectionKcal, 1.0e-12);
        assertEquals(6.0,
                summary.maximumAbsoluteCorrectionKcal, 1.0e-12);

        double scale = 0.4;
        EnergyMatrix corrected = new EnergyMatrix(base);
        corrections.applyResidualTo(corrected, scale, pairEta);
        assertTrue(corrected.hasHigherOrderTerms());

        for (int rc3 = 0; rc3 < cards[3]; rc3++) {
            for (int rc2 = 0; rc2 < cards[2]; rc2++) {
                for (int rc1 = 0; rc1 < cards[1]; rc1++) {
                    for (int rc0 = 0; rc0 < cards[0]; rc0++) {
                        int[] conf = {rc0, rc1, rc2, rc3};
                        double jointScore =
                                corrections.scoreJointCorrection(conf);
                        double residualScore =
                                corrections.scoreResidual(conf, pairEta);
                        assertEquals(rc0 + 2.0 * rc1
                                        + 4.0 * rc2 - 1.0,
                                jointScore, 1.0e-12);
                        // Pair eta is already outside this fitted residual.
                        // The lookup must not be allocated or subtracted.
                        assertEquals(jointScore,
                                residualScore, 1.0e-12);
                        double matrixDelta = corrected.getInternalEnergy(
                                new RCTuple(conf))
                                - base.getInternalEnergy(new RCTuple(conf));
                        assertEquals(scale * residualScore,
                                matrixDelta, 1.0e-12,
                                "conf=" + java.util.Arrays.toString(conf));
                    }
                }
            }
        }
    }

    @Test
    public void residualFitSelectsAStableNonCliqueTripleAndDeclaresFillEdges() {
        int[] cards = {2, 2, 2, 2};
        RCs rcs = makeRCs(cards);
        EnergyMatrix base = makeTriangleBase(cards);
        InteractionGraph graph = InteractionGraph.buildFromEnergyMatrix(
                base, base, rcs, 0.5);
        PackStarTripleEtaCorrections localPrior =
                new PackStarTripleEtaCorrections(
                        cards.length, 0.5, 0L,
                        List.of(), 0.0);

        int[][] confs = new int[16][4];
        double[] residuals = new double[16];
        int index = 0;
        for (int rc3 = 0; rc3 < 2; rc3++) {
            for (int rc2 = 0; rc2 < 2; rc2++) {
                for (int rc1 = 0; rc1 < 2; rc1++) {
                    for (int rc0 = 0; rc0 < 2; rc0++) {
                        confs[index] = new int[]{rc0, rc1, rc2, rc3};
                        residuals[index] = ((rc0 ^ rc1 ^ rc3) == 0)
                                ? -2.0 : 2.0;
                        index++;
                    }
                }
            }
        }

        PackStarTripleEtaCorrections fitted =
                localPrior.fitSelectedResidual(
                        rcs, graph, confs, residuals,
                        1, 3, 2,
                        0.0, 0.0, 3.0, 100L);
        assertEquals(1, fitted.positionTriples().size());
        assertEquals("0,1,3",
                fitted.positionTriples().get(0).toString());
        assertEquals(8L, fitted.factorAssignments);
        assertEquals(2, fitted.requiredFillEdges(graph).size());

        PackStarTripleEtaCorrections.PairEtaLookup ignoredPairEta =
                (pos1, rc1, pos2, rc2) -> 1000.0;
        for (int sample = 0; sample < confs.length; sample++) {
            assertEquals(residuals[sample],
                    fitted.scoreResidual(
                            confs[sample], ignoredPairEta),
                    1.0e-12);
        }
    }
}
