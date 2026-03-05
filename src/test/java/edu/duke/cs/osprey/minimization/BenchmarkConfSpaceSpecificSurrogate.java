package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.ConfSpaceSpecificSurrogateFactory;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.markstar.framework.ApproximatedGridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.branch.BranchDecomposition;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeNode;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.io.File;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Manual benchmark/demo for task-specific surrogate workflow.
 * Disabled by default because this is a heavyweight benchmark.
 */
@Disabled("Manual benchmark: enable explicitly when validating task-specific surrogate flow")
public class BenchmarkConfSpaceSpecificSurrogate {

    @Test
    public void demoTaskSpecificSurrogateOnGridDPAndCCD() {

        ForcefieldParams ffparams = new ForcefieldParams();
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
        strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
        strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism).build();

        SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, minimizingEcalc)
                .build().calcReferenceEnergies();
        ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, minimizingEcalc)
                .setReferenceEnergies(eref).build();

        EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                .build().calcEnergyMatrix();
        EnergyMatrix ematRigid = ematMinimized;

        RCs rcs = new RCs(confSpace);
        InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                confSpace, ematRigid, ematMinimized, rcs, 8.0, 0.1
        );
        BranchDecomposition bd = new BranchDecomposition(ig);
        bd.compute();
        RootedTreeNode rootedRoot = bd.rootBranchTree(rcs);
        RootedTreeEdge.postOrderCompLlambda(rootedRoot);
        RootedTreeEdge rootEdge = rootedRoot.getLeftChild().getChildOfEdge();
        rootEdge.compactTree();

        ConfAStarTree astar = new ConfAStarTree.Builder(ematMinimized, rcs)
                .setTraditional().build();
        ConfSearch.ScoredConf scoredConf = astar.nextConf();
        int[] conf = scoredConf.getAssignments();
        RCTuple tuple = new RCTuple(conf);

        File cacheRoot = new File("/tmp/osprey-confspace-surrogate-cache");
        ConfEnergyCalculator surrogateConfEcalc = ConfSpaceSpecificSurrogateFactory.withTaskSpecificApproximator(
                confEcalc, cacheRoot, "demo-1cc8-3res", 10, 1e6
        );

        ApproximatedGridDPMinimizer gridDP = new ApproximatedGridDPMinimizer(
                confSpace, ig, rootEdge, 2, ffparams, eref, true, surrogateConfEcalc.amat, true
        );
        var gridResult = gridDP.minimize(conf);

        var approxObjective = ConfSpaceSpecificSurrogateFactory.makeApproximationObjective(
                surrogateConfEcalc, tuple, false
        );
        SimpleCCDMinimizer ccd = new SimpleCCDMinimizer((context) -> new QuadraticLineSearcher());
        ccd.init(approxObjective);
        ccd.setMaxIterations(3);

        DoubleMatrix1D startx = DoubleFactory1D.dense.make(approxObjective.getNumDOFs());
        if (gridResult.bestDOFValues != null) {
            for (int d = 0; d < startx.size(); d++) {
                startx.set(d, gridResult.bestDOFValues[d]);
            }
        } else {
            DoubleMatrix1D lo = approxObjective.getConstraints()[0];
            DoubleMatrix1D hi = approxObjective.getConstraints()[1];
            for (int d = 0; d < startx.size(); d++) {
                startx.set(d, (lo.get(d) + hi.get(d)) / 2.0);
            }
        }

        Minimizer.Result approxCCDResult = ccd.minimizeFrom(startx);

        System.out.println(String.format(
                "Task-specific surrogate demo: GridDP=%.6f  CCD(surr)=%.6f",
                gridResult.energy, approxCCDResult.energy
        ));
        assertTrue(Double.isFinite(gridResult.energy));
        assertTrue(Double.isFinite(approxCCDResult.energy));

        ccd.close();
        minimizingEcalc.close();
    }
}
