package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar.ConfSpaces;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * Test to track MARK* behavior and answer two questions:
 * 1. Does error bound ever increase during execution?
 * 2. Are nodes that are expanded first the ones that contribute most to the final partition function?
 */
public class TestMARKStarTracking {

    public static class BoundSnapshot {
        public int iteration;
        public double epsilonBound;
        public BigDecimal upperBound;
        public BigDecimal lowerBound;
        public int numConfsEnergied;
        public int numConfsScored;

        public BoundSnapshot(int iteration, double epsilonBound, BigDecimal upperBound,
                           BigDecimal lowerBound, int numConfsEnergied, int numConfsScored) {
            this.iteration = iteration;
            this.epsilonBound = epsilonBound;
            this.upperBound = upperBound;
            this.lowerBound = lowerBound;
            this.numConfsEnergied = numConfsEnergied;
            this.numConfsScored = numConfsScored;
        }

        @Override
        public String toString() {
            return String.format("Iteration %4d: epsilon=%8.6f, bounds=[%12.6e, %12.6e], confs_energied=%4d, confs_scored=%4d",
                iteration, epsilonBound, lowerBound, upperBound, numConfsEnergied, numConfsScored);
        }
    }

    public static class NodeExpansionRecord {
        public int iteration;
        public int[] assignments;
        public double confLowerBound;
        public double confUpperBound;
        public BigDecimal subtreeUpperBound;
        public BigDecimal subtreeLowerBound;
        public BigDecimal finalContribution;  // Will be set at the end

        public NodeExpansionRecord(int iteration, int[] assignments, double confLowerBound,
                                  double confUpperBound, BigDecimal subtreeUpperBound,
                                  BigDecimal subtreeLowerBound) {
            this.iteration = iteration;
            this.assignments = assignments.clone();
            this.confLowerBound = confLowerBound;
            this.confUpperBound = confUpperBound;
            this.subtreeUpperBound = subtreeUpperBound;
            this.subtreeLowerBound = subtreeLowerBound;
        }

        public String assignmentsStr() {
            StringBuilder sb = new StringBuilder("[");
            for (int i = 0; i < assignments.length; i++) {
                if (i > 0) sb.append(",");
                sb.append(assignments[i] == -1 ? "_" : String.valueOf(assignments[i]));
            }
            sb.append("]");
            return sb.toString();
        }
    }

    @Test
    public void testTinyMARKStarWithTracking() throws IOException {
        // Create a very small conf space
        ConfSpaces confSpaces = makeSimple1GUA(2);

        // Track snapshots
        List<BoundSnapshot> snapshots = new ArrayList<>();
        List<NodeExpansionRecord> expansions = new ArrayList<>();

        // Run MARK* with custom tracking
        runMARKStarWithTracking(confSpaces, 0.68, snapshots, expansions);

        // Analysis 1: Check if error bound ever increases
        System.out.println("\n========== ERROR BOUND ANALYSIS ==========");
        boolean errorBoundIncreased = false;
        for (int i = 1; i < snapshots.size(); i++) {
            BoundSnapshot prev = snapshots.get(i-1);
            BoundSnapshot curr = snapshots.get(i);

            System.out.println(curr.toString());

            if (curr.epsilonBound > prev.epsilonBound &&
                curr.epsilonBound - prev.epsilonBound > 1e-6) {
                System.out.println("  *** ERROR BOUND INCREASED! ***");
                System.out.println("      Previous: " + prev.epsilonBound);
                System.out.println("      Current:  " + curr.epsilonBound);
                System.out.println("      Increase: " + (curr.epsilonBound - prev.epsilonBound));
                errorBoundIncreased = true;
            }
        }

        if (!errorBoundIncreased) {
            System.out.println("\n✓ Error bound NEVER increased - monotonically decreasing!");
        } else {
            System.out.println("\n✗ Error bound DID increase at some point!");
        }

        // Analysis 2: Check correlation between expansion order and final contribution
        System.out.println("\n========== NODE EXPANSION vs CONTRIBUTION ANALYSIS ==========");
        System.out.println("Note: Nodes expanded early should ideally have large contributions to final Z");

        // Write detailed log to file
        try (FileWriter writer = new FileWriter("markstar_tracking_log.txt")) {
            writer.write("MARK* Execution Tracking Log\n");
            writer.write("============================\n\n");

            writer.write("Error Bound Progression:\n");
            for (BoundSnapshot snap : snapshots) {
                writer.write(snap.toString() + "\n");
            }

            writer.write("\n\nNode Expansion Order and Contributions:\n");
            for (int i = 0; i < Math.min(50, expansions.size()); i++) {
                NodeExpansionRecord rec = expansions.get(i);
                writer.write(String.format("Expansion %4d: %s conf_bounds=[%8.3f, %8.3f] subtree_upper=%12.6e\n",
                    rec.iteration, rec.assignmentsStr(), rec.confLowerBound, rec.confUpperBound,
                    rec.subtreeUpperBound));
            }

            writer.write("\n\nSummary:\n");
            writer.write("Error bound increased? " + (errorBoundIncreased ? "YES" : "NO") + "\n");
            writer.write("Total iterations: " + snapshots.size() + "\n");
            writer.write("Total nodes expanded: " + expansions.size() + "\n");
        }

        System.out.println("\nDetailed log written to: markstar_tracking_log.txt");
    }

    private static void runMARKStarWithTracking(ConfSpaces confSpaces, double epsilon,
                                                List<BoundSnapshot> snapshots,
                                                List<NodeExpansionRecord> expansions) {
        Parallelism parallelism = Parallelism.makeCpu(2);

        // Define energy calculators
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
            .setParallelism(parallelism)
            .setIsMinimizing(true)
            .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
            .setIsMinimizing(false)
            .build();

        // Build conf energy calculator
        ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpaces.complex, minimizingEcalc)
            .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, minimizingEcalc)
                .build()
                .calcReferenceEnergies())
            .setEnergyPartition(EnergyPartition.Traditional)
            .build();

        // Calculate energy matrices
        EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
            .build()
            .calcEnergyMatrix();

        EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(
            new ConfEnergyCalculator.Builder(confSpaces.complex, rigidEcalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, rigidEcalc)
                    .build()
                    .calcReferenceEnergies())
                .build())
            .build()
            .calcEnergyMatrix();

        // Create RCs
        edu.duke.cs.osprey.astar.conf.RCs rcs = new edu.duke.cs.osprey.astar.conf.RCs(confSpaces.complex);

        // Create MARK* bound calculator with tracking
        MARKStarBoundWithTracking pfunc = new MARKStarBoundWithTracking(
            confSpaces.complex, rigidEmat, minimizingEmat, confEcalc, rcs, parallelism,
            snapshots, expansions);

        pfunc.init(epsilon);
        pfunc.compute();

        System.out.println("\nFinal Result:");
        System.out.println("Status: " + pfunc.getStatus());
        System.out.println("Values: " + pfunc.getValues());
    }

    /**
     * Extended MARKStarBound that tracks bound progression and node expansions
     */
    public static class MARKStarBoundWithTracking extends MARKStarBound {
        private List<BoundSnapshot> snapshots;
        private List<NodeExpansionRecord> expansions;
        private int iterationCount = 0;

        public MARKStarBoundWithTracking(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                        EnergyMatrix minimizingEmat, ConfEnergyCalculator confEcalc,
                                        edu.duke.cs.osprey.astar.conf.RCs rcs, Parallelism parallelism,
                                        List<BoundSnapshot> snapshots, List<NodeExpansionRecord> expansions) {
            super(confSpace, rigidEmat, minimizingEmat, confEcalc, rcs, parallelism);
            this.snapshots = snapshots;
            this.expansions = expansions;
            this.debug = true; // Enable debug output
        }

        @Override
        protected void updateBound() {
            super.updateBound();

            // Record snapshot
            snapshots.add(new BoundSnapshot(
                iterationCount++,
                epsilonBound,
                rootNode.getUpperBound(),
                rootNode.getLowerBound(),
                getNumConfsEvaluated(),
                getNumConfsScored()
            ));
        }

        // Note: We would need to override methods that expand nodes to track expansions
        // This is more complex and would require deeper hooks into the expansion logic
    }

    private static ConfSpaces makeSimple1GUA(int numFlex) {
        ConfSpaces confSpaces = new ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .addMoleculeForWildTypeRotamers(mol)
            .build();

        // Protein strand - very small
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("1", "180")
            .build();

        for (int i = 21; i < 21 + numFlex; i++) {
            protein.flexibility.get(i + "").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        // Ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("181", "215")
            .build();
        ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }
}
