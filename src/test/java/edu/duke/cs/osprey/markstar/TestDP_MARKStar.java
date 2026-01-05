package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Test DP-enhanced MARK* on real proteins
 * Compares DP correction selection vs original greedy approach
 */
public class TestDP_MARKStar {

    private static final int NUM_CPUs = 4;
    private static final boolean REDUCE_MINIMIZATIONS = true;

    @Test
    public void test2RL0_GreedyVsDP() {
        System.out.println("\n================================================================");
        System.out.println(" DP-MARKStar vs Original MARK* on 2RL0 Protein");
        System.out.println("================================================================\n");

        // Create conf space
        TestKStar.ConfSpaces confSpaces = make2RL0();

        // Test with greedy first
        System.out.println("=== Running MARK* with GREEDY correction selection ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        TestMARKStar.Result greedyResult = runMARKStarWithDP(confSpaces, 0.68, "GREEDY");

        System.out.println("\n\n=== Running MARK* with DP correction selection ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        TestMARKStar.Result dpResult = runMARKStarWithDP(confSpaces, 0.68, "DP");

        // Print comparison
        printComparison(greedyResult, dpResult);
    }

    private TestKStar.ConfSpaces make2RL0() {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

        // Configure forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

        // Template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .build();

        // Define protein strand - LARGE: 5 flexible residues
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Define ligand strand - LARGE: 4 flexible residues
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Make conf spaces
        confSpaces.protein = new SimpleConfSpace.Builder()
            .addStrand(protein)
            .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
            .addStrand(ligand)
            .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
            .addStrands(protein, ligand)
            .build();

        return confSpaces;
    }

    public static TestMARKStar.Result runMARKStarWithDP(TestKStar.ConfSpaces confSpaces, double epsilon, String label) {

        Stopwatch totalWatch = new Stopwatch().start();

        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        // Define the minimizing energy calculator
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        // Define the rigid energy calculator
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        // Conf energy calculator factory
        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
            return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
        };

        TestMARKStar.Result result = new TestMARKStar.Result();

        MARKStar.Settings settings = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setEnergyMatrixCachePattern("*.dp.mark.emat")
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .setReduceMinimizations(REDUCE_MINIMIZATIONS)
                .build();

        result.markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        System.out.println("Pre-calculating energy matrices...");
        result.markstar.precalcEmats();

        System.out.println("Running MARK* search (" + label + ")...");
        result.scores = result.markstar.run();

        totalWatch.stop();

        System.out.println("\n=== " + label + " Results ===");
        System.out.println("Total runtime: " + totalWatch.getTime(2));
        if (result.scores != null && result.scores.size() > 0) {
            MARKStar.ScoredSequence best = result.scores.get(0);
            System.out.println("Best sequence: " + best.sequence);
            System.out.println("K* score: " + best.score);
        }

        return result;
    }

    private void printComparison(TestMARKStar.Result greedy, TestMARKStar.Result dp) {
        System.out.println("\n================================================================");
        System.out.println(" FINAL COMPARISON: Greedy vs DP on 2RL0 Protein");
        System.out.println("================================================================\n");

        if (greedy.scores != null && greedy.scores.size() > 0 &&
            dp.scores != null && dp.scores.size() > 0) {

            MARKStar.ScoredSequence greedyBest = greedy.scores.get(0);
            MARKStar.ScoredSequence dpBest = dp.scores.get(0);

            System.out.println("Best Sequences:");
            System.out.println("  Greedy: " + greedyBest.sequence);
            System.out.println("  DP:     " + dpBest.sequence);
            System.out.println("  Match:  " + greedyBest.sequence.toString().equals(dpBest.sequence.toString()));

            System.out.println("\nK* Scores:");
            System.out.println("  Greedy: " + greedyBest.score);
            System.out.println("  DP:     " + dpBest.score);

            // Performance metrics
            System.out.println("\nPerformance Metrics:");
            System.out.println("  Protein: 2RL0 (9 flexible residues)");
            System.out.println("  Protein flex: 5 residues (G649 with 4 mutations, G650-653)");
            System.out.println("  Ligand flex: 4 residues (A172, A192-194)");
            System.out.println("  Epsilon: 0.68");

            System.out.println("\n================================================================");
            System.out.println("CONCLUSION:");
            System.out.println("DP algorithm successfully ran on real protein (2RL0)");
            System.out.println("Ready for larger-scale testing with more memory");
            System.out.println("================================================================");
        }
    }
}
