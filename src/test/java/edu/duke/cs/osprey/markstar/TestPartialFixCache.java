package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.PartialFixIntegration;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
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

/**
 * Performance comparison: PartialFixCache (Phase 4)
 *
 * Tests two configurations:
 * 1. Original (no PartialFixCache)
 * 2. With PartialFixCache (BWM*-inspired L-set/M-set caching)
 */
public class TestPartialFixCache {

    private static final int NUM_CPUs = 20;

    @Test
    public void testPartialFixCachePerformance() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIALFIX CACHE PERFORMANCE TEST (Phase 4)");
        System.out.println("Comparing: Original vs PartialFixCache (BWM*-inspired L-set/M-set caching)");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no PartialFixCache)
        System.out.println("=== TEST 1: Original (No PartialFixCache) ===\n");
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = false;
        long timeOriginal = testScale(scale, "Original-NoCache");

        // Test 2: With PartialFixCache
        System.out.println("\n=== TEST 2: With PartialFixCache ===\n");
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;

        long timePartialFix = testScale(scale, "PartialFixCache");

        // Print PartialFixCache statistics
        System.out.println("\n=== PartialFixCache Detailed Statistics ===");
        PartialFixIntegration.printAllStatistics();

        // Summary
        printComparison(timeOriginal, timePartialFix);
    }

    @Test
    public void testPartialFixCacheLargerScale() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIALFIX CACHE PERFORMANCE TEST - LARGER SCALE");
        System.out.println("Comparing: Original vs PartialFixCache");
        System.out.println("Configuration: 9 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 9;

        // Test 1: Original (no PartialFixCache)
        System.out.println("=== TEST 1: Original (No PartialFixCache) - 9 residues ===\n");
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = false;
        long timeOriginal = testScale(scale, "Original-NoCache-9res");

        // Test 2: With PartialFixCache
        System.out.println("\n=== TEST 2: With PartialFixCache - 9 residues ===\n");
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;
        long timePartialFix = testScale(scale, "PartialFixCache-9res");

        // Print PartialFixCache statistics
        System.out.println("\n=== PartialFixCache Detailed Statistics (9 residues) ===");
        PartialFixIntegration.printAllStatistics();

        // Summary
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE COMPARISON - 9 Residues");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s", "With PartialFixCache", formatTime(timePartialFix)));
        System.out.println();

        double ratio = (double) timeOriginal / timePartialFix;
        double improvement = ((double) timeOriginal - timePartialFix) / timeOriginal * 100;

        if (ratio > 1.0) {
            System.out.println(String.format("PartialFixCache provides %.2fx SPEEDUP!", ratio));
            System.out.println(String.format("Time reduction: %.1f%%", improvement));
        } else {
            System.out.println(String.format("PartialFixCache has %.2fx overhead", 1.0/ratio));
            System.out.println(String.format("Time increase: %.1f%%", -improvement));
        }

        System.out.println();
        System.out.println("ANALYSIS:");
        System.out.println("- PartialFixCache uses BWM*-style branch decomposition");
        System.out.println("- L-set (internal positions): Cached and reused");
        System.out.println("- M-set (separator, 1-3 positions): Quick optimized (5 CCD iterations)");
        System.out.println("- Expected benefit: Faster upper bound tightening + skip more minimizations");
        System.out.println("- For 9 residues: More complex tree → More caching opportunities");
        System.out.println("=".repeat(100) + "\n");
    }

    private void printComparison(long timeOriginal, long timePartialFix) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE COMPARISON");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No PartialFixCache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "With PartialFixCache", formatTime(timePartialFix),
            (double) timeOriginal / timePartialFix));
        System.out.println();

        // Analysis
        double improvement = ((double) timeOriginal - timePartialFix) / timeOriginal * 100;

        System.out.println("IMPROVEMENT ANALYSIS:");
        System.out.println(String.format("- PartialFixCache vs Original:  %+.1f%%", improvement));
        System.out.println();

        System.out.println("KEY INSIGHTS:");
        System.out.println("- PartialFixCache: BWM*-inspired L-set/M-set separation strategy");
        System.out.println("- L-set caching:   Reuses minimized DOF values for internal positions");
        System.out.println("- M-set quick opt: 5 CCD iterations for separator (1-3 positions)");
        System.out.println("- Skip decisions:  3 conditions (lower tight, upper tight, gap tight)");
        System.out.println("- Upper bound:     Tightened symmetrically to triple correction");

        if (improvement > 30) {
            System.out.println("- ✓✓✓ PartialFixCache provides EXCELLENT speedup!");
        } else if (improvement > 20) {
            System.out.println("- ✓✓ PartialFixCache provides SIGNIFICANT speedup!");
        } else if (improvement > 10) {
            System.out.println("- ✓ PartialFixCache provides measurable speedup");
        } else if (improvement > 0) {
            System.out.println("- PartialFixCache provides modest improvement");
        } else {
            System.out.println("- PartialFixCache has overhead (may benefit larger problems)");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    private long testScale(int numFlexible, String label) {
        return testScale(numFlexible, label, null);
    }

    private long testScale(int numFlexible, String label, TestKStar.ConfSpaces confSpaces) {
        if (confSpaces == null) {
            confSpaces = buildConfSpace(numFlexible);
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions:  " + confSpaces.complex.positions.size());
        System.out.println("  PartialFixCache:     " + (PartialFixIntegration.ENABLE_PARTIALFIX_CACHE ? "ON" : "OFF"));
        System.out.println();

        Stopwatch watch = new Stopwatch().start();
        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
            return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
        };

        MARKStar.Settings settings = new MARKStar.Settings.Builder()
                .setEpsilon(0.99)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        System.out.println("Running...");
        markstar.precalcEmats();
        var scores = markstar.run();

        watch.stop();
        long timeMs = (long) watch.getTimeMs();

        System.out.println("Runtime: " + formatTime(timeMs) + "\n");

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        return timeMs;
    }

    private TestKStar.ConfSpaces buildConfSpace(int numFlexible) {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();

        // Add flexible residues
        if (numFlexible >= 1) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 2) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 3) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 4) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 5) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 6) {
            ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 7) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 8) {
            ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 9) {
            protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Test with skip decision statistics
     * Shows how many minimizations were skipped and for what reasons
     */
    @Test
    public void testPartialFixCacheSkipStatistics() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIALFIX CACHE SKIP STATISTICS TEST");
        System.out.println("Shows skip decisions: triple_correction_lower, partialfix_upper, tight_gap");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        System.out.println("=== With PartialFixCache - Detailed Skip Statistics ===\n");
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;

        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);

        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(false)
                .build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
            return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
        };

        MARKStar.Settings settings = new MARKStar.Settings.Builder()
                .setEpsilon(0.99)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .build();

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

        System.out.println("Running MARK* with PartialFixCache...");
        System.out.println("Watch for [MINIMIZATION_DECISION] and [SKIP_MINIMIZATION] log entries\n");

        markstar.precalcEmats();
        var scores = markstar.run();

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        System.out.println("\n" + "=".repeat(100));
        System.out.println("SKIP STATISTICS SUMMARY");
        System.out.println("=".repeat(100));
        System.out.println("Review the log output above to see:");
        System.out.println("- [MINIMIZATION_DECISION]: Shows decision for each conformation");
        System.out.println("  * reason=triple_correction_lower: Skipped due to triple correction");
        System.out.println("  * reason=partialfix_upper: Skipped due to PartialFixCache tightening upper");
        System.out.println("  * reason=tight_gap: Skipped due to tight bounds gap");
        System.out.println("  * reason=both_bounds_tightened: Both lower and upper tightened");
        System.out.println("  * reason=no_sufficient_correction: Full minimization executed");
        System.out.println();
        System.out.println("- [SKIP_MINIMIZATION]: Details of skip (lower/upper improvements)");
        System.out.println("- [PARTIALFIX_UPPER]: PartialFixCache cache hit/miss statistics");
        System.out.println("=".repeat(100) + "\n");
    }

    private String formatTime(long ms) {
        if (ms < 1000) {
            return ms + " ms";
        } else if (ms < 60000) {
            return String.format("%.2f s", ms / 1000.0);
        } else {
            long minutes = ms / 60000;
            long seconds = (ms % 60000) / 1000;
            return String.format("%d min %d s", minutes, seconds);
        }
    }
}
