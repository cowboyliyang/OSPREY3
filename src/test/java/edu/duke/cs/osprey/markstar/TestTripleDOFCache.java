package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.CachedMinimizer;
import edu.duke.cs.osprey.ematrix.PartialFixIntegration;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SubtreeDOFCache;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.jupiter.api.Test;

/**
 * Performance comparison: Triple DOF Cache
 *
 * Tests three configurations:
 * 1. Original (no cache)
 * 2. Subtree Cache only
 * 3. Subtree Cache + Triple DOF Cache
 */
public class TestTripleDOFCache {

    private static final int NUM_CPUs = 20;

    @Test
    public void testTripleCachePerformance() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("TRIPLE DOF CACHE PERFORMANCE TEST");
        System.out.println("Comparing: Original vs Subtree Cache vs Subtree+Triple Cache");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeOriginal = testScale(scale, "Original-NoCache");

        // Test 2: Subtree Cache only
        System.out.println("\n=== TEST 2: Subtree Cache Only ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;

        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeOnly = testScale(scale, "SubtreeOnly", confSpaces);
        System.out.println("\n=== Subtree Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Test 3: Subtree Cache + Triple DOF Cache
        System.out.println("\n=== TEST 3: Subtree Cache + Triple DOF Cache ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true; // NEW: Enable triple cache

        confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeAndTriple = testScale(scale, "Subtree+Triple", confSpaces);
        System.out.println("\n=== Subtree + Triple Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Summary
        printComparison(timeOriginal, timeSubtreeOnly, timeSubtreeAndTriple);
    }

    @Test
    public void testTripleCacheLargerScale() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("TRIPLE DOF CACHE PERFORMANCE TEST - LARGER SCALE");
        System.out.println("Comparing: Subtree Cache vs Subtree+Triple Cache");
        System.out.println("Configuration: 9 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 9;

        // Test 1: Subtree Cache only
        System.out.println("=== TEST 1: Subtree Cache Only (9 residues) ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;

        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeOnly = testScale(scale, "SubtreeOnly-9res", confSpaces);
        System.out.println("\n=== Subtree Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Test 2: Subtree Cache + Triple DOF Cache
        System.out.println("\n=== TEST 2: Subtree Cache + Triple DOF Cache (9 residues) ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true; // NEW: Enable triple cache

        confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeAndTriple = testScale(scale, "Subtree+Triple-9res", confSpaces);
        System.out.println("\n=== Subtree + Triple Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Summary
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE COMPARISON - 9 Residues");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s", "Subtree Cache Only", formatTime(timeSubtreeOnly)));
        System.out.println(String.format("%-35s: %s", "Subtree + Triple Cache", formatTime(timeSubtreeAndTriple)));
        System.out.println();

        double ratio = (double) timeSubtreeOnly / timeSubtreeAndTriple;
        double improvement = ((double) timeSubtreeOnly - timeSubtreeAndTriple) / timeSubtreeOnly * 100;

        if (ratio > 1.0) {
            System.out.println(String.format("Triple Cache provides %.2fx SPEEDUP!", ratio));
            System.out.println(String.format("Time reduction: %.1f%%", improvement));
        } else {
            System.out.println(String.format("Triple Cache has %.2fx overhead", 1.0/ratio));
            System.out.println(String.format("Time increase: %.1f%%", -improvement));
        }

        System.out.println();
        System.out.println("ANALYSIS:");
        System.out.println("- Triple cache stores minimized DOF values from triple corrections");
        System.out.println("- When minimizing larger conformations, triple DOFs are reused as initial guess");
        System.out.println("- Expected benefit: Better initial guess → Faster convergence");
        System.out.println("- For 9 residues: C(9,3) = 84 potential triple matches per conformation");
        System.out.println("=".repeat(100) + "\n");
    }

    private void printComparison(long timeOriginal, long timeSubtreeOnly, long timeSubtreeAndTriple) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE COMPARISON");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "Subtree Cache Only", formatTime(timeSubtreeOnly),
            (double) timeOriginal / timeSubtreeOnly));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original, %.2fx vs Subtree)",
            "Subtree + Triple Cache", formatTime(timeSubtreeAndTriple),
            (double) timeOriginal / timeSubtreeAndTriple,
            (double) timeSubtreeOnly / timeSubtreeAndTriple));
        System.out.println();

        // Analysis
        double subtreeImprovement = ((double) timeOriginal - timeSubtreeOnly) / timeOriginal * 100;
        double tripleImprovement = ((double) timeSubtreeOnly - timeSubtreeAndTriple) / timeSubtreeOnly * 100;
        double totalImprovement = ((double) timeOriginal - timeSubtreeAndTriple) / timeOriginal * 100;

        System.out.println("IMPROVEMENT ANALYSIS:");
        System.out.println(String.format("- Subtree Cache vs Original:           %+.1f%%", subtreeImprovement));
        System.out.println(String.format("- Triple Cache vs Subtree:              %+.1f%%", tripleImprovement));
        System.out.println(String.format("- Total improvement vs Original:        %+.1f%%", totalImprovement));
        System.out.println();

        System.out.println("KEY INSIGHTS:");
        System.out.println("- Subtree Cache: Reuses minimized DOF values for entire subtrees");
        System.out.println("- Triple Cache:  Provides better initial guess from cached triple DOFs");
        System.out.println("- Combined:      Two-level caching strategy for maximum benefit");

        if (tripleImprovement > 10) {
            System.out.println("- ✓✓ Triple cache provides SIGNIFICANT additional speedup!");
        } else if (tripleImprovement > 5) {
            System.out.println("- ✓ Triple cache provides measurable additional speedup");
        } else if (tripleImprovement > 0) {
            System.out.println("- Triple cache provides modest improvement");
        } else {
            System.out.println("- Triple cache has overhead on this scale (may benefit larger problems)");
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
        System.out.println("  Subtree cache:       " + (CachedMinimizer.ENABLE_SUBTREE_CACHE ? "ON" : "OFF"));
        System.out.println("  Triple cache:        " + (SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE ? "ON" : "OFF"));
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

        // Stop timer and print minimization statistics
        edu.duke.cs.osprey.markstar.MinimizationTimer.stopProgram();
        edu.duke.cs.osprey.markstar.MinimizationTimer.printStatistics();

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

        // Add mutable residues with multiple amino acid options
        // 7 mutables ALL ON PROTEIN
        if (numFlexible >= 1) {
            protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU").addWildTypeRotamers().setContinuous(); // MUTABLE 1
        }
        if (numFlexible >= 2) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "PHE").addWildTypeRotamers().setContinuous(); // MUTABLE 2
        }
        if (numFlexible >= 3) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU", "ASP").addWildTypeRotamers().setContinuous(); // MUTABLE 3
        }
        if (numFlexible >= 4) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASN", "GLN").addWildTypeRotamers().setContinuous(); // MUTABLE 4
        }
        if (numFlexible >= 5) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType, "SER", "THR").addWildTypeRotamers().setContinuous(); // MUTABLE 5
        }
        if (numFlexible >= 6) {
            protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType, "ALA", "GLY").addWildTypeRotamers().setContinuous(); // MUTABLE 6
        }
        if (numFlexible >= 7) {
            protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "LYS", "ARG").addWildTypeRotamers().setContinuous(); // MUTABLE 7
        }
        if (numFlexible >= 8) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // Flexible only
        }
        if (numFlexible >= 9) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // Flexible only
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Test minimization timing with 6 mutable + 3 flexible
     * Quick test to verify timing instrumentation works correctly
     */
    @Test
    public void testMinimizationTiming6Mutable3Flexible() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("MINIMIZATION TIMING TEST - 6 Mutable + 3 Flexible");
        System.out.println("Testing comprehensive timing instrumentation across all phases");
        System.out.println("Comparing THREE versions to verify consistency");
        System.out.println("=".repeat(100) + "\n");

        // Enable detailed minimization logging
        CachedMinimizer.ENABLE_MINIMIZATION_LOGGING = true;

        // Test 1: Original (No Cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        CachedMinimizer.clearGlobalCache();

        long timeOriginal = testScale6Mutable3Flexible("Original-NoCache", null);

        // Test 2: Subtree Cache Only
        System.out.println("\n=== TEST 2: Subtree Cache Only ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;

        TestKStar.ConfSpaces confSpaces = buildConfSpace6Mutable3Flexible();
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeOnly = testScale6Mutable3Flexible("SubtreeOnly", confSpaces);

        System.out.println("\n=== Subtree Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Test 3: Subtree + Triple Cache
        System.out.println("\n=== TEST 3: Subtree + Triple Cache ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;

        confSpaces = buildConfSpace6Mutable3Flexible();
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeAndTriple = testScale6Mutable3Flexible("Subtree+Triple", confSpaces);

        System.out.println("\n=== Subtree + Triple Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Disable logging
        CachedMinimizer.ENABLE_MINIMIZATION_LOGGING = false;

        // Print comparison
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE COMPARISON - 6 Mutable + 3 Flexible");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "Subtree Cache Only", formatTime(timeSubtreeOnly),
            (double) timeOriginal / timeSubtreeOnly));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original, %.2fx vs Subtree)",
            "Subtree + Triple Cache", formatTime(timeSubtreeAndTriple),
            (double) timeOriginal / timeSubtreeAndTriple,
            (double) timeSubtreeOnly / timeSubtreeAndTriple));
        System.out.println("=".repeat(100) + "\n");
    }

    /**
     * Quick consistency test - runs same conformations through all three versions
     * to verify they produce identical minimization values
     */
    @Test
    public void testConsistencyAcrossVersions() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("CONSISTENCY TEST - Verifying identical minimization results");
        System.out.println("Testing: Same conformations through Original, Subtree, and Subtree+Triple versions");
        System.out.println("=".repeat(100) + "\n");

        // Use smaller epsilon for faster test
        testConsistency(0.95);
    }

    private void testConsistency(double epsilon) {
        System.out.println("Building configuration space...");
        TestKStar.ConfSpaces confSpaces = buildConfSpace6Mutable3Flexible();

        // Enable detailed logging to see all minimization results
        CachedMinimizer.ENABLE_MINIMIZATION_LOGGING = true;

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
                .setEpsilon(epsilon)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .build();

        // Run all three versions on same problem
        for (int version = 1; version <= 3; version++) {
            String versionName;
            if (version == 1) {
                versionName = "Original (No Cache)";
                CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
                SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
                CachedMinimizer.clearGlobalCache();
            } else if (version == 2) {
                versionName = "Subtree Cache Only";
                CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
                SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
                CachedMinimizer.clearGlobalCache();
                CachedMinimizer.initializeGlobalCache(confSpaces.protein);
                CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
                CachedMinimizer.initializeGlobalCache(confSpaces.complex);
            } else {
                versionName = "Subtree + Triple Cache";
                CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
                SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;
                CachedMinimizer.clearGlobalCache();
                CachedMinimizer.initializeGlobalCache(confSpaces.protein);
                CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
                CachedMinimizer.initializeGlobalCache(confSpaces.complex);
            }

            System.out.println("\n" + "=".repeat(100));
            System.out.println("VERSION " + version + ": " + versionName);
            System.out.println("=".repeat(100) + "\n");

            MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                    rigidEcalc, minimizingEcalc, confEcalcFactory, settings);

            markstar.precalcEmats();
            var scores = markstar.run();

            System.out.println("\nVersion " + version + " completed.");
        }

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        CachedMinimizer.ENABLE_MINIMIZATION_LOGGING = false;

        System.out.println("\n" + "=".repeat(100));
        System.out.println("CONSISTENCY TEST COMPLETE");
        System.out.println("Review the [Minimize-*] log entries above to compare energies across versions");
        System.out.println("=".repeat(100) + "\n");
    }

    private long testScale6Mutable3Flexible(String label, TestKStar.ConfSpaces confSpaces) {
        if (confSpaces == null) {
            confSpaces = buildConfSpace6Mutable3Flexible();
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Total positions:     " + confSpaces.complex.positions.size());
        System.out.println("  Mutable positions:   6");
        System.out.println("  Flexible positions:  3");
        System.out.println("  Subtree cache:       " + (CachedMinimizer.ENABLE_SUBTREE_CACHE ? "ON" : "OFF"));
        System.out.println("  Triple cache:        " + (SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE ? "ON" : "OFF"));
        System.out.println();

        // Reset and start global minimization timer
        edu.duke.cs.osprey.markstar.MinimizationTimer.reset();
        edu.duke.cs.osprey.markstar.MinimizationTimer.startProgram();

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

        // Stop timer and print minimization statistics
        edu.duke.cs.osprey.markstar.MinimizationTimer.stopProgram();
        edu.duke.cs.osprey.markstar.MinimizationTimer.printStatistics();

        return timeMs;
    }

    private TestKStar.ConfSpaces buildConfSpace6Mutable3Flexible() {
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

        // 6 Mutable residues: Wild + 5 types each
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();

        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers().setContinuous();

        // 3 Flexible residues: Wild-type only, continuous minimization
        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Larger scale test: 7 mutable + 4 flexible positions
     * More amino acid types per mutable position (10+ types)
     * This provides a more realistic and computationally intensive test case
     */
    @Test
    public void testTripleCache7Mutable4Flexible() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("TRIPLE DOF CACHE PERFORMANCE TEST - LARGE SCALE");
        System.out.println("Comparing: Original vs Subtree Cache vs Subtree+Triple Cache");
        System.out.println("Configuration: 7 mutable + 4 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        // Test 1: Original (no cache)
        System.out.println("=== TEST 1: Original (No Cache) - 7 mutable + 4 flexible ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeOriginal = testScale7Mutable4Flexible("Original-NoCache");

        // Test 2: Subtree Cache only
        System.out.println("\n=== TEST 2: Subtree Cache Only - 7 mutable + 4 flexible ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;

        TestKStar.ConfSpaces confSpaces = buildConfSpace7Mutable4Flexible();
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeOnly = testScale7Mutable4Flexible("SubtreeOnly", confSpaces);
        System.out.println("\n=== Subtree Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Test 3: Subtree Cache + Triple DOF Cache
        System.out.println("\n=== TEST 3: Subtree Cache + Triple DOF Cache - 7 mutable + 4 flexible ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true; // NEW: Enable triple cache

        confSpaces = buildConfSpace7Mutable4Flexible();
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timeSubtreeAndTriple = testScale7Mutable4Flexible("Subtree+Triple", confSpaces);
        System.out.println("\n=== Subtree + Triple Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Summary
        printComparison7Mutable4Flexible(timeOriginal, timeSubtreeOnly, timeSubtreeAndTriple);
    }

    private long testScale7Mutable4Flexible(String label) {
        return testScale7Mutable4Flexible(label, null);
    }

    private long testScale7Mutable4Flexible(String label, TestKStar.ConfSpaces confSpaces) {
        if (confSpaces == null) {
            confSpaces = buildConfSpace7Mutable4Flexible();
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Total positions:     " + confSpaces.complex.positions.size());
        System.out.println("  Mutable positions:   7");
        System.out.println("  Flexible positions:  4");
        System.out.println("  Subtree cache:       " + (CachedMinimizer.ENABLE_SUBTREE_CACHE ? "ON" : "OFF"));
        System.out.println("  Triple cache:        " + (SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE ? "ON" : "OFF"));
        System.out.println();

        // Reset and start global minimization timer
        edu.duke.cs.osprey.markstar.MinimizationTimer.reset();
        edu.duke.cs.osprey.markstar.MinimizationTimer.startProgram();

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

        // Stop timer and print minimization statistics
        edu.duke.cs.osprey.markstar.MinimizationTimer.stopProgram();
        edu.duke.cs.osprey.markstar.MinimizationTimer.printStatistics();

        return timeMs;
    }

    /**
     * Build conf space with 7 mutable + 4 flexible residues
     * - Mutable: 10+ amino acid types per position for more realistic design space
     * - Flexible: Wild-type only with continuous minimization (.setContinuous())
     */
    private TestKStar.ConfSpaces buildConfSpace7Mutable4Flexible() {
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

        // 7 Mutable residues with expanded amino acid library
        // Each position has Wild + 10 types = 11 RCs per position
        // Common design includes: hydrophobic (ALA, VAL, ILE, LEU, PHE, TRP),
        // aromatic (TYR, PHE, TRP), polar (SER, THR), charged (ASP, GLU, LYS, ARG)

        // Protein mutable positions (4 positions)
        protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "TRP", "SER", "THR", "ASN")
            .addWildTypeRotamers().setContinuous();

        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "TRP", "SER", "THR", "GLN")
            .addWildTypeRotamers().setContinuous();

        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "TRP", "ASP", "GLU", "LYS")
            .addWildTypeRotamers().setContinuous();

        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "TRP", "SER", "THR", "ARG")
            .addWildTypeRotamers().setContinuous();

        // Ligand mutable positions (3 positions)
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "LEU", "PHE", "TYR", "TRP", "SER", "THR", "ASN", "GLN")
            .addWildTypeRotamers().setContinuous();

        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "SER", "ASP", "GLU")
            .addWildTypeRotamers().setContinuous();

        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType,
            "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "LYS", "ARG", "HIS")
            .addWildTypeRotamers().setContinuous();

        // 4 Flexible residues (wild-type only, continuous minimization)
        // C(4,3) = 4 different triple combinations
        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    private void printComparison7Mutable4Flexible(long timeOriginal, long timeSubtreeOnly, long timeSubtreeAndTriple) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE COMPARISON - 7 Mutable + 4 Flexible Residues");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "Subtree Cache Only", formatTime(timeSubtreeOnly),
            (double) timeOriginal / timeSubtreeOnly));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original, %.2fx vs Subtree)",
            "Subtree + Triple Cache", formatTime(timeSubtreeAndTriple),
            (double) timeOriginal / timeSubtreeAndTriple,
            (double) timeSubtreeOnly / timeSubtreeAndTriple));
        System.out.println();

        // Analysis
        double subtreeImprovement = ((double) timeOriginal - timeSubtreeOnly) / timeOriginal * 100;
        double tripleImprovement = ((double) timeSubtreeOnly - timeSubtreeAndTriple) / timeSubtreeOnly * 100;
        double totalImprovement = ((double) timeOriginal - timeSubtreeAndTriple) / timeOriginal * 100;

        System.out.println("IMPROVEMENT ANALYSIS:");
        System.out.println(String.format("- Subtree Cache vs Original:           %+.1f%%", subtreeImprovement));
        System.out.println(String.format("- Triple Cache vs Subtree:              %+.1f%%", tripleImprovement));
        System.out.println(String.format("- Total improvement vs Original:        %+.1f%%", totalImprovement));
        System.out.println();

        System.out.println("CONFIGURATION DETAILS:");
        System.out.println("- 7 Mutable residues:  11 amino acid types per position (Wild + 10 mutations)");
        System.out.println("- 4 Flexible residues: Wild-type only, continuous DOF minimization");
        System.out.println("- Total design space:  11^7 = 19.5 million sequence combinations");
        System.out.println("- Triple combinations: C(4,3) = 4 different triples per conformation");
        System.out.println();

        System.out.println("KEY INSIGHTS:");
        System.out.println("- Subtree Cache: Reuses minimized DOF values for entire subtrees (4+ positions)");
        System.out.println("- Triple Cache:  Provides better initial guess from cached triple DOFs");
        System.out.println("- Combined:      Two-level caching strategy for maximum benefit");
        System.out.println("- More flexible residues (4 vs 3) = more triple diversity = better cache reuse");

        if (tripleImprovement > 15) {
            System.out.println("- ✓✓✓ Triple cache provides EXCELLENT additional speedup!");
        } else if (tripleImprovement > 10) {
            System.out.println("- ✓✓ Triple cache provides SIGNIFICANT additional speedup!");
        } else if (tripleImprovement > 5) {
            System.out.println("- ✓ Triple cache provides measurable additional speedup");
        } else if (tripleImprovement > 0) {
            System.out.println("- Triple cache provides modest improvement");
        } else {
            System.out.println("- Triple cache has overhead on this configuration");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
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

    /**
     * Test that captures DOF values before and after minimization for each conformation.
     * Based on test1 from testTripleCachePerformance (Original - No Cache).
     * Saves initial and minimized DOF values to track how conformations change during minimization.
     */
    @Test
    public void testSaveDOFValuesBeforeAndAfterMinimization() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("DOF VALUE TRACKING TEST - Test 1 (No Cache)");
        System.out.println("Saving initial and minimized DOF values for each conformation");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no cache) - same configuration as first test
        System.out.println("=== TEST 1: Original (No Cache) - Tracking DOF Changes ===\n");

        // Disable ALL caches to match the true "Original" configuration
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = false; // Disable PartialFixCache (Phase 4)
        CachedMinimizer.clearGlobalCache();

        System.out.println("All caches disabled:");
        System.out.println("  - Subtree Cache: OFF");
        System.out.println("  - Triple DOF Cache: OFF");
        System.out.println("  - PartialFix Cache: OFF");
        System.out.println();

        // Enable DOF value logging to capture before/after DOF values
        // Use SimpleCCDMinimizer since that's what's actually used by the EnergyCalculator
        SimpleCCDMinimizer.ENABLE_DOF_VALUE_LOGGING = true;
        CCDMinimizer.ENABLE_DOF_VALUE_LOGGING = true; // Also enable for CCDMinimizer just in case
        EnergyCalculator.ENABLE_RCTUPLE_LOGGING = true; // Enable RC tuple logging

        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);

        System.out.println("Configuration:");
        System.out.println("  Flexible positions:  " + confSpaces.complex.positions.size());
        System.out.println("  Subtree cache:       OFF");
        System.out.println("  Triple cache:        OFF");
        System.out.println("  DOF value logging:   ON");
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

        System.out.println("Running MARKStar with DOF tracking...");
        markstar.precalcEmats();
        var scores = markstar.run();

        watch.stop();
        long timeMs = (long) watch.getTimeMs();

        System.out.println("Runtime: " + formatTime(timeMs) + "\n");

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        // Disable DOF value logging
        SimpleCCDMinimizer.ENABLE_DOF_VALUE_LOGGING = false;
        CCDMinimizer.ENABLE_DOF_VALUE_LOGGING = false;
        EnergyCalculator.ENABLE_RCTUPLE_LOGGING = false;

        // Re-enable PartialFixCache for other tests
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;

        System.out.println("\n" + "=".repeat(100));
        System.out.println("DOF VALUE TRACKING COMPLETE");
        System.out.println("Check the log output above for [DOF-BEFORE] and [DOF-AFTER] entries");
        System.out.println("Each conformation's DOF values are saved before and after minimization");
        System.out.println("=".repeat(100) + "\n");
    }
}
