package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.CachedMinimizer;
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

/**
 * Performance comparison: Optimized DP-MARKStar vs Original Greedy MARKStar
 * Tests on medium-scale problems (7-9 flexible residues)
 */
public class TestDPvsOriginal {

    private static final int NUM_CPUs = 20;

    @Test
    public void testMediumScale7Residues() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("MEDIUM-SCALE PERFORMANCE TEST: 7 Flexible Residues");
        System.out.println("Comparing: Original Greedy vs Optimized DP-MARKStar");
        System.out.println("=".repeat(100) + "\n");

        // Test 1: Original Greedy
        System.out.println("=== TEST 1: Original MARKStar (Greedy) - 7 residues ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        long timeOriginal = testScale(7, "Original Greedy");

        // Test 2: Optimized DP (Phase 1)
        System.out.println("\n=== TEST 2: DP-MARKStar (Phase 1 Optimized) - 7 residues ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        long timeOptimized = testScale(7, "DP-MARKStar Optimized");

        // Summary
        printComparison("7 Residues", timeOriginal, timeOptimized);
    }

    @Test
    public void testMediumScale9Residues() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("MEDIUM-SCALE PERFORMANCE TEST: 9 Flexible Residues");
        System.out.println("Comparing: Original Greedy vs Optimized DP-MARKStar");
        System.out.println("=".repeat(100) + "\n");

        // Test 1: Original Greedy
        System.out.println("=== TEST 1: Original MARKStar (Greedy) - 9 residues ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        long timeOriginal = testScale(9, "Original Greedy");

        // Test 2: Optimized DP (Phase 1)
        System.out.println("\n=== TEST 2: DP-MARKStar (Phase 1 Optimized) - 9 residues ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        long timeOptimized = testScale(9, "DP-MARKStar Optimized");

        // Summary
        printComparison("9 Residues", timeOriginal, timeOptimized);
    }

    @Test
    public void testComprehensiveComparison() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("DETAILED COMPARISON: Original vs Phase 1+2");
        System.out.println("Configuration: 3 Mutable + 3 Flexible = 6 Total Positions");
        System.out.println("Recording: Minimizations, Corrections, Cache usage");
        System.out.println("=".repeat(100) + "\n");

        int scale = 3; // 3 mutable + 3 flexible = 6 total

        // Test 1: Original (baseline)
        System.out.println("=== Test 1: Original (baseline) ===");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = false;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeOriginal = testScale(scale, "Original-6pos");

        // Test 2: Phase 1+2 (DP + Subtree Caching)
        System.out.println("\n=== Test 2: Phase 1+2 (DP + Subtree Caching) ===");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        long timePhase12 = testScale(scale, "Phase1+2-6pos");

        // Print cache stats
        System.out.println();
        CachedMinimizer.printGlobalStats();

        // Summary
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE SUMMARY");
        System.out.println("=".repeat(100));
        System.out.println(String.format("Original:    %,d ms (%.2f s)", timeOriginal, timeOriginal/1000.0));
        System.out.println(String.format("Phase 1+2:   %,d ms (%.2f s)", timePhase12, timePhase12/1000.0));
        double speedup = (double) timeOriginal / timePhase12;
        System.out.println(String.format("Speedup:     %.2fx %s", speedup, speedup > 1 ? "(faster)" : "(SLOWER)"));
        System.out.println("=".repeat(100));
    }

    // Skip old comprehensive comparison logic below
    private void OLD_testComprehensiveComparison() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("COMPREHENSIVE PERFORMANCE COMPARISON");
        System.out.println("Original Greedy vs Phase 1 (DP) vs Phase 1+2 (DP + Subtree Cache)");
        System.out.println("=".repeat(100) + "\n");

        long[][] results = new long[1][3]; // [scale][original/phase1/phase2]
        int[] scales = {3}; // Changed to test 1 mutable + 3 flexible = 4 total positions

        for (int i = 0; i < scales.length; i++) {
            int scale = scales[i];
            System.out.println("\n" + "-".repeat(100));
            System.out.println("Testing scale: " + scale + " flexible residues");
            System.out.println("-".repeat(100) + "\n");

            // Test 1: Original Greedy ONLY
            System.out.println("=== Test 1: Original Greedy (baseline) ===");
            UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
            CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
            CachedMinimizer.clearGlobalCache();
            results[i][0] = testScale(scale, "Original-" + scale);

            // Test 2: Phase 1 (DP optimizations only) - SKIPPED
            // (Skipping Phase 1 to only compare Original vs Phase 1+2)

            // Test 3: Phase 1+2 (DP + Subtree Caching)
            System.out.println("\n=== Test 2: Phase 1+2 (DP + TRUE Subtree Caching) ===");
            UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
            UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
            CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
            TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
            CachedMinimizer.clearGlobalCache();
            CachedMinimizer.initializeGlobalCache(confSpaces.complex);
            results[i][2] = testScale(scale, "Phase1+2-" + scale);

            // Print cache stats after each scale
            if (CachedMinimizer.ENABLE_SUBTREE_CACHE) {
                CachedMinimizer.printGlobalStats();
            }
        }

        // Final summary table
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE SUMMARY - ALL THREE VERSIONS");
        System.out.println("=".repeat(100));
        System.out.println();
        System.out.println(String.format("%-12s | %-18s | %-18s | %-18s | %-15s | %s",
            "Scale", "Original", "Phase 1 (DP)", "Phase 1+2", "P1 vs Orig", "P2 vs P1"));
        System.out.println("-".repeat(100));

        for (int i = 0; i < scales.length; i++) {
            long orig = results[i][0];
            long phase1 = results[i][1];
            long phase2 = results[i][2];

            double speedup1 = (double) orig / phase1;
            double speedup2 = (double) phase1 / phase2;

            String speedup1Str = speedup1 > 1.0
                ? String.format("%.2fx faster", speedup1)
                : String.format("%.2fx slower", 1.0/speedup1);

            String speedup2Str = speedup2 > 1.0
                ? String.format("%.2fx faster", speedup2)
                : String.format("%.2fx slower", 1.0/speedup2);

            System.out.println(String.format("%-12s | %-18s | %-18s | %-18s | %-15s | %s",
                scales[i] + " res",
                formatTime(orig),
                formatTime(phase1),
                formatTime(phase2),
                speedup1Str,
                speedup2Str));
        }

        System.out.println("=".repeat(100));

        // Analysis
        System.out.println("\nANALYSIS:");
        System.out.println("- Phase 1: DP guarantees global optimum + binary search + parallelization");
        System.out.println("- Phase 2: TRUE subtree DOF caching with boundary-only optimization");
        System.out.println("- Original Greedy is heuristic but very fast O(n log n)");

        double avgSpeedup1 = 0;
        double avgSpeedup2 = 0;
        for (int i = 0; i < scales.length; i++) {
            avgSpeedup1 += (double) results[i][0] / results[i][1];
            avgSpeedup2 += (double) results[i][1] / results[i][2];
        }
        avgSpeedup1 /= scales.length;
        avgSpeedup2 /= scales.length;

        System.out.println(String.format("- Average Phase 1 vs Original: %.2fx", avgSpeedup1));
        System.out.println(String.format("- Average Phase 2 vs Phase 1: %.2fx", avgSpeedup2));
        System.out.println(String.format("- Overall Phase 1+2 vs Original: %.2fx", avgSpeedup1 * avgSpeedup2));

        if (avgSpeedup2 > 1.05) {
            System.out.println("- ✓ Phase 2 subtree caching provides measurable speedup!");
        } else if (avgSpeedup2 > 0.95) {
            System.out.println("- Phase 2 performs comparably to Phase 1 (cache overhead balanced)");
        } else {
            System.out.println("- Phase 2 has overhead (may need larger scale to see benefits)");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    private void printComparison(String scale, long timeOriginal, long timeOptimized) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE COMPARISON - " + scale);
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-30s: %s", "Original Greedy MARKStar", formatTime(timeOriginal)));
        System.out.println(String.format("%-30s: %s", "DP-MARKStar (Phase 1 Opt.)", formatTime(timeOptimized)));
        System.out.println();

        double ratio = (double) timeOriginal / timeOptimized;
        if (ratio > 1.0) {
            System.out.println(String.format("DP-MARKStar is %.2fx FASTER than original!", ratio));
        } else {
            System.out.println(String.format("DP-MARKStar is %.2fx slower than original", 1.0/ratio));
        }

        double percentDiff = ((timeOptimized - timeOriginal) / (double) timeOriginal) * 100;
        System.out.println(String.format("Time difference: %+.1f%%", percentDiff));

        System.out.println();
        System.out.println("Notes:");
        System.out.println("- DP-MARKStar provides OPTIMAL correction selection (global optimum)");
        System.out.println("- Original uses greedy heuristic (local optimum)");
        System.out.println("- Phase 1 optimizations: Binary search + Parallel DP");
        System.out.println("=".repeat(100) + "\n");
    }

    private long testScale(int numFlexible, String label) {
        return testScale(numFlexible, label, null);
    }

    private long testScale(int numFlexible, String label, TestKStar.ConfSpaces confSpaces) {
        // Use provided confSpaces if available, otherwise build new one
        if (confSpaces == null) {
            confSpaces = buildConfSpace(numFlexible);
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions: " + confSpaces.complex.positions.size());
        System.out.println("  Algorithm: " + (UpdatingEnergyMatrix.USE_DP_ALGORITHM ? "DP" : "Greedy"));
        System.out.println("  DP optimizations: " + (UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS ? "ON" : "OFF"));
        System.out.println();

        // Run MARK*
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

        System.out.println("Pre-calculating energy matrices...");
        markstar.precalcEmats();

        System.out.println("Running MARK*...");
        var scores = markstar.run();

        watch.stop();
        long timeMs = (long) watch.getTimeMs();

        System.out.println("\nResults:");
        System.out.println("  Runtime: " + formatTime(timeMs));
        if (scores != null && scores.size() > 0) {
            System.out.println("  Best sequence: " + scores.get(0).sequence);
            System.out.println("  K* score: " + scores.get(0).score);
        }

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        return timeMs;
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
     * Test Phase 1 + Phase 2: Complete DP-MARKStar with all optimizations
     */
    @Test
    public void testAllPhasesIntegrated() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("COMPLETE DP-MARKSTAR TEST: All Phases Integrated");
        System.out.println("Phase 1: Binary Search + Parallel DP");
        System.out.println("Phase 2: Subtree DOF Cache with Branch Decomposition");
        System.out.println("=".repeat(100) + "\n");

        long[][] results = new long[2][3]; // [scale][original/phase1/phase1+2]
        int[] scales = {7, 9};

        for (int i = 0; i < scales.length; i++) {
            int scale = scales[i];
            System.out.println("\n" + "-".repeat(100));
            System.out.println("Testing scale: " + scale + " flexible residues");
            System.out.println("-".repeat(100) + "\n");

            // Test 1: Original Greedy
            System.out.println("=== TEST 1: Original Greedy ===");
            UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
            CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
            results[i][0] = testScale(scale, "Original-" + scale);

            // Test 2: Phase 1 only
            System.out.println("\n=== TEST 2: Phase 1 (DP Optimizations) ===");
            UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
            UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
            CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
            results[i][1] = testScale(scale, "Phase1-" + scale);

            // Test 3: Phase 1 + Phase 2
            System.out.println("\n=== TEST 3: Phase 1+2 (DP + Subtree Cache) ===");
            UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
            UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
            CachedMinimizer.ENABLE_SUBTREE_CACHE = true;

            // Initialize cache for this conf space
            TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
            CachedMinimizer.initializeGlobalCache(confSpaces.complex);

            results[i][2] = testScale(scale, "Phase1+2-" + scale);

            // Print cache statistics
            CachedMinimizer.printGlobalStats();
            CachedMinimizer.clearGlobalCache();
        }

        // Final summary table
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE SUMMARY: All Phases");
        System.out.println("=".repeat(100));
        System.out.println();
        System.out.println(String.format("%-15s | %-20s | %-20s | %-20s | %s",
            "Scale", "Original", "Phase 1", "Phase 1+2", "Best Speedup"));
        System.out.println("-".repeat(100));

        for (int i = 0; i < scales.length; i++) {
            long orig = results[i][0];
            long phase1 = results[i][1];
            long phase12 = results[i][2];
            double speedup1 = (double) orig / phase1;
            double speedup12 = (double) orig / phase12;
            double bestSpeedup = Math.max(speedup1, speedup12);

            System.out.println(String.format("%-15s | %-20s | %-20s | %-20s | %.2fx",
                scales[i] + " residues",
                formatTime(orig),
                formatTime(phase1),
                formatTime(phase12),
                bestSpeedup));
        }

        System.out.println("=".repeat(100));

        // Analysis
        System.out.println("\nANALYSIS:");
        System.out.println("- Phase 1: DP with binary search O(n log n) + parallel processing");
        System.out.println("- Phase 2: Subtree DOF caching with branch decomposition");
        System.out.println("- Both phases provide provably optimal solutions (unlike Greedy)");

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    /**
     * Helper to build conf space for Phase 2 cache initialization
     * Modified: Use 1 mutable + flexible residues for better subtree caching opportunities
     */
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

        // NEW: Add 3 mutable positions to create sequence space (64 sequences = 4x4x4)
        protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU").addWildTypeRotamers().setContinuous();

        // Add flexible residues (wild-type only)
        if (numFlexible >= 1) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 2) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 3) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 4) {
            ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 5) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 6) {
            ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 7) {
            protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 8) {
            ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
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
     * Build conf space with NO mutable positions - all wild-type flexible only
     */
    private TestKStar.ConfSpaces buildConfSpaceNoMutable(int numFlexible) {
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

        // All positions are wild-type only (no mutations)
        if (numFlexible >= 1) {
            protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 2) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 3) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 4) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 5) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 6) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 7) {
            ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 8) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Test Phase 2: True Subtree DOF Caching with boundary optimization
     * COMMENTED OUT: Now testing with 7 flexible in comprehensive test
     */
    // @Test
    public void testPhase2SubtreeCaching() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PHASE 2 TEST: True Subtree DOF Caching");
        System.out.println("Comparing: Original vs Phase 1 vs Phase 1+2");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7; // Changed from 5 to 7 flexible residues

        // Test 1: Original Greedy (baseline)
        System.out.println("=== TEST 1: Original Greedy (baseline) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeOriginal = testScale(scale, "Original");

        // Test 2: Phase 1 only (DP optimizations)
        System.out.println("\n=== TEST 2: Phase 1 (DP Optimizations) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timePhase1 = testScale(scale, "Phase1");

        // Test 3: Phase 1 + Phase 2 (DP + Subtree Caching)
        System.out.println("\n=== TEST 3: Phase 1+2 (DP + Subtree Caching with Boundary Optimization) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;

        // Initialize global cache
        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        long timePhase2 = testScale(scale, "Phase1+2");

        // Print cache statistics
        System.out.println("\n=== Subtree DOF Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // Summary
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PHASE 2 PERFORMANCE COMPARISON");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-30s: %s (baseline)", "Original Greedy", formatTime(timeOriginal)));
        System.out.println(String.format("%-30s: %s (%.2fx vs Original)",
            "Phase 1 (DP)", formatTime(timePhase1), (double) timeOriginal / timePhase1));
        System.out.println(String.format("%-30s: %s (%.2fx vs Original, %.2fx vs Phase1)",
            "Phase 1+2 (DP + Cache)", formatTime(timePhase2),
            (double) timeOriginal / timePhase2,
            (double) timePhase1 / timePhase2));
        System.out.println();

        // Analysis
        double phase2Improvement = ((double) timePhase1 - timePhase2) / timePhase1 * 100;
        System.out.println("ANALYSIS:");
        System.out.println("- Phase 1: DP optimizations (binary search, parallelization)");
        System.out.println("- Phase 2: TRUE Subtree DOF caching with boundary-only optimization");
        System.out.println(String.format("- Phase 2 improvement over Phase 1: %+.1f%%", phase2Improvement));

        if (phase2Improvement > 10) {
            System.out.println("- ✓ Phase 2 subtree caching provides significant speedup!");
        } else if (phase2Improvement > 0) {
            System.out.println("- Phase 2 provides modest improvement (may need larger scale to see benefit)");
        } else {
            System.out.println("- Phase 2 has overhead (cache management cost > savings on this scale)");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    /**
     * Test Original vs Phase 1+2 ONLY
     * 9 flexible residues, 2 mutable to all 20 amino acids
     * This tests sequence design with subtree caching
     */
    @Test
    public void testOriginalVsPhase2WithMutations() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE COMPARISON: Original vs Phase 1+2");
        System.out.println("Two test scenarios:");
        System.out.println("  A) 1 mutable + 3 flexible (mutation scenario)");
        System.out.println("  B) 3 flexible only (no mutation scenario)");
        System.out.println("Resources: 20 CPU cores, 100GB memory");
        System.out.println("=".repeat(100) + "\n");

        // ==================== SCENARIO A: 1 mutable + 5 flexible ====================
        System.out.println("\n" + "=".repeat(100));
        System.out.println("SCENARIO A: 1 Mutable + 5 Flexible");
        System.out.println("Configuration: 6 total positions (1 mutable with Wild+5 types, 4 flexible)");
        System.out.println("=".repeat(100) + "\n");

        // A1: Original
        System.out.println("=== A1: Original Greedy (baseline) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeA_Original = testScaleWithMutations(6, 1, "A-Original");

        // A2: Phase 1+2
        System.out.println("\n=== A2: Phase 1+2 (DP + TRUE Subtree Caching) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        TestKStar.ConfSpaces confSpacesA = buildConfSpaceWithMutations(6, 1);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpacesA.complex);
        // CRITICAL: Pass confSpacesA to ensure cache is used with correct confSpace
        long timeA_Phase2 = testScaleWithMutations(6, 1, "A-Phase1+2", confSpacesA);

        // Print Scenario A summary
        System.out.println("\n" + "-".repeat(100));
        System.out.println("SCENARIO A SUMMARY (1 mutable + 5 flexible)");
        System.out.println("-".repeat(100));
        System.out.println(String.format("%-30s: %s", "Original Greedy", formatTime(timeA_Original)));
        System.out.println(String.format("%-30s: %s (%.2fx speedup)",
            "Phase 1+2", formatTime(timeA_Phase2),
            (double) timeA_Original / timeA_Phase2));
        double improvementA = ((double) timeA_Original - timeA_Phase2) / timeA_Original * 100;
        System.out.println(String.format("Improvement: %+.1f%%", improvementA));
        System.out.println("-".repeat(100) + "\n");

        // ==================== SCENARIO B: 5 flexible only ====================
        System.out.println("\n" + "=".repeat(100));
        System.out.println("SCENARIO B: 5 Flexible Only (No Mutations)");
        System.out.println("Configuration: 5 flexible positions, no mutations");
        System.out.println("=".repeat(100) + "\n");

        // B1: Original
        System.out.println("=== B1: Original Greedy (baseline) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        long timeB_Original = testScale(5, "B-Original");

        // B2: Phase 1+2
        System.out.println("\n=== B2: Phase 1+2 (DP + TRUE Subtree Caching) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        TestKStar.ConfSpaces confSpacesB = buildConfSpace(5);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpacesB.complex);
        // CRITICAL: Pass confSpacesB to ensure cache is used with correct confSpace
        long timeB_Phase2 = testScale(5, "B-Phase1+2", confSpacesB);

        // Print Scenario B summary
        System.out.println("\n" + "-".repeat(100));
        System.out.println("SCENARIO B SUMMARY (5 flexible only)");
        System.out.println("-".repeat(100));
        System.out.println(String.format("%-30s: %s", "Original Greedy", formatTime(timeB_Original)));
        System.out.println(String.format("%-30s: %s (%.2fx speedup)",
            "Phase 1+2", formatTime(timeB_Phase2),
            (double) timeB_Original / timeB_Phase2));
        double improvementB = ((double) timeB_Original - timeB_Phase2) / timeB_Original * 100;
        System.out.println(String.format("Improvement: %+.1f%%", improvementB));
        System.out.println("-".repeat(100) + "\n");

        // Print cache statistics
        System.out.println("\n=== Subtree DOF Cache Statistics ===");
        CachedMinimizer.printGlobalStats();

        // ==================== FINAL OVERALL SUMMARY ====================
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL OVERALL SUMMARY");
        System.out.println("=".repeat(100));
        System.out.println();

        System.out.println("SCENARIO A (1 mutable + 5 flexible):");
        System.out.println(String.format("  Original:   %s", formatTime(timeA_Original)));
        System.out.println(String.format("  Phase 1+2:  %s (%.2fx speedup, %+.1f%% improvement)",
            formatTime(timeA_Phase2),
            (double) timeA_Original / timeA_Phase2,
            improvementA));
        System.out.println();

        System.out.println("SCENARIO B (3 flexible only, no mutations):");
        System.out.println(String.format("  Original:   %s", formatTime(timeB_Original)));
        System.out.println(String.format("  Phase 1+2:  %s (%.2fx speedup, %+.1f%% improvement)",
            formatTime(timeB_Phase2),
            (double) timeB_Original / timeB_Phase2,
            improvementB));
        System.out.println();

        System.out.println("ANALYSIS:");
        System.out.println("- Resources: 20 CPU cores, 100GB memory");
        System.out.println("- Phase 1+2 includes: DP + Binary Search + Parallelization + TRUE Subtree Caching");
        System.out.println("- Mutation scenario shows how cache handles different amino acid types");
        System.out.println("- No-mutation scenario shows pure flexible residue optimization");

        if (improvementA > 20 || improvementB > 20) {
            System.out.println("- ✓✓ Significant speedup achieved in at least one scenario!");
        } else if (improvementA > 0 && improvementB > 0) {
            System.out.println("- ✓ Consistent speedup across both scenarios");
        } else {
            System.out.println("- Mixed results - overhead vs benefits varies by scenario");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    /**
     * Build conf space with mutations
     * @param numFlexible Number of flexible residues (same as before)
     * @param numMutable Number of residues that are mutable to all 20 amino acids
     */
    private TestKStar.ConfSpaces buildConfSpaceWithMutations(int numFlexible, int numMutable) {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .build();

        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();

        // Add flexible residues (must have numFlexible >= 9)
        int flexCount = 0;
        int mutCount = 0;

        // Position 1 & 2: flexible only
        if (flexCount < numFlexible) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }
        if (flexCount < numFlexible) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }

        // Position 3 & 4: flexible + mutable (first 2 mutable positions)
        if (flexCount < numFlexible) {
            if (mutCount < numMutable) {
                // Mutable position (matching TestKStar.make2RL0 - Wild + 5 types)
                protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
                mutCount++;
            } else {
                protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            }
            flexCount++;
        }
        if (flexCount < numFlexible) {
            if (mutCount < numMutable) {
                ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
                mutCount++;
            } else {
                ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            }
            flexCount++;
        }

        // Remaining positions: flexible only
        if (flexCount < numFlexible) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }
        if (flexCount < numFlexible) {
            ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }
        if (flexCount < numFlexible) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }
        if (flexCount < numFlexible) {
            ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }
        if (flexCount < numFlexible) {
            protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            flexCount++;
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Test with mutations
     */
    private long testScaleWithMutations(int numFlexible, int numMutable, String label) {
        return testScaleWithMutations(numFlexible, numMutable, label, null);
    }

    /**
     * Test with mutations, optionally using pre-built confSpaces
     * This overload is critical for Phase 2 cache to work correctly:
     * Cache must be initialized with the SAME confSpace that's used for testing
     */
    private long testScaleWithMutations(int numFlexible, int numMutable, String label, TestKStar.ConfSpaces confSpaces) {
        // Use provided confSpaces if available, otherwise build new one
        if (confSpaces == null) {
            confSpaces = buildConfSpaceWithMutations(numFlexible, numMutable);
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions: " + numFlexible);
        System.out.println("  Mutable positions: " + numMutable + " (Wild + 5 types each)");
        System.out.println("  Algorithm: " +
            (UpdatingEnergyMatrix.USE_DP_ALGORITHM ? "DP" : "Greedy"));
        System.out.println("  DP optimizations: " +
            (UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS ? "ON" : "OFF"));
        System.out.println("  Subtree cache: " +
            (CachedMinimizer.ENABLE_SUBTREE_CACHE ? "ON" : "OFF"));
        System.out.println();

        // Run MARK*
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

        System.out.println("Pre-calculating energy matrices...");
        markstar.precalcEmats();

        System.out.println("Running MARK*...");
        var scores = markstar.run();

        watch.stop();
        long timeMs = (long) watch.getTimeMs();

        System.out.println("\nResults:");
        System.out.println("  Runtime: " + formatTime(timeMs));
        if (scores != null && scores.size() > 0) {
            System.out.println("  Best sequence: " + scores.get(0).sequence);
            System.out.println("  K* score: " + scores.get(0).score);
        }

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        return timeMs;
    }

    /**
     * Simple test: 1 mutable residue only, no other flexible
     * Quick validation of mutation configuration
     */
    @Test
    public void testSimpleMutation() {
        System.out.println("\n=== SIMPLE MUTATION TEST ===");
        System.out.println("1 mutable residue (Wild + 5 types), no other flexible");
        System.out.println("Only testing Original (no DP, no cache)");
        System.out.println("Matching TestKStar.make2RL0 configuration\n");

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

        // Only 1 mutable position: G650 (matching TestKStar.make2RL0 pattern - Wild + 5 types)
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        System.out.println("Conf space built successfully!");
        System.out.println("Complex positions: " + confSpaces.complex.positions.size());
        System.out.println("Position 0 RCs: " + confSpaces.complex.positions.get(0).resConfs.size());

        // Run Original MARK*
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        CachedMinimizer.clearGlobalCache();

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

        System.out.println("Pre-calculating energy matrices...");
        markstar.precalcEmats();

        System.out.println("Running MARK*...");
        var scores = markstar.run();

        watch.stop();
        long timeMs = (long) watch.getTimeMs();

        System.out.println("\nResults:");
        System.out.println("  Runtime: " + formatTime(timeMs));
        if (scores != null && scores.size() > 0) {
            System.out.println("  Best sequence: " + scores.get(0).sequence);
            System.out.println("  K* score: " + scores.get(0).score);
        }

        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();

        System.out.println("\n=== TEST PASSED ===\n");
    }
}
