package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.CachedMinimizer;
import edu.duke.cs.osprey.ematrix.PartialFixIntegration;
import edu.duke.cs.osprey.ematrix.PartialStartIntegration;
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
 * Phase 5: Partial Start Cache Performance Test
 *
 * Tests three configurations:
 * 1. Original (no cache) - baseline
 * 2. Phase 5: Partial Start Cache only (warm start from 5/6-subset matches)
 * 3. All caches enabled (Subtree + Triple + PartialFix + PartialStart)
 *
 * Partial Start Cache provides warm start DOF values by:
 * - Storing 7-tuple minimization results as 5-subset and 6-subset entries
 * - When querying, find matching subsets and use cached DOF values
 * - For overlapping positions: use cached DOF values (warm start)
 * - For non-overlapping positions: use default DOF values
 */
public class TestPartialStartCache {

    private static final int NUM_CPUs = 20;

    /**
     * Main test: Compare Original vs Phase 5 (Partial Start Cache)
     */
    @Test
    public void testPartialStartCachePerformance() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PHASE 5: PARTIAL START CACHE PERFORMANCE TEST");
        System.out.println("Comparing: Original vs Partial Start Cache");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        disableAllCaches();
        long timeOriginal = runTest(scale, "Original-NoCache");

        // Test 2: Partial Start Cache only
        System.out.println("\n=== TEST 2: Partial Start Cache Only ===\n");
        disableAllCaches();
        PartialStartIntegration.ENABLE_PARTIALSTART_CACHE = true;
        long timePartialStart = runTest(scale, "PartialStartOnly");

        // Print comparison
        printComparison(timeOriginal, timePartialStart);
    }

    /**
     * Extended test: Compare all cache configurations
     */
    @Test
    public void testAllCacheConfigurations() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FULL CACHE COMPARISON TEST");
        System.out.println("Comparing: Original vs Subtree vs Subtree+Triple vs Subtree+Triple+PartialFix vs All");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        disableAllCaches();
        long timeOriginal = runTest(scale, "Original-NoCache");

        // Test 2: Subtree Cache only
        System.out.println("\n=== TEST 2: Subtree Cache Only ===\n");
        disableAllCaches();
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        long timeSubtree = runTest(scale, "SubtreeOnly", confSpaces);

        // Test 3: Subtree + Triple Cache
        System.out.println("\n=== TEST 3: Subtree + Triple Cache ===\n");
        disableAllCaches();
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;
        confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        long timeSubtreeTriple = runTest(scale, "Subtree+Triple", confSpaces);

        // Test 4: Subtree + Triple + PartialFix Cache
        System.out.println("\n=== TEST 4: Subtree + Triple + PartialFix Cache ===\n");
        disableAllCaches();
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;
        confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        long timeSubtreeTriplePartialFix = runTest(scale, "Subtree+Triple+PartialFix", confSpaces);

        // Test 5: All Caches (Subtree + Triple + PartialFix + PartialStart)
        System.out.println("\n=== TEST 5: All Caches (+ Partial Start) ===\n");
        disableAllCaches();
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = true;
        PartialStartIntegration.ENABLE_PARTIALSTART_CACHE = true;
        confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        long timeAll = runTest(scale, "AllCaches", confSpaces);

        // Print full comparison
        printFullComparison(timeOriginal, timeSubtree, timeSubtreeTriple,
                           timeSubtreeTriplePartialFix, timeAll);
    }

    /**
     * Test Partial Start Cache only (isolated from other caches)
     */
    @Test
    public void testPartialStartCacheIsolated() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIAL START CACHE ISOLATED TEST");
        System.out.println("Comparing: Original vs Partial Start Cache (no other caches)");
        System.out.println("Configuration: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (no cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        disableAllCaches();
        long timeOriginal = runTest(scale, "Original-NoCache");

        // Test 2: Partial Start Cache only (isolated)
        System.out.println("\n=== TEST 2: Partial Start Cache Only (Isolated) ===\n");
        disableAllCaches();
        PartialStartIntegration.ENABLE_PARTIALSTART_CACHE = true;
        long timePartialStart = runTest(scale, "PartialStartOnly");

        // Print comparison
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PARTIAL START CACHE ISOLATED COMPARISON");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "Partial Start Cache Only", formatTime(timePartialStart),
            (double) timeOriginal / timePartialStart));

        double improvement = ((double) timeOriginal - timePartialStart) / timeOriginal * 100;
        System.out.println();
        System.out.println(String.format("Improvement: %+.1f%%", improvement));

        if (improvement > 20) {
            System.out.println("✓✓ Partial Start Cache provides SIGNIFICANT speedup!");
        } else if (improvement > 10) {
            System.out.println("✓ Partial Start Cache provides good speedup");
        } else if (improvement > 0) {
            System.out.println("Partial Start Cache provides modest improvement");
        } else {
            System.out.println("Partial Start Cache has overhead (may benefit larger problems)");
        }

        // Print cache statistics
        PartialStartIntegration.printGlobalStats();

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    /**
     * Disable all caches for clean baseline
     */
    private void disableAllCaches() {
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        PartialFixIntegration.ENABLE_PARTIALFIX_CACHE = false;
        PartialStartIntegration.ENABLE_PARTIALSTART_CACHE = false;
        CachedMinimizer.clearGlobalCache();
    }

    /**
     * Run a single test configuration
     */
    private long runTest(int numFlexible, String label) {
        return runTest(numFlexible, label, null);
    }

    private long runTest(int numFlexible, String label, TestKStar.ConfSpaces confSpaces) {
        if (confSpaces == null) {
            confSpaces = buildConfSpace(numFlexible);
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions:    " + confSpaces.complex.positions.size());
        System.out.println("  Subtree cache:         " + (CachedMinimizer.ENABLE_SUBTREE_CACHE ? "ON" : "OFF"));
        System.out.println("  Triple cache:          " + (SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE ? "ON" : "OFF"));
        System.out.println("  PartialFix cache:      " + (PartialFixIntegration.ENABLE_PARTIALFIX_CACHE ? "ON" : "OFF"));
        System.out.println("  PartialStart cache:    " + (PartialStartIntegration.ENABLE_PARTIALSTART_CACHE ? "ON" : "OFF"));
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
                .setShowPfuncProgress(false)
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

        // Print cache statistics if enabled
        if (PartialStartIntegration.ENABLE_PARTIALSTART_CACHE) {
            PartialStartIntegration.printGlobalStats();
        }

        return timeMs;
    }

    /**
     * Build conf space with specified number of flexible residues
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

        // Add mutable residues with expanded amino acid options for more diversity
        if (numFlexible >= 1) {
            protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "SER", "THR").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 2) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "ASN", "GLN").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 3) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "GLU", "ASP", "LYS").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 4) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "ASN", "GLN", "ARG").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 5) {
            protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "SER", "THR", "HIS").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 6) {
            protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "GLY", "PRO", "MET").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 7) {
            protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "LYS", "ARG", "GLU").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 8) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE").addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 9) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE").addWildTypeRotamers().setContinuous();
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Print comparison between Original and Partial Start Cache
     */
    private void printComparison(long timeOriginal, long timePartialStart) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("PERFORMANCE COMPARISON: Original vs Partial Start Cache");
        System.out.println("=".repeat(100));
        System.out.println(String.format("%-35s: %s (baseline)", "Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-35s: %s (%.2fx vs Original)",
            "Partial Start Cache", formatTime(timePartialStart),
            (double) timeOriginal / timePartialStart));
        System.out.println();

        double improvement = ((double) timeOriginal - timePartialStart) / timeOriginal * 100;
        System.out.println(String.format("Improvement: %+.1f%%", improvement));
        System.out.println();

        System.out.println("PHASE 5 DETAILS:");
        System.out.println("- Stores 7-tuple results as 5-subset (21 indices) and 6-subset (7 indices)");
        System.out.println("- Query: O(28) hash lookups instead of O(N) linear scan");
        System.out.println("- 6-subset match = 6 positions overlap (better warm start)");
        System.out.println("- 5-subset match = 5 positions overlap (good warm start)");
        System.out.println("- Non-overlapping positions use default DOF values");

        if (improvement > 30) {
            System.out.println("\n✓✓✓ Partial Start Cache provides EXCELLENT speedup!");
        } else if (improvement > 20) {
            System.out.println("\n✓✓ Partial Start Cache provides SIGNIFICANT speedup!");
        } else if (improvement > 10) {
            System.out.println("\n✓ Partial Start Cache provides good speedup");
        } else if (improvement > 0) {
            System.out.println("\nPartial Start Cache provides modest improvement");
        } else {
            System.out.println("\nPartial Start Cache has overhead (may benefit larger problems)");
        }

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    /**
     * Print full comparison of all cache configurations
     */
    private void printFullComparison(long timeOriginal, long timeSubtree, long timeSubtreeTriple,
                                     long timeSubtreeTriplePartialFix, long timeAll) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FULL CACHE PERFORMANCE COMPARISON");
        System.out.println("=".repeat(100));

        System.out.println(String.format("%-45s: %s (baseline)",
            "1. Original (No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-45s: %s (%.2fx)",
            "2. Subtree Cache", formatTime(timeSubtree),
            (double) timeOriginal / timeSubtree));
        System.out.println(String.format("%-45s: %s (%.2fx)",
            "3. Subtree + Triple Cache", formatTime(timeSubtreeTriple),
            (double) timeOriginal / timeSubtreeTriple));
        System.out.println(String.format("%-45s: %s (%.2fx)",
            "4. Subtree + Triple + PartialFix", formatTime(timeSubtreeTriplePartialFix),
            (double) timeOriginal / timeSubtreeTriplePartialFix));
        System.out.println(String.format("%-45s: %s (%.2fx)",
            "5. All Caches (+ PartialStart)", formatTime(timeAll),
            (double) timeOriginal / timeAll));

        System.out.println();
        System.out.println("INCREMENTAL IMPROVEMENTS:");
        System.out.println(String.format("  Subtree vs Original:              %+.1f%%",
            improvement(timeOriginal, timeSubtree)));
        System.out.println(String.format("  Triple vs Subtree:                %+.1f%%",
            improvement(timeSubtree, timeSubtreeTriple)));
        System.out.println(String.format("  PartialFix vs Triple:             %+.1f%%",
            improvement(timeSubtreeTriple, timeSubtreeTriplePartialFix)));
        System.out.println(String.format("  PartialStart vs PartialFix:       %+.1f%%",
            improvement(timeSubtreeTriplePartialFix, timeAll)));
        System.out.println();
        System.out.println(String.format("  TOTAL vs Original:                %+.1f%%",
            improvement(timeOriginal, timeAll)));

        System.out.println("\n" + "=".repeat(100) + "\n");
    }

    private double improvement(long baseline, long improved) {
        return ((double) baseline - improved) / baseline * 100;
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
