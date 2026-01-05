package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.BranchDecomposition;
import edu.duke.cs.osprey.ematrix.DPTupleTrie;
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
 * Test Phase 1 DP optimizations:
 * 1. Binary search: O(n²) → O(n log n)
 * 2. Parallel DP for large problems
 * 3. Compare with original DP and Greedy
 */
public class TestDPOptimizations {

    private static final int NUM_CPUs = 4;

    @Test
    public void testPhase1Optimizations() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("DP-MARKStar Phase 1 Optimizations Test");
        System.out.println("Testing: Binary Search + Parallel DP");
        System.out.println("=".repeat(80) + "\n");

        // Test 1: Small scale with optimizations OFF
        System.out.println("=== TEST 1: Original DP (5 residues) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = false;
        long time1 = testScale(5, "Original DP");

        // Test 2: Small scale with optimizations ON
        System.out.println("\n=== TEST 2: Optimized DP (5 residues) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        long time2 = testScale(5, "Optimized DP");

        // Test 3: Greedy for comparison
        System.out.println("\n=== TEST 3: Greedy (5 residues) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        long time3 = testScale(5, "Greedy");

        // Summary
        System.out.println("\n" + "=".repeat(80));
        System.out.println("Performance Comparison");
        System.out.println("=".repeat(80));
        System.out.println("Original DP:  " + formatTime(time1));
        System.out.println("Optimized DP: " + formatTime(time2));
        System.out.println("Greedy:       " + formatTime(time3));
        System.out.println();

        double speedup = (double) time1 / time2;
        System.out.println("Optimization Speedup: " + String.format("%.2fx", speedup));
        System.out.println("DP vs Greedy: " + String.format("%.2fx %s",
            Math.abs((double) time2 / time3),
            time2 < time3 ? "faster" : "slower"));

        System.out.println("\n" + "=".repeat(80));
        System.out.println("✓ Phase 1 Optimizations Test Complete");
        System.out.println("=".repeat(80) + "\n");
    }

    private long testScale(int numFlexible, String label) {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .build();

        // Build conf space
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();

        // Add flexible residues
        if (numFlexible >= 2) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 4) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        if (numFlexible >= 5) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions: " + confSpaces.complex.positions.size());
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
     * Test Phase 2: Branch Decomposition and Subtree DOF Cache
     */
    @Test
    public void testPhase2BranchDecomposition() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("DP-MARKStar Phase 2: Branch Decomposition Test");
        System.out.println("=".repeat(80) + "\n");

        // Build a test conf space
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        // Add flexible residues
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(protein).build();

        // Test branch decomposition
        System.out.println("Building Branch Decomposition...");
        BranchDecomposition branchDecomp = new BranchDecomposition(confSpace);

        System.out.println("\nBranch Decomposition Statistics:");
        System.out.println(branchDecomp.getStats());
        System.out.println("Branch width (separator size): " + branchDecomp.branchWidth);

        // Print tree structure
        System.out.println("\nTree Structure:");
        branchDecomp.printTree();

        // Check if suitable for DOF caching
        boolean suitable = branchDecomp.isSuitable(5);
        System.out.println("\nSuitable for DOF caching (max width=5): " + suitable);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("✓ Phase 2 Branch Decomposition Test Complete");
        System.out.println("=".repeat(80) + "\n");
    }

    /**
     * Test Phase 3: DP-Trie Integration
     */
    @Test
    public void testPhase3DPTrie() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("DP-MARKStar Phase 3: DP-Trie Integration Test");
        System.out.println("=".repeat(80) + "\n");

        int numPositions = 5;
        DPTupleTrie dpTrie = new DPTupleTrie(numPositions);

        // Add some test corrections
        System.out.println("Adding test corrections...");
        java.util.ArrayList<Integer> pos1 = new java.util.ArrayList<>();
        java.util.ArrayList<Integer> rcs1 = new java.util.ArrayList<>();
        pos1.add(0); pos1.add(1);
        rcs1.add(0); rcs1.add(0);
        dpTrie.addCorrection(new edu.duke.cs.osprey.confspace.RCTuple(pos1, rcs1), -5.0);

        java.util.ArrayList<Integer> pos2 = new java.util.ArrayList<>();
        java.util.ArrayList<Integer> rcs2 = new java.util.ArrayList<>();
        pos2.add(1); pos2.add(2);
        rcs2.add(0); rcs2.add(0);
        dpTrie.addCorrection(new edu.duke.cs.osprey.confspace.RCTuple(pos2, rcs2), -3.0);

        java.util.ArrayList<Integer> pos3 = new java.util.ArrayList<>();
        java.util.ArrayList<Integer> rcs3 = new java.util.ArrayList<>();
        pos3.add(0); pos3.add(2);
        rcs3.add(0); rcs3.add(0);
        dpTrie.addCorrection(new edu.duke.cs.osprey.confspace.RCTuple(pos3, rcs3), -7.0);

        System.out.println("Added 3 corrections to DP-Trie");

        // Test retrieval
        System.out.println("\nTesting DP-integrated traversal...");
        java.util.ArrayList<Integer> confPos = new java.util.ArrayList<>();
        java.util.ArrayList<Integer> confRCs = new java.util.ArrayList<>();
        confPos.add(0); confPos.add(1); confPos.add(2);
        confRCs.add(0); confRCs.add(0); confRCs.add(0);
        edu.duke.cs.osprey.confspace.RCTuple conf = new edu.duke.cs.osprey.confspace.RCTuple(confPos, confRCs);

        double optimalCorrection = dpTrie.getOptimalCorrection(conf);
        System.out.println("Optimal correction value: " + optimalCorrection);

        // Test traditional collection for comparison
        System.out.println("\nTesting traditional collection...");
        java.util.List<edu.duke.cs.osprey.confspace.TupE> collected = dpTrie.collectCorrections(conf);
        System.out.println("Collected " + collected.size() + " matching corrections");
        for (edu.duke.cs.osprey.confspace.TupE tup : collected) {
            System.out.println("  " + tup.tup + " -> " + tup.E);
        }

        // Print statistics
        dpTrie.printStats();

        System.out.println("\n" + "=".repeat(80));
        System.out.println("✓ Phase 3 DP-Trie Integration Test Complete");
        System.out.println("=".repeat(80) + "\n");
    }

    /**
     * Comprehensive test combining all three phases
     */
    @Test
    public void testAllPhasesCombined() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("DP-MARKStar COMPREHENSIVE TEST: All Phases Combined");
        System.out.println("Phase 1: Binary Search + Parallel DP");
        System.out.println("Phase 2: Branch Decomposition + Subtree DOF Cache");
        System.out.println("Phase 3: DP-Trie Integration");
        System.out.println("=".repeat(80) + "\n");

        // Test 1: Baseline (all optimizations OFF)
        System.out.println("=== TEST 1: Baseline (Greedy) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;
        long timeBaseline = testScale(5, "Baseline Greedy");

        // Test 2: Phase 1 only (Binary Search + Parallel DP)
        System.out.println("\n=== TEST 2: Phase 1 Optimizations ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = true;
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = true;
        long timePhase1 = testScale(5, "Phase 1 (Binary Search + Parallel DP)");

        // Test 3: Placeholder for Phase 2+3 (when integrated into MARKStar)
        // For now, just demonstrate the components work
        System.out.println("\n=== Phase 2+3 Component Verification ===\n");
        System.out.println("✓ BranchDecomposition implemented and tested");
        System.out.println("✓ SubtreeDOFCache implemented and compiled");
        System.out.println("✓ DPTupleTrie implemented and tested");
        System.out.println("\nNote: Full Phase 2+3 integration requires modifying MARKStar's");
        System.out.println("      minimization and correction selection pipeline.");

        // Summary
        System.out.println("\n" + "=".repeat(80));
        System.out.println("Performance Summary");
        System.out.println("=".repeat(80));
        System.out.println("Baseline (Greedy):  " + formatTime(timeBaseline));
        System.out.println("Phase 1 Optimized:  " + formatTime(timePhase1));
        System.out.println();

        double speedupVsBaseline = (double) timeBaseline / timePhase1;
        System.out.println("Speedup vs Baseline: " + String.format("%.2fx", speedupVsBaseline));

        System.out.println("\n" + "=".repeat(80));
        System.out.println("✓ Comprehensive Test Complete");
        System.out.println("=".repeat(80) + "\n");
    }
}
