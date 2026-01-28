package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.CachedMinimizer;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SubtreeDOFCache;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Test to compare CCD Minimizer energy values across three configurations:
 * 1. Original (no cache) - Pure CCD minimization
 * 2. Subtree Cache Only - CCD with subtree DOF caching
 * 3. Subtree + Triple Cache - CCD with both subtree and triple DOF caching
 *
 * The goal is to verify that all three methods produce identical (or nearly identical)
 * minimized energy values for the same conformations.
 */
public class TestCCDMinimizerEnergyComparison {

    private static final int NUM_CPUs = 4;  // Use fewer CPUs for deterministic results
    private static final double ENERGY_TOLERANCE = 0.001;  // Tolerance for energy comparison

    /**
     * Main test: Compare CCD minimizer energy across three configurations
     * for a set of test conformations.
     */
    @Test
    public void testCCDMinimizerEnergyConsistency() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("CCD MINIMIZER ENERGY COMPARISON TEST");
        System.out.println("Comparing minimized energies across:");
        System.out.println("  1. Original (No Cache) - Pure CCD minimization");
        System.out.println("  2. Subtree Cache Only - CCD with subtree DOF caching");
        System.out.println("  3. Subtree + Triple Cache - CCD with both caching methods");
        System.out.println("=".repeat(100) + "\n");

        // Build conf space
        TestKStar.ConfSpaces confSpaces = buildTestConfSpace();

        // Generate test conformations
        List<int[]> testConformations = generateTestConformations(confSpaces.complex);
        System.out.println("Generated " + testConformations.size() + " test conformations\n");

        // Store energy results for comparison
        List<EnergyResult> results = new ArrayList<>();

        // Test 1: Original (No Cache)
        System.out.println("=== TEST 1: Original (No Cache) ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        CachedMinimizer.clearGlobalCache();

        List<Double> energiesOriginal = minimizeConformations(confSpaces, testConformations, "Original");

        // Test 2: Subtree Cache Only
        System.out.println("\n=== TEST 2: Subtree Cache Only ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = false;
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        List<Double> energiesSubtree = minimizeConformations(confSpaces, testConformations, "Subtree");

        // Test 3: Subtree + Triple Cache
        System.out.println("\n=== TEST 3: Subtree + Triple Cache ===\n");
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
        SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE = true;
        CachedMinimizer.clearGlobalCache();
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);

        List<Double> energiesTriple = minimizeConformations(confSpaces, testConformations, "Triple");

        // Compare results
        compareAndReportResults(testConformations, energiesOriginal, energiesSubtree, energiesTriple);
    }

    /**
     * Minimize conformations and return energies
     */
    private List<Double> minimizeConformations(TestKStar.ConfSpaces confSpaces,
                                                List<int[]> conformations,
                                                String label) {
        List<Double> energies = new ArrayList<>();

        Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {

            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpaces.complex, ecalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, ecalc)
                            .build()
                            .calcReferenceEnergies())
                    .build();

            for (int i = 0; i < conformations.size(); i++) {
                int[] conf = conformations.get(i);
                RCTuple tuple = new RCTuple(conf);

                // Minimize the conformation
                double energy = confEcalc.calcEnergy(tuple).energy;
                energies.add(energy);

                if (i < 10 || i == conformations.size() - 1) {
                    System.out.println(String.format("[%s] Conf %d: %s -> E=%.6f",
                        label, i, formatConf(conf), energy));
                } else if (i == 10) {
                    System.out.println("  ... (showing first 10 and last conformation)");
                }
            }

            ecalc.tasks.waitForFinish();
        }

        return energies;
    }

    /**
     * Compare energy results and report findings
     */
    private void compareAndReportResults(List<int[]> conformations,
                                          List<Double> energiesOriginal,
                                          List<Double> energiesSubtree,
                                          List<Double> energiesTriple) {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("ENERGY COMPARISON RESULTS");
        System.out.println("=".repeat(100) + "\n");

        int exactMatch = 0;
        int closeMatch = 0;  // < 0.01 difference
        int moderateDiff = 0;  // < 0.1 difference
        int largeDiff = 0;  // >= 0.1 difference

        double maxDiff = 0;
        int maxDiffIndex = -1;
        double sumDiff = 0;
        double sumRelDiff = 0;

        System.out.println("Detailed comparison for each conformation:");
        System.out.println("-".repeat(120));
        System.out.println(String.format("%-30s | %-15s | %-15s | %-15s | %-10s | %-10s",
            "Conformation", "Original", "Subtree", "Triple", "Max Diff", "Status"));
        System.out.println("-".repeat(120));

        for (int i = 0; i < conformations.size(); i++) {
            double e1 = energiesOriginal.get(i);
            double e2 = energiesSubtree.get(i);
            double e3 = energiesTriple.get(i);

            double diff12 = Math.abs(e1 - e2);
            double diff13 = Math.abs(e1 - e3);
            double diff23 = Math.abs(e2 - e3);
            double maxDiffForConf = Math.max(Math.max(diff12, diff13), diff23);

            sumDiff += maxDiffForConf;
            if (Math.abs(e1) > 0.001) {
                sumRelDiff += maxDiffForConf / Math.abs(e1);
            }

            if (maxDiffForConf > maxDiff) {
                maxDiff = maxDiffForConf;
                maxDiffIndex = i;
            }

            String status;
            if (maxDiffForConf < 0.001) {
                exactMatch++;
                status = "EXACT";
            } else if (maxDiffForConf < 0.01) {
                closeMatch++;
                status = "CLOSE";
            } else if (maxDiffForConf < 0.1) {
                moderateDiff++;
                status = "MODERATE";
            } else {
                largeDiff++;
                status = "LARGE DIFF!";
            }

            System.out.println(String.format("%-30s | %15.6f | %15.6f | %15.6f | %10.6f | %s",
                formatConf(conformations.get(i)), e1, e2, e3, maxDiffForConf, status));
        }

        System.out.println("-".repeat(120));

        // Summary statistics
        int total = conformations.size();
        double avgDiff = sumDiff / total;
        double avgRelDiff = sumRelDiff / total * 100;

        System.out.println("\n" + "=".repeat(100));
        System.out.println("SUMMARY STATISTICS");
        System.out.println("=".repeat(100));
        System.out.println(String.format("Total conformations tested: %d", total));
        System.out.println();
        System.out.println("Energy difference categories:");
        System.out.println(String.format("  - Exact match (diff < 0.001):     %d (%.1f%%)", exactMatch, 100.0 * exactMatch / total));
        System.out.println(String.format("  - Close match (diff < 0.01):      %d (%.1f%%)", closeMatch, 100.0 * closeMatch / total));
        System.out.println(String.format("  - Moderate diff (diff < 0.1):     %d (%.1f%%)", moderateDiff, 100.0 * moderateDiff / total));
        System.out.println(String.format("  - Large diff (diff >= 0.1):       %d (%.1f%%)", largeDiff, 100.0 * largeDiff / total));
        System.out.println();
        System.out.println(String.format("Average max difference: %.6f", avgDiff));
        System.out.println(String.format("Average relative difference: %.3f%%", avgRelDiff));
        System.out.println(String.format("Maximum difference: %.6f (conf %d: %s)",
            maxDiff, maxDiffIndex, maxDiffIndex >= 0 ? formatConf(conformations.get(maxDiffIndex)) : "N/A"));

        // Final verdict
        System.out.println("\n" + "=".repeat(100));
        if (largeDiff == 0 && moderateDiff == 0) {
            System.out.println("VERDICT: ALL THREE METHODS PRODUCE CONSISTENT RESULTS");
            System.out.println("All energy differences are within acceptable tolerance (< 0.01)");
        } else if (largeDiff == 0) {
            System.out.println("VERDICT: METHODS ARE MOSTLY CONSISTENT");
            System.out.println("Some moderate differences exist but no large discrepancies");
        } else {
            System.out.println("VERDICT: SIGNIFICANT DIFFERENCES DETECTED");
            System.out.println("Review the conformations with large differences above");
        }
        System.out.println("=".repeat(100) + "\n");
    }

    /**
     * Build a test conf space with 5 flexible residues
     */
    private TestKStar.ConfSpaces buildTestConfSpace() {
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

        // 5 flexible residues for testing
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }

    /**
     * Generate a set of test conformations
     */
    private List<int[]> generateTestConformations(SimpleConfSpace confSpace) {
        List<int[]> conformations = new ArrayList<>();
        int numPositions = confSpace.positions.size();

        // Get the number of rotamers at each position
        int[] numRotamers = new int[numPositions];
        for (int i = 0; i < numPositions; i++) {
            numRotamers[i] = confSpace.positions.get(i).resConfs.size();
        }

        // Generate some test conformations
        // Start with the all-0 conformation
        int[] conf0 = new int[numPositions];
        conformations.add(conf0);

        // Generate conformations varying one position at a time
        for (int pos = 0; pos < numPositions; pos++) {
            for (int rc = 0; rc < Math.min(numRotamers[pos], 3); rc++) {
                int[] conf = new int[numPositions];
                conf[pos] = rc;
                if (!containsConf(conformations, conf)) {
                    conformations.add(conf);
                }
            }
        }

        // Add some random-ish conformations
        for (int i = 0; i < 10; i++) {
            int[] conf = new int[numPositions];
            for (int pos = 0; pos < numPositions; pos++) {
                conf[pos] = (i + pos) % numRotamers[pos];
            }
            if (!containsConf(conformations, conf)) {
                conformations.add(conf);
            }
        }

        return conformations;
    }

    private boolean containsConf(List<int[]> list, int[] conf) {
        for (int[] existing : list) {
            boolean same = true;
            for (int i = 0; i < conf.length; i++) {
                if (conf[i] != existing[i]) {
                    same = false;
                    break;
                }
            }
            if (same) return true;
        }
        return false;
    }

    private String formatConf(int[] conf) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < conf.length; i++) {
            if (i > 0) sb.append(" ");
            sb.append(conf[i]);
        }
        return sb.toString();
    }

    /**
     * Simple result container
     */
    private static class EnergyResult {
        int[] conf;
        double energyOriginal;
        double energySubtree;
        double energyTriple;
    }
}
