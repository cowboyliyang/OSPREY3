#!/bin/bash
#SBATCH --job-name=osprey_cache_only
#SBATCH --output=test_cache_only_%j.out
#SBATCH --error=test_cache_only_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --partition=compsci

# Load Java module
module load Java/17.0.4

# Set working directory
cd /home/users/lz280/IdeaProjects/OSPREY3

echo "========================================"
echo "OSPREY3 Cache-Only Performance Test"
echo "Comparing: Original vs Subtree Cache ONLY (with bugfix)"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "========================================"
echo ""

# Create a temporary test file that only compares Original vs Cache
cat > src/test/java/edu/duke/cs/osprey/markstar/TestCacheOnly.java << 'EOF'
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

public class TestCacheOnly {

    private static final int NUM_CPUs = 20;

    @Test
    public void testOriginalVsSubtreeCacheOnly() {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("CACHE-ONLY PERFORMANCE TEST");
        System.out.println("Comparing: Original (Greedy, No Cache) vs Subtree Cache ONLY");
        System.out.println("Scale: 7 flexible residues");
        System.out.println("=".repeat(100) + "\n");

        int scale = 7;

        // Test 1: Original (NO DP, NO Cache)
        System.out.println("=== TEST 1: Original (Greedy + No Cache) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;  // Disable DP
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = false;  // Disable DP optimizations
        CachedMinimizer.ENABLE_SUBTREE_CACHE = false;  // Disable cache
        CachedMinimizer.clearGlobalCache();
        long timeOriginal = runTest(scale, "Original");

        // Test 2: Subtree Cache ONLY (NO DP, WITH Cache + BUGFIX)
        System.out.println("\n=== TEST 2: Subtree Cache ONLY (Greedy + Cache with Bugfix) ===\n");
        UpdatingEnergyMatrix.USE_DP_ALGORITHM = false;  // Disable DP (use Greedy)
        UpdatingEnergyMatrix.USE_DP_OPTIMIZATIONS = false;  // Disable DP optimizations
        CachedMinimizer.ENABLE_SUBTREE_CACHE = true;  // Enable cache

        TestKStar.ConfSpaces confSpaces = buildConfSpace(scale);
        CachedMinimizer.clearGlobalCache();

        // BUGFIX: Initialize cache for ALL ConfSpaces
        System.out.println("[BUGFIX] Initializing caches for all ConfSpaces...");
        CachedMinimizer.initializeGlobalCache(confSpaces.protein);
        CachedMinimizer.initializeGlobalCache(confSpaces.ligand);
        CachedMinimizer.initializeGlobalCache(confSpaces.complex);
        System.out.println();

        long timeCacheOnly = runTest(scale, "CacheOnly", confSpaces);

        // Print cache statistics
        System.out.println("\n" + "=".repeat(100));
        System.out.println("CACHE STATISTICS");
        System.out.println("=".repeat(100));
        CachedMinimizer.printGlobalStats();

        // Final summary
        System.out.println("\n" + "=".repeat(100));
        System.out.println("FINAL PERFORMANCE SUMMARY - Cache Impact");
        System.out.println("=".repeat(100));
        System.out.println();
        System.out.println(String.format("%-30s: %s", "Original (Greedy, No Cache)", formatTime(timeOriginal)));
        System.out.println(String.format("%-30s: %s", "Subtree Cache ONLY", formatTime(timeCacheOnly)));
        System.out.println();

        double speedup = (double) timeOriginal / timeCacheOnly;
        double improvement = ((double) timeOriginal - timeCacheOnly) / timeOriginal * 100;

        System.out.println(String.format("Speedup: %.2fx", speedup));
        System.out.println(String.format("Performance improvement: %+.1f%%", improvement));
        System.out.println();

        if (speedup > 1.1) {
            System.out.println("✓✓ Subtree cache provides significant speedup!");
        } else if (speedup > 1.0) {
            System.out.println("✓ Subtree cache provides modest speedup");
        } else {
            System.out.println("Cache overhead exceeds benefits on this scale");
        }

        System.out.println();
        System.out.println("NOTES:");
        System.out.println("- Both tests use Greedy algorithm (NO DP optimizations)");
        System.out.println("- Only difference is Subtree DOF Cache enabled/disabled");
        System.out.println("- This isolates the cache performance impact");
        System.out.println("- Bugfix ensures proper cache for each ConfSpace");
        System.out.println("=".repeat(100) + "\n");
    }

    private long runTest(int numFlexible, String label) {
        return runTest(numFlexible, label, null);
    }

    private long runTest(int numFlexible, String label, TestKStar.ConfSpaces confSpaces) {
        if (confSpaces == null) {
            confSpaces = buildConfSpace(numFlexible);
        }

        System.out.println("Configuration (" + label + "):");
        System.out.println("  Flexible positions: " + confSpaces.complex.positions.size());
        System.out.println("  Algorithm: " + (UpdatingEnergyMatrix.USE_DP_ALGORITHM ? "DP" : "Greedy"));
        System.out.println("  Cache enabled: " + CachedMinimizer.ENABLE_SUBTREE_CACHE);
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

    private TestKStar.ConfSpaces buildConfSpace(int numMutable) {
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

        // Add mutable positions
        if (numMutable >= 1) {
            protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 2) {
            protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 3) {
            protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 4) {
            protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 5) {
            ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 6) {
            ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
        }
        if (numMutable >= 7) {
            ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
        }

        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
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
EOF

echo "Temporary test file created"
echo ""

# Compile
echo "Compiling..."
./gradlew compileTestJava

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo ""
echo "========================================"
echo "Running Cache-Only Performance Test"
echo "========================================"
echo ""

# Run the test
./gradlew test --tests "edu.duke.cs.osprey.markstar.TestCacheOnly.testOriginalVsSubtreeCacheOnly"

TEST_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Test completed with exit code: $TEST_EXIT_CODE"
echo "End time: $(date)"
echo "========================================"

# Copy results
if [ -f "build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestCacheOnly.xml" ]; then
    cp build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestCacheOnly.xml \
       /home/users/lz280/"260114 triple correction analysis"/test_cache_only_${SLURM_JOB_ID}.xml
    echo "Test XML results copied"
fi

# Clean up temporary test file
rm -f src/test/java/edu/duke/cs/osprey/markstar/TestCacheOnly.java
echo "Temporary test file removed"

echo ""
echo "========================================"
echo "Job complete!"
echo "========================================"

exit $TEST_EXIT_CODE
