/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
**
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
**
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
**
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
**
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.tools.Log.log;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.util.*;

/**
 * Test to compare MARK* vs K* partition function calculations.
 * This test checks:
 * 1. Whether MARK* produces smaller (tighter) partition function values than K*
 * 2. The relationship between score, minimized energy, lower bound, and upper bound
 * 3. Tests on multiple small systems with 3-5 flexible residues
 */
public class TestMARKStarVsKStarPartitionFunction {

    private static ForcefieldParams ffparams = new ForcefieldParams();

    /**
     * Result holder for partition function comparison
     */
    private static class PartitionFunctionResult {
        String algorithm;
        BigDecimal qstar;  // lower bound
        BigDecimal pstar;  // upper bound (or qprime for K*)
        BigDecimal freeEnergyLower;
        BigDecimal freeEnergyUpper;
        int numConfsEvaluated;
        long timeMs;

        // Energy bound information for individual conformations
        List<ConformationBounds> confBounds = new ArrayList<>();

        PartitionFunctionResult(String algorithm) {
            this.algorithm = algorithm;
        }

        @Override
        public String toString() {
            return String.format(
                "%s:\n" +
                "  Q* (lower):           %12.6e\n" +
                "  P* (upper):           %12.6e\n" +
                "  Free Energy Lower:    %12.6f\n" +
                "  Free Energy Upper:    %12.6f\n" +
                "  Confs Evaluated:      %d\n" +
                "  Time (ms):            %d",
                algorithm,
                qstar, pstar,
                freeEnergyLower.doubleValue(), freeEnergyUpper.doubleValue(),
                numConfsEvaluated, timeMs
            );
        }
    }

    /**
     * Holds bounds information for a single conformation
     */
    private static class ConformationBounds {
        int[] conf;
        double score;           // A* score (lower bound before minimization)
        double minimizedEnergy; // Actual minimized energy
        double lowerBound;      // Lower bound used in pfunc
        double upperBound;      // Upper bound used in pfunc

        ConformationBounds(int[] conf, double score, double minimizedEnergy,
                          double lowerBound, double upperBound) {
            this.conf = conf.clone();
            this.score = score;
            this.minimizedEnergy = minimizedEnergy;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }

        /**
         * Validate that bounds are consistent:
         * lowerBound <= minimizedEnergy <= upperBound
         * score <= minimizedEnergy (score is pre-minimization lower bound)
         */
        boolean areBoundsConsistent() {
            boolean consistent = true;
            StringBuilder errors = new StringBuilder();

            if (lowerBound > minimizedEnergy + 1e-6) {
                consistent = false;
                errors.append(String.format("  Lower bound (%f) > minimized energy (%f)\n",
                    lowerBound, minimizedEnergy));
            }

            if (minimizedEnergy > upperBound + 1e-6) {
                consistent = false;
                errors.append(String.format("  Minimized energy (%f) > upper bound (%f)\n",
                    minimizedEnergy, upperBound));
            }

            if (score > minimizedEnergy + 1e-6) {
                consistent = false;
                errors.append(String.format("  Score (%f) > minimized energy (%f)\n",
                    score, minimizedEnergy));
            }

            if (!consistent) {
                log("Bounds inconsistency for conf %s:\n%s",
                    SimpleConfSpace.formatConfRCs(conf), errors.toString());
            }

            return consistent;
        }
    }

    /**
     * Create a small test system with 3 flexible residues
     */
    private static SimpleConfSpace makeSmallSystem3Flex() {
        Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

        Strand protein = new Strand.Builder(mol)
            .setResidues("G648", "G654")
            .build();

        // Only 3 flexible positions
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA", "VAL").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        return new SimpleConfSpace.Builder()
            .addStrand(protein)
            .build();
    }

    /**
     * Create a small test system with 4 flexible residues
     */
    private static SimpleConfSpace makeSmallSystem4Flex() {
        Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

        Strand protein = new Strand.Builder(mol)
            .setResidues("G648", "G654")
            .build();

        // 4 flexible positions
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER").addWildTypeRotamers().setContinuous();

        return new SimpleConfSpace.Builder()
            .addStrand(protein)
            .build();
    }

    /**
     * Create a small test system with 5 flexible residues
     */
    private static SimpleConfSpace makeSmallSystem5Flex() {
        Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

        Strand ligand = new Strand.Builder(mol)
            .setResidues("A155", "A194")
            .build();

        // 5 flexible positions
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "VAL").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        return new SimpleConfSpace.Builder()
            .addStrand(ligand)
            .build();
    }

    /**
     * Calculate partition function using standard K* gradient descent
     */
    private static PartitionFunctionResult calculateWithKStar(
            SimpleConfSpace confSpace,
            ConfEnergyCalculator confEcalc,
            EnergyMatrix emat,
            double epsilon) {

        PartitionFunctionResult result = new PartitionFunctionResult("K* GradientDescent");

        RCs rcs = new RCs(confSpace);

        GradientDescentPfunc pfunc = new GradientDescentPfunc(
            confEcalc,
            new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
            new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
            rcs.getNumConformations()
        );

        try (TaskExecutor.ContextGroup contexts = confEcalc.ecalc.tasks.contextGroup()) {
            pfunc.setInstanceId(0);
            pfunc.putTaskContexts(contexts);
            pfunc.setReportProgress(true);

            long startTime = System.currentTimeMillis();
            pfunc.init(epsilon);
            pfunc.compute();
            result.timeMs = System.currentTimeMillis() - startTime;

            PartitionFunction.Result pfuncResult = pfunc.makeResult();
            result.qstar = pfuncResult.values.qstar;
            result.pstar = pfuncResult.values.qprime.add(pfuncResult.values.qstar); // K* uses qprime as the gap
            result.numConfsEvaluated = pfuncResult.numConfs;

            MathTools.DoubleBounds feBounds = pfuncResult.values.calcFreeEnergyBounds();
            result.freeEnergyLower = new BigDecimal(feBounds.lower);
            result.freeEnergyUpper = new BigDecimal(feBounds.upper);
        }

        return result;
    }

    /**
     * Calculate partition function using MARK*
     */
    private static PartitionFunctionResult calculateWithMARKStar(
            SimpleConfSpace confSpace,
            ConfEnergyCalculator confEcalc,
            EnergyMatrix rigidEmat,
            EnergyMatrix minimizingEmat,
            double epsilon) {

        PartitionFunctionResult result = new PartitionFunctionResult("MARK*");

        RCs rcs = new RCs(confSpace);

        MARKStarBound markstar = new MARKStarBound(
            confSpace,
            rigidEmat,
            minimizingEmat,
            confEcalc,
            rcs,
            Parallelism.makeCpu(4)
        );

        // Initialize correction matrix for MARK*
        UpdatingEnergyMatrix correctionMatrix = new UpdatingEnergyMatrix(confSpace, rigidEmat);
        markstar.setCorrections(correctionMatrix);

        long startTime = System.currentTimeMillis();
        markstar.init(epsilon);
        markstar.setReportProgress(true);
        markstar.compute();
        result.timeMs = System.currentTimeMillis() - startTime;

        PartitionFunction.Result pfuncResult = markstar.makeResult();
        result.qstar = pfuncResult.values.qstar;
        result.pstar = pfuncResult.values.pstar;
        result.numConfsEvaluated = pfuncResult.numConfs;

        MathTools.DoubleBounds feBounds = pfuncResult.values.calcFreeEnergyBounds();
        result.freeEnergyLower = new BigDecimal(feBounds.lower);
        result.freeEnergyUpper = new BigDecimal(feBounds.upper);

        return result;
    }

    /**
     * Test partition function comparison on a small system
     */
    private static void testSystemComparison(String systemName, SimpleConfSpace confSpace, double epsilon) {
        log("\n========================================");
        log("Testing: %s", systemName);
        log("Number of flexible residues: %d", confSpace.positions.size());
        log("Epsilon: %.4f", epsilon);
        log("========================================\n");

        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
            .setParallelism(Parallelism.makeCpu(8))
            .build()) {

            // Setup energy calculator
            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                    .build()
                    .calcReferenceEnergies()
                )
                .build();

            // Calculate energy matrices
            log("Calculating minimizing energy matrix...");
            EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                .build()
                .calcEnergyMatrix();

            log("Calculating rigid energy matrix...");
            EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                .build()
                .calcEnergyMatrix();

            // Test with K*
            log("\nRunning K* GradientDescent...");
            PartitionFunctionResult kstarResult = calculateWithKStar(confSpace, confEcalc, minimizingEmat, epsilon);
            log(kstarResult.toString());

            // Test with MARK*
            log("\nRunning MARK*...");
            PartitionFunctionResult markstarResult = calculateWithMARKStar(
                confSpace, confEcalc, rigidEmat, minimizingEmat, epsilon);
            log(markstarResult.toString());

            // Compare results
            log("\n========================================");
            log("COMPARISON RESULTS:");
            log("========================================");

            // Compare partition function bounds
            double qstarRatio = markstarResult.qstar.divide(
                kstarResult.qstar, PartitionFunction.decimalPrecision).doubleValue();
            double pstarRatio = markstarResult.pstar.divide(
                kstarResult.pstar, PartitionFunction.decimalPrecision).doubleValue();

            log("Q* (lower bound) ratio (MARK*/K*): %.6f", qstarRatio);
            log("P* (upper bound) ratio (MARK*/K*): %.6f", pstarRatio);

            // Check if MARK* produces tighter bounds
            BigDecimal kstarGap = kstarResult.pstar.subtract(kstarResult.qstar);
            BigDecimal markstarGap = markstarResult.pstar.subtract(markstarResult.qstar);
            double gapRatio = markstarGap.divide(kstarGap, PartitionFunction.decimalPrecision).doubleValue();

            log("Gap between bounds (P* - Q*):");
            log("  K*:     %12.6e", kstarGap);
            log("  MARK*:  %12.6e", markstarGap);
            log("  Ratio (MARK*/K*): %.6f", gapRatio);

            if (gapRatio < 1.0) {
                log("✓ MARK* produces tighter bounds than K* (%.2f%% of K* gap)", gapRatio * 100);
            } else {
                log("✗ MARK* gap is larger than K* gap");
            }

            // Compare free energies
            log("\nFree Energy Comparison:");
            log("  K* range:     [%.6f, %.6f]",
                kstarResult.freeEnergyLower.doubleValue(),
                kstarResult.freeEnergyUpper.doubleValue());
            log("  MARK* range:  [%.6f, %.6f]",
                markstarResult.freeEnergyLower.doubleValue(),
                markstarResult.freeEnergyUpper.doubleValue());

            // Compare performance
            log("\nPerformance:");
            log("  K* evaluated %d conformations in %d ms",
                kstarResult.numConfsEvaluated, kstarResult.timeMs);
            log("  MARK* evaluated %d conformations in %d ms",
                markstarResult.numConfsEvaluated, markstarResult.timeMs);
            double speedup = (double) kstarResult.timeMs / markstarResult.timeMs;
            log("  Speedup: %.2fx", speedup);

            log("========================================\n");
        }
    }

    @Test
    public void testSmallSystem3Flex() {
        SimpleConfSpace confSpace = makeSmallSystem3Flex();
        testSystemComparison("Small System (3 flexible residues)", confSpace, 0.10);
    }

    @Test
    public void testSmallSystem4Flex() {
        SimpleConfSpace confSpace = makeSmallSystem4Flex();
        testSystemComparison("Small System (4 flexible residues)", confSpace, 0.10);
    }

    @Test
    public void testSmallSystem5Flex() {
        SimpleConfSpace confSpace = makeSmallSystem5Flex();
        testSystemComparison("Small System (5 flexible residues)", confSpace, 0.10);
    }

    /**
     * Test all three systems and summarize results
     */
    @Test
    public void testAllSmallSystems() {
        log("\n");
        log("================================================================================");
        log("MARK* vs K* PARTITION FUNCTION COMPARISON TEST");
        log("Testing on multiple small systems (3-5 flexible residues)");
        log("================================================================================\n");

        testSystemComparison("System 1: 3 flexible residues", makeSmallSystem3Flex(), 0.10);
        testSystemComparison("System 2: 4 flexible residues", makeSmallSystem4Flex(), 0.10);
        testSystemComparison("System 3: 5 flexible residues", makeSmallSystem5Flex(), 0.10);

        log("\n");
        log("================================================================================");
        log("ALL TESTS COMPLETED");
        log("================================================================================\n");
    }
}
