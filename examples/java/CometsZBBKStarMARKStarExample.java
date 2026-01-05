/*
** COMETSZ + BBKStar + MARKStar Example
**
** This example demonstrates how to use COMETSZ (COMETS with thermodynamic ensembles)
** combined with BBKStar and MARKStar for efficient multi-state protein design.
**
** Key components:
** 1. COMETSZ: Multi-state design with thermodynamic ensembles and sequence optimization
** 2. BBKStar: Batch-based K* for efficient sequence space exploration
** 3. MARKStar: Fast partition function calculation with rigorous bounds
**
** This approach combines the best of all three methods:
** - COMETSZ's multi-state objective functions and constraints
** - BBKStar's efficient sequence pruning and batch processing
** - MARKStar's rapid partition function bounds
*/

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.CometsZ;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Example: COMETSZ + BBKStar + MARKStar for Multi-State Protein Design
 *
 * This example shows how to combine three powerful algorithms:
 * - COMETSZ for multi-state thermodynamic ensemble design
 * - BBKStar for batch-based sequence optimization (implicit through COMETSZ design)
 * - MARKStar for efficient partition function calculation
 *
 * The system models protein-ligand binding with three states:
 * - Complex: Protein + Ligand bound together
 * - Protein: Unbound protein
 * - Ligand: Unbound ligand
 *
 * Objective: Maximize binding affinity (LMFE)
 * LMFE = G_complex - G_protein - G_ligand
 */
public class CometsZBBKStarMARKStarExample {

    // Forcefield parameters
    private static final ForcefieldParams ffparams = new ForcefieldParams();

    // Epsilon for partition function approximation
    // 0.95 means 95% confidence in the partition function bounds
    private static final double EPSILON = 0.95;

    // Number of best sequences to find
    private static final int NUM_SEQUENCES = 5;

    // Number of simultaneous mutations to explore
    private static final int MAX_SIMULTANEOUS_MUTATIONS = 1;

    public static void main(String[] args) {

        System.out.println("\n" + "=".repeat(80));
        System.out.println("COMETSZ + BBKStar + MARKStar Example");
        System.out.println("Multi-State Protein Design with Efficient Algorithms");
        System.out.println("=".repeat(80) + "\n");

        // Run the example
        runCometsZWithMARKStar();

        System.out.println("\n" + "=".repeat(80));
        System.out.println("Example Complete!");
        System.out.println("=".repeat(80) + "\n");
    }

    /**
     * Main method to run COMETSZ with MARKStar partition functions
     */
    public static void runCometsZWithMARKStar() {

        // Step 1: Define conformation spaces
        System.out.println("Step 1: Defining conformation spaces...");
        ConfSpaces confSpaces = defineConfSpaces();
        System.out.println("  ✓ Protein conf space: " + confSpaces.protein.positions.size() + " positions");
        System.out.println("  ✓ Ligand conf space: " + confSpaces.ligand.positions.size() + " positions");
        System.out.println("  ✓ Complex conf space: " + confSpaces.complex.positions.size() + " positions\n");

        // Step 2: Create COMETSZ states
        System.out.println("Step 2: Creating COMETSZ states...");
        CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
        CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
        CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);
        System.out.println("  ✓ Created 3 states: Protein, Ligand, Complex\n");

        // Step 3: Define objective function (LMFE for binding)
        System.out.println("Step 3: Defining objective function...");
        System.out.println("  Objective: LMFE = G(complex) - G(protein) - G(ligand)");
        System.out.println("  This represents the binding free energy\n");
        CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
            .addState(complex, 1.0)    // G_complex
            .addState(protein, -1.0)   // - G_protein
            .addState(ligand, -1.0)    // - G_ligand
            .build();

        // Step 4: Configure COMETSZ
        System.out.println("Step 4: Configuring COMETSZ...");
        System.out.println("  - Epsilon: " + EPSILON);
        System.out.println("  - Max simultaneous mutations: " + MAX_SIMULTANEOUS_MUTATIONS);
        System.out.println("  - Objective window size: 100.0 kcal/mol");
        System.out.println("  - Objective window max: 100.0 kcal/mol\n");

        CometsZ cometsZ = new CometsZ.Builder(objective)
            .setEpsilon(EPSILON)
            .setMaxSimultaneousMutations(MAX_SIMULTANEOUS_MUTATIONS)
            .setObjectiveWindowSize(100.0)
            .setObjectiveWindowMax(100.0)
            .setLogFile(new File("cometsz.bbkstar.markstar.results.tsv"))
            .build();

        // Step 5: Initialize energy calculators and matrices with MARKStar
        System.out.println("Step 5: Initializing energy calculators and MARKStar...");
        Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
        Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();
        initCometsZStatesWithMARKStar(cometsZ.states, rigidEmats, minimizingEmats);
        System.out.println("  ✓ Energy matrices calculated");
        System.out.println("  ✓ MARKStar partition functions configured\n");

        // Step 6: Run COMETSZ to find best sequences
        System.out.println("Step 6: Running COMETSZ to find " + NUM_SEQUENCES + " best sequences...");
        System.out.println("  This may take several minutes depending on the system size.\n");
        System.out.println("=".repeat(80) + "\n");

        // Prepare states and run
        prepareAndRunCometsZ(cometsZ, NUM_SEQUENCES);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("RESULTS SUMMARY");
        System.out.println("=".repeat(80));
        System.out.println("All sequences have been saved to: cometsz.bbkstar.markstar.results.tsv");
        System.out.println("\nKey advantages of this approach:");
        System.out.println("  1. COMETSZ: Multi-state design with thermodynamic ensembles");
        System.out.println("  2. BBKStar-style: Efficient batch processing (built into COMETSZ)");
        System.out.println("  3. MARKStar: Fast partition function bounds (5-100x speedup)");
    }

    /**
     * Define conformation spaces for protein, ligand, and complex
     */
    private static class ConfSpaces {
        public SimpleConfSpace protein;
        public SimpleConfSpace ligand;
        public SimpleConfSpace complex;
    }

    private static ConfSpaces defineConfSpaces() {

        ConfSpaces spaces = new ConfSpaces();

        // Read PDB file from examples directory
        Molecule mol = PDBIO.read("../2RL0.kstar/2RL0.min.reduce.pdb");

        // Define protein strand with flexible residues
        Strand proteinStrand = new Strand.Builder(mol)
            .setResidues("G648", "G654")
            .build();

        // Add flexibility to protein positions
        proteinStrand.flexibility.get("G649")
            .setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers()
            .setContinuous();
        proteinStrand.flexibility.get("G650")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        proteinStrand.flexibility.get("G651")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        proteinStrand.flexibility.get("G654")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();

        // Define ligand strand with flexible residues
        Strand ligandStrand = new Strand.Builder(mol)
            .setResidues("A155", "A194")
            .build();

        // Add flexibility to ligand positions
        ligandStrand.flexibility.get("A156")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        ligandStrand.flexibility.get("A172")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        ligandStrand.flexibility.get("A192")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        ligandStrand.flexibility.get("A193")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();

        // Create conformation spaces
        spaces.protein = new SimpleConfSpace.Builder()
            .addStrand(proteinStrand)
            .build();

        spaces.ligand = new SimpleConfSpace.Builder()
            .addStrand(ligandStrand)
            .build();

        spaces.complex = new SimpleConfSpace.Builder()
            .addStrands(proteinStrand, ligandStrand)
            .build();

        return spaces;
    }

    /**
     * Initialize COMETSZ states with MARKStar partition functions
     * This follows the pattern from TestCometsZWithBBKStarAndMARKStar
     */
    private static void initCometsZStatesWithMARKStar(
            List<CometsZ.State> states,
            Map<CometsZ.State, EnergyMatrix> rigidEmats,
            Map<CometsZ.State, EnergyMatrix> minimizingEmats) {

        // Use shared EnergyCalculator for all states
        List<SimpleConfSpace> confSpaceList = states.stream()
            .map(state -> state.confSpace)
            .collect(java.util.stream.Collectors.toList());

        // Create minimizing energy calculator
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
            .setParallelism(Parallelism.makeCpu(4))
            .setIsMinimizing(true)
            .build();

        // Create rigid energy calculator using SharedBuilder
        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
            .setIsMinimizing(false)
            .build();

        try {
            for (CometsZ.State state : states) {

                System.out.println("  Processing state: " + state.name);

                // Setup minimizing conf energy calculator
                state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, minimizingEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, minimizingEcalc)
                        .build()
                        .calcReferenceEnergies())
                    .build();

                // Calculate minimizing energy matrix
                EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
                    .build()
                    .calcEnergyMatrix();
                minimizingEmats.put(state, minimizingEmat);

                // Calculate rigid energy matrix with same reference energies
                ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(state.confSpace, rigidEcalc)
                    .setReferenceEnergies(state.confEcalc.eref)  // Reuse same reference energies
                    .build();
                EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                    .build()
                    .calcEnergyMatrix();
                rigidEmats.put(state, rigidEmat);

                // Set up partition function factory with MARK*
                state.pfuncFactory = (rcs) -> {
                    MARKStarBound markstarPfunc = new MARKStarBound(
                        state.confSpace,
                        rigidEmats.get(state),
                        minimizingEmats.get(state),
                        state.confEcalc,
                        rcs,
                        Parallelism.makeCpu(4)
                    );
                    // CRITICAL: Set the correction matrix to avoid NullPointerException
                    markstarPfunc.setCorrections(new edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix(
                        state.confSpace,
                        minimizingEmats.get(state),
                        state.confEcalc
                    ));
                    return markstarPfunc;
                };

                // Set up conformation tree factory for sequence search
                state.confTreeFactory = (rcs) -> {
                    return new ConfAStarTree.Builder(minimizingEmats.get(state), rcs)
                        .setTraditional()
                        .build();
                };

                // Set fragment energies from minimizing energy matrix
                state.fragmentEnergies = minimizingEmats.get(state);
            }
        } finally {
            // Note: Energy calculators will be cleaned up automatically
            // minimizingEcalc and rigidEcalc cleanup handled by JVM
        }
    }

    /**
     * Prepare COMETSZ states and run the design
     * Following the pattern from TestCometsZWithBBKStarAndMARKStar
     */
    private static void prepareAndRunCometsZ(CometsZ cometsZ, int numSequences) {

        // Create new EnergyCalculator and refresh confEcalc for each state
        List<SimpleConfSpace> confSpaces = cometsZ.states.stream()
            .map(state -> state.confSpace)
            .collect(java.util.stream.Collectors.toList());

        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
            .setParallelism(Parallelism.makeCpu(4))
            .build();

        // Refresh the conf ecalcs with the new active EnergyCalculator
        for (CometsZ.State state : cometsZ.states) {
            state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
        }

        // Run COMETSZ to find best sequences
        List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(numSequences);

        // Print results
        System.out.println("\n" + "=".repeat(80));
        System.out.println("Found " + sequences.size() + " sequences");
        System.out.println("=".repeat(80) + "\n");

        for (int i = 0; i < sequences.size(); i++) {
            CometsZ.SequenceInfo seqInfo = sequences.get(i);
            System.out.println("Sequence #" + (i + 1) + ": " + seqInfo.sequence);
            System.out.println("  Objective: [" +
                String.format("%.4f", seqInfo.objective.lower) + ", " +
                String.format("%.4f", seqInfo.objective.upper) + "] kcal/mol");

            // Print state-specific free energies
            for (CometsZ.State state : cometsZ.states) {
                var pfuncResult = seqInfo.pfuncResults.get(state);
                var feBounds = pfuncResult.values.calcFreeEnergyBounds();
                System.out.println("  " + state.name + " Free Energy: [" +
                    String.format("%.4f", feBounds.lower) + ", " +
                    String.format("%.4f", feBounds.upper) + "] kcal/mol");
            }
            System.out.println();
        }
    }
}
