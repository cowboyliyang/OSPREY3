/*
** TRUE BBK* + MARK* Example
**
** This example correctly demonstrates BBK* (not COMETSZ) with MARK* partition functions.
** Based on TestBBKStar.java pattern.
**
** Key differences from previous example:
** 1. Uses BBKStar class directly (not CometsZ)
** 2. Configures BOTH rigid and minimizing energy matrices
** 3. Properly integrates MARK* as the partition function method
*/

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.List;

/**
 * TRUE Example: BBK* + MARK* Integration
 *
 * This correctly shows how to use BBKStar with MARK* partition functions.
 * Pattern follows TestBBKStar.runBBKStar() method.
 */
public class TrueBBKStarMARKStarExample {

    private static final ForcefieldParams ffparams = new ForcefieldParams();

    public static void main(String[] args) {

        System.out.println("\n" + "=".repeat(80));
        System.out.println("TRUE BBK* + MARK* Example");
        System.out.println("Following TestBBKStar.java pattern");
        System.out.println("=".repeat(80) + "\n");

        // Define the system (protein-ligand binding)
        System.out.println("Step 1: Loading molecule and defining conformation spaces...");
        ConfSpaces confSpaces = defineConfSpaces();
        System.out.println("  ✓ Protein positions: " + confSpaces.protein.positions.size());
        System.out.println("  ✓ Ligand positions: " + confSpaces.ligand.positions.size());
        System.out.println("  ✓ Complex positions: " + confSpaces.complex.positions.size() + "\n");

        // Run BBK* with MARK*
        System.out.println("Step 2: Running BBK* with MARK* partition functions...\n");
        int numSequences = 5;
        double epsilon = 0.95;
        List<KStar.ScoredSequence> results = runBBKStarWithMARKStar(confSpaces, numSequences, epsilon);

        // Display results
        System.out.println("\n" + "=".repeat(80));
        System.out.println("Results: Top " + results.size() + " sequences found");
        System.out.println("=".repeat(80) + "\n");

        for (int i = 0; i < results.size(); i++) {
            KStar.ScoredSequence scoredSeq = results.get(i);
            System.out.println(String.format("#%d: %s", i + 1,
                scoredSeq.sequence.toString(Sequence.Renderer.ResType)));
            if (scoredSeq.score != null && scoredSeq.score.score != null) {
                System.out.println(String.format("    K* score: %.6e", scoredSeq.score.score));
                System.out.println(String.format("    K* bounds: [%.2f, %.2f] log10",
                    scoredSeq.score.lowerBoundLog10(), scoredSeq.score.upperBoundLog10()));
            }
            System.out.println();
        }

        System.out.println("=".repeat(80));
        System.out.println("✓ Example Complete!");
        System.out.println("=".repeat(80) + "\n");
    }

    /**
     * Run BBK* with MARK* partition functions
     * Pattern: TestBBKStar.runBBKStar() lines 67-169
     */
    public static List<KStar.ScoredSequence> runBBKStarWithMARKStar(
            ConfSpaces confSpaces, int numSequences, double epsilon) {

        Parallelism parallelism = Parallelism.makeCpu(4);

        // Create minimizing energy calculator (shared across all states)
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, ffparams)
            .setParallelism(parallelism)
            .build();

        try {
            // ===== CRITICAL: Configure K* and BBK* Settings =====
            // This is what makes it BBK*, not regular K*

            KStar.Settings kstarSettings = new KStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setStabilityThreshold(null)
                .setMaxSimultaneousMutations(1)
                .addScoreConsoleWriter(new KStarScoreWriter.Formatter() {
                    @Override
                    public String format(KStarScoreWriter.ScoreInfo info) {
                        return String.format("Sequence: %s  K* score: %.6e",
                            info.sequence.toString(Sequence.Renderer.ResType),
                            info.kstarScore.score);
                    }
                })
                .build();

            // BBK* specific settings - this enables batch processing!
            BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
                .setNumBestSequences(numSequences)
                .setNumConfsPerBatch(8)  // Process 8 conformations per batch
                .build();

            // ===== Create BBKStar instance (NOT KStar, NOT CometsZ) =====
            BBKStar bbkstar = new BBKStar(
                confSpaces.protein,
                confSpaces.ligand,
                confSpaces.complex,
                kstarSettings,
                bbkstarSettings
            );

            // ===== Configure each conformation space =====
            // This is the critical part that integrates MARK*

            for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {
                SimpleConfSpace confSpace = (SimpleConfSpace)info.confSpace;

                System.out.println("Configuring " + info.type + " state...");

                // 1. Setup minimizing conf energy calculator with reference energies
                info.confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                    ).build();

                // 2. Calculate minimizing energy matrix
                System.out.println("  - Calculating minimizing energy matrix...");
                EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(info.confEcalcMinimized)
                    .build()
                    .calcEnergyMatrix();

                info.confSearchFactoryMinimized = (rcs) ->
                    new ConfAStarTree.Builder(ematMinimized, rcs)
                        .setTraditional()
                        .build();

                // 3. BBK* CRITICAL: Setup rigid energies (BBK* needs both rigid and minimizing!)
                System.out.println("  - Calculating rigid energy matrix...");
                EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false)
                    .build();

                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(info.confEcalcMinimized, ecalcRigid);

                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                    .build()
                    .calcEnergyMatrix();

                info.confSearchFactoryRigid = (rcs) ->
                    new ConfAStarTree.Builder(ematRigid, rcs)
                        .setTraditional()
                        .build();

                // 4. MARK* Integration: Use MARK* as the partition function method
                System.out.println("  - Setting up MARK* partition function...");

                // Store the matrices for use in the lambda
                final EnergyMatrix finalEmatRigid = ematRigid;
                final EnergyMatrix finalEmatMinimized = ematMinimized;

                info.pfuncFactory = rcs -> {
                    // Create MARK* partition function with fast queues
                    MARKStarBound pfunc = new MARKStarBoundFastQueues(
                        confSpace,
                        finalEmatRigid,          // Rigid matrix for upper bounds
                        finalEmatMinimized,      // Minimizing matrix for lower bounds
                        info.confEcalcMinimized,
                        rcs,
                        parallelism
                    );

                    // CRITICAL: Set correction matrix to avoid NullPointerException
                    pfunc.setCorrections(new UpdatingEnergyMatrix(
                        confSpace,
                        finalEmatMinimized,
                        info.confEcalcMinimized
                    ));

                    // Initialize with epsilon
                    pfunc.init(epsilon);

                    return pfunc;
                };

                // No ConfDB for this example
                info.confDBFile = null;

                System.out.println("  ✓ " + info.type + " configuration complete\n");
            }

            // ===== RUN BBK* =====
            System.out.println("Starting BBK* + MARK* computation...");
            System.out.println("This may take several minutes...\n");

            List<KStar.ScoredSequence> sequences = bbkstar.run(ecalcMinimized.tasks);

            System.out.println("\n✓ BBK* computation complete!");

            return sequences;

        } finally {
            ecalcMinimized.close();
        }
    }

    /**
     * Define conformation spaces for protein-ligand binding
     * Using the 2RL0 test system
     */
    public static class ConfSpaces {
        public SimpleConfSpace protein;
        public SimpleConfSpace ligand;
        public SimpleConfSpace complex;
    }

    public static ConfSpaces defineConfSpaces() {
        ConfSpaces spaces = new ConfSpaces();

        // Load PDB structure
        Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

        // Define protein strand
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(ffparams.templLib)
            .setResidues("G648", "G654")
            .build();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Define ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(ffparams.templLib)
            .setResidues("A155", "A194")
            .build();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Create conformation spaces
        spaces.protein = new SimpleConfSpace.Builder()
            .addStrand(protein)
            .build();

        spaces.ligand = new SimpleConfSpace.Builder()
            .addStrand(ligand)
            .build();

        spaces.complex = new SimpleConfSpace.Builder()
            .addStrands(protein, ligand)
            .build();

        return spaces;
    }
}
