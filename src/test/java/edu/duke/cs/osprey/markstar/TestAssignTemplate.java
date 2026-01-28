package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

/**
 * Test to verify if assignTemplate() modifies coordinates
 */
public class TestAssignTemplate {

    @Test
    public void testAssignTemplateModifiesCoordinates() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("TEST: Does assignTemplate() modify coordinates?");
        System.out.println("=".repeat(80) + "\n");

        // Read PDB
        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        System.out.println("Loaded molecule with " + mol.residues.size() + " residues");

        // Save original coordinates of first few residues
        int numResiduesToCheck = Math.min(5, mol.residues.size());
        double[][] originalCoords = new double[numResiduesToCheck][];
        String[] resNames = new String[numResiduesToCheck];

        System.out.println("\n--- Original coordinates (before any Strand creation) ---");
        for (int i = 0; i < numResiduesToCheck; i++) {
            Residue res = mol.residues.get(i);
            if (res.coords != null) {
                originalCoords[i] = Arrays.copyOf(res.coords, res.coords.length);
                resNames[i] = res.getPDBResNumber();

                // Print first 9 values (3 atoms)
                int numToPrint = Math.min(9, res.coords.length);
                System.out.print("Residue " + resNames[i] + " coords[0-" + (numToPrint-1) + "]: ");
                for (int j = 0; j < numToPrint; j++) {
                    System.out.printf("%.3f ", res.coords[j]);
                }
                System.out.println();
            }
        }

        // Create template library
        ForcefieldParams ffparams = new ForcefieldParams();
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // Create protein strand (this calls assignTemplate internally)
        System.out.println("\n--- Creating protein strand ---");
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();
        System.out.println("Protein strand created with " + protein.mol.residues.size() + " residues");

        // Check if original mol was modified
        System.out.println("\n--- Checking if original mol was modified ---");
        boolean originalModified = false;
        for (int i = 0; i < numResiduesToCheck; i++) {
            Residue res = mol.residues.get(i);
            if (res.coords != null && originalCoords[i] != null) {
                if (!Arrays.equals(originalCoords[i], res.coords)) {
                    System.err.println("⚠️  Original mol residue " + resNames[i] + " WAS MODIFIED!");
                    originalModified = true;

                    // Show what changed
                    int numToPrint = Math.min(9, res.coords.length);
                    System.out.print("  Original: ");
                    for (int j = 0; j < numToPrint; j++) {
                        System.out.printf("%.3f ", originalCoords[i][j]);
                    }
                    System.out.println();
                    System.out.print("  Current:  ");
                    for (int j = 0; j < numToPrint; j++) {
                        System.out.printf("%.3f ", res.coords[j]);
                    }
                    System.out.println();
                } else {
                    System.out.println("✓ Original mol residue " + resNames[i] + " unchanged");
                }
            }
        }

        if (!originalModified) {
            System.out.println("\n✓ GOOD: Original mol was NOT modified");
        } else {
            System.err.println("\n❌ BAD: Original mol WAS modified!");
        }

        // Now compare protein.mol coordinates with original
        System.out.println("\n--- Comparing protein.mol with original mol ---");
        for (int i = 0; i < protein.mol.residues.size(); i++) {
            Residue proteinRes = protein.mol.residues.get(i);
            String proteinResNum = proteinRes.getPDBResNumber();

            // Find corresponding residue in original mol
            Residue originalRes = mol.getResByPDBResNumber(proteinResNum);

            if (originalRes != null && originalRes.coords != null && proteinRes.coords != null) {
                // Check if same length
                if (originalRes.coords.length != proteinRes.coords.length) {
                    System.out.println("⚠️  Residue " + proteinResNum + " has different number of coords");
                    System.out.println("    Original: " + originalRes.coords.length + ", Protein: " + proteinRes.coords.length);
                    continue;
                }

                // Check if coordinates are the same (but possibly reordered)
                boolean identical = Arrays.equals(originalRes.coords, proteinRes.coords);

                if (identical) {
                    System.out.println("✓ Residue " + proteinResNum + " coords IDENTICAL (same order)");
                } else {
                    // Check if values are same but reordered
                    boolean sameValues = haveSameValues(originalRes.coords, proteinRes.coords);
                    if (sameValues) {
                        System.out.println("⚠️  Residue " + proteinResNum + " coords REORDERED (same values, different order)");
                    } else {
                        System.err.println("❌ Residue " + proteinResNum + " coords VALUES CHANGED!");

                        // Show differences
                        int numToPrint = Math.min(9, originalRes.coords.length);
                        System.out.print("  Original: ");
                        for (int j = 0; j < numToPrint; j++) {
                            System.out.printf("%.3f ", originalRes.coords[j]);
                        }
                        System.out.println();
                        System.out.print("  Protein:  ");
                        for (int j = 0; j < numToPrint; j++) {
                            System.out.printf("%.3f ", proteinRes.coords[j]);
                        }
                        System.out.println();
                    }
                }
            }
        }

        // Create ligand strand
        System.out.println("\n--- Creating ligand strand ---");
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();
        System.out.println("Ligand strand created with " + ligand.mol.residues.size() + " residues");

        // Compare a few ligand residues
        System.out.println("\n--- Comparing ligand.mol with original mol ---");
        int numLigandResToCheck = Math.min(3, ligand.mol.residues.size());
        for (int i = 0; i < numLigandResToCheck; i++) {
            Residue ligandRes = ligand.mol.residues.get(i);
            String ligandResNum = ligandRes.getPDBResNumber();

            Residue originalRes = mol.getResByPDBResNumber(ligandResNum);

            if (originalRes != null && originalRes.coords != null && ligandRes.coords != null) {
                boolean identical = Arrays.equals(originalRes.coords, ligandRes.coords);

                if (identical) {
                    System.out.println("✓ Residue " + ligandResNum + " coords IDENTICAL");
                } else {
                    boolean sameValues = haveSameValues(originalRes.coords, ligandRes.coords);
                    if (sameValues) {
                        System.out.println("⚠️  Residue " + ligandResNum + " coords REORDERED");
                    } else {
                        System.err.println("❌ Residue " + ligandResNum + " coords VALUES CHANGED!");
                    }
                }
            }
        }

        System.out.println("\n" + "=".repeat(80));
        System.out.println("TEST COMPLETE");
        System.out.println("=".repeat(80) + "\n");
    }

    /**
     * Check if two arrays have the same values (possibly in different order)
     */
    private boolean haveSameValues(double[] arr1, double[] arr2) {
        if (arr1.length != arr2.length) return false;

        double[] sorted1 = Arrays.copyOf(arr1, arr1.length);
        double[] sorted2 = Arrays.copyOf(arr2, arr2.length);
        Arrays.sort(sorted1);
        Arrays.sort(sorted2);

        for (int i = 0; i < sorted1.length; i++) {
            if (Math.abs(sorted1[i] - sorted2[i]) > 1e-6) {
                return false;
            }
        }
        return true;
    }

    @Test
    public void testTemplateHashConsistency() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("TEST: Template Hash Consistency");
        System.out.println("=".repeat(80) + "\n");

        // Read PDB
        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

        ForcefieldParams ffparams = new ForcefieldParams();
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // Create protein strand
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        // Create ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();

        // Compute hash for protein.mol
        int proteinHash = computeMoleculeHash(protein.mol);
        System.out.println("Protein.mol hash: " + proteinHash);

        // Compute hash for ligand.mol
        int ligandHash = computeMoleculeHash(ligand.mol);
        System.out.println("Ligand.mol hash:  " + ligandHash);

        // Create another protein strand from fresh mol
        Molecule mol2 = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        Strand protein2 = new Strand.Builder(mol2)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        int protein2Hash = computeMoleculeHash(protein2.mol);
        System.out.println("Protein2.mol hash: " + protein2Hash);

        if (proteinHash == protein2Hash) {
            System.out.println("\n✓ GOOD: Creating protein strand twice gives same hash");
        } else {
            System.err.println("\n❌ BAD: Creating protein strand twice gives DIFFERENT hash!");
        }

        System.out.println("\n" + "=".repeat(80) + "\n");
    }

    private int computeMoleculeHash(Molecule mol) {
        int hash = 0;
        for (Residue res : mol.residues) {
            if (res.coords != null) {
                int n = Math.min(9, res.coords.length);
                for (int i = 0; i < n; i++) {
                    hash = 31 * hash + Double.hashCode(res.coords[i]);
                }
            }
        }
        return hash;
    }
}
