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
 * Simplified test to verify if assignTemplate() modifies coordinates
 */
public class TestAssignTemplateSimple {

    @Test
    public void testOriginalMolNotModified() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("TEST: Does Strand.Builder modify original molecule?");
        System.out.println("=".repeat(80) + "\n");

        // Read PDB
        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        System.out.println("Loaded molecule with " + mol.residues.size() + " residues");

        // Save original coordinates of first residue
        Residue firstRes = mol.residues.get(0);
        double[] originalCoords = Arrays.copyOf(firstRes.coords, firstRes.coords.length);
        String resName = firstRes.getPDBResNumber();

        System.out.println("\nOriginal first residue " + resName + " coords[0-8]:");
        for (int i = 0; i < Math.min(9, originalCoords.length); i++) {
            System.out.printf("%.3f ", originalCoords[i]);
        }
        System.out.println();

        // Create template library
        ForcefieldParams ffparams = new ForcefieldParams();
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // Create protein strand (this calls assignTemplate internally)
        System.out.println("\nCreating protein strand...");
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();
        System.out.println("Protein strand created");

        // Check if original mol was modified
        System.out.println("\nChecking if original mol first residue was modified...");
        if (Arrays.equals(originalCoords, firstRes.coords)) {
            System.out.println("✓ GOOD: Original mol was NOT modified");
        } else {
            System.err.println("❌ BAD: Original mol WAS modified!");
            System.out.print("Current coords[0-8]: ");
            for (int i = 0; i < Math.min(9, firstRes.coords.length); i++) {
                System.out.printf("%.3f ", firstRes.coords[i]);
            }
            System.out.println();
        }

        System.out.println("\n" + "=".repeat(80) + "\n");
    }

    @Test
    public void testStrandConsistency() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("TEST: Are strands created consistently?");
        System.out.println("=".repeat(80) + "\n");

        ForcefieldParams ffparams = new ForcefieldParams();
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // Create protein strand from mol1
        Molecule mol1 = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        Strand protein1 = new Strand.Builder(mol1)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        // Create protein strand from mol2 (fresh copy)
        Molecule mol2 = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));
        Strand protein2 = new Strand.Builder(mol2)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();

        // Compute hashes for both
        int hash1 = computeMoleculeHash(protein1.mol);
        int hash2 = computeMoleculeHash(protein2.mol);

        System.out.println("Protein1 hash: " + hash1);
        System.out.println("Protein2 hash: " + hash2);

        if (hash1 == hash2) {
            System.out.println("\n✓ GOOD: Same strand created twice has same hash");
        } else {
            System.err.println("\n❌ BAD: Same strand created twice has DIFFERENT hash!");
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
