package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

/**
 * Test to identify the template hashes for protein, ligand, and complex ConfSpaces
 */
public class TestThreeConfSpaceHashes {

    @Test
    public void testIdentifyConfSpaceHashes() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("IDENTIFYING TEMPLATE HASHES FOR PROTEIN/LIGAND/COMPLEX");
        System.out.println("=".repeat(80) + "\n");

        // Load molecule
        Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

        // Create template library
        ForcefieldParams ffparams = new ForcefieldParams();
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // Create protein strand
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("G648", "G654")
            .build();
        protein.flexibility.get("G648").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Create ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();
        ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // Build ConfSpaces
        SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder().addStrand(protein).build();
        SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder().addStrand(ligand).build();
        SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        // Compute hashes
        int proteinHash = computeTemplateHash(proteinConfSpace);
        int ligandHash = computeTemplateHash(ligandConfSpace);
        int complexHash = computeTemplateHash(complexConfSpace);

        System.out.println("Protein  ConfSpace template hash: " + proteinHash);
        System.out.println("Ligand   ConfSpace template hash: " + ligandHash);
        System.out.println("Complex  ConfSpace template hash: " + complexHash);

        System.out.println("\nExpected hashes from previous test:");
        System.out.println("  - 1310777396  (appeared 3 times, energy=-1.889)");
        System.out.println("  - 1320578138  (unknown frequency)");
        System.out.println("  - -451549042  (appeared 3 times, energy=-11.463)");

        System.out.println("\nHash matching:");
        if (proteinHash == 1310777396) System.out.println("  ✓ Protein matches 1310777396");
        if (ligandHash == 1310777396) System.out.println("  ✓ Ligand matches 1310777396");
        if (complexHash == 1310777396) System.out.println("  ✓ Complex matches 1310777396");

        if (proteinHash == 1320578138) System.out.println("  ✓ Protein matches 1320578138");
        if (ligandHash == 1320578138) System.out.println("  ✓ Ligand matches 1320578138");
        if (complexHash == 1320578138) System.out.println("  ✓ Complex matches 1320578138");

        if (proteinHash == -451549042) System.out.println("  ✓ Protein matches -451549042");
        if (ligandHash == -451549042) System.out.println("  ✓ Ligand matches -451549042");
        if (complexHash == -451549042) System.out.println("  ✓ Complex matches -451549042");

        System.out.println("\n" + "=".repeat(80) + "\n");
    }

    private int computeTemplateHash(SimpleConfSpace confSpace) {
        // Use makeMolecule with empty assignment to get template molecule
        ParametricMolecule pmol = confSpace.makeMolecule(new int[confSpace.numPos()]);
        Molecule template = pmol.mol;
        int hash = 0;
        for (Residue res : template.residues) {
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
