package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

/**
 * Test to understand position mapping in protein/ligand/complex ConfSpaces
 */
public class TestPositionMapping {

    @Test
    public void testPositionIndices() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("POSITION MAPPING IN PROTEIN/LIGAND/COMPLEX CONFSPACES");
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

        // Print position information
        System.out.println("PROTEIN ConfSpace:");
        System.out.println("  Total positions: " + proteinConfSpace.numPos());
        for (int i = 0; i < proteinConfSpace.numPos(); i++) {
            SimpleConfSpace.Position pos = proteinConfSpace.positions.get(i);
            System.out.println("    pos" + i + ": " + pos.resNum + " (" + pos.resFlex.wildType + "), RCs=" + pos.resConfs.size());
        }
        System.out.println();

        System.out.println("LIGAND ConfSpace:");
        System.out.println("  Total positions: " + ligandConfSpace.numPos());
        for (int i = 0; i < ligandConfSpace.numPos(); i++) {
            SimpleConfSpace.Position pos = ligandConfSpace.positions.get(i);
            System.out.println("    pos" + i + ": " + pos.resNum + " (" + pos.resFlex.wildType + "), RCs=" + pos.resConfs.size());
        }
        System.out.println();

        System.out.println("COMPLEX ConfSpace:");
        System.out.println("  Total positions: " + complexConfSpace.numPos());
        for (int i = 0; i < complexConfSpace.numPos(); i++) {
            SimpleConfSpace.Position pos = complexConfSpace.positions.get(i);
            System.out.println("    pos" + i + ": " + pos.resNum + " (" + pos.resFlex.wildType + "), RCs=" + pos.resConfs.size());
        }
        System.out.println();

        System.out.println("=".repeat(80));
        System.out.println("ANALYSIS:");
        System.out.println("If we see a triple like [pos2=RC29, pos3=RC23, pos4=RC8], we need to know:");
        System.out.println("1. Which ConfSpace was this triple calculated in?");
        System.out.println("2. What residues do pos2, pos3, pos4 correspond to in that ConfSpace?");
        System.out.println("3. Why would the same residue triple be calculated in different ConfSpaces?");
        System.out.println("=".repeat(80) + "\n");
    }
}
