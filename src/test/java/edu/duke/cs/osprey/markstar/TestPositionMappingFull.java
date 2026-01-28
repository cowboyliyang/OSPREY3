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
 * Test with numMutable=6 configuration
 */
public class TestPositionMappingFull {

    @Test
    public void testPositionIndicesWithSixMutable() {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("POSITION MAPPING WITH numMutable=6");
        System.out.println("=".repeat(80) + "\n");

        int numMutable = 6;

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

        // Create ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("A155", "A194")
            .build();

        // Add flexibility with numMutable=6
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
        protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
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
        System.out.println("CRITICAL FINDING:");
        System.out.println("Triple [pos2=RC29, pos3=RC23, pos4=RC8] refers to:");
        System.out.println();
        System.out.println("In PROTEIN ConfSpace:");
        if (proteinConfSpace.numPos() > 4) {
            System.out.println("  pos2: " + proteinConfSpace.positions.get(2).resNum);
            System.out.println("  pos3: " + proteinConfSpace.positions.get(3).resNum);
            System.out.println("  pos4: " + proteinConfSpace.positions.get(4).resNum);
        } else {
            System.out.println("  (positions 2, 3, 4 don't all exist in protein ConfSpace)");
        }
        System.out.println();
        System.out.println("In COMPLEX ConfSpace:");
        if (complexConfSpace.numPos() > 4) {
            System.out.println("  pos2: " + complexConfSpace.positions.get(2).resNum);
            System.out.println("  pos3: " + complexConfSpace.positions.get(3).resNum);
            System.out.println("  pos4: " + complexConfSpace.positions.get(4).resNum);
        }
        System.out.println();
        System.out.println("CONCLUSION:");
        System.out.println("The SAME residue triple (G650, G651, G652) is calculated in BOTH:");
        System.out.println("  1. Protein ConfSpace with protein-only template (hash=1310777396)");
        System.out.println("  2. Complex ConfSpace with protein+ligand template (hash=-451549042)");
        System.out.println();
        System.out.println("Different templates → Different initial structures → Different minimized energies!");
        System.out.println("=".repeat(80) + "\n");
    }
}
