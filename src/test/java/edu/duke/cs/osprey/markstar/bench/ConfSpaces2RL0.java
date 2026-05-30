package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;

/**
 * 2RL0 protein–ligand conf-space builders used by the BranchMARK* benchmark suite.
 * Extracted from TestBranchMARKStar to keep that file focused on tests.
 *
 * All builders use {@code examples/python.KStar/2RL0.min.reduce.pdb}:
 *   - Protein strand: chain A residues A155–A194
 *   - Ligand strand : chain G residues G648–G654 (rigid)
 */
public final class ConfSpaces2RL0 {
    private ConfSpaces2RL0() {}

    private static final String PDB_PATH = "examples/python.KStar/2RL0.min.reduce.pdb";

    private static Molecule readMol() {
        return PDBIO.read(FileTools.readFile(PDB_PATH));
    }

    private static Strand makeRigidLigand(Molecule mol, ResidueTemplateLibrary lib) {
        return new Strand.Builder(mol)
                .setTemplateLibrary(lib)
                .setResidues("G648", "G654")
                .build();
    }

    private static Strand makeProteinStrand(Molecule mol, ResidueTemplateLibrary lib) {
        return new Strand.Builder(mol)
                .setTemplateLibrary(lib)
                .setResidues("A155", "A194")
                .build();
    }

    private static TestKStar.ConfSpaces assemble(Strand protein, Strand ligand,
                                                 ForcefieldParams ffparams) {
        TestKStar.ConfSpaces cs = new TestKStar.ConfSpaces();
        cs.ffparams = ffparams;
        cs.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        cs.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        cs.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();
        return cs;
    }

    /**
     * Wild-type-only confspace with N flexible residues on protein, rigid ligand.
     * Used by smoke tests (testBranchMARKStarOnly, testFallbackOnDenseGraph) and
     * the parallelization / edge-selection ablations.
     */
    public static TestKStar.ConfSpaces buildWildTypeConfSpace(int numFlexible) {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        Strand ligand = makeRigidLigand(mol, lib);

        String[] proteinResidues = {"A156","A157","A158","A159","A160","A161","A162","A163","A164","A165",
                "A166","A167","A168","A169","A170","A171","A172","A173","A174","A175"};
        for (int i = 0; i < Math.min(numFlexible, proteinResidues.length); i++) {
            protein.flexibility.get(proteinResidues[i])
                    .setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }
        return assemble(protein, ligand, ffp);
    }

    /**
     * 4 mutable positions, high-rotamer AA palette (ARG, LYS, MET, GLU + wt) = 5^4 = 625 sequences.
     * Positions A156, A164, A172, A192. Used by GNN benchmark on highrot conf space.
     */
    public static TestKStar.ConfSpaces buildHighRotamerConfSpace() {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        protein.flexibility.get("A156").setLibraryRotamers("ARG","LYS","MET","GLU","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG","LYS","MET","GLU","TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG","LYS","MET","GLU","LEU").addWildTypeRotamers().setContinuous();

        Strand ligand = makeRigidLigand(mol, lib);
        return assemble(protein, ligand, ffp);
    }

    /**
     * 8 mutable highrot positions = 5^8 = 390,625 sequences. Primary scaling benchmark
     * conf space; matches the 8-pos GNN training set.
     */
    public static TestKStar.ConfSpaces buildHighRotamerConfSpace8pos() {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        protein.flexibility.get("A156").setLibraryRotamers("ARG","LYS","MET","GLU","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG","LYS","MET","GLU","HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG","LYS","MET","GLU","THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG","LYS","MET","GLU","ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG","LYS","MET","GLU","TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG","LYS","MET","GLU","ILE").addWildTypeRotamers().setContinuous();

        Strand ligand = makeRigidLigand(mol, lib);
        return assemble(protein, ligand, ffp);
    }

    /**
     * 8 mutable highrot (same as 8pos) + 4 flexible-WT (A157,A159,A165,A167) = 12 positions.
     * Same sequence count as 8pos benchmark at maxMut=1, but larger conformational space.
     */
    public static TestKStar.ConfSpaces buildHighRotamerConfSpace12pos() {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        protein.flexibility.get("A156").setLibraryRotamers("ARG","LYS","MET","GLU","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A158").setLibraryRotamers("ARG","LYS","MET","GLU","HIS").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A162").setLibraryRotamers("ARG","LYS","MET","GLU","THR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A168").setLibraryRotamers("ARG","LYS","MET","GLU","ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("ARG","LYS","MET","GLU","TRP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A174").setLibraryRotamers("ARG","LYS","MET","GLU","TYR").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ARG","LYS","MET","GLU","ILE").addWildTypeRotamers().setContinuous();

        for (String res : new String[]{"A157","A159","A165","A167"}) {
            protein.flexibility.get(res).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        Strand ligand = makeRigidLigand(mol, lib);
        return assemble(protein, ligand, ffp);
    }

    /**
     * 4 mutable positions, limited AA palette = 4*3*3*3 = 108 sequences.
     * Same 4 positions as buildAllMutableConfSpace; matches GNN model's 4-pos input.
     */
    public static TestKStar.ConfSpaces buildMediumMutableConfSpace() {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        protein.flexibility.get("A156").setLibraryRotamers("ALA","VAL","ILE","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers("ALA","TYR","PHE").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers("LYS","ARG","ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers("ILE","VAL","LEU").addWildTypeRotamers().setContinuous();

        Strand ligand = makeRigidLigand(mol, lib);
        return assemble(protein, ligand, ffp);
    }

    /**
     * 4 mutable positions × 20 AA = 160,000 sequences. ~300 RCs/pos → ~300^4 conformations.
     */
    public static TestKStar.ConfSpaces buildAllMutableConfSpace() {
        ForcefieldParams ffp = new ForcefieldParams();
        Molecule mol = readMol();
        ResidueTemplateLibrary lib = new ResidueTemplateLibrary.Builder(ffp.forcefld).build();

        Strand protein = makeProteinStrand(mol, lib);
        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY",
                "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
                "THR","TRP","TYR","VAL"};
        protein.flexibility.get("A156").setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A164").setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A172").setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A192").setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();

        Strand ligand = makeRigidLigand(mol, lib);
        return assemble(protein, ligand, ffp);
    }
}
