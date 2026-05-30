package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.markstar.TestBranchMARKStar;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;

/**
 * EGFR GNN training-data export with a small-molecule ligand (AQ4/Erlotinib).
 *
 * Pulls GAFF parameters (prepi + frcmod + template coords + rotatable dihedrals)
 * from prepped files, builds a TestKStar.ConfSpaces with N hotspot residues
 * flexible+mutable to all 20 AA, then exports GNN data for protein + complex strands.
 */
public class EGFRBench {

    @Test
    public void exportGNNDataEGFR5pos() throws Exception {
        String proteinPdb = System.getProperty("osprey.gnn.egfr.proteinPdb",
                "gnn_data/egfr_erlotinib/protein.pdb");
        String ligandPdb = System.getProperty("osprey.gnn.egfr.ligandPdb",
                "gnn_data/egfr_erlotinib/ligand.h.pdb");
        String prepiPath = System.getProperty("osprey.gnn.egfr.prepFile",
                "gnn_data/egfr_erlotinib/ligand.prepi");
        String frcmodPath = System.getProperty("osprey.gnn.egfr.frcmodFile",
                "gnn_data/egfr_erlotinib/ligand.frcmod");
        String tcPath = System.getProperty("osprey.gnn.egfr.tcFile",
                "gnn_data/egfr_erlotinib/ligand.tc");
        String rotPath = System.getProperty("osprey.gnn.egfr.rotFile",
                "gnn_data/egfr_erlotinib/ligand.rot");
        String ligandResname = System.getProperty("osprey.gnn.egfr.ligandResname", "AQ4");
        String hotspotStr = System.getProperty("osprey.gnn.egfr.hotspotResNums", "694,766,768,772,773");
        String chain = System.getProperty("osprey.gnn.egfr.chain", "A");
        String baseDir = System.getProperty("osprey.gnn.egfr.outputDir", "gnn_data/egfr_5pos");
        int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
        int numCPUs = Integer.getInteger("osprey.gnn.numCPUs", TestBranchMARKStar.NUM_CPUs);

        int[] hotspotResNums = Arrays.stream(hotspotStr.split(","))
                .mapToInt(s -> Integer.parseInt(s.trim())).toArray();

        System.out.println("=== EGFR 5pos GNN Data Export ===");
        System.out.println("Protein PDB: " + proteinPdb);
        System.out.println("Ligand PDB: " + ligandPdb);
        System.out.println("Hotspot residues: " + hotspotStr);
        System.out.println("Ligand resname: " + ligandResname);
        System.out.println("Output: " + baseDir);

        TestKStar.ConfSpaces confSpaces = buildConfSpaceEGFR(
                new File(proteinPdb), new File(ligandPdb),
                new File(prepiPath), new File(frcmodPath),
                new File(tcPath), new File(rotPath),
                ligandResname, chain, hotspotResNums);

        System.out.println("\nProtein positions: " + confSpaces.protein.positions.size());
        for (var pos : confSpaces.protein.positions) {
            System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }
        System.out.println("Complex positions: " + confSpaces.complex.positions.size());
        for (var pos : confSpaces.complex.positions) {
            System.out.println("  " + pos.resNum + ": " + pos.resTypes.size() + " AA, " + pos.resConfs.size() + " RCs");
        }

        Parallelism parallelism = Parallelism.makeCpu(numCPUs);

        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .build()) {
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false).build();

            for (String strand : new String[]{"protein", "complex"}) {
                System.out.println("\n========== Exporting " + strand.toUpperCase()
                        + " data (EGFR 5pos) ==========");
                SimpleConfSpace cs = strand.equals("protein") ? confSpaces.protein : confSpaces.complex;
                ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, ecalcMinimized)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMinimized)
                                .build().calcReferenceEnergies())
                        .build();
                EnergyMatrix ematMin = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                        .build().calcEnergyMatrix();
                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalc, ecalcRigid);
                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .build().calcEnergyMatrix();
                RCs rcs = new RCs(cs);
                InteractionGraph ig = InteractionGraph.buildWithDualCutoff(
                        cs, ematRigid, ematMin, rcs, 8.0, 0.1);
                new GNNDataExporter(confEcalc, ematMin, ematRigid, ig, rcs)
                        .export(numSamples, new File(baseDir + "/" + strand));
            }
            System.out.println("\n========== EGFR 5pos GNN data export complete ==========");
        }
    }

    /**
     * Build TestKStar.ConfSpaces for an EGFR system with small-molecule ligand.
     * Ligand uses GAFF parameters (prepi + frcmod). Hotspot residues are flexible
     * and mutable to all 20 AA.
     */
    private TestKStar.ConfSpaces buildConfSpaceEGFR(
            File proteinPdb, File ligandPdb,
            File prepiFile, File frcmodFile, File tcFile, File rotFile,
            String ligandResname, String chain, int[] hotspotResNums
    ) throws Exception {
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

        ForcefieldParams.Forcefield ff = ForcefieldParams.Forcefield.AMBER;
        ForcefieldFileParser parser;
        if (frcmodFile.exists()) {
            parser = new ForcefieldFileParser(
                    ForcefieldParams.class.getResourceAsStream(ff.paramsPath),
                    frcmodFile.toPath());
        } else {
            parser = new ForcefieldFileParser(
                    ForcefieldParams.class.getResourceAsStream(ff.paramsPath));
        }
        confSpaces.ffparams = new ForcefieldParams(ff, parser);

        String prepiContent = Files.readString(prepiFile.toPath());
        String tcContent = tcFile.exists() ? Files.readString(tcFile.toPath()) : "";
        String rotContent = rotFile.exists() ? Files.readString(rotFile.toPath()) : "";

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addTemplates(prepiContent)
                .addTemplateCoords(tcContent)
                .addRotamers(rotContent)
                .build();

        Molecule proteinMol = PDBIO.readFile(proteinPdb);
        Molecule ligandMol = PDBIO.readFile(ligandPdb);

        int firstRes = Integer.MAX_VALUE, lastRes = Integer.MIN_VALUE;
        for (Residue res : proteinMol.residues) {
            String resNum = Residues.normalizeResNum(res.getPDBResNumber());
            int i = 0;
            while (i < resNum.length() && !Character.isDigit(resNum.charAt(i))) i++;
            int rn = Integer.parseInt(resNum.substring(i).replaceAll("[^0-9-].*$", ""));
            if (rn < firstRes) firstRes = rn;
            if (rn > lastRes) lastRes = rn;
        }
        System.out.println("Protein residues: " + chain + firstRes + " to " + chain + lastRes);

        Strand protein = new Strand.Builder(proteinMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(chain + firstRes, chain + lastRes)
                .build();

        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
                "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
        for (int resNum : hotspotResNums) {
            String id = chain + resNum;
            Strand.ResidueFlex rf = protein.flexibility.get(id);
            if (rf != null) {
                rf.setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
                System.out.println("  Flexible: " + id + " (wt=" + rf.wildType + ")");
            } else {
                System.out.println("  WARNING: " + id + " not found in protein flexibility");
            }
        }

        int ligFirstRes = Integer.MAX_VALUE, ligLastRes = Integer.MIN_VALUE;
        String ligChain = chain;
        for (Residue res : ligandMol.residues) {
            try {
                String rnStr = Residues.normalizeResNum(res.getPDBResNumber());
                int i = 0;
                while (i < rnStr.length() && !Character.isDigit(rnStr.charAt(i))) i++;
                int rn = Integer.parseInt(rnStr.substring(i).replaceAll("[^0-9-].*$", ""));
                if (rn < ligFirstRes) { ligFirstRes = rn; ligChain = String.valueOf(res.getChainId()); }
                if (rn > ligLastRes) ligLastRes = rn;
            } catch (NumberFormatException e) {
                // skip non-numeric residue IDs (e.g., waters)
            }
        }
        Strand ligand = new Strand.Builder(ligandMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(ligChain + ligFirstRes, ligChain + ligLastRes)
                .build();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();
        return confSpaces;
    }
}
