package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.lute.ConfSampler;
import edu.duke.cs.osprey.lute.LUTE;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEIO;
import edu.duke.cs.osprey.lute.LUTEPfunc;
import edu.duke.cs.osprey.lute.LUTEState;
import edu.duke.cs.osprey.lute.RandomizedDFSConfSampler;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.packstar.PackStarPartitionFunction;
import edu.duke.cs.osprey.packstar.PackStarCasePreflight;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.wmb.MeanFieldBound;
import edu.duke.cs.osprey.wmb.WeightedMiniBucket;

import java.math.BigDecimal;
import java.math.BigInteger;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.*;
import java.util.*;

/**
 * Generic PDB benchmark runner for comparing K*, MARK*, and MARK*+GNN S9.
 *
 * Reads design specification from system properties:
 *   osprey.bench.pdbPath       — path to prepped PDB
 *   osprey.bench.proteinChains — comma-separated chain IDs for protein (e.g. "A,B")
 *   osprey.bench.ligandChains  — comma-separated chain IDs for ligand (e.g. "C")
 *   osprey.bench.mutable       — semicolon-separated mutable residue IDs (e.g. "A96;B85")
 *   osprey.bench.flexible      — semicolon-separated flexible residue IDs (e.g. "C7;C6;C5")
 *   osprey.bench.method        — kstar | markstar | packstar | pac | kstar_lute | dp_profile | gnn_s9 | gnn_s10 | gnn_s11
 *   osprey.bench.epsilon       — approximation ratio (default 0.683)
 *   osprey.bench.numCPUs       — number of CPUs (default 8)
 *   osprey.bench.designId      — design identifier for output
 *   osprey.bench.outputDir     — output directory
 *
 * For GNN S9:
 *   osprey.gnn.eval.proteinModelPath
 *   osprey.gnn.eval.complexModelPath
 *   osprey.gnn.eval.proteinSubtreeModelPath
 *   osprey.gnn.eval.complexSubtreeModelPath
 *   osprey.gnn.gpuBatchSize
 */
public class GenericPDBBench {

    public static void main(String[] args) throws Exception {
        // Read design spec from system properties
        String pdbPath = System.getProperty("osprey.bench.pdbPath");
        String proteinChainsProp = System.getProperty("osprey.bench.proteinChains", "");
        String ligandChainsProp = System.getProperty("osprey.bench.ligandChains", "");
        String mutableStr = System.getProperty("osprey.bench.mutable", "");
        String flexibleStr = System.getProperty("osprey.bench.flexible", "");
        String method = System.getProperty("osprey.bench.method", "markstar");
        double epsilon = Double.parseDouble(System.getProperty("osprey.bench.epsilon", "0.683"));
        int cpus = Integer.getInteger("osprey.bench.numCPUs", 8);
        String designId = System.getProperty("osprey.bench.designId", "unknown");
        String outputDir = System.getProperty("osprey.bench.outputDir", "/usr/xtmp/lz280/bench_comparison/results");

        String[] mutableResidues = mutableStr.split(";");
        String[] flexibleResidues = flexibleStr.isEmpty() ? new String[0] : flexibleStr.split(";");

        System.out.println("==============================================");
        System.out.println("  Generic PDB Benchmark");
        System.out.println("  Design: " + designId);
        System.out.println("  PDB: " + pdbPath);
        System.out.println("  Method: " + method);
        System.out.println("  Mutable: " + mutableStr);
        System.out.println("  Flexible: " + flexibleStr);
        System.out.println("  Epsilon: " + epsilon + ", CPUs: " + cpus);
        System.out.println("==============================================");

        if (pdbPath == null || pdbPath.isEmpty()) {
            System.err.println("ERROR: osprey.bench.pdbPath is required");
            System.exit(1);
        }

        // Build conf spaces
        ForcefieldParams ffparams = new ForcefieldParams();
        Molecule mol = PDBIO.readFile(pdbPath);
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

        // All 20 AA types (matching MARK* paper: 13-19 other amino acids per mutable)
        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY",
                "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
                "THR","TRP","TYR","VAL"};

        // Detect all chain ranges
        Map<String, String[]> chainRanges = getChainRanges(mol);
        System.out.println("  Chains detected: " + chainRanges.keySet());

        Set<String> proteinChainSet = parseChainSet(proteinChainsProp);
        if (proteinChainSet.isEmpty()) {
            for (String res : mutableResidues) {
                res = res.trim();
                if (!res.isEmpty() && res.length() >= 2) proteinChainSet.add(String.valueOf(res.charAt(0)));
            }
        }
        Set<String> ligandChainSet = parseChainSet(ligandChainsProp);
        if (ligandChainSet.isEmpty()) {
            for (String ch : chainRanges.keySet()) {
                if (!proteinChainSet.contains(ch) && !ch.trim().isEmpty()) ligandChainSet.add(ch);
            }
        }
        System.out.println("  Protein chains: " + proteinChainSet);
        System.out.println("  Ligand chains: " + ligandChainSet);

        // Build one protein strand per chain so multi-chain proteins do not
        // accidentally swallow ligand/intermediate chains in PDB order.
        List<Strand> proteinStrands = new ArrayList<>();
        for (String ch : proteinChainSet) {
            String[] r = chainRanges.get(ch);
            if (r == null) {
                System.err.println("  WARNING: requested protein chain " + ch + " not detected in PDB");
                continue;
            }
            System.out.println("  Protein strand: " + r[0] + " - " + r[1]);
            proteinStrands.add(new Strand.Builder(mol).setTemplateLibrary(templateLib)
                    .setResidues(r[0], r[1]).build());
        }
        if (proteinStrands.isEmpty()) {
            throw new IllegalArgumentException("no protein strands could be built from protein chains " + proteinChainSet);
        }

        // Apply mutable residues to protein strand.
        // Token syntax: "<chain><resnum>" => all 20 AAs (back-compat);
        //               "<chain><resnum>=A,S,D" => only the listed AAs (1-letter) + wild type.
        for (String tok : mutableResidues) {
            tok = tok.trim(); if (tok.isEmpty()) continue;
            String res; String[] rotamers;
            int eq = tok.indexOf('=');
            if (eq >= 0) {
                res = tok.substring(0, eq).trim();
                rotamers = parseAAList(tok.substring(eq + 1));
            } else {
                res = tok; rotamers = all20;
            }
            if (applyResidueFlex(proteinStrands, res, rotamers)) {
                System.out.println("  Protein mutable: " + res + " -> " + String.join(",", rotamers));
            } else {
                System.err.println("  WARNING: mutable residue " + res + " not on protein strand");
            }
        }
        // Apply protein-side flexible residues
        for (String res : flexibleResidues) {
            res = res.trim(); if (res.isEmpty()) continue;
            if (!proteinChainSet.contains(String.valueOf(res.charAt(0)))) continue;
            if (applyResidueFlex(proteinStrands, res, Strand.WildType)) {
                System.out.println("  Protein flexible: " + res);
            } else {
                System.err.println("  WARNING: flexible residue " + res + " not on protein strand");
            }
        }

        // Assemble conf spaces
        TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();
        confSpaces.ffparams = ffparams;

        if (!ligandChainSet.isEmpty()) {
            // Build one ligand strand per chain so we don't accidentally include
            // protein residues that lie between ligand chains in PDB ordering.
            List<Strand> ligandStrands = new ArrayList<>();
            for (String ch : ligandChainSet) {
                String[] r = chainRanges.get(ch);
                if (r == null) continue;
                System.out.println("  Ligand strand: " + r[0] + " - " + r[1]);
                Strand ligStrand = new Strand.Builder(mol).setTemplateLibrary(templateLib)
                        .setResidues(r[0], r[1]).build();

                // Apply ligand-side flexible residues for this chain only
                for (String res : flexibleResidues) {
                    res = res.trim(); if (res.isEmpty()) continue;
                    if (!ch.equals(String.valueOf(res.charAt(0)))) continue;
                    if (applyResidueFlex(Collections.singletonList(ligStrand), res, Strand.WildType)) {
                        System.out.println("  Ligand flexible: " + res);
                    } else {
                        System.err.println("  WARNING: flexible residue " + res + " not on ligand strand");
                    }
                }
                ligandStrands.add(ligStrand);
            }

            confSpaces.protein = new SimpleConfSpace.Builder().addStrands(proteinStrands).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().addStrands(ligandStrands).build();
            SimpleConfSpace.Builder complexBuilder = new SimpleConfSpace.Builder().addStrands(proteinStrands);
            complexBuilder.addStrands(ligandStrands);
            confSpaces.complex = complexBuilder.build();
        } else {
            System.out.println("  WARNING: no ligand chains — K* will be trivial");
            confSpaces.protein = new SimpleConfSpace.Builder().addStrands(proteinStrands).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().build();
            confSpaces.complex = confSpaces.protein;
        }

        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        System.out.println("  Ligand positions: " + confSpaces.ligand.positions.size());
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());

        int wmbGpus = Integer.getInteger("osprey.wmb.numGpus", 0);
        Parallelism parallelism = wmbGpus > 0
                ? Parallelism.make(cpus, wmbGpus, Integer.getInteger("osprey.wmb.streamsPerGpu", 64))
                : Parallelism.makeCpu(cpus);
        String ematDir = outputDir + "/emat_cache/" + designId;
        new File(ematDir).mkdirs();

        long t0 = System.currentTimeMillis();

        switch (method) {
            case "kstar":
                runKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir);
                break;
            case "kstar_lute":
                runKStarLute(confSpaces, epsilon, parallelism, ematDir, designId, outputDir);
                break;
            case "markstar":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, false, false);
                break;
            case "wmb_decoupled":
                runWmbDecoupled(confSpaces, epsilon, designId, outputDir);
                break;
            case "branch":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, true, false);
                break;
            case "packstar":
            case "pac":
                runPackStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, method);
                break;
            case "dp_profile":
                runDPProfile(confSpaces, parallelism, ematDir, designId);
                break;
            case "gnn_s9":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false, false);
                break;
            case "gnn_s10":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false, false);
                break;
            case "gnn_s11":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false, false);
                break;
            case "export_gnn":
                String gnnOutputDir = System.getProperty("osprey.gnn.outputDir", outputDir + "/gnn_models/" + designId);
                int numSamples = Integer.getInteger("osprey.gnn.numSamples", 200000);
                exportGNNData(confSpaces, parallelism, ematDir, gnnOutputDir, numSamples);
                break;
            default:
                System.err.println("Unknown method: " + method);
                System.exit(1);
        }

        long elapsed = System.currentTimeMillis() - t0;
        System.out.println("\n=== TOTAL TIME: " + String.format("%.1f", elapsed / 1000.0) + " s ===");
    }

    private static void exportGNNData(TestKStar.ConfSpaces confSpaces,
                                       Parallelism parallelism, String ematDir,
                                       String gnnOutputDir, int numSamples) {
        System.out.println("\n=== Exporting GNN training data ===");
        System.out.println("  Output: " + gnnOutputDir);
        System.out.println("  Samples: " + numSamples);

        try (EnergyCalculator ecalcMin = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
             EnergyCalculator ecalcRigid = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build()) {

            for (String spaceName : new String[]{"protein", "ligand", "complex"}) {
                SimpleConfSpace cs = spaceName.equals("protein") ? confSpaces.protein
                        : spaceName.equals("ligand") ? confSpaces.ligand : confSpaces.complex;

                if (cs.positions.isEmpty()) {
                    System.out.println("  SKIP " + spaceName + ": no positions");
                    continue;
                }

                System.out.println("\n--- " + spaceName + " (" + cs.positions.size() + " positions) ---");

                ConfEnergyCalculator confEcalcMin = new ConfEnergyCalculator.Builder(cs, ecalcMin)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcMin)
                                .build().calcReferenceEnergies())
                        .build();

                EnergyMatrix ematMin = new SimplerEnergyMatrixCalculator.Builder(confEcalcMin)
                        .setCacheFile(new File(ematDir + "/export." + spaceName + ".min.dat"))
                        .build().calcEnergyMatrix();

                ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator.Builder(cs, ecalcRigid)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalcRigid)
                                .build().calcReferenceEnergies())
                        .build();

                EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                        .setCacheFile(new File(ematDir + "/export." + spaceName + ".rigid.dat"))
                        .build().calcEnergyMatrix();

                RCs rcs = new RCs(cs);
                InteractionGraph ig = InteractionGraph.buildFromEnergyMatrix(ematRigid, ematMin, rcs, 0.1);

                GNNDataExporter exporter = new GNNDataExporter(
                        confEcalcMin, ematMin, ematRigid, ig, rcs);

                File outDir = new File(gnnOutputDir + "/" + spaceName);
                exporter.export(numSamples, outDir);
                System.out.println("  Exported " + spaceName + " to " + outDir);
            }
        } catch (Exception e) {
            System.err.println("Export failed: " + e.getMessage());
            e.printStackTrace();
        }
    }

    /**
     * Parse a comma-separated list of AA codes into OSPREY residue template names.
     * Accepts standard 1-letter codes (A,S,D,...) AND already-3-letter OSPREY
     * template names (HID,HIE,HIP,...) passed through unchanged/uppercased --
     * this lets callers request alternate-protonation-state templates (which
     * exist in all_amino94.in + LovellRotamer.dat for HIS tautomers) without
     * being silently collapsed to their 1-letter parent by aa1to3().
     */
    private static String[] parseAAList(String csv) {
        String[] parts = csv.split(",");
        List<String> out = new ArrayList<>();
        for (String p : parts) {
            p = p.trim();
            if (p.isEmpty()) continue;
            out.add(p.length() > 1 ? p.toUpperCase() : aa1to3(p.charAt(0)));
        }
        return out.toArray(new String[0]);
    }

    private static String aa1to3(char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return "ALA"; case 'R': return "ARG"; case 'N': return "ASN";
            case 'D': return "ASP"; case 'C': return "CYS"; case 'E': return "GLU";
            case 'Q': return "GLN"; case 'G': return "GLY"; case 'H': return "HIS";
            case 'I': return "ILE"; case 'L': return "LEU"; case 'K': return "LYS";
            case 'M': return "MET"; case 'F': return "PHE"; case 'P': return "PRO";
            case 'S': return "SER"; case 'T': return "THR"; case 'W': return "TRP";
            case 'Y': return "TYR"; case 'V': return "VAL";
            default: throw new IllegalArgumentException("Unknown 1-letter AA code: " + c);
        }
    }

    private static void runKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                  Parallelism parallelism, String ematDir,
                                  String designId, String outputDir) {
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build()) {

            KStar.Settings settings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(1)
                    .setShowPfuncProgress(true)
                    .build();

            KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand,
                    confSpaces.complex, settings);

            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                SimpleConfSpace cs = (SimpleConfSpace) info.confSpace;
                info.confEcalc = new ConfEnergyCalculator.Builder(cs, ecalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalc)
                                .build().calcReferenceEnergies())
                        .build();
                EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
                        .setCacheFile(new File(ematDir + "/kstar." + info.type.name().toLowerCase() + ".dat"))
                        .build().calcEnergyMatrix();
                info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                        info.confEcalc,
                        new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
                        new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
                        rcs.getNumConformations());
                info.confDBFile = null;
            }

            long kstarT0 = System.currentTimeMillis();
            List<KStar.ScoredSequence> scores = kstar.run(ecalc.tasks);
            double kstarElapsed = (System.currentTimeMillis() - kstarT0) / 1000.0;
            writeKStarResults(scores, designId, "kstar", outputDir, epsilon, kstarElapsed);
        }
    }

    /**
     * K* baseline using LUTE (Hallen 2017) in place of per-conformation CCD
     * minimization: fit a pairwise-decomposable tuple expansion once per state
     * (protein/ligand/complex, after SimpleDEE pruning), then run classic K*
     * enumeration (LUTEPfunc) against the fitted energy — no PAC/deterministic
     * guarantee, a point-estimate comparison for the "learned local correction"
     * family (LUTE/EPIC) discussed in the paper's Discussion section.
     * Method name for osprey.bench.method: "kstar_lute".
     */
    private static void runKStarLute(TestKStar.ConfSpaces confSpaces, double epsilon,
                                      Parallelism parallelism, String ematDir,
                                      String designId, String outputDir) {
        try (EnergyCalculator ecalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build()) {

            KStar.Settings settings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(1)
                    .setShowPfuncProgress(true)
                    .build();

            KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand,
                    confSpaces.complex, settings);

            double luteMaxRMSE = Double.parseDouble(System.getProperty("osprey.lute.maxRMSE", "0.1"));
            double luteMaxOverfit = Double.parseDouble(System.getProperty("osprey.lute.maxOverfittingScore", "1.5"));
            int luteSeed = Integer.getInteger("osprey.lute.randomSeed", 12345);
            String luteDir = ematDir + "/lute";
            new File(luteDir).mkdirs();

            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                SimpleConfSpace cs = (SimpleConfSpace) info.confSpace;
                String stateName = info.type.name().toLowerCase();

                info.confEcalc = new ConfEnergyCalculator.Builder(cs, ecalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalc)
                                .build().calcReferenceEnergies())
                        .build();
                EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
                        .setCacheFile(new File(ematDir + "/kstar_lute." + stateName + ".dat"))
                        .build().calcEnergyMatrix();

                if (cs.positions.isEmpty()) {
                    // trivial state (e.g. rigid ligand with no flexible/mutable positions):
                    // nothing for LUTE to fit, fall back to classic pfunc (matches runKStar)
                    info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                            info.confEcalc,
                            new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
                            new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
                            rcs.getNumConformations());
                    info.confDBFile = null;
                    continue;
                }

                PruningMatrix pmat = new SimpleDEE.Runner()
                        .setSinglesThreshold(100.0)
                        .setPairsThreshold(100.0)
                        .setGoldsteinDiffThreshold(10.0)
                        .setShowProgress(true)
                        .setCacheFile(new File(luteDir + "/" + stateName + ".pmat.dat"))
                        .setParallelism(parallelism)
                        .run(cs, emat);

                File luteFile = new File(luteDir + "/" + stateName + ".dat");
                LUTEConfEnergyCalculator luteEcalc;
                if (luteFile.exists()) {
                    System.out.println("  Loading cached LUTE " + stateName + " from " + luteFile);
                    luteEcalc = new LUTEConfEnergyCalculator(cs, LUTEIO.read(luteFile));
                } else {
                    System.out.println("  Training LUTE " + stateName + " (maxRMSE=" + luteMaxRMSE + ")...");
                    File confDBFile = new File(luteDir + "/" + stateName + ".conf.db");
                    try (ConfDB confdb = new ConfDB(cs, confDBFile)) {
                        ConfDB.ConfTable confTable = confdb.new ConfTable("lute");
                        LUTE lute = new LUTE(cs);
                        ConfSampler sampler = new RandomizedDFSConfSampler(cs, pmat, luteSeed);
                        lute.sampleTuplesAndFit(info.confEcalc, emat, pmat, confTable, sampler,
                                LUTE.Fitter.OLSCG, luteMaxOverfit, luteMaxRMSE);
                        lute.reportConfSpaceSize(pmat);
                        lute.save(luteFile);
                        luteEcalc = new LUTEConfEnergyCalculator(cs, new LUTEState(lute.getTrainingSystem()));
                    }
                }

                final PruningMatrix finalPmat = pmat;
                final LUTEConfEnergyCalculator finalLuteEcalc = luteEcalc;
                info.pfuncFactory = (rcs) -> {
                    RCs prunedRcs = new RCs(rcs, finalPmat);
                    ConfAStarTree astar = new ConfAStarTree.Builder(null, prunedRcs)
                            .setLUTE(finalLuteEcalc)
                            .build();
                    return new LUTEPfunc(finalLuteEcalc, astar, prunedRcs.getNumConformations());
                };
                info.confDBFile = null;
            }

            long kstarT0 = System.currentTimeMillis();
            List<KStar.ScoredSequence> scores = kstar.run(ecalc.tasks);
            double kstarElapsed = (System.currentTimeMillis() - kstarT0) / 1000.0;
            writeKStarResults(scores, designId, "kstar_lute", outputDir, epsilon, kstarElapsed);
        }
    }

    private static void runPackStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                    Parallelism parallelism, String ematDir,
                                    String designId, String outputDir,
                                    String outputMethodName) {

        EnergyCalculator.Type ecalcType = Integer.getInteger("osprey.wmb.numGpus", 0) > 0
                ? EnergyCalculator.Type.ResidueCudaCCD : EnergyCalculator.Type.Cpu;
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType).setIsMinimizing(false).build();

        try {
            KStar.Settings settings = new KStar.Settings.Builder()
                    .setEpsilon(epsilon)
                    .setStabilityThreshold(null)
                    .setMaxSimultaneousMutations(1)
                    .setShowPfuncProgress(true)
                    .build();
            KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand,
                    confSpaces.complex, settings);

            int maxNumConfs = Integer.getInteger("osprey.packstar.maxNumConfs", -1);
            int leafMinimizationBatchSize = Integer.getInteger(
                    "osprey.packstar.leafMinimizationBatchSize", 0);
            boolean fullParallelLeafBatch = Boolean.parseBoolean(
                    System.getProperty("osprey.packstar.fullParallelLeafBatch", "false"));
            boolean reduceMinimizations = Boolean.parseBoolean(
                    System.getProperty("osprey.packstar.reduceMinimizations", "true"));
            boolean correctionTighteningEnabled = Boolean.parseBoolean(
                    System.getProperty("osprey.packstar.correctionTightening", "true"));

            for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
                SimpleConfSpace cs = (SimpleConfSpace) info.confSpace;
                String stateName = info.type.name();
                String cachePrefix = ematDir + "/packstar." + stateName.toLowerCase(Locale.ROOT);

                ConfEnergyCalculator minimizingConfEcalc = new ConfEnergyCalculator.Builder(cs, minimizingEcalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, minimizingEcalc)
                                .build().calcReferenceEnergies())
                        .build();
                ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(cs, rigidEcalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, rigidEcalc)
                                .build().calcReferenceEnergies())
                        .build();
                info.confEcalc = minimizingConfEcalc;

                EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                        .setCacheFile(new File(cachePrefix + ".rigid.dat"))
                        .build().calcEnergyMatrix();
                EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
                        .setCacheFile(new File(cachePrefix + ".minimizing.dat"))
                        .build().calcEnergyMatrix();
                UpdatingEnergyMatrix corrections = new UpdatingEnergyMatrix(cs,
                        minimizingEmat, minimizingConfEcalc);

                info.pfuncFactory = (rcs) -> {
                    PackStarPartitionFunction pfunc = new PackStarPartitionFunction(
                            cs, rigidEmat, minimizingEmat, minimizingConfEcalc,
                            rcs, parallelism, stateName);
                    pfunc.setCorrections(corrections);
                    pfunc.setReduceMinimizations(reduceMinimizations);
                    pfunc.setCorrectionTighteningEnabled(correctionTighteningEnabled);
                    if (maxNumConfs > 0) {
                        pfunc.setMaxNumConfs(maxNumConfs);
                    }
                    if (leafMinimizationBatchSize > 0) {
                        pfunc.setLeafMinimizationBatchSize(leafMinimizationBatchSize);
                    } else if (fullParallelLeafBatch) {
                        pfunc.useFullParallelLeafBatch();
                    }
                    return pfunc;
                };
                info.confDBFile = null;
            }

            // Production admission enumerates exactly the WT + mutant state set
            // that KStar.run() will request, de-duplicates filtered unbound
            // sequences, and performs only allocation-free branch-DP previews.
            // A positive caseSlaHours is the opt-in switch; an over-budget case
            // is rejected here before K* can materialize a DP table.
            BranchDpAdmission.CaseSummary admission =
                    PackStarCasePreflight.runIfConfigured(kstar);
            if (admission != null && Boolean.getBoolean(
                    "osprey.bench.packstarPreflightOnly")) {
                System.out.println(String.format(Locale.ROOT,
                        "PACK* preflight-only complete: predictedCaseHours=%.4f caseSlaHours=%.4f; formal K* run and DP-table materialization were not started.",
                        admission.totalHours(), admission.slaHours));
                return;
            }

            long packT0 = System.currentTimeMillis();
            List<KStar.ScoredSequence> scores = kstar.run(minimizingEcalc.tasks);
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
            double packElapsed = (System.currentTimeMillis() - packT0) / 1000.0;
            writeKStarResults(scores, designId, outputMethodName, outputDir, epsilon, packElapsed);
        } finally {
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
        }
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useGNN, boolean useBranch,
                                     boolean fullParallelLeafBatch) {
        runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir,
                useGNN, useBranch, fullParallelLeafBatch, 0, true);
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useGNN, boolean useBranch,
                                     boolean fullParallelLeafBatch, int leafMinimizationBatchSize) {
        runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir,
                useGNN, useBranch, fullParallelLeafBatch, leafMinimizationBatchSize, true);
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useGNN, boolean useBranch,
                                     boolean fullParallelLeafBatch, int leafMinimizationBatchSize,
                                     boolean correctionTighteningEnabled) {
        EnergyCalculator.Type ecalcType = Integer.getInteger("osprey.wmb.numGpus", 0) > 0
                ? EnergyCalculator.Type.ResidueCudaCCD : EnergyCalculator.Type.Cpu;
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType).setIsMinimizing(false).build();

        MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (cs, ecalc) ->
                new ConfEnergyCalculator.Builder(cs, ecalc)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, ecalc)
                                .build().calcReferenceEnergies())
                        .build();

        MARKStar.Settings.Builder sb = new MARKStar.Settings.Builder()
                .setEpsilon(epsilon)
                .setMaxSimultaneousMutations(1)
                .setShowPfuncProgress(true)
                .setParallelism(parallelism)
                .setEnergyMatrixCachePattern(ematDir + "/markstar.*.dat")
                .setFullParallelLeafBatch(fullParallelLeafBatch)
                .setLeafMinimizationBatchSize(leafMinimizationBatchSize)
                .setCorrectionTighteningEnabled(correctionTighteningEnabled);
        if (useBranch) sb.setUseBranchDecomposition(true);

        MARKStar markstar = new MARKStar(confSpaces.protein, confSpaces.ligand,
                confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, sb.build());
        markstar.precalcEmats();

        if (useGNN) {
            String pModel = System.getProperty("osprey.gnn.eval.proteinModelPath");
            String cModel = System.getProperty("osprey.gnn.eval.complexModelPath");
            String lModel = System.getProperty("osprey.gnn.eval.ligandModelPath");
            String pSub = System.getProperty("osprey.gnn.eval.proteinSubtreeModelPath");
            String cSub = System.getProperty("osprey.gnn.eval.complexSubtreeModelPath");
            String lSub = System.getProperty("osprey.gnn.eval.ligandSubtreeModelPath");
            int gpuBatch = Integer.getInteger("osprey.gnn.gpuBatchSize", 1000);

            // Pick strategy: s9 (default), s10, or s11.
            String gnnMethod = System.getProperty("osprey.bench.method", "gnn_s9");
            boolean useS10 = "gnn_s10".equals(gnnMethod);
            boolean useS11 = "gnn_s11".equals(gnnMethod);
            String s11LandscapeMode = System.getProperty("osprey.gnn.s11.landscapeMode", "mix")
                    .trim().toLowerCase();
            boolean s11NavigatorOn = !useS11
                    || (Boolean.parseBoolean(System.getProperty("osprey.gnn.s11.subtreeNavigator", "true"))
                    && !s11LandscapeMode.equals("off")
                    && !s11LandscapeMode.equals("none")
                    && !s11LandscapeMode.equals("false"));
            if (useS11 && !s11NavigatorOn) {
                System.out.println("  S11 subtree navigator disabled: skipping subtree model load");
            }

            // Per-confspace independent GNN setup: each (protein/ligand/complex) is enabled
            // independently if its leaf model exists; subtree is also independent.
            boolean anyGNN = false;
            String pState = "off", cState = "off", lState = "off";

            // Protein leaf GNN — enable both S9 (preferred) and S7 (fallback if no subtree)
            if (pModel != null && !pModel.isEmpty() && new File(pModel).isFile()
                    && confSpaces.protein.positions.size() > 0) {
                markstar.protein.gnnCalc = new GNNConfEnergyCalculator(
                        new File(pModel), markstar.protein.minimizingEmat,
                        confSpaces.protein.positions.size());
                if (useS11) {
                    markstar.protein.useStrategy11 = true;
                } else if (useS10) {
                    markstar.protein.useStrategy10 = true;
                } else {
                    markstar.protein.useStrategy9 = true;
                    markstar.protein.useStrategy7 = true;
                }
                markstar.protein.s7GPUBatchSize = gpuBatch;
                pState = "leaf";
                anyGNN = true;
            }
            // Complex leaf GNN
            if (cModel != null && !cModel.isEmpty() && new File(cModel).isFile()
                    && confSpaces.complex.positions.size() > 0) {
                markstar.complex.gnnCalc = new GNNConfEnergyCalculator(
                        new File(cModel), markstar.complex.minimizingEmat,
                        confSpaces.complex.positions.size());
                if (useS11) {
                    markstar.complex.useStrategy11 = true;
                } else if (useS10) {
                    markstar.complex.useStrategy10 = true;
                } else {
                    markstar.complex.useStrategy9 = true;
                    markstar.complex.useStrategy7 = true;
                }
                markstar.complex.s7GPUBatchSize = gpuBatch;
                cState = "leaf";
                anyGNN = true;
            }
            // Ligand leaf GNN
            if (lModel != null && !lModel.isEmpty() && new File(lModel).isFile()
                    && confSpaces.ligand.positions.size() > 0) {
                markstar.ligand.gnnCalc = new GNNConfEnergyCalculator(
                        new File(lModel), markstar.ligand.minimizingEmat,
                        confSpaces.ligand.positions.size());
                if (useS11) {
                    markstar.ligand.useStrategy11 = true;
                } else if (useS10) {
                    markstar.ligand.useStrategy10 = true;
                } else {
                    markstar.ligand.useStrategy9 = true;
                    markstar.ligand.useStrategy7 = true;
                }
                markstar.ligand.s7GPUBatchSize = gpuBatch;
                lState = "leaf";
                anyGNN = true;
            }
            // Subtree GNN (independent per confspace)
            if (s11NavigatorOn && pSub != null && !pSub.isEmpty() && new File(pSub).isFile()
                    && confSpaces.protein.positions.size() > 0) {
                markstar.protein.subtreeGnnCalc = new GNNSubtreeEnergyCalculator(
                        new File(pSub), markstar.protein.minimizingEmat,
                        confSpaces.protein.positions.size());
                pState = (pState.equals("leaf") ? "leaf+sub" : "sub");
                anyGNN = true;
            }
            if (s11NavigatorOn && cSub != null && !cSub.isEmpty() && new File(cSub).isFile()
                    && confSpaces.complex.positions.size() > 0) {
                markstar.complex.subtreeGnnCalc = new GNNSubtreeEnergyCalculator(
                        new File(cSub), markstar.complex.minimizingEmat,
                        confSpaces.complex.positions.size());
                cState = (cState.equals("leaf") ? "leaf+sub" : "sub");
                anyGNN = true;
            }
            if (s11NavigatorOn && lSub != null && !lSub.isEmpty() && new File(lSub).isFile()
                    && confSpaces.ligand.positions.size() > 0) {
                markstar.ligand.subtreeGnnCalc = new GNNSubtreeEnergyCalculator(
                        new File(lSub), markstar.ligand.minimizingEmat,
                        confSpaces.ligand.positions.size());
                lState = (lState.equals("leaf") ? "leaf+sub" : "sub");
                anyGNN = true;
            }

            if (anyGNN) {
                String strategyName = useS11 ? "S11" : (useS10 ? "S10" : "S9");
                System.out.println("  GNN " + strategyName + " enabled: gpuBatch=" + gpuBatch
                        + "  protein=" + pState + "  ligand=" + lState + "  complex=" + cState);
            } else {
                System.err.println("WARNING: no GNN models available, falling back to MARK* CCD");
            }
        }

        long markT0 = System.currentTimeMillis();
        List<MARKStar.ScoredSequence> scores = markstar.run();
        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();
        double markElapsed = (System.currentTimeMillis() - markT0) / 1000.0;
        String outputMethodName = System.getProperty("osprey.bench.method",
                useGNN ? "gnn_s9" : "markstar");
        writeMARKStarResults(scores, designId, outputMethodName, outputDir, epsilon, markElapsed);
    }

    private static final class WmbRun {
        String label;
        double seconds;
        double wmbSeconds;
        double meanFieldSeconds;
        double astarSeconds;
        double ccdSeconds;
        double residualSeconds;
        double initialEpsilon;
        int minimizations;
        int scored;
        String status;
        double logZLower;
        double logZUpper;
        double delta;
    }

    private static void writeWmbCsv(String outputDir, String designId, String state,
                                    int numPos, double epsilon, boolean gpu, int numGpus,
                                    int streams, int iBound, int sweeps, double ematSec,
                                    List<WmbRun> runs) {
        double totSec = 0;
        double totWmbSec = 0;
        double totMfSec = 0;
        double totAstarSec = 0;
        double totCcdSec = 0;
        double totResidualSec = 0;
        long totMin = 0, totScored = 0;
        for (WmbRun r : runs) {
            totSec += r.seconds;
            totWmbSec += r.wmbSeconds;
            totMfSec += r.meanFieldSeconds;
            totAstarSec += r.astarSeconds;
            totCcdSec += r.ccdSeconds;
            totResidualSec += r.residualSeconds;
            totMin += r.minimizations;
            totScored += r.scored;
        }
        System.out.println("\n  --- WMB+MF per-sequence totals (" + state + ", " + runs.size()
                + " sequences) ---");
        System.out.println("  total minimizations=" + totMin + "  total scored=" + totScored
                + "  total pfunc time=" + String.format("%.1f", totSec) + "s"
                + "  (+ one-time emat " + String.format("%.1f", ematSec) + "s)");
        if (totWmbSec + totMfSec + totAstarSec + totCcdSec + totResidualSec > 0) {
            System.out.println("  breakdown: wmb=" + String.format("%.3f", totWmbSec)
                    + "s  mean_field=" + String.format("%.3f", totMfSec)
                    + "s  astar=" + String.format("%.3f", totAstarSec)
                    + "s  gpu_ccd=" + String.format("%.3f", totCcdSec)
                    + "s  residual=" + String.format("%.3f", totResidualSec) + "s");
        }

        try {
            File dir = new File(outputDir, "wmb_decoupled");
            dir.mkdirs();
            File out = new File(dir, designId + "_" + state + "_perseq.csv");
            String dev = gpu ? "gpu" : "cpu";
            try (PrintWriter w = new PrintWriter(new FileWriter(out))) {
                w.println("design_id,state,num_pos,epsilon,device,num_gpus,streams_per_gpu,i_bound,"
                        + "mean_field_sweeps,emat_build_s,sequence,seconds,minimizations,scored,"
                        + "wmb_s,mean_field_s,astar_s,gpu_ccd_s,residual_s,initial_epsilon,"
                        + "log10_z_lower,log10_z_upper,eff_epsilon,status");
                for (WmbRun run : runs) {
                    w.println(String.join(",", designId, state, String.valueOf(numPos),
                            String.valueOf(epsilon), dev, String.valueOf(numGpus),
                            String.valueOf(streams), String.valueOf(iBound), String.valueOf(sweeps),
                            String.format("%.1f", ematSec), "\"" + run.label + "\"",
                            String.format("%.2f", run.seconds), String.valueOf(run.minimizations),
                            String.valueOf(run.scored), String.format("%.4f", run.wmbSeconds),
                            String.format("%.4f", run.meanFieldSeconds),
                            String.format("%.4f", run.astarSeconds),
                            String.format("%.4f", run.ccdSeconds),
                            String.format("%.4f", run.residualSeconds),
                            String.format("%.5f", run.initialEpsilon),
                            String.format("%.4f", run.logZLower),
                            String.format("%.4f", run.logZUpper), String.format("%.5f", run.delta),
                            run.status));
                }
            }
            System.out.println("  wrote " + out.getPath());
        } catch (IOException e) {
            System.err.println("  failed to write WMB csv: " + e.getMessage());
        }
    }

    /** log10 of a positive BigDecimal that may overflow double range. */
    private static double log10(BigDecimal z) {
        if (z == null || z.signum() <= 0) {
            return Double.NEGATIVE_INFINITY;
        }
        int digitsBeforePoint = z.precision() - z.scale();
        double mantissa = z.movePointLeft(digitsBeforePoint).doubleValue(); // in [0.1, 1)
        return digitsBeforePoint + Math.log10(mantissa);
    }

    /**
     * Decoupled WMB audit -- the architecture the plan actually calls for, with NO MARK* search.
     * For each sequence's pfunc: the WMB upper bound is computed once on the minimized emat and a
     * mean-field lower once on the rigid emat (both cheap, zero conformations); conformations are
     * then enumerated in emat-energy order by A* and minimized in big parallel GPU batches; the
     * un-audited residual mass is bounded deterministically by
     * min(WMB_upper - audited_emat_mass, num_remaining * exp(-threshold/RT)).  Iterates batches
     * until the bracket [max(Z_audited_exact, Z_meanfield), Z_audited_exact + residual] reaches
     * epsilon.  All the expensive CCD lands in dependency-free parallel batches that fill the GPU.
     */
    private static void runWmbDecoupled(TestKStar.ConfSpaces confSpaces, double epsilon,
                                        String designId, String outputDir) {
        int cpus = Integer.getInteger("osprey.bench.numCPUs", 8);
        int numGpus = Integer.getInteger("osprey.wmb.numGpus", 1);
        int streams = Integer.getInteger("osprey.wmb.streamsPerGpu", 64);
        int iBound = Integer.getInteger("osprey.wmb.iBound", 3);
        int sweeps = Integer.getInteger("osprey.wmb.meanFieldSweeps", 100);
        int batchSize = Math.max(1, Integer.getInteger("osprey.wmb.batchSize",
                Math.max(1, numGpus * streams) * 2));
        int initialBatchSize = Math.max(1, Integer.getInteger("osprey.wmb.initialBatchSize",
                Math.min(16, batchSize)));
        initialBatchSize = Math.min(initialBatchSize, batchSize);
        long maxAudit = Long.getLong("osprey.wmb.maxAudit", 3_000_000L);
        String stateName = System.getProperty("osprey.wmb.state", "complex").trim().toLowerCase();
        SimpleConfSpace cs = stateName.equals("protein") ? confSpaces.protein
                : stateName.equals("ligand") ? confSpaces.ligand : confSpaces.complex;
        if (!stateName.equals("protein") && !stateName.equals("ligand")) {
            stateName = "complex";
        }

        boolean gpu = numGpus > 0;
        Parallelism parallelism = gpu
                ? Parallelism.make(cpus, numGpus, streams) : Parallelism.makeCpu(cpus);
        EnergyCalculator.Type type = gpu ? EnergyCalculator.Type.ResidueCudaCCD : EnergyCalculator.Type.Cpu;
        System.out.println("\n=== WMB DECOUPLED (no search): state=" + stateName + " pos=" + cs.positions.size()
                + " device=" + (gpu ? ("gpu x" + numGpus + " streams/gpu=" + streams) : ("cpu x" + cpus))
                + " iBound=" + iBound + " batch=" + initialBatchSize + "->" + batchSize
                + " eps=" + epsilon + " ===");

        EnergyCalculator minEcalc = new EnergyCalculator.Builder(cs, confSpaces.ffparams)
                .setParallelism(parallelism).setType(type).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minEcalc)
                .setIsMinimizing(false).build();
        try {
            ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(cs, minEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, minEcalc)
                            .build().calcReferenceEnergies())
                    .build();
            ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(cs, rigidEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, rigidEcalc)
                            .build().calcReferenceEnergies())
                    .build();

            long ematT0 = System.currentTimeMillis();
            EnergyMatrix minEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                    .build().calcEnergyMatrix();
            EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                    .build().calcEnergyMatrix();
            double ematSec = (System.currentTimeMillis() - ematT0) / 1000.0;
            System.out.println("  emat build: " + String.format("%.1f", ematSec) + " s");

            int maxMut = Integer.getInteger("osprey.wmb.maxSimultaneousMutations", 1);
            List<Sequence> sequences = new ArrayList<>();
            if (cs.seqSpace.containsWildTypeSequence()) {
                sequences.add(cs.seqSpace.makeWildTypeSequence());
            }
            sequences.addAll(cs.seqSpace.getMutants(maxMut, true));
            System.out.println("  sequences: " + sequences.size());

            List<WmbRun> runs = new ArrayList<>();
            for (Sequence seq : sequences) {
                RCs seqRcs = seq.makeRCs(cs);
                String label = seq.toString(Sequence.Renderer.ResType);
                runs.add(decoupledPfunc(label, minEmat, rigidEmat, confEcalc, seqRcs,
                        iBound, sweeps, epsilon, initialBatchSize, batchSize, maxAudit));
            }
            writeWmbCsv(outputDir, designId, stateName + "_decoupled", cs.positions.size(),
                    epsilon, gpu, numGpus, streams, iBound, sweeps, ematSec, runs);
        } finally {
            rigidEcalc.clean();
            minEcalc.clean();
        }
    }

    private static WmbRun decoupledPfunc(String label, EnergyMatrix minEmat, EnergyMatrix rigidEmat,
                                         ConfEnergyCalculator confEcalc, RCs rcs, int iBound,
                                         int sweeps, double targetEps, int initialBatchSize,
                                         int maxBatchSize, long maxAudit) {
        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        double rt = bc.R * bc.T;
        ExpFunction ef = new ExpFunction(PartitionFunction.decimalPrecision);

        long t0 = System.currentTimeMillis();
        BigInteger totalConfs = rcs.getNumConformations();
        int[] noAssign = new int[rcs.getNumPos()];
        java.util.Arrays.fill(noAssign, -1);
        double wmbSec = 0;
        double mfSec = 0;
        double astarSec = 0;
        double ccdSec = 0;
        double residualSec = 0;

        // WMB upper (minimized emat) and mean-field lower (rigid emat): once, no conformations.
        BigDecimal zWmbUpper;
        long timer = System.nanoTime();
        try {
            zWmbUpper = ef.exp(WeightedMiniBucket.upperLogZ(minEmat, rcs, noAssign, iBound, rt));
        } catch (RuntimeException e) {
            zWmbUpper = null;
        } finally {
            wmbSec += elapsedSeconds(timer);
        }
        BigDecimal zMeanField = BigDecimal.ZERO;
        timer = System.nanoTime();
        try {
            double mf = MeanFieldBound.lowerLogZ(rigidEmat, rcs, noAssign,
                    sweeps, MeanFieldBound.DEFAULT_TOLERANCE, rt);
            if (!Double.isNaN(mf) && Math.abs(mf) < 1.0e4) {
                zMeanField = ef.exp(mf);
            }
        } catch (RuntimeException e) {
            // keep zero
        } finally {
            mfSec += elapsedSeconds(timer);
        }

        BigDecimal zExact = BigDecimal.ZERO;   // audited true Boltzmann mass
        BigDecimal zEmat = BigDecimal.ZERO;    // audited minimized-emat Boltzmann mass
        BigInteger audited = BigInteger.ZERO;
        BigDecimal zUpper = (zWmbUpper != null) ? zWmbUpper : null;
        BigDecimal zLower = zMeanField;
        double eps = effectiveEpsilon(zLower, zUpper);
        double initialEps = eps;
        long minimizations = 0;
        boolean exhausted = false;

        if (eps > targetEps && maxAudit > 0) {
            ConfAStarTree tree = new ConfAStarTree.Builder(minEmat, rcs).setTraditional().build();
            int batchSize = Math.max(1, Math.min(initialBatchSize, maxBatchSize));
            while (eps > targetEps && minimizations < maxAudit) {
                int thisBatchSize = (int) Math.min((long) batchSize, maxAudit - minimizations);
                List<ConfSearch.ScoredConf> batch = new ArrayList<>(thisBatchSize);
                timer = System.nanoTime();
                for (int i = 0; i < thisBatchSize; i++) {
                    ConfSearch.ScoredConf c = tree.nextConf();
                    if (c == null) { exhausted = true; break; }
                    batch.add(c);
                }
                astarSec += elapsedSeconds(timer);
                if (batch.isEmpty()) { exhausted = true; break; }

                // one dependency-free parallel GPU minimization batch
                List<Double> energies = java.util.Collections.synchronizedList(new ArrayList<>());
                timer = System.nanoTime();
                for (ConfSearch.ScoredConf conf : batch) {
                    confEcalc.tasks.submit(() -> confEcalc.calcEnergy(conf),
                            (econf) -> energies.add(econf.getEnergy()));
                }
                confEcalc.tasks.waitForFinish();
                ccdSec += elapsedSeconds(timer);
                minimizations += batch.size();

                timer = System.nanoTime();
                double thresholdEmatE = Double.NEGATIVE_INFINITY;
                for (ConfSearch.ScoredConf conf : batch) {
                    double ematE = ematEnergy(minEmat, conf.getAssignments());
                    zEmat = zEmat.add(bc.calc(ematE));
                    thresholdEmatE = Math.max(thresholdEmatE, ematE);
                }
                for (double e : energies) {
                    zExact = zExact.add(bc.calc(e));
                }
                audited = audited.add(BigInteger.valueOf(batch.size()));

                // residual upper bound on the un-audited true mass: tighter of WMB and count*exp
                BigInteger remaining = totalConfs.subtract(audited);
                BigDecimal residualCount = remaining.signum() <= 0 ? BigDecimal.ZERO
                        : new BigDecimal(remaining).multiply(bc.calc(thresholdEmatE));
                BigDecimal residual = residualCount;
                if (zWmbUpper != null) {
                    BigDecimal residualWmb = zWmbUpper.subtract(zEmat).max(BigDecimal.ZERO);
                    residual = residual.min(residualWmb);
                }

                zUpper = zExact.add(residual);
                zLower = zExact.max(zMeanField);
                eps = effectiveEpsilon(zLower, zUpper);
                residualSec += elapsedSeconds(timer);

                if (batchSize < maxBatchSize) {
                    batchSize = Math.min(maxBatchSize, Math.max(batchSize + 1, batchSize * 2));
                }
            }
        }
        double sec = (System.currentTimeMillis() - t0) / 1000.0;

        WmbRun run = new WmbRun();
        run.label = label;
        run.seconds = sec;
        run.wmbSeconds = wmbSec;
        run.meanFieldSeconds = mfSec;
        run.astarSeconds = astarSec;
        run.ccdSeconds = ccdSec;
        run.residualSeconds = residualSec;
        run.initialEpsilon = initialEps;
        run.minimizations = (int) Math.min(minimizations, Integer.MAX_VALUE);
        run.scored = run.minimizations;
        run.status = (eps <= targetEps || exhausted) ? "Estimated" : "Budget";
        run.logZLower = log10(zLower);
        run.logZUpper = (zUpper == null) ? Double.POSITIVE_INFINITY : log10(zUpper);
        run.delta = eps;
        System.out.println("  [" + label + "] time=" + String.format("%.1f", sec)
                + "s  minimizations=" + minimizations + "/" + totalConfs
                + "  logZ=[" + String.format("%.3f", run.logZLower) + ", "
                + String.format("%.3f", run.logZUpper) + "]  eps=" + String.format("%.3f", eps)
                + "  breakdown=[wmb " + String.format("%.3f", wmbSec)
                + "s, mf " + String.format("%.3f", mfSec)
                + "s, astar " + String.format("%.3f", astarSec)
                + "s, gpu_ccd " + String.format("%.3f", ccdSec)
                + "s, residual " + String.format("%.3f", residualSec) + "s]"
                + "  " + run.status);
        return run;
    }

    private static double elapsedSeconds(long startNs) {
        return (System.nanoTime() - startNs) / 1.0e9;
    }

    private static double effectiveEpsilon(BigDecimal zLower, BigDecimal zUpper) {
        if (zUpper == null) {
            return 1.0;
        }
        if (zUpper.signum() == 0) {
            return 0.0;
        }
        BigDecimal gap = zUpper.subtract(zLower).max(BigDecimal.ZERO);
        return gap.divide(zUpper, PartitionFunction.decimalPrecision).doubleValue();
    }

    /** Full minimized-emat energy of a conformation (matches WmbModel: includes the const term). */
    private static double ematEnergy(EnergyMatrix emat, int[] conf) {
        double e = emat.getConstTerm();
        for (int i = 0; i < conf.length; i++) {
            if (conf[i] < 0) {
                continue;
            }
            e += emat.getEnergy(i, conf[i]);
            for (int j = i + 1; j < conf.length; j++) {
                if (conf[j] >= 0) {
                    e += emat.getEnergy(i, conf[i], j, conf[j]);
                }
            }
        }
        return e;
    }

    private static void runDPProfile(TestKStar.ConfSpaces confSpaces,
                                      Parallelism parallelism, String ematDir,
                                      String designId) {
        if (System.getProperty("branchdp.dp.cache") == null) {
            System.setProperty("branchdp.dp.cache", "false");
        }

        String stateProp = System.getProperty("osprey.dpProfile.state", "complex")
                .trim().toLowerCase(Locale.ROOT);
        int seqIndex = Integer.getInteger("osprey.dpProfile.seqIndex", 0);
        int maxMut = Integer.getInteger("osprey.dpProfile.maxMut", 1);

        SimpleConfSpace cs;
        String stateName;
        switch (stateProp) {
            case "protein":
                cs = confSpaces.protein;
                stateName = "Protein";
                break;
            case "ligand":
                cs = confSpaces.ligand;
                stateName = "Ligand";
                break;
            case "complex":
                cs = confSpaces.complex;
                stateName = "Complex";
                break;
            default:
                throw new IllegalArgumentException("Unknown osprey.dpProfile.state: " + stateProp);
        }

        List<Sequence> sequences = new ArrayList<>();
        if (confSpaces.complex.seqSpace.containsWildTypeSequence()) {
            sequences.add(confSpaces.complex.seqSpace.makeWildTypeSequence());
        }
        sequences.addAll(confSpaces.complex.seqSpace.getMutants(maxMut, true));
        if (seqIndex < 0 || seqIndex >= sequences.size()) {
            throw new IllegalArgumentException("osprey.dpProfile.seqIndex=" + seqIndex
                    + " outside [0," + sequences.size() + ")");
        }

        Sequence globalSequence = sequences.get(seqIndex);
        Sequence stateSequence = globalSequence.filter(cs.seqSpace);
        System.out.println("\n=== DP Profile ===");
        System.out.println("  Design: " + designId);
        System.out.println("  State: " + stateName);
        System.out.println("  Sequence index: " + seqIndex + " / " + (sequences.size() - 1));
        System.out.println("  Sequence: " + globalSequence.toString(Sequence.Renderer.ResType));
        System.out.println("  State sequence: " + stateSequence.toString(Sequence.Renderer.ResType));
        System.out.println("  Positions: " + cs.positions.size());
        System.out.println("  DP cache: " + System.getProperty("branchdp.dp.cache"));

        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

        try {
            ConfEnergyCalculator minimizingConfEcalc = new ConfEnergyCalculator.Builder(cs, minimizingEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, minimizingEcalc)
                            .build().calcReferenceEnergies())
                    .build();
            ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(cs, rigidEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(cs, rigidEcalc)
                            .build().calcReferenceEnergies())
                    .build();

            String cacheFamily = System.getProperty("osprey.dpProfile.cacheFamily", "markstar")
                    .trim().toLowerCase(Locale.ROOT);
            if (!cacheFamily.equals("markstar") && !cacheFamily.equals("packstar")) {
                throw new IllegalArgumentException("Unknown osprey.dpProfile.cacheFamily: "
                        + cacheFamily);
            }
            String cachePrefix = ematDir + "/" + cacheFamily + "."
                    + stateName.toLowerCase(Locale.ROOT);
            EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                    .setCacheFile(new File(cachePrefix + ".rigid.dat"))
                    .build().calcEnergyMatrix();
            EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
                    .setCacheFile(new File(cachePrefix + ".minimizing.dat"))
                    .build().calcEnergyMatrix();

            long t0 = System.currentTimeMillis();
            new BranchMARKStarBound(cs, rigidEmat, minimizingEmat, minimizingConfEcalc,
                    stateSequence.makeRCs(cs), parallelism, stateName);
            long elapsedMs = System.currentTimeMillis() - t0;
            System.out.println("[DP_PROFILE] done design=" + designId
                    + " state=" + stateName
                    + " seqIndex=" + seqIndex
                    + " elapsedMs=" + elapsedMs
                    + " threads=" + System.getProperty("branchdp.dp.parallel.threads", "auto"));
        } finally {
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
        }
    }

    /** CSV columns shared by K* / MARK* / GNN-S9 outputs. */
    private static final String RESULTS_HEADER =
        "rank,sequence,method,target_eps,score_log10,lb_log10,ub_log10," +
        "prot_qstar_lb_log10,prot_qstar_ub_log10,prot_status,prot_eps,prot_nconf," +
        "prot_nscored,prot_npartial,prot_s9_leafGNN,prot_s9_subtreeGNN,prot_s9_ccdFromGNN,prot_s9_onnxCalls," +
        "lig_qstar_lb_log10,lig_qstar_ub_log10,lig_status,lig_eps,lig_nconf," +
        "lig_nscored,lig_npartial,lig_s9_leafGNN,lig_s9_subtreeGNN,lig_s9_ccdFromGNN,lig_s9_onnxCalls," +
        "comp_qstar_lb_log10,comp_qstar_ub_log10,comp_status,comp_eps,comp_nconf," +
        "comp_nscored,comp_npartial,comp_s9_leafGNN,comp_s9_subtreeGNN,comp_s9_ccdFromGNN,comp_s9_onnxCalls," +
        "total_time_s";

    /** Format one row from a KStarScore + sequence index. */
    private static String formatScoreRow(int rank, Sequence sequence,
                                          edu.duke.cs.osprey.kstar.KStarScore score,
                                          String method, double targetEps, double totalTimeS) {
        Double scoreLog = score.scoreLog10();  // null when not Estimated
        String scoreStr = (scoreLog == null || scoreLog.isNaN()) ? "" : String.format("%.6f", scoreLog);
        Double lbLog = score.lowerBoundLog10();
        Double ubLog = score.upperBoundLog10();
        String lbStr = (lbLog == null || lbLog.isNaN()) ? "" : String.format("%.6f", lbLog);
        String ubStr = (ubLog == null || ubLog.isNaN()) ? "" : String.format("%.6f", ubLog);
        return String.format("%d,%s,%s,%.6f,%s,%s,%s,%s,%s,%s,%.1f",
            rank,
            sequence.toString(Sequence.Renderer.ResType),
            method,
            targetEps,
            scoreStr, lbStr, ubStr,
            formatPfunc(score.protein),
            formatPfunc(score.ligand),
            formatPfunc(score.complex),
            totalTimeS);
    }

    /**
     * Returns 11 columns for one pfunc result:
     * "qstar_lb_log10,qstar_ub_log10,status,eps,numConfs,
     *  nscored,npartial,s9_leafGNN,s9_subtreeGNN,s9_ccdFromGNN,s9_onnxCalls"
     */
    private static String formatPfunc(edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Result r) {
        if (r == null) return ",,N/A,,0,0,0,0,0,0,0";
        Double lb = edu.duke.cs.osprey.kstar.KStarScore.scoreToLog10(r.values.calcLowerBound());
        Double ub = edu.duke.cs.osprey.kstar.KStarScore.scoreToLog10(r.values.calcUpperBound());
        String lbStr = (lb == null || lb.isNaN()) ? "" : String.format("%.6f", lb);
        String ubStr = (ub == null || ub.isNaN()) ? "" : String.format("%.6f", ub);
        String epsStr = "";
        try {
            double eps = r.values.getEffectiveEpsilon();
            epsStr = Double.isNaN(eps) || Double.isInfinite(eps) ? "" : String.format("%.6f", eps);
        } catch (RuntimeException e) {
            // MARK*/PackStar bounds can use MagicBigDecimal infinities; keep the
            // CSV row writable even when epsilon is not numerically meaningful.
        }
        long leafGNN  = r.getStat("s9LeafGNNBounded");
        long subGNN   = r.getStat("s9SubtreeGNNBounded");
        long ccdFromGNN = r.getStat("s9CCDFromGNN");
        long onnxCalls = r.getStat("s9LeafOnnxCalls") + r.getStat("s9SubtreeOnnxCalls");
        return String.format("%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d",
            lbStr, ubStr, r.status.name(), epsStr, r.numConfs,
            r.getStat("numConfsScored"), r.getStat("numPartialMinimizations"),
            leafGNN, subGNN, ccdFromGNN, onnxCalls);
    }

    private static void writeKStarResults(List<KStar.ScoredSequence> scores,
                                           String designId, String method, String outputDir,
                                           double targetEps, double totalTimeS) {
        String csvPath = outputDir + "/" + designId + "_" + method + ".csv";
        try (PrintWriter pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println(RESULTS_HEADER);
            scores.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
            for (int i = 0; i < scores.size(); i++) {
                KStar.ScoredSequence s = scores.get(i);
                pw.println(formatScoreRow(i + 1, s.sequence, s.score, method, targetEps, totalTimeS));
            }
            System.out.println("Results written to " + csvPath);
        } catch (IOException e) {
            System.err.println("Error writing CSV: " + e.getMessage());
        }
    }

    private static void writeMARKStarResults(List<MARKStar.ScoredSequence> scores,
                                              String designId, String method, String outputDir,
                                              double targetEps, double totalTimeS) {
        String csvPath = outputDir + "/" + designId + "_" + method + ".csv";
        try (PrintWriter pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println(RESULTS_HEADER);
            scores.sort((a, b) -> Double.compare(b.score.lowerBoundLog10(), a.score.lowerBoundLog10()));
            for (int i = 0; i < scores.size(); i++) {
                MARKStar.ScoredSequence s = scores.get(i);
                pw.println(formatScoreRow(i + 1, s.sequence, s.score, method, targetEps, totalTimeS));
            }
            System.out.println("Results written to " + csvPath);
        } catch (IOException e) {
            System.err.println("Error writing CSV: " + e.getMessage());
        }
    }

    private static Set<String> parseChainSet(String chains) {
        Set<String> out = new LinkedHashSet<>();
        for (String token : chains.split(",")) {
            token = token.trim();
            if (!token.isEmpty()) out.add(token);
        }
        return out;
    }

    private static boolean applyResidueFlex(List<Strand> strands, String res, String... rotamers) {
        for (Strand strand : strands) {
            var flex = strand.flexibility.get(res);
            if (flex == null) continue;
            flex.setLibraryRotamers(rotamers).addWildTypeRotamers().setContinuous();
            return true;
        }
        return false;
    }

    private static Map<String, String[]> getChainRanges(Molecule mol) {
        Map<String, String[]> ranges = new LinkedHashMap<>();
        Map<String, String> firstRes = new LinkedHashMap<>();
        Map<String, String> lastRes = new LinkedHashMap<>();
        for (var res : mol.residues) {
            String chain = String.valueOf(res.getChainId());
            String fullId = res.getPDBResNumber().trim();
            if (!firstRes.containsKey(chain)) firstRes.put(chain, fullId);
            lastRes.put(chain, fullId);
        }
        for (String ch : firstRes.keySet()) {
            ranges.put(ch, new String[]{firstRes.get(ch), lastRes.get(ch)});
        }
        return ranges;
    }
}
