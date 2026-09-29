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
 * Generic PDB benchmark runner for comparing K*, MARK*, and PACK*.
 *
 * Reads design specification from system properties:
 *   osprey.bench.pdbPath       — path to prepped PDB
 *   osprey.bench.proteinChains — comma-separated chain IDs for protein (e.g. "A,B")
 *   osprey.bench.ligandChains  — comma-separated chain IDs for ligand (e.g. "C")
 *   osprey.bench.mutable       — semicolon-separated mutable residue IDs (e.g. "A96;B85")
 *   osprey.bench.flexible      — semicolon-separated flexible residue IDs (e.g. "C7;C6;C5")
 *   osprey.bench.method        — kstar | markstar | packstar | pac | kstar_lute | dp_profile
 *   osprey.bench.epsilon       — approximation ratio (default 0.683)
 *   osprey.bench.numCPUs       — number of CPUs (default 8)
 *   osprey.bench.designId      — design identifier for output
 *   osprey.bench.outputDir     — output directory
 *
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

        int numGpus = Integer.getInteger("osprey.bench.numGpus", 0);
        Parallelism parallelism = numGpus > 0
                ? Parallelism.make(cpus, numGpus, Integer.getInteger("osprey.bench.streamsPerGpu", 64))
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
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, false);
                break;
            case "branch":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false);
                break;
            case "packstar":
            case "pac":
                runPackStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, method);
                break;
            case "packstar_pfunc":
                runPackStarPfunc(confSpaces, epsilon, parallelism, ematDir, designId, outputDir);
                break;
            case "sequence_dump":
                runSequenceDump(confSpaces, outputDir, designId);
                break;
            case "dp_profile":
                runDPProfile(confSpaces, parallelism, ematDir, designId);
                break;
            default:
                System.err.println("Unknown method: " + method);
                System.exit(1);
        }

        long elapsed = System.currentTimeMillis() - t0;
        System.out.println("\n=== TOTAL TIME: " + String.format("%.1f", elapsed / 1000.0) + " s ===");
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

        EnergyCalculator.Type ecalcType = Integer.getInteger("osprey.bench.numGpus", 0) > 0
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
                if ("auto".equalsIgnoreCase(edu.duke.cs.osprey.packstar.PackStarConfig
                        .getProperty("packstar.admission.mode", "sla").trim())) {
                    System.out.println("PACK* automatic preflight complete: structural optimization only; "
                            + "formal K* run and DP-table materialization were not started.");
                } else System.out.println(String.format(Locale.ROOT,
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

    /**
     * Run exactly one PACK* partition function for one state and one global
     * sequence. This lets rescue-array tasks retry only the pfunc that failed
     * in a previous full K* run.
     */
    private static void runPackStarPfunc(TestKStar.ConfSpaces confSpaces,
                                         double epsilon,
                                         Parallelism parallelism,
                                         String ematDir,
                                         String designId,
                                         String outputDir) {
        String stateProp = System.getProperty(
                "osprey.packstarPfunc.state", "complex")
                .trim().toLowerCase(Locale.ROOT);
        int seqIndex = Integer.getInteger("osprey.packstarPfunc.seqIndex", 0);
        int maxMut = Integer.getInteger("osprey.packstarPfunc.maxMut", 1);

        SimpleConfSpace cs;
        String stateName;
        int stateInstanceId;
        switch (stateProp) {
            case "protein":
                cs = confSpaces.protein;
                stateName = "Protein";
                stateInstanceId = 0;
                break;
            case "ligand":
                cs = confSpaces.ligand;
                stateName = "Ligand";
                stateInstanceId = 1;
                break;
            case "complex":
                cs = confSpaces.complex;
                stateName = "Complex";
                stateInstanceId = 2;
                break;
            default:
                throw new IllegalArgumentException(
                        "Unknown osprey.packstarPfunc.state: " + stateProp);
        }

        List<Sequence> sequences = new ArrayList<>();
        if (confSpaces.complex.seqSpace.containsWildTypeSequence()) {
            sequences.add(confSpaces.complex.seqSpace.makeWildTypeSequence());
        }
        sequences.addAll(confSpaces.complex.seqSpace.getMutants(maxMut, true));
        if (seqIndex < 0 || seqIndex >= sequences.size()) {
            throw new IllegalArgumentException("osprey.packstarPfunc.seqIndex="
                    + seqIndex + " outside [0," + sequences.size() + ")");
        }
        Sequence globalSequence = sequences.get(seqIndex);
        Sequence stateSequence = globalSequence.filter(cs.seqSpace);

        System.out.println("\n=== PACK* single-pfunc run ===");
        System.out.println("  Design: " + designId);
        System.out.println("  State: " + stateName);
        System.out.println("  Sequence index: " + seqIndex + " / "
                + (sequences.size() - 1));
        System.out.println("  Sequence: "
                + globalSequence.toString(Sequence.Renderer.ResType));

        EnergyCalculator.Type ecalcType = Integer.getInteger("osprey.bench.numGpus", 0) > 0
                ? EnergyCalculator.Type.ResidueCudaCCD : EnergyCalculator.Type.Cpu;
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setType(ecalcType)
                .setIsMinimizing(false).build();

        try {
            ConfEnergyCalculator minimizingConfEcalc =
                    new ConfEnergyCalculator.Builder(cs, minimizingEcalc)
                    .setReferenceEnergies(
                            new SimplerEnergyMatrixCalculator.Builder(
                                    cs, minimizingEcalc)
                                    .build().calcReferenceEnergies())
                    .build();
            ConfEnergyCalculator rigidConfEcalc =
                    new ConfEnergyCalculator.Builder(cs, rigidEcalc)
                    .setReferenceEnergies(
                            new SimplerEnergyMatrixCalculator.Builder(
                                    cs, rigidEcalc)
                                    .build().calcReferenceEnergies())
                    .build();

            String cachePrefix = ematDir + "/packstar."
                    + stateName.toLowerCase(Locale.ROOT);
            EnergyMatrix rigidEmat =
                    new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                    .setCacheFile(new File(cachePrefix + ".rigid.dat"))
                    .build().calcEnergyMatrix();
            EnergyMatrix minimizingEmat =
                    new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
                    .setCacheFile(new File(cachePrefix + ".minimizing.dat"))
                    .build().calcEnergyMatrix();
            UpdatingEnergyMatrix corrections = new UpdatingEnergyMatrix(
                    cs, minimizingEmat, minimizingConfEcalc);

            PartitionFunction.Result result;
            long start = System.currentTimeMillis();
            try (PackStarPartitionFunction pfunc =
                         new PackStarPartitionFunction(
                                 cs, rigidEmat, minimizingEmat,
                                 minimizingConfEcalc,
                                 stateSequence.makeRCs(cs), parallelism,
                                 stateName)) {
                pfunc.setCorrections(corrections);
                pfunc.setReduceMinimizations(Boolean.parseBoolean(
                        System.getProperty(
                                "osprey.packstar.reduceMinimizations", "true")));
                pfunc.setCorrectionTighteningEnabled(Boolean.parseBoolean(
                        System.getProperty(
                                "osprey.packstar.correctionTightening", "true")));
                pfunc.setInstanceId(stateInstanceId);
                pfunc.setReportProgress(true);
                pfunc.init(epsilon);
                pfunc.compute();
                result = pfunc.makeResult();
            }
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
            double elapsedS = (System.currentTimeMillis() - start) / 1000.0;
            double effectiveEps;
            try {
                effectiveEps = result.values.getEffectiveEpsilon();
            } catch (RuntimeException ex) {
                if (result.status == PartitionFunction.Status.Estimated) {
                    throw ex;
                }
                effectiveEps = 1.0;
            }
            Double lowerLog10 = edu.duke.cs.osprey.kstar.KStarScore.scoreToLog10(
                    result.values.calcLowerBound());
            Double upperLog10 = edu.duke.cs.osprey.kstar.KStarScore.scoreToLog10(
                    result.values.calcUpperBound());
            String lowerText = lowerLog10 == null || !Double.isFinite(lowerLog10)
                    ? "" : String.format(Locale.ROOT, "%.9f", lowerLog10);
            String upperText = upperLog10 == null || !Double.isFinite(upperLog10)
                    ? "" : String.format(Locale.ROOT, "%.9f", upperLog10);

            String summaryPath = System.getProperty(
                    "osprey.packstarPfunc.outputTsv",
                    outputDir + "/" + designId + "_packstar_pfunc.tsv");
            File summaryFile = new File(summaryPath);
            if (summaryFile.getParentFile() != null) {
                summaryFile.getParentFile().mkdirs();
            }
            boolean writeHeader = !summaryFile.exists();
            try (PrintWriter writer = new PrintWriter(
                    new FileWriter(summaryFile, true))) {
                if (writeHeader) {
                    writer.println("design_id\tstate\tseq_index\tsequence\tstatus\tepsilon\tlower_log10\tupper_log10\tnum_confs\telapsed_s");
                }
                writer.printf(Locale.ROOT,
                        "%s\t%s\t%d\t%s\t%s\t%.9f\t%s\t%s\t%d\t%.3f%n",
                        designId, stateName, seqIndex,
                        globalSequence.toString(Sequence.Renderer.ResType),
                        result.status.name(), effectiveEps,
                        lowerText, upperText, result.numConfs, elapsedS);
            } catch (IOException ex) {
                throw new RuntimeException("PACK* pfunc summary write failed", ex);
            }
            System.out.println(String.format(Locale.ROOT,
                    "[PACKSTAR_PFUNC_DONE] design=%s state=%s seqIndex=%d status=%s epsilon=%.9f elapsed=%.3fs output=%s",
                    designId, stateName, seqIndex, result.status.name(),
                    effectiveEps, elapsedS, summaryPath));
        } finally {
            minimizingEcalc.tasks.waitForFinish();
            rigidEcalc.tasks.waitForFinish();
        }
    }

    /** Dump the exact global WT + single-mutant order used by KStar. */
    private static void runSequenceDump(TestKStar.ConfSpaces confSpaces,
                                        String outputDir,
                                        String designId) {
        int maxMut = Integer.getInteger("osprey.sequenceDump.maxMut", 1);
        List<Sequence> sequences = new ArrayList<>();
        if (confSpaces.complex.seqSpace.containsWildTypeSequence()) {
            sequences.add(confSpaces.complex.seqSpace.makeWildTypeSequence());
        }
        sequences.addAll(confSpaces.complex.seqSpace.getMutants(maxMut, true));
        String outputPath = System.getProperty("osprey.sequenceDump.output",
                outputDir + "/" + designId + "_sequences.tsv");
        File output = new File(outputPath);
        if (output.getParentFile() != null) {
            output.getParentFile().mkdirs();
        }
        try (PrintWriter writer = new PrintWriter(new FileWriter(output))) {
            writer.println("design_id\tseq_index\tsequence");
            for (int i = 0; i < sequences.size(); i++) {
                writer.printf(Locale.ROOT, "%s\t%d\t%s%n", designId, i,
                        sequences.get(i).toString(Sequence.Renderer.ResType));
            }
        } catch (IOException ex) {
            throw new RuntimeException("sequence dump write failed", ex);
        }
        System.out.println("[SEQUENCE_DUMP_DONE] design=" + designId
                + " count=" + sequences.size() + " output=" + outputPath);
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useBranch,
                                     boolean fullParallelLeafBatch) {
        runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir,
                useBranch, fullParallelLeafBatch, 0, true);
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useBranch,
                                     boolean fullParallelLeafBatch, int leafMinimizationBatchSize) {
        runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir,
                useBranch, fullParallelLeafBatch, leafMinimizationBatchSize, true);
    }

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useBranch,
                                     boolean fullParallelLeafBatch, int leafMinimizationBatchSize,
                                     boolean correctionTighteningEnabled) {
        EnergyCalculator.Type ecalcType = Integer.getInteger("osprey.bench.numGpus", 0) > 0
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

        long markT0 = System.currentTimeMillis();
        List<MARKStar.ScoredSequence> scores = markstar.run();
        minimizingEcalc.tasks.waitForFinish();
        rigidEcalc.tasks.waitForFinish();
        double markElapsed = (System.currentTimeMillis() - markT0) / 1000.0;
        String outputMethodName = System.getProperty("osprey.bench.method",
                "markstar");
        writeMARKStarResults(scores, designId, outputMethodName, outputDir, epsilon, markElapsed);
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

    /** CSV schema shared by K*, MARK*, and PACK*. */
    private static final String RESULTS_HEADER =
        "rank,sequence,method,target_eps,score_log10,lb_log10,ub_log10," +
        "prot_qstar_lb_log10,prot_qstar_ub_log10,prot_status,prot_eps,prot_nconf," +
        "prot_nscored,prot_npartial," +
        "lig_qstar_lb_log10,lig_qstar_ub_log10,lig_status,lig_eps,lig_nconf," +
        "lig_nscored,lig_npartial," +
        "comp_qstar_lb_log10,comp_qstar_ub_log10,comp_status,comp_eps,comp_nconf," +
        "comp_nscored,comp_npartial," +
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
     * Returns seven columns for one pfunc result:
     * "qstar_lb_log10,qstar_ub_log10,status,eps,numConfs,
     *  nscored,npartial"
     */
    private static String formatPfunc(edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Result r) {
        if (r == null) return ",,N/A,,0,0,0";
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
        return String.format("%s,%s,%s,%s,%d,%d,%d",
            lbStr, ubStr, r.status.name(), epsStr, r.numConfs,
            r.getStat("numConfsScored"), r.getStat("numPartialMinimizations"));
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
