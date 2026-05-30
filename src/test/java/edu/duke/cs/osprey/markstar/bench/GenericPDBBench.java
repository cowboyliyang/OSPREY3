package edu.duke.cs.osprey.markstar.bench;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNDataExporter;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
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
 *   osprey.bench.method        — kstar | markstar | gnn_s9 | gnn_s10 | gnn_s11
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
        String[] proteinChains = System.getProperty("osprey.bench.proteinChains", "A").split(",");
        String[] ligandChains = System.getProperty("osprey.bench.ligandChains", "B").split(",");
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

        // Protein chains = chains with MUTABLE residues
        // Ligand chains = all other non-empty chains
        Set<String> proteinChainSet = new LinkedHashSet<>();
        for (String res : mutableResidues) {
            res = res.trim();
            if (!res.isEmpty() && res.length() >= 2) proteinChainSet.add(String.valueOf(res.charAt(0)));
        }
        Set<String> ligandChainSet = new LinkedHashSet<>();
        for (String ch : chainRanges.keySet()) {
            if (!proteinChainSet.contains(ch) && !ch.trim().isEmpty()) ligandChainSet.add(ch);
        }
        System.out.println("  Protein chains: " + proteinChainSet);
        System.out.println("  Ligand chains: " + ligandChainSet);

        // Build protein strand
        String pFirst = null, pLast = null;
        for (String ch : proteinChainSet) {
            String[] r = chainRanges.get(ch);
            if (r != null) { if (pFirst == null) pFirst = r[0]; pLast = r[1]; }
        }
        System.out.println("  Protein strand: " + pFirst + " - " + pLast);

        Strand proteinStrand = new Strand.Builder(mol).setTemplateLibrary(templateLib)
                .setResidues(pFirst, pLast).build();

        // Apply mutable residues to protein strand
        for (String res : mutableResidues) {
            res = res.trim(); if (res.isEmpty()) continue;
            try {
                proteinStrand.flexibility.get(res).setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
                System.out.println("  Protein mutable: " + res);
            } catch (Exception e) {
                System.err.println("  WARNING: mutable residue " + res + " not on protein strand: " + e.getMessage());
            }
        }
        // Apply protein-side flexible residues
        for (String res : flexibleResidues) {
            res = res.trim(); if (res.isEmpty()) continue;
            if (!proteinChainSet.contains(String.valueOf(res.charAt(0)))) continue;
            try {
                proteinStrand.flexibility.get(res).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
                System.out.println("  Protein flexible: " + res);
            } catch (Exception e) {
                System.err.println("  WARNING: flexible residue " + res + " not on protein strand: " + e.getMessage());
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
                    try {
                        ligStrand.flexibility.get(res).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
                        System.out.println("  Ligand flexible: " + res);
                    } catch (Exception e) {
                        System.err.println("  WARNING: flexible residue " + res + " not on ligand strand: " + e.getMessage());
                    }
                }
                ligandStrands.add(ligStrand);
            }

            confSpaces.protein = new SimpleConfSpace.Builder().addStrand(proteinStrand).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().addStrands(ligandStrands).build();
            SimpleConfSpace.Builder complexBuilder = new SimpleConfSpace.Builder().addStrand(proteinStrand);
            for (Strand ls : ligandStrands) complexBuilder.addStrand(ls);
            confSpaces.complex = complexBuilder.build();
        } else {
            System.out.println("  WARNING: no ligand chains — K* will be trivial");
            confSpaces.protein = new SimpleConfSpace.Builder().addStrand(proteinStrand).build();
            confSpaces.ligand = new SimpleConfSpace.Builder().build();
            confSpaces.complex = confSpaces.protein;
        }

        System.out.println("  Protein positions: " + confSpaces.protein.positions.size());
        System.out.println("  Ligand positions: " + confSpaces.ligand.positions.size());
        System.out.println("  Complex positions: " + confSpaces.complex.positions.size());

        Parallelism parallelism = Parallelism.makeCpu(cpus);
        String ematDir = outputDir + "/emat_cache/" + designId;
        new File(ematDir).mkdirs();

        long t0 = System.currentTimeMillis();

        switch (method) {
            case "kstar":
                runKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir);
                break;
            case "markstar":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, false);
                break;
            case "branch":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, true);
                break;
            case "pac":
                // PAC (Probably Approximately Correct) estimation runs on top of the
                // BranchMARK* tree. Force the PAC flag on so BranchMARKStarBound bypasses
                // the exact search loop and uses Rao-Blackwellized importance sampling.
                System.setProperty("branchmarkstar.usePAC", "true");
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, false, true);
                break;
            case "gnn_s9":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false);
                break;
            case "gnn_s10":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false);
                break;
            case "gnn_s11":
                runMARKStar(confSpaces, epsilon, parallelism, ematDir, designId, outputDir, true, false);
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

    private static void runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon,
                                     Parallelism parallelism, String ematDir,
                                     String designId, String outputDir,
                                     boolean useGNN, boolean useBranch) {
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(
                confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism).setIsMinimizing(false).build();

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
                .setEnergyMatrixCachePattern(ematDir + "/markstar.*.dat");
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
        double eps = r.values.getEffectiveEpsilon();
        String epsStr = Double.isNaN(eps) || Double.isInfinite(eps) ? "" : String.format("%.6f", eps);
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
