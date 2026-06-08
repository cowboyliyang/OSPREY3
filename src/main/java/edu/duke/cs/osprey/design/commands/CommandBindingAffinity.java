package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.StrandFlex;
import edu.duke.cs.osprey.design.*;
import edu.duke.cs.osprey.design.models.AffinityDesign;
import edu.duke.cs.osprey.design.models.MoleculeDto;
import edu.duke.cs.osprey.design.models.ResidueModifier;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import edu.duke.cs.osprey.energy.approximation.branch.GNNConfEnergyCalculator;
import edu.duke.cs.osprey.energy.approximation.branch.GNNSubtreeEnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.markstar.framework.BranchMARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import one.util.streamex.IntStreamEx;
import one.util.streamex.StreamEx;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Parameters(commandDescription = CommandBindingAffinity.CommandDescription)
public class CommandBindingAffinity extends RunnableCommand {

    public static final String CommandName = "affinity";
    public static final String CommandDescription = "Compute an epsilon approximation to binding affinity (K*).";

    @Parameter(names = "--do-scan", description = "Runs a scan using the scan settings specified in the design.")
    public boolean doScan;

    @Parameter(names = "--scan-flex-distance", description = "Distance (in angstroms) around mutable residues in which residues are flexible.")
    public double scanFlexDistance = 2.0;

    @Parameter(names = "--scan-distance", description = "Distance (in angstroms) around a scan's target residue in which residues are made mutable.")
    public double scanDistance = 5.0;

    @Parameter(names = "--scan-output", description = "Specifies the output directory to save the scan designs in.")
    public String scanOutput;

    @Parameter(names = "--max-num-confs", description = "Sets an upper bound on the number of conformations evaluated.")
    private int maxNumberConfs = -1;

    @Parameter(names = "--use-markstar", description = "Use MARK* instead of Gradient Descent for the partition function calculation.")
    public boolean useMarkstar;

    @Parameter(names = "--use-branchmarkstar", description = "Use BranchMARK* (branch-decomposition variant of MARK*) for the partition function calculation. Faster than MARK* on kinase-sized active sites.")
    public boolean useBranchMarkstar;

    @Parameter(names = "--use-branchmarkstar-pac", description = "Use BranchMARK* with PAC partition-function sampling.")
    public boolean useBranchMarkstarPac;

    @Parameter(names = "--branchmarkstar-pac-samples", description = "Number of PAC samples for BranchMARK* PAC mode. Uses the BranchMARK* default when omitted.")
    public int branchMarkstarPacSamples = -1;

    @Parameter(names = "--branchmarkstar-pac-confidence", description = "PAC confidence level for BranchMARK* PAC mode. Uses the BranchMARK* default when omitted.")
    public double branchMarkstarPacConfidence = Double.NaN;

    @Parameter(names = "--branchmarkstar-pac-residual-bound", description = "Deterministic absolute residual-energy bound |xi| in kcal/mol. Required for a strict PAC certificate.")
    public double branchMarkstarPacResidualBound = Double.NaN;

    // === GNN Strategy 8 (leaf + subtree GNN) — only used with --use-branchmarkstar ===
    @Parameter(names = "--gnn-strategy", description = "GNN strategy ID (0=off, 8=leaf+subtree). Only effective with --use-branchmarkstar.")
    public int gnnStrategy = 0;

    @Parameter(names = "--gnn-protein-leaf", description = "Path to protein leaf-GNN ONNX model.")
    public String gnnProteinLeaf;

    @Parameter(names = "--gnn-protein-subtree", description = "Path to protein subtree-GNN ONNX model.")
    public String gnnProteinSubtree;

    @Parameter(names = "--gnn-complex-leaf", description = "Path to complex leaf-GNN ONNX model.")
    public String gnnComplexLeaf;

    @Parameter(names = "--gnn-complex-subtree", description = "Path to complex subtree-GNN ONNX model.")
    public String gnnComplexSubtree;

    // Training conf-space spec (for RC mapping)
    @Parameter(names = "--gnn-train-protein-pdb", description = "Training protein PDB.")
    public String gnnTrainProteinPdb;

    @Parameter(names = "--gnn-train-ligand-pdb", description = "Training ligand PDB (HETATM).")
    public String gnnTrainLigandPdb;

    @Parameter(names = "--gnn-train-prepi", description = "Training ligand prepi file.")
    public String gnnTrainPrepi;

    @Parameter(names = "--gnn-train-frcmod", description = "Training ligand frcmod file.")
    public String gnnTrainFrcmod;

    @Parameter(names = "--gnn-train-tc", description = "Training ligand tc file.")
    public String gnnTrainTc;

    @Parameter(names = "--gnn-train-rot", description = "Training ligand rot file.")
    public String gnnTrainRot;

    @Parameter(names = "--gnn-train-chain", description = "Training protein chain ID (default: A).")
    public String gnnTrainChain = "A";

    @Parameter(names = "--gnn-train-hotspots", description = "Comma-separated training hotspot residue numbers, e.g. 694,766,768,772,773.")
    public String gnnTrainHotspots;

    @Parameter(names = "--gnn-gpu-batch", description = "GPU batch size for leaf GNN (default 1000).")
    public int gnnGpuBatch = 1000;

    @Parameter(names = "--gnn-leaf-gpu-protein", description = "GPU device id for protein leaf GNN (default 0).")
    public int gnnLeafGpuProtein = 0;

    @Parameter(names = "--gnn-subtree-gpu-protein", description = "GPU device id for protein subtree GNN (default 0).")
    public int gnnSubtreeGpuProtein = 0;

    @Parameter(names = "--gnn-leaf-gpu-complex", description = "GPU device id for complex leaf GNN (default 0).")
    public int gnnLeafGpuComplex = 0;

    @Parameter(names = "--gnn-subtree-gpu-complex", description = "GPU device id for complex subtree GNN (default 0).")
    public int gnnSubtreeGpuComplex = 0;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    private String[] args;

    private Map<String, String> dbSettings = new HashMap<>();
    private static final String dbHostnameKey = "dbHostname";
    private static final String dbPortKey = "dbPort";
    private static final String dbNameKey = "dbName";
    private static final String dbUserNameKey = "dbUserName";
    private static final String dbPasswordKey = "dbPassword";
    private static final List<String> allMutableAaTypes = List.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIE", "HID", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL");

    @Override
    public int run(JCommander commander, String[] args) {
        this.args = args;

        return processHelpAndNoArgs(commander, args)
                .orElseGet(() -> parseAndValidate(delegate.design)
                .map(this::runAffinityDesign)
                .orElse(Main.Failure));
    }

    private int runAffinityDesign(AffinityDesign design) {
        if (useBranchMarkstarPac) {
            useBranchMarkstar = true;
            System.setProperty("branchmarkstar.usePAC", "true");
            if (branchMarkstarPacSamples > 0) {
                System.setProperty("branchmarkstar.pac.samples", Integer.toString(branchMarkstarPacSamples));
            }
            if (!Double.isNaN(branchMarkstarPacConfidence) && branchMarkstarPacConfidence > 0) {
                System.setProperty("branchmarkstar.pac.confidence", Double.toString(branchMarkstarPacConfidence));
            }
            if (!Double.isNaN(branchMarkstarPacResidualBound) && branchMarkstarPacResidualBound >= 0) {
                System.setProperty("branchmarkstar.pac.residualBound",
                        Double.toString(branchMarkstarPacResidualBound));
            }
        }

        var paramsAndStrands = new ForceFieldParamsAndStrands(delegate, design);

        if (doScan && design.scanSettings != null) {
            return makeScanDesigns(design, paramsAndStrands.makeAllResidues());
        }

        // Exit early if just trying to validate input
        if (delegate.verifyInput) {
            System.out.println("Design file validated.");
            return Main.Success;
        }

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        var minimizingECalc = new EnergyCalculator.Builder(paramsAndStrands.complexConfSpace, paramsAndStrands.forcefieldParams)
                .setParallelism(delegate.getParallelism())
                .build();
        var rigidECalc = new EnergyCalculator.SharedBuilder(minimizingECalc)
                .setIsMinimizing(false)
                .build();

        var epsilon = delegate.epsilon > 0 ? delegate.epsilon : 0.999999;
        var kstar = new KStar(paramsAndStrands.proteinConfSpace, paramsAndStrands.ligandConfSpace, paramsAndStrands.complexConfSpace, makeKStarSettings(epsilon));

        // ----- Lazy-build training conf spaces for GNN RC mapping (Strategy 8) -----
        SimpleConfSpace trainProteinCS = null;
        SimpleConfSpace trainComplexCS = null;
        boolean useGnnS8 = useBranchMarkstar && gnnStrategy == 8;
        if (useGnnS8) {
            try {
                int[] hotspots = Arrays.stream(gnnTrainHotspots.split(","))
                        .mapToInt(s -> Integer.parseInt(s.trim())).toArray();
                var train = buildTrainingConfSpaces(
                        new File(gnnTrainProteinPdb), new File(gnnTrainLigandPdb),
                        new File(gnnTrainPrepi), new File(gnnTrainFrcmod),
                        new File(gnnTrainTc), new File(gnnTrainRot),
                        gnnTrainChain, hotspots);
                trainProteinCS = train[0];
                trainComplexCS = train[1];
                System.out.println("[GNN-S8] Training conf spaces built: protein="
                        + trainProteinCS.positions.size() + "pos, complex="
                        + trainComplexCS.positions.size() + "pos");
            } catch (Exception e) {
                throw new RuntimeException("[GNN-S8] Failed to build training conf spaces: " + e.getMessage(), e);
            }
        }
        final SimpleConfSpace finalTrainProteinCS = trainProteinCS;
        final SimpleConfSpace finalTrainComplexCS = trainComplexCS;

        for (var info : kstar.confSpaceInfos()) {
            var referenceEnergies = new SimpleReferenceEnergies.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .build();
            var minimizingConfECalc = new ConfEnergyCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .setReferenceEnergies(referenceEnergies)
                    .build();
            info.confEcalc = minimizingConfECalc;

            var minimizedEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(minimizingConfECalc)
                    .build()
                    .calcEnergyMatrix();

            if (useBranchMarkstar) {
                var rigidConfECalc = new ConfEnergyCalculator(info.confEcalc, rigidECalc);
                var rigidEnergymatrix = new SimplerEnergyMatrixCalculator.Builder(rigidConfECalc)
                        .build()
                        .calcEnergyMatrix();

                final boolean isProtein = info.type == KStar.ConfSpaceType.Protein;
                final boolean isComplex = info.type == KStar.ConfSpaceType.Complex;
                final EnergyMatrix capturedEmatMin = minimizedEnergyMatrix;

                info.pfuncFactory = (rcs) -> {
                    var pfn = new BranchMARKStarBound(minimizingConfECalc.confSpace, rigidEnergymatrix, minimizedEnergyMatrix, info.confEcalc, rcs, minimizingECalc.parallelism);
                    pfn.setCorrections(new UpdatingEnergyMatrix(info.confEcalc.confSpace, minimizedEnergyMatrix, info.confEcalc));

                    if (useGnnS8 && (isProtein || isComplex)) {
                        String leafModelPath  = isProtein ? gnnProteinLeaf    : gnnComplexLeaf;
                        String subtreeModelPath = isProtein ? gnnProteinSubtree : gnnComplexSubtree;
                        int leafGpu     = isProtein ? gnnLeafGpuProtein     : gnnLeafGpuComplex;
                        int subtreeGpu  = isProtein ? gnnSubtreeGpuProtein  : gnnSubtreeGpuComplex;
                        SimpleConfSpace inferenceCS = (SimpleConfSpace) info.confSpace;
                        SimpleConfSpace trainCS = isProtein ? finalTrainProteinCS : finalTrainComplexCS;

                        if (leafModelPath != null) {
                            var leafGNN = new GNNConfEnergyCalculator(
                                    new File(leafModelPath), capturedEmatMin,
                                    inferenceCS.positions.size(), leafGpu);
                            leafGNN.setRCMapping(GNNConfEnergyCalculator.buildRCMapping(inferenceCS, trainCS));
                            pfn.setGNNBatchCalculator(leafGNN);
                            pfn.setGPUBatchSize(gnnGpuBatch);
                            pfn.setCPParams(0.001, 0.10, 0.06);
                            System.out.println("[GNN-S8] " + (isProtein ? "protein" : "complex")
                                    + " leaf GNN attached: " + leafModelPath + " (gpu " + leafGpu + ")");
                        }
                        if (subtreeModelPath != null) {
                            var subtreeGNN = new GNNSubtreeEnergyCalculator(
                                    new File(subtreeModelPath), capturedEmatMin,
                                    inferenceCS.positions.size(), subtreeGpu);
                            subtreeGNN.setRCMapping(GNNConfEnergyCalculator.buildRCMapping(inferenceCS, trainCS));
                            pfn.setSubtreeGNN(subtreeGNN);
                            System.out.println("[GNN-S8] " + (isProtein ? "protein" : "complex")
                                    + " subtree GNN attached: " + subtreeModelPath + " (gpu " + subtreeGpu + ")");
                        }
                    }
                    return pfn;
                };
            } else if (useMarkstar) {
                var rigidConfECalc = new ConfEnergyCalculator(info.confEcalc, rigidECalc);
                var rigidEnergymatrix = new SimplerEnergyMatrixCalculator.Builder(rigidConfECalc)
                        .build()
                        .calcEnergyMatrix();

                info.pfuncFactory = (rcs) -> {
                    var pfn = new MARKStarBoundFastQueues(minimizingConfECalc.confSpace, rigidEnergymatrix, minimizedEnergyMatrix, info.confEcalc, rcs, minimizingECalc.parallelism);
                    pfn.setCorrections(new UpdatingEnergyMatrix(info.confEcalc.confSpace, minimizedEnergyMatrix, info.confEcalc));
                    return pfn;
                };
            } else {
                info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                        info.confEcalc,
                        new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                        new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                        rcs.getNumConformations()
                );
            }
        }

        printResults(kstar.run(minimizingECalc.tasks));
        return Main.Success;
    }

    // I don't like this, think of a better way to encapsulate this
    public static class ForceFieldParamsAndStrands {
        private final SimpleConfSpace complexConfSpace;
        public ForcefieldParams forcefieldParams;
        public SimpleConfSpace proteinConfSpace;
        public SimpleConfSpace ligandConfSpace;
        public List<Strand> allStrands;

        private static void addStrandsToComplexBuilder(MoleculeDto dto, SimpleConfSpace.Builder builder, List<Strand> strands) {
            if (dto.translateRotate != null) {
                var tr = new StrandFlex.TranslateRotate(
                        dto.translateRotate.rotateDegrees,
                        dto.translateRotate.translateAngstroms
                );

                StreamEx.of(strands).forEach(cs -> builder.addStrand(cs, tr));
            } else {
                StreamEx.of(strands).forEach(builder::addStrand);
            }
        }

        public ForceFieldParamsAndStrands(DesignFileDelegate delegate, AffinityDesign design) {

            forcefieldParams = delegate.frcmodPath == null
                    ? new ForcefieldParams()
                    : new ForcefieldParams(ForcefieldParams.Forcefield.AMBER, new ForcefieldFileParser(
                            getClass().getResourceAsStream(ForcefieldParams.Forcefield.AMBER.paramsPath), Paths.get(delegate.frcmodPath))
            );

            proteinConfSpace = delegate.createConfSpace(design.protein, forcefieldParams);
            ligandConfSpace = delegate.createConfSpace(design.ligand, forcefieldParams);

            allStrands = StreamEx.of(proteinConfSpace.strands, ligandConfSpace.strands)
                    .flatMap(List::stream)
                    .toList();

            var builder = new SimpleConfSpace.Builder();
            addStrandsToComplexBuilder(design.protein, builder, proteinConfSpace.strands);
            addStrandsToComplexBuilder(design.ligand, builder, ligandConfSpace.strands);
            complexConfSpace = builder.build();
        }

        List<Residue> makeAllResidues() {
            return StreamEx.of(allStrands)
                    .flatMap(x -> x.mol.residues.stream())
                    .toList();
        }
    }

    private int makeScanDesigns(AffinityDesign design, List<Residue> allResidues) {
        var target = design.scanSettings.target;
        var residues = design.scanSettings.residues;

        if (target.isEmpty() && residues.isEmpty() || !target.isEmpty() && !residues.isEmpty()) {
            System.err.println("Either target or residues must be specified, but not both");
            return Main.Failure;
        }

        residues.forEach(r -> r.mutability = r.mutability.isEmpty() ? allMutableAaTypes : r.mutability);
        var mutableTargets = residues.isEmpty()
                ? findMutableResiduesAroundTarget(design, scanDistance, target, allResidues)
                : residues;

        return createScanDesigns(design, mutableTargets, scanFlexDistance);
    }

    @NotNull
    private List<ResidueModifier> findMutableResiduesAroundTarget(AffinityDesign design, double dist, String target, List<Residue> allResidues) {
        var targetRes = allResidues.stream()
                .filter(a -> a.getPDBResNumber().equals(target))
                .findFirst()
                .orElseThrow();

        var modifiers = StreamEx.of(allResidues)
                .filter(x -> x != targetRes)
                .filter(x -> !design.scanSettings.excluding.contains(x.getPDBResNumber()))
                .filter(x -> x.distanceTo(targetRes) <= dist)
                .map(CommandBindingAffinity::makeFlexibleResidueModifier)
                .toList();

        modifiers.forEach(a -> a.mutability = allMutableAaTypes);
        return modifiers;
    }

    public static Stream<Residue> nearbyResidues(Residue target, Molecule molecule, double withinDistance) {
        return molecule.residues.stream()
                .filter(x -> x.distanceTo(target) <= withinDistance)
                .filter(x -> target.getChainId() != x.getChainId()
                        || !target.getPDBResNumber().equals(x.getPDBResNumber())
                        || !target.getType().equals(x.getType())
                );
    }

    public static Stream<ResidueModifier> nearbyResidueModifiers(Residue target, Molecule molecule, double withinDistance) {
        return nearbyResidues(target, molecule, withinDistance).map(CommandBindingAffinity::makeFlexibleResidueModifier);
    }

    private int createScanDesigns(AffinityDesign designTemplate, List<ResidueModifier> mutableTargets, double flexDist) {

        var proteinMol = designTemplate.makeProteinMolecule();
        var ligandMol = designTemplate.makeLigandMolecule();

        for (var comboIndices : StreamEx.ofCombinations(mutableTargets.size(), delegate.maxSimultaneousMutations)) {

            var mutableResidues = IntStreamEx.of(comboIndices)
                    .mapToObj(mutableTargets::get)
                    .toList();
            var mutableTargetNums = StreamEx.of(mutableResidues)
                    .map(x -> x.identity.positionIdentifier())
                    .toSet();
            var mutableTargetRes = StreamEx.of(proteinMol.residues)
                    .append(ligandMol.residues)
                    .filter(res -> mutableTargetNums.contains(res.getPDBResNumber()))
                    .toList();

            List<ResidueModifier> proteinResMods = StreamEx.of(mutableTargetRes)
                    .flatMap(residue -> nearbyResidueModifiers(residue, proteinMol, flexDist))
                    .distinct()
                    .toList();

            var ligandResMods = StreamEx.of(mutableTargetRes)
                    .flatMap(residue -> nearbyResidueModifiers(residue, ligandMol, flexDist))
                    .distinct()
                    .toList();

            StreamEx.of(mutableResidues)
                    .zipWith(StreamEx.of(mutableTargetRes))
                    .forKeyValue((residueModifier, residue) -> {
                        if (proteinMol.residues.contains(residue)) {
                            proteinResMods.add(residueModifier);
                        } else {
                            ligandResMods.add(residueModifier);
                        }
                    });

            try {
                var nameTemp = delegate.design.getName().substring(0, delegate.design.getName().lastIndexOf("."));
                var comboName = String.join("-", StreamEx.of(mutableTargetRes)
                        .map(Residue::getPDBResNumber)
                        .toList());

                var path = Path.of(String.format("%s.%s.yaml", nameTemp, comboName));
                Files.deleteIfExists(path);
                var outFile = Files.createFile(path);

                var designCopy = designTemplate.copy();
                designCopy.protein.residueModifiers = StreamEx.of(designCopy.protein.residueModifiers).append(proteinResMods).toList();
                designCopy.ligand.residueModifiers = StreamEx.of(designCopy.ligand.residueModifiers).append(ligandResMods).toList();
                designCopy.scanSettings = null; // no need to copy this into the newly-created design
                designCopy.write(outFile);
            } catch (IOException e) {
                e.printStackTrace();
                return Main.Failure;
            }
        }

        return Main.Success;
    }

    @NotNull
    private static ResidueModifier makeFlexibleResidueModifier(Residue x) {
        var res = new edu.duke.cs.osprey.design.models.Residue();
        res.aminoAcidType = x.getType();
        res.chain = String.valueOf(x.getChainId());
        res.residueNumber = Integer.parseInt(x.getPDBResNumber().substring(1));

        var modifier = new ResidueModifier();
        modifier.identity = res;
        modifier.mutability = List.of();

        return modifier;
    }

    private void printResults(List<KStar.ScoredSequence> results) {
        // nop
    }

    private KStar.Settings makeKStarSettings(double epsilon) {
        var builder = new KStar.Settings.Builder();
        builder.setEpsilon(epsilon);
        builder.addScoreConsoleWriter();
        if (maxNumberConfs > 0) {
            builder.setMaxNumConf(maxNumberConfs);
        }
        builder.setMaxSimultaneousMutations(delegate.maxSimultaneousMutations);
        builder.setStabilityThreshold(stabilityThreshold < 0 ? null : stabilityThreshold);

        if (delegate.saveResultsToDb) { // assuming the validation has already been done at this point.
            String designFile = "";
            try {
                designFile = Files.readString(delegate.design.toPath());
            } catch (IOException e) {
                e.printStackTrace();
            }

            var commandLineArgs = String.join(" ", args);

            var connectionString = String.format("jdbc:postgresql://%s:%s/%s", dbSettings.get(dbHostnameKey), dbSettings.get(dbPortKey), dbSettings.get(dbNameKey));
            var pgsqlConnInfo = new PostgresConnectionInfo(dbSettings.get(dbUserNameKey), dbSettings.get(dbPasswordKey), connectionString);
            var s3ConnInfo = new S3Settings("us-east-1", "duke-osprey"); // TODO: get this from config
            var dbScoreWriter = new PostgresScoreWriter(pgsqlConnInfo, s3ConnInfo, delegate.design.getName(), List.of(designFile, commandLineArgs), delegate.numConfs);

            builder.addScoreWriter(dbScoreWriter);
        } else if (delegate.numConfs > 0) {
            var saveDir = delegate.saveDir;
            var scoreWriter = new StructureFileScoreWriter(saveDir, delegate.numConfs);
            builder.addScoreWriter(scoreWriter);
        }

        return builder.build();
    }

    private Optional<AffinityDesign> parseAndValidate(File designSpec) {

        if (delegate.saveResultsToDb && delegate.propertiesFile == null) {
            System.err.println("You requested to save the design to a database, but a properties file containing the database settings was not specified. Exiting.");
            return Optional.empty();
        }

        if (delegate.saveResultsToDb) {
            try(var is = new FileInputStream(delegate.propertiesFile)) {
                var props = new Properties();
                props.load(is);

                var expecting = List.of(dbHostnameKey, dbPortKey, dbNameKey, dbUserNameKey, dbPasswordKey);

                expecting.stream()
                        .filter(key -> !props.containsKey(key))
                        .forEach(key -> System.out.printf("The properties file must have a key `%s` for configuring the database but does not.%n", key));

                if (expecting.stream().anyMatch(key -> !props.containsKey(key))) {
                    System.err.println("Properties file does not have all required properties, exiting.");
                    return Optional.empty();
                }

                dbSettings = expecting.stream().collect(Collectors.toMap(key -> key, props::getProperty));

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return getAffinityDesignFromFile(designSpec);
    }

    @NotNull
    public static Optional<AffinityDesign> getAffinityDesignFromFile(File designSpec) {
        AffinityDesign design;
        try {
            design = AffinityDesign.parse(designSpec);
        } catch (IOException e) {
            e.printStackTrace();
            return Optional.empty();
        }

        var specErrors = design.validate();
        if (!specErrors.isEmpty()) {
            System.err.println("Invalid design specification. The following validations failed:");
            specErrors.stream().map(s -> String.format("- %s", s)).forEach(System.err::println);
            return Optional.empty();
        }

        return Optional.of(design);
    }

    /**
     * Build training {protein, complex} conf spaces matching the GNN training spec.
     * Replicates TestBranchMARKStar.buildConfSpaceEGFR — same hotspot residues, all20, continuous.
     * Returns array [protein, complex].
     */
    private static SimpleConfSpace[] buildTrainingConfSpaces(
            File proteinPdb, File ligandPdb,
            File prepiFile, File frcmodFile, File tcFile, File rotFile,
            String chain, int[] hotspotResNums) throws IOException {

        ForcefieldParams.Forcefield ff = ForcefieldParams.Forcefield.AMBER;
        ForcefieldFileParser parser = (frcmodFile != null && frcmodFile.exists())
                ? new ForcefieldFileParser(ForcefieldParams.class.getResourceAsStream(ff.paramsPath), frcmodFile.toPath())
                : new ForcefieldFileParser(ForcefieldParams.class.getResourceAsStream(ff.paramsPath));
        ForcefieldParams ffparams = new ForcefieldParams(ff, parser);

        String prepiContent = Files.readString(prepiFile.toPath());
        String tcContent  = (tcFile != null && tcFile.exists())  ? Files.readString(tcFile.toPath())  : "";
        String rotContent = (rotFile != null && rotFile.exists()) ? Files.readString(rotFile.toPath()) : "";

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
                .addTemplates(prepiContent)
                .addTemplateCoords(tcContent)
                .addRotamers(rotContent)
                .build();

        Molecule proteinMol = PDBIO.readFile(proteinPdb);
        Molecule ligandMol  = PDBIO.readFile(ligandPdb);

        // Protein residue range
        int firstRes = Integer.MAX_VALUE, lastRes = Integer.MIN_VALUE;
        for (Residue r : proteinMol.residues) {
            String s = Residues.normalizeResNum(r.getPDBResNumber());
            int i = 0;
            while (i < s.length() && !Character.isDigit(s.charAt(i))) i++;
            int rn = Integer.parseInt(s.substring(i).replaceAll("[^0-9-].*$", ""));
            if (rn < firstRes) firstRes = rn;
            if (rn > lastRes)  lastRes  = rn;
        }
        Strand protein = new Strand.Builder(proteinMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(chain + firstRes, chain + lastRes)
                .build();

        String[] all20 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
                          "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
        for (int rn : hotspotResNums) {
            String id = chain + rn;
            Strand.ResidueFlex rf = protein.flexibility.get(id);
            if (rf != null) rf.setLibraryRotamers(all20).addWildTypeRotamers().setContinuous();
        }

        // Ligand strand (rigid)
        int ligFirst = Integer.MAX_VALUE, ligLast = Integer.MIN_VALUE;
        String ligChain = chain;
        for (Residue r : ligandMol.residues) {
            try {
                String s = Residues.normalizeResNum(r.getPDBResNumber());
                int i = 0;
                while (i < s.length() && !Character.isDigit(s.charAt(i))) i++;
                int rn = Integer.parseInt(s.substring(i).replaceAll("[^0-9-].*$", ""));
                if (rn < ligFirst) { ligFirst = rn; ligChain = String.valueOf(r.getChainId()); }
                if (rn > ligLast)  ligLast  = rn;
            } catch (NumberFormatException ignored) {}
        }
        Strand ligand = new Strand.Builder(ligandMol)
                .setTemplateLibrary(templateLib)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
                .setResidues(ligChain + ligFirst, ligChain + ligLast)
                .build();

        SimpleConfSpace proteinCS = new SimpleConfSpace.Builder().addStrand(protein).build();
        SimpleConfSpace complexCS = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();
        return new SimpleConfSpace[] { proteinCS, complexCS };
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }
}
