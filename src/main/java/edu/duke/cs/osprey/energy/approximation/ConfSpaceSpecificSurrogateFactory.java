package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.branch.BranchResidualDatasetBuilder;
import edu.duke.cs.osprey.energy.approximation.branch.BranchResidualFeatureEncoder;
import edu.duke.cs.osprey.energy.approximation.branch.BranchSharedRidgeModel;
import edu.duke.cs.osprey.energy.approximation.branch.BranchSparseSampler;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.lute.BranchResidualConfEnergyCalculator;
import edu.duke.cs.osprey.lute.ConfSampler;
import edu.duke.cs.osprey.lute.LUTE;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEIO;
import edu.duke.cs.osprey.lute.LUTEState;
import edu.duke.cs.osprey.lute.RandomizedDFSConfSampler;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.branchdp.RootedTreeEdge;
import edu.duke.cs.osprey.minimization.ApproximationObjectiveFunction;
import edu.duke.cs.osprey.minimization.GlobalMLPSurrogateObjectiveFunction;
import edu.duke.cs.osprey.minimization.MLPSurrogateObjectiveFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.io.File;

/**
 * Factory helpers for per-confspace surrogate training/loading and objective creation.
 */
public class ConfSpaceSpecificSurrogateFactory {

    public static class BranchResidualRidgeConfig {
        public File cacheRoot;
        public String taskTag = "default";
        public boolean forceRetrain = false;

        public int featureHashDim = 65536;
        public String samplerMode = "sparse";
        public int pretrainGridDPSamples = 20000;
        public int finetuneCCDSamples = 2000;
        public double lambdaPre = 1e-2;
        public double lambdaFine = 1e-1;

        public boolean pretrainUseForcefieldFinalEval = true;
        public int gridSize = 2;
        public int gridThreads = 1;
        public int lowEnergyPoolLimit = 4096;

        public int randomSeed = 12345;
        public int ridgeMaxIterations = 256;
        public double ridgeTolerance = 1e-6;
    }

    /**
     * Build a ConfEnergyCalculator with a task-specific approximator matrix attached.
     */
    public static ConfEnergyCalculator withTaskSpecificApproximator(
            ConfEnergyCalculator baseConfEcalc,
            File cacheRoot,
            String taskTag,
            int numSamplesPerParam,
            double approximationErrorBudget
    ) {

        ConfSpaceSpecificApproximatorCache.Config cfg = new ConfSpaceSpecificApproximatorCache.Config();
        cfg.cacheRoot = cacheRoot;
        cfg.taskTag = taskTag;
        cfg.numSamplesPerParam = numSamplesPerParam;
        cfg.type = ApproximatorMatrixCalculator.ApproximatorType.Quadratic;

        ApproximatorMatrix amat = ConfSpaceSpecificApproximatorCache.loadOrTrain(baseConfEcalc, cfg);

        return new ConfEnergyCalculator.Builder(baseConfEcalc.confSpace, baseConfEcalc.ecalc)
                .setEnergyPartition(baseConfEcalc.epart)
                .setReferenceEnergies(baseConfEcalc.eref)
                .addResEntropy(baseConfEcalc.addResEntropy)
                .setAddShellInters(baseConfEcalc.addShellInters)
                .setApproximatorMatrix(amat)
                .setApproximationErrorBudget(approximationErrorBudget)
                .build();
    }

    /**
     * Create a CCD objective function that uses only the task-specific surrogate (no forcefield calls).
     */
    public static ApproximationObjectiveFunction makeApproximationObjective(
            ConfEnergyCalculator confEcalcWithApproximator,
            RCTuple tuple,
            boolean requireFullApproximation
    ) {
        if (confEcalcWithApproximator.amat == null) {
            throw new IllegalArgumentException("ConfEnergyCalculator has no approximator matrix");
        }

        ResidueInteractions inters = confEcalcWithApproximator.makeFragInters(tuple);
        ResidueInteractionsApproximator approximator = confEcalcWithApproximator.amat.get(
                tuple, inters, confEcalcWithApproximator.approximationErrorBudget
        );

        if (requireFullApproximation && approximator.ffInters.size() > 0) {
            throw new IllegalStateException(String.format(
                    "Approximation left %d forcefield interactions for tuple %s",
                    approximator.ffInters.size(), tuple.stringListing()
            ));
        }

        ParametricMolecule pmol = confEcalcWithApproximator.confSpace.makeMolecule(tuple);
        return new ApproximationObjectiveFunction(pmol, approximator.approximator);
    }

    /**
     * Train or load a confspace-specific MLP surrogate matrix for one-body and pairwise energies.
     */
    public static MLPSurrogateMatrix loadOrTrainTaskSpecificMLPSurrogate(
            ConfEnergyCalculator baseConfEcalc,
            File cacheRoot,
            String taskTag,
            int numSamplesPerParam,
            int minSamplesPerModel,
            int maxSamplesPerModel,
            int hidden1,
            int hidden2,
            int epochs,
            int batchSize,
            double learningRate
    ) {
        ConfSpaceSpecificMLPSurrogateCache.Config cfg = new ConfSpaceSpecificMLPSurrogateCache.Config();
        cfg.cacheRoot = cacheRoot;
        cfg.taskTag = taskTag;
        cfg.numSamplesPerParam = numSamplesPerParam;
        cfg.minSamplesPerModel = minSamplesPerModel;
        cfg.maxSamplesPerModel = maxSamplesPerModel;
        cfg.hidden1 = hidden1;
        cfg.hidden2 = hidden2;
        cfg.epochs = epochs;
        cfg.batchSize = batchSize;
        cfg.learningRate = learningRate;
        return ConfSpaceSpecificMLPSurrogateCache.loadOrTrain(baseConfEcalc, cfg);
    }

    /**
     * Train or load a single confspace-specific global MLP surrogate for full tuple energies.
     */
    public static TaskGlobalMLPSurrogate loadOrTrainTaskSpecificGlobalMLPSurrogate(
            ConfEnergyCalculator baseConfEcalc,
            File cacheRoot,
            String taskTag,
            int numTupleSamples,
            int samplesPerTuple,
            int hidden1,
            int hidden2,
            int epochs,
            int batchSize,
            double learningRate
    ) {
        ConfSpaceSpecificGlobalMLPSurrogateCache.Config cfg = new ConfSpaceSpecificGlobalMLPSurrogateCache.Config();
        cfg.cacheRoot = cacheRoot;
        cfg.taskTag = taskTag;
        cfg.numTupleSamples = numTupleSamples;
        cfg.samplesPerTuple = samplesPerTuple;
        cfg.hidden1 = hidden1;
        cfg.hidden2 = hidden2;
        cfg.epochs = epochs;
        cfg.batchSize = batchSize;
        cfg.learningRate = learningRate;
        return ConfSpaceSpecificGlobalMLPSurrogateCache.loadOrTrain(baseConfEcalc, cfg);
    }

    /**
     * Train or load a task-specific branch-shared ridge residual model and return a conf-level calculator.
     */
    public static BranchResidualConfEnergyCalculator loadOrTrainTaskSpecificBranchResidualRidge(
            ConfEnergyCalculator confEcalcMinimized,
            EnergyMatrix ematMinimized,
            InteractionGraph interactionGraph,
            RootedTreeEdge rootEdge,
            ForcefieldParams ffparams,
            SimpleReferenceEnergies eref,
            BranchResidualRidgeConfig config
    ) {

        if (config == null) {
            throw new IllegalArgumentException("config must not be null");
        }
        if (config.cacheRoot == null) {
            throw new IllegalArgumentException("config.cacheRoot must not be null");
        }

        SimpleConfSpace confSpace = confEcalcMinimized.confSpace;
        String confHash = ConfSpaceSpecificApproximatorCache.fingerprint(confSpace);

        File taskDir = new File(config.cacheRoot, sanitize(config.taskTag));
        if (!taskDir.exists() && !taskDir.mkdirs()) {
            throw new RuntimeException("Failed to create branch residual cache directory: " + taskDir.getAbsolutePath());
        }

        File luteFile = new File(taskDir, confHash + ".branchResidual.lute");
        File ridgeFile = new File(taskDir, confHash + ".branchResidual.ridge");

        if (config.forceRetrain) {
            if (luteFile.exists() && !luteFile.delete()) {
                throw new RuntimeException("Failed to delete old cached LUTE file: " + luteFile.getAbsolutePath());
            }
            if (ridgeFile.exists() && !ridgeFile.delete()) {
                throw new RuntimeException("Failed to delete old cached ridge file: " + ridgeFile.getAbsolutePath());
            }
        }

        LUTEState luteState;
        if (luteFile.exists()) {
            luteState = LUTEIO.read(luteFile);
        } else {
            luteState = trainAndCacheLuteState(confEcalcMinimized, ematMinimized, luteFile, config.randomSeed);
        }

        LUTEConfEnergyCalculator baseLute = new LUTEConfEnergyCalculator(confSpace, luteState);
        BranchSharedRidgeModel ridgeModel;

        if (ridgeFile.exists()) {
            ridgeModel = BranchSharedRidgeModel.readFromFile(ridgeFile);
            if (ridgeModel.hashDim != config.featureHashDim) {
                ridgeModel = trainAndCacheBranchResidualModel(
                        confEcalcMinimized,
                        baseLute,
                        ematMinimized,
                        interactionGraph,
                        rootEdge,
                        ffparams,
                        eref,
                        config,
                        ridgeFile
                );
            }
        } else {
            ridgeModel = trainAndCacheBranchResidualModel(
                    confEcalcMinimized,
                    baseLute,
                    ematMinimized,
                    interactionGraph,
                    rootEdge,
                    ffparams,
                    eref,
                    config,
                    ridgeFile
            );
        }

        return makeBranchResidualConfEnergyCalculator(confSpace, luteState, ridgeModel, interactionGraph, rootEdge);
    }

    public static BranchResidualConfEnergyCalculator makeBranchResidualConfEnergyCalculator(
            SimpleConfSpace confSpace,
            LUTEState luteState,
            BranchSharedRidgeModel ridgeModel,
            InteractionGraph interactionGraph,
            RootedTreeEdge rootEdge
    ) {
        LUTEConfEnergyCalculator baseLute = new LUTEConfEnergyCalculator(confSpace, luteState);
        BranchResidualFeatureEncoder encoder = new BranchResidualFeatureEncoder(
                confSpace,
                interactionGraph,
                rootEdge,
                ridgeModel.hashDim
        );
        return new BranchResidualConfEnergyCalculator(confSpace, baseLute, encoder, ridgeModel);
    }

    private static LUTEState trainAndCacheLuteState(
            ConfEnergyCalculator confEcalcMinimized,
            EnergyMatrix ematMinimized,
            File luteFile,
            int randomSeed
    ) {

        SimpleConfSpace confSpace = confEcalcMinimized.confSpace;
        PruningMatrix pmat = new PruningMatrix(confSpace);

        try (ConfDB confDB = new ConfDB(confSpace)) {
            ConfDB.ConfTable confTable = confDB.new ConfTable("branchResidualLUTE");

            LUTE lute = new LUTE(confSpace);
            ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);

            lute.sampleTuplesAndFit(
                    confEcalcMinimized,
                    ematMinimized,
                    pmat,
                    confTable,
                    sampler,
                    LUTE.Fitter.OLSCG,
                    1.5,
                    0.1
            );

            if (lute.getTrainingSystem() == null) {
                throw new RuntimeException("LUTE training did not produce a training system");
            }

            LUTEState state = new LUTEState(lute.getTrainingSystem());
            LUTEIO.write(state, luteFile);
            return state;
        }
    }

    private static BranchSharedRidgeModel trainAndCacheBranchResidualModel(
            ConfEnergyCalculator confEcalcMinimized,
            LUTEConfEnergyCalculator baseLute,
            EnergyMatrix ematMinimized,
            InteractionGraph interactionGraph,
            RootedTreeEdge rootEdge,
            ForcefieldParams ffparams,
            SimpleReferenceEnergies eref,
            BranchResidualRidgeConfig config,
            File ridgeFile
    ) {

        PruningMatrix pmat = new PruningMatrix(confEcalcMinimized.confSpace);

        BranchResidualFeatureEncoder encoder = new BranchResidualFeatureEncoder(
                confEcalcMinimized.confSpace,
                interactionGraph,
                rootEdge,
                config.featureHashDim
        );

        BranchSparseSampler sampler = new BranchSparseSampler(
                confEcalcMinimized.confSpace,
                interactionGraph,
                pmat,
                config.randomSeed
        );

        BranchResidualDatasetBuilder datasetBuilder = new BranchResidualDatasetBuilder(
                confEcalcMinimized,
                baseLute,
                encoder,
                sampler,
                interactionGraph,
                rootEdge,
                ffparams,
                eref,
                ematMinimized,
                pmat
        )
                .setSamplerMode(config.samplerMode)
                .setGridSize(config.gridSize)
                .setGridThreads(config.gridThreads)
                .setLowEnergyPoolLimit(config.lowEnergyPoolLimit)
                .setPretrainUseForcefieldFinalEval(config.pretrainUseForcefieldFinalEval);

        BranchSharedRidgeModel.TrainingData pretrain = datasetBuilder.buildPretrainData(config.pretrainGridDPSamples);
        BranchSharedRidgeModel preModel = BranchSharedRidgeModel.fit(
                pretrain,
                config.lambdaPre,
                config.ridgeMaxIterations,
                config.ridgeTolerance
        );

        BranchSharedRidgeModel.TrainingData fineTune = datasetBuilder.buildFineTuneData(config.finetuneCCDSamples);
        BranchSharedRidgeModel finalModel = BranchSharedRidgeModel.fitWithPrior(
                fineTune,
                config.lambdaFine,
                preModel,
                config.ridgeMaxIterations,
                config.ridgeTolerance
        );

        finalModel.writeToFile(ridgeFile);
        return finalModel;
    }

    /**
     * Create a CCD objective function backed only by MLP one-body/pair surrogates.
     */
    public static MLPSurrogateObjectiveFunction makeMLPSurrogateObjective(
            ConfEnergyCalculator confEcalc,
            MLPSurrogateMatrix surrogate,
            RCTuple tuple,
            boolean requireFullApproximation
    ) {
        double erefOffset = 0.0;
        if (confEcalc.eref != null) {
            for (int i = 0; i < tuple.pos.size(); i++) {
                erefOffset += confEcalc.eref.getOffset(confEcalc.confSpace, tuple.pos.get(i), tuple.RCs.get(i));
            }
        }
        return new MLPSurrogateObjectiveFunction(
                confEcalc.confSpace, tuple, surrogate, requireFullApproximation, erefOffset
        );
    }

    /**
     * Create a CCD objective function backed by a single global MLP surrogate.
     */
    public static GlobalMLPSurrogateObjectiveFunction makeGlobalMLPSurrogateObjective(
            ConfEnergyCalculator confEcalc,
            TaskGlobalMLPSurrogate surrogate,
            RCTuple tuple
    ) {
        return new GlobalMLPSurrogateObjectiveFunction(confEcalc.confSpace, tuple, surrogate);
    }

    private static String sanitize(String s) {
        if (s == null || s.trim().isEmpty()) {
            return "default";
        }
        return s.replaceAll("[^a-zA-Z0-9._-]", "_");
    }
}
