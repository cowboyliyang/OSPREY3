package edu.duke.cs.osprey.energy.approximation.branch;

import cern.colt.matrix.DoubleFactory1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.markstar.framework.GridDPMinimizer;
import edu.duke.cs.osprey.markstar.framework.branch.InteractionGraph;
import edu.duke.cs.osprey.markstar.framework.branch.RootedTreeEdge;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Builds residual training datasets:
 *   y = E_target(conf) - E_base(conf)
 *
 * Pretrain labels: GridDP (optionally FF final eval at GridDP point)
 * Fine-tune labels: CCD minimized forcefield energy
 */
public class BranchResidualDatasetBuilder {

    private final SimpleConfSpace confSpace;
    private final ConfEnergyCalculator confEcalcMinimized;
    private final LUTEConfEnergyCalculator baseConfEcalc;
    private final BranchResidualFeatureEncoder encoder;
    private final BranchSparseSampler sampler;

    private final InteractionGraph interactionGraph;
    private final RootedTreeEdge rootEdge;
    private final ForcefieldParams ffparams;
    private final SimpleReferenceEnergies eref;
    private final EnergyMatrix ematMinimized;

    private String samplerMode = "sparse";
    private BranchSparseSampler.Split sparseSplit = BranchSparseSampler.Split.defaults();
    private int lowEnergyPoolLimit = 4096;

    private int gridSize = 2;
    private int gridThreads = 1;
    private boolean pretrainUseForcefieldFinalEval = true;

    public BranchResidualDatasetBuilder(
            ConfEnergyCalculator confEcalcMinimized,
            LUTEConfEnergyCalculator baseConfEcalc,
            BranchResidualFeatureEncoder encoder,
            BranchSparseSampler sampler,
            InteractionGraph interactionGraph,
            RootedTreeEdge rootEdge,
            ForcefieldParams ffparams,
            SimpleReferenceEnergies eref,
            EnergyMatrix ematMinimized,
            PruningMatrix pmat
    ) {
        this.confSpace = confEcalcMinimized.confSpace;
        this.confEcalcMinimized = confEcalcMinimized;
        this.baseConfEcalc = baseConfEcalc;
        this.encoder = encoder;
        this.sampler = sampler;

        this.interactionGraph = interactionGraph;
        this.rootEdge = rootEdge;
        this.ffparams = ffparams;
        this.eref = eref;
        this.ematMinimized = ematMinimized;
    }

    public BranchResidualDatasetBuilder setSamplerMode(String mode) {
        if (mode != null && !mode.trim().isEmpty()) {
            this.samplerMode = mode.trim().toLowerCase(Locale.ROOT);
        }
        return this;
    }

    public BranchResidualDatasetBuilder setSparseSplit(BranchSparseSampler.Split split) {
        if (split != null) {
            this.sparseSplit = split;
        }
        return this;
    }

    public BranchResidualDatasetBuilder setLowEnergyPoolLimit(int val) {
        this.lowEnergyPoolLimit = Math.max(1, val);
        return this;
    }

    public BranchResidualDatasetBuilder setGridSize(int gridSize) {
        this.gridSize = Math.max(2, gridSize);
        return this;
    }

    public BranchResidualDatasetBuilder setGridThreads(int gridThreads) {
        this.gridThreads = Math.max(1, gridThreads);
        return this;
    }

    public BranchResidualDatasetBuilder setPretrainUseForcefieldFinalEval(boolean val) {
        this.pretrainUseForcefieldFinalEval = val;
        return this;
    }

    public BranchSharedRidgeModel.TrainingData buildPretrainData(int numSamples) {

        List<int[]> confs = sampleConfs(numSamples);
        List<BranchResidualFeatureEncoder.SparseFeatures> features = new ArrayList<>(confs.size());
        List<Double> targets = new ArrayList<>(confs.size());

        GridDPMinimizer gridDP = new GridDPMinimizer(
                confSpace,
                interactionGraph,
                rootEdge,
                gridSize,
                ffparams,
                eref,
                true,
                gridThreads
        );

        for (int[] conf : confs) {
            double baseEnergy = baseConfEcalc.calcEnergy(conf);
            GridDPMinimizer.Result result = gridDP.minimize(conf);

            double target = result.energy;
            if (pretrainUseForcefieldFinalEval && result.bestDOFValues != null) {
                target = evaluateForcefieldEnergyAt(conf, result.bestDOFValues);
            }

            double residual = target - baseEnergy;
            features.add(encoder.encode(conf));
            targets.add(residual);
        }

        return BranchSharedRidgeModel.TrainingData.fromLists(encoder.getHashDim(), features, targets);
    }

    public BranchSharedRidgeModel.TrainingData buildFineTuneData(int numSamples) {

        List<int[]> confs = sampleConfs(numSamples);
        List<BranchResidualFeatureEncoder.SparseFeatures> features = new ArrayList<>(confs.size());
        List<Double> targets = new ArrayList<>(confs.size());

        for (int[] conf : confs) {
            double baseEnergy = baseConfEcalc.calcEnergy(conf);
            double target = confEcalcMinimized.calcEnergy(new RCTuple(conf)).energy;
            double residual = target - baseEnergy;

            features.add(encoder.encode(conf));
            targets.add(residual);
        }

        return BranchSharedRidgeModel.TrainingData.fromLists(encoder.getHashDim(), features, targets);
    }

    private List<int[]> sampleConfs(int numSamples) {
        int n = Math.max(1, numSamples);

        if ("lute".equalsIgnoreCase(samplerMode)) {
            return sampler.sample(
                    n,
                    new BranchSparseSampler.Split(0.0, 0.0, 1.0),
                    null,
                    lowEnergyPoolLimit
            );
        }

        return sampler.sample(n, sparseSplit, ematMinimized, lowEnergyPoolLimit);
    }

    private double evaluateForcefieldEnergyAt(int[] conf, double[] dofValues) {
        RCTuple tuple = new RCTuple(conf);
        ParametricMolecule pmol = confSpace.makeMolecule(tuple);

        for (int d = 0; d < dofValues.length && d < pmol.dofs.size(); d++) {
            pmol.dofs.get(d).apply(dofValues[d]);
        }

        ResidueInteractions inters = confEcalcMinimized.makeFragInters(tuple);
        try (EnergyFunction efunc = confEcalcMinimized.ecalc.makeEnergyFunction(pmol, inters)) {
            return efunc.getEnergy();
        }
    }
}
