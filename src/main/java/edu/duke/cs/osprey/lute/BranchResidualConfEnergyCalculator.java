package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.branch.BranchResidualFeatureEncoder;
import edu.duke.cs.osprey.energy.approximation.branch.BranchSharedRidgeModel;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;

/**
 * Full-conf calculator:
 *   E_pred(conf) = E_base_lute(conf) + residual_ridge(conf)
 */
public class BranchResidualConfEnergyCalculator extends ConfEnergyCalculator {

    public final LUTEConfEnergyCalculator base;
    public final BranchResidualFeatureEncoder encoder;
    public final BranchSharedRidgeModel model;

    public BranchResidualConfEnergyCalculator(
            SimpleConfSpace confSpace,
            LUTEConfEnergyCalculator base,
            BranchResidualFeatureEncoder encoder,
            BranchSharedRidgeModel model
    ) {
        super(confSpace, new TaskExecutor());
        this.base = base;
        this.encoder = encoder;
        this.model = model;
    }

    private static class NotSupportedByBranchResidualException extends RuntimeException {
        private NotSupportedByBranchResidualException() {
            super("BranchResidualConfEnergyCalculator only supports full-conformation energies");
        }
    }

    @Override
    public ResidueInteractions makeSingleInters(int pos, int rc) {
        throw new NotSupportedByBranchResidualException();
    }

    @Override
    public ResidueInteractions makePairInters(int pos1, int rc1, int pos2, int rc2) {
        throw new NotSupportedByBranchResidualException();
    }

    @Override
    public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf) {
        return new ConfSearch.EnergiedConf(conf, calcEnergy(conf.getAssignments()));
    }

    @Override
    public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf, ResidueInteractions inters) {
        throw new NotSupportedByBranchResidualException();
    }

    @Override
    public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {
        throw new NotSupportedByBranchResidualException();
    }

    @Override
    public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag) {
        int[] conf = Conf.make(confSpace, frag);
        if (!Conf.isCompletelyAssigned(conf)) {
            throw new NotSupportedByBranchResidualException();
        }
        return new EnergyCalculator.EnergiedParametricMolecule((ParametricMolecule) null, null, calcEnergy(conf));
    }

    @Override
    public MoleculeObjectiveFunction makeIntraShellObjFcn(int pos, int rc) {
        throw new NotSupportedByBranchResidualException();
    }

    @Override
    public MoleculeObjectiveFunction makePairwiseObjFcn(int pos1, int rc1, int pos2, int rc2) {
        throw new NotSupportedByBranchResidualException();
    }

    public double calcEnergy(int[] conf) {

        numCalculations.incrementAndGet();

        double baseEnergy = base.calcEnergy(conf);
        BranchResidualFeatureEncoder.SparseFeatures feat = encoder.encode(conf);
        double residual = model.predict(feat);

        return baseEnergy + residual;
    }
}
