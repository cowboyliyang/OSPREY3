package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.minimization.ApproximationObjectiveFunction;

import java.io.File;

/**
 * Factory helpers for per-confspace surrogate training/loading and objective creation.
 */
public class ConfSpaceSpecificSurrogateFactory {

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
}
