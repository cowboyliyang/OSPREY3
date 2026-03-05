package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.minimization.ApproximationObjectiveFunction;
import edu.duke.cs.osprey.minimization.MLPSurrogateObjectiveFunction;

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

    /**
     * Train or load a confspace-specific MLP surrogate matrix for one-body and pairwise energies.
     */
    public static MLPSurrogateMatrix loadOrTrainTaskSpecificMLPSurrogate(
            ConfEnergyCalculator baseConfEcalc,
            File cacheRoot,
            String taskTag,
            int numSamplesPerParam,
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
        cfg.hidden1 = hidden1;
        cfg.hidden2 = hidden2;
        cfg.epochs = epochs;
        cfg.batchSize = batchSize;
        cfg.learningRate = learningRate;
        return ConfSpaceSpecificMLPSurrogateCache.loadOrTrain(baseConfEcalc, cfg);
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
}
