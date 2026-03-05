package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.Random;

import static edu.duke.cs.osprey.tools.Log.log;

/**
 * Train or load per-task MLP surrogates for one-body and pairwise energies.
 */
public class ConfSpaceSpecificMLPSurrogateCalculator {

    public final ConfEnergyCalculator confEcalc;
    public final MLPEnergyModel.TrainingConfig trainingConfig = new MLPEnergyModel.TrainingConfig();

    private int numSamplesPerParam = 8;
    private int minSamplesPerModel = 64;
    private int maxSamplesPerModel = 4096;
    private File cacheFile = null;

    public ConfSpaceSpecificMLPSurrogateCalculator(ConfEnergyCalculator confEcalc) {
        this.confEcalc = confEcalc;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setNumSamplesPerParam(int val) {
        numSamplesPerParam = val;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setMinSamplesPerModel(int val) {
        minSamplesPerModel = val;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setMaxSamplesPerModel(int val) {
        maxSamplesPerModel = val;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setCacheFile(File val) {
        cacheFile = val;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setHiddenSizes(int h1, int h2) {
        trainingConfig.hidden1 = h1;
        trainingConfig.hidden2 = h2;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setEpochs(int epochs) {
        trainingConfig.epochs = epochs;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setBatchSize(int batchSize) {
        trainingConfig.batchSize = batchSize;
        return this;
    }

    public ConfSpaceSpecificMLPSurrogateCalculator setLearningRate(double lr) {
        trainingConfig.learningRate = lr;
        return this;
    }

    public MLPSurrogateMatrix calc() {

        MLPSurrogateMatrix surrogate = new MLPSurrogateMatrix(confEcalc.confSpace);

        if (cacheFile != null && cacheFile.exists()) {
            surrogate.readFrom(cacheFile);
            log("read MLP surrogate matrix from file: %s", cacheFile.getAbsolutePath());
            return surrogate;
        }

        int numSingles = confEcalc.confSpace.getNumResConfs();
        int numPairs = countPairs(confEcalc.confSpace);
        Progress progress = new Progress(numSingles + numPairs);
        log("training %d MLP surrogate models (%d singles + %d pairs) ...", progress.getTotalWork(), numSingles, numPairs);

        for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

                confEcalc.tasks.submit(
                        () -> calcOneBody(pos1, rc1),
                        (model) -> {
                            surrogate.setOneBody(pos1.index, rc1.index, model);
                            progress.incrementProgress();
                        }
                );
            }
        }

        for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                for (SimpleConfSpace.Position pos2 : confEcalc.confSpace.positions.subList(0, pos1.index)) {
                    for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {

                        confEcalc.tasks.submit(
                                () -> calcPair(pos1, rc1, pos2, rc2),
                                (model) -> {
                                    surrogate.setPair(pos1.index, rc1.index, pos2.index, rc2.index, model);
                                    progress.incrementProgress();
                                }
                        );
                    }
                }
            }
        }

        confEcalc.tasks.waitForFinish();

        if (cacheFile != null) {
            surrogate.writeTo(cacheFile);
            log("wrote MLP surrogate matrix to file: %s", cacheFile.getAbsolutePath());
        }

        return surrogate;
    }

    private MLPEnergyModel calcOneBody(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc) {
        RCTuple tuple = new RCTuple(pos.index, rc.index);
        ResidueInteractions inters = new ResidueInteractions();
        inters.addPair(pos.resNum, pos.resNum);
        for (String shellResNum : confEcalc.confSpace.shellResNumbers) {
            inters.addPair(pos.resNum, shellResNum);
        }
        return calc(tuple, inters);
    }

    private MLPEnergyModel calcPair(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1,
                                    SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2) {
        RCTuple tuple = orderedPairTuple(pos1.index, rc1.index, pos2.index, rc2.index);
        String resNum1 = confEcalc.confSpace.positions.get(tuple.pos.get(0)).resNum;
        String resNum2 = confEcalc.confSpace.positions.get(tuple.pos.get(1)).resNum;
        ResidueInteractions inters = new ResidueInteractions();
        inters.addPair(resNum1, resNum2);
        return calc(tuple, inters);
    }

    private MLPEnergyModel calc(RCTuple tuple, ResidueInteractions inters) {

        ParametricMolecule pmol = confEcalc.confSpace.makeMolecule(tuple);
        try (EnergyFunction ff = confEcalc.ecalc.makeEnergyFunction(pmol, inters)) {

            MoleculeObjectiveFunction f = new MoleculeObjectiveFunction(pmol, ff);
            int numDofs = pmol.dofs.size();
            double[] mins = new double[numDofs];
            double[] maxs = new double[numDofs];
            for (int d = 0; d < numDofs; d++) {
                mins[d] = pmol.dofBounds.getMin(d);
                maxs[d] = pmol.dofBounds.getMax(d);
            }

            if (numDofs == 0) {
                double e = ff.getEnergy();
                return MLPEnergyModel.train(
                        new double[][]{new double[0]},
                        new double[]{e},
                        mins,
                        maxs,
                        trainingConfig,
                        tuple.hashCode()
                );
            }

            int paramCount = MLPEnergyModel.numParams(numDofs, trainingConfig.hidden1, trainingConfig.hidden2);
            int numSamples = 1 + numSamplesPerParam * paramCount;
            numSamples = Math.max(minSamplesPerModel, numSamples);
            numSamples = Math.min(maxSamplesPerModel, numSamples);

            double[][] x = new double[numSamples][numDofs];
            double[] y = new double[numSamples];

            // Sample 0: minimized center point to anchor the local minimum region
            SimpleCCDMinimizer centerMin = new SimpleCCDMinimizer(f);
            DoubleMatrix1D center = centerMin.minimizeFromCenter().dofValues;
            for (int d = 0; d < numDofs; d++) {
                x[0][d] = center.get(d);
            }
            y[0] = f.getValue(center);
            centerMin.close();

            Random rand = new Random(tuple.hashCode());
            for (int i = 1; i < numSamples; i++) {
                DoubleMatrix1D xi = DoubleFactory1D.dense.make(numDofs);
                for (int d = 0; d < numDofs; d++) {
                    double val = mins[d] + rand.nextDouble() * (maxs[d] - mins[d]);
                    x[i][d] = val;
                    xi.set(d, val);
                }
                y[i] = f.getValue(xi);
            }

            return MLPEnergyModel.train(x, y, mins, maxs, trainingConfig, tuple.hashCode());
        }
    }

    private static int countPairs(SimpleConfSpace confSpace) {
        int pairs = 0;
        for (SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
                    pairs += pos2.resConfs.size();
                }
            }
        }
        return pairs;
    }

    private static RCTuple orderedPairTuple(int pos1, int rc1, int pos2, int rc2) {
        if (pos1 < pos2) {
            return new RCTuple(pos1, rc1, pos2, rc2);
        }
        return new RCTuple(pos2, rc2, pos1, rc1);
    }
}

