package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.math.BigDecimal;

/**
 * PACK* partition-function entry point.
 *
 * <p>This is currently a thin facade over the branch-decomposition PAC path.
 * The facade gives callers a PACK*-named API and keeps PAC mode mandatory for
 * this path.</p>
 */
public class PackStarPartitionFunction implements PartitionFunction.WithConfDB, PackStarSampleTraceable, AutoCloseable {

    private final PackStarBound delegate;

    public PackStarPartitionFunction(SimpleConfSpace confSpace,
                                     EnergyMatrix rigidEmat,
                                     EnergyMatrix minimizingEmat,
                                     ConfEnergyCalculator minimizingConfEcalc,
                                     RCs rcs,
                                     Parallelism parallelism) {
        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism, null);
    }

    public PackStarPartitionFunction(SimpleConfSpace confSpace,
                                     EnergyMatrix rigidEmat,
                                     EnergyMatrix minimizingEmat,
                                     ConfEnergyCalculator minimizingConfEcalc,
                                     RCs rcs,
                                     Parallelism parallelism,
                                     String stateNameOverride) {
        this.delegate = makeDelegate(confSpace, rigidEmat, minimizingEmat,
                minimizingConfEcalc, rcs, parallelism, stateNameOverride);
    }

    private static PackStarBound makeDelegate(SimpleConfSpace confSpace,
                                              EnergyMatrix rigidEmat,
                                              EnergyMatrix minimizingEmat,
                                              ConfEnergyCalculator minimizingConfEcalc,
                                              RCs rcs,
                                              Parallelism parallelism,
                                              String stateNameOverride) {
        return new PackStarBound(confSpace, rigidEmat, minimizingEmat,
                minimizingConfEcalc, rcs, parallelism, stateNameOverride);
    }

    public PackStarBound getDelegate() {
        return delegate;
    }

    public BranchDpAdmission.Prediction getAdmissionPrediction() {
        return delegate.getAdmissionPrediction();
    }

    @Override
    public void setInstanceId(int val) {
        delegate.setInstanceId(val);
    }

    public void setCorrections(UpdatingEnergyMatrix corrections) {
        delegate.setCorrections(corrections);
    }

    public void setMaxNumConfs(int maxNumConfs) {
        delegate.setMaxNumConfs(maxNumConfs);
    }

    public void setLeafMinimizationBatchSize(int batchSize) {
        delegate.setLeafMinimizationBatchSize(batchSize);
    }

    public void useFullParallelLeafBatch() {
        delegate.useFullParallelLeafBatch();
    }

    public void setCorrectionTighteningEnabled(boolean enabled) {
        delegate.setCorrectionTighteningEnabled(enabled);
    }

    public void setReduceMinimizations(boolean enabled) {
        delegate.setReduceMinimizations(enabled);
    }

    @Override
    public void setSampleListener(PackStarSampleListener listener) {
        delegate.setSampleListener(listener);
    }

    @Override
    public void setReportProgress(boolean val) {
        delegate.setReportProgress(val);
    }

    @Override
    public void setConfListener(ConfListener val) {
        delegate.setConfListener(val);
    }

    @Override
    public void init(double targetEpsilon) {
        delegate.init(targetEpsilon);
    }

    @Override
    public void setStabilityThreshold(BigDecimal stabilityThreshold) {
        delegate.setStabilityThreshold(stabilityThreshold);
    }

    @Override
    public Status getStatus() {
        return delegate.getStatus();
    }

    @Override
    public Values getValues() {
        return delegate.getValues();
    }

    @Override
    public int getParallelism() {
        return delegate.getParallelism();
    }

    @Override
    public int getNumConfsEvaluated() {
        return delegate.getNumConfsEvaluated();
    }

    @Override
    public void compute(int maxNumConfs) {
        delegate.compute(maxNumConfs);
    }

    @Override
    public Result makeResult() {
        return delegate.makeResult();
    }

    @Override
    public void setConfDB(ConfDB confDB, ConfDB.Key key) {
        delegate.setConfDB(confDB, key);
    }

    @Override
    public void close() {
        delegate.close();
    }
}
