package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.math.BigDecimal;
import java.util.function.Supplier;

/**
 * Adapter that keeps PACK* branch-DP lifecycle calls inside the PACK* runtime scope.
 *
 * <p>{@link PackStarBranchDpBackend} owns PACK* backend policy; this wrapper
 * ensures shared branch primitives see the correct runtime scope.</p>
 */
final class BranchPackStarBackend implements PackStarBackend {

    private final PackStarBranchDpBackend delegate;

    BranchPackStarBackend(SimpleConfSpace confSpace,
                          EnergyMatrix rigidEmat,
                          EnergyMatrix minimizingEmat,
                          ConfEnergyCalculator minimizingConfEcalc,
                          RCs rcs,
                          Parallelism parallelism,
                          String stateNameOverride) {
        try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
            this.delegate = new PackStarBranchDpBackend(confSpace, rigidEmat, minimizingEmat,
                    minimizingConfEcalc, rcs, parallelism, stateNameOverride);
        }
    }

    @Override
    public void setCorrections(UpdatingEnergyMatrix corrections) {
        runInPackStarBackend(() -> delegate.setCorrections(corrections));
    }

    @Override
    public void setMaxNumConfs(int maxNumConfs) {
        runInPackStarBackend(() -> delegate.setMaxNumConfs(maxNumConfs));
    }

    @Override
    public void setLeafMinimizationBatchSize(int batchSize) {
        runInPackStarBackend(() -> delegate.setLeafMinimizationBatchSize(batchSize));
    }

    @Override
    public void useFullParallelLeafBatch() {
        runInPackStarBackend(delegate::useFullParallelLeafBatch);
    }

    @Override
    public void setCorrectionTighteningEnabled(boolean enabled) {
        runInPackStarBackend(() -> delegate.setCorrectionTighteningEnabled(enabled));
    }

    @Override
    public void setReduceMinimizations(boolean enabled) {
        runInPackStarBackend(() -> delegate.setReduceMinimizations(enabled));
    }

    @Override
    public void setReportProgress(boolean val) {
        runInPackStarBackend(() -> delegate.setReportProgress(val));
    }

    @Override
    public void setConfListener(ConfListener val) {
        runInPackStarBackend(() -> delegate.setConfListener(val));
    }

    @Override
    public void init(double targetEpsilon) {
        runInPackStarBackend(() -> delegate.init(targetEpsilon));
    }

    @Override
    public void setStabilityThreshold(BigDecimal stabilityThreshold) {
        runInPackStarBackend(() -> delegate.setStabilityThreshold(stabilityThreshold));
    }

    @Override
    public Status getStatus() {
        return callInPackStarBackend(delegate::getStatus);
    }

    @Override
    public Values getValues() {
        return callInPackStarBackend(delegate::getValues);
    }

    @Override
    public int getParallelism() {
        return callInPackStarBackend(delegate::getParallelism);
    }

    @Override
    public int getNumConfsEvaluated() {
        return callInPackStarBackend(delegate::getNumConfsEvaluated);
    }

    @Override
    public void compute(int maxNumConfs) {
        runInPackStarBackend(() -> delegate.compute(maxNumConfs));
    }

    @Override
    public Result makeResult() {
        return callInPackStarBackend(delegate::makeResult);
    }

    @Override
    public void setConfDB(ConfDB confDB, ConfDB.Key key) {
        runInPackStarBackend(() -> delegate.setConfDB(confDB, key));
    }

    private void runInPackStarBackend(Runnable action) {
        try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
            action.run();
        }
    }

    private <T> T callInPackStarBackend(Supplier<T> supplier) {
        try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
            return supplier.get();
        }
    }
}
