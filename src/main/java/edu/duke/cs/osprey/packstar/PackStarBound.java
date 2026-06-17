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

/**
 * PACK* branch-decomposition partition-function bound.
 *
 * <p>This class owns the PACK* API while reusing the legacy branch-DP backend
 * during migration.  PACK* estimator mode is forced for the backend; callers do
 * not need to set the old BranchMARK* PAC switch.</p>
 */
public class PackStarBound implements PartitionFunction.WithConfDB {

    private final PackStarBackend backend;

    public PackStarBound(SimpleConfSpace confSpace,
                         EnergyMatrix rigidEmat,
                         EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc,
                         RCs rcs,
                         Parallelism parallelism) {
        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism, null);
    }

    public PackStarBound(SimpleConfSpace confSpace,
                         EnergyMatrix rigidEmat,
                         EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc,
                         RCs rcs,
                         Parallelism parallelism,
                         String stateNameOverride) {
        this(new BranchPackStarBackend(confSpace, rigidEmat, minimizingEmat,
                minimizingConfEcalc, rcs, parallelism, stateNameOverride));
    }

    PackStarBound(PackStarBackend backend) {
        this.backend = backend;
    }

    public void setCorrections(UpdatingEnergyMatrix corrections) {
        backend.setCorrections(corrections);
    }

    public void setMaxNumConfs(int maxNumConfs) {
        backend.setMaxNumConfs(maxNumConfs);
    }

    public void setLeafMinimizationBatchSize(int batchSize) {
        backend.setLeafMinimizationBatchSize(batchSize);
    }

    public void useFullParallelLeafBatch() {
        backend.useFullParallelLeafBatch();
    }

    public void setCorrectionTighteningEnabled(boolean enabled) {
        backend.setCorrectionTighteningEnabled(enabled);
    }

    public void setReduceMinimizations(boolean enabled) {
        backend.setReduceMinimizations(enabled);
    }

    @Override
    public void setReportProgress(boolean val) {
        backend.setReportProgress(val);
    }

    @Override
    public void setConfListener(ConfListener val) {
        backend.setConfListener(val);
    }

    @Override
    public void init(double targetEpsilon) {
        backend.init(targetEpsilon);
    }

    @Override
    public void setStabilityThreshold(BigDecimal stabilityThreshold) {
        backend.setStabilityThreshold(stabilityThreshold);
    }

    @Override
    public Status getStatus() {
        return backend.getStatus();
    }

    @Override
    public Values getValues() {
        return backend.getValues();
    }

    @Override
    public int getParallelism() {
        return backend.getParallelism();
    }

    @Override
    public int getNumConfsEvaluated() {
        return backend.getNumConfsEvaluated();
    }

    @Override
    public void compute(int maxNumConfs) {
        backend.compute(maxNumConfs);
    }

    @Override
    public Result makeResult() {
        return backend.makeResult();
    }

    @Override
    public void setConfDB(ConfDB confDB, ConfDB.Key key) {
        backend.setConfDB(confDB, key);
    }
}
