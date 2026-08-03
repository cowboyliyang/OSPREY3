package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDpAdmission;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

/**
 * Internal PACK* backend contract.
 *
 * <p>This keeps the public PACK* facade independent while PACK* shares the
 * neutral branch-DP engine.</p>
 */
public interface PackStarBackend extends PartitionFunction.WithConfDB, PackStarSampleTraceable, AutoCloseable {

    void setCorrections(UpdatingEnergyMatrix corrections);

    void setMaxNumConfs(int maxNumConfs);

    void setLeafMinimizationBatchSize(int batchSize);

    void useFullParallelLeafBatch();

    void setCorrectionTighteningEnabled(boolean enabled);

    void setReduceMinimizations(boolean enabled);

    BranchDpAdmission.Prediction getAdmissionPrediction();

    /** Stable K-star state role used to derive the estimator random stream. */
    @Override
    void setInstanceId(int val);

    /** Release rooted DP/search storage after the scalar result has been copied. */
    @Override
    void close();
}
