package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

/**
 * Internal PACK* backend contract.
 *
 * <p>This keeps the public PACK* facade independent from the legacy
 * BranchMARK* entry point while PACK* shares the neutral branch-DP engine.</p>
 */
public interface PackStarBackend extends PartitionFunction.WithConfDB {

    void setCorrections(UpdatingEnergyMatrix corrections);

    void setMaxNumConfs(int maxNumConfs);

    void setLeafMinimizationBatchSize(int batchSize);

    void useFullParallelLeafBatch();

    void setCorrectionTighteningEnabled(boolean enabled);

    void setReduceMinimizations(boolean enabled);
}
