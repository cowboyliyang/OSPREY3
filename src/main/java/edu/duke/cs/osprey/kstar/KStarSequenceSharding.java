package edu.duke.cs.osprey.kstar;

/**
 * Canonical sequence-bundle ownership for process-level K-star sharding.
 *
 * <p>The wild-type bundle is replicated because every rank needs the same
 * local stability baseline. Every other bundle has exactly one owner. This
 * class is shared by formal execution and PACK* preflight so their manifests
 * and workload predictions describe the work that is actually run.</p>
 */
public final class KStarSequenceSharding {

    public static final int NO_REPLICATED_BUNDLE = -1;
    public static final int WILD_TYPE_BUNDLE_ORDINAL = 0;

    private KStarSequenceSharding() {
    }

    public static int replicatedBundleOrdinal(boolean hasWildType) {
        return hasWildType
                ? WILD_TYPE_BUNDLE_ORDINAL : NO_REPLICATED_BUNDLE;
    }

    public static boolean isAssignedToShard(int bundleOrdinal,
                                             boolean hasWildType,
                                             int shardIndex,
                                             int shardCount) {
        validate(bundleOrdinal, shardIndex, shardCount);
        if (shardCount == 1) {
            return true;
        }
        if (hasWildType && bundleOrdinal == WILD_TYPE_BUNDLE_ORDINAL) {
            return true;
        }
        return bundleOrdinal % shardCount == shardIndex;
    }

    private static void validate(int bundleOrdinal, int shardIndex,
                                 int shardCount) {
        if (bundleOrdinal < 0) {
            throw new IllegalArgumentException(
                    "bundle ordinal must be non-negative: " + bundleOrdinal);
        }
        if (shardCount < 1 || shardIndex < 0 || shardIndex >= shardCount) {
            throw new IllegalArgumentException(
                    "invalid K* sequence shard index=" + shardIndex
                            + " for shardCount=" + shardCount);
        }
    }
}
