package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDpConfig;

/**
 * Thread-local marker for calls that are executing through the PACK* backend.
 *
 * <p>The PACK* entry point owns when branch-DP primitives may prefer
 * {@code packstar.*} aliases over neutral {@code branchdp.*} keys.</p>
 */
final class PackStarBackendRuntime {

    private PackStarBackendRuntime() {
    }

    static boolean isPackStarBackendActive() {
        return BranchDpConfig.isPackStarAliasScopeActive();
    }

    static Scope enter() {
        return new Scope(BranchDpConfig.enterPackStarAliasScope());
    }

    static final class Scope implements AutoCloseable {

        private final BranchDpConfig.Scope branchDpScope;
        private boolean closed = false;

        private Scope(BranchDpConfig.Scope branchDpScope) {
            this.branchDpScope = branchDpScope;
        }

        @Override
        public void close() {
            if (closed) {
                return;
            }
            branchDpScope.close();
            closed = true;
        }
    }
}
