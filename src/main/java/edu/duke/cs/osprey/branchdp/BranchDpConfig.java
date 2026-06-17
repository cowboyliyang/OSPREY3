package edu.duke.cs.osprey.branchdp;

/**
 * Shared branch-DP system-property resolver.
 *
 * <p>The branch-DP primitives own the legacy {@code branchmarkstar.*} keys.
 * PACK* callers may enter a scoped alias mode where matching {@code packstar.*}
 * keys are preferred, while plain BranchMARK* callers continue to see only the
 * legacy keys.</p>
 */
public final class BranchDpConfig {

    private static final ThreadLocal<Integer> PACKSTAR_ALIAS_DEPTH =
            ThreadLocal.withInitial(() -> 0);

    private BranchDpConfig() {
    }

    public static boolean isPackStarAliasScopeActive() {
        return PACKSTAR_ALIAS_DEPTH.get() > 0;
    }

    public static String getBackendName() {
        return isPackStarAliasScopeActive() ? "PACK*" : "BranchMARK*";
    }

    public static String getBackendLogPrefix() {
        return getBackendName() + ":";
    }

    public static String getBackendThreadNamePrefix() {
        return isPackStarAliasScopeActive() ? "packstar" : "branchmarkstar";
    }

    public static Scope enterPackStarAliasScope() {
        PACKSTAR_ALIAS_DEPTH.set(PACKSTAR_ALIAS_DEPTH.get() + 1);
        return new Scope();
    }

    public static String getBackendProperty(String key, String defaultValue) {
        return getProperty(key, defaultValue, isPackStarAliasScopeActive());
    }

    public static boolean getBackendBoolean(String key, boolean defaultValue) {
        String value = getBackendProperty(key, null);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    public static int getBackendInteger(String key, int defaultValue, String warningPrefix) {
        String value = getBackendProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value.trim());
        } catch (NumberFormatException e) {
            warn(warningPrefix, "Invalid integer", key, value, Integer.toString(defaultValue));
            return defaultValue;
        }
    }

    public static long getBackendLong(String key, long defaultValue, String warningPrefix) {
        String value = getBackendProperty(key, null);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Long.parseLong(value.trim());
        } catch (NumberFormatException e) {
            warn(warningPrefix, "Invalid long", key, value, Long.toString(defaultValue));
            return defaultValue;
        }
    }

    static String getProperty(String key, String defaultValue, boolean preferPackStarAliases) {
        if (!preferPackStarAliases) {
            String value = directProperty(key);
            return value != null ? value : defaultValue;
        }
        String preferredKey = preferPackStarKey(key);
        String value = directProperty(preferredKey);
        if (value == null) {
            String legacyKey = legacyBranchMarkStarKey(preferredKey);
            if (!legacyKey.equals(preferredKey)) {
                value = directProperty(legacyKey);
            }
        }
        return value != null ? value : defaultValue;
    }

    private static String directProperty(String key) {
        String value = System.getProperty(key);
        return value != null ? value : System.getProperty("osprey." + key);
    }

    private static String preferPackStarKey(String key) {
        String prefix = "branchmarkstar.";
        if (key.startsWith(prefix)) {
            return "packstar." + key.substring(prefix.length());
        }
        return key;
    }

    private static String legacyBranchMarkStarKey(String key) {
        String prefix = "packstar.";
        if (key.startsWith(prefix)) {
            return "branchmarkstar." + key.substring(prefix.length());
        }
        return key;
    }

    private static void warn(String warningPrefix, String kind, String key, String value, String defaultValue) {
        System.err.println(warningPrefix + " " + kind + " for '" + key
                + "': '" + value + "', using " + defaultValue + ".");
    }

    public static final class Scope implements AutoCloseable {

        private boolean closed = false;

        private Scope() {
        }

        @Override
        public void close() {
            if (closed) {
                return;
            }
            int depth = PACKSTAR_ALIAS_DEPTH.get();
            if (depth <= 1) {
                PACKSTAR_ALIAS_DEPTH.remove();
            } else {
                PACKSTAR_ALIAS_DEPTH.set(depth - 1);
            }
            closed = true;
        }
    }
}
