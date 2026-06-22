package edu.duke.cs.osprey.packstar;

import java.util.Locale;

/**
 * Shared PACK* system-property resolver.
 *
 * <p>PACK* owns the {@code packstar.*} keys. When a shared branch-DP key is
 * passed in, the matching PACK* key is preferred and the neutral
 * {@code branchdp.*} key remains available as a fallback.</p>
 */
public final class PackStarConfig {

    private PackStarConfig() {
    }

    public static String getProperty(String key, String defaultValue) {
        return getProperty(key, defaultValue, true);
    }

    public static String getProperty(String key, String defaultValue, boolean preferPackStarAliases) {
        if (preferPackStarAliases) {
            String preferredKey = preferPackStarKey(key);
            String value = directProperty(preferredKey);
            if (value == null) {
                String neutralKey = neutralBranchDpKey(preferredKey);
                if (!neutralKey.equals(preferredKey)) {
                    value = directProperty(neutralKey);
                }
            }
            return value != null ? value : defaultValue;
        }

        String value = directProperty(key);
        return value != null ? value : defaultValue;
    }

    public static String getBackendProperty(String key, String defaultValue) {
        return getProperty(key, defaultValue, PackStarBackendRuntime.isPackStarBackendActive());
    }

    public static boolean getBoolean(String key, boolean defaultValue) {
        return getBoolean(key, defaultValue, true);
    }

    public static boolean getBoolean(String key, boolean defaultValue, boolean preferPackStarAliases) {
        String value = getProperty(key, null, preferPackStarAliases);
        return value == null ? defaultValue : Boolean.parseBoolean(value);
    }

    public static boolean getBackendBoolean(String key, boolean defaultValue) {
        return getBoolean(key, defaultValue, PackStarBackendRuntime.isPackStarBackendActive());
    }

    public static int getInteger(String key, int defaultValue, String warningPrefix) {
        return getInteger(key, defaultValue, warningPrefix, true);
    }

    public static int getInteger(String key, int defaultValue, String warningPrefix,
                                 boolean preferPackStarAliases) {
        String value = getProperty(key, null, preferPackStarAliases);
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

    public static int getBackendInteger(String key, int defaultValue, String warningPrefix) {
        return getInteger(key, defaultValue, warningPrefix,
                PackStarBackendRuntime.isPackStarBackendActive());
    }

    public static long getLong(String key, long defaultValue, String warningPrefix) {
        return getLong(key, defaultValue, warningPrefix, true);
    }

    public static long getLong(String key, long defaultValue, String warningPrefix,
                               boolean preferPackStarAliases) {
        String value = getProperty(key, null, preferPackStarAliases);
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

    public static long getBackendLong(String key, long defaultValue, String warningPrefix) {
        return getLong(key, defaultValue, warningPrefix,
                PackStarBackendRuntime.isPackStarBackendActive());
    }

    public static double getDouble(String key, double defaultValue, String warningPrefix) {
        return getDouble(key, defaultValue, warningPrefix, true);
    }

    public static double getDouble(String key, double defaultValue, String warningPrefix,
                                   boolean preferPackStarAliases) {
        String value = getProperty(key, null, preferPackStarAliases);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            warn(warningPrefix, "Invalid double", key, value, Double.toString(defaultValue));
            return defaultValue;
        }
    }

    public static long getBytes(String key, long defaultValue, String warningPrefix) {
        return getBytes(key, defaultValue, warningPrefix, true);
    }

    public static long getBytes(String key, long defaultValue, String warningPrefix,
                                boolean preferPackStarAliases) {
        String value = getProperty(key, null, preferPackStarAliases);
        if (value == null || value.trim().isEmpty()) {
            return defaultValue;
        }
        try {
            return parseByteCount(value.trim());
        } catch (NumberFormatException e) {
            warn(warningPrefix, "Invalid byte count", key, value, Long.toString(defaultValue));
            return defaultValue;
        }
    }

    public static long parseByteCount(String value) {
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace("_", "");
        long multiplier = 1L;
        if (normalized.endsWith("kib")) {
            multiplier = 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("kb") || normalized.endsWith("k")) {
            multiplier = 1024L;
            normalized = normalized.replaceAll("kb?$", "");
        } else if (normalized.endsWith("mib")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("mb") || normalized.endsWith("m")) {
            multiplier = 1024L * 1024L;
            normalized = normalized.replaceAll("mb?$", "");
        } else if (normalized.endsWith("gib")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("gb") || normalized.endsWith("g")) {
            multiplier = 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("gb?$", "");
        } else if (normalized.endsWith("tib")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.substring(0, normalized.length() - 3);
        } else if (normalized.endsWith("tb") || normalized.endsWith("t")) {
            multiplier = 1024L * 1024L * 1024L * 1024L;
            normalized = normalized.replaceAll("tb?$", "");
        } else if (normalized.endsWith("b")) {
            normalized = normalized.substring(0, normalized.length() - 1);
        }

        double amount = Double.parseDouble(normalized.trim());
        if (amount < 0 || amount > Long.MAX_VALUE / (double) multiplier) {
            throw new NumberFormatException(value);
        }
        return (long) (amount * multiplier);
    }

    private static String directProperty(String key) {
        String value = System.getProperty(key);
        return value != null ? value : System.getProperty("osprey." + key);
    }

    private static String preferPackStarKey(String key) {
        String prefix = "branchdp.";
        if (key.startsWith(prefix)) {
            return "packstar." + key.substring(prefix.length());
        }
        return key;
    }

    private static String neutralBranchDpKey(String key) {
        String prefix = "packstar.";
        if (key.startsWith(prefix)) {
            return "branchdp." + key.substring(prefix.length());
        }
        return key;
    }

    private static void warn(String warningPrefix, String kind, String key, String value, String defaultValue) {
        System.err.println(warningPrefix + " " + kind + " for '" + key
                + "': '" + value + "', using " + defaultValue + ".");
    }
}
