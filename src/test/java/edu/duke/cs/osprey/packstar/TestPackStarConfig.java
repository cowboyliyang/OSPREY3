package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDpConfig;
import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestPackStarConfig {

    @Test
    public void packStarAliasWinsWhenEnabled() {
        withProperties(Map.of(
                "packstar.example.int", "7",
                "branchmarkstar.example.int", "3"
        ), () -> assertEquals(7,
                PackStarConfig.getInteger("branchmarkstar.example.int", 0, "[test]", true)));
    }

    @Test
    public void legacyFallbackIsAcceptedWhenAliasEnabled() {
        withProperties(Map.of(
                "branchmarkstar.example.int", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getInteger("packstar.example.int", 0, "[test]", true)));
    }

    @Test
    public void branchMarkStarModeIgnoresPackStarAlias() {
        withProperties(Map.of(
                "packstar.example.int", "7",
                "branchmarkstar.example.int", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getInteger("branchmarkstar.example.int", 0, "[test]", false)));
    }

    @Test
    public void branchMarkStarModeUsesDefaultWithoutLegacyKey() {
        withProperties(Map.of(
                "packstar.example.int", "7"
        ), () -> assertEquals(0,
                PackStarConfig.getInteger("branchmarkstar.example.int", 0, "[test]", false)));
    }

    @Test
    public void backendConfigUsesPackStarAliasInsidePackStarScope() {
        withProperties(Map.of(
                "packstar.example.backend", "7",
                "branchmarkstar.example.backend", "3"
        ), () -> {
            try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
                assertEquals(7,
                        PackStarConfig.getBackendInteger("branchmarkstar.example.backend", 0, "[test]"));
            }
        });
    }

    @Test
    public void backendConfigIgnoresPackStarAliasOutsidePackStarScope() {
        withProperties(Map.of(
                "packstar.example.backend", "7",
                "branchmarkstar.example.backend", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getBackendInteger("branchmarkstar.example.backend", 0, "[test]")));
    }

    @Test
    public void backendIdentityTracksPackStarScope() {
        assertEquals("BranchMARK*", BranchDpConfig.getBackendName());
        assertEquals("BranchMARK*:", BranchDpConfig.getBackendLogPrefix());
        assertEquals("branchmarkstar", BranchDpConfig.getBackendThreadNamePrefix());

        try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
            assertEquals("PACK*", BranchDpConfig.getBackendName());
            assertEquals("PACK*:", BranchDpConfig.getBackendLogPrefix());
            assertEquals("packstar", BranchDpConfig.getBackendThreadNamePrefix());
        }

        assertEquals("BranchMARK*", BranchDpConfig.getBackendName());
        assertEquals("BranchMARK*:", BranchDpConfig.getBackendLogPrefix());
        assertEquals("branchmarkstar", BranchDpConfig.getBackendThreadNamePrefix());
    }

    private static void withProperties(Map<String, String> values, Runnable test) {
        Map<String, String> oldValues = new LinkedHashMap<>();
        for (String key : values.keySet()) {
            oldValues.put(key, System.getProperty(key));
            System.clearProperty(key);
            System.clearProperty("osprey." + key);
        }
        try {
            for (Map.Entry<String, String> entry : values.entrySet()) {
                System.setProperty(entry.getKey(), entry.getValue());
            }
            test.run();
        } finally {
            for (String key : values.keySet()) {
                restoreProperty(key, oldValues.get(key));
                System.clearProperty("osprey." + key);
            }
        }
    }

    private static void restoreProperty(String key, String value) {
        if (value == null) {
            System.clearProperty(key);
        } else {
            System.setProperty(key, value);
        }
    }
}
