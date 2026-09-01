package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.branchdp.BranchDpConfig;
import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestPackStarConfig {

    @Test
    public void packStarPropertyWinsWhenEnabled() {
        withProperties(Map.of(
                "packstar.example.int", "7",
                "branchdp.example.int", "3"
        ), () -> assertEquals(7,
                PackStarConfig.getInteger("branchdp.example.int", 0, "[test]", true)));
    }

    @Test
    public void branchDpFallbackIsAcceptedWhenPackStarKeyIsMissing() {
        withProperties(Map.of(
                "branchdp.example.int", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getInteger("packstar.example.int", 0, "[test]", true)));
    }

    @Test
    public void directModeIgnoresPackStarAlias() {
        withProperties(Map.of(
                "packstar.example.int", "7",
                "branchdp.example.int", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getInteger("branchdp.example.int", 0, "[test]", false)));
    }

    @Test
    public void directModeUsesDefaultWithoutBranchDpKey() {
        withProperties(Map.of(
                "packstar.example.int", "7"
        ), () -> assertEquals(0,
                PackStarConfig.getInteger("branchdp.example.int", 0, "[test]", false)));
    }

    @Test
    public void backendConfigUsesPackStarAliasInsidePackStarScope() {
        withProperties(Map.of(
                "packstar.example.backend", "7",
                "branchdp.example.backend", "3"
        ), () -> {
            try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
                assertEquals(7,
                        PackStarConfig.getBackendInteger("branchdp.example.backend", 0, "[test]"));
            }
        });
    }

    @Test
    public void backendConfigIgnoresPackStarAliasOutsidePackStarScope() {
        withProperties(Map.of(
                "packstar.example.backend", "7",
                "branchdp.example.backend", "3"
        ), () -> assertEquals(3,
                PackStarConfig.getBackendInteger("branchdp.example.backend", 0, "[test]")));
    }

    @Test
    public void backendIdentityTracksPackStarScope() {
        assertEquals("Branch-DP", BranchDpConfig.getBackendName());
        assertEquals("Branch-DP:", BranchDpConfig.getBackendLogPrefix());
        assertEquals("branchdp", BranchDpConfig.getBackendThreadNamePrefix());

        try (PackStarBackendRuntime.Scope scope = PackStarBackendRuntime.enter()) {
            assertEquals("PACK*", BranchDpConfig.getBackendName());
            assertEquals("PACK*:", BranchDpConfig.getBackendLogPrefix());
            assertEquals("packstar", BranchDpConfig.getBackendThreadNamePrefix());
        }

        assertEquals("Branch-DP", BranchDpConfig.getBackendName());
        assertEquals("Branch-DP:", BranchDpConfig.getBackendLogPrefix());
        assertEquals("branchdp", BranchDpConfig.getBackendThreadNamePrefix());
    }

    @Test
    public void packStarAdmissionCountsWorstCaseDpSweeps() {
        assertEquals(6, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4, false));
        assertEquals(2, PackStarBranchDpBackend.conservativeDpSweeps(
                true, false, 4, false));
        assertEquals(2, PackStarBranchDpBackend.conservativeDpSweeps(
                false, true, 4, false));
        assertEquals(7, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4, true));
        assertEquals(12, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4, 6, false));
        assertEquals(13, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4, 6, true));
        assertEquals(13, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4, 7, false));
        assertEquals(16,
                PackStarBranchDpBackend
                        .conservativeEtaV4AdditionalDpSweeps(2, 4));
        assertEquals(22, PackStarBranchDpBackend.conservativeDpSweeps(
                true, true, 4,
                PackStarBranchDpBackend
                        .conservativeEtaV4AdditionalDpSweeps(2, 4),
                false));
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
