package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;

/**
 * Loads or trains a confspace-specific approximator matrix and stores it in a deterministic cache path.
 *
 * The cache key is a fingerprint of the task's confspace, so each OSPREY task
 * gets its own approximation model.
 */
public class ConfSpaceSpecificApproximatorCache {

    public static class Config {
        public File cacheRoot;
        public String taskTag = "default";
        public int numSamplesPerParam = 10;
        public ApproximatorMatrixCalculator.ApproximatorType type = ApproximatorMatrixCalculator.ApproximatorType.Quadratic;
        public boolean forceRetrain = false;
    }

    public static ApproximatorMatrix loadOrTrain(ConfEnergyCalculator confEcalc, Config config) {

        if (config == null) {
            throw new IllegalArgumentException("config must not be null");
        }
        if (config.cacheRoot == null) {
            throw new IllegalArgumentException("config.cacheRoot must not be null");
        }

        String confHash = fingerprint(confEcalc.confSpace);
        File taskDir = new File(config.cacheRoot, sanitize(config.taskTag));
        File cacheFile = new File(taskDir, String.format("%s.%s.amat", confHash, config.type.name().toLowerCase(Locale.ROOT)));

        if (!taskDir.exists() && !taskDir.mkdirs()) {
            throw new RuntimeException("Failed to create approximator cache directory: " + taskDir.getAbsolutePath());
        }

        if (config.forceRetrain && cacheFile.exists() && !cacheFile.delete()) {
            throw new RuntimeException("Failed to delete existing approximator cache: " + cacheFile.getAbsolutePath());
        }

        return new ApproximatorMatrixCalculator(confEcalc)
                .setApproximatorType(config.type)
                .setNumSamplesPerParam(config.numSamplesPerParam)
                .setCacheFile(cacheFile)
                .calc();
    }

    public static String fingerprint(SimpleConfSpace confSpace) {
        StringBuilder sb = new StringBuilder();

        sb.append("numPositions=").append(confSpace.positions.size()).append('\n');
        for (SimpleConfSpace.Position pos : confSpace.positions) {
            sb.append("pos=").append(pos.index)
                    .append("|mindex=").append(pos.mindex)
                    .append("|resNum=").append(pos.resNum)
                    .append("|numRC=").append(pos.resConfs.size())
                    .append('\n');

            for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
                sb.append("rc=").append(rc.index)
                        .append("|template=").append(rc.template.name)
                        .append("|type=").append(rc.type)
                        .append("|rotamer=").append(rc.rotamerIndex == null ? "null" : rc.rotamerIndex)
                        .append('\n');

                List<String> dofNames = new ArrayList<>(rc.dofBounds.keySet());
                Collections.sort(dofNames);
                for (String dofName : dofNames) {
                    double[] bounds = rc.dofBounds.get(dofName);
                    sb.append("dof=").append(dofName)
                            .append("|min=").append(bounds[0])
                            .append("|max=").append(bounds[1])
                            .append('\n');
                }
            }
        }

        List<String> shellNums = new ArrayList<>(confSpace.shellResNumbers);
        Collections.sort(shellNums);
        sb.append("shell=").append(shellNums).append('\n');

        return sha256Hex(sb.toString());
    }

    private static String sanitize(String s) {
        if (s == null || s.trim().isEmpty()) {
            return "default";
        }
        return s.replaceAll("[^a-zA-Z0-9._-]", "_");
    }

    private static String sha256Hex(String s) {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException ex) {
            throw new RuntimeException("SHA-256 unavailable", ex);
        }
        byte[] bytes = digest.digest(s.getBytes(StandardCharsets.UTF_8));
        StringBuilder hex = new StringBuilder(bytes.length * 2);
        for (byte b : bytes) {
            hex.append(String.format("%02x", b));
        }
        return hex.toString();
    }
}
