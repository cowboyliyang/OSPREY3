package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.energy.ConfEnergyCalculator;

import java.io.File;
import java.util.Locale;

/**
 * Load or train per-task MLP surrogate matrices.
 */
public class ConfSpaceSpecificMLPSurrogateCache {

    public static class Config {
        public File cacheRoot;
        public String taskTag = "default";
        public boolean forceRetrain = false;

        public int numSamplesPerParam = 8;
        public int minSamplesPerModel = 64;
        public int maxSamplesPerModel = 4096;
        public int hidden1 = 32;
        public int hidden2 = 32;
        public int epochs = 300;
        public int batchSize = 64;
        public double learningRate = 1e-3;
        public double l2 = 1e-8;
    }

    public static MLPSurrogateMatrix loadOrTrain(ConfEnergyCalculator confEcalc, Config config) {

        if (config == null) {
            throw new IllegalArgumentException("config must not be null");
        }
        if (config.cacheRoot == null) {
            throw new IllegalArgumentException("config.cacheRoot must not be null");
        }

        String confHash = ConfSpaceSpecificApproximatorCache.fingerprint(confEcalc.confSpace);
        File taskDir = new File(config.cacheRoot, sanitize(config.taskTag));
        File cacheFile = new File(taskDir, String.format(Locale.ROOT, "%s.mlp.matrix", confHash));

        if (!taskDir.exists() && !taskDir.mkdirs()) {
            throw new RuntimeException("Failed to create MLP surrogate cache directory: " + taskDir.getAbsolutePath());
        }

        if (config.forceRetrain && cacheFile.exists() && !cacheFile.delete()) {
            throw new RuntimeException("Failed to delete existing MLP surrogate cache: " + cacheFile.getAbsolutePath());
        }

        ConfSpaceSpecificMLPSurrogateCalculator calc = new ConfSpaceSpecificMLPSurrogateCalculator(confEcalc)
                .setNumSamplesPerParam(config.numSamplesPerParam)
                .setMinSamplesPerModel(config.minSamplesPerModel)
                .setMaxSamplesPerModel(config.maxSamplesPerModel)
                .setHiddenSizes(config.hidden1, config.hidden2)
                .setEpochs(config.epochs)
                .setBatchSize(config.batchSize)
                .setLearningRate(config.learningRate)
                .setCacheFile(cacheFile);
        calc.trainingConfig.l2 = config.l2;

        return calc.calc();
    }

    private static String sanitize(String s) {
        if (s == null || s.trim().isEmpty()) {
            return "default";
        }
        return s.replaceAll("[^a-zA-Z0-9._-]", "_");
    }
}

