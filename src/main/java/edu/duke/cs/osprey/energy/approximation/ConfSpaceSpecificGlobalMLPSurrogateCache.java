package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.energy.ConfEnergyCalculator;

import java.io.File;
import java.util.Locale;

/**
 * Load or train one global MLP surrogate per conf space task.
 */
public class ConfSpaceSpecificGlobalMLPSurrogateCache {

    public static class Config {
        public File cacheRoot;
        public String taskTag = "default";
        public boolean forceRetrain = false;

        public int numTupleSamples = 256;
        public int samplesPerTuple = 32;
        public int hidden1 = 64;
        public int hidden2 = 64;
        public int epochs = 400;
        public int batchSize = 128;
        public double learningRate = 1e-3;
        public double l2 = 1e-8;
    }

    public static TaskGlobalMLPSurrogate loadOrTrain(ConfEnergyCalculator confEcalc, Config config) {

        if (config == null) {
            throw new IllegalArgumentException("config must not be null");
        }
        if (config.cacheRoot == null) {
            throw new IllegalArgumentException("config.cacheRoot must not be null");
        }

        String confHash = ConfSpaceSpecificApproximatorCache.fingerprint(confEcalc.confSpace);
        File taskDir = new File(config.cacheRoot, sanitize(config.taskTag));
        File cacheFile = new File(taskDir, String.format(Locale.ROOT, "%s.mlp.global", confHash));

        if (!taskDir.exists() && !taskDir.mkdirs()) {
            throw new RuntimeException("Failed to create global MLP surrogate cache directory: " + taskDir.getAbsolutePath());
        }

        if (config.forceRetrain && cacheFile.exists() && !cacheFile.delete()) {
            throw new RuntimeException("Failed to delete existing global MLP cache: " + cacheFile.getAbsolutePath());
        }

        ConfSpaceSpecificGlobalMLPSurrogateCalculator calc = new ConfSpaceSpecificGlobalMLPSurrogateCalculator(confEcalc)
                .setNumTupleSamples(config.numTupleSamples)
                .setSamplesPerTuple(config.samplesPerTuple)
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
