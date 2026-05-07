import ai.onnxruntime.*;
import java.nio.LongBuffer;
import java.util.Collections;
import java.util.Arrays;

public class TestOnnx {
    public static void main(String[] args) throws Exception {
        String modelPath = "/usr/xtmp/lz280/osprey_data/gnn_data/2RL0_all20_4pos/protein/model/gnn_model.onnx";
        
        OrtEnvironment env = OrtEnvironment.getEnvironment();
        OrtSession.SessionOptions opts = new OrtSession.SessionOptions();
        opts.setOptimizationLevel(OrtSession.SessionOptions.OptLevel.ALL_OPT);
        OrtSession session = env.createSession(modelPath, opts);
        
        System.out.println("Model loaded successfully");
        
        for (var entry : session.getInputInfo().entrySet()) {
            System.out.println("Input: " + entry.getKey() + " -> " + entry.getValue().getInfo());
        }
        for (var entry : session.getOutputInfo().entrySet()) {
            System.out.println("Output: " + entry.getKey() + " -> " + entry.getValue().getInfo());
        }
        
        int[][] testConfs = {
            {86, 90, 89, 87},
            {33, 30, 16, 34},
            {72, 188, 72, 188}
        };
        
        int batchSize = testConfs.length;
        int numPos = 4;
        
        // Method 1: LongBuffer (same as GNNConfEnergyCalculator)
        long[] flatConfs = new long[batchSize * numPos];
        for (int i = 0; i < batchSize; i++) {
            for (int j = 0; j < numPos; j++) {
                flatConfs[i * numPos + j] = testConfs[i][j];
            }
        }
        System.out.println("\nflat array: " + Arrays.toString(flatConfs));
        
        OnnxTensor inputTensor = OnnxTensor.createTensor(
            env,
            LongBuffer.wrap(flatConfs),
            new long[]{batchSize, numPos}
        );
        System.out.println("Input tensor shape: " + Arrays.toString(inputTensor.getInfo().getShape()));
        
        OrtSession.Result result = session.run(Collections.singletonMap("confs", inputTensor));
        Object rawOutput = result.get(0).getValue();
        System.out.println("Output type: " + rawOutput.getClass().getName());
        
        if (rawOutput instanceof float[][]) {
            float[][] out2d = (float[][]) rawOutput;
            for (int i = 0; i < batchSize; i++) {
                System.out.printf("conf=%s residual=%.4f (float[][])%n",
                    Arrays.toString(testConfs[i]), out2d[i][0]);
            }
        } else if (rawOutput instanceof float[]) {
            float[] out1d = (float[]) rawOutput;
            for (int i = 0; i < batchSize; i++) {
                System.out.printf("conf=%s residual=%.4f (float[])%n",
                    Arrays.toString(testConfs[i]), out1d[i]);
            }
        } else {
            System.out.println("Unexpected: " + rawOutput.getClass().getName() + " = " + rawOutput);
        }
        
        inputTensor.close();
        result.close();
        
        // Method 2: long[][] directly
        System.out.println("\n--- Method 2: long[][] ---");
        long[][] confs2d = new long[batchSize][numPos];
        for (int i = 0; i < batchSize; i++) {
            for (int j = 0; j < numPos; j++) {
                confs2d[i][j] = testConfs[i][j];
            }
        }
        OnnxTensor inputTensor2 = OnnxTensor.createTensor(env, confs2d);
        OrtSession.Result result2 = session.run(Collections.singletonMap("confs", inputTensor2));
        Object rawOutput2 = result2.get(0).getValue();
        
        if (rawOutput2 instanceof float[][]) {
            float[][] out2d = (float[][]) rawOutput2;
            for (int i = 0; i < batchSize; i++) {
                System.out.printf("conf=%s residual=%.4f%n", Arrays.toString(testConfs[i]), out2d[i][0]);
            }
        } else if (rawOutput2 instanceof float[]) {
            float[] out1d = (float[]) rawOutput2;
            for (int i = 0; i < batchSize; i++) {
                System.out.printf("conf=%s residual=%.4f%n", Arrays.toString(testConfs[i]), out1d[i]);
            }
        }
        
        inputTensor2.close();
        result2.close();
        session.close();
    }
}
