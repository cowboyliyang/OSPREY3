import ai.onnxruntime.*;
import java.nio.LongBuffer;
import java.util.Collections;
import java.util.Arrays;

public class TestOnnxBoth {
    public static void main(String[] args) throws Exception {
        String[] models = {
            "/usr/xtmp/lz280/osprey_data/gnn_data/2RL0_all20_4pos/protein/model/gnn_model.onnx",
            "/usr/xtmp/lz280/osprey_data/gnn_data/2RL0_all20_4pos/complex/model/gnn_model.onnx"
        };
        String[] names = {"protein", "complex"};

        int[][] testConfs = {
            {15, 13, 66, 15},
            {15, 13, 66, 4},
            {15, 13, 66, 5},
            {86, 90, 89, 87},
            {33, 30, 16, 34},
            {72, 188, 72, 188}
        };

        OrtEnvironment env = OrtEnvironment.getEnvironment();

        for (int m = 0; m < models.length; m++) {
            System.out.println("=== " + names[m] + " ===");
            OrtSession.SessionOptions opts = new OrtSession.SessionOptions();
            opts.setOptimizationLevel(OrtSession.SessionOptions.OptLevel.ALL_OPT);
            OrtSession session = env.createSession(models[m], opts);

            int batchSize = testConfs.length;
            int numPos = 4;
            long[] flat = new long[batchSize * numPos];
            for (int i = 0; i < batchSize; i++)
                for (int j = 0; j < numPos; j++)
                    flat[i * numPos + j] = testConfs[i][j];

            OnnxTensor input = OnnxTensor.createTensor(env, LongBuffer.wrap(flat), new long[]{batchSize, numPos});
            OrtSession.Result result = session.run(Collections.singletonMap("confs", input));
            Object raw = result.get(0).getValue();

            float[] out;
            if (raw instanceof float[]) {
                out = (float[]) raw;
            } else if (raw instanceof float[][]) {
                float[][] out2d = (float[][]) raw;
                out = new float[batchSize];
                for (int i = 0; i < batchSize; i++) out[i] = out2d[i][0];
            } else {
                System.out.println("Unexpected type: " + raw.getClass()); continue;
            }

            for (int i = 0; i < batchSize; i++) {
                System.out.printf("  conf=%-25s residual=%.4f%n", Arrays.toString(testConfs[i]), out[i]);
            }

            input.close(); result.close(); session.close();
        }
    }
}
