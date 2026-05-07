package edu.duke.cs.osprey.markstar;

/**
 * Standalone entry point for BranchMARK* benchmarks (with and without GNN).
 *
 * Strategies:
 *   branch_ccd  — BranchMARK* with CCD-only minimization (no GNN)
 *   branch_gnn  — BranchMARK* + GNN decoupled pool (Strategy 7-style)
 *
 * System properties:
 *   osprey.gnn.benchStrategy   = branch_ccd | branch_gnn (required)
 *   osprey.gnn.confSpace       = medium | highrot (default: medium)
 *   osprey.gnn.numSeqs         = number of sequences (default: 20)
 *   osprey.gnn.gpuBatchSize    = GPU batch size for GNN pool (default: 1000)
 *   osprey.gnn.eval.proteinModelPath = path to protein ONNX model
 *   osprey.gnn.eval.complexModelPath = path to complex ONNX model
 *   osprey.gnn.trainingConfSpace     = all20 | highrot | medium (for RC mapping)
 *   osprey.gnn.numCPUs         = number of CPUs (default: 4)
 *   osprey.gnn.ematCache       = emat cache pattern
 *
 * Usage:
 *   java -Dosprey.gnn.benchStrategy=branch_gnn \
 *        -Dosprey.gnn.confSpace=highrot \
 *        -Dosprey.gnn.gpuBatchSize=1000 \
 *        -Dosprey.gnn.eval.proteinModelPath=gnn_data/.../protein/model/gnn_model.onnx \
 *        -Dosprey.gnn.eval.complexModelPath=gnn_data/.../complex/model/gnn_model.onnx \
 *        -cp ... edu.duke.cs.osprey.markstar.RunBranchBenchmark
 */
public class RunBranchBenchmark {
    public static void main(String[] args) throws Exception {
        String which = System.getProperty("osprey.gnn.benchStrategy", "branch_ccd");
        if (!which.startsWith("branch_")) {
            System.err.println("RunBranchBenchmark only supports branch_ccd and branch_gnn strategies.");
            System.err.println("Got: " + which);
            System.err.println("Falling back to branch_ccd");
            System.setProperty("osprey.gnn.benchStrategy", "branch_ccd");
        }
        TestBranchMARKStar test = new TestBranchMARKStar();
        test.benchmarkGNNStrategies();
    }
}
