package edu.duke.cs.osprey.markstar;

/**
 * Simple main() entry point to run benchmarkGNNStrategies without JUnit.
 * Usage: java -cp ... edu.duke.cs.osprey.markstar.RunBenchmark
 * System properties are read inside benchmarkGNNStrategies().
 */
public class RunBenchmark {
    public static void main(String[] args) throws Exception {
        TestBranchMARKStar test = new TestBranchMARKStar();
        test.benchmarkGNNStrategies();
    }
}
