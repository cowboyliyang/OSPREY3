package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.markstar.bench.KStarBaselineBench;

/**
 * Main entry point for K* baseline benchmark (no JUnit).
 * Usage: java -cp ... edu.duke.cs.osprey.markstar.RunKStarBaseline
 */
public class RunKStarBaseline {
    public static void main(String[] args) {
        new KStarBaselineBench().benchmarkKStar();
    }
}
