package edu.duke.cs.osprey.packstar;

import java.util.function.IntFunction;

/** Shared final-sample event evaluation; failures never become negative hits. */
final class PackStarFunctionalObservableEvaluator {

    private PackStarFunctionalObservableEvaluator() {}

    static PackStarFunctionalObservableResult evaluate(
            String name, PackStarFunctionalEvent event, double rt,
            double[] logWeights,
            IntFunction<PackStarFunctionalEvent.Sample> samples) {
        PackStarFunctionalObservableAccumulator accumulator =
                new PackStarFunctionalObservableAccumulator(name, rt);
        int hits = 0;
        try {
            for (int i = 0; i < logWeights.length; i++) {
                boolean hit = event.test(samples.apply(i));
                if (hit) hits++;
                accumulator.add(hit, logWeights[i]);
            }
            return accumulator.finish();
        } catch (RuntimeException ex) {
            return PackStarFunctionalObservableResult.eventError(
                    name, logWeights.length, hits, rt,
                    "event evaluation failed: " + ex.getClass().getSimpleName()
                            + ": " + String.valueOf(ex.getMessage()).replace('\t', ' ')
                            .replace('\n', ' ').replace('\r', ' '));
        }
    }
}
