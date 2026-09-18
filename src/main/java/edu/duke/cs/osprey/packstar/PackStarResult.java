package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import java.util.Objects;

/**
 * PACK* result retained by ordinary K* result caches and score writers.
 * The immutable functional estimate remains available after the pfunc closes.
 * No event configuration is required: disabled events report NOT_CONFIGURED.
 */
public final class PackStarResult extends PartitionFunction.Result {

    private final PackStarFunctionalObservableResult functionalObservable;

    PackStarResult(PartitionFunction.Result source,
                   PackStarFunctionalObservableResult functionalObservable) {
        super(source.status, source.values, source.numConfs);
        source.getStats().forEach(this::setStat);
        this.functionalObservable = Objects.requireNonNull(functionalObservable);
    }

    public PackStarFunctionalObservableResult getFunctionalObservableResult() {
        return functionalObservable;
    }
}
