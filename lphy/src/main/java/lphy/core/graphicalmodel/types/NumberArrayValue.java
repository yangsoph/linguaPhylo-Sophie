package lphy.core.graphicalmodel.types;

import lphy.core.graphicalmodel.components.DeterministicFunction;
import lphy.core.vectorization.VectorValue;

public class NumberArrayValue extends VectorValue<Number> {

    public NumberArrayValue(String id, Number[] value) {
        super(id, value);
    }

    public NumberArrayValue(String id, Number[] value, DeterministicFunction function) {
        super(id, value, function);
    }
}
