package lphy.core.functions;

import lphy.graphicalModel.*;

import java.util.Arrays;
import java.util.List;

import static lphy.core.ParameterNames.NoParamName0;
import static lphy.core.ParameterNames.NoParamName1;

/**
 * @author Walter Xie
 */
public class Intersect<T> extends DeterministicFunction<T[]> {

    public Intersect(@ParameterInfo(name = NoParamName0, description = "set 1.") Value<T[]> a,
                     @ParameterInfo(name = NoParamName1, description = "set 2.") Value<T[]> b) {
        setInput(NoParamName0, a);
        setInput(NoParamName1, b);
    }

    @Override
    @GeneratorInfo(name = "intersect", verbClause = "is", narrativeName = "intersection",
            description = "A function to get intersection between two sets.")
    public Value<T[]> apply() {
        Value<T[]> a = (Value<T[]>)paramMap.get(NoParamName0);
        Class<?> aTy = a.value().getClass().getComponentType();
        Value<T[]> b = (Value<T[]>)paramMap.get(NoParamName1);
        Class<?> bTy = b.value().getClass().getComponentType();

        if (!aTy.equals(bTy))
            throw new IllegalArgumentException("Must use the same type !");

//        Set<String> intersection = set1.stream()
//                .filter(set2::contains)
//                .collect(Collectors.toSet());
        List<T> intersect = Arrays.stream(a.value()).distinct()
                .filter(x -> Arrays.asList(b.value()).contains(x)).toList();

        System.out.println("Intersect a vector "+ a.getCanonicalId() + " (" + a.value().length +
                ") with " + b.getCanonicalId() + " (" + b.value().length +
                "), the result length = " + intersect.size());
        return ValueUtils.createValue(intersect, this);
    }

}
