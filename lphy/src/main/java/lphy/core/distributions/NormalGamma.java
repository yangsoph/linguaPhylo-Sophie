package lphy.core.distributions;

import lphy.graphicalModel.*;
import lphy.util.RandomUtils;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Map;
import java.util.TreeMap;

import static lphy.core.distributions.DistributionConstants.*;
import static lphy.graphicalModel.ValueUtils.doubleValue;

/**
 * Normal-gamma distribution
 */
public class NormalGamma implements GenerativeDistribution<Double[]> {

    private Value<Number> shape;
    private Value<Number> scale;
    private Value<Number> mean;
    private Value<Number> precision;

    GammaDistribution gammaDistribution;
    NormalDistribution normalDistribution;

    static final String precisionParamName = "precision";

    public NormalGamma(@ParameterInfo(name = shapeParamName, description = "the shape of the distribution.") Value<Number> shape,
                       @ParameterInfo(name = scaleParamName, description = "the scale of the distribution.") Value<Number> scale,
                       @ParameterInfo(name = meanParamName, description = "the mean of the distribution.") Value<Number> mean,
                       @ParameterInfo(name = precisionParamName, narrativeName = "precision", description = "the standard deviation of the distribution.") Value<Number> precision) {

        this.shape = shape;
        this.scale = scale;

        this.mean = mean;
        this.precision = precision;

        constructDistribution();
    }

    @GeneratorInfo(name = "NormalGamma", verbClause = "has", narrativeName = "normal-gamma prior",
            category = GeneratorCategory.PROB_DIST, examples = {"simplePhyloBrownian.lphy","simplePhyloOU.lphy"},
            description = "The normal-gamma probability distribution.")
    public RandomVariable<Double[]> sample() {
        double T = gammaDistribution.sample();
        double x = normalDistribution.sample();
        return new RandomVariable<>(null, new Double[] {x, T}, this);
    }

    private void constructDistribution() {
        if (mean == null) throw new IllegalArgumentException("The mean value can't be null!");
        if (precision == null) throw new IllegalArgumentException("The precision value can't be null!");

        double m = doubleValue(mean);
        double sh = doubleValue(shape);
        double sc = doubleValue(scale);
        double lambda = doubleValue(precision);

        RandomGenerator random = RandomUtils.getRandom();
        gammaDistribution = new GammaDistribution(random, sh, sc);

        double T = gammaDistribution.sample();
        normalDistribution = new NormalDistribution(random, m, lambda * T);
    }

    @Override
    public double density(Double[] x) {
        return gammaDistribution.density(x[0]) * normalDistribution.density(x[1]);
    }

    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(shapeParamName, shape);
            put(scaleParamName, scale);
            put(meanParamName, mean);
            put(precisionParamName, precision);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case meanParamName -> mean = value;
            case precisionParamName -> precision = value;
            case shapeParamName -> shape = value;
            case scaleParamName -> scale = value;
            default -> throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
        constructDistribution();
    }

    public String toString() {
        return getName();
    }

    public Value<Number> getMean() {
        return mean;
    }

    public Value<Number> getPrecision() {
        return precision;
    }
}