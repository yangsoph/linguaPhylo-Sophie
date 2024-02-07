package lphy.base.evolution.substitutionmodel;

import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

/**
 * Created by Alexei Drummond on 2/02/20.
 */
public class GeneralTimeReversible extends RateMatrix {

    protected final String ratesParamName = SubstModelParamNames.RatesParamName;
    protected final String freqParamName = SubstModelParamNames.FreqParamName;

    int numStates;
    int ratesDim;

    public GeneralTimeReversible(@ParameterInfo(name = ratesParamName, description = "the relative rates of the GTR process.") Value<Double[]> rates,
                                 @ParameterInfo(name = freqParamName, description = "the base frequencies.") Value<Double[]> freq,
                                 @ParameterInfo(name = RateMatrix.meanRateParamName, description = "the base frequencies.", optional = true) Value<Number> meanRate) {

        super(meanRate);

        setParam(ratesParamName, rates);
        setParam(freqParamName, freq);
        update(rates, freq);
    }

    // cannot have verbClause and narrativeName here, see tutorials/h5n1.lphy
    @GeneratorInfo(name = "generalTimeReversible",
            category = GeneratorCategory.RATE_MATRIX,
            examples = {"https://linguaphylo.github.io/tutorials/discrete-phylogeography/"},
            description = "The general time reversible instantaneous rate matrix. Takes relative rates and base frequencies and produces an general time reversible rate matrix.")
    public Value<Double[][]> apply() {
        Value<Double[]> rates = getRates();
        Value<Double[]> freq = getFreq();
        update (rates, freq);

        return new DoubleArray2DValue(generalTimeReversible(rates.value(), freq.value()), this);
    }

    // symmetric rate matrix: rates and indicators dimension = n(n-1)/2
    private void update(Value<Double[]> rates, Value<Double[]> freq) {
        numStates = freq.value().length;
        ratesDim = numStates * (numStates - 1) / 2;
        if (rates.value().length != ratesDim)
            throw new RuntimeException("Expected dimension of " + ratesDim + " for the rates of a " + numStates + " state model.");
    }

    private Double[][] generalTimeReversible(Double[] rates, Double[] freqs) {

        Double[][] Q = new Double[numStates][numStates];

        double[] totalRates = new double[numStates];

        // construct off-diagonals
        int upper = 0;
        for (int i = 0; i < numStates; i++) {
            for (int j = i + 1; j < numStates; j++) {
                Q[i][j] = rates[upper] * freqs[j];
                Q[j][i] = rates[upper] * freqs[i];
                upper += 1;
            }
        }

        // construct diagonals
        for (int i = 0; i < numStates; i++) {
            double totalRate = 0.0;
            for (int j = 0; j < numStates; j++) {
                if (j != i) {
                    totalRate += Q[i][j];
                }
            }
            Q[i][i] = -totalRate;
        }
        // normalise rate matrix to one expected substitution per unit time
        normalize(freqs, Q, totalRateDefault1());

        return Q;
    }

    public Value<Double[]> getRates() {
        return getParams().get(ratesParamName);
    }

    public Value<Double[]> getFreq() {
        return getParams().get(freqParamName);
    }
}
