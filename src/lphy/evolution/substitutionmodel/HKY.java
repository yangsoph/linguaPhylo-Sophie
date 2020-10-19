package lphy.evolution.substitutionmodel;

import lphy.graphicalModel.*;
import lphy.graphicalModel.types.DoubleArray2DValue;

import java.util.Map;

import static lphy.graphicalModel.ValueUtils.doubleValue;

/**
 * Created by adru001 on 2/02/20.
 */
public class HKY extends RateMatrix {

    public static final String kappaParamName = "kappa";
    public static final String freqParamName =  "freq";
    public static final String rateParamName =  "rate";


    public HKY(@ParameterInfo(name = "kappa", description = "the kappa of the HKY process.") Value<Number> kappa,
               @ParameterInfo(name = "freq", description = "the base frequencies.") Value<Double[]> freq,
               @ParameterInfo(name = "rate", description = "the total rate of substitution per unit time. Default 1.0.", optional = true) Value<Number> rate) {
        setParam(kappaParamName, kappa);
        setParam(freqParamName, freq);
        if (rate != null) setParam(rateParamName, rate);
    }


    @GeneratorInfo(name = "hky", description = "The HKY instantaneous rate matrix. Takes a kappa and base frequencies (and optionally a total rate) and produces an HKY85 rate matrix.")
    public Value<Double[][]> apply() {

        Map<String, Value> params = getParams();
        double kappa = doubleValue((Value<Number>)params.get(kappaParamName));
        Double[] freq = ((Value<Double[]>)params.get(freqParamName)).value();
        Value<Number> rate = params.getOrDefault(rateParamName, Value.Double_1);

        return new DoubleArray2DValue(hky(kappa, freq, doubleValue(rate)), this);
    }

    public Value<Double> getKappa() {
        return getParams().get(kappaParamName);
    }

    public Value<Double[]> getFreq() {
        return getParams().get(freqParamName);
    }

    private Double[][] hky(double kappa, Double[] freqs, double rate) {

        int numStates = 4;
        
        Double[][] Q = new Double[numStates][numStates];

        double[] totalRates = new double[numStates];

        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (i != j) {
                    if (Math.abs(i-j) == 2) {
                        Q[i][j] = kappa * freqs[j];
                    } else {
                        Q[i][j] = freqs[j];
                    }
                } else Q[i][i] = 0.0;
                totalRates[i] += Q[i][j];
            }
            Q[i][i] = -totalRates[i];
        }

        // normalise rate matrix to rate
        normalize(freqs, Q, rate);

        return Q;
    }
}