package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static lphy.base.evolution.coalescent.ExactMethod.computeProbability;

public class TestHalfCladeSize {

    Random random = new Random();

    // blockSize <= numOfTaxa / 2 ; leftCladeSize -- random or half
    public static List<Number> nonRandomTestCase(double thetaExact, int simReps, int numOfTaxa, int blockSize) {

        List<Number> results = new ArrayList<>();

        double err = 0;

        List<Double> tauList = new ArrayList<>();
        List<Taxon> taxaList = new ArrayList<>(); // all taxa
        List<Taxon> taxaListCopy = new ArrayList<>(); // a copy of all taxa, for removing later
        List<Taxon> leftCladeList = new ArrayList<>();
        List<Taxon> rightCladeList = new ArrayList<>();

        for (int i = numOfTaxa; i >= 1; i--) { // create all taxa
            String name = "taxon" + String.valueOf(i);
            double time = (double) (i - 1) / (numOfTaxa - 1); // total time = 1
            taxaList.add(new Taxon(name, time));
            taxaListCopy.add(new Taxon(name, time));
            tauList.add((double) 1 / (numOfTaxa - 1));
        }
        tauList.remove(tauList.size() - 1);

        // e.g. taxa [taxa1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; blockSize = 4
        for (int i = 0; i < (Math.ceil((double) numOfTaxa/blockSize)); i++) { // i = [0, 1, 2]
            if (i != Math.ceil((double) numOfTaxa/blockSize) - 1) { // if it is not the last block
                for (int j = i * blockSize; j < (i * blockSize + blockSize); j++) { // j = [0, 1, 2, 3], [4, 5, 6, 7]
//                  System.out.println("i = " + i + ",  j = " + j);
                    if (i % 2 == 0) { leftCladeList.add(taxaListCopy.remove(taxaListCopy.size()-1)); }
                    else { rightCladeList.add(taxaListCopy.remove(taxaListCopy.size()-1)); }
                }
            }
            else { // if it is the last block, the taxa left may not be enough to fill a whole block
                for (int k = 0; k < taxaListCopy.size(); k++) { // k get taxa 9 & 10
//                  System.out.println("i = " + i + ",  k = " + k);
                    if (i % 2 == 0) { leftCladeList.add(taxaListCopy.get(k)); }
                    else { rightCladeList.add(taxaListCopy.get(k)); }
                }
            }
        }

//      System.out.println("leftCladeList = " + leftCladeList);
//      System.out.println("rightCladeList = " + rightCladeList);

            // convert list to array
        Taxon[] leftCladeArray = leftCladeList.toArray(new Taxon[]{});
        Taxon[] rightCladeArray = rightCladeList.toArray(new Taxon[]{});
        Taxon[] taxaArray = taxaList.toArray(new Taxon[]{});

        // exact method input
        double exactResults = computeProbability(thetaExact, leftCladeList, rightCladeList, tauList);
        System.out.println(" ");
        System.out.println("Exact clade split prob = " + exactResults);

        // simulation input
        Value<Number> thetaSim = new Value<>("theta", thetaExact);
        Value<Taxa> taxaSim = new Value<>("taxa", new Taxa.Simple(taxaArray));
        Value<Integer> nSim = new Value<>("n", numOfTaxa);
        Value<Taxa> leftCladeSim = new Value<>("taxa", new Taxa.Simple(leftCladeArray));
        Value<Taxa> rightCladeSim = new Value<>("taxa", new Taxa.Simple(rightCladeArray));

        long starting = System.currentTimeMillis();

        // for each number of simulation, run that many times, average
        double[] weights = new double[simReps];
        for (int i = 0; i < simReps; i++) {
            SerialCoalClade coalescent = new SerialCoalClade(thetaSim, nSim, taxaSim, null, leftCladeSim, rightCladeSim);
            RandomVariable<TimeTree> sample = coalescent.sample();
            weights[i] = sample.value().getRoot().getWeight();
        }
        Mean mean = new Mean();
        double meanWeight = mean.evaluate(weights);
        StandardDeviation standardDeviation = new StandardDeviation();
        double stderr = standardDeviation.evaluate(weights) / Math.sqrt(simReps);

        long ending = System.currentTimeMillis();

        System.out.println("num of sim = " + simReps);
        System.out.println("Simulation mean weight = " + meanWeight + " +/- " + stderr);
        double time = (ending - starting) / 1000;
        System.out.println("Time took: " + time + " seconds");
        err = ((meanWeight - exactResults) / exactResults) * 100;
        System.out.println("error between exact and sim results= " + err + " %");

        // numOfTaxa, simReps, exactResult, simResult, error, simTime
        results.add(numOfTaxa);
        results.add(simReps);
        results.add(exactResults);
        results.add(meanWeight);
        results.add(err);
        results.add(time);

        return results;
    }


    public static void main(String[] args) throws IOException {
        double thetaExact = 1.0;
        int simRepsExp = 10;
        int numOfTaxa = 10;
        int blockSize = 4;

        int eachSimRepRepeat = 10;
        // double[][] eachSimRepErr = new double[simRepsExp][eachSimRepRepeat];

        for (int expNum = 0; expNum < simRepsExp; expNum++) { // for different number of simulations
            int simReps = (int) Math.pow(2, expNum + 3);
            double[] eachSimRepErr = new double[eachSimRepRepeat];
            for (int i = 0; i < eachSimRepRepeat; i++) {
                List<Number> results = nonRandomTestCase(thetaExact, simReps, numOfTaxa, blockSize);
                eachSimRepErr[i] = (double) results.get(4);
                Mean mean = new Mean();
                double eachSimRepErrMean = mean.evaluate(eachSimRepErr);
            }


//            for (int i = 0; i < eachSimRepRepeat; i++) {
//                double err = nonRandomTestCase(thetaExact, simReps, numOfTaxa, blockSize);
//            }
        }
    }
}

