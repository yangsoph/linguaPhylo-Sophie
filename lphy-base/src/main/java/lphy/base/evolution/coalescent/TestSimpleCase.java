package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import java.io.File;
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;

import static lphy.base.evolution.coalescent.ExactMethod.computeProbability;

public class TestSimpleCase {

    public static void main(String[] args) throws IOException {

        Random random = new Random();

        Value<Number> thetaSim = new Value<>("theta", 1);
        double thetaExact = ValueUtils.doubleValue(thetaSim);
        // int numOfTaxa = 20;

        int cladeReps = 100;
        int simRepsExp = 15;

        double[][] exactResults = new double[simRepsExp][cladeReps];
        double[][] simResults = new double[simRepsExp][cladeReps];

        File outputFile = new File("/Users/zyan598/Documents/GitHub/CCD_prior/testing/output_error_time_varyTaxa.csv");
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);

        // header: numOfTaxa, simReps, exactResult, simResult, error, simTime
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("numOfTaxa").append(separator);
        sb.append("simReps").append(separator);
        sb.append("exactResult").append(separator);
        sb.append("simResult").append(separator);
        sb.append("error").append(separator);
        sb.append("exactTime").append(separator);
        sb.append("simTime").append(separator);
        writer.println(sb.toString());

        for (int numOfTaxa = 5; numOfTaxa <= 20; numOfTaxa += 5) {
            for (int expNum = 4; expNum < simRepsExp; expNum++) { // for different number of simulations
                int simReps = (int) Math.pow(2, expNum);
                //int simReps = 10000;

                for (int cReps = 0; cReps < cladeReps; cReps++) { // for different clade splits

                    List<Double> tauList = new ArrayList<>();
                    List<Taxon> taxaList = new ArrayList<>(); // all taxa
                    List<Taxon> rightCladeList = new ArrayList<>(); // a copy of all taxa, for removing later
                    List<Taxon> leftCladeList = new ArrayList<>();

                    // Draw a number for left clade size
                    int leftSize = (int) (random.nextDouble() * (numOfTaxa / 2)) + 1; // random integer in range [1, numOfTaxa/2]

                    for (int i = 1; i <= numOfTaxa; i++) { // create all taxa
                        String name = "taxon" + String.valueOf(i);
                        double time = (double) (i - 1) / (numOfTaxa - 1); // total time = 1
                        taxaList.add(new Taxon(name, time));
                        rightCladeList.add(new Taxon(name, time));
                        tauList.add((double) 1 / (numOfTaxa - 1));
                    }
                    tauList.remove(tauList.size() - 1);

                    for (int i = 0; i < leftSize; i++) { // randomly pick half of the taxa to be in the left clade
                        int randomIndex = random.nextInt(rightCladeList.size());
                        Taxon randomTaxon = rightCladeList.remove(randomIndex);
                        leftCladeList.add(randomTaxon);
                    }

                    // convert list to array
                    Taxon[] leftCladeArray = leftCladeList.toArray(new Taxon[]{});
                    Taxon[] rightCladeArray = rightCladeList.toArray(new Taxon[]{});
                    Taxon[] taxaArray = taxaList.toArray(new Taxon[]{});

                    long startingExact = System.currentTimeMillis();
                    exactResults[expNum][cReps] = computeProbability(thetaExact, leftCladeList, rightCladeList, tauList);
                    long endingExact = System.currentTimeMillis();
                    long timeExact = (endingExact - startingExact);

                    //                System.out.println(" ");
                    //                System.out.println("Exact clade split prob = " + exactResults[expNum][cReps]);

                    Value<Taxa> taxaSim = new Value<>("taxa", new Taxa.Simple(taxaArray));
                    Value<Integer> nSim = new Value<>("n", numOfTaxa);
                    Value<Taxa> leftCladeSim = new Value<>("taxa", new Taxa.Simple(leftCladeArray));
                    Value<Taxa> rightCladeSim = new Value<>("taxa", new Taxa.Simple(rightCladeArray));

                    long startingSim = System.currentTimeMillis();

                    // for each number of simulation, run that many times, average
                    double[] weights = new double[simReps];
                    for (int i = 0; i < simReps; i++) {
                        SerialCoalClade coalescent = new SerialCoalClade(thetaSim, nSim, taxaSim, null, leftCladeSim, rightCladeSim);
                        RandomVariable<TimeTree> sample = coalescent.sample();
                        weights[i] = sample.value().getRoot().getWeight();
                    }
                    Mean mean = new Mean();
                    simResults[expNum][cReps] = mean.evaluate(weights);

                    long endingSim = System.currentTimeMillis();

                    StandardDeviation standardDeviation = new StandardDeviation();
                    double stderr = standardDeviation.evaluate(weights) / Math.sqrt(simReps);

                    //                System.out.println("num of sim = " + simReps);
                    //                System.out.println("Simulation mean weight = " + simResults[expNum][cReps] + " +/- " + stderr);

                    long timeSim = (endingSim - startingSim);
                    //                System.out.println("Time took: " + timeSim + " ms");

                    double err = ((simResults[expNum][cReps] - exactResults[expNum][cReps]) / exactResults[expNum][cReps]) * 100;
                    //                System.out.println("error = " + err + " %");

                    // numOfTaxa, simReps, exactResult, simResult, error, simTime
                    sb = new StringBuilder();
                    sb.append(numOfTaxa).append(separator);
                    sb.append(simReps).append(separator);
                    sb.append(exactResults[expNum][cReps]).append(separator);
                    sb.append(simResults[expNum][cReps]).append(separator);
                    sb.append(err).append(separator);
                    sb.append(timeExact).append(separator);
                    sb.append(timeSim).append(separator);
                    writer.println(sb.toString());
                    writer.flush();

                }
            }
        }
        writer.flush();
        writer.close();
    }
}