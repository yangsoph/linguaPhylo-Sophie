package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.util.ArrayList;
import java.util.List;

import java.io.File;
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;

import static lphy.base.evolution.coalescent.ExactMethod.computeProbability;

public class TestCombined {

    // public final static int SEED = 177;

    public static void main(String[] args) {

        // RandomUtils.setSeed(SEED);
        // RandomGenerator random = RandomUtils.getRandom();
        Random random = new Random();

        Value<Number> theta = new Value<>("theta", 1);
        double theta_1 = 1;
        int numOfTaxa = 20;

        int cladeReps = 10;
        int numOfSimReps = 10;

        double[][] simulationResults = new double[numOfSimReps][cladeReps];
        double[][] exactResults = new double[numOfSimReps][cladeReps];

        // reps for simulation progressively increase [8, 16, 32, 64, 128, 256, 512]
        for (int S = 0; S < numOfSimReps; S++) {
            int simReps = (int) Math.pow(2, S + 3);

            for (int C = 0; C < cladeReps; C++) {

                List<Taxon> taxaList = new ArrayList<>(); // all taxa
                List<Taxon> rightClade = new ArrayList<>(); // a copy of all taxa, for removing later
                List<Double> tau = new ArrayList<>();

                for (int i = 1; i <= numOfTaxa; i++) { // create all taxa
                    String name = "taxon" + String.valueOf(i);
                    double time = (double) (i - 1) / (numOfTaxa - 1); // total time = 1
                    taxaList.add(new Taxon(name, time));
                    rightClade.add(new Taxon(name, time));
                    tau.add((double) 1 / (numOfTaxa - 1));
                }
                tau.remove(tau.size() - 1);

                // How many taxa are there in the left clade
                int leftSize = (int) (random.nextDouble() * (numOfTaxa / 2)) + 1; // random integer in range [1, numOfTaxa/2]

                System.out.println("leftSize = " + leftSize);

                List<Taxon> leftClade = new ArrayList<>();

                for (int i = 0; i < leftSize; i++) { // randomly pick half of the taxa to be in the left clade
                    int randomIndex = random.nextInt(rightClade.size());
                    Taxon randomTaxon = rightClade.remove(randomIndex);
                    leftClade.add(randomTaxon);
                }

                // convert list to array
                Taxon[] leftCladeArray = leftClade.toArray(new Taxon[]{});
                Taxon[] rightCladeArray = rightClade.toArray(new Taxon[]{});
                Taxon[] taxaArray = taxaList.toArray(new Taxon[]{});

                exactResults[S][C] = computeProbability(theta_1, leftClade, rightClade, tau);
                System.out.println("Exact clade split prob = " + exactResults[S][C]);

                Value<Taxa> taxa = new Value<>("taxa", new Taxa.Simple(taxaArray));
                Value<Integer> n = new Value<>("n", numOfTaxa);
                Value<Taxa> leftCladeInput = new Value<>("taxa", new Taxa.Simple(leftCladeArray));
                Value<Taxa> rightCladeInput = new Value<>("taxa", new Taxa.Simple(rightCladeArray));

                long starting = System.currentTimeMillis();

                double[] weights = new double[simReps];

                for (int i = 0; i < simReps; i++) {
                    SerialCoalClade coalescent = new SerialCoalClade(theta, n, taxa, null, leftCladeInput, rightCladeInput);
                    RandomVariable<TimeTree> sample = coalescent.sample();
                    weights[i] = sample.value().getRoot().getWeight();
                }

                Mean mean = new Mean();
                simulationResults[S][C] = mean.evaluate(weights);

                StandardDeviation standardDeviation = new StandardDeviation();
                double stderr = standardDeviation.evaluate(weights) / Math.sqrt(simReps);

                System.out.println("Simulation mean weight = " + simulationResults[S][C] + " +/- " + stderr);

                long ending = System.currentTimeMillis();
                System.out.println("Time took: " + ((ending - starting) / 1000) + " seconds");

                double err = ((simulationResults[S][C] - exactResults[S][C]) / exactResults[S][C]) * 100;
                System.out.println("error = " + err + " %");

                /*
                try {
                    File myObj = new File("C:\\Users\\MyName\\filename.txt");
                    FileWriter myWriter = new FileWriter("filename.txt");
                    myWriter.write("Files in Java might be tricky, but it is fun enough!");
                    myWriter.close();
                    System.out.println("Successfully wrote to the file.");
                } catch (IOException e) {
                    System.out.println("An error occurred.");
                    e.printStackTrace();
                }

                 */
            }
        }
    }
}