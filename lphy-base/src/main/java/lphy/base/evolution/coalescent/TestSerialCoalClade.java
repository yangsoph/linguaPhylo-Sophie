package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class TestSerialCoalClade {
    public static void main(String[] args) {

        Value<Number> theta = new Value<>("theta", 1);

        int numOfTaxa = 50;

        List<Taxon> taxaList = new ArrayList<>(); // all taxa
        List<Taxon> taxaCopy = new ArrayList<>(); // a copy of all taxa

        for (int i = 1; i <= numOfTaxa; i++) { // create all taxa
            String name = "taxon" + String.valueOf(i);
            double time = (double) (i - 1) / (numOfTaxa - 1); // total time = 1
            taxaList.add(new Taxon(name, time));
            taxaCopy.add(new Taxon(name, time));
        }

        // How many taxa are there in the left clade
        int leftSize = (int) (Math.random() * (numOfTaxa / 2)) + 1; // random integer in range [1, numOfTaxa/2]

        Taxon[] leftCladeArray = new Taxon[leftSize];

        Random rand = new Random();
        for (int i = 0; i < leftSize; i++) { // randomly pick half of the taxa to be in the left clade
            int randomIndex = rand.nextInt(taxaCopy.size());
            Taxon randomTaxon = taxaCopy.remove(randomIndex);
            leftCladeArray[i] = randomTaxon;
        }

        // convert list to array
        Taxon[] rightCladeArray = taxaCopy.toArray(new Taxon[]{});
        Taxon[] taxaArray = taxaList.toArray(new Taxon[]{});

        Value<Taxa> taxa = new Value<>("taxa", new Taxa.Simple(taxaArray));
        Value<Integer> n = new Value<>("n", numOfTaxa);
        Value<Taxa> leftClade = new Value<>("taxa", new Taxa.Simple(leftCladeArray));
        Value<Taxa> rightClade = new Value<>("taxa", new Taxa.Simple(rightCladeArray));

        long starting = System.currentTimeMillis();

        int reps = 100000;
        double[] weights = new double[reps];

        for (int i = 0; i < reps; i++) {
            SerialCoalClade coalescent = new SerialCoalClade(theta, n, taxa, null, leftClade, rightClade);
            RandomVariable<TimeTree> sample = coalescent.sample();
            weights[i] = sample.value().getRoot().getWeight();
        }

        Mean mean = new Mean();
        double meanWeight = mean.evaluate(weights);

        StandardDeviation standardDeviation = new StandardDeviation();
        double stderr = standardDeviation.evaluate(weights) / Math.sqrt(reps);

        System.out.println("Mean weight = " + meanWeight + " +/- " + stderr);

        long ending = System.currentTimeMillis();

        System.out.println((ending - starting) / 1000);
    }
}
