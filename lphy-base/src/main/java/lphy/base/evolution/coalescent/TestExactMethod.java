package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxon;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static lphy.base.evolution.coalescent.ExactMethod.computeProbability;

public class TestExactMethod {

    public final static int SEED = 123;
    // SEED = 123, Clade split prob = 3.021760104493697E-7

    public static void main(String[] args) {

        int x = 1;

        for (int j = 0; j < x; j++) {

            Random random = new Random(SEED);
            /// Random random = new Random();
            double theta = 1;

            int numOfTaxa = 20;

            List<Taxon> rightClade = new ArrayList<>();
            List<Double> tau = new ArrayList<>();

            for (int i = 1; i <= numOfTaxa; i++) { // create all taxa
                String name = "taxon" + String.valueOf(i);
                Double time = (double) (i - 1) / (numOfTaxa - 1); // so that total time add up to 1, regardless of numOfTaxa
                rightClade.add(new Taxon(name, time));
                tau.add((double) 1 / (numOfTaxa - 1)); // but this will not work if tau is randomly drawn from an interval
            }
            tau.remove(tau.size() - 1);
            //System.out.println(tau);

            // How many taxa are there in the left clade
            int leftSize = (int) (random.nextDouble() * (numOfTaxa / 2)) + 1; // random integer in range [1, numOfTaxa/2]

            List<Taxon> leftClade = new ArrayList<>();

            for (int i = 0; i < leftSize; i++) { // randomly pick half of the taxa to be in the left clade
                int randomIndex = random.nextInt(rightClade.size());
                Taxon randomTaxon = rightClade.remove(randomIndex);
                leftClade.add(randomTaxon);
            }

            double cladeSplitProb = computeProbability(theta, leftClade, rightClade, tau);

            System.out.println("Clade split prob = " + cladeSplitProb);

            /*
            if (cladeSplitProb <= 0) {

                System.out.println("cladeSplitProb <= 0, num of taxa = " + numOfTaxa);
                System.out.println("Clade split prob = " + cladeSplitProb);
                System.out.println("leftClade = " + leftClade);
                System.out.println("rightClade = " + rightClade);
                System.out.println("tau = " + tau);
                System.out.println();

            }

             */
        }
    }
}
