package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static lphy.base.evolution.coalescent.ExactMethod.computeProbability;

public class TestSimpleCase {
    public static void main(String[] args) {

        Random random = new Random();

        Value<Number> theta = new Value<>("theta", 1);
        double theta_1 = 1;

        // case 1
        Taxon A = new Taxon("A");
        Taxon B = new Taxon("B", 0.5);
        Taxon C = new Taxon("C", 1.0);

        List<Double> tau = Arrays.asList(0.5, 0.5);

        List<Taxon> left1 = Arrays.asList(A, C);
        List<Taxon> right1 = Arrays.asList(B);

        double exact1 = computeProbability(theta_1, left1, right1, tau);
        System.out.println("Case1: Exact clade split prob = " + exact1);

        Value<Taxa> taxa = new Value<>("taxa", new Taxa.Simple(new Taxon[]{A, B, C}));
        Value<Integer> n = new Value<>("n", taxa.value().ntaxa());
        Value<Taxa> leftClade1 = new Value<>("taxa", new Taxa.Simple(new Taxon[]{A, C}));
        Value<Taxa> rightClade1 = new Value<>("taxa", new Taxa.Simple(new Taxon[]{B}));

        int simRep = 100000;
        double[] weights1 = new double[simRep];
        for (int i = 0; i < simRep; i++) {
            SerialCoalClade coalescent = new SerialCoalClade(theta, n, taxa, null, leftClade1, rightClade1);
            RandomVariable<TimeTree> sample = coalescent.sample();
            weights1[i] = sample.value().getRoot().getWeight();
        }
        Mean mean1 = new Mean();
        double meanWeight1 = mean1.evaluate(weights1);
        StandardDeviation standardDeviation1 = new StandardDeviation();
        double stderr1 = standardDeviation1.evaluate(weights1) / Math.sqrt(simRep);
        System.out.println("Mean weight = " + meanWeight1 + " +/- " + stderr1);

        double err1 = ((meanWeight1 - exact1) / exact1) * 100;
        System.out.println("error = " + err1 + " %");

        // case 2
        List<Taxon> left2 = Arrays.asList(A, B);
        List<Taxon> right2 = Arrays.asList(C);

        double exact2 = computeProbability(theta_1, left2, right2, tau);
        System.out.println("Case2: Exact clade split prob = " + exact2);

        Value<Taxa> leftClade2 = new Value<>("taxa", new Taxa.Simple(new Taxon[]{A, B}));
        Value<Taxa> rightClade2 = new Value<>("taxa", new Taxa.Simple(new Taxon[]{C}));

        double[] weights2 = new double[simRep];
        for (int i = 0; i < simRep; i++) {
            SerialCoalClade coalescent = new SerialCoalClade(theta, n, taxa, null, leftClade2, rightClade2);
            RandomVariable<TimeTree> sample = coalescent.sample();
            weights2[i] = sample.value().getRoot().getWeight();
        }
        Mean mean2 = new Mean();
        double meanWeight2 = mean2.evaluate(weights2);
        StandardDeviation standardDeviation2 = new StandardDeviation();
        double stderr2 = standardDeviation2.evaluate(weights2) / Math.sqrt(simRep);
        System.out.println("Mean weight = " + meanWeight2 + " +/- " + stderr2);

        double err2 = ((meanWeight2 - exact2) / exact2) * 100;
        System.out.println("error = " + err2 + " %");

    }
}
