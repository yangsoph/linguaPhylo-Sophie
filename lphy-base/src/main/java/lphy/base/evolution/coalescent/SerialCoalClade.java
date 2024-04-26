package lphy.base.evolution.coalescent;

import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TaxaConditionedTreeGenerator;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.util.*;

import static lphy.base.evolution.coalescent.CoalParamNames.thetaParamName;

public class SerialCoalClade extends TaxaConditionedTreeGenerator {

    private Value<Number> theta;
    private Value<Taxa> leftClade;
    private Value<Taxa> rightClade;

    public SerialCoalClade
            (@ParameterInfo(name = thetaParamName, narrativeName = "coalescent parameter", description = "effective population size, possibly scaled to mutations or calendar units.") Value<Number> theta,
             @ParameterInfo(name = DistributionConstants.nParamName, description = "number of taxa.", optional = true) Value<Integer> n,
             @ParameterInfo(name = TaxaConditionedTreeGenerator.taxaParamName, description = "Taxa object, (e.g. Taxa or TimeTree or Object[])", optional = true) Value<Taxa> taxa,
             @ParameterInfo(name = TaxaConditionedTreeGenerator.agesParamName, description = "an array of leaf node ages.", optional = true) Value<Double[]> ages,
             @ParameterInfo(name = "leftClade", description = "left clade of clade split constraint") Value<Taxa> leftClade,
             @ParameterInfo(name = "rightClade", description = "right clade of clade split constraint.") Value<Taxa> rightClade) {

        super(n, taxa, ages);

        this.theta = theta;
        this.ages = ages;
        this.leftClade = leftClade;
        this.rightClade = rightClade;

        super.checkTaxaParameters(true);
    }

    public RandomVariable<TimeTree> sample() {

        double totalWeight = 1;

        TimeTree tree = new TimeTree();

        List<TimeTreeNode> leafNodes = createLeafTaxa(tree);

        // active nodes in the three groups
        List<TimeTreeNode> activeLeft = new ArrayList<>();
        List<TimeTreeNode> activeRight = new ArrayList<>();

        // leaves to be added, all and in the three groups
        List<TimeTreeNode> leavesToBeAdded = new ArrayList<>();
        List<TimeTreeNode> toBeAddedLeft = new ArrayList<>();
        List<TimeTreeNode> toBeAddedRight = new ArrayList<>();

        double time = 0.0;
        double popSize = ValueUtils.doubleValue(theta);

        // initialise: mark leaves with clade
        for (TimeTreeNode leaf : leafNodes) {
            if (leftClade.value().indexOfTaxon(leaf.getId()) != -1) {
                leaf.setMetaData("clade", "left");
            } else {
                leaf.setMetaData("clade", "right");
            }
        }

        // initialise: divide leaves into activeLeft, activeRight, activeRest, leavesToBeAdded
        // initialise: divide leavesToBeAdded into left, right and rest
        for (TimeTreeNode leaf : leafNodes) {
            if (leaf.getAge() <= time) {
                if (leaf.getMetaData("clade").equals("left")) {
                    activeLeft.add(leaf);
                } else {
                    activeRight.add(leaf);
                }
            } else {
                leavesToBeAdded.add(leaf);
                if (leaf.getMetaData("clade").equals("left")) {
                    toBeAddedLeft.add(leaf);
                } else {
                    toBeAddedRight.add(leaf);
                }
            }
        }

        // REVERSE ORDER - youngest age at end of list
        leavesToBeAdded.sort((o1, o2) -> Double.compare(o2.getAge(), o1.getAge()));
        toBeAddedLeft.sort((o1, o2) -> Double.compare(o2.getAge(), o1.getAge()));
        toBeAddedRight.sort((o1, o2) -> Double.compare(o2.getAge(), o1.getAge()));

        int validPairCount = 0;
        int totalPairCount;

        // initialise: number of valid pairs in each group
        int validPairCountLeft = 0;
        int validPairCountRight = 0;

        // initialise: canCreateC, existC. C is the parent node of the left and right clade
        boolean canCreateC = false;

        //------------------------------------------------------------------------------------------------------------------------------
        while ((activeLeft.size() + activeRight.size() + leavesToBeAdded.size()) > 1) {

            // update totalPairCount
            int activeNodesCount = activeLeft.size() + activeRight.size();
            totalPairCount = activeNodesCount * (activeNodesCount - 1) / 2; // activeNodesCount choose 2

            // check can create C or not
            if (activeLeft.size() == 1 && activeRight.size() == 1
                    && toBeAddedLeft.size() == 0 && toBeAddedRight.size() == 0) {
                canCreateC = true;
            }

            // update the number of valid pairs in each group
            validPairCountLeft = activeLeft.size() * (activeLeft.size() - 1) / 2; // activeLeft choose 2
            validPairCountRight = activeRight.size() * (activeRight.size() - 1) / 2; // activeRight choose 2

            // update validPairCount
            if (canCreateC == false) {
                validPairCount = validPairCountLeft + validPairCountRight;
            } else { // canCreateC == true
                validPairCount = 1;
            }

            // if there's only 1 active node left, cannot coalesce. Go to the next sampling time
            if (activeNodesCount == 1) {
                time = leavesToBeAdded.get(leavesToBeAdded.size() - 1).getAge();
            } else { // there are more than 1 active node, can coalesce

                // if there is valid pair, draw time
                if (validPairCount != 0) {
                    double rate = (activeNodesCount * (activeNodesCount - 1.0)) / (popSize * 2.0);
                    double x = -Math.log(random.nextDouble()) / rate;
                    time += x;

                    // if time too large, go to next sampling time
                    if (leavesToBeAdded.size() > 0 && time > leavesToBeAdded.get(leavesToBeAdded.size() - 1).getAge()) {
                        time = leavesToBeAdded.get(leavesToBeAdded.size() - 1).getAge();
                    }
                    else {
                        // At least one group has >= 2 active lineages that can coalesce
                        if (canCreateC == false) {
                            int randomNum = random.nextInt(validPairCount);
                            if (randomNum < validPairCountLeft) { // left pair coalescence
                                coalescentEvent(activeLeft, time);
                            } else { // right
                                coalescentEvent(activeRight, time);
                            }
                            totalWeight *= (double) validPairCount / totalPairCount;
                        } else { // canCreateC == true
                            coalescentEvent(activeLeft, activeRight, activeLeft, time);
                        }
                    }
                } else { // no valid pari, go to next time, p_nc
                    // the probability that no coalescent before next sample, exp( -tau * (k choose 2) / theta )
                    double nextTime = leavesToBeAdded.get(leavesToBeAdded.size() - 1).getAge();
                    double tau = nextTime - time;
                    double noCoalescentProb = Math.exp(-tau * ((double) totalPairCount / theta.value().doubleValue()));
                    totalWeight *= noCoalescentProb;
                    time = nextTime;
                }
            }

            // At the next sampling time, add nodes to activeNodes
            while (leavesToBeAdded.size() > 0 && leavesToBeAdded.get(leavesToBeAdded.size() - 1).getAge() == time) {
                TimeTreeNode youngest = leavesToBeAdded.remove(leavesToBeAdded.size() - 1);
                if (youngest.getMetaData("clade").equals("left")) {
                    activeLeft.add(youngest);
                    toBeAddedLeft.remove(toBeAddedLeft.size() - 1);
                } else {
                    activeRight.add(youngest);
                    toBeAddedRight.remove(toBeAddedRight.size() - 1);
                }
            }
        }

        tree.setRoot(activeLeft.get(0));
        TimeTreeNode root = tree.getRoot();
        root.setWeight(totalWeight);

        return new RandomVariable<>("\u03C8", tree, this);
    }

    /* Helper method */
    private void coalescentEvent(List<TimeTreeNode> activeClade, double time) {
        coalescentEvent(activeClade, activeClade, activeClade, time);
    }

    /* Helper method */
    private void coalescentEvent(List<TimeTreeNode> firstChildClade, List<TimeTreeNode> secondChildClade, List<TimeTreeNode> parentClade, double time) {
        TimeTreeNode parent;
        TimeTreeNode firstChild;
        TimeTreeNode secondChild;
        firstChild = firstChildClade.remove(random.nextInt(firstChildClade.size()));
        secondChild = secondChildClade.remove(random.nextInt(secondChildClade.size()));
        parent = new TimeTreeNode(time, new TimeTreeNode[]{firstChild, secondChild});
        parentClade.add(parent);
    }

    public Value<Number> getTheta() {
        return theta;
    }

    public static void main(String[] args) {
        Value<Number> theta = new Value<>("theta", 1);

        Taxon A = new Taxon("A");
        Taxon B = new Taxon("B", 0.5);
        Taxon C = new Taxon("C", 1);

        Value<Taxa> taxa = new Value<>("taxa", new Taxa.Simple(new Taxon[]{A, B, C}));
        Value<Integer> n = new Value<>("n", taxa.value().ntaxa());
        Value<Taxa> leftClade = new Value<>("taxa", new Taxa.Simple(new Taxon[]{A, C}));
        Value<Taxa> rightClade = new Value<>("taxa", new Taxa.Simple(new Taxon[]{B}));

        int reps = 100;
        double[] weights = new double[reps];

        for (int i = 0; i < reps; i++) {
            //System.out.println(" ");
            SerialCoalClade coalescent = new SerialCoalClade(theta, n, taxa, null, leftClade, rightClade);
            RandomVariable<TimeTree> sample = coalescent.sample();
            weights[i] = sample.value().getRoot().getWeight();
        }

        Mean mean = new Mean();
        double meanWeight = mean.evaluate(weights);
        StandardDeviation standardDeviation = new StandardDeviation();
        double stderr = standardDeviation.evaluate(weights) / Math.sqrt(reps);

        System.out.println("Mean weight = " + meanWeight + " +/- " + stderr);

        // Exact clade split prob = 0.2021768865708778

        double exactResult = 0.2021768865708778;
        double err = ((meanWeight - exactResult) / exactResult) * 100;
        System.out.println("error = " + err + " %");
    }
}







