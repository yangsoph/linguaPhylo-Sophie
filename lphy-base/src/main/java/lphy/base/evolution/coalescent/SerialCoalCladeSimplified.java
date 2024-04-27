package lphy.base.evolution.coalescent;

import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.RandomVariable;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.DoubleStream;


public class SerialCoalCladeSimplified {
    // double theta, double[] timesLeft, double[] timesRight

    Random random = new Random();

    public double sample(double theta, List<Double> timesLeft, List<Double> timesRight) {

        double totalWeight = 1;

        double popSize = theta;

        // nLeft, nRight: total number of taxa in the left and right clade
        int nLeft = timesLeft.size();
        int nRight = timesRight.size();
        // sort in reverse order, at coal event, remove the youngest, add the parent to the end of the list
        Collections.sort(timesLeft, Collections.reverseOrder());
        Collections.sort(timesRight, Collections.reverseOrder());

        // samplingTimes: distinct sampling times
        double[] samplingTimes = DoubleStream.concat((DoubleStream) Arrays.stream(timesLeft), (DoubleStream) Arrays.stream(timesRight))
                .distinct().sorted().toArray();

        double time = samplingTimes[0]; // initial time: the first sampling time
        double lastSamplingTime = samplingTimes[samplingTimes.length - 1];
        // no good? sort the time, remove it along the way
        int currentTimeIndex = 0; // for jumping to the next sampling time, initially start from 0

        // sampledLeftAtTime, sampledRightAtTime: how many samples are added at a sampling time
        int[] sampledLeftAtTime = new int[samplingTimes.length];
        int[] sampledRightAtTime = new int[samplingTimes.length];
        // taus: sampling intervals
        double[] taus = new double[samplingTimes.length - 1];

        // count and assign sampledLeftAtTime, sampledRightAtTime
        int iLeft = 0;
        int iRight = 0;
        for (int i = 0; i < samplingTimes.length; i++) {
            while ((iLeft < timesLeft.length) && (timesLeft[iLeft] == samplingTimes[i])) {
                sampledLeftAtTime[i]++;

                iLeft++;
            }
            while ((iRight < timesRight.length) && (timesRight[iRight] == samplingTimes[i])) {
                sampledRightAtTime[i]++;
                iRight++;
            }
            if (i != (samplingTimes.length - 1)) {
                taus[i] = samplingTimes[i + 1] - samplingTimes[i];
            }
        }

        // initially, the taxa sampled at the first sampling time are active, the rest are in toBeAdded
        int activeLeftSize = sampledLeftAtTime[0];
        int activeRightSize = sampledRightAtTime[0];
        int toBeAddedLeftSize = timesLeft.length - activeLeftSize;
        int toBeAddedRightSize = timesRight.length - activeRightSize;

        int validPairCount;
        int totalPairCount;

        // initialise: number of valid pairs in each group
        int validPairCountLeft;
        int validPairCountRight;

        // initialise: canCreateC, existC. C is the parent node of the left and right clade
        boolean canCreateC = false;

        //------------------------------------------------------------------------------------------------------------------------------
        while (time != lastSamplingTime) {

            // update totalPairCount
            int activeNodesCount = activeLeftSize + activeRightSize;
            totalPairCount = activeNodesCount * (activeNodesCount - 1) / 2; // activeNodesCount choose 2

            // check can create C or not
            if (activeLeftSize == 1 && activeRightSize == 1
                    && toBeAddedLeftSize == 0 && toBeAddedRightSize == 0) {
                canCreateC = true;
            }

            // update the number of valid pairs in each group
            validPairCountLeft = activeLeftSize * (activeLeftSize - 1) / 2; // activeLeft choose 2
            validPairCountRight = activeRightSize * (activeRightSize - 1) / 2; // activeRight choose 2

            // update validPairCount
            if (canCreateC == false) {
                validPairCount = validPairCountLeft + validPairCountRight;
            } else { // canCreateC == true
                validPairCount = 1;
            }

            // if there's only 1 active node left, cannot coalesce. Go to the next sampling time
            if (activeNodesCount == 1) {
                time = samplingTimes[currentTimeIndex+1];
            } else { // there are more than 1 active node, can coalesce

                // if there is valid pair, draw time
                if (validPairCount != 0) {
                    double rate = (activeNodesCount * (activeNodesCount - 1.0)) / (popSize * 2.0);
                    double x = -Math.log(random.nextDouble()) / rate;
                    time += x;

                    // if time too large, go to next sampling time
                    if ((toBeAddedLeftSize + toBeAddedRightSize) > 0 &&
                            time > samplingTimes[currentTimeIndex+1]) {
                        time = samplingTimes[currentTimeIndex+1];
                    }
                    else {
                        // At least one group has >= 2 active lineages that can coalesce
                        if (canCreateC == false) {
                            int randomNum = random.nextInt(validPairCount);
                            if (randomNum < validPairCountLeft) { // left pair coalescence
                                coalescentEvent(timesLeft);
                            } else { // right
                                coalescentEvent(timesRight);
                            }
                            totalWeight *= (double) validPairCount / totalPairCount;
                        } else { // canCreateC == true
                            coalescentEvent(timesLeft, timesRight, timesLeft);
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
    private void coalescentEvent(Double[] activeClade, double time) {
        coalescentEvent(activeClade, activeClade, activeClade, time);
    }

    /* Helper method */
    private void coalescentEvent(Double[] firstChildClade, Double[] secondChildClade, Double[] parentClade, double time) {
        Double parent;
        Double firstChild;
        Double secondChild;
        firstChild = firstChildClade.remove(random.nextInt(firstChildClade.size()));
        secondChild = secondChildClade.remove(random.nextInt(secondChildClade.size()));
        parent = new TimeTreeNode(time, new TimeTreeNode[]{firstChild, secondChild});
        parentClade.add(parent);
    }

}
