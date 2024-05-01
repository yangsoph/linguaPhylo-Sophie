package lphy.base.evolution.coalescent;

import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.*;
import java.util.stream.DoubleStream;


public class SerialCoalCladeSimplified {
    // double theta, double[] timesLeft, double[] timesRight

    Random random = new Random();

    public double sample(double theta, List<Double> timesLeft, List<Double> timesRight) {

        double totalWeight = 1;

        double popSize = theta;

        // sort input in reverse order, at coal event, remove the youngest, add the parent to the end of the list
        Collections.sort(timesLeft, Collections.reverseOrder());
        Collections.sort(timesRight, Collections.reverseOrder());

        // samplingTimes: distinct sampling times
        double[] samplingTimes = DoubleStream.concat((DoubleStream) timesLeft, (DoubleStream) timesRight)
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
            while ((iLeft < timesLeft.size()) && (timesLeft.get(iLeft) == samplingTimes[i])) {
                sampledLeftAtTime[i]++;
                iLeft++;
            }
            while ((iRight < timesRight.size()) && (timesRight.get(iRight) == samplingTimes[i])) {
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
        int toBeAddedLeftSize = timesLeft.size() - activeLeftSize;
        int toBeAddedRightSize = timesRight.size() - activeRightSize;

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
                                coalescentEvent(timesLeft, time);
                                activeLeftSize -= 1; // remove two children, add one parent
                            } else { // right
                                coalescentEvent(timesRight, time);
                                activeRightSize -= 1; // remove two children, add one parent
                            }
                            totalWeight *= (double) validPairCount / totalPairCount;
                        }
                        // Now this case might never happen? Bcs they will be contemp case
                        else { // canCreateC == true.
                            coalescentEvent(timesLeft, timesRight, timesLeft, time);
                            activeRightSize -= 1;
                        }
                    }
                } else { // no valid pari, go to next time, p_nc
                    // the probability that no coalescent before next sample, exp( -tau * (k choose 2) / theta )
                    double nextTime = samplingTimes[currentTimeIndex+1];
                    double tau = nextTime - time;
                    double noCoalescentProb = Math.exp(-tau * ((double) totalPairCount / theta));
                    totalWeight *= noCoalescentProb;
                    time = nextTime;
                }
            }

            // The end of this loop, update time index, now it points to the next sampling time
            currentTimeIndex += 1;

            // At the next sampling time, add nodes to activeNodes
            while ((toBeAddedLeftSize + toBeAddedRightSize) > 0 && samplingTimes[currentTimeIndex] == time) {
                activeLeftSize += sampledLeftAtTime[currentTimeIndex];
                activeRightSize += sampledRightAtTime[currentTimeIndex];
                toBeAddedLeftSize -= sampledLeftAtTime[currentTimeIndex];
                toBeAddedRightSize -= sampledRightAtTime[currentTimeIndex];
            }
        }

//        totalWeight *= contempFormula((activeLeftSize + activeRightSize), activeLeftSize);

        return totalWeight;
    }

    /* Helper method */
    private void coalescentEvent(List<Double> activeClade, double time) {
        coalescentEvent(activeClade, activeClade, activeClade, time);
    }

    /* Helper method */
    private void coalescentEvent(List<Double> firstChildClade, List<Double> secondChildClade, List<Double> parentClade, double time) {
        firstChildClade.remove(firstChildClade.size()-1);
        secondChildClade.remove(secondChildClade.size()-1);
        parentClade.add(time); // use the coalescent time to represent the parent
    }

//    public static double contempFormula(int parentSize, int leftSize) {
//        return (2.0 / (parentSize - 1)) * (1.0 / binomialCoeff(parentSize, leftSize));
//    }

//    public static long binomialCoeff(int n, int k) {
//        if (n < binomialCoefficientThreshold) {
//            return combination[n][k];
//        } else {
//            long x = 1;
//            for (int i = 0; i < (k - 1); i++) {
//                x *= (n - i) / (k - i);
//            }
//            return x;
//        }
//    }

    public static void main(String[] args) { // ------------------------------------------------------------------------

        int parentCladeSize = 3;
//        SerialCoalCladeSimplified testCase = new SerialCoalCladeSimplified(parentCladeSize);

        double theta = 1.0;

        List<Double> timesLeft = new ArrayList<>(Arrays.asList(0.0, 0.0));
        List<Double> timesRight = new ArrayList<>(Arrays.asList(1.0));

//        System.out.println(testCase.sample(theta, timesLeft, timesRight));
    }

}
