package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxon;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.*;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

/**
 * TODO
 *
 * @author Sophie Yang, Remco Bouckaert, Jonathan Klawitter
 */
public class ExactMethodV2 {

    Random random = new Random();

    // static cache tables
    static int binomialCoefficientThreshold = 50;
    private static long[][] combination;
    private static double[][] squareBracketsDivParentheses;
    private static double[][] parenthesesDivFactorials;
    private static double[][][][] pairPickingCache;

    /**
     * Construct new CladePartitionProbabilityComputer with caches set up for fast computation.
     *
     * @param n size of largest parent clade; used to set up caches
     */

    public ExactMethodV2(int n) {
        setUpCache(n);
    }

    private void setUpCache(int n) {
        squareBracketsDivParentheses = new double[n][n];
        for (int x = 1; x < n; x++) {
            double output3 = 1;
            for (int i = 0; i < x; i++) {
                output3 *= (x - i) / (double) (x + i);
                squareBracketsDivParentheses[x][i + 1] = (double) ((i + 1) * 2 - 1.0) * output3;
            }
        }

        parenthesesDivFactorials = new double[n][n + 1];
        for (int y = 0; y < n; y++) {
            for (int n_min_1 = Math.max(y - 1, 0); n_min_1 <= n; n_min_1++) {
                double output = 1;
                int upper = Math.min(y, n_min_1);
                for (int i = 0; i < upper; i++) {
                    output *= (y + i) / (i + 1.0);
                }
                for (int i = y; i < n_min_1; i++) {
                    output *= (y + i) / (n_min_1 + 1.0 - i);
                }
                if (y > n_min_1) {
                    output /= y;
                }
                parenthesesDivFactorials[y][n_min_1] = output;
            }
        }

//	    for (int i = 0; i < n; i++) {
//	    	for (int j = Math.max(i-1,0); j <= n;j++) {
//	    		double p = parenthesesDivFactorials[i][j];
//	    		double t = parentheses[i][j]/(factorial[i]*factorial[j+1-i]);
//	    		System.out.println(i+ " " + j + " " +  p +" " + t + " " + (p-t) + " " + p/t);
//	    		// parenthesesDivFactorials[i][j] = t;
//	    	}
//	    }

        combination = new long[binomialCoefficientThreshold + 1][binomialCoefficientThreshold + 1];
        for (int i = 0; i <= binomialCoefficientThreshold; i++) {
            for (int j = 0; j <= i; j++) {
                combination[i][j] = CombinatoricsUtils.binomialCoefficient(i, j);
            }
        }

        pairPickingCache = new double[n][n][n][n];
    }

    public double getProbability(double theta, double[] timesLeft, double[] timesRight) {

        // nLeft, nRight: total number of taxa in the left and right clade
        int nLeft = timesLeft.length;
        int nRight = timesRight.length;
        Arrays.sort(timesLeft);
        Arrays.sort(timesRight);

        // samplingTimes: distinct sampling times
        double[] samplingTimes = DoubleStream.concat(Arrays.stream(timesLeft), Arrays.stream(timesRight))
                .distinct().sorted().toArray();

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

        // dimensions of the table of the DP
        // since we can have 0 to n1 active left, we have nLeft + 1 rows
        // since we can have 0 to n2 active right, we have nRight + 1 columns
        int rows = nLeft + 1;
        int cols = nRight + 1;
        // use flattened table for speedup
        // thus T[i][j] = T[i * cols + j]
        double[] T = new double[rows * cols];
        // furthermore, the DP loops over several tables in theory,
        // but we can use a single table and overwrite previous values
        // since any entry T[i,j] at time i only depends
        // on entries T[k,l] at time i-1 with index k >= i and l >=j

        // initializing the table with initially active lineages
        int leftSampledSoFar = sampledLeftAtTime[0];
        int rightSampledSoFar = sampledRightAtTime[0];
        T[leftSampledSoFar * cols + rightSampledSoFar] = 1;

        // main loop of the DP
        for (int i = 0; i < taus.length; i++) {
            double tau = taus[i];

            // rows
            for (int j = 0; j <= nLeft; j++) {
                // if at least one left lineage has been sampled,
                // we can skip row 0, since we cannot go back to zero left lineages
                if (leftSampledSoFar != 0 && j == 0) {
                    continue;
                }

                // columns
                for (int k = 0; k <= rightSampledSoFar; k++) {
                    // if at least one right lineage has been sampled,
                    // we can skip column 0, since we cannot go back to zero right lineages
                    if (rightSampledSoFar != 0 && k == 0) {
                        continue;
                    }

                    double prob = 0;
                    // recall that we can overwrite old entries,
                    // because we only use values to the right and below
                    for (int jPrev = j; jPrev <= leftSampledSoFar; jPrev++) {
                        for (int kPrev = k; kPrev <= rightSampledSoFar; kPrev++) {
                            double prevProb = T[jPrev * cols + kPrev];

                            int x = jPrev + kPrev;
                            int y = j + k;
                            double coalesProb1 = coalescentFormula(x, y, tau, theta);
                            double coalesProb2 = pairPickingFormula(jPrev, kPrev, j, k);

                            prob += prevProb * coalesProb1 * coalesProb2;
                        }
                    }
                    T[j * cols + k] = prob;
                }
            }

            // since we sample new lineages,
            // we shift the whole table correspondingly from bottom right to top left
            int addedLeft = sampledLeftAtTime[i + 1];
            int addedRight = sampledRightAtTime[i + 1];
            leftSampledSoFar += addedLeft;
            rightSampledSoFar += addedRight;
            for (int j = (leftSampledSoFar - addedLeft); j >= 0; j--) {
                if (((leftSampledSoFar - addedLeft) != 0) && (j == 0)) {
                    continue;
                }
                for (int k = (rightSampledSoFar - addedRight); k >= 0; k--) {
                    if (((rightSampledSoFar - addedRight) != 0) && (k == 0)) {
                        continue;
                    }
                    double value = T[j * cols + k];
                    T[j * cols + k] = 0;
                    T[(j + addedLeft) * cols + k + addedRight] = value;
                }
            }
        }

        // we have now reached the contemporaneous case,
        // so we combine all entries  with the probability of combining into the right clades now
        double finalProbability = 0;
        for (int prev_j = 1; prev_j <= leftSampledSoFar; prev_j++) {
            for (int prev_k = 1; prev_k <= rightSampledSoFar; prev_k++) {
                double prev_prob = T[prev_j * cols + prev_k];
                double contemp_prob = contempFormula(prev_j + prev_k, prev_j);
                finalProbability += prev_prob * contemp_prob;
            }
        }

        return finalProbability;
    }

    private static double coalescentFormula(int prevActiveLineages, int nowActiveLineages, double tau, double theta) {
        double prob = 0;
        final double g = tau / theta * 0.5;
        int k = 1;
        if (nowActiveLineages != 1) {
            for (int n = nowActiveLineages; n <= prevActiveLineages; n++) {
                double expPart = Math.exp(-n * (n - 1) * g);
                double numerator = k * squareBracketsDivParentheses[prevActiveLineages][n]
                        * parenthesesDivFactorials[nowActiveLineages][n - 1];
                k = -k;
                prob += expPart * numerator;
            }
        } else {
            for (int n = 2; n <= prevActiveLineages; n++) {
                double expPart = Math.exp(-n * (n - 1) * g);
                double numerator = k * squareBracketsDivParentheses[prevActiveLineages][n];
                k = -k;
                prob += expPart * numerator;
            }
            prob = 1 - prob;
        }
        return prob;
    }

    private static double pairPickingFormula(int jPrev, int kPrev, int j, int k) {
        int dj = jPrev - j;
        int dk = kPrev - k;

        double r = pairPickingCache[jPrev][kPrev][dj][dk];
        if (r == 0) {
            final int prev = jPrev + kPrev;
            r = 1.0;
            for (int m = 0; m < dj; m++) {
                r *= (jPrev - m) * (jPrev - m - 1.0) / ((prev - m) * (prev - m - 1));
                // r *= (dj + dk - m)/(m +1.0);
            }
            for (int n = 0; n < dk; n++) {
                r *= (kPrev - n) * (kPrev - n - 1.0) / ((prev - dj - n) * (prev - dj - n - 1));
                // r *= (dk - n)/(dj + n + 1.0);
            }
            r *= binomialCoeff((dj + dk), dj);
            pairPickingCache[jPrev][kPrev][dj][dk] = r;
        }

        return r;
    }

    public static double contempFormula(int parentSize, int leftSize) {
        return (2.0 / (parentSize - 1)) * (1.0 / binomialCoeff(parentSize, leftSize));
    }

    public static long binomialCoeff(int n, int k) {
        if (n < binomialCoefficientThreshold) {
            return combination[n][k];
        } else {
            long x = 1;
            for (int i = 0; i < (k - 1); i++) {
                x *= (n - i) / (k - i);
            }
            return x;
        }
    }

    public double sample(double theta, List<Double> timesLeft, List<Double> timesRight) {

        double totalWeight = 1;

        double popSize = theta;

        List<Double> samplingTimes = new ArrayList<>();
        for (int i = 0; i < timesLeft.size(); i++) {
            if (samplingTimes.contains(timesLeft.get(i)) == false) {
                samplingTimes.add(timesLeft.get(i));
            }
        }
        for (int i = 0; i < timesRight.size(); i++) {
            if (samplingTimes.contains(timesRight.get(i)) == false) {
                samplingTimes.add(timesRight.get(i));
            }
        }
        samplingTimes.sort(null);

        System.out.println("initial samplingTimes = " + samplingTimes);

        // sort input in reverse order, at coal event, remove the youngest, add the parent to the end of the list
        Collections.sort(timesLeft, Collections.reverseOrder());
        Collections.sort(timesRight, Collections.reverseOrder());

        double time = samplingTimes.get(0); // initial time: the first sampling time
        double lastSamplingTime = samplingTimes.get(samplingTimes.size()-1);
        // no good? sort the time, remove it along the way
        int currentTimeIndex = 0; // for jumping to the next sampling time, initially start from 0

        // sampledLeftAtTime, sampledRightAtTime: how many samples are added at a sampling time
//        int[] sampledLeftAtTime = new int[samplingTimes.size()];
//        int[] sampledRightAtTime = new int[samplingTimes.size()];
//
//        // count and assign sampledLeftAtTime, sampledRightAtTime

//        int iLeft = 0;
////        int iRight = 0;
//        for (int i = 0; i < samplingTimes.size(); i++) {
////            System.out.println("i = " + i + ",  samplingTimes.get(i) = " + samplingTimes.get(i));
//            for (int j = 0; j < timesLeft.size(); j++) {
////                System.out.println("j = "+ j + ",  timesLeft.get(j) = " + timesLeft.get(j));
//                System.out.println("i = " + i + ",  samplingTimes.get(i) = " + samplingTimes.get(i) + ",  j = "+ j + ",  timesLeft.get(j) = " + timesLeft.get(j) + ",  timesLeft.get(j) == samplingTimes.get(i) " + (timesLeft.get(j) == samplingTimes.get(i)));
//                if (timesLeft.get(j) == samplingTimes.get(i)) {
////                    System.out.println("sampledLeftAtTime[i] = " + sampledLeftAtTime[i]);
//                    sampledLeftAtTime[i]++;
//                }
//            }
//            for (int k = 0; k < timesRight.size(); k++) {
//                System.out.println("k = "+ k);
//                if (timesRight.get(k) == samplingTimes.get(i)) {
//                    sampledRightAtTime[i]++;
//                }
//            }
//        }
////            while ((iLeft < timesLeft.size()) && (timesLeft.get(iLeft) == samplingTimes.get(i))) {
//                System.out.println("sampledLeftAtTime[i] = " + sampledLeftAtTime[i]);
//                sampledLeftAtTime[i]++;
//                iLeft++;
//            }
//            while ((iRight < timesRight.size()) && (timesRight.get(iRight) == samplingTimes.get(i))) {
//                sampledRightAtTime[i]++;
//                iRight++;
//            }
//        }

//        System.out.println("sampledLeftAtTime = " + Arrays.toString(sampledLeftAtTime) + ",  sampledRightAtTime = " + Arrays.toString(sampledRightAtTime));

        int activeLeftSize = 0;
        int activeRightSize = 0;
        // initially, the taxa sampled at the first sampling time are active, the rest are in toBeAdded
        int iLeft = timesLeft.size() - 1;
        int iRight = 0;
        while (timesLeft.get(iLeft) == samplingTimes.get(currentTimeIndex)) {
            activeLeftSize ++;
            iLeft --;
        }
        while (timesRight.get(iRight) == samplingTimes.get(currentTimeIndex)) {
            activeRightSize ++;
            iRight --;
        }
        System.out.println("activeLeftSize = " + activeLeftSize + ",  activeRightSize = " + activeRightSize);


//        int activeLeftSize = sampledLeftAtTime[0];
//        int activeRightSize = sampledRightAtTime[0];
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

            System.out.println(" ");

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

//            System.out.println("activeLeftSize = " + activeLeftSize + ",  activeRightSize = " + activeRightSize);
            System.out.println("totalPairCount = " + totalPairCount + ",  validPairCount = " + validPairCount);
            System.out.println("canCreateC = " + canCreateC);

            // if there's only 1 active node left, cannot coalesce. Go to the next sampling time
            if (activeNodesCount == 1) {
                time = samplingTimes.get(currentTimeIndex+1);
            } else { // there are more than 1 active node, can coalesce

                // if there is valid pair, draw time
                if (validPairCount != 0) {
                    double rate = (activeNodesCount * (activeNodesCount - 1.0)) / (popSize * 2.0);
                    double x = -Math.log(random.nextDouble()) / rate;
                    time += x;

                    // if time too large, go to next sampling time
                    if ((toBeAddedLeftSize + toBeAddedRightSize) > 0 &&
                            time > samplingTimes.get(currentTimeIndex+1)) {
                        time = samplingTimes.get(currentTimeIndex+1);
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
                    double nextTime = samplingTimes.get(currentTimeIndex+1);
                    double tau = nextTime - time;
                    double noCoalescentProb = Math.exp(-tau * ((double) totalPairCount / theta));
                    totalWeight *= noCoalescentProb;
                    time = nextTime;
                }
            }

            // The end of this loop, update time index, now it points to the next sampling time
            currentTimeIndex += 1;

            // At the next sampling time, add nodes to activeNodes
            while ((toBeAddedLeftSize + toBeAddedRightSize) > 0 && samplingTimes.get(currentTimeIndex) == time) {
//                activeLeftSize += sampledLeftAtTime[currentTimeIndex];
//                activeRightSize += sampledRightAtTime[currentTimeIndex];
//                toBeAddedLeftSize -= sampledLeftAtTime[currentTimeIndex];
//                toBeAddedRightSize -= sampledRightAtTime[currentTimeIndex];
            }
        }

        totalWeight *= contempFormula((activeLeftSize + activeRightSize), activeLeftSize);

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

    public static void main(String[] args) { // ------------------------------------------------------------------------

        int parentCladeSize = 3;
        ExactMethodV2 testCase = new ExactMethodV2(parentCladeSize);

        double theta = 1.0;
        double[] timesLeft = new double[]{0, 0};
        double[] timesRight = new double[]{0};
        List<Double> timesLeftList = new ArrayList<>(Arrays.asList(0.0, 0.0));
        List<Double> timesRightList = new ArrayList<>(Arrays.asList(1.0));

        // System.out.println(testCase.getProbability(theta, timesLeft, timesRight));
        System.out.println(testCase.sample(theta, timesLeftList, timesRightList));

    }

}


