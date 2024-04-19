package lphy.base.evolution.coalescent;

import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTreeNode;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.*;

public class ExactMethod {

    public static void main(String[] args) {

        double theta = 1;

        /*
        Contemp case, correct. Got 0.05 i.e. 1/ 20
        Taxon A = new Taxon("A");
        Taxon B = new Taxon("B");
        Taxon C = new Taxon("C");
        Taxon D = new Taxon("D");
        Taxon E = new Taxon("E");
        List<Taxon> leftClade = Arrays.asList(A, B);
        List<Taxon> rightClade = Arrays.asList(C, D, E);
        List<Double> tau = new ArrayList<>();
         */

        Taxon A = new Taxon("A");
        Taxon B = new Taxon("B");
        Taxon C = new Taxon("C", 1);
        List<Taxon> leftClade = Arrays.asList(A, B);
        List<Taxon> rightClade = Arrays.asList(C);
        List<Double> tau = Arrays.asList(1.0);
        // Correct. Got 0.7547470392190384
        /*

         */

        // System.out.println(computeProbability(theta, leftClade, rightClade, tau));

        /*
        long startTime = System.currentTimeMillis();
        for (int i = 0; i < 10000; i++) {
            computeProbability(C1, C2, tau, S, theta);
        }
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("took " + totalTime + " millis");
        // took 2647 millis
        */
    }

    public static double computeProbability(double theta, List<Taxon> leftClade, List<Taxon> rightClade, List<Double> tau) {

        double inf = Double.POSITIVE_INFINITY;

        double time = 0.0;

        int n_C1 = 0; // the number of leaves activated from the beginning till now, in the left clade
        int n_C2 = 0; // ... in the right clade

        List<Taxon> toBeAddedLeft = new ArrayList<>();
        List<Taxon> toBeAddedRight = new ArrayList<>();

        // Initialise: count how many leaves are sampled at time 0 in each clade, add the others to toBeAdded
        for (Taxon leaf : leftClade) {
            if (leaf.getAge() <= time) {
                n_C1 += 1;
            } else {
                toBeAddedLeft.add(leaf);
            }
        }
        for (Taxon leaf : rightClade) {
            if (leaf.getAge() <= time) {
                n_C2 += 1;
            } else {
                toBeAddedRight.add(leaf);
            }
        }
        // REVERSE ORDER - youngest age at end of list
        toBeAddedLeft.sort((o1, o2) -> Double.compare(o2.getAge(), o1.getAge()));
        toBeAddedRight.sort((o1, o2) -> Double.compare(o2.getAge(), o1.getAge()));

        // The number of rows and columns in the table
        int rows = leftClade.size() + 1; // +1 for j = 0
        int cols = rightClade.size() + 1;
        double[][] table = new double[rows][cols]; // initialise the table
        table[n_C1][n_C2] = 1; // that entry has initial probability of 1

        if (tau.size() != 0) {
            for (int idx = 0; idx < tau.size(); idx++) { // loop through tau

                // System.out.println("loop through tau");

                //System.out.println();

                //System.out.println("time = " + time);

                Double tau_i = tau.get(idx);

                // Access each cell (j, k) that needs to be filled
                for (int j = 0; j <= n_C1; j++) {
                    if (n_C1 != 0 && j == 0) continue;
                    for (int k = 0; k <= n_C2; k++) {
                        if (n_C2 != 0 && k == 0) continue;

                        double prob = 0;
                        // The range to look for from previous table
                        for (int prev_j = j; prev_j <= n_C1; prev_j++) {
                            for (int prev_k = k; prev_k <= n_C2; prev_k++) {
                                double prev_prob = table[prev_j][prev_k];
                                double coales_prob = function1(prev_j, prev_k, j, k, tau_i, theta);
                                prob += prev_prob * coales_prob;
                                //System.out.println("prev_prob = " + prev_prob + "  coales_prob = " + coales_prob);
                                // bug: coales_prob is -ve
                            }
                        }
                        table[j][k] = prob;
                        //System.out.println("table[" + j + "][" + k + "] = " + table[j][k]);
                    }
                }

                // updates at the next sampling time
                double timeLeft;
                double timeRight;

                if (toBeAddedLeft.size() > 0) {
                    timeLeft = toBeAddedLeft.get(toBeAddedLeft.size() - 1).getAge();
                } else {
                    timeLeft = inf;
                }
                if (toBeAddedRight.size() > 0) {
                    timeRight = toBeAddedRight.get(toBeAddedRight.size() - 1).getAge();
                } else {
                    timeRight = inf;
                }
                time = Math.min(timeLeft, timeRight);

                //System.out.println("timeLeft = " + timeLeft);
                //System.out.println("timeRight = " + timeRight);

                //System.out.println("update time = " + time);

                List<Taxon> newLeftSamples = new ArrayList<>();
                List<Taxon> newRightSamples = new ArrayList<>();

                // get the new leaves sampled at the new time in each clade
                while (toBeAddedLeft.size() > 0 && toBeAddedLeft.get(toBeAddedLeft.size() - 1).getAge() == time) {
                    Taxon youngest = toBeAddedLeft.remove(toBeAddedLeft.size() - 1);
                    newLeftSamples.add(youngest);
                }
                while (toBeAddedRight.size() > 0 && toBeAddedRight.get(toBeAddedRight.size() - 1).getAge() == time) {
                    Taxon youngest = toBeAddedRight.remove(toBeAddedRight.size() - 1);
                    newRightSamples.add(youngest);
                }

                int n_C1_add = newLeftSamples.size();
                int n_C2_add = newRightSamples.size();

                n_C1 += n_C1_add;
                n_C2 += n_C2_add;

                // shift in decreasing index order, so won't overwrite the values
                // newLeftCount: shift down, (j + newLeftCount, k). newRightCount: shift right, (j, k + newRightCount)
                for (int j = n_C1 - n_C1_add; j >= 0; j--) {
                    if (n_C1 - n_C1_add != 0 && j == 0) continue;
                    for (int k = n_C2 - n_C2_add; k >= 0; k--) {
                        if (n_C2 - n_C2_add != 0 && k == 0) continue;
                        double value = table[j][k];
                        table[j][k] = 0;
                        table[j + n_C1_add][k + n_C2_add] = value;
                    }
                }
            }
        }

        // All contemporaneous, or after the last sampling time. Only need the cell (j=1, k=1)
        double finalProb = 0;
        for (int prev_j = 1; prev_j <= n_C1; prev_j++) {
            for (int prev_k = 1; prev_k <= n_C2; prev_k++) {
                double prev_prob = table[prev_j][prev_k];
                double contemp_prob = contempFormula(prev_j + prev_k, prev_j);
                finalProb += prev_prob * contemp_prob;
            }
        }
        return finalProb;
    }

    public static double function1(int prev_j, int prev_k, int j, int k, Double tau_i, double theta) {
        int dj = prev_j - j;
        int dk = prev_k - k;

        int x = prev_j + prev_k;
        int y = j + k;
        double result1 = coalesFormula(x, y, tau_i, theta);

        long orders = combination(dj + dk, dj);

        double result2 = 1;
        for (int m = 0; m < dj; m++) {
            result2 *= (prev_j - m) * (prev_j - m - 1);
        }
        for (int n = 0; n < dk; n++) {
            result2 *= (prev_k - n) * (prev_k - n - 1);
        }
        for (int r = 0; r < dj + dk; r++) {
            result2 /= (prev_j + prev_k - r) * (prev_j + prev_k - r - 1);
        }

        // System.out.println("result1 = " + result1 + ";;  orders = " + orders + ";;  result2 = " + result2);
        // bug: result1 is -ve. I.e. coalesFormula gives -ve probability

        return result1 * orders * result2;
    }

    public static double contempFormula(int C_size, int C1_size) {
        return (2.0 / (C_size - 1)) * (1.0 / combination(C_size, C1_size));
    }

    public static long combination(int n, int k) {
        return CombinatoricsUtils.binomialCoefficient(n, k);
    }

    // coalesFormula: the probability that x lineages coalesce to y lineages in time T = tau / theta
    public static double coalesFormula(int x, int y, Double tau, double theta) {
        double output = 0;
        if (y != 1) {
            // System.out.println("y != 1");
            for (int k = y; k <= x; k++) {

                double expPart = Math.exp(-k * (k - 1) * (tau / theta) * 0.5);

                double term1 = (2 * k - 1);
                double term2 = Math.pow(-1, k - y);
                double term3 = squareBrackets(x, k);
                double term4 = parentheses(y, (k - 1));
                // System.out.println("term1 = " + term1 + "  term2 = " + term2 + "  term3 = " + term3 + "  term4 = " + term4);
                // double numerator = (2 * n - 1) * Math.pow(-1, n - y) * squareBrackets(x, n) * parentheses(y, n - 1);
                double numerator = term1 * term2 * term3 * term4;

                // double denominator = CombinatoricsUtils.factorial(y) * CombinatoricsUtils.factorial((n - y)) * parentheses(x, n);
                // double denominator = term5 * term6 * term7;
                // System.out.println("expPart * numerator / denominator = " + (expPart * numerator / denominator));

                double part1 = divideFactorial(1, y);
                double part2 = divideFactorial(1, (k - y));
                double part3 = (expPart * numerator) / parentheses(x, k);

                output += (part1 * part2 * part3);
                // System.out.println("inside y != 1, output = " + output);
                // System.out.println("added = " + (part1 * part2 * part3));
                // System.out.println("expPart = " + expPart + "  numerator = " + numerator + "  denominator = " + denominator);
            }
        } else {
            // System.out.println("y == 1");
            for (int k = 2; k <= x; k++) {
                double expPart = Math.exp(-k * (k - 1) * tau / theta * 0.5);
                double numerator = (2 * k - 1) * Math.pow(-1, k) * squareBrackets(x, k);
                double denominator = parentheses(x, k);
                output += expPart * numerator / denominator;
            }
            output = 1 - output;
        }
        // System.out.println("coalesFormula output = " + output);
        // coalesFormula output = -216783.46933198246 (y != 1)
        return output;
    }

    private static long squareBrackets(int x, int k) {
        long output = 1;
        for (int i = 0; i < k; i++) {
            output *= (x - i);
        }
        return output;
    }

    private static long parentheses(int x, int k) {
        long output = 1;
        for (int i = 0; i < k; i++) {
            output *= (x + i);
        }
        return output;
    }

    public static double divideFactorial(double x, double y) {
        double factorial = 1;
        for (int i = 2; i <= y; i++) {
            factorial *= i;
        }
        return (x / factorial);
    }

}
