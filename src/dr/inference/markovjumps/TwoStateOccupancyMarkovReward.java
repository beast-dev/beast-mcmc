/*
 * TwoStateOccupancyMarkovReward.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.markovjumps;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.math.Binomial;
import dr.math.GammaFunction;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.GeneralizedIntegerGammaDistribution;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc Suchard
 * @author Christophe Fraiser
 * @author Andrew Rambaut
 *
 *
 * https://en.wikipedia.org/wiki/Generalized_integer_gamma_distribution
 *
 */


// This class only works for 2-state models with equal rates, starting and ending in state 0

public class TwoStateOccupancyMarkovReward implements MarkovReward {


    private static final boolean DEBUG = true;
    private static final boolean DEBUG2 = false;

    public TwoStateOccupancyMarkovReward(double[] Q) {
        this(Q, 1E-10);
    }

    public TwoStateOccupancyMarkovReward(double[] Q, double epsilon) {
        this.Q = Q;
        this.maxTime = 0;
        this.epsilon = epsilon;

        this.maxK = 10;  // TODO How to determine?

        eigenSystem = new DefaultEigenSystem(2);

//        if (Q[idx(0, 0)] != Q[idx(1, 1)]) {
//            throw new IllegalArgumentException("Only currently implemented for equal rates models");
//        }
    }

    private int idx(int i, int j) {
        return i * 2 + j; // row-major
    }

    private double[] jumpProbabilities = null;

//    private double getAki(double lambda1, double lambda2, int k, int i) {
////        double logA = Binomial.logChoose(k + i - 1, i) +
////                k * Math.log(lambda1) + (k - 1) * Math.log(lambda2)
////        int sign = 1;
////        return sign * Math.exp(logA);
//        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
//                * Math.pow(lambda2, k - 1) / Math.pow(lambda2 - lambda1, k + i);
//    }
//
//    private double getBki(double lambda1, double lambda2, int k, int i) {
//        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
//                * Math.pow(lambda2, k - 1) / Math.pow(lambda1 - lambda2, k + i);
//    }
//
//    private double getCki(double lambda1, double lambda2, int k, int i) {
//        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
//                * Math.pow(lambda2, k) / Math.pow(lambda2 - lambda1, k + i);
//    }
//
//    private double getDki(double lambda1, double lambda2, int k, int i) {
//        return Binomial.choose(k + i, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
//                * Math.pow(lambda2, k) / Math.pow(lambda1 - lambda2, k + i + 1);
//    }


    private void computeCDForJumpProbabilities(double lambda1, double lambda2, double[][] C, double[][] D) {

        if (lambda1 == lambda2) return;

        for (int k = 0; k < maxK / 2; ++k) {
            final double l1l2k = Math.pow(lambda1 * lambda2, k);

            for (int i = 0; i <= k + 1; ++i) {
                final double sign = Math.pow(-1, i);

                if (k == 0 && i == 0) {
                    C[k][i] = 1.0;
                } else {
                    C[k][i] = sign * l1l2k * Binomial.choose(k + i - 1, i)
                            / Math.pow(lambda2 - lambda1, k + i);
                }
                D[k][i] = sign * l1l2k * Binomial.choose(k + i, i)
                        / Math.pow(lambda1 - lambda2, k + i + 1);

//                C[k][i] = C1(lambda1, lambda2, k, i);
//                D[k][i] = D1(lambda1, lambda2, k, i);
            }
        }
    }

//    private static double C1(double lambda1, double lambda2, int k, int index) {
//        if (k == 0 && index == 0) return 1.0;
//
//        return Math.pow(-1, index) * Math.pow(lambda1 * lambda2, k) * Binomial.choose(k + index - 1, index) /
//                Math.pow(lambda2 - lambda1, k + index);
//    }
//
//    private static double D1(double lambda1, double lambda2, int k, int index) {
//        return Math.pow(-1, index) * Math.pow(lambda1 * lambda2, k) * Binomial.choose(k + index, index) /
//                Math.pow(lambda1 - lambda2, k + index + 1);
//    }

    private void computeJumpProbabilities(double lambda1, double lambda2, double time,
                                          final double[][] C, final double[][] D, double[] jumpProbabilities) {

//        jumpProbabilities = new double[maxK];

        final double expLambda1Time = Math.exp(-lambda1 * time);
        final double expLambda2Time = Math.exp(-lambda2 * time);

        if (lambda1 == lambda2) {
            for (int m = 1; m < maxK / 2 + 1; ++m) {
                final int k = 2 * m;
                jumpProbabilities[k] = expLambda1Time * Math.pow(lambda1 * time, k)
                        / Math.exp(GammaFunction.lnGamma(k + 1));
            }
        } else {

            for (int k = 1; k < maxK / 2; ++k) {
                double sum = 0.0;
                double multiplicativeFactor = 1.0;
                for (int i = 1; i <= k + 1; ++i) {

                    if (i > 1) {
                        multiplicativeFactor *= time / (i - 1);
                    }

                    sum += C[k][k - i + 1]
                            * multiplicativeFactor
                            * expLambda1Time;

                    if (i <= k) {
                        sum += D[k][k - i]
                                * multiplicativeFactor
                                * expLambda2Time;
                    }
                }
                jumpProbabilities[2 * k] = sum;
            }
        }
    }

    public double[][] getC() { return C; }

    public double[][] getD() { return D; }

    public double[] getJumpProbabilities() { return jumpProbabilities; }

    public double computeCdf(double x, double time, int i, int j) {
        throw new RuntimeException("Not yet implemented");
    }

    public double computePdf(double x, double time, int i, int j) {

        if (i != 0 || j != 0) throw new RuntimeException("Not yet implemented");

        final double lambda0 = -Q[idx(0,0)];
        final double lambda1 = -Q[idx(1,1)];

        final boolean symmetric = (lambda0 == lambda1);

        if (!symmetric && C == null) {
            C = new double[(maxK / 2) + 1][(maxK / 2) + 1];
            D = new double[(maxK / 2) + 1][(maxK / 2) + 1];
            computeCDForJumpProbabilities(lambda0, lambda1, C, D);
        }

        if (jumpProbabilities == null) {
            jumpProbabilities = new double[maxK + 1];
//            computeJumpProbabilities(lambda0, lambda1, time, C, D, jumpProbabilities);   // Error: probs are function of time
        }
        computeJumpProbabilities(lambda0, lambda1, time, C, D, jumpProbabilities); // are function of time.
        // TODO Could cache computeJumpProbabilities(key = time) in HashMap

        if (symmetric) {
            // Single rate (symmetric)
            final double scale = 1.0 / lambda0;
            double sum = 0.0;

            // if time - x > 0, then there must have been at least k = 2 jumps
            for (int m = 1; m <= maxK / 2; ++m) {
                final int k = 2 * m;
                sum +=  jumpProbabilities[k] *
                        Math.exp(
                                + GammaDistribution.logPdf(x, m, scale)
                                + GammaDistribution.logPdf(time - x, m + 1, scale)
                                - GammaDistribution.logPdf(time, k + 1, scale)
                );
            }
            return sum;
        } else {
            // Two rate model
            double sum = 0.0;
            for (int m = 1; m <= maxK / 2; ++m) {
                final int k = 2 * m;
                sum += jumpProbabilities[k] *
                        GammaDistribution.pdf(x, m, 1.0 / lambda1) *
                        GammaDistribution.pdf(time - x, m + 1, 1.0 / lambda0) /
                        GeneralizedIntegerGammaDistribution.pdf(time, m, m + 1, lambda1, lambda0); // TODO Cache
            }
            // TODO Remove code duplication in if (symmetric) { } else { }
            return sum;
        }
    }

    public double[] computePdf(double x, double time) {
        throw new RuntimeException("Not yet implemented");
    }

    private double[][] squareMatrix(final double[] mat) {
        final int dim = 2;
        double[][] rtn = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                rtn[i][j] = mat[idx(i, j)];
            }
        }
        return rtn;
    }

    private int determineNumberOfSteps(double time, double lambda) {
//        final double tolerance = (1.0 - epsilon) / Math.exp(-lambda * time);
//        final double logTolerance = Math.log(1.0 - epsilon);

        int i = -1;
//        double sum = 0.0;
//        int factorialI = 1;
//
//        while (sum < tolerance) {
//            i++;
//            sum += Math.pow(lambda * time, i) / factorialI;
//            factorialI *= (i + 1); // + 1 because used on next iterate
//        }
//
//        int firstN = i;
//
//                            LogTricks.logSum()
//        i = -1;
        final double tolerance2 =
//                -epsilon;
                (1.0 - epsilon);
        double sum2 = 0.0;
//        double sum2 = Double.NEGATIVE_INFINITY;
        while (Math.abs(sum2 - tolerance2) > epsilon && sum2 < 1.0) {
//        while (sum2 < tolerance2) {
            i++;
            double logDensity = -lambda * time + i * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(i + 1);
            sum2 += Math.exp(logDensity);
//            sum2 = LogTricks.logSum(sum2, logDensity);
            if (DEBUG2) {
                System.err.println(sum2 + " " + tolerance2 + " " + Math.abs(sum2 - tolerance2) + " " + epsilon * 0.01);
//            if (i > 500) System.exit(-1);
            }
        }

//        System.err.println("First: " + firstN);
//        System.err.println("Second:" + i);
//        System.exit(-1);

        return i;
    }

//    private static int getKMax(int n) {
//        if (n % 2 == 0) {
//            return n / 2;
//        } else {
//            return (n + 1) / 2;
//        }
//    }

//    private static void hypgeoF11_k_2k(final double x, double[] result, double[] work) {
//        final int len = result.length;
//        if (len == 1 || x <= 0.0) return;
//
//        final int nmax = len - 1;
//        final int kmax = getKMax(nmax);
//        final double y = x / 2.0;
//        final double lx = Math.log(0.25 * x);
//        // TODO finish
//
//
//
//    }

//    private static void pclt(double t, double lambda1, double lambda2, double[] result) {
//
//        final int len = result.length;
//        final int nmax = len - 1;
//        final int kmax = getKMax(nmax);
//
//        double[] work = new double[kmax + 1];
//
//        final double x = t * (lambda2 - lambda1);
//        final double abs_x = Math.abs(x);
//        final double log_l1 = Math.log(lambda1), log_l2 = Math.log(lambda2), log_t = Math.log(t);
//        final double max_rate = Math.max(lambda1, lambda2);
//        double log_gamma = 0.0;
//
//        hypgeoF11_k_2k(abs_x, result, work);
//        // TODO finish
//
//    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Q: " + new Vector(Q) + "\n");
//        sb.append("r: " + new Vector(r) + "\n");
//        sb.append("lambda: " + lambda + "\n");
//        sb.append("N: " + getNfromC() + "\n");
//        sb.append("maxTime: " + maxTime + "\n");
//        sb.append("cprob at maxTime: " + new Vector(computeConditionalProbabilities(maxTime)) + "\n");
        return sb.toString();
    }

    private EigenDecomposition getEigenDecomposition() {
        if (eigenDecomposition == null) {
            eigenDecomposition = eigenSystem.decomposeMatrix(squareMatrix(Q));
        }
        return eigenDecomposition;
    }

    private EigenDecomposition eigenDecomposition;

    public double[] computeConditionalProbabilities(double distance) {

        double[] matrix = new double[4];
        eigenSystem.computeExponential(getEigenDecomposition(), distance, matrix);

        return matrix;
    }

    public double computeConditionalProbability(double distance, int i, int j) {
        return eigenSystem.computeExponential(getEigenDecomposition(), distance, i, j);
    }

    private final double[] Q;
    private final int maxK;

    private final double epsilon;

    private final EigenSystem eigenSystem;

    private double maxTime;

    private double[][] C = null;
    private double[][] D = null;
}
