/*
 * TwoStateOccupancyMarkovReward.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beagle.evomodel.substmodel.DefaultEigenSystem;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.EigenSystem;
import dr.math.Binomial;
import dr.math.GammaFunction;
import dr.math.distributions.GammaDistribution;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc Suchard
 * @author Christophe Fraiser
 * @author Andrew Rambaut
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

        if (Q[idx(0, 0)] != Q[idx(1, 1)]) {
            throw new IllegalArgumentException("Only currently implemented for equal rates models");
        }
    }

    private int idx(int i, int j) {
        return i * 2 + j; // row-major
    }

    private double[] jumpProbabilities = null;

    private double getAki(double lambda1, double lambda2, int k, int i) {
//        double logA = Binomial.logChoose(k + i - 1, i) +
//                k * Math.log(lambda1) + (k - 1) * Math.log(lambda2)
//        int sign = 1;
//        return sign * Math.exp(logA);
        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
                * Math.pow(lambda2, k - 1) / Math.pow(lambda2 - lambda1, k + i);
    }

    private double getBki(double lambda1, double lambda2, int k, int i) {
        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
                * Math.pow(lambda2, k - 1) / Math.pow(lambda1 - lambda2, k + i);
    }

    private double getCki(double lambda1, double lambda2, int k, int i) {
        return Binomial.choose(k + i - 1, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
                * Math.pow(lambda2, k) / Math.pow(lambda2 - lambda1, k + i);
    }

    private double getDki(double lambda1, double lambda2, int k, int i) {
        return Binomial.choose(k + i, i) * Math.pow(-1, i) * Math.pow(lambda1, k)
                * Math.pow(lambda2, k) / Math.pow(lambda1 - lambda2, k + i + 1);
    }


    private void computeJumpProbabilities(double lambda1, double lambda2, double[] jumpProbabilities) {
        jumpProbabilities = new double[maxK];

//        if (lambda1 == lambda2) { // Poisson process
//
//        } else {
//            for (int k = 0; k < maxK / 2; ++k) {
//                double sumC = 0.0;
//                double sumD = 0.0;
//                for (int i = 1; i <= k + 1; ++i) {
//                    sumC += getCki(lambda1, lambda2, k, k - i + 1) * Math
//                }
//                jumpProbabilities[2 * k] = sumC + sumD;
//            }
//        }
    }

    public double computeCdf(double x, double time, int i, int j) {
        return 0.0;
    }

    public double computePdf(double x, double time, int i, int j) {

//        if (jumpProbabilities == null) {
//            computeJumpProbabilities(Q[idx(0,0)], Q[idx(1,1)], jumpProbabilities);
//        }

        final double lambda = -Q[idx(0, 0)];
        final double rate = 1.0 / lambda;
        final double logLambdaTime = Math.log(lambda) + Math.log(time);

        final double time2 = time - x;

        final double multiplier = Math.exp(-lambda * time);

        double sum = 0.0;

        // if time - x > 0, then there must have been at least k = 2 jumps
        for (int m = 1; m <= maxK / 2; ++m) {
            final int k = 2 * m;
            sum += Math.exp(k * logLambdaTime - GammaFunction.lnGamma(k + 1)
                    + GammaDistribution.logPdf(x, m, rate)
                    + GammaDistribution.logPdf(time - x, m + 1, rate)
                    - GammaDistribution.logPdf(time, k + 1, rate)
            );
        }
        return multiplier * sum;
    }

    public double[] computePdf(double x, double time) {
//        return computePdf(new double[]{x}, time)[0];
        return null;
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
}
