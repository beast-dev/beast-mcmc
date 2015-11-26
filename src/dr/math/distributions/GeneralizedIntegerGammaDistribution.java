/*
 * Hyperexponential.java
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

package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * GeneralizedIntegerGammaDistribution distribution.
 * @author Marc A. Suchard
 */
public class GeneralizedIntegerGammaDistribution implements Distribution {

    private int shape1, shape2;
    private double rate1, rate2;

    private double[] A = null;
    private double[] B = null;

    public GeneralizedIntegerGammaDistribution(int shape1, int shape2, double rate1, double rate2) {
        this.shape1 = shape1;
        this.shape2 = shape2;
        this.rate1 = rate1;
        this.rate2 = rate2;
    }

    @Override
    public double logPdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    public double generatingFunction(double s) {
        return Math.pow(rate1 / (rate1 + s), shape1) * Math.pow(rate2 / (rate2 + s), shape2);
    }

    //        http://www.ism.ac.jp/editsec/aism/pdf/034_3_0591.pdf
    public double generatingFunctionPartialFraction(double s) {
        if (A == null) {
            computeCoefficients();
        }
        double sum = 0.0;

        for (int i = 1; i <= shape1; ++i) {
            sum += A[i] / Math.pow(rate1 + s, i);
        }

        for (int i = 1; i <= shape2; ++i) {
            sum += B[i] / Math.pow(rate2 + s, i);
        }

        return sum;
    }

    private void computeCoefficients() {
        A = new double[shape1 + 1];
        B = new double[shape2 + 1];

        final double lambdaFactor = Math.pow(rate1, shape1) * Math.pow(rate2, shape2);

        int sign = 1;
        double factorial = 1.0;
        for (int i = 1; i <= shape1; ++i) {
            if (i > 1 && (shape2 + i - 2) > 1) {
                factorial *= shape2 + i - 2;
                factorial /= i - 1;
            }
            //                System.err.println("A[" + (shape1 - i + 1) + "]: " + factorial);
            A[shape1 - i + 1] = factorial * sign * lambdaFactor / Math.pow(rate2 - rate1, shape2 + i - 1);
            sign *= -1;
        }

        sign = 1;
        factorial = 1.0;
        for (int i = 1; i <= shape2; ++i) {
            if (i > 1 && (shape1 + i - 2) > 1) {
                factorial *= shape1 + i - 2;
                factorial /= i - 1;
            }

            //                System.err.println("B[" + (shape2 - i + 1) + "]: " + factorial);
            B[shape2 - i + 1] = factorial * sign * lambdaFactor / Math.pow(rate1 - rate2, shape1 + i - 1);
            sign *= -1;
        }
    }

    public double pdf(double x) {
        if (A == null) {
            computeCoefficients();
        }

        final double expRate1X = Math.exp(-rate1 * x);
        final double expRate2X = Math.exp(-rate2 * x);

        double sum = 0.0;
        double power = 1.0;
        int factorial = 1;

        for (int i = 1; i <= shape1; ++i) {
            sum += A[i] * power * expRate1X / factorial;

            power *= x; // x^{i - 1}
            factorial *= i; // (i - 1)!
        }

        power = 1.0;
        factorial = 1;
        for (int i = 1; i <= shape2; ++i) {
            sum += B[i] * power * expRate2X / factorial;

            power *= x; // x^{i - 1}
            factorial *= i; // (i - 1)!
        }

        return sum;
    }

    public static double pdf(double x, int shape1, int shape2, double rate1, double rate2) {
        return new GeneralizedIntegerGammaDistribution(shape1, shape2, rate1, rate2).pdf(x);
    }


// https://en.wikipedia.org/wiki/Generalized_integer_gamma_distribution

// https://en.wikipedia.org/wiki/Generalized_integer_gamma_distribution

//        http://arxiv.org/pdf/math/0408189v1.pdf

//        http://www.math.utep.edu/Faculty/moschopoulos/Publications/1985-The_Distribution_of_the_Sum_of_Independent_Gamma_Random_Variables.pdf

//          partial fractions, inverse Fourier transform
//        http://stats.stackexchange.com/questions/72479/general-sum-of-gamma-distributions


}

