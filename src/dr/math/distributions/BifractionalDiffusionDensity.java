/*
 * BifractionalDiffusionDensity.java
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
import cern.jet.stat.Gamma;


/**
 * @author Marc Suchard
 *         <p/>
 *         **
 *         Most of the following functions come from:
 *         <p/>
 *         R Gorenflo, A Iskenderov and Y Luchko (2000) Mapping between solutions of fractional diffusion-wave equations.
 *         Fractional Calculus and Applied Analysis, 3, 75 - 86
 *         <p/>
 *         Also see:
 *         F. Mainardi and G. Pagnini (2003) The Wright functions as solutions of the time-fractional diffusion equation.
 *         Applied Mathematics and Computation, 141, pages?
 *
 *
 * Need to read: http://www.hindawi.com/journals/ijde/2010/104505.html and http://www1.tfh-berlin.de/~luchko/papers/Luchko_111.pdf
 */

public class BifractionalDiffusionDensity implements Distribution {

    public BifractionalDiffusionDensity(double v, double alpha, double beta) {
        this.v = v;
        this.alpha = alpha;
        this.beta = beta;
        coefficients = constructBifractionalDiffusionCoefficients(alpha, beta);
    }

    public BifractionalDiffusionDensity(double alpha, double beta) {
        this(1.0, alpha, beta);
    }

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return pdf(x, v);
    }

    public double pdf(double x, double v) {
        return pdf(x, v, alpha, beta, coefficients);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return logPdf(x, v, alpha, beta, coefficients);
    }

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * @return a probability density function representing this distribution
     */
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    public static double logPdf(double x, double v, double alpha, double beta) {
        return Math.log(pdf(x, v, alpha, beta));
    }

    public static double logPdf(double x, double v, double alpha, double beta, double[][][] coefficients) {
        return Math.log(pdf(x, v, alpha, beta, coefficients));
    }

//    /*
//     * Taken from:  Saichev AI and Zaslavsky GM (1997) Fractional kinetic equations: solutions and applications.
//     *              Chaos, 7, 753-764
//     * @param x evaluation point
//     * @param t evaluation time
//     * @param alpha coefficient
//     * @param beta coefficient
//     * @return probability density
//     */
//
//    public static double logPdf(double x, double t, double alpha, double beta) {
//        final double mu = beta / alpha;
//        final double absX = Math.abs(x);
//        final double absY = absX / Math.pow(t,mu);
//        double density = 0;
//        double incr = Double.MAX_VALUE;
//        int m = 1; // 0th term = 0 \propto cos(pi/2)
//        int sign = -1;
//
//        while (Math.abs(incr) > eps && m < maxK) {
//            incr =  sign / Math.pow(absY, m * alpha)
//                         * Gamma.gamma(m * alpha + 1)
//                         / Gamma.gamma(m * beta + 1)
//                         * Math.cos(halfPI * (m * alpha + 1));
//            density += incr;
//            sign *= -1;
//            m++;
//        }
//
//        return Math.log(density / (Math.PI * absX));
//    }


    static class SignedDouble {
        double x;
        boolean positive;

        SignedDouble(double x, boolean positive) {
            this.x = x;
            this.positive = positive;
        }
    }

    public static SignedDouble logGamma(double z) {

        // To extend the gamma function to negative (non-integer) numbers, apply the relationship
        // \Gamma(z) = \frac{ \Gamma(z+n) }{ z(z+1)\cdots(z+n-1),
        // by choosing n such that z+n is positive
        if (z > 0) {
            return new SignedDouble(Gamma.logGamma(z), true);
        }
        int n = ((int) -z);
        if (z + n == 0.0) {
            return new SignedDouble(Double.NaN, true); // Zero and the negative integers are diverging poles
        }
        n++;
        boolean positive = (n % 2 == 0);
        return new SignedDouble(Gamma.logGamma(z + n) - Gamma.logGamma(-z + 1) + Gamma.logGamma(-z - n + 1), positive);
    }


    public static double gamma(double z) {

        // To extend the gamma function to negative (non-integer) numbers, apply the relationship
        // \Gamma(z) = \frac{ \Gamma(z+n) }{ z(z+1)\cdots(z+n-1),
        // by choosing n such that z+n is positive
        if (z > 0.0) {
            return Gamma.gamma(z);
        }
        int n = ((int) -z);
        if (z+n == 0.0) {
            return Double.NaN; // Zero and the negative integers are diverging poles
        }
        n++;
        boolean positive = (n % 2 == 0);
        double result = Gamma.gamma(z + n) / Gamma.gamma(-z + 1) * Gamma.gamma(-z - n + 1);
        if (!positive) result *= -1;
        return result;
    }


// The following comments out functions may be faster than using the generalized Wright function.
// Keep for comparison
//
//    public static double infiniteSumAlphaGreaterThanBeta1(double z, double alpha, double beta) {
//        double sum = 0;
//        int k = 0;
//        boolean isPositive = true;
//
//        double incr = Double.MAX_VALUE;
//        while (//(Math.abs(incr) > 1E-20) &&
//                (k < maxK)) {
//
//            double x1, x2, x3;
//            x1 = gamma(0.5 - 0.5 * alpha - 0.5 * alpha * k);
//            if (!Double.isNaN(x1)) {
//                incr = x1;
//            } else {
//                System.err.println("Big problem!");
//                System.exit(-1);
//            }
//            x2 = gamma(1.0 - beta - beta * k);
//            if (!Double.isNaN(x2)) {
//                incr /= x2;
//            } else {
//                 incr = 0;
//            }
//            x3 = gamma(0.5 * alpha + 0.5 * alpha * k);
//            if (!Double.isNaN(x3)) {
//                incr /= x3;
//            } else {
//                incr = 0;
//            }
//            incr *= Math.pow(z, k);
//            if (isPositive) {
//                sum += incr;
//            } else {
//                sum -= incr;
//            }
//            isPositive = !isPositive;
//            k++;
//        }
//        return sum;
//    }
//
//    public static double infiniteSumAlphaGreaterThanBeta2(double z, double alpha, double beta) {
//        double sum = 0;
//        int m = 0;
//        boolean isPositive = true;
//
//        double incr = Double.MAX_VALUE;
//        while (// (Math.abs(incr) > 1E-20) &&
//                (m < maxK)) {
//
//            double x1, x2, x3, x4;
//
//            x1 = gamma(1.0 / alpha + 2.0 / alpha * m);
//            if (!Double.isNaN(x1)) {
//                incr = x1;
//            } else {
//                System.err.println("Big problem!");
//                System.exit(-1);
//            }
//            x2 = gamma(1.0 - 1.0 / alpha - 2.0 / alpha * m);
//            if (!Double.isNaN(x2)) {
//                incr *= x2;
//            } else {
//                System.err.println("Big problem!");
//                System.exit(-1);
//            }
//            x3 = gamma(0.5 + m);
//            if (!Double.isNaN(x3)) {
//                incr /= x3;
//            } else {
//                incr = 0;
//            }
//            x4 = gamma(1.0 - beta / alpha - 2 * beta / alpha * m);
//            if (!Double.isNaN(x4)) {
//                incr /= x4;
//            } else {
//                incr = 0;
//            }
//            incr /= gamma(m + 1);
//            incr *= Math.pow(z, m);
//            if (isPositive) {
//                sum += incr;
//            } else {
//                sum -= incr;
//            }
//            isPositive = !isPositive;
//            m++;
//
//        }
//        return sum;
//    }

    /*
     * Taken from Equation (17) in Gorenflo, Iskenderov and Luchko (2000)
     */
    private static double evaluateGreensFunctionAtZero(double t, double alpha, double beta) {
        if (beta == 1) {
            final double oneOverAlpha = 1.0 / alpha;
            return gamma(oneOverAlpha) / (Math.PI * alpha * Math.pow(t, oneOverAlpha));
        } else {
            final double betaOverAlpha = beta / alpha;
            return 1.0 / (alpha * Math.pow(t, betaOverAlpha) * Math.sin(Math.PI / alpha) * gamma(1.0 - betaOverAlpha));
        }
    }

    /*
     * Taken from Equation (20) in Gorenflo, Iskenderov and Luchko (2000)
     */
    private static double evaluateGreensFunctionAlphaEqualsBeta(double x, double t, double alpha) {

        final double absX = Math.abs(x);
        final double twoAlpha = 2.0 * alpha;
        final double tPowAlpha = Math.pow(t, alpha);
        final double piHalfAlpha = 0.5 * Math.PI * alpha;

        double green = Math.pow(absX, alpha - 1.0) * tPowAlpha * Math.sin(piHalfAlpha) /
                       (Math.pow(t, twoAlpha) + 2 * Math.pow(absX, alpha) * tPowAlpha * Math.cos(piHalfAlpha)
                                              + Math.pow(absX, twoAlpha));
        return oneOverPi * green;
    }


    /*
     * Taken from Equation (18) in Gorenflo, Iskenderov and Luchko (2000)
     */
    private static double evaluateGreensFunctionBetaGreaterThanAlpha(double x, double t, double alpha, double beta,
                                                                     double[][][] coefficients) {
        double z = Math.pow(2.0, alpha) * Math.pow(t, beta) / Math.pow(Math.abs(x), alpha);
        return oneOverSqrtPi / Math.sqrt(Math.abs(x)) * generalizedWrightFunction(-z, coefficients[0], coefficients[1]);
    }


    /*
     * Taken from Equation (19) in Gorenflo, Iskenderov and Luchko (2000)
     */
    private static double evaluateGreensFunctionAlphaGreaterThanBeta(double x, double t, double alpha, double beta,
                                                                     double[][][] coefficients) {

         double z1 = Math.pow(Math.abs(x),alpha) / (Math.pow(2,alpha) * Math.pow(t, beta));
         double z2 = x * x / (4.0 * Math.pow(t, 2 * beta / alpha));

         double green1 = oneOverSqrtPi * Math.pow(Math.abs(x), alpha - 1.0) /
                 (Math.pow(2,alpha) * Math.pow(t, beta)) *
                 generalizedWrightFunction(-z1, coefficients[2], coefficients[3]);
        
         double green2 = oneOverSqrtPi /  (alpha * Math.pow(t, beta/alpha)) *
                 generalizedWrightFunction(-z2, coefficients[4], coefficients[5]);
   
         return green1 + green2;
    }

    public static double pdf(double x, double v, double alpha, double beta) {
        return pdf(x, v, alpha, beta, null);
    }

    public static double pdf(double x, double v, double alpha, double beta, double[][][] coefficients) {
        final double t = 0.5 * v;
        if (x == 0) {
            return evaluateGreensFunctionAtZero(t, alpha, beta);
        }
        if (alpha == beta) {
            return evaluateGreensFunctionAlphaEqualsBeta(x, t, alpha);
        }
        if (coefficients == null) {
            coefficients = constructBifractionalDiffusionCoefficients(alpha, beta);
        }
        if (alpha > beta) {
            return evaluateGreensFunctionAlphaGreaterThanBeta(x, t, alpha, beta, coefficients);
        } else {
            return evaluateGreensFunctionBetaGreaterThanAlpha(x, t, alpha, beta, coefficients);
        }
    }

    /*
     * Helper function to construct coefficients for generalized Wright function
     */
    public static double[][][] constructBifractionalDiffusionCoefficients(double alpha, double beta) {
        double[][][] coefficients = new double[6][][];

        // coefficients[0:1][][] : Greens function beta > alpha
        coefficients[0] = new double[][]{{0.5, alpha / 2.0}, {1.0, 1.0}};
        coefficients[1] = new double[][]{{1.0, beta}, {0.0, -alpha / 2.0}};

        // coefficients[2:3][][] : Greens function #1 alpha > beta
        coefficients[2] = new double[][]{{0.5 - alpha / 2.0, -alpha / 2.0}, {1.0, 1.0}};
        coefficients[3] = new double[][]{{1.0 - beta, -beta}, {alpha / 2.0, alpha / 2.0}};

        // coefficients[4:5][][] : Greens function #2 alpha > beta
        coefficients[4] = new double[][]{{1.0 / alpha, 2.0 / alpha}, {1.0 - 1.0 / alpha, -2.0 / alpha}};
        coefficients[5] = new double[][]{{0.5, 1.0}, {1.0 - beta / alpha, -2.0 * beta / alpha}};

        return coefficients;
    }

    public static double generalizedWrightFunction(double z, double[][] aAp, double[][] bBq) {
        final int p = aAp.length;
        final int q = bBq.length;
        double sum = 0.0;
        double incr;
        double zPowK = 1.0;
        int k = 0;

        while (// incr > eps &&
                k < maxK) {
            incr = 1;
            for (int i = 0; i < p; i++) {
                final double[] aAi = aAp[i];

                double x = gamma(aAi[0] + aAi[1] * k); // TODO Precompute these factors
                if (!Double.isNaN(x)) {
                    incr *= x;
                } else {
                    incr = Double.NaN;
                }
            }
            for (int i = 0; i < q; i++) {
                final double[] bBi = bBq[i];
                double x = gamma(bBi[0] + bBi[1] * k); // TODO Precompute these factors
                if (!Double.isNaN(x)) {
                    incr /= x;
                } else {
                    incr = 0.0;
                }
            }
            incr /= gamma(k+1); // k! TODO Precompute these factors
            incr *= zPowK;
            sum += incr;

            // Get ready for next loop
            zPowK *= z;
            k++;
        }
        return sum;
    }

    public static void main(String[] arg) {

        double alpha = 2.0;
        double beta = 0.8;
        double z1 = -2.34;
        SignedDouble result = logGamma(z1);
        System.err.println("logGamma("+z1+") = "+result.x+" "+(result.positive ? "(+)" : "(-)"));
        System.err.println("gamma("+z1+") = "+ gamma(z1));
        System.err.println("gamma(-2.0) = "+gamma(-2.0));
        System.err.println("");

        double var = 4.0;
        double t = 0.5 * var;
        double x = 1.0;

        double[][][] coefficients = constructBifractionalDiffusionCoefficients(alpha, beta);
        System.err.println("p(x = "+x+", v = "+var+") = " + evaluateGreensFunctionAlphaGreaterThanBeta(x, t, alpha,
                beta, coefficients));


        alpha = 0.7;
        beta = 1.4;
        coefficients = constructBifractionalDiffusionCoefficients(alpha, beta);
//        System.err.println("p(x = "+x+", v = "+var+") = " + evaluateGreensFunctionBetaGreaterThanAlpha(x, t, alpha, beta));
        System.err.println("p(x = "+x+", v = "+var+") = " + evaluateGreensFunctionBetaGreaterThanAlpha(x, t, alpha,
                beta, coefficients));        



    }

    private static double oneOverSqrtPi = 1.0 / Math.sqrt(Math.PI);
    private static double oneOverPi = 1.0 / Math.PI;
    public static final int maxK = 50;

    private double alpha;
    private double beta;
    private double v;
    private double[][][] coefficients;
}