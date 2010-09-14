/*
 * InverseGaussianDistribution.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
import dr.math.interfaces.OneVariableFunction;
import dr.math.iterations.BisectionZeroFinder;
import dr.math.iterations.NewtonZeroFinder;

/**
 * normal distribution (pdf, cdf, quantile)
 *
 * @author Wai Lok Sibon Li
 * @version $Id: InverseGaussianDistribution.java,v 1.7 2008/04/24 20:26:01 rambaut Exp $
 *
 * Reading: Chhikara, R. S., and Folks, J. Leroy, (1989). The
 * inverse Gaussian distribution: Theory, methodology, and applications. Marcel Dekker, New York.
 */
public class InverseGaussianDistribution implements Distribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     *
     * @param mean  mean
     * @param shape shape
     */
    public InverseGaussianDistribution(double mean, double shape) {
        this.m = mean;
        this.shape = shape;
        this.sd = calculateSD(mean, shape);
    }

    public double getMean() {
        return m;
    }

    public void setMean(double value) {
        m = value;
    }

    public double getShape() {
        return shape;
    }

    public void setShape(double value) {
        shape = value;
        sd = calculateSD(m, shape);
    }

    public static double calculateSD(double mean, double shape) {
        return Math.sqrt((mean * mean * mean) / shape);
    }
    //public double getSD() {
    //return sd;
    //}

    //public void setSD(double value) {
    //sd = value;
    //}

    public double pdf(double x) {
        return pdf(x, m, shape);
    }

    public double logPdf(double x) {
        return logPdf(x, m, shape);
    }

    public double cdf(double x) {
        return cdf(x, m, shape);
    }

    public double quantile(double y) {
        return quantile(y, m, shape);
    }

    public double mean() {
        return mean(m, shape);
    }

    public double variance() {
        return variance(m, shape);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
            //return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * probability density function
     *
     * @param x     argument
     * @param m     mean
     * @param shape shape parameter
     * @return pdf at x
     */
    public static double pdf(double x, double m, double shape) {
        double a = Math.sqrt(shape / (2.0 * Math.PI * x * x * x));
        double b = ((-shape) * (x - m) * (x - m)) / (2.0 * m * m * x);
        return a * Math.exp(b);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x     argument
     * @param m     mean
     * @param shape shape parameter
     * @return log pdf at x
     */
    public static double logPdf(double x, double m, double shape) {
        double a = Math.sqrt(shape / (2.0 * Math.PI * x * x * x));
        double b = ((-shape) * (x - m) * (x - m)) / (2.0 * m * m * x);
        return Math.log(a) + b;
    }

    /**
     * cumulative density function
     *
     * @param x     argument
     * @param m     mean
     * @param shape shape parameter
     * @return cdf at x
     */
    public static double cdf(double x, double m, double shape) {
        if (x <= 0 || m <= 0 || shape <= 0) {
            return Double.NaN;
        }
        /* Taken from R code, package SuppDists */
        double a = Math.sqrt(shape / x);
        double b = x / m;
        //double p1 = NormalDistribution.cdf(a*(b - 1.0),0.0,1.0);
        double p1 = NormalDistribution.cdf(a * (b - 1.0), 0.0, 1.0, false);
        //double p2 = NormalDistribution.cdf(-a*(b + 1.0),0.0,1.0);
        double p2 = NormalDistribution.cdf(-a * (b + 1.0), 0.0, 1.0, false);
        if (p2 == 0.0) {
            return p1;
        }
        else {
            double c=2.0 * shape / m;
            if (c>=0x1.fffffffffffffP+1023) {// Double.MAX_EXPONENT is Java 1.6 feature
                return Double.POSITIVE_INFINITY;
            }
            return p1 + Math.exp(c) * p2;
        }

        /* Another implementation of the inverse Gaussian cdf that doesn't use the Normal distribution function
         * Is not as accurate (due to error function issues)
         */
//        double a = Math.sqrt(shape / (2.0 * x)) * ((x / m) - 1);
//        double b = (1.0 + ErrorFunction.erf(a));
//        double c = Math.sqrt(shape / (2.0 * x)) * ((x / m) + 1);
//        double d = ((2.0 * shape) / m) + Math.log(1 - ErrorFunction.erf(c));
//        return 0.5*b + 0.5*Math.exp(d);
    }

    /**
     * quantiles (=inverse cumulative density function)
     * <p/>
     *
     * Same implementation as SuppleDists in R.
     *
     * Using Whitmore and Yalovsky for an initial guess. Works well for
	 * large t=lambda/mu > 2 perhaps
	 * Whitmore, G.A. and Yalovsky, M. (1978). A normalizing logarithmic
	 * transformation for inverse Gaussian random variables,
	 * Technometrics 20-2, 207-208
	 * For small t, with x<0.5 mu, use gamma approx to 1/x -- alpha=1/2 and beta =2/lambda and 1-p
	 * When x>0.5mu, approx x with gamma for p, and exponentiate -- don't know why this works.
     *
     * There are cases  which even this method produces inaccurate results (e.g. when shape = 351). Therefore,
     * we have a catch that will determine whether or not the result is accurate enough and if not,
     * then a zerofinder will be used to find a more accurate approximation. 
     *
     * @param z     argument
     * @param m     mean
     * @param shape shape parameter
     * @return icdf at z
     */
    public static double quantile(double z, double m, double shape) {
        if(z < 0.01 || z > 0.99) {
            throw new RuntimeException("Quantile is too low/high to calculate (numerical estimation for extreme values is incomplete");
        }

        /* Approximation method used by Mudholkar GS, Natarajan R (1999)
         * Approximations for the inverse gaussian probabilities and percentiles.
         * Communications in Statistics - Simulation and Computation 28: 1051 - 1071.
         */
        double initialGuess;
        if (shape / m > 2.0) {
            initialGuess=(NormalDistribution.quantile(z,0.0,1.0)-0.5*Math.sqrt(m/shape))/Math.sqrt(shape/m);
            initialGuess=m*Math.exp(initialGuess);
        }
        else {
            initialGuess=shape/(GammaDistribution.quantile(1.0-z,0.5,1.0)*2.0);
            if (initialGuess > m / 2.0) {		// too large for the gamma approx
                initialGuess=m*Math.exp(GammaDistribution.quantile(z,0.5,1.0)*0.1);  // this seems to work for the upper tail ???
            }
        }
//        double phi = shape / m;
//        if(phi>50.0) {
            // Use Normal Distribution
//            initialGuess = (NormalDistribution.quantile(z, m,Math.sqrt(m*m*m/shape)));//-0.5*Math.sqrt(m/shape))/Math.sqrt(m*m*m/shape);
//        }

        final InverseGaussianDistribution f = new InverseGaussianDistribution(m, shape);
        final double y = z;
        NewtonZeroFinder zeroFinder = new NewtonZeroFinder(new OneVariableFunction() {
            public double value (double x) {
                return f.cdf(x) - y;
            }
        }, initialGuess);
        zeroFinder.evaluate();

        if(Double.isNaN(zeroFinder.getResult()) || zeroFinder.getPrecision() > 0.000005) {
            zeroFinder = new NewtonZeroFinder(new OneVariableFunction() {
                public double value (double x) {
                    return f.cdf(x) - y;
                }
            }, initialGuess);
            zeroFinder.initializeIterations();
            int i;
            double previousPrecision = 0.0, previousResult = Double.NaN;
            double max = 10000.0, min = 0.00001;
            for(i=0; i < 50; i++) {
                zeroFinder.evaluateIteration();
                double precision = f.cdf(zeroFinder.getResult()) - z;
                if((previousPrecision > 0 && precision < 0) || (previousPrecision < 0 && precision > 0))  {
                    max = Math.max(previousResult, zeroFinder.getResult());
                    min = Math.min(previousResult, zeroFinder.getResult());
                    max = Math.min(10000.0, max);
                    break;
                }

                previousPrecision = precision;
                previousResult = zeroFinder.getResult();

            }
            return calculateZeroFinderApproximation(z, m, shape, min, max, initialGuess);
        }
        return zeroFinder.getResult();
    }

    /** Finds the approximation of the inverse Gaussian quantile using a zero finder
     * until it converges
     *
     * @param z            quantile
     * @param m            mean
     * @param shape        shape
     * @param min          min search value
     * @param max          max search value
     * @param initialGuess first guess of the quantile
     * @return estimated x value at quantile z
     */
    private static double calculateZeroFinderApproximation(double z, double m, double shape, double min, double max, double initialGuess) {
        final InverseGaussianDistribution f = new InverseGaussianDistribution(m, shape);
        final double y = z;
        BisectionZeroFinder bisectionZeroFinder = new BisectionZeroFinder(new OneVariableFunction() {
            public double value(double x) {
                return f.cdf(x) - y;
            }
        }, min, max);
        bisectionZeroFinder.setInitialValue(initialGuess);
        bisectionZeroFinder.initializeIterations();

        double bestValue = Double.NaN; /* I found that the converged value is not necesssarily the best */
        double bestPrecision = 10;
        double precision = 10;
        double previousPrecision = 10;
        int count = 0;
        while(precision > 0.001 &&  count < 10) {
            bisectionZeroFinder.evaluateIteration();
            precision = Math.abs(f.cdf(bisectionZeroFinder.getResult()) - z);
            if(precision < bestPrecision) {
                bestPrecision = precision;
                bestValue = bisectionZeroFinder.getResult();
            }
            else if(previousPrecision == precision) {
                count++;
            }
            previousPrecision = precision;
        }
        bisectionZeroFinder.finalizeIterations();
        /* Turns out the final answer is not necessarily the most accurate */
        //return bisectionZeroFinder.getResult();
        return bestValue;
    }


    /** Calculates the gamma approximation of the quantile estimate of Inverse Gaussian
     * Shifted Gamma
     * (see Mudholkar GS, Natarajan R (1999))
     * UNUSED METHOD
     *
     * @param z     quantile
     * @param m     mean
     * @param shape shape
     * @return approximation of x
     */
    private static double calculateShiftedGammaApproximation(double z, double m, double shape) {
        double a = (3 * m * m) / (4 * shape);
        double b = (m / 3);
        double nu  = (8 * shape) / (9 * m);
        return a * ChiSquareDistribution.quantile(z, nu) + b;
    }

    /** Calculates the gamma approximation of the quantile estimate of Inverse Gaussian
     * Shifted Gamma adapted to reciprocal Inverse Gaussian
     * (see Mudholkar GS, Natarajan R (1999))
     * UNUSED METHOD
     *
     * @param z     quantile
     * @param m     mean
     * @param shape shape
     * @return approximation of x
     */
    private static double calculateShiftedGammaApproximationWithRIG(double z, double m, double shape) {
        double a = (3 * shape + 8 * m)/(4 * shape * (shape + 2 * m));
        double b = (shape + 3 * m)/(m * (3 * shape + 8 * m));
        double nu = (8 * Math.pow((shape + 2 * m), 3)) / (m * Math.pow((8 * m + 3 * shape), 2));
        double y_hat = a * ChiSquareDistribution.quantile(z, nu) + b;
        return 1 / y_hat;
    }

    /** Finds the approximation of the inverse Gaussian quantile using a zero finder
     * given a set number of maximum iterations
     * UNUSED METHOD
     *
     * @param z             quantile
     * @param m             mean
     * @param shape         shape
     * @param numIterations number of iterations used to approximate value
     * @param min           min search value
     * @param max           max search value
     * @return estimated x value at quantile z
     */
    private static double calculateZeroFinderApproximation(double z, double m, double shape, int numIterations, double min, double max) {
        final InverseGaussianDistribution f = new InverseGaussianDistribution(m, shape);
        final double y = z;
        BisectionZeroFinder bisectionZeroFinder = new BisectionZeroFinder(new OneVariableFunction() {
            public double value(double x) {
                return f.cdf(x) - y;
            }
        }, min, max);
        //}, 0.0001, 100000);

        bisectionZeroFinder.setMaximumIterations(numIterations);
        bisectionZeroFinder.evaluate();
        return bisectionZeroFinder.getResult();
    }
    
    /**
     * mean
     *
     * @param m     mean
     * @param shape shape parameter
     * @return mean
     */
    public static double mean(double m, double shape) {
        return m;
    }

    /**
     * variance
     *
     * @param m     mean
     * @param shape shape parameter
     * @return variance
     */
    public static double variance(double m, double shape) {
        double sd = calculateSD(m, shape);
        return sd * sd;
    }

    // Private

    protected double m, sd, shape;

}
