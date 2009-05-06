package dr.math.distributions;

import dr.math.UnivariateFunction;
import dr.math.ErrorFunction;

/*
 * NormalDistribution.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
/**
 * normal distribution (pdf, cdf, quantile)
 *
 * @author Wai Lok Sibon Li
 * @version $Id: InverseGaussianDistribution.java,v 1.7 2008/04/24 20:26:01 rambaut Exp $
 */
public class InverseGaussianDistribution implements Distribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     * @param mean  mean
     * @param shape  shape
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
        return Math.sqrt((mean*mean*mean)/shape);
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
     * @param x  argument
     * @param m  mean
     * @param shape  shape parameter
     * @return pdf at x
     */
    public static double pdf(double x, double m, double shape) {
        // For normal
        //double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * sd);
        //double b = -(x - m) * (x - m) / (2.0 * sd * sd);

        double a = Math.sqrt(shape/(2.0 * Math.PI * x * x * x));
        double b = ((-shape) * (x - m) * (x - m))/(2.0 * m * m * x);

        return a * Math.exp(b);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x  argument
     * @param m  mean
     * @param shape  shape parameter
     * @return log pdf at x
     */
    public static double logPdf(double x, double m, double shape) {
        // For normal
        //double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * sd);
        //double b = -(x - m) * (x - m) / (2.0 * sd * sd);

        double a = Math.sqrt(shape/(2.0 * Math.PI * x * x * x));
        double b = ((-shape) * (x - m) * (x - m))/(2.0 * m * m * x);

        return Math.log(a) + b;
    }

    // ============================================HAVE NOT COMPLETED FROM HERE DOWNWARDS=========================================================================================

    /**
     * cumulative density function
     *
     * @param x  argument
     * @param m  mean
     * @param shape  shape parameter
     * @return cdf at x
     */
    public static double cdf(double x, double m, double shape) {
        // For normal
        //double a = (x - m) / (Math.sqrt(2.0) * sd);

        //return 0.5 * (1.0 + ErrorFunction.erf(a));
        double a = Math.sqrt(shape / (2.0 * x)) * ((x / m) - 1);
        double b = (1.0 + ErrorFunction.erf(a));
        double c = Math.sqrt(shape / (2.0 * x)) * ((x / m) + 1);
        double d = Math.exp((2.0 * shape)/m) * (1 - ErrorFunction.erf(c));

        return 0.5*b + 0.5*d;
    }

    /**
     * quantiles (=inverse cumulative density function)
     *
     * CURRENTLY NOT IMPLEMENTED PROPERLY. Can be implemented later using a Zero finder function
     * (See JUnit test for LogNormal distribution for zero finder). Alternatively find out
     * how they do it with SuppleDists in R (download source code and open the C function which
     * contains the implementation. Reading: Chhikara, R. S., and Folks, J. Leroy, (1989). The
     * inverse Gaussian distribution: Theory, methodology, and applications. Marcel Dekker, New York.
     *
     * @param z  argument
     * @param m  mean
     * @param shape  shape parameter
     * @return icdf at z
     */
    public static double quantile(double z, double m, double shape) {
        double x=0;

        throw new RuntimeException("Quantile function for Inverse Gaussian Distribution is not yet implemented");

        //System.out.println(cdf(x, m, shape));

        // For normal
        //return m + Math.sqrt(2.0) * sd * ErrorFunction.inverseErf(2.0 * z - 1.0);
        //return 0.0;
    }

    /**
     * mean
     *
     * @param m  mean
     * @param shape  shape parameter
     * @return mean
     */
    public static double mean(double m, double shape) {
        return m;
    }

    /**
     * variance
     *
     * @param m  mean
     * @param shape  shape parameter
     * @return variance
     */
    public static double variance(double m, double shape) {
        double sd = calculateSD(m, shape);
        return sd * sd;
    }

    // Private

    protected double m, sd, shape;

}
