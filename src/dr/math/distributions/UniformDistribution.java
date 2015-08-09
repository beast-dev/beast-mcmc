/*
 * UniformDistribution.java
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
import dr.util.DataTable;

/**
 * uniform distribution.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: UniformDistribution.java,v 1.3 2005/07/11 14:06:25 rambaut Exp $
 */
public class UniformDistribution implements Distribution {
    //
    // Public stuff
    //

    /*
     * Constructor
     */

    public UniformDistribution(double lower, double upper) {
        this.lower = lower;
        this.upper = upper;
    }

    public double pdf(double x) {
        return pdf(x, lower, upper);
    }

    public double logPdf(double x) {
        return logPdf(x, lower, upper);
    }

    public double cdf(double x) {
        return cdf(x, lower, upper);
    }

    public double quantile(double y) {
        return quantile(y, lower, upper);
    }

    public double mean() {
        return mean(lower, upper);
    }

    public double variance() {
        return variance(lower, upper);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return lower;
        }

        public final double getUpperBound() {
            return upper;
        }
    };


    /**
     * probability density function of the uniform distribution
     *
     * @param x     argument
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return pdf value
     */
    public static double pdf(double x, double lower, double upper) {
        return (x >= lower && x <= upper ? 1.0 / (upper - lower) : 0.0);
    }

    /**
     * the natural log of the probability density function of the uniform distribution
     *
     * @param x     argument
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return log pdf value
     */
    public static double logPdf(double x, double lower, double upper) {
        if (x < lower || x > upper) return Double.NEGATIVE_INFINITY;

        // improve numerical stability:
        return - Math.log(upper - lower);
//        return Math.log(pdf(x, lower, upper));
    }

    /**
     * cumulative density function of the uniform distribution
     *
     * @param x     argument
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return cdf value
     */
    public static double cdf(double x, double lower, double upper) {
        if (x < lower) return 0.0;
        if (x > upper) return 1.0;
        return (x - lower) / (upper - lower);
    }


    /**
     * quantile (inverse cumulative density function) of the uniform distribution
     *
     * @param y     argument
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return icdf value
     */
    public static double quantile(double y, double lower, double upper) {
        if (!(y >= 0.0 && y <= 1.0)) throw new IllegalArgumentException("y must in range [0,1]");
        return (y * (upper - lower)) + lower;
    }

    /**
     * mean of the uniform distribution
     *
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return mean
     */
    public static double mean(double lower, double upper) {
        return (upper + lower) / 2;
    }

    /**
     * variance of the uniform distribution
     *
     * @param lower the lower bound of the uniform distribution
     * @param upper the upper bound of the uniform distribution
     * @return variance
     */
    public static double variance(double lower, double upper) {
        return (upper - lower) * (upper - lower) / 12;
    }

    // Private

    private final double upper;
    private final double lower;
}
