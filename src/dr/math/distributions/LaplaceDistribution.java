/*
 * LaplaceDistribution.java
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
 * @author Alexei Drummond
 */
public class LaplaceDistribution implements Distribution {

    // the mean parameter
    double mu;

    // the scale parameter
    double beta;

    // the maximum density
    double c;

    /**
     * This general constructor creates a new laplace distribution with a
     * specified mean and scale
     *
     * @param mu   mean
     * @param beta scale
     */
    public LaplaceDistribution(double mu, double beta) {
        setParameters(mu, beta);
    }

    /**
     * This default constructor creates a new laplace distribution with rate 1
     */
    public LaplaceDistribution() {
        this(0, 1);
    }

    /**
     * Set parameters of this distribution
     *
     * @param k mean
     * @param b scale
     */
    public void setParameters(double k, double b) {
        if (b <= 0) b = 1;
        mu = k;
        beta = b;

        //Normalizing constant
        c = 1 / (2 * beta);
    }

    /**
     * @return mean parameter
     */
    public double getMu() {
        return mu;
    }

    /**
     * @return scale parameter
     */
    public double getBeta() {
        return beta;
    }

    /**
     * @return maximum value of pdf function
     */
    public double getMaxDensity() {
        return c;
    }

    /**
     * Cumulative distribution function
     */
    public double cdf(double x) {
        if (x == mu) return 0.5;
        else return (0.5) * (1 + ((x - mu) / Math.abs(x - mu))
                * (1 - Math.exp(-Math.abs(x - mu) / beta)));
    }

    /**
     * Density function
     */
    public double pdf(double x) {
        return c * Math.exp(-Math.abs(x - mu) / beta);
    }

    public double logPdf(double x) {
        return Math.log(c) - (Math.abs(x - mu) / beta);
    }

    public double quantile(double y) {
        double sign = 1;
        if (y < 0.5) sign = -1;

        return mu - beta * sign * Math.log(1 - 2 * Math.abs(y - 0.5));
    }

    public double mean() {
        return mu;
    }

    public double variance() {
        return 2 * beta * beta;
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };
}
