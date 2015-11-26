/*
 * ChiSquareDistribution.java
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

/**
 * chi-square distribution
 * (distribution of sum of squares of n N(0,1) random variables)
 * <p/>
 * (Parameter: n; mean: n; variance: 2*n)
 * <p/>
 * The chi-square distribution is a special case of the Gamma distribution
 * (shape parameter = n/2.0, scale = 2.0).
 *
 * @author Korbinian Strimmer
 * @version $Id: ChiSquareDistribution.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public class ChiSquareDistribution extends GammaDistribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     */
    public ChiSquareDistribution(double n) {
        super(n / 2.0, 2.0);
        this.n = n;
    }

    public double pdf(double x) {
        return pdf(x, n);
    }

    public double cdf(double x) {
        return cdf(x, n);
    }

    public double quantile(double y) {
        return quantile(y, n);
    }

    public double mean() {
        return mean(n);
    }

    public double variance() {
        return variance(n);
    }

    /**
     * probability density function of the chi-square distribution
     *
     * @param x argument
     * @param n degrees of freedom
     * @return pdf value
     */
    public static double pdf(double x, double n) {
        return pdf(x, n / 2.0, 2.0);
    }

    /**
     * cumulative density function of the chi-square distribution
     *
     * @param x argument
     * @param n degrees of freedom
     * @return cdf value
     */
    public static double cdf(double x, double n) {
        return cdf(x, n / 2.0, 2.0);
    }

    /**
     * quantile (inverse cumulative density function) of the chi-square distribution
     *
     * @param y argument
     * @param n degrees of freedom
     * @return icdf value
     */
    public static double quantile(double y, double n) {
        return quantile(y, n / 2.0, 2.0);
    }

    /**
     * mean of the chi-square distribution
     *
     * @param n degrees of freedom
     * @return mean
     */
    public static double mean(double n) {
        return n;
    }

    /**
     * variance of the chi-square distribution
     *
     * @param n degrees of freedom
     * @return variance
     */
    public static double variance(double n) {
        return 2.0 * n;
    }

    // Private

    protected double n;
}
