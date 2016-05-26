/*
 * OffsetPositiveDistribution.java
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
 * A distribution that is offset, so that the origin is greater than 0
 *
 * @author Alexei Drummond
 */
public class OffsetPositiveDistribution implements Distribution {

    /**
     * Constructor
     *
     * @param distribution distribution to offset
     * @param offset       a (positive) location parameter that allows this distribution to start
     *                     at a non-zero location
     */
    public OffsetPositiveDistribution(Distribution distribution, double offset) {

        if (offset < 0.0) throw new IllegalArgumentException();
        this.offset = offset;
        this.distribution = distribution;
    }

    /**
     * probability density function of the offset distribution
     *
     * @param x argument
     * @return pdf value
     */
    public final double pdf(double x) {
        if (offset < 0) return 0.0;
        return distribution.pdf(x - offset);
    }

    /**
     * log probability density function of the offset distribution
     *
     * @param x argument
     * @return pdf value
     */
    public final double logPdf(double x) {
        if (offset < 0) return Math.log(0.0);
        return distribution.logPdf(x - offset);
    }

    /**
     * cumulative density function of the offset distribution
     *
     * @param x argument
     * @return cdf value
     */
    public final double cdf(double x) {
        if (offset < 0) return 0.0;
        return distribution.cdf(x - offset);
    }

    /**
     * quantile (inverse cumulative density function) of the (offset) distribution
     *
     * @param y the p-value
     * @return icdf value
     */
    public final double quantile(double y) {
        return distribution.quantile(y) + offset;
    }

    /**
     * mean of the offset distribution
     *
     * @return mean
     */
    public final double mean() {
        return distribution.mean() + offset;
    }

    /**
     * variance of the offset distribution
     *
     * @return variance
     */
    public final double variance() {
        throw new UnsupportedOperationException();
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return offset + distribution.getProbabilityDensityFunction().getLowerBound();
        }

        public final double getUpperBound() {
            return offset + distribution.getProbabilityDensityFunction().getUpperBound();
        }
    };

    // the location parameter of the start of the positive distribution
    private double offset = 0.0;

    // the distribution to offset
    private Distribution distribution;
}
