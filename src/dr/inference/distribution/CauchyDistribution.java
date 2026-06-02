/*
 * CauchyDistribution.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.distribution;



import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;
import org.apache.commons.math.distribution.CauchyDistributionImpl;

public class CauchyDistribution extends AbstractContinuousDistribution implements Distribution {
    CauchyDistributionImpl distribution;
    double median;
    double scale;

    public CauchyDistribution(double median, double scale){
        distribution = new CauchyDistributionImpl(median, scale);
        this.median = median;
        this.scale = scale;
    }


    @Override
    public double pdf(double x) {
        return distribution.density(x);
    }

    @Override
    public double logPdf(double x) {
        return Math.log(distribution.density(x));
    }

    @Override
    public double cdf(double x) {
        return distribution.cumulativeProbability(x);
    }

    @Override
    public double quantile(double y) {
        return distribution.inverseCumulativeProbability(y);
    }

    @Override
    public double mean() {
        return Double.NaN;
    }

    @Override
    public double variance() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public final UnivariateFunction getProbabilityDensityFunction() {
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

    @Override
    protected double getInitialDomain(double v) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    protected double getDomainLowerBound(double v) {
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    protected double getDomainUpperBound(double v) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public double cumulativeProbability(double v) throws MathException {
        return 0;
    }
}
