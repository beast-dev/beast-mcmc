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

import dr.inference.model.GradientProvider;
import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import org.apache.commons.math.distribution.CauchyDistributionImpl;

public class CauchyDistribution implements Distribution, GradientProvider {
    private final CauchyDistributionImpl distribution;
    private final double median;
    private final double scale;

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
        public double evaluate(double x) {
            return pdf(x);
        }

        public double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object input) {
        double x = (double) input;

        double dev = x - median;
        double dx = 2 * dev / (dev * dev  + scale * scale);
        return new double[] { dx };
    }
}
