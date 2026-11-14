/*
 * NegativeBinomialDistribution.java
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

package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * @author Marc Suchard
 */
public class GeometricDistribution implements Distribution {

    private final double p;

    public GeometricDistribution(double p) {
        this.p = p;
    }

    public double pdf(double x) {
        return pdf((int) x, p);
    }

    public double logPdf(double x) {
        return logPdf((int) x, p);
    }

    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    public double quantile(double y) { throw new RuntimeException("Not yet implemented"); }

    public double mean() {
        return 1.0 / p;
    }

    public double variance() { return (1.0 - p) / (p * p); }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    public static double pdf(int k, double p) {
        if (k < 1)  return 0;
        return Math.pow(1.0 - p, k - 1) * p;
    }

    public static double logPdf(int k, double p) {
        if (k < 1)  return Double.NEGATIVE_INFINITY;
        return (k - 1) * Math.log(1.0 - p) + Math.log(p);
    }
}
