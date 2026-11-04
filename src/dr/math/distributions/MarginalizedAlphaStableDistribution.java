/*
 * MarginalizedAlphaStableDistribution.java
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
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class MarginalizedAlphaStableDistribution implements Distribution {

    private final double scale;
    private final double alpha;

    public MarginalizedAlphaStableDistribution(double scale, double alpha) {
        this.scale = scale;
        this.alpha = alpha;
    }

    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    @Override
    public double logPdf(double x) {
        return logPdf(x, scale, alpha);
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    public static double logPdf(double x, double scale, double alpha) {
        return -Math.log(scale) - Math.pow(Math.abs(x) / scale, alpha);
    }

    public static double gradLogPdf(double x, double scale, double alpha) {
        return -alpha * Math.pow(Math.abs(x) / scale, alpha - 1.0) * gradAbsX(x) / scale;
    }

    private static double gradAbsX(double x) {
        if (x < 0) {
            return -1;
        } else if (x > 0) {
            return 1;
        } else {
            return Double.NaN;
        }
    }
}
