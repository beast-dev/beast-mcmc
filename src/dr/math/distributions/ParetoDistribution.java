/*
 * ParetoDistribution.java
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

public class ParetoDistribution implements Distribution {

    private double scale;
    private double shape;

    public ParetoDistribution(double scale, double shape) {
        if (scale <= 0 || shape <= 0) {
            throw new RuntimeException("Shape and scale must be positive.");
        }
        this.scale = scale;
        this.shape = shape;
    }

    @Override
    public double pdf(double x) {
        return x > scale ? shape * Math.pow(scale, shape) / Math.pow(x, shape + 1) : 0.0;
    }

    @Override
    public double logPdf(double x) {
        if (x < scale) return Double.NEGATIVE_INFINITY;
        return Math.log(shape) + shape * Math.log(scale) - (shape + 1) * Math.log(x);
    }

    @Override
    public double cdf(double x) {
        if (x < scale) return 0.0;
        return 1 - Math.pow(scale / x, shape);
    }

    @Override
    public double quantile(double p) {
        return 1.0 / Math.pow(1.0 - p, 1 / shape);
    }//todo

    @Override
    public double mean() {
        if (shape <= 2) return Double.POSITIVE_INFINITY;
        else return scale * shape / (shape - 1);
    }

    @Override
    public double variance() {
        if (shape <= 2) return Double.POSITIVE_INFINITY;
        else return scale * scale * shape / (shape - 1) * (shape - 1) * (shape - 2);
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        return null;
    }
}
