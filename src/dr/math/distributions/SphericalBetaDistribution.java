/*
 * SphericalBetaDistribution.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.GradientProvider;
import dr.math.GammaFunction;
import dr.util.EuclideanToInfiniteNormUnitBallTransform;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class SphericalBetaDistribution implements MultivariateDistribution, GradientProvider {

    private double shape;
    private int dim;
    private double logNormalizationConstant;

    public static final String TYPE = "SphericalBetaDistribution";

    public SphericalBetaDistribution(int dim, double shape) {

        assert (shape > 0);

        this.shape = shape;
        this.dim = dim;
        this.logNormalizationConstant = computeLogNormalizationConstant();

    }

    SphericalBetaDistribution(int dim) { // returns a non-informative (uniform) density
        this(dim, 1.0);
    }

    private double computeLogNormalizationConstant() {
        // Lewandowski, Kurowicka, and Joe (2009) Lemma 7
        return GammaFunction.lnGamma(shape + 0.5 * dim)
                - 0.5 * dim * Math.log(Math.PI)
                - GammaFunction.lnGamma(shape);
    }

    @Override
    public double logPdf(double[] x) {
        assert (x.length == dim);
        return (shape == 1.0) ? logNormalizationConstant
                : logNormalizationConstant + (shape - 1) * Math.log(1 - squaredNorm(x));
    }

    private double squaredNorm(double[] x) {
        double norm = EuclideanToInfiniteNormUnitBallTransform.squaredNorm(x);
        assert (norm <= 1.0);
        return norm;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    private double[] gradLogPdf(double[] x) {
        assert (x.length == dim);

        if (shape == 1.0) { // Uniform
            return new double[x.length];
        } else {
            double factor = -2 * (shape - 1) / (1 - squaredNorm(x));
            double[] gradient = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                gradient[i] = factor * x[i];
            }
            return gradient;
        }
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }

}
