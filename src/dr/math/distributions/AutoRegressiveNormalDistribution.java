/*
 * AutoRegressiveNormalDistribution.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.model.HessianProvider;
import dr.inference.model.Likelihood;

/**
 * @author Marc Suchard
 */
public class AutoRegressiveNormalDistribution implements MultivariateDistribution, GaussianProcessRandomGenerator,
        GradientProvider, HessianProvider {

    public static final String TYPE = "AutoRegressiveNormal";

    private final int dim;
    private final double marginal;
    private final double decay;

    public AutoRegressiveNormalDistribution(int dim, double marginal, double decay) {
        this.dim = dim;
        this.marginal = marginal;
        this.decay = decay;

        if (marginal != 1.0) {
            throw new IllegalArgumentException("Not yet implemented");
        }
    }

    public String getType() {
        return TYPE;
    }

    private double getLogDet() {
        return  (1 - dim) * Math.log(1.0 - decay * decay);
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        return new double[dim];
    }

    @Override
    public double logPdf(double[] x) {

        double SSE = 0.0;
        for (int i = 1; i < dim; ++i) {
            final double delta =  x[i] - decay * x[i - 1];
            SSE += delta * delta;
        }

        return dim * logNormalize + 0.5 * (getLogDet() - SSE);
    }

    public double[] gradLogPdf(double[] x) {

        double[] gradient = new double[dim];

        gradient[0] = -(x[0] - decay * x[1]);

        for (int i = 1; i < dim - 1; ++i) {
            gradient[i] = -(-decay * x[i - 1] + x[i] - decay * x[i + 1]);
        }

        gradient[dim - 1] = -(-decay * x[dim - 2] + x[dim - 1]);

        return gradient;
    }

    private double[][] hessianLogPdf(double[] x) {
        throw new RuntimeException("Not yet implemented");
    }

    private double[] diagonalHessianLogPdf(double[] x) {
        throw new RuntimeException("Not yet implemented");
    }

    private static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

    @Override
    public Object nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    @Override
    public double[][] getPrecisionMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        return diagonalHessianLogPdf((double[]) x);
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        return hessianLogPdf((double[]) x);
    }
}
