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
    private final double precisionScale; //scale term in the the formula to invert a AR1 matrix.

    public AutoRegressiveNormalDistribution(int dim, double marginal, double decay) {
        this.dim = dim;
        this.marginal = marginal;
        this.decay = decay;
        this.precisionScale = 1 / (1 - decay * decay);

        if (marginal != 1.0) {
            throw new IllegalArgumentException("Not yet implemented");
        }

        if (Math.abs(decay) >= 1.0) {
            throw new IllegalArgumentException("|Rho| must be < 1.0");
        }
    }

    public String getType() {
        return TYPE;
    }

    private double getLogDet() {
        return  (dim - 1) * Math.log(1.0 - decay * decay);
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        return new double[dim];
    }

    public double[] getPrecisionColumn(int index) {
        double[] column = new double[dim];

        if (index == 0) {
            column[0] = 1.0 * precisionScale;
            column[1] = -decay * precisionScale;
        } else if (index == dim - 1) {
            column[dim - 2] = -decay * precisionScale;
            column[dim - 1] = 1.0 * precisionScale;
        } else {
            column[index - 1] = -decay * precisionScale;
            column[index] = (1.0 + decay * decay) * precisionScale;
            column[index + 1] = -decay * precisionScale;
        }

        return column;
    }

    @Override
    public double logPdf(double[] x) {

        double SSE = x[0] * x[0] + x[dim - 1] * x[dim - 1];

        for (int i = 1; i < dim - 1; ++i) {
            SSE += (1.0 + decay * decay) * x[i] * x[i];
        }

        for (int i = 1; i < dim; ++i) {
            SSE -= 2 * decay * x[i - 1] * x[i];
        }

        SSE = SSE * precisionScale;

        final double logDet = getLogDet();

        return dim * logNormalize + 0.5 * (-logDet - SSE);
    }

    private double[] scaledPrecisionVectorProduct(double[] x, double scale) {

        assert (x.length == dim);

        double[] product = new double[dim];

        product[0] = scale * (x[0] - decay * x[1]) * precisionScale;

        for (int i = 1; i < dim - 1; ++i) {
            product[i] = scale * (-decay * x[i - 1] + (1 + decay * decay) * x[i] - decay * x[i + 1]) * precisionScale;
        }

        product[dim - 1] = scale * (-decay * x[dim - 2] + x[dim - 1]) * precisionScale;

        return product;
    }

    public double[] getDiagonal() {

        double[] diagonal = new double[dim];
        diagonal[0] = precisionScale;

        for (int i = 1; i < dim - 1; ++i) {
            diagonal[i] = (1 + decay * decay) * precisionScale;
        }

        diagonal[dim - 1] = precisionScale;

        return diagonal;
    }

    public double[] gradLogPdf(double[] x) {
        return scaledPrecisionVectorProduct(x, -1.0);
    }

    public double[] getPrecisionVectorProduct(double[] x) {
        return  scaledPrecisionVectorProduct(x, 1.0);
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
