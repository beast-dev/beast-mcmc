/*
 * MultivariateNormalDistribution.java
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

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.matrixAlgebra.*;

/**
 * @author Marc Suchard
 */
public class GaussianMarkovRandomField extends RandomFieldDistribution {

    public static final String TYPE = "GaussianMarkovRandomField";

    private final double[] mean;
    private final double[][] precision;
    private final int dim;
    private final Parameter start;
    private final Parameter incrementPrecision;
    private double[][] variance = null;
    private double[][] cholesky = null;
    private Double logDet = null;

    public GaussianMarkovRandomField(int dim, Parameter incrementPrecision, Parameter start) {

        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;
        this.start = start;
        this.mean = new double[dim];
        final double x0 = start.getParameterValue(0);
        for(int i=0; i<dim; ++i) {
            this.mean[i] = x0;
        }
        this.incrementPrecision = incrementPrecision;
        final double k = incrementPrecision.getParameterValue(0);
        this.precision = new double[dim][dim];
        this.precision[0][0] = k;
        this.precision[0][1] = -1*k;
        this.precision[dim-1][dim-1] = k;
        this.precision[dim-1][dim-2] = -1*k;
        for (int i = 1; i < dim-1; ++i) {
            this.precision[i][i] = 2*k;
            this.precision[i][i-1] = -1*k;
            this.precision[i][i+1] = -1*k;
        }
    }

    public GradientProvider getGradientWrt(Parameter parameter) {
        // TODO Should return a GradientProvider for the specified the hyper-parameter
        throw new RuntimeException("Not yet implemented");
    }

    public String getType() {
        return TYPE;
    }

    public double[][] getVariance() {
        final double k = incrementPrecision.getParameterValue(0);
        if (variance == null) {

            for (int i=0; i<dim; ++i) {
               for (int j=0; j<dim; ++j) {
                   if(j == i) variance[j][j] = j / k;
                   else {
                       variance[i][j] = Math.abs(j-i)/k;
                   }
               }
            }
        }
        return variance;
    }

    public double[][] getCholeskyDecomposition() {
        if (cholesky == null) {
            cholesky = getCholeskyDecomposition(getVariance());
        }
        return cholesky;
    }

    public double getLogDet() {
        final double k = incrementPrecision.getParameterValue(0);
        if (logDet == null) {
            double det = Math.pow(k, dim);
            for(int i=2; i<=dim; ++i) {
                det = det * (2 - 2 * Math.cos((i-1)*(Math.PI/dim)));
            }
            logDet = Math.log(det);
        }

        return logDet;
    }

    private boolean isDiagonal(double x[][]) {
        for (int i = 0; i < x.length; ++i) {
            for (int j = i + 1; j < x.length; ++j) {
                if (x[i][j] != 0.0) {
                    return false;
                }
            }
        }
        return true;
    }

    private double logDetForDiagonal(double x[][]) {
        double logDet = 0;
        for (int i = 0; i < x.length; ++i) {
            logDet += Math.log(x[i][i]);
        }
        return logDet;
    }


    public double[][] getScaleMatrix() {
        return precision;
    }

    public double[] getMean() {
        return mean;
    }


    public static double calculatePrecisionMatrixDeterminate(double[][] precision) {
        try {
            return new Matrix(precision).determinant();
        } catch (IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }


    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    public double logPdf(double[] x) {

            return logPdf(x, mean, precision, getLogDet());

    }

    public double[] gradLogPdf(double[] x) {

            return gradLogPdf(x, mean, precision);

    }



    public static double[] gradLogPdf(double[] x, double[] mean, double[][] precision) {
        final int dim = x.length;
        final double[] gradient = new double[dim];
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; ++i) {
            delta[i] = mean[i] - x[i];
        }
        gradient[0] = precision[0][0] * delta[0] + precision[0][1] * delta[1];
        gradient[dim-1] = precision[dim-1][dim-2] * delta[dim-2] + precision[dim-1][dim-1] * delta[dim-1];
        for (int i = 1; i < dim-1; ++i) {
            gradient[i] = precision[i][i-1] * delta[i-1] + precision[i][i] * delta[i] + precision[i][i+1] * delta[i+1];
        }
        return gradient;
    }

    public double[][] hessianLogPdf(double[] x) {

            return hessianLogPdf(x, mean, precision);
    }



    public static double[][] hessianLogPdf(double[] x, double[] mean, double[][] precision) {

        final int dim = x .length;
        final double[][] hessian = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                hessian[i][j] = -precision[i][j];
            }
        }
        return hessian;
    }

    public double[] diagonalHessianLogPdf(double[] x) {

            return diagonalHessianLogPdf(x, mean, precision);

    }



    public static double[] diagonalHessianLogPdf(double[] x, double[] mean, double[][] precision) {
        final int dim = x.length;
        final double[] hessian = new double[dim];

        for (int i = 0; i < dim; ++i) {
            hessian[i] = -precision[i][i];
        }

        return hessian;
    }

    // scale only modifies precision
    // in one dimension, this is equivalent to:
    // PDF[NormalDistribution[mean, Sqrt[scale]*Sqrt[1/precison]], x]
    public double logPdf(double[] x, double[] mean, double[][] precision,
                                double logDet) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        final int dim = x.length;
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        double SSE = precision[dim-1][dim-1] * delta[dim-1] * delta[dim-1];

        for (int i = 0; i < dim-1; i++) {
            SSE += precision[i][i] * delta[i] * delta[i] + 2 * precision[i][i + 1] * delta[i] * delta[i + 1];
        }
        return (dim-1) * logNormalize + 0.5 * (logDet - SSE);   // There was an error here.
        // Variance = (scale * Precision^{-1})
    }



    private static double[][] getInverse(double[][] x) {
        return new SymmetricMatrix(x).inverse().toComponents();
    }

    private static double[][] getCholeskyDecomposition(double[][] variance) {
        double[][] cholesky;
        try {
            cholesky = (new CholeskyDecomposition(variance)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
        }
        return cholesky;
    }



    private static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);


    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }



    @Override
    public int getDimension() { return mean.length; }

    public Parameter getincrementPrecision() { return incrementPrecision; }

    public Parameter getstart() { return start; }


    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }



    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        return diagonalHessianLogPdf((double[]) x);
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        return hessianLogPdf((double[]) x);
    }


    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Should be no sub-models");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // TODO
    }

    @Override
    protected void storeState() {
        // TODO
    }

    @Override
    protected void restoreState() {
        // TODO
    }

    @Override
    protected void acceptState() { }
}
