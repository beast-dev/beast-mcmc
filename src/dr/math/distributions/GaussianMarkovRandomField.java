/*
 * GaussianMarkovRandomField.java
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

import dr.inference.distribution.RandomField;
import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;

import java.util.Arrays;

/**
 * @author Marc Suchard
 */
public class GaussianMarkovRandomField extends RandomFieldDistribution {

    public static final String TYPE = "GaussianMarkovRandomField";

    private final int dim;
    private final Parameter meanParameter;
    private final Parameter precisionParameter;
    private final RandomField.WeightProvider weightProvider;

    private final double[] mean;
    private final double[][] precision; // TODO Use a sparse matrix, like in GmrfSkyrideLikelihood
    private double logDet;

//    private double[][] variance = null;
//    private double[][] cholesky = null;

    private boolean meanKnown;
    private boolean precisionKnown;
    private boolean determinantKnown;

    public GaussianMarkovRandomField(int dim,
                                     Parameter incrementPrecision,
                                     Parameter start) {
        this(dim, incrementPrecision, start, null);
    }

    public GaussianMarkovRandomField(int dim,
                                     Parameter incrementPrecision,
                                     Parameter start,
                                     RandomField.WeightProvider weightProvider) {

        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;
        this.meanParameter = start;
        this.precisionParameter = incrementPrecision;
        this.weightProvider = weightProvider;

        this.mean = new double[dim];
        this.precision = new double[dim][dim];

//        populateMean(this.mean);
//        populatePrecision(this.precision);

        meanKnown = false;
        precisionKnown = false;
        determinantKnown = false;
    }

//    private void check() {
//        if (!meanKnown) {
//            populateMean(mean);
//            meanKnown = true;
//        }
//        if (!precisionKnown) {
//            populatePrecision(precision);
//            precisionKnown = true;
//        }
//    }

    public double[] getMean() {
        if (!meanKnown) {
            if (meanParameter.getDimension() == 1) {
                Arrays.fill(mean, meanParameter.getParameterValue(0));
            } else {
                for (int i = 0; i < mean.length; ++i) {
                    mean[i] = meanParameter.getParameterValue(i);
                }
            }

            meanKnown = true;
        }
        return mean;
    }

    private double[][] getPrecision() {

        if (!precisionKnown) {
            final double k = precisionParameter.getParameterValue(0);

            precision[0][0] = k;
            precision[0][1] = -1 * k;
            precision[dim - 1][dim - 1] = k;
            precision[dim - 1][dim - 2] = -1 * k;
            for (int i = 1; i < dim - 1; ++i) {
                precision[i][i] = 2 * k;
                precision[i][i - 1] = -1 * k;
                precision[i][i + 1] = -1 * k;
            }

            precisionKnown = true;
        }
        return precision;
    }

    public GradientProvider getGradientWrt(Parameter parameter) {
        // TODO Should return a GradientProvider for the specified the hyper-parameter
        throw new RuntimeException("Not yet implemented");
    }

    public String getType() {
        return TYPE;
    }

//    public double[][] getVariance() {
//        final double k = precisionParameter.getParameterValue(0);
//        if (variance == null) {
//
//            for (int i=0; i<dim; ++i) {
//               for (int j=0; j<dim; ++j) {
//                   if(j == i) variance[j][j] = j / k;
//                   else {
//                       variance[i][j] = Math.abs(j-i)/k;
//                   }
//               }
//            }
//        }
//        return variance;
//    }
//
//    public double[][] getCholeskyDecomposition() {
//        if (cholesky == null) {
//            cholesky = getCholeskyDecomposition(getVariance());
//        }
//        return cholesky;
//    }

    public double getLogDet() {

        if (!determinantKnown) {
            final double k = precisionParameter.getParameterValue(0);
            double det = Math.pow(k, dim);
            for(int i=2; i<=dim; ++i) {
                det = det * (2 - 2 * Math.cos((i-1)*(Math.PI/dim)));
            }
            logDet = Math.log(det);

            determinantKnown = true;
        }

        return logDet;
    }

//    private boolean isDiagonal(double x[][]) {
//        for (int i = 0; i < x.length; ++i) {
//            for (int j = i + 1; j < x.length; ++j) {
//                if (x[i][j] != 0.0) {
//                    return false;
//                }
//            }
//        }
//        return true;
//    }
//
//    private double logDetForDiagonal(double x[][]) {
//        double logDet = 0;
//        for (int i = 0; i < x.length; ++i) {
//            logDet += Math.log(x[i][i]);
//        }
//        return logDet;
//    }


    public double[][] getScaleMatrix() {
        return getPrecision();
    }

//    public static double calculatePrecisionMatrixDeterminate(double[][] precision) {
//        try {
//            return new Matrix(precision).determinant();
//        } catch (IllegalDimension e) {
//            throw new RuntimeException(e.getMessage());
//        }
//    }


    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    public double logPdf(double[] x) {
        return logPdf(x, getMean(), getPrecision(), getLogDet());
    }

    public double[] gradLogPdf(double[] x) {
        return gradLogPdf(x, getMean(), getPrecision());
    }

    public double[][] hessianLogPdf(double[] x) {
        return hessianLogPdf(x, getPrecision());
    }

    public double[] diagonalHessianLogPdf(double[] x) {
        return diagonalHessianLogPdf(x, getPrecision());
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

    public static double[][] hessianLogPdf(double[] x, double[][] precision) {
        final int dim = x .length;
        final double[][] hessian = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                hessian[i][j] = -precision[i][j];
            }
        }
        return hessian;
    }

    public static double[] diagonalHessianLogPdf(double[] x, double[][] precision) {
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
    public static double logPdf(double[] x, double[] mean, double[][] precision,
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

//    private static double[][] getInverse(double[][] x) {
//        return new SymmetricMatrix(x).inverse().toComponents();
//    }

//    private static double[][] getCholeskyDecomposition(double[][] variance) {
//        double[][] cholesky;
//        try {
//            cholesky = (new CholeskyDecomposition(variance)).getL();
//        } catch (IllegalDimension illegalDimension) {
//            throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
//        }
//        return cholesky;
//    }

    private static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

//    public double logPdf(Object x) {
//        double[] v = (double[]) x;
//        return logPdf(v);
//    }


    @Override
    public int getDimension() { return mean.length; }

//    public Parameter getincrementPrecision() { return precisionParameter; }

//    public Parameter getstart() { return meanParameter; }


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
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == meanParameter) {
            meanKnown = false;
        } else if (variable == precisionParameter) {
            precisionKnown = false;
            determinantKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
        // TODO
    }

    @Override
    protected void restoreState() { // TODO with caching
        meanKnown = false;
        precisionKnown = false;
        determinantKnown = false;
    }

    @Override
    protected void acceptState() { }
}
