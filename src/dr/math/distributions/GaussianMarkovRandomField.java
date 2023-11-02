/*
 * GaussianMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

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
                                     Parameter precision,
                                     Parameter start) {
        this(dim, precision, start, null);
    }

    public GaussianMarkovRandomField(int dim,
                                     Parameter precision,
                                     Parameter start,
                                     RandomField.WeightProvider weightProvider) {

        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;
        this.meanParameter = start;
        this.precisionParameter = precision;
        this.weightProvider = weightProvider;

        this.mean = new double[dim];
        this.precision = new double[dim][dim];

        meanKnown = false;
        precisionKnown = false;
        determinantKnown = false; // TODO No need to be computed separately
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
            if (meanParameter == null) {
                Arrays.fill(mean, 0.0);
            } else if (meanParameter.getDimension() == 1) {
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

    @Override
    public double[][] getScaleMatrix() {
        return getPrecision();
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    @Override
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

    // TODO Below is the relevant code from GMRFMultilocusSkyrideLikelihood for building a `SymmTridiagMatrix`
    // TODO `getFieldScalar` rescaling should be handled by `WeightsProvider`

//    protected double getFieldScalar() {
//        final double rootHeight;
//        if (rescaleByRootHeight) {
//            rootHeight = tree.getNodeHeight(tree.getRoot());
//        } else {
//            rootHeight = 1.0;
//        }
//        return rootHeight;
//    }
//
//    protected void setupGMRFWeights() {
//
//        setupSufficientStatistics();
//
//        //Set up the weight Matrix
//        double[] offdiag = new double[fieldLength - 1];
//        double[] diag = new double[fieldLength];
//
//        //First set up the offdiagonal entries;
//
//        if (!timeAwareSmoothing) {
//            for (int i = 0; i < fieldLength - 1; i++) {
//                offdiag[i] = -1.0;
//            }
//        } else {
//            for (int i = 0; i < fieldLength - 1; i++) {
//                offdiag[i] = -2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1]) * getFieldScalar();
//            }
//        }
//
//        //Then set up the diagonal entries;
//        for (int i = 1; i < fieldLength - 1; i++)
//            diag[i] = -(offdiag[i] + offdiag[i - 1]);
//
//        //Take care of the endpoints
//        diag[0] = -offdiag[0];
//        diag[fieldLength - 1] = -offdiag[fieldLength - 2];
//
//        weightMatrix = new SymmTridiagMatrix(diag, offdiag);
//    }
//
//    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
//        SymmTridiagMatrix a = weightMatrix.copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, a.get(i, i) * precision);
//            a.set(i + 1, i, a.get(i + 1, i) * precision);
//        }
//        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
//        return a;
//    }
//
//    public SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda) {
//        if (lambda == 1)
//            return getScaledWeightMatrix(precision);
//
//        SymmTridiagMatrix a = weightMatrix.copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
//            a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
//        }
//
//        a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
//        return a;
//    }
//
//    private DenseVector getMeanAdjustedGamma() {
//        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
//        updateGammaWithCovariates(currentGamma);
//        return currentGamma;
//    }
//
//    double getLogFieldLikelihood() {
//
//        DenseVector diagonal1 = new DenseVector(fieldLength);
//        DenseVector currentGamma = getMeanAdjustedGamma();
//
//        double currentLike = handleMissingValues();
//
//        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
//        currentQ.mult(currentGamma, diagonal1);
//
//        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
//            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }
//
//        return currentLike;
//    }

    private static double logPdf(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q,
                                 double precision, double lambda) {

        final int dim = x.length;
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; ++i) {
            delta[i] = x[i] - mean[i];
        }

        double SSE = 0.0;
        for (int i = 0; i < dim - 1; i++) {
            SSE += Q.diagonal[i] * delta[i] * delta[i] + 2 * Q.offDiagonal[i] * delta[i] * delta[i + 1];
        }
        SSE += Q.diagonal[dim - 1] * delta[dim - 1] * delta[dim - 1];

        double logLikelihood = 0.5 * (dim - 1) * Math.log(precision) - 0.5 * SSE;
        if (lambda == 1.0) {
            logLikelihood -= (dim - 1) * logNormalize;
        } else {
            logLikelihood -= dim * logNormalize;
        }

        return logLikelihood;
    }


    class SymmetricTriDiagonalMatrix {

        double[] diagonal;
        double[] offDiagonal;

        SymmetricTriDiagonalMatrix(double[] diagonal, double[] offDiagonal) {
            this.diagonal = diagonal;
            this.offDiagonal = offDiagonal;
        }

        void copy(SymmetricTriDiagonalMatrix copy) {
            System.arraycopy(diagonal, 0, copy.diagonal, 0, diagonal.length);
            System.arraycopy(offDiagonal, 0, copy.offDiagonal, 0, offDiagonal.length);
        }

        void swap(SymmetricTriDiagonalMatrix swap) {
            double[] tmp1 = diagonal;
            diagonal = swap.diagonal;
            swap.diagonal = tmp1;

            double[] tmp2 = offDiagonal;
            offDiagonal = swap.offDiagonal;
            swap.offDiagonal = tmp2;
        }
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
    public int getDimension() { return dim; }

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
