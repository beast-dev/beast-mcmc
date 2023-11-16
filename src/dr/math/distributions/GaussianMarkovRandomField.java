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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.inference.distribution.RandomField;
import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.matrixAlgebra.RobustEigenDecomposition;

import java.util.Arrays;

/**
 * @author Marc Suchard
 * @author Pratyusa Data
 * @author Xiang Ji
 */
public class GaussianMarkovRandomField extends RandomFieldDistribution {

    public static final String TYPE = "GaussianMarkovRandomField";

    protected final int dim;
    private final Parameter meanParameter;
    private final Parameter precisionParameter;
    private final RandomField.WeightProvider weightProvider;

    private final double[] mean;
    private final double[][] precision; // TODO Use a sparse matrix, like in GmrfSkyrideLikelihood

    SymmetricTriDiagonalMatrix Q;
    private SymmetricTriDiagonalMatrix savedQ;

    private boolean meanKnown;
    private boolean precisionKnown;
    boolean qKnown;

    private final double logMatchTerm;

    public GaussianMarkovRandomField(String name,
                                     int dim,
                                     Parameter precision,
                                     Parameter start) {
        this(name, dim, precision, start, null, true);
    }

    public GaussianMarkovRandomField(String name,
                                     int dim,
                                     Parameter precision,
                                     Parameter start,
                                     RandomField.WeightProvider weightProvider,
                                     boolean matchPseudoDeterminant) {
        super(name);

        this.dim = dim;
        this.meanParameter = start;
        this.precisionParameter = precision;
        this.weightProvider = weightProvider;

        addVariable(meanParameter);
        addVariable(precisionParameter);

        if (weightProvider != null) {
            addModel(weightProvider);
        }

        this.mean = new double[dim];
        this.precision = new double[dim][dim];

        this.Q = new SymmetricTriDiagonalMatrix(dim);

        this.logMatchTerm = matchPseudoDeterminant ? matchPseudoDeterminantTerm(dim) : 0.0;

        meanKnown = false;
        precisionKnown = false;
        qKnown = false;
    }

    @Override
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

    protected SymmetricTriDiagonalMatrix getQ() {
        if (!qKnown) {
            double precision = precisionParameter.getParameterValue(0);
            Q.diagonal[0] = precision;
            for (int i = 1; i < dim - 1; ++i) {
                Q.diagonal[i] = 2 * precision;
            }
            Q.diagonal[dim - 1] = precision;

            for (int i = 0; i < dim - 1; ++i) {
                Q.offDiagonal[i] = -precision;
            }
            // TODO Update for lambda != 1 and for weights

            qKnown = true;
        }
        return Q;
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

    @Override
    public double getIncrement(int i, Parameter field) {

        double[] mean = getMean();
        return (field.getParameterValue(i) - mean[i]) - (field.getParameterValue(i + 1) - mean[i + 1]);
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {

        if (parameter == precisionParameter) {
            return new GradientProvider() {
                @Override
                public int getDimension() {
                    return 1;
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    double gradient =  gradLogPdfWrtPrecision((double[]) x, getMean(), getQ(),
                            precisionParameter.getParameterValue(0));
                    return new double[]{gradient};
                }
            };
        } else if (parameter == meanParameter) {
            return new GradientProvider() {
                @Override
                public int getDimension() {
                    return meanParameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {

                    double[] gradient = gradLogPdf((double[]) x, getMean(), getQ());

                    if (meanParameter.getDimension() == dim) {
                        for (int i = 0; i < dim; ++i) {
                            gradient[i] *= -1;
                        }
                        return gradient;
                    } else if (meanParameter.getDimension() == 1) {
                        double sum = 0.0;
                        for (int i = 0; i < dim; ++i) {
                            sum += gradient[i];
                        }
                        return new double[]{sum};
                    }

                    throw new IllegalArgumentException("Unknown mean parameter structure");
                }
            };
        } else {
            throw new RuntimeException("Unknown parameter");
        }
    }

    public String getType() {
        return TYPE;
    }

    private static double matchPseudoDeterminantTerm(int dim) {
        double term = 0.0;
        for (int i = 1; i < dim; ++i) {
            double x = (2 - 2 * Math.cos(i * Math.PI / dim));
            term += Math.log(x);
        }
        return term;
    }

    private double getLogDeterminant() {

        double logDet = (dim - 1) * Math.log(precisionParameter.getParameterValue(0)) + logMatchTerm;

        if (CHECK_DETERMINANT) {

            RobustEigenDecomposition ed = new RobustEigenDecomposition(new DenseDoubleMatrix2D(getPrecision()));
            DoubleMatrix1D values = ed.getRealEigenvalues();
            double sum = 0.0;
            for (int i = 0; i < values.size(); ++i) {
                double v = values.get(i);
                if (Math.abs(v) > 1E-6) {
                    sum += Math.log(v);
                }
            }

            if (Math.abs(sum - logDet) > 1E-6) {
                throw new RuntimeException("Incorrect pseudo-determinant");
            }
        }

        return logDet;
    }

    private static final boolean CHECK_DETERMINANT = false;

    @Override
    public double[][] getScaleMatrix() {
        return getPrecision();
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return meanParameter;
    }

    @Override
    public double logPdf(double[] x) {
//        double x1 = logPdf(x, getMean(), getPrecision(), getLogDet());
//        double x2 = logPdf(x, getMean(), getQ(), precisionParameter.getParameterValue(0),
//                1.0, getLogDeterminant());
//
//        System.err.println(x1 + " ?= " + x2);
//
        return logPdf(x, getMean(), getQ(), precisionParameter.getParameterValue(0),
                1.0, getLogDeterminant());
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

    public static double gradLogPdfWrtPrecision(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q,
                                                double precision) {
        final int dim = x.length;
        return 0.5 * (dim - 1 - getSSE(x, mean, Q)) / precision; // TODO Not correct with lambda != 1.0
    }

    public static double[] gradLogPdf(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q) {

        final int dim = x.length;

        final double[] gradient = new double[dim];
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; ++i) {
            delta[i] = mean[i] - x[i];
        }

        gradient[0] = Q.diagonal[0] * delta[0] + Q.offDiagonal[1] * delta[1];
        for (int i = 1; i < dim - 1; ++i) {
            gradient[i] = Q.offDiagonal[i - 1] * delta[i - 1] + Q.diagonal[i] * delta[i] + Q.offDiagonal[i] * delta[i + 1];
        }
        gradient[dim - 1] = Q.offDiagonal[dim - 2] * delta[dim - 2] + Q.diagonal[dim - 1] * delta[dim - 1];

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

    public static double[][] hessianLogPdf(double[] x, SymmetricTriDiagonalMatrix Q) { // TODO test
        final int dim = x .length;
        final double[][] hessian = new double[dim][dim];

        hessian[0][0] = -Q.diagonal[0];
        hessian[0][1] = -Q.offDiagonal[0];

        for (int i = 1; i < dim - 1; ++i) {
            hessian[i][i - 1] = -Q.offDiagonal[i - 1];
            hessian[i][i]     = -Q.diagonal[i];
            hessian[i][i + 1] = -Q.offDiagonal[i];
        }

        hessian[dim - 1][dim - 2] = -Q.offDiagonal[dim - 2];
        hessian[dim - 1][dim - 1] = -Q.diagonal[dim - 1];

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

    public static double[] diagonalHessianLogPdf(double[] x, SymmetricTriDiagonalMatrix Q) {
        final int dim = x.length;
        final double[] hessian = new double[dim];

        System.arraycopy(Q.diagonal, 0, hessian, 0, dim);

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
                                 double precision, double lambda, double logDeterminant) {
        return getLogNormalization(x.length, precision, lambda, logDeterminant) - 0.5 * getSSE(x, mean, Q);
    }

    private static double getSSE(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q) {

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

        return SSE;
    }

    private static double getLogNormalization(int dim, double precision, double lambda, double logDeterminant) {

        double logNorm = 0.5 * logDeterminant;

        if (lambda == 1.0) {
            logNorm -= (dim - 1) * HALF_LOG_TWO_PI;
        } else {
            logNorm -= dim * HALF_LOG_TWO_PI;
        }

        return logNorm;
    }

    static class SymmetricTriDiagonalMatrix {

        double[] diagonal;
        double[] offDiagonal;

        SymmetricTriDiagonalMatrix(int dim) {
            this.diagonal = new double[dim];
            this.offDiagonal = new double[dim - 1];
        }

        SymmetricTriDiagonalMatrix(double[] diagonal, double[] offDiagonal) {
            this.diagonal = diagonal;
            this.offDiagonal = offDiagonal;
        }

        void copyTo(SymmetricTriDiagonalMatrix copy) {
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

    public static double logPdf(double[] x, double[] mean, double[][] precision,
                                double logDet) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        return getLogNormalization(x.length, logDet) - 0.5 * getSSE(x, mean, precision);
    }

    public static double getLogNormalization(int dim, double logDet) {
        return -(dim - 1) * HALF_LOG_TWO_PI + 0.5 * logDet;
    }

    public static double getSSE(double[] x, double[] mean, double[][] precision) {

        final int dim = x.length;

        final double[] delta = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        double SSE = precision[dim-1][dim-1] * delta[dim-1] * delta[dim-1];

        for (int i = 0; i < dim-1; i++) {
            SSE += precision[i][i] * delta[i] * delta[i] + 2 * precision[i][i + 1] * delta[i] * delta[i + 1];
        }

        return SSE;
    }

    private static final double HALF_LOG_TWO_PI = Math.log(2.0 * Math.PI) / 2;

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
//        double[] g1 = gradLogPdf((double[]) x, getMean(), getPrecision());
//        double[] g2 = gradLogPdf((double[]) x, getMean(), getQ());
//
//        System.err.println(new WrappedVector.Raw(g1));
//        System.err.println(new WrappedVector.Raw(g2));

        return gradLogPdf((double[]) x, getMean(), getQ());
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        // TODO update for Q
        return diagonalHessianLogPdf((double[]) x, getPrecision());
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        // TODO update for Q
        return hessianLogPdf((double[]) x, getPrecision());
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
            qKnown = false;
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
        qKnown = false;
    }

    @Override
    protected void acceptState() { }
}
