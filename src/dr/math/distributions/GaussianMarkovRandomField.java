/*
 * GaussianMarkovRandomField.java
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
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.inference.distribution.RandomField;
import dr.inference.model.*;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import no.uib.cipr.matrix.*;

import java.util.Arrays;
import java.util.function.BinaryOperator;

/**
 * @author Marc Suchard
 * @author Pratyusa Data
 * @author Xiang Ji
 */
public class GaussianMarkovRandomField extends RandomFieldDistribution {

    public static final String TYPE = "GaussianMarkovRandomField";

    protected final int dim;
    private final Parameter meanParameter;
    protected final Parameter precisionParameter;
    protected final Parameter lambdaParameter;
    protected final RandomField.WeightProvider weightProvider;

    private final double[] mean;

    final SymmetricTriDiagonalMatrix Q;
    private final SymmetricTriDiagonalMatrix savedQ;

    private boolean meanKnown;
    protected boolean qKnown;
    private boolean savedQKnown;

    private final double logMatchTerm;

    private final int bandWidth;
    private int nonZeroEntryCount = -1;

    private boolean fieldDeterminantKnown;
    private boolean savedFieldDeterminantKnown;
    private double logFieldDeterminant;
    private double savedLogFieldDeterminant;

    public GaussianMarkovRandomField(String name,
                                     int dim,
                                     Parameter precision,
                                     Parameter mean,
                                     Parameter lambda,
                                     RandomField.WeightProvider weightProvider,
                                     boolean matchPseudoDeterminant) {
        super(name);

        this.dim = dim;
        this.bandWidth = 1;
        this.meanParameter = mean;
        this.precisionParameter = precision;
        this.lambdaParameter = lambda;
        this.weightProvider = weightProvider;

        if (meanParameter != null) {
            addVariable(meanParameter);
        }

        addVariable(precisionParameter);

        if (lambda != null) {
            addVariable(lambdaParameter);
        }

        if (weightProvider != null) {
            addModel(weightProvider);
        }

        this.mean = new double[dim];

        this.Q = new SymmetricTriDiagonalMatrix(dim);
        this.savedQ = new SymmetricTriDiagonalMatrix(dim);

        this.logMatchTerm = matchPseudoDeterminant ? matchPseudoDeterminantTerm(dim) : 0.0;

        meanKnown = false;
        qKnown = false;
        fieldDeterminantKnown = false;
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
            if (weightProvider == null) {
                Q.diagonal[0] = 1.0;
                for (int i = 1; i < dim - 1; ++i) {
                    Q.diagonal[i] = 2.0;
                }
                Q.diagonal[dim - 1] = 1.0;
                for (int i = 0; i < dim - 1; ++i) {
                    Q.offDiagonal[i] = -1.0;
                }
            } else {

                Q.diagonal[0] = weightProvider.weight(0, 1);
                for (int i = 1; i < dim - 1; ++i) {
                    Q.diagonal[i] = weightProvider.weight(i - 1, i) + weightProvider.weight(i, i + 1);
                }
                Q.diagonal[dim - 1] = weightProvider.weight(dim - 2, dim - 1);
                for (int i = 0; i < dim - 1; ++i) {
                    Q.offDiagonal[i] = -weightProvider.weight(i, i + 1);
                }
            }
            if (lambdaParameter != null) {
                double lambda = lambdaParameter.getParameterValue(0);
                for (int i = 0; i < dim - 1; ++i) {
                    Q.offDiagonal[i] = Q.offDiagonal[i] * lambda;
                    // TODO what about correction to Q.diagonal?
                }
            }

            qKnown = true;
        }
        return Q;
    }

    private static double[][] makePrecisionMatrix(SymmetricTriDiagonalMatrix Q, double scalar) {

        final int dim = Q.diagonal.length;
        double[][] precision = new double[dim][dim];

        for (int i = 0; i < dim; ++i) {
            precision[i][i] = Q.diagonal[i] * scalar;
        }

        for (int i = 0; i < dim - 1; ++i) {
            precision[i][i + 1] = Q.offDiagonal[i] * scalar;
            precision[i + 1][i] = Q.offDiagonal[i] * scalar;
        }

        return precision;
    }

    public double getFieldValue(int i, int j) {

        SymmetricTriDiagonalMatrix q = getQ();

        int whichDiagonal = Math.abs(i - j);
        if (whichDiagonal == 0) {
            return q.diagonal[i];
        } else if (whichDiagonal == 1) {
            return q.offDiagonal[Math.min(i, j)];
        } else {
            return 0.0;
        }
    }

    public int getNonZeroEntryCount() {
        if (nonZeroEntryCount == -1) {
            int[] count = new int[]{ 0 };
            internalMapAll((i, j, value) -> ++count[0]);
            nonZeroEntryCount = count[0];
        }
        return nonZeroEntryCount;
    }

    @FunctionalInterface
    public interface ReduceBlockFunction<R> {
        R apply(int i, int j, double fieldValue);
    }

    @FunctionalInterface
    public interface BlockFunction {
        void apply(int i, int j, double fieldValue);
    }

    @SuppressWarnings("unused")
    public void mapOverAllNonZeroEntries(BlockFunction map) {
        internalMapAll(map::apply);
    }

    @SuppressWarnings("unused")
    public void mapOverNonZeroEntriesInRow(BlockFunction map, int row) {
        internalMapRow(map::apply, row);
    }

    @SuppressWarnings({"unused", "unchecked"})
    public <R> R mapReduceOverAllNonZeroEntries(ReduceBlockFunction<R> map,
                                                BinaryOperator<R> reduce, R initial) {

        final R[] sum = (R[]) new Object[]{initial};
        internalMapAll((i, j, fieldValue) -> sum[0] = reduce.apply(sum[0], map.apply(i, j, fieldValue)));
        return sum[0];
    }

    @FunctionalInterface
    public interface InternalFunction {
        void apply(int i, int j, double value);
    }

    private void internalMapAll(InternalFunction function) {
        for (int i = 0; i < dim; ++i) {
            internalMapRow(function, i);
        }
    }

    private void internalMapRow(InternalFunction function, int row) {
        final int begin = Math.max(0, row - bandWidth);
        final int end = Math.min(dim, row + bandWidth + 1);
        for (int j = begin; j < end; ++j) {
            function.apply(row, j, getFieldValue(row, j));
        }
    }

    public boolean isImproper() {
        return lambdaParameter == null || lambdaParameter.getParameterValue(0) == 1.0;
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
                            precisionParameter.getParameterValue(0), isImproper());
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

                    double[] gradient = gradLogPdf((double[]) x, getMean(), precisionParameter.getParameterValue(0), getQ());

                    if (meanParameter.getDimension() == dim) {
                        for (int i = 0; i < dim; ++i) {
                            gradient[i] *= -1;
                        }
                        return gradient;
                    } else if (meanParameter.getDimension() == 1) {
                        double sum = 0.0;
                        for (int i = 0; i < dim; ++i) {
                            sum -= gradient[i];
                        }
                        return new double[]{sum};
                    }

                    throw new IllegalArgumentException("Unknown mean parameter structure");
                }
            };
        } else if (parameter == lambdaParameter) {
            throw new RuntimeException("Not yet implemented"); // TODO
        } else {
            throw new RuntimeException("Unknown parameter");
        }
    }

    public String getType() {
        return TYPE;
    }

    private double matchPseudoDeterminantTerm(int dim) {
        double term = 0.0;
        if (isImproper() && weightProvider == null) {
            for (int i = 1; i < dim; ++i) {
                double x = (2 - 2 * Math.cos(i * Math.PI / dim));
                term += Math.log(x);
            }
        }
        return term;
    }

    public double getLogDeterminant() {

        int effectiveDim = isImproper() ? dim - 1 : dim;
        double logDet = effectiveDim * Math.log(precisionParameter.getParameterValue(0)) + logMatchTerm;

        if (!fieldDeterminantKnown) {
            double logFieldDet = 0.0;
            if (!isImproper() || weightProvider != null) {
                double[][] precision = makePrecisionMatrix(Q, 1.0);
                RobustEigenDecomposition ed = new RobustEigenDecomposition(new DenseDoubleMatrix2D(precision));
                DoubleMatrix1D values = ed.getRealEigenvalues();
                for (int i = 0; i < values.size(); ++i) {
                    double v = values.get(i);
                    if (Math.abs(v) > 1E-6) {
                        logFieldDet += Math.log(v);
                    }
                }
            }

            logFieldDeterminant = logFieldDet;
            fieldDeterminantKnown = true;
        }

        logDet += logFieldDeterminant;

        if (CHECK_DETERMINANT) {

            double[][] precision = makePrecisionMatrix(Q, 1.0);
            RobustEigenDecomposition ed = new RobustEigenDecomposition(new DenseDoubleMatrix2D(precision));
            DoubleMatrix1D values = ed.getRealEigenvalues();
            double sum = 0.0;
            for (int i = 0; i < values.size(); ++i) {
                double v = values.get(i);
                if (Math.abs(v) > 1E-6) {
                    sum += Math.log(v);
                }
            }

            if (Math.abs(sum - logDet) > 1E-6) {
                throw new RuntimeException("Incorrect (pseudo-) determinant");
            }
        }

        return logDet;
    }

    public UpperTriangBandMatrix getCholeskyDecomposition() {

        SymmetricTriDiagonalMatrix Q = getQ();

        UpperSPDBandMatrix A = new UpperSPDBandMatrix(dim, 1);
        double precision = precisionParameter.getParameterValue(0);
        for (int i = 0; i < dim; ++i) {
            A.set(i, i, precision * Q.diagonal[i]);
        }

        for (int i = 0; i < dim - 1; ++i) {
            A.set(i, i + 1, precision * Q.offDiagonal[i]);
        }

        BandCholesky chol = BandCholesky.factorize(A);
        return chol.getU();
    }
    
    public static double[][] testCholeskyUpper(Matrix band, int dim) {
        double[][] result = new double[dim][dim];

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += band.get(k, i) * band.get(k, j); // U^t U
                }
                result[i][j] = sum;
            }
        }

        return result;
    }

    private static final boolean CHECK_DETERMINANT = false;

    @Override
    public double[][] getScaleMatrix() {
        return makePrecisionMatrix(getQ(), precisionParameter.getParameterValue(0));
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return meanParameter;
    }

    @Override
    public double logPdf(double[] x) {
        return logPdf(x, getMean(), precisionParameter.getParameterValue(0), getQ(), isImproper(), getLogDeterminant());
    }

    public static double gradLogPdfWrtPrecision(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q,
                                                double precision, boolean isImproper) {
        final int effectiveDim = isImproper ? x.length - 1 : x.length;
        return 0.5 * ((effectiveDim - getSSE(x, mean, precision, Q))/precision);
    }

    public static double[] gradLogPdf(double[] x, double[] mean, double precision, SymmetricTriDiagonalMatrix Q) {

        final int dim = x.length;

        final double[] gradient = new double[dim];
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; ++i) {
            delta[i] = mean[i] - x[i];
        }

        gradient[0] = precision * (Q.diagonal[0] * delta[0] + Q.offDiagonal[0] * delta[1]);
        for (int i = 1; i < dim - 1; ++i) {
            gradient[i] = precision * (Q.offDiagonal[i - 1] * delta[i - 1] + Q.diagonal[i] * delta[i] + Q.offDiagonal[i] * delta[i + 1]);
        }
        gradient[dim - 1] = precision * (Q.offDiagonal[dim - 2] * delta[dim - 2] + Q.diagonal[dim - 1] * delta[dim - 1]);

        return gradient;
    }

    public static double[][] hessianLogPdf(double[] x, double precision, SymmetricTriDiagonalMatrix Q) { // TODO test
        final int dim = x .length;
        final double[][] hessian = new double[dim][dim];

        hessian[0][0] = -precision * Q.diagonal[0];
        hessian[0][1] = -precision * Q.offDiagonal[0];

        for (int i = 1; i < dim - 1; ++i) {
            hessian[i][i - 1] = -precision * Q.offDiagonal[i - 1];
            hessian[i][i]     = -precision * Q.diagonal[i];
            hessian[i][i + 1] = -precision * Q.offDiagonal[i];
        }

        hessian[dim - 1][dim - 2] = -precision * Q.offDiagonal[dim - 2];
        hessian[dim - 1][dim - 1] = -precision * Q.diagonal[dim - 1];

        return hessian;
    }

    public static double[] diagonalHessianLogPdf(double[] x, double precision, SymmetricTriDiagonalMatrix Q) {
        final int dim = x.length;
        final double[] hessian = new double[dim];

        System.arraycopy(Q.diagonal, 0, hessian, 0, dim);
        // TODO do we not need to negate each element of hessian?

        for (int i = 0; i < dim; i++) {
            hessian[i] = hessian[i] * precision;
        }

        return hessian;
    }

    protected static double logPdf(double[] x, double[] mean, double precision, SymmetricTriDiagonalMatrix Q,
                                 boolean isImproper, double logDeterminant) {
        return getLogNormalization(x.length, isImproper, logDeterminant) - 0.5 * getSSE(x, mean, precision, Q);
    }

    private static double getSSE(double[] x, double[] mean, double precision, SymmetricTriDiagonalMatrix Q) {

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

        return SSE * precision;
    }

    protected static class SymmetricTriDiagonalMatrix {

        public double[] diagonal;
        public double[] offDiagonal;

        SymmetricTriDiagonalMatrix(int dim) {
            this(new double[dim], new double[dim - 1]);
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

    protected static double getLogNormalization(int dim, boolean isImproper, double logDeterminant) {
        final int effectiveDim = isImproper ? dim - 1 : dim;
        return -effectiveDim * HALF_LOG_TWO_PI + 0.5 * logDeterminant;
    }

    private static final double HALF_LOG_TWO_PI = Math.log(2.0 * Math.PI) / 2;

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x, getMean(), precisionParameter.getParameterValue(0), getQ());
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        return diagonalHessianLogPdf((double[]) x, precisionParameter.getParameterValue(0), getQ());
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        return hessianLogPdf((double[]) x, precisionParameter.getParameterValue(0), getQ());
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == weightProvider) {
            qKnown = false;
            fieldDeterminantKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == meanParameter) {
            meanKnown = false;
        } else if (variable == precisionParameter || variable == lambdaParameter) {
            qKnown = false;
            if (variable == lambdaParameter) {
                fieldDeterminantKnown = false;
            }
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
        if (qKnown) {
            Q.copyTo(savedQ);
        }
        savedQKnown = qKnown;

        savedLogFieldDeterminant = logFieldDeterminant;
        savedFieldDeterminantKnown = fieldDeterminantKnown;
    }

    @Override
    protected void restoreState() { // TODO cache mean
        meanKnown = false;

        qKnown = savedQKnown;
        if (qKnown) {
            savedQ.swap(Q);
        }

        logFieldDeterminant = savedLogFieldDeterminant;
        fieldDeterminantKnown = savedFieldDeterminantKnown;
    }

    @Override
    protected void acceptState() { }
}
