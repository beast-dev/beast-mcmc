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

import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.*;

/**
 * @author Marc Suchard
 */
public class MultivariateNormalDistribution implements MultivariateDistribution, GaussianProcessRandomGenerator,
        GradientProvider {

    public static final String TYPE = "MultivariateNormal";

    private final double[] mean;
    private final double[][] precision;
    private double[][] variance = null;
    private double[][] cholesky = null;
    private Double logDet = null;

    private final boolean hasSinglePrecision;
    private final double singlePrecision;

    public MultivariateNormalDistribution(double[] mean, double[][] precision) {
        this.mean = mean;
        this.precision = precision;
        this.hasSinglePrecision = false;
        this.singlePrecision = 1.0;
    }

    public MultivariateNormalDistribution(double[] mean, double singlePrecision) {
        this.mean = mean;
        this.hasSinglePrecision = true;
        this.singlePrecision = singlePrecision;

        final int dim = mean.length;
        this.precision = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            this.precision[i][i] = singlePrecision;
        }
    }

    public String getType() {
        return TYPE;
    }

    public double[][] getVariance() {
        if (variance == null) {
            variance = new SymmetricMatrix(precision).inverse().toComponents();
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
        if (logDet == null) {
            logDet = Math.log(calculatePrecisionMatrixDeterminate(precision));
        }
        if (Double.isInfinite(logDet)) {
            if (isDiagonal(precision)) {
                logDet = logDetForDiagonal(precision);
            }
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

    public double[] nextMultivariateNormal() {
        return nextMultivariateNormalCholesky(mean, getCholeskyDecomposition(), 1.0);
    }

    public double[] nextMultivariateNormal(double[] x) {
        return nextMultivariateNormalCholesky(x, getCholeskyDecomposition(), 1.0);
    }

    // Scale lives in variance-space
    public double[] nextScaledMultivariateNormal(double[] mean, double scale) {
        return nextMultivariateNormalCholesky(mean, getCholeskyDecomposition(), Math.sqrt(scale));
    }

    // Scale lives in variance-space
    public void nextScaledMultivariateNormal(double[] mean, double scale, double[] result) {
        nextMultivariateNormalCholesky(mean, getCholeskyDecomposition(), Math.sqrt(scale), result);
    }


    public static double calculatePrecisionMatrixDeterminate(double[][] precision) {
        try {
            return new Matrix(precision).determinant();
        } catch (IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public double logPdf(double[] x) {
        if (hasSinglePrecision) {
            return logPdf(x, mean, singlePrecision, 1.0);
        } else {
            return logPdf(x, mean, precision, getLogDet(), 1.0);
        }
    }

    public double[] gradLogPdf(double[] x) {
        if (hasSinglePrecision) {
            return gradLogPdf(x, mean, singlePrecision);
        } else {
            return gradLogPdf(x, mean, precision);
        }
    }

    public static double[] gradLogPdf(double[] x, double[] mean, double singlePrecision) {

        final int dim = x .length;
        final double[] gradient = new double[dim];

        for (int i = 0; i < dim; ++i) {
            gradient[i] = singlePrecision * (mean[i] - x[i]);
        }

        return gradient;
    }

    public static double[] gradLogPdf(double[] x, double[] mean, double[][] precision) {
        final int dim = x.length;
        final double[] gradient = new double[dim];
        final double[] delta = new double[dim];

        for (int i = 0; i < dim; ++i) {
            delta[i] = mean[i] - x[i];
        }

        for (int i = 0; i < dim; ++i) {
            double sum = 0;
            for (int j = 0; j <dim; ++j) {
                sum += precision[i][j] * delta[j];
            }
            gradient[i] = sum;
        }

        return gradient;
    }

    // scale only modifies precision
    // in one dimension, this is equivalent to:
    // PDF[NormalDistribution[mean, Sqrt[scale]*Sqrt[1/precison]], x]
    public static double logPdf(double[] x, double[] mean, double[][] precision,
                                double logDet, double scale) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        final int dim = x.length;
        final double[] delta = new double[dim];
        final double[] tmp = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                tmp[i] += delta[j] * precision[j][i];
            }
        }

        double SSE = 0;

        for (int i = 0; i < dim; i++)
            SSE += tmp[i] * delta[i];

        return dim * logNormalize + 0.5 * (logDet - dim * Math.log(scale) - SSE / scale);   // There was an error here.
        // Variance = (scale * Precision^{-1})
    }

    /* Equal precision, independent dimensions */
    public static double logPdf(double[] x, double[] mean, double precision, double scale) {

        final int dim = x.length;

        double SSE = 0;
        for (int i = 0; i < dim; i++) {
            double delta = x[i] - mean[i];
            SSE += delta * delta;
        }

        return dim * logNormalize + 0.5 * (dim * (Math.log(precision) - Math.log(scale)) - SSE * precision / scale);
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

    public static double[] nextMultivariateNormalPrecision(double[] mean, double[][] precision) {
        return nextMultivariateNormalVariance(mean, getInverse(precision));
    }

    public static double[] nextMultivariateNormalVariance(double[] mean, double[][] variance) {
        return nextMultivariateNormalVariance(mean, variance, 1.0);
    }

    public static double[] nextMultivariateNormalVariance(double[] mean, double[][] variance, double scale) {
        return nextMultivariateNormalCholesky(mean, getCholeskyDecomposition(variance), Math.sqrt(scale));
    }

    public static double[] nextMultivariateNormalCholesky(double[] mean, double[][] cholesky) {
        return nextMultivariateNormalCholesky(mean, cholesky, 1.0);
    }

    public static double[] nextMultivariateNormalCholesky(double[] mean, double[][] cholesky, double sqrtScale) {

        double[] result = new double[mean.length];
        nextMultivariateNormalCholesky(mean, cholesky, sqrtScale, result);
        return result;
    }

    public static void nextMultivariateNormalCholesky(double[] mean, double[][] cholesky, double sqrtScale, double[] result) {

        final int dim = mean.length;

        System.arraycopy(mean, 0, result, 0, dim);

        double[] epsilon = new double[dim];
        for (int i = 0; i < dim; i++)
            epsilon[i] = MathUtils.nextGaussian() * sqrtScale;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
                result[i] += cholesky[i][j] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
        }
    }

    public static void nextMultivariateNormalCholesky(final WrappedVector mean, final double[][] cholesky,
                                                      final double sqrtScale, final WrappedVector result,
                                                      final double[] epsilon) {

        final int dim = mean.getDim();
        for (int i = 0; i < dim; i++) {
            epsilon[i] = MathUtils.nextGaussian() * sqrtScale;
        }

        for (int i = 0; i < dim; i++) {
            double x = mean.get(i);
            for (int j = 0; j <= i; j++) {
                x += cholesky[i][j] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
            result.set(i, x);
        }
    }

    public static void nextMultivariateNormalCholesky(final double[] mean, final int meanOffset, final double[][] cholesky,
                                                      final double sqrtScale, final double[] result, final int resultOffset,
                                                      final double[] epsilon) {

        final int dim = epsilon.length;

        System.arraycopy(mean, meanOffset, result, resultOffset, dim);

        for (int i = 0; i < dim; i++)
            epsilon[i] = MathUtils.nextGaussian() * sqrtScale;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
                result[resultOffset + i] += cholesky[i][j] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
        }
    }

    // TODO should be a junit test
    public static void main(String[] args) {
        testPdf();
        testRandomDraws();
    }

    public static void testPdf() {
        double[] start = {1, 2};
        double[] stop = {0, 0};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        double scale = 0.2;
        System.err.println("logPDF = " + logPdf(start, stop, precision, Math.log(calculatePrecisionMatrixDeterminate(precision)), scale));
        System.err.println("Should = -19.94863\n");

        System.err.println("logPDF = " + logPdf(start, stop, 2, 0.2));
        System.err.println("Should = -24.53529\n");
    }

    public static void testRandomDraws() {

        double[] start = {1, 2};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        int length = 100000;


        System.err.println("Random draws (via precision) ...");
        double[] mean = new double[2];
        double[] SS = new double[2];
        double[] var = new double[2];
        double ZZ = 0;
        for (int i = 0; i < length; i++) {
            double[] draw = nextMultivariateNormalPrecision(start, precision);
            for (int j = 0; j < 2; j++) {
                mean[j] += draw[j];
                SS[j] += draw[j] * draw[j];
            }
            ZZ += draw[0] * draw[1];
        }

        for (int j = 0; j < 2; j++) {
            mean[j] /= length;
            SS[j] /= length;
            var[j] = SS[j] - mean[j] * mean[j];
        }
        ZZ /= length;
        ZZ -= mean[0] * mean[1];

        System.err.println("Mean: " + new Vector(mean));
        System.err.println("TRUE: [ 1 2 ]\n");
        System.err.println("MVar: " + new Vector(var));
        System.err.println("TRUE: [ 0.571 1.14 ]\n");
        System.err.println("Covv: " + ZZ);
        System.err.println("TRUE: -0.286");
    }

    public static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

    // RandomGenerator interface
    public Object nextRandom() {
        return nextMultivariateNormal();
    }

    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() { return mean.length; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    @Override
    public double[][] getPrecisionMatrix() {
        return precision;
    }
}
