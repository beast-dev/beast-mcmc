/*
 * LKJCholeskyTransformConstrained.java
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

package dr.util;

import dr.math.matrixAlgebra.WrappedMatrix;

/**
 * @author Paul Bastide
 */

public class LKJCholeskyTransformConstrained extends Transform.MultivariableTransform {

    // LKJ transform with CPCs constrained between -1 and 1.
    // transform: from cholesky of correlation matrix to constrained CPCs

    protected int dim;

    public LKJCholeskyTransformConstrained(int dim) {
        this.dim = dim;
    }

    public double[] transform(double[] values) {
        return transform(values, 0, values.length);
    }

    public double[] inverse(double[] values) {
        return inverse(values, 0, values.length);
    }

    // values = CPCs
    @Override
    public double[] inverse(double[] values, int from, int to) {
        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";
        assert dim * (dim - 1) / 2 == values.length : "The transform function can only be applied to the whole array of values.";
        for (int k = 0; k < dim; k++) {
            assert values[k] <= 1.0 && values[k] >= -1.0 : "CPCs must be between -1.0 and 1.0";
        }

        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix L
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(dim);
        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix Z
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(values, dim, 1.0);

        double acc;
        double temp;
//        L.set(0, 0, 1.0);
        for (int j = 1; j < dim; j++) {
            acc = 1.0;
            for (int i = 0; i < j; i++) {
                temp = Z.get(i, j);
                L.set(i, j, temp * acc);
                acc *= Math.sqrt(1 - Math.pow(temp, 2));
            }
        }

        return L.getBuffer();
    }

    // values = Cholesky of correlation
    @Override
    public double[] transform(double[] values, int from, int to) {
        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";

        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix L
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(values, dim);

        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix Z
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(dim);

        double acc;
        double temp;
        for (int j = 1; j < dim; j++) {
            acc = 1.0;
            for (int i = 0; i < j; i++) {
                temp = L.get(i, j) / acc;
                Z.set(i, j, temp);
                acc *= Math.sqrt(1 - Math.pow(temp, 2));
            }
        }

        return Z.getBuffer();
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant for the LKJ transform.");
    }

    public String getTransformName() {
        return "LKJTransform";
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        // values = untransformed (R)
        double[] transformedValues = transform(value);
        // Jacobian of inverse
        double[][] jacobianInverse = computeJacobianMatrixInverse(transformedValues);
        // gradient of log jacobian of the inverse
        double[] gradientLogJacobianInverse = getGradientLogJacobianInverse(transformedValues);
        // Matrix multiplication (upper triangular) + updated gradient
        double[] updatedGradient = new double[gradient.length];
        for (int i = 0; i < gradient.length; i++) {
            for (int j = i; j < gradient.length; j++) {
                updatedGradient[i] += jacobianInverse[i][j] * gradient[j];
            }
            updatedGradient[i] += gradientLogJacobianInverse[i];
        }
        return updatedGradient;
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getLogJacobian(double[] values, int from, int to) {
        assert from == 0 && to == values.length
                : "The logJacobian function can only be applied to the whole array of values.";
        double[] transformedValues = transform(values);
        double logJacobian = 0;
        int k = 0;
        for (int i = 0; i < dim - 2; i++) {
            k++;
            for (int j = i + 2; j < dim; j++) {
                logJacobian += (j - i - 1) * Math.log(1.0 - Math.pow(transformedValues[k], 2));
                k++;
            }
        }
        return -0.5 * logJacobian;
    }

    public double[] getGradientLogJacobianInverse(double[] values) {
        double[] gradientLogJacobian = new double[values.length];
        int k = 0;
        for (int i = 0; i < dim - 2; i++) { // Sizes of conditioning sets
            for (int j = i + 1; j < dim; j++) {
                gradientLogJacobian[k] = -(j - i - 1) * values[k] / (1.0 - Math.pow(values[k], 2));
                k++;
            }
        }
        return gradientLogJacobian;
    }
    // ************************************************************************* //
    // Computation of the jacobian matrix
    // ************************************************************************* //

    // Returns the *transpose* of the Jacobian matrix: jacobian[posStrict(k, l)][posStrict(i, j)] = d W_{ij} / d Z_{kl}
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim * (dim - 1) / 2][dim * (dim - 1) / 2];

        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix Z
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(values, dim, 1.0);

        for (int j = 1; j < dim; j++) {
            for (int k = 0; k < j + 1; k++) {
                recursionJacobian(jacobian, Z, k, j);
            }
        }

        return jacobian;
    }

    private void recursionJacobian(double[][] jacobian, WrappedMatrix.WrappedStrictlyUpperTriangularMatrix Z,
                                   int k, int j) {
        // 0 <= k <= j

        WrappedMatrix.WrappedStrictlyUpperTriangularMatrix L
                = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(jacobian[posStrict(k, j)], dim);

        double acc;
        double temp;
        acc = 1.0;
        // before k
        for (int i = 0; i < k; i++) {
            acc *= Math.sqrt(1 - Math.pow(Z.get(i, j), 2));
        }
        // k
        if (k < j) L.set(k, j, acc);
        temp = Z.get(k, j);
        acc *= -temp / Math.sqrt(1 - Math.pow(temp, 2));
        // after k
        for (int i = k + 1; i < j; i++) {
            temp = Z.get(i, j);
            L.set(i, j, temp * acc);
            acc *= Math.sqrt(1 - Math.pow(temp, 2));
        }
    }

    // ************************************************************************* //
    // Helper functions to deal with upper triangular matrices
    // ************************************************************************* //

    private int posStrict(int i, int j) {
        return i * (2 * dim - i - 1) / 2 + (j - i - 1);
    }
}
