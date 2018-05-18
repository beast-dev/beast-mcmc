/*
 * CorrelationFromCholesky.java
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

import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.extractUpperTriangular;
import static dr.math.matrixAlgebra.WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal;

/**
 * @author Paul Bastide
 */

public class CorrelationToCholesky extends Transform.MultivariateTransform {

    // Transform a correlation matrix into a Cholesky matrix

    protected int dim;

    public CorrelationToCholesky(int dim) {
        this.dim = dim;
    }

    // values = cholesky
    @Override
    public double[] inverse(double[] values, int from, int to) {
        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";
        assert dim * (dim - 1) / 2 == values.length : "The transform function can only be applied to the whole array of values.";

        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(values, dim);

        SymmetricMatrix R = L.transposedProduct();

        return extractUpperTriangular(R);
    }

    // values = correlation
    @Override
    public double[] transform(double[] values, int from, int to) {
        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";

        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(values, dim);
        double[] L;
        try {
            L = (new CholeskyDecomposition(R)).getStrictlyUpperTriangular();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in LKJ inverse transform.");
        }

        return L;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant for the correlation to Cholesky transform.");
    }

    public String getTransformName() {
        return "CorrelationToCholeskyTransform";
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getLogJacobian(double[] values, int from, int to) {
        assert from == 0 && to == values.length
                : "The logJacobian function can only be applied to the whole array of values.";
        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(transform(values), dim);
        double logJacobian = 0;
        for (int i = 0; i < dim - 1; i++) {
            logJacobian += (dim - i - 1) * Math.log(L.get(i, i));
        }
        return -logJacobian;
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(values, dim);
        double[] gradientLogJacobian = new double[values.length];
        int k = 0;
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                gradientLogJacobian[k] = -(dim - j - 1) * L.get(i, j) / Math.pow(L.get(j, j), 2);
                k++;
            }
        }
        return gradientLogJacobian;
    }

    // ************************************************************************* //
    // Computation of the jacobian matrix
    // ************************************************************************* //

    // Returns the *transpose* of the Jacobian matrix: jacobian[posStrict(k, l)][posStrict(i, j)] = d R_{ij} / d V_{kl}
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim * (dim - 1) / 2][dim * (dim - 1) / 2];

        WrappedMatrix.WrappedUpperTriangularMatrix W = fillDiagonal(values, dim);

        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                double temp = W.get(i, j) / W.get(i, i);
                for (int k = 0; k < i; k++) {
                    jacobian[posStrict(k, i)][posStrict(i, j)] = W.get(k, j) - W.get(k, i) * temp;
                    jacobian[posStrict(k, j)][posStrict(i, j)] = W.get(k, i);
                }
                jacobian[posStrict(i, j)][posStrict(i, j)] = W.get(i, i);
            }
        }

        return jacobian;
    }

    // ************************************************************************* //
    // Helper functions to deal with upper triangular matrices
    // ************************************************************************* //

    private int posStrict(int i, int j) {
        return i * (2 * dim - i - 1) / 2 + (j - i - 1);
    }

}

