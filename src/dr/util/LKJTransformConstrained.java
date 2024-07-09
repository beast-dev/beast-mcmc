/*
 * LKJTransformConstrained.java
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

public class LKJTransformConstrained extends LKJCholeskyTransformConstrained {

    private static boolean DEBUG = false;

    // LKJ transform with CPCs constrained between -1 and 1.
    // transform: from correlation matrix to constrained CPCs

    public LKJTransformConstrained(int dimVector) {
        super(dimVector);
    }

    // "Cholesky" transform from Stan manual
    @Override
    protected double[] inverse(double[] values) {
        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(super.inverse(values), dimVector);

        SymmetricMatrix R = L.transposedProduct();

        if (DEBUG) {
            System.err.println("Z: " + compoundCorrelationSymmetricMatrix(values, dimVector));
            System.err.println("R: " + R);
            try {
                if (!R.isPD()) {
                    throw new RuntimeException("The LKJ transform should produce a Positive Definite matrix.");
                }
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
        }

        return extractUpperTriangular(R);
    }

    @Override
    protected double[] transform(double[] values) {
        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(values, dimVector);
        double[] L;
        try {
            L = (new CholeskyDecomposition(R)).getStrictlyUpperTriangular();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in LKJ inverse transform.");
        }

        double[] results = super.transform(L);

        if (DEBUG) {
            System.err.println("R: " + compoundCorrelationSymmetricMatrix(values, dimVector));
            System.err.println("L: " + new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(L, dimVector));
            System.err.println("Z: " + compoundCorrelationSymmetricMatrix(results, dimVector));
        }

        return results;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant for the LKJ transform.");
    }

    public String getTransformName() {
        return "LKJTransform";
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        double[] transformedValues = transform(values);
        double logJacobian = 0;
        int k = 0;
        for (int i = 0; i < dimVector - 2; i++) { // Sizes of conditioning sets
            for (int j = i + 1; j < dimVector; j++) {
                logJacobian += (dimVector - i - 2) * Math.log(1.0 - Math.pow(transformedValues[k], 2));
                k++;
            }
        }
        return -0.5 * logJacobian;
    }

    public double[] getGradientLogJacobianInverse(double[] values) {
        double[] gradientLogJacobian = new double[values.length];
        int k = 0;
        for (int i = 0; i < dimVector - 2; i++) { // Sizes of conditioning sets
            for (int j = i + 1; j < dimVector; j++) {
                gradientLogJacobian[k] = -(dimVector - i - 2) * values[k] / (1.0 - Math.pow(values[k], 2));
                k++;
            }
        }
        return gradientLogJacobian;
    }
    // ************************************************************************* //
    // Computation of the jacobian matrix
    // ************************************************************************* //

    // Returns the *transpose* of the Jacobian matrix: jacobian[pos(k, l)][pos(i, j)] = d R_{ij} / d Z_{kl}
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];

        // Initialization
        for (int j = 1; j < dimVector; j++) {
            jacobian[pos(0, j)][j - 1] = 1.0;
        }

        recursionJacobian(jacobian, values);

        return jacobian;
    }

    private void recursionJacobian(double[][] jacobian, double[] values) {
        for (int i = 1; i < dimVector - 1; i++) {
            for (int j = i + 1; j < dimVector; j++) {
                jacobian[pos(i, j)][pos(i, j)] = values[pos(i, j)];
                for (int iota = 1; iota < i + 1; iota++) {
                    setUpperTriangular(jacobian[pos(i, j)], i, j, recursionFormulaJacobian(jacobian[pos(i, j)], values, i, j, iota, i, j));
                }
                // jacobian[pos(k, l)][pos(i, j)] = d R_{ij} / d Z_{kl}
                for (int k = 0; k < i; k++) {
                    jacobian[pos(k, i)][pos(i, j)] = values[pos(i, j)];
                    jacobian[pos(k, j)][pos(i, j)] = values[pos(i, j)];
                    for (int iota = 1; iota < i + 1; iota++) {
                        setUpperTriangular(jacobian[pos(k, i)], i, j, recursionFormulaJacobian(jacobian[pos(k, i)], values, i, j, iota, k, i));
                        setUpperTriangular(jacobian[pos(k, j)], i, j, recursionFormulaJacobian(jacobian[pos(k, j)], values, i, j, iota, k, j));
                    }
                }
            }
        }
    }

    private double recursionFormulaJacobian(double[] trans, double[] values, int i, int j, int iota, int k, int l) {
        double Rimi = getUpperTriangular(values, i - iota, i);
        double Rimj = getUpperTriangular(values, i - iota, j);

        if ((i == k) && (j == l) && (iota == 1)) {
            return Math.sqrt((1 - Rimi * Rimi) * (1 - Rimj * Rimj));
        }
        if ((i - iota == k) && (i == l)) {
            return getUpperTriangular(trans, i, j) * (-Rimi / Math.sqrt((1 - Rimi * Rimi))) * Math.sqrt((1 - Rimj * Rimj)) + Rimj;
        }
        if ((i - iota == k) && (j == l)) {
            return getUpperTriangular(trans, i, j) * (-Rimj / Math.sqrt((1 - Rimj * Rimj))) * Math.sqrt((1 - Rimi * Rimi)) + Rimi;
        }
        if ( i - iota < k ) {
            return getUpperTriangular(trans, i, j) * Math.sqrt((1 - Rimi * Rimi) * (1 - Rimj * Rimj));
        }
        return getUpperTriangular(trans, i, j) * Math.sqrt((1 - Rimi * Rimi) * (1 - Rimj * Rimj)) + Rimi * Rimj;
    }

    // ************************************************************************* //
    // Helper functions to deal with upper triangular matrices
    // ************************************************************************* //

    private double getUpperTriangular(double[] values, int i, int j) {
        assert i <= j;
        if (i == j) return 1.0; // Z_{ii} = 1 by convention.
        return values[pos(i, j)];
    }

    private void setUpperTriangular(double[] values, int i, int j, double val) {
        assert i < j;
        values[pos(i, j)] = val;
    }

    private int pos(int i, int j) {
        return i * (2 * dimVector - i - 1) / 2 + (j - i - 1);
    }

    // ************************************************************************* //
    // Unused recursive transformation (kept for test purposes)
    // ************************************************************************* //

    public double[] inverseRecursion(double[] values, int from, int to) {
        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";
        assert dimVector * (dimVector - 1) / 2 == values.length : "The transform function can only be applied to the whole array of values.";
        for (int k = 0; k < dim; k++) {
            assert values[k] <= 1.0 && values[k] >= -1.0 : "CPCs must be between -1.0 and 1.0";
        }

        double[] trans = new double[values.length];
        System.arraycopy(values, 0, trans, 0, values.length);

        recursionInverse(trans, values);

        return trans;
    }

    public double[] transformRecursion(double[] values, int from, int to) {
        assert from == 0 && to == values.length
                : "The transform function can only be applied to the whole array of values.";

        double[] inv = new double[values.length];
        System.arraycopy(values, 0, inv, 0, values.length);

        recursion(inv);

        if (DEBUG) {
            SymmetricMatrix Z = compoundCorrelationSymmetricMatrix(values, dimVector);
            try {
                if (!Z.isPD()) {
                    throw new RuntimeException("The LKJ transform should produce a Positive Definite matrix.");
                }
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
        }

        return inv;
    }

    private void recursionInverse(double[] trans, double[] values) {
        for (int i = 1; i < dimVector; i++) {
            for (int j = i + 1; j < dimVector; j++) {
                for (int iota = 1; iota < i + 1; iota++) {
                    setUpperTriangular(trans, i, j, recursionInverseFormula(trans, values, i, j, iota));
                }
            }
        }
    }

    private double recursionInverseFormula(double[] trans, double[] values, int i, int j, int iota) {
        double Rimi = getUpperTriangular(values, i - iota, i);
        double Rimj = getUpperTriangular(values, i - iota, j);

        return getUpperTriangular(trans, i, j) * Math.sqrt((1 - Rimi * Rimi) * (1 - Rimj * Rimj)) + Rimi * Rimj;
    }

    private void recursion(double[] inv) {
        for (int i = 1; i < dimVector; i++) {
            for (int j = i + 1; j < dimVector; j++) {
                for (int iota = 1; iota < i + 1; iota++) {
                    setUpperTriangular(inv, i, j, recursionFormula(inv, i, j, iota));
                }
            }
        }
    }

    private double recursionFormula(double[] inv, int i, int j, int iota) {
        double Zimi = getUpperTriangular(inv, iota - 1, i);
        double Zimj = getUpperTriangular(inv, iota - 1, j);

        return (getUpperTriangular(inv, i, j) - Zimi * Zimj) / Math.sqrt((1 - Zimi * Zimi) * (1 - Zimj * Zimj));
    }
}