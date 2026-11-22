/*
 * CorrelationToCholesky.java
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
import dr.xml.*;

import java.util.Arrays;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.extractUpperTriangular;

/**
 * @author Paul Bastide
 */

public class CorrelationToCholesky extends Transform.MultivariateTransform {

    // Transform a correlation matrix into a Cholesky matrix

    private int dimVector;
    private final ThreadLocal<double[]> diagonalScratch;

    public CorrelationToCholesky(int dimVector) {
        super(dimVector * (dimVector - 1) / 2);
        this.dimVector = dimVector;
        this.diagonalScratch = ThreadLocal.withInitial(() -> new double[dimVector]);
    }

    // values = cholesky
    @Override
    protected double[] inverse(double[] values) {
        WrappedMatrix.WrappedUpperTriangularMatrix L = WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal(values, dimVector);

        SymmetricMatrix R = L.transposedProduct();

        return extractUpperTriangular(R);
    }

    // values = correlation
    @Override
    protected double[] transform(double[] values) {
        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(values, dimVector);
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

    @Override
    public boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    public String getTransformName() {
        return "CorrelationToCholeskyTransform";
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        WrappedMatrix.WrappedUpperTriangularMatrix L = WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal(transform(values), dimVector);
        double logJacobian = 0;
        for (int i = 0; i < dimVector - 1; i++) {
            logJacobian += (dimVector - i - 1) * Math.log(L.get(i, i));
        }
        return -logJacobian;
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        double[] diagonal = diagonalScratch.get();
        fillDiagonal(values, diagonal);
        double[] gradientLogJacobian = new double[values.length];
        int k = 0;
        for (int i = 0; i < dimVector - 1; i++) {
            for (int j = i + 1; j < dimVector; j++) {
                gradientLogJacobian[k] = -(dimVector - j - 1) * getStrictUpper(values, i, j) / Math.pow(diagonal[j], 2);
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
        double[][] jacobian = new double[dim][dim];

        WrappedMatrix.WrappedUpperTriangularMatrix W = WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal(values, dimVector);

        for (int i = 0; i < dimVector - 1; i++) {
            for (int j = i + 1; j < dimVector; j++) {
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

    @Override
    protected double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value) {
        double[] updatedGradient = new double[gradient.length];
        updateGradientInverseUnWeightedLogDensity(gradient, value, updatedGradient);
        return updatedGradient;
    }

    public void updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, double[] result) {
        updateGradientInverseUnWeightedLogDensity(gradient, value, result, diagonalScratch.get());
    }

    public void updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, double[] result, double[] diagonal) {
        assert gradient.length == dim && value.length == dim && result.length == dim;
        assert diagonal.length >= dimVector;

        Arrays.fill(result, 0.0);
        fillDiagonal(value, diagonal);

        for (int i = 0; i < dimVector - 1; i++) {
            double wii = diagonal[i];
            for (int j = i + 1; j < dimVector; j++) {
                int ij = posStrict(i, j);
                double gradientIJ = gradient[ij];
                double temp = getStrictUpper(value, i, j) / wii;
                for (int k = 0; k < i; k++) {
                    double wki = getStrictUpper(value, k, i);
                    result[posStrict(k, i)] += (getStrictUpper(value, k, j) - wki * temp) * gradientIJ;
                    result[posStrict(k, j)] += wki * gradientIJ;
                }
                result[ij] += wii * gradientIJ;
            }
        }
    }

    // ************************************************************************* //
    // Helper functions to deal with upper triangular matrices
    // ************************************************************************* //

    private int posStrict(int i, int j) {
        return i * (2 * dimVector - i - 1) / 2 + (j - i - 1);
    }

    private double getStrictUpper(double[] values, int i, int j) {
        return values[posStrict(i, j)];
    }

    private void fillDiagonal(double[] values, double[] diagonal) {
        for (int j = 0; j < dimVector; j++) {
            double sum = 0.0;
            for (int i = 0; i < j; i++) {
                double value = getStrictUpper(values, i, j);
                sum += value * value;
            }
            if (sum > 1.0) {
                if (Math.abs(sum - 1.0) > 1E-6) {
                    throw new RuntimeException("Values are not consistent with the cholesky decomposition of " +
                            "a correlation matrix. Sum of squared values must be less than 1 (got " + sum + ")");
                }
                sum = 1.0;
            }
            diagonal[j] = Math.sqrt(1 - sum);
        }
    }


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String DIMENSION = "dimension";
        private static final String CORRELATION_TO_CHOLESKY = "correlationToCholeskyTransform";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int dim = xo.getIntegerAttribute(DIMENSION);
            return new CorrelationToCholesky(dim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION)
            };
        }

        @Override
        public String getParserDescription() {
            return "transforms the off-diagonal elements of a correlation to the off-diagonal elements of its" +
                    " Cholesky decomposition";
        }

        @Override
        public Class getReturnType() {
            return CorrelationToCholesky.class;
        }

        @Override
        public String getParserName() {
            return CORRELATION_TO_CHOLESKY;
        }
    };

}
