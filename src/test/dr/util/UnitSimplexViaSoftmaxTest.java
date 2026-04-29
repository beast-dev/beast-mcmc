/*
 * LKJTransformTest.java
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

package test.dr.util;

import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.util.UnitSimplexToRealsTransform;
import dr.util.UnitSimplexViaSoftmaxToRealsTransform;
import test.dr.math.MathTestCase;

/**
 * @author Marc Suchard
 */

public class UnitSimplexViaSoftmaxTest extends MathTestCase {

    public void testGetGradientLogDetJacobian() {

        System.out.println("\nTest gradient of log-det Jacobian");
        double[] valuesOnReals = new double[] { -1.0, 1.0, 3.0, 2.0 };
        int dim = valuesOnReals.length;

        UnitSimplexViaSoftmaxToRealsTransform transform = new UnitSimplexViaSoftmaxToRealsTransform(dim);
        double[] numericalGradient = NumericalDerivative.gradient(new MultivariateFunction() {
            @Override
            public double evaluate(double[] argument) {
                return -transform.getLogJacobian(transform.inverse(valuesOnReals, 0, dim));
            }

            @Override
            public int getNumArguments() { return dim; }

            @Override
            public double getLowerBound(int n) { return Double.NEGATIVE_INFINITY; }

            @Override
            public double getUpperBound(int n) { return Double.POSITIVE_INFINITY; }
        }, valuesOnReals);
        double[] gradient = transform.getGradientLogJacobianInverse(valuesOnReals);

        assertEquals(gradient, numericalGradient, 1E-6);

        System.out.println("Success");
    }

    public void testTransformation() {

        System.out.println("\nTest Unit-Simplex transform.");

        double[] valuesOnReals = new double[]{ -1.0, 1.0, Double.NaN };
        double[] valuesOnSimplex = new double[]{ 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

        UnitSimplexViaSoftmaxToRealsTransform transform = new UnitSimplexViaSoftmaxToRealsTransform(valuesOnReals.length);
        double[] result1 = transform.inverse(valuesOnReals, 0, valuesOnReals.length);

        assertEquals(result1, valuesOnSimplex, 1E-10);

        System.out.println("Success");
    }

//    public void testGetLogDetJacobian() {
//
//        System.out.println("\nTest log-det Jacobian");
//        double[] valuesOnReals = new double[] { 1.0, 1.0, 2.0, 2.0 };
//        int dim = valuesOnReals.length;
//
//        UnitSimplexViaSoftmaxToRealsTransform transform = new UnitSimplexViaSoftmaxToRealsTransform(dim);
//        double[] valuesOSimplex = transform.inverse(valuesOnReals, 0, dim);
//
//        double logDet1 = transform.getLogJacobian(valuesOSimplex);
//        double[][] tmp = transform.computeJacobianMatrixInverse(valuesOnReals);
////        double[][] jacobian = new double[tmp.length + 1][];
////        System.arraycopy(tmp, 0, jacobian, 0, tmp.length);
////        jacobian[tmp.length] = new double[tmp.length + 1];
////        jacobian[tmp.length][tmp.length] = 1.0;
////        double[][] jacobian = tmp;
//        double[][] jacobian = new double[tmp.length][tmp.length];
//        for (int i = 0; i < tmp.length; ++i) {
//            System.arraycopy(tmp[i], 0, jacobian[i], 0, tmp.length);
//        }
//
//        double logDet2 = Double.NaN;
//        double logDet3 = 0.0;
//
//
//
//        try {
//            CholeskyDecomposition decomp = new CholeskyDecomposition(jacobian);
//            double[][] L = decomp.getL();
//            for (int i = 0; i < 3; ++i) {
//                logDet3 += Math.log(L[i][i]);
//            }
//            logDet3 *= 2;
//
//
//            logDet2 = new Matrix(jacobian).logDeterminant();
//        } catch (IllegalDimension e) {
//            throw new RuntimeException("Illegal dimensions");
//        }
//
//        assertEquals(logDet1, logDet2, 1E-5);
//
//        System.out.println("Success");
//    }

    public void testJacobianMatrix() {

        System.out.println("\nTest Jacobian matrix");

        double[] valuesOnReals = new double[]{ -1.0, 0.0, 0.0, 1.0 };
        int dim = valuesOnReals.length;

        UnitSimplexViaSoftmaxToRealsTransform transform = new UnitSimplexViaSoftmaxToRealsTransform(dim);

        double[][] numericalJacobian = new double[dim][];
        for (int i = 0; i < dim; ++i) {

            final int index = i;
            numericalJacobian[i] = NumericalDerivative.gradient(new MultivariateFunction() {
                @Override
                public double evaluate(double[] argument) {
                    return transform.inverse(argument, 0, argument.length)[index];
                }

                @Override
                public int getNumArguments() {
                    return valuesOnReals.length;
                }

                @Override
                public double getLowerBound(int n) {
                    return Double.NEGATIVE_INFINITY;
                }

                @Override
                public double getUpperBound(int n) {
                    return Double.POSITIVE_INFINITY;
                }
            }, valuesOnReals);
        }
        transpose(numericalJacobian);

        double[][] jacobian = transform.computeJacobianMatrixInverse(valuesOnReals);

        for (int i = 0; i < dim - 1; ++i) {
            assertEquals(jacobian[i], numericalJacobian[i], 1E-5);
        }

        System.out.println("Success");
    }

    private static void transpose(double[][] matrix) {
        assert matrix.length == matrix[0].length;

        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < i; ++j) {
                double tmp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = tmp;
            }
        }
    }
}
