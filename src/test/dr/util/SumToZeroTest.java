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
import dr.util.SumToZeroTransform;
import dr.util.UnitSimplexViaSoftmaxToRealsTransform;
import test.dr.math.MathTestCase;

/**
 * @author Marc Suchard
 */

public class SumToZeroTest extends MathTestCase {

    public void testTransformation() {

        System.out.println("\nTest sum-to-zero transform.");

        double[] unconstrainedValues = new double[]{ 1.0, 1.0, 1.0 };
        double[] sumToZeroValues = new double[]{ 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

        SumToZeroTransform transform = new SumToZeroTransform(unconstrainedValues.length + 1);
        double[] result1 = transform.inverse(unconstrainedValues, 0, unconstrainedValues.length);

//        assertEquals(result1, sumToZeroValues, 1E-10);

        System.out.println("Success");
    }

    public void testJacobianDeterminant() {
        System.out.println("\nTest determinant");

        double[] unconstrainedValues = new double[]{ 1, -1, 1, 1 };
        int outDim = unconstrainedValues.length;
        int inDim = outDim + 1;

        SumToZeroTransform transform = new SumToZeroTransform(inDim);

        double logDet = transform.getLogJacobian(null);

        // TODO

        System.out.println("Success");
    }

    public void testJacobianMatrix() {

        System.out.println("\nTest Jacobian matrix");

        double[] unconstrainedValues = new double[]{ 1, -1, 1, 1 };
        int outDim = unconstrainedValues.length;
        int inDim = outDim + 1;

        SumToZeroTransform transform = new SumToZeroTransform(inDim);

        double[][] numericalJacobian = new double[inDim][outDim];
        for (int i = 0; i < inDim; ++i) {

            final int index = i;
            double[] tmp = NumericalDerivative.gradient(new MultivariateFunction() {
                @Override
                public double evaluate(double[] argument) {
                    return transform.inverse(argument, 0, argument.length)[index];
                }

                @Override
                public int getNumArguments() {
                    return outDim;
                }

                @Override
                public double getLowerBound(int n) {
                    return Double.NEGATIVE_INFINITY;
                }

                @Override
                public double getUpperBound(int n) {
                    return Double.POSITIVE_INFINITY;
                }
            }, unconstrainedValues);

            numericalJacobian[i] = new double[outDim];
            for (int j = 0; j < tmp.length; ++j) {
                numericalJacobian[i][j] = tmp[j];
            }
        }
        numericalJacobian = transpose(numericalJacobian);

        double[][] jacobian = transform.computeJacobianMatrixInverse(unconstrainedValues);

        for (int i = 0; i < jacobian.length; ++i) {
            assertEquals(jacobian[i], numericalJacobian[i], 1E-5);
        }

        System.out.println("Success");
    }

    private static double[][] transpose(double[][] matrix) {

        if (matrix.length == matrix[0].length) {
            for (int i = 0; i < matrix.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    double tmp = matrix[i][j];
                    matrix[i][j] = matrix[j][i];
                    matrix[j][i] = tmp;
                }
            }
            return matrix;
        } else {
            double[][] result = new double[matrix[0].length][matrix.length];
            for (int i = 0; i < result.length; ++i) {
                for (int j = 0; j < result[i].length; ++j) {
                    result[i][j] = matrix[j][i];
                }
            }
            return result;
        }
    }
}
