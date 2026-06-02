/*
 * LKJCholeskyTransformConstrainedWithDiag.java
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


/**
 * @author Zhenyu Zhang
 */
public class LKJCholeskyTransformConstrainedWithDiag extends LKJCholeskyTransformConstrained {

    private int dimCPC;
    private final int totalDimension;

    LKJCholeskyTransformConstrainedWithDiag(int dim) {
        super(dim);
        this.dimCPC = dim * (dim - 1) / 2;
        this.totalDimension = dim * (dim + 1) / 2;
    }

    @Override
    protected double[] transform(double[] values) {

        double[] choleskyFactor = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);
        double[] CPCs = super.transform(choleskyFactor, 0, dimCPC);

        return pasteTogether(CPCs, diagonals);
    }

    @Override
    //values: CPCs appended with log-transformed diagonals
    protected double[] inverse(double[] values) {

        assert values.length == totalDimension : "The transform function can only be applied to CPCs appended " +
                "with diagonals";
        for (int k = 0; k < dimCPC; k++) {
            assert values[k] <= 1.0 && values[k] >= -1.0 : "CPCs must be between -1.0 and 1.0";
        }

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);
        double[] choleskyFactor = super.inverse(CPCs, 0, CPCs.length);
        return pasteTogether(choleskyFactor, diagonals);
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[][] jacobian = super.computeJacobianMatrixInverse(CPCs);
        return appendIdentityMatrix(jacobian);
    }

    public double[] getGradientLogJacobianInverse(double[] values) {

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[] gradientLogJacobianInverse = super.getGradientLogJacobianInverse(CPCs);
        return pasteTogether(gradientLogJacobianInverse, new double[dimVector]);
    }

    private double[] subsetCholeskyOrCPCs(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == totalDimension;
        double[] choleskyOrCPC = new double[dimCPC];
        System.arraycopy(values, dimVector, choleskyOrCPC, 0, dimCPC);
        return choleskyOrCPC;
    }

    private double[] subsetDiagonals(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == totalDimension;
        double[] diagonals = new double[dimVector];
        System.arraycopy(values, 0, diagonals, 0, dimVector);
        return diagonals;
    }

    private double[][] appendIdentityMatrix(double[][] jacobian) {

        assert jacobian.length == dimCPC;
        int length = dimCPC + dimVector;
        double[][] appendedJacobian = new double[length][length];

        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i >= dimCPC || j >= dimCPC) {
                    if (i == j) {
                        appendedJacobian[i][j] = 1;
                    }
                } else {
                    appendedJacobian[i][j] = jacobian[i][j];
                }
            }
        }

        return appendedJacobian;
    }

    private double[] pasteTogether(double[] choleskyOrCPCs, double[] diagonals) {

        double[] concatenatedArray = new double[dimCPC + dimVector];
        System.arraycopy(diagonals, 0, concatenatedArray, 0, dimVector);
        System.arraycopy(choleskyOrCPCs, 0, concatenatedArray, dimVector, dimCPC);
        return concatenatedArray;
    }

    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }
}
