/*
 * DiagonalMatrix.java
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

package dr.inference.model;

/**
 * @author Marc Suchard
 */
public class DiagonalMatrix extends MatrixParameter
        implements MatrixParameterInterface, VariableListener {

    private Parameter diagonalParameter;

    public DiagonalMatrix(Parameter param) {
        super(MATRIX_PARAMETER);
        addParameter(param);
        diagonalParameter = param;
    }

//    public static DiagonalMatrix buildIdentityTimesElementMatrix(int dim, double value) {
//        Parameter param = new Parameter.Default(dim, value);
//        param.addBounds(new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));
//        return new DiagonalMatrix(param);
//    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == diagonalParameter) {
            fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }

    }

    public String getDimensionName(int dim) {
        int row = dim / diagonalParameter.getDimension();
        int col = dim % diagonalParameter.getDimension();
        return getId() + "." + row + "." + col;
    }

    public Parameter getDiagonalParameter() { return diagonalParameter; }

    public double getParameterValue(int row, int col) {
        if (row != col)
            return 0.0;
        return diagonalParameter.getParameterValue(row);
    }

    public double[][] getParameterAsMatrix() {
        final int I = getRowDimension();
        double[][] parameterAsMatrix = new double[I][I];
        for (int i = 0; i < I; i++) {
            parameterAsMatrix[i][i] = diagonalParameter.getParameterValue(i);
        }
        return parameterAsMatrix;
    }

    public void inverse() {
        for (int i = 0; i < diagonalParameter.getDimension(); i++) {
            {
                if (diagonalParameter.getValue(i) == 0) {
                    throw new RuntimeException("Diagonal matrix is not of full rank");
                } else {
                    setParameterValue(i, 1 / diagonalParameter.getParameterValue(i));
                }
            }
        }
    }

    public int getDimension() {
        return diagonalParameter.getDimension() * diagonalParameter.getDimension();
    }

    public double getParameterValue(int i) {
        return getParameterValue(i / diagonalParameter.getDimension(), i % diagonalParameter.getDimension());
    }

//    public MatrixParameter preMultiply(MatrixParameter right) {
//        if (right.getRowDimension() != this.getColumnDimension()) {
//            throw new RuntimeException("Incompatible Dimensions: " + right.getRowDimension() + " does not equal " + this.getColumnDimension() + ".\n");
//        }
//        MatrixParameter answer = new MatrixParameter(null);
//        answer.setDimensions(right.getRowDimension(), right.getColumnDimension());
//        for (int i = 0; i < right.getRowDimension(); i++) {
//            for (int j = 0; j < right.getColumnDimension(); j++) {
//                answer.setParameterValueQuietly(i, j, right.getParameterValue(i, j) * getParameterValue(i));
//            }
//
//        }
//        return answer;
//    }
//
//    public MatrixParameter preMultiplyInPlace(MatrixParameter right, MatrixParameter answer) {
//        if (right.getRowDimension() != this.getColumnDimension()) {
//            throw new RuntimeException("Incompatible Dimensions: " + right.getRowDimension() + " does not equal " + this.getColumnDimension() + ".\n");
//        }
////        MatrixParameter answer=new MatrixParameter(null);
////        answer.setDimensions(right.getRowDimension(), right.getColumnDimension());
//        for (int i = 0; i < right.getRowDimension(); i++) {
//            for (int j = 0; j < right.getColumnDimension(); j++) {
//                answer.setParameterValueQuietly(i, j, right.getParameterValue(i, j) * getParameterValue(i));
//            }
//
//        }
//        return answer;
//    }
//
//
//    public MatrixParameter postMultiply(MatrixParameter left) {
//        if (left.getColumnDimension() != this.getRowDimension()) {
//            throw new RuntimeException("Incompatible Dimensions: " + this.getColumnDimension() + " does not equal " + left.getRowDimension() + ".\n");
//        }
//        MatrixParameter answer = new MatrixParameter(null);
//        answer.setDimensions(left.getRowDimension(), left.getColumnDimension());
//        for (int i = 0; i < left.getRowDimension(); i++) {
//            for (int j = 0; j < left.getColumnDimension(); j++) {
//                answer.setParameterValueQuietly(i, j, left.getParameterValue(i, j) * getParameterValue(j));
//            }
//
//        }
//        return answer;
//    }
//
//    public MatrixParameter postMultiplyInPlace(MatrixParameter left, MatrixParameter answer) {
//        if (left.getColumnDimension() != this.getRowDimension()) {
//            throw new RuntimeException("Incompatible Dimensions: " + this.getColumnDimension() + " does not equal " + left.getRowDimension() + ".\n");
//        }
////        MatrixParameter answer=new MatrixParameter(null);
////        answer.setDimensions(left.getRowDimension(), left.getColumnDimension());
//        for (int i = 0; i < left.getRowDimension(); i++) {
//            for (int j = 0; j < left.getColumnDimension(); j++) {
//                answer.setParameterValueQuietly(i, j, left.getParameterValue(i, j) * getParameterValue(j));
//            }
//
//        }
//        return answer;
//    }

    public MatrixParameter add(MatrixParameter Right) {

        MatrixParameter answer = new MatrixParameter(null);
        answer.setDimensions(this.getRowDimension(), this.getColumnDimension());
        return addInPlace(Right, answer);
    }

    public MatrixParameter addInPlace(MatrixParameter Right, MatrixParameter answer) {
        if (this.getColumnDimension() != Right.getColumnDimension() || this.getRowDimension() != Right.getRowDimension()) {
            throw new RuntimeException("You cannot add a " + getRowDimension() + " by " + getColumnDimension() + " matrix to a " + Right.getRowDimension() + " by " + Right.getColumnDimension() + " matrix.");
        }

        for (int i = 0; i < this.getRowDimension(); i++) {
            for (int j = 0; j < this.getColumnDimension(); j++) {
                if (i == j) {
                    answer.setParameterValueQuietly(i, j, this.getParameterValue(i, j) + Right.getParameterValue(i, j));
                } else {
                    answer.setParameterValueQuietly(i, j, Right.getParameterValue(i, j));
                }
            }
        }
        return answer;
    }

    public void product(double a) {
        for (int i = 0; i < this.getRowDimension(); i++) {
            this.setParameterValueQuietly(i, i, a * this.getParameterValue(i, i));

        }
        this.fireParameterChangedEvent();
    }

    public MatrixParameter productInPlace(double a, MatrixParameter answer) {
        for (int i = 0; i < this.getRowDimension(); i++) {
            answer.setParameterValueQuietly(i, i, a * this.getParameterValue(i, i));

        }
        return answer;
    }

    public double getDeterminant() {
        double product = 1;
        for (int i = 0; i < diagonalParameter.getDimension(); i++) {
            product *= diagonalParameter.getParameterValue(i);
        }
        return product;
    }

    public void setParameterValue(int i, double j) {
        diagonalParameter.setParameterValue(i, j);
    }

    public void setParameterValueQuietly(int i, double j) {
        diagonalParameter.setParameterValueQuietly(i, j);
    }

    public int getColumnDimension() {
        return diagonalParameter.getDimension();
    }

    public int getRowDimension() {
        return diagonalParameter.getDimension();
    }

}
