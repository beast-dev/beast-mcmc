/*
 * CorrelationDiagonalMatrix.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 * @author Marc A. Suchard
 */

public class DiagonalConstrainedMatrixView extends CompoundParameter implements MatrixParameterInterface {

    private final MatrixParameterInterface matrix;
    private final Parameter mask;
    private final double constraintValue;
    private final double sqrtConstraintValue;


    public DiagonalConstrainedMatrixView(String id,
                                         MatrixParameterInterface matrix,
                                         Parameter mask,
                                         double constraintValue) {
        super(id);
        this.matrix = matrix;
        this.mask = mask;
        this.constraintValue = constraintValue;
        this.sqrtConstraintValue = Math.sqrt(constraintValue);

        matrix.addVariableListener(this);
    }

    private boolean constrained(int i) {
        return mask.getParameterValue(i) == 1.0;
    }

    @Override
    public double getParameterValue(int row, int col) {

        if (row == col && constrained(row)) {
            return constraintValue;
        }

        double entry = matrix.getParameterValue(row, col);

        if (constrained(row)) {
            entry /= Math.sqrt(matrix.getParameterValue(row, row));
            entry *= sqrtConstraintValue;
        }

        if (constrained(col)) {
            entry /= Math.sqrt(matrix.getParameterValue(col, col));
            entry *= sqrtConstraintValue;
        }

        return entry;
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {
        final int rowDim = getRowDimension();
        final int colDim = getColumnDimension();

        double[][] matrix = new double[rowDim][colDim];
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                matrix[i][j] = getParameterValue(i, j);
            }
        }
        return matrix;
    }

    @Override
    public int getDimension() {
        return matrix.getDimension();
    }

    public double getParameterValue(int dim) {
        final int row = dim / getRowDimension();
        final int col = dim % getRowDimension();
        return getParameterValue(row, col);
    }

    @Override
    public String getDimensionName(int i) {
        return matrix.getDimensionName(i);
    }

    @Override
    public int getColumnDimension() {
        return matrix.getColumnDimension();
    }

    @Override
    public int getRowDimension() {
        return matrix.getRowDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == matrix) {
            fireParameterChangedEvent(index, type);
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }
}
