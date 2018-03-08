/*
 * CompoundSymmetricMatrix.java
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
 * @author Zhenyu Zhang
 */
public class CorrelationSymmetricMatrix extends MatrixParameter {

    public enum Type {
        AS_CORRELATION,
        AS_IS
    }

    private final Parameter diagonals;
    private final Parameter offDiagonals;
    private final Type type;
    private final int dim;

    public CorrelationSymmetricMatrix(Parameter diagonals, Parameter offDiagonals) {
        this(diagonals, offDiagonals, Type.AS_CORRELATION);
    }

    public CorrelationSymmetricMatrix(Parameter diagonals, Parameter offDiagonals, Type type) {
        super(MATRIX_PARAMETER);
        this.diagonals = diagonals;
        this.offDiagonals = offDiagonals;
        addParameter(diagonals);
        addParameter(offDiagonals);
        this.type = type;
        this.dim = diagonals.getDimension();

        if (offDiagonals.getDimension() !=
                (diagonals.getDimension() * (diagonals.getDimension() - 1) / 2)) {
            throw new IllegalArgumentException("Invalid parameter dimensions");
        }
    }

    public double[] getAttributeValue() {
        double[] stats = new double[dim * dim];
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                stats[index] = getParameterValue(i, j);
                index++;
            }
        }        
        return stats;
    }

    @Override
    public double getParameterValue(int row, int col) {
        if (row != col) {
            double value = offDiagonals.getParameterValue(getOffDiagonalIndex(row, col, dim));
            if (type == Type.AS_CORRELATION) {
                value *= Math.sqrt(diagonals.getParameterValue(row) * diagonals.getParameterValue(col));
            }
            return value;
        } else {
            return diagonals.getParameterValue(row);
        }
    }

    @Override
    public void setParameterValue(int row, int column, double a) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int column, double a){
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int column, double val){
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {
        double[][] parameterAsMatrix = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            parameterAsMatrix[i][i] = getParameterValue(i, i);
            for (int j = i + 1; j < dim; j++) {
                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
            }
        }
        return parameterAsMatrix;
    }

    @Override
    public int getColumnDimension() {
        return diagonals.getDimension();
    }

    @Override
    public int getRowDimension() {
        return diagonals.getDimension();
    }

    private static int getOffDiagonalIndex(int row, int col, final int dim) {
        if (col < row) {
            int tmp = row;
            row = col;
            col = tmp;
        }
        int index = 0;
        for (int r = 0; r < row; ++r) {
            index += dim - 1 - r;
        }
        index += col - row - 1;

        return index;
    }
}
