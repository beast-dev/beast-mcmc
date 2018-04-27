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
 */
public class CompoundSymmetricMatrix extends MatrixParameter {

    private final Parameter diagonalParameter;
    private final Parameter offDiagonalParameter;

    private boolean asCorrelation = false;
    private int dim;

    public CompoundSymmetricMatrix(Parameter diagonals, Parameter offDiagonal, boolean asCorrelation) {
        super(MATRIX_PARAMETER);
        diagonalParameter = diagonals;
        offDiagonalParameter = offDiagonal;
        addParameter(diagonalParameter);
        addParameter(offDiagonal);
        dim = diagonalParameter.getDimension();
        this.asCorrelation = asCorrelation;
    }

    @Override
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
    public int getDimension() {
        return getColumnDimension() * getRowDimension();
    }
    
    @Override
    public double getParameterValue(int index) {
        final int dim = getColumnDimension();
        return getParameterValue(index / dim, index % dim);
    }

    @Override
    public String getDimensionName(int index) {
        int dim = getColumnDimension();
        String row = Integer.toString(index / dim);
        String col = Integer.toString(index % dim);

        return getId() + row + col;
    }

    @Override
    public double getParameterValue(int row, int col) {
        if (row != col) {
            if (asCorrelation) {
                return offDiagonalParameter.getParameterValue(getUpperTriangularIndex(row, col)) *
                        Math.sqrt(diagonalParameter.getParameterValue(row) * diagonalParameter.getParameterValue(col));
            }
            return offDiagonalParameter.getParameterValue(0);
        }
        return diagonalParameter.getParameterValue(row);
    }

    private int getUpperTriangularIndex(int i, int j) {
        assert i != j;
        if (i < j) {
            return upperTriangularTransformation(i, j);
        } else {
            return upperTriangularTransformation(j, i);
        }
    }

    private int upperTriangularTransformation(int i, int j) {
        return i * (2 * dim - i - 1) / 2 + (j - i - 1);
    }

    @Override
    public double[][] getParameterAsMatrix() {
        final int I = dim;
        double[][] parameterAsMatrix = new double[I][I];
        for (int i = 0; i < I; i++) {
            parameterAsMatrix[i][i] = getParameterValue(i, i);
            for (int j = i + 1; j < I; j++) {
                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
            }
        }
        return parameterAsMatrix;
    }

    @Override
    public int getColumnDimension() {
        return diagonalParameter.getDimension();
    }

    @Override
    public int getRowDimension() {
        return diagonalParameter.getDimension();
    }

    @Override
    public void setParameterValue(int index, double a) {
        throw new RuntimeException("Do not set entries of a CompoundSymmetricMatrix directly");
    }

    @Override
    public void setParameterValue(int row, int column, double a) {
        throw new RuntimeException("Do not set entries of a CompoundSymmetricMatrix directly");
    }
}
