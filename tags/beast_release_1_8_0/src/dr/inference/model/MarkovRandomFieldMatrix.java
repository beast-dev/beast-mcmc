/*
 * MarkovRandomFieldMatrix.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
public class MarkovRandomFieldMatrix extends MatrixParameter {

    private Parameter diagonalParameter;
    private Parameter offDiagonalParameter;

    private boolean asCorrelation = false;

    private int dim;

    public MarkovRandomFieldMatrix(Parameter diagonals, Parameter offDiagonal, boolean asCorrelation) {
        super(MATRIX_PARAMETER);
        diagonalParameter = diagonals;
        offDiagonalParameter = offDiagonal;
        addParameter(diagonalParameter);
        addParameter(offDiagonalParameter);
        dim = diagonalParameter.getDimension();
        this.asCorrelation = asCorrelation;
    }

//    public double[] getAttributeValue() {
//        double[] stats = new double[dim * dim];
//        int index = 0;
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                stats[index] = getParameterValue(i, j);
//                index++;
//            }
//        }
//        return stats;
//    }
//
//    public int getDimension() {
//        return dim * dim;
//    }

    public double getParameterValue(int i) {
        int row = i / dim;
        int col = i % dim;
        return getParameterValue(row, col);
    }

    public double getParameterValue(int row, int col) {
        if (row == col) {
            return diagonalParameter.getParameterValue(row);
        } else if (row == (col - 1) || row == (col + 1)) {
            if (asCorrelation) {
                return -offDiagonalParameter.getParameterValue(0) *
                        Math.sqrt(diagonalParameter.getParameterValue(row) * diagonalParameter.getParameterValue(col));
            }
            return offDiagonalParameter.getParameterValue(0);
        }
        return 0.0;
    }

//    public double[][] getParameterAsMatrix() {
//        final int I = dim;
//        double[][] parameterAsMatrix = new double[I][I];
//        for (int i = 0; i < I; i++) {
//            parameterAsMatrix[i][i] = getParameterValue(i, i);
//            for (int j = i + 1; j < I; j++) {
//                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
//            }
//        }
//        return parameterAsMatrix;
//    }

    public int getColumnDimension() {
        return diagonalParameter.getDimension();
    }

    public int getRowDimension() {
        return diagonalParameter.getDimension();
    }
}
