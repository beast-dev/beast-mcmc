/*
 * DiagonalMatrix.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
public class DiagonalMatrix extends MatrixParameter {

    private Parameter diagonalParameter;

    public DiagonalMatrix(Parameter param) {
        super(MATRIX_PARAMETER);
        addParameter(param);
        diagonalParameter = param;
    }

//	public DiagonalMatrix(String name, Parameter parameter) {
//		Parameter.Default(name, parameters);
//	}

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

    public int getDimension() {
        return diagonalParameter.getDimension() * diagonalParameter.getDimension();
    }

    public double getParameterValue(int i) {
        return getParameterValue(i / diagonalParameter.getDimension(), i % diagonalParameter.getDimension());
    }

    public int getColumnDimension() {
        return diagonalParameter.getDimension();
    }

    public int getRowDimension() {
        return diagonalParameter.getDimension();
    }

}
