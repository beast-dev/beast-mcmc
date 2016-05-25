/*
 * TransposedMatrixParameter.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class TransposedMatrixParameter extends MatrixParameter {

    public TransposedMatrixParameter(String name) {
        super(name + ".transpose");
    }

    public TransposedMatrixParameter(String name, Parameter[] parameters) {
        super(name + ".transpose", parameters);
    }

    public static TransposedMatrixParameter recast(String name, CompoundParameter compoundParameter) {
        final int count = compoundParameter.getParameterCount();
        Parameter[] parameters = new Parameter[count];
        for (int i = 0; i < count; ++i) {
            parameters[i] = compoundParameter.getParameter(i);
        }
        return new TransposedMatrixParameter(name, parameters);
    }

//    public double getParameterValue(int row, int col) {
//        // transposed
//        return super.getParameterValue(col, row);
//    }

    public int getColumnDimension() {
        // transposed
        return super.getParameter(0).getDimension();
    }

    public int getRowDimension() {
        // transposed
        return super.getParameterCount();
    }

    public int getParameterCount() {
        // MatrixParameter.getParamaterCount examines unique parameters
        // and not column dimension, as it probably should
        return getColumnDimension();
    }

    public double[][] getParameterAsMatrix() {
        final int I = getColumnDimension();
        final int J = getRowDimension();
        double[][] parameterAsMatrix = new double[J][I];
        for (int i = 0; i < I; i++) {
            for (int j = 0; j < J; j++)
                parameterAsMatrix[j][i] = getParameterValue(i, j);

        }
        return parameterAsMatrix;
    }

    public double getParameterValue(int dim) {
        // TODO Map to transposed dimension
        int transposedDim = dim;
        return super.getParameterValue(transposedDim);
    }

    public Parameter getParameter(int index) {
        if (slices == null) {
            // construct vector_slices
            slices = new ArrayList<Parameter>();
            for (int i = 0; i < getColumnDimension(); ++i) {
                VectorSliceParameter thisSlice = new VectorSliceParameter(getParameterName() + "." + i, i);
                for (int j = 0; j < getRowDimension(); ++j) {
                    thisSlice.addParameter(super.getParameter(j));
                }
                slices.add(thisSlice);
            }
        }
        return slices.get(index);
    }

//    @Override
//    public void setParameterValueQuietly(int row, int column, double a) {
//        super.setParameterValueQuietly(column,row, a);
//    }

    MatrixParameter transposeBack(){
        return MatrixParameter.recast(null, this);
    }

    private List<Parameter> slices = null;
}
