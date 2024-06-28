/*
 * TransformedMatrixParameter.java
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

package dr.inference.model;

import dr.util.Transform;


public class TransformedMatrixParameter extends TransformedMultivariateParameter implements MatrixParameterInterface {

    private final int rowDim;
    private final int colDim;

    public TransformedMatrixParameter(Parameter parameter, Transform.MultivariableTransform transform,
                                      Boolean inverse, int rowDim, int colDim) {
        super(parameter, transform, inverse);
        this.rowDim = rowDim;
        this.colDim = colDim;
    }


    public TransformedMatrixParameter(MatrixParameterInterface parameter, Transform.MultivariableTransform transform,
                                      Boolean inverse) {
        this(parameter, transform, inverse, parameter.getRowDimension(), parameter.getColumnDimension());
    }


    @Override
    public double getParameterValue(int row, int col) {
        return getParameterValue(index(row, col));
    }

    @Override
    public Parameter getParameter(int column) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        setParameterValue(index(row, col), value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        setParameterValueQuietly(index(row, col), value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        setParameterValueNotifyChangedAll(index(row, col), value);
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {
        return MatrixParameterInterface.getParameterAsMatrix(this);
    }

    @Override
    public int getColumnDimension() {
        return colDim;
    }

    @Override
    public int getRowDimension() {
        return rowDim;
    }

    @Override
    public int getUniqueParameterCount() {
        return 1;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return parameter;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

}


