/*
 * DiagonalFloorMatrixView.java
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

/**
 * Read-only view of a square matrix parameter with a constant added to every diagonal entry:
 * {@code view[i][j] = matrix[i][j] + (i == j ? floor : 0)}.
 *
 * <p>The view forwards change events of the wrapped matrix and never modifies it.  Because the
 * added term is constant, the derivative of the view with respect to the wrapped matrix is the
 * identity, so gradients pulled back through the wrapped matrix's own parameterization remain
 * valid.  Used by {@code OUProcessModel} to apply a declared diffusion floor
 * {@code Q + floor * I} consistently to every likelihood and gradient engine.</p>
 */
public class DiagonalFloorMatrixView extends CompoundParameter implements MatrixParameterInterface {

    private final MatrixParameterInterface matrix;
    private final double floor;

    public DiagonalFloorMatrixView(final String id,
                                   final MatrixParameterInterface matrix,
                                   final double floor) {
        super(id);
        if (matrix == null) {
            throw new IllegalArgumentException("matrix must not be null");
        }
        if (matrix.getRowDimension() != matrix.getColumnDimension()) {
            throw new IllegalArgumentException("DiagonalFloorMatrixView requires a square matrix");
        }
        if (!Double.isFinite(floor) || floor < 0.0) {
            throw new IllegalArgumentException("floor must be finite and non-negative, found " + floor);
        }
        this.matrix = matrix;
        this.floor = floor;
        matrix.addVariableListener(this);
    }

    public MatrixParameterInterface getWrappedMatrix() {
        return matrix;
    }

    public double getFloor() {
        return floor;
    }

    @Override
    public double getParameterValue(final int row, final int col) {
        final double entry = matrix.getParameterValue(row, col);
        return row == col ? entry + floor : entry;
    }

    @Override
    public Parameter getParameter(final int column) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setParameterValue(final int row, final int col, final double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(final int row, final int col, final double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(final int row, final int col, final double value) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public double[] getColumnValues(final int col) {
        final double[] out = new double[getRowDimension()];
        for (int row = 0; row < out.length; ++row) {
            out[row] = getParameterValue(row, col);
        }
        return out;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        return MatrixParameterInterface.getParameterAsMatrix(this);
    }

    @Override
    public int getDimension() {
        return matrix.getDimension();
    }

    @Override
    public double getParameterValue(final int dim) {
        // Flat storage of the wrapped matrix; for a square matrix the diagonal entries sit at
        // multiples of (rows + 1) in both row-major and column-major layouts.
        final double entry = matrix.getParameterValue(dim);
        return dim % (getRowDimension() + 1) == 0 ? entry + floor : entry;
    }

    @Override
    public double[] getParameterValues() {
        final double[] out = new double[getDimension()];
        for (int i = 0; i < out.length; ++i) {
            out[i] = getParameterValue(i);
        }
        return out;
    }

    @Override
    public String getDimensionName(final int i) {
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
        return matrix.getUniqueParameterCount();
    }

    @Override
    public Parameter getUniqueParameter(final int index) {
        return matrix.getUniqueParameter(index);
    }

    @Override
    public void copyParameterValues(final double[] destination, final int offset) {
        final double[] values = getParameterValues();
        System.arraycopy(values, 0, destination, offset, values.length);
    }

    @Override
    public void setAllParameterValuesQuietly(final double[] values, final int offset) {
        throw new IllegalArgumentException("Not implemented");
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return matrix.isConstrainedSymmetric();
    }

    @Override
    public void variableChangedEvent(final Variable variable, final int index, final Parameter.ChangeType type) {
        if (variable == matrix) {
            fireParameterChangedEvent(index, type);
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }
}
