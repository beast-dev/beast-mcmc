/*
 * TransformedCompoundMatrix.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.Transform;

/**
 * @author Marc Suchard
 * @author Paul Bastide
 */
abstract public class AbstractTransformedCompoundMatrix extends MatrixParameter {

    final Parameter diagonalParameter;
    final Parameter offDiagonalParameter;

    protected boolean isStrictlyUpperTriangular = true;

    final CompoundParameter untransformedCompoundParameter;

    protected final int dim;

    AbstractTransformedCompoundMatrix(Parameter diagonals, Parameter offDiagonal) {
        super(MATRIX_PARAMETER);
        diagonalParameter = diagonals;
        dim = diagonalParameter.getDimension();
        offDiagonalParameter = offDiagonal;
        addParameter(diagonalParameter);
        addParameter(offDiagonalParameter);

        untransformedCompoundParameter = new CompoundParameter(this.getParameterName());
        setUntransformedCompoundParameter();
    }

    AbstractTransformedCompoundMatrix(Parameter diagonals, Parameter offDiagonal,
                                      Transform.MultivariableTransform transform, Boolean inverse) {
        super(MATRIX_PARAMETER);
        diagonalParameter = diagonals;
        dim = diagonalParameter.getDimension();
        offDiagonalParameter =
                (transform == null) ? offDiagonal : new TransformedMultivariateParameter(offDiagonal, transform, inverse);
        addParameter(diagonalParameter);
        addParameter(offDiagonalParameter);

        untransformedCompoundParameter = new CompoundParameter(this.getParameterName());
        setUntransformedCompoundParameter();
    }

    private void setUntransformedCompoundParameter() {
        untransformedCompoundParameter.addParameter(getDiagonalParameter());
        untransformedCompoundParameter.addParameter(getUntransformedOffDiagonalParameter());
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
    public double getParameterUntransformedValue(int index) {
        final int dim = getColumnDimension();
        return getParameterUntransformedValue(index / dim, index % dim);
    }

    private double getParameterUntransformedValue(int row, int col) {
        if (row != col) {
            return offDiagonalParameter.getParameterUntransformedValue(getUpperTriangularIndex(row, col));
        }
        return diagonalParameter.getParameterUntransformedValue(row);
    }

    @Override
    public void setParameterUntransformedValue(int index, double a) {
        final int dim = getColumnDimension();
        setParameterUntransformedValue(index / dim, index % dim, a);
    }

    private void setParameterUntransformedValue(int row, int col, double a) {
        if (row != col) {
            offDiagonalParameter.setParameterUntransformedValue(getUpperTriangularIndex(row, col), a);
        } else {
            diagonalParameter.setParameterUntransformedValue(row, a);
        }
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
        throw new RuntimeException("Do not set entries of a TransformedCompoundMatrix directly");
    }

    @Override
    public void setParameterValue(int row, int column, double a) {
        throw new RuntimeException("Do not set entries of a TransformedCompoundMatrix directly");
    }

    public double[] getDiagonal() {
        return diagonalParameter.getParameterValues();
    }

    public Parameter getDiagonalParameter() {
        return diagonalParameter;
    }

    public Parameter getOffDiagonalParameter() {
        return offDiagonalParameter;
    }

    public Parameter getUntransformedOffDiagonalParameter() {
        if (offDiagonalParameter instanceof TransformedMultivariateParameter) {
            return ((TransformedMultivariateParameter) offDiagonalParameter).getUntransformedParameter();
        }
        return offDiagonalParameter;
    }

    public CompoundParameter getUntransformedCompoundParameter() {
        return untransformedCompoundParameter;
    }

    abstract double[] updateGradientDiagonal(double[] gradient);

    abstract double[] updateGradientOffDiagonal(double[] gradient);

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
    public String getDimensionName(int index) {
        int dim = getColumnDimension();
        String row = Integer.toString(index / dim);
        String col = Integer.toString(index % dim);

        return getId() + row + col;
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }

    int getUpperTriangularIndex(int i, int j) {
        if (i < j) {
            return upperTriangularTransformation(i, j);
        } else {
            return upperTriangularTransformation(j, i);
        }
    }

    private int upperTriangularTransformation(int i, int j) {
        if (isStrictlyUpperTriangular) {
            assert i != j;
            return strictlyUpperTriangularTransformation(i, j);
        }
        return weaklyUpperTriangularTransformatino(i, j);
    }

    private int strictlyUpperTriangularTransformation(int i, int j) {
        return i * (2 * dim - i - 1) / 2 + (j - i - 1);
    }

    private int weaklyUpperTriangularTransformatino(int i, int j) {
        return i * (2 * dim - i + 1) / 2 + (j - i);
    }

    public void setStrictlyUpperTriangular(boolean b) {
        this.isStrictlyUpperTriangular = b;
    }
}
