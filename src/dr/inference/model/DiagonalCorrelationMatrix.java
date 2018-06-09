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

import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.util.CorrelationToCholesky;

/**
 * @author Marc Suchard
 * @author Paul Bastide
 */
public class DiagonalCorrelationMatrix extends CompoundParameter implements MatrixParameterInterface {
    
    // TODO This is going to be a cleaner, more generic version of CompoundSymmetricMatrix

    enum SpaceTransform {
        NONE,
        INVERSE
    }

    enum CorrelationParametrization {
        AS_CORRELATION,
        AS_CHOLESKY
    }

    private final Parameter diagonalParameter;
    private final Parameter offDiagonalParameter;

    private final boolean asCorrelation;
    private final boolean isCholesky;
    private final int dim;

    public DiagonalCorrelationMatrix(String name,
                                     Parameter diagonals, Parameter offDiagonal,
                                     boolean asCorrelation, boolean isCholesky) {
        super(name, new Parameter[] { diagonals, offDiagonal});

        assert asCorrelation || !isCholesky; // cholesky only allowed when used as correlation.
        diagonalParameter = diagonals;
        dim = diagonalParameter.getDimension();
        if (!isCholesky) {
            offDiagonalParameter = offDiagonal;
        } else {
            offDiagonalParameter
                    = new TransformedMultivariateParameter(offDiagonal, new CorrelationToCholesky(dim), true);
        }
        addParameter(diagonalParameter);
        addParameter(offDiagonal);
        this.asCorrelation = asCorrelation;
        this.isCholesky = isCholesky;
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
    public String toString() {
        return toStringCompoundParameter(dim * (dim + 1) / 2);
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
            return offDiagonalParameter.getParameterValue(getUpperTriangularIndex(row, col));
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
    public int getUniqueParameterCount() {
        return 0;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return null;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {

    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {

    }

    @Override
    public String toSymmetricString() {
        return null;
    }

    @Override
    public void setParameterValue(int index, double a) {
        throw new RuntimeException("Do not set entries of a DiagonalCorrelationMatrix directly");
    }

    @Override
    public void setParameterValue(int row, int column, double a) {
        throw new RuntimeException("Do not set entries of a DiagonalCorrelationMatrix directly");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {

    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {

    }

    @Override
    public double[] getColumnValues(int col) {
        return new double[0];
    }

    public boolean isCholesky() {
        return isCholesky;
    }

    public boolean asCorrelation() {
        return asCorrelation;
    }

    public double[] getDiagonal() {
        return diagonalParameter.getParameterValues();
    }

    public Parameter getDiagonalParameter() { return diagonalParameter; }

    public Parameter getOffDiagonalParameter() { return offDiagonalParameter; }

    public double[][] getCorrelationMatrix() {
        SymmetricMatrix correlation
                = SymmetricMatrix.compoundCorrelationSymmetricMatrix(offDiagonalParameter.getParameterValues(), dim);
        if (!asCorrelation) {
            for (int i = 0; i < dim; i++) {
                for (int j = i + 1; j < dim; j++) {
                    correlation.setSymmetric(i, j,
                            correlation.component(i, j) / Math.sqrt(diagonalParameter.getParameterValue(i) * diagonalParameter.getParameterValue(j)));
                }
            }
        }
        return correlation.toComponents();
    }

    public double[] updateGradientCorrelation(double[] gradient) {
        if (!isCholesky) {
            return gradient;
        } else {
            CorrelationToCholesky transform = new CorrelationToCholesky(dim);
            double[][] jacobian
                    = transform.computeJacobianMatrixInverse(((TransformedMultivariateParameter) offDiagonalParameter).getParameterUntransformedValues());
            // Matrix multiplication (upper triangular)
            double[] updatedGradient = new double[gradient.length];
            for (int i = 0; i < gradient.length; i++) {
                for (int j = i; j < gradient.length; j++) {
                    updatedGradient[i] += jacobian[i][j] * gradient[j];
                }
            }
            return updatedGradient;
        }
    }
}
