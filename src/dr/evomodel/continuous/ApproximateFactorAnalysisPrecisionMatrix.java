/*
 * ApproximateFactorAnalysisPrecisionMatrix.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;

/**
 * Created by msuchard on 11/8/16.
 */
public class ApproximateFactorAnalysisPrecisionMatrix extends CompoundParameter implements MatrixParameterInterface {

    private double[] values;
    private double[] storedValues;

    private boolean valuesKnown;
    private boolean storedValuesKnow;

    final private MatrixParameterInterface L;
    final private Parameter gamma;
    final int dim;

    public ApproximateFactorAnalysisPrecisionMatrix(String name,
                                                    MatrixParameterInterface L,
                                                    Parameter gamma) {
        super(name, new Parameter[] {L, gamma});
        this.L = L;
        this.gamma = gamma;
        this.dim = L.getColumnDimension();
        if (L.getRowDimension() != dim) {
            throw new IllegalArgumentException("Can only use square matrices");
        }
    }

    private void computeValues() {
        if (!valuesKnown) {
            computeValuesImp();
            valuesKnown = true;
        }
    }

    private void computeValuesImp() {
        double[][] matrix = new double[dim][dim];

        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                matrix[row][col] = L.getParameterValue(row, col) * L.getParameterValue(col, row);
            }

            matrix[row][row] += gamma.getParameterValue(row);
        }

        Matrix inverse = new Matrix(matrix).inverse();

        int index = 0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                values[index] = inverse.component(row, col);
                ++index;
            }
        }
    }

    @Override
    public double getParameterValue(int dim) {
        computeValues();
        throw new RuntimeException("Not yet implemented");
//        return values[dim];
    }

    @Override
    public double[][] getParameterAsMatrix() {
        computeValues();
        double[][] matrix = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(values, i * dim, matrix[i], 0, dim);
        }
        return matrix;
    }

    @Override
    public double getParameterValue(int row, int col) {
        computeValues();
        throw new RuntimeException("Not yet implemented");
//        return values[col * dim + row];
    }

    @Override
    public double[] getParameterValues() {
        computeValues();
        throw new RuntimeException("Not yet implemented");
//        double[] matrix = new double[values.length];
//        System.arraycopy(value, 0, matrix, 0, values.length);
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getColumnDimension() {
        return 0;
    }

    @Override
    public int getRowDimension() {
        return 0;
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
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void storeValues() {
        super.storeValues();

        System.arraycopy(values, 0, storedValues, 0, values.length);
        storedValuesKnow = valuesKnown;
    }

    @Override
    public void restoreValues() {
        super.restoreValues();

        double[] tmp = values;
        values = storedValues;
        storedValues = tmp;
        valuesKnown = storedValuesKnow;
    }
}
