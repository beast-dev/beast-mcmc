/*
 * MaskedMatrixParameter.java
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

import java.util.ArrayList;

public class MaskedMatrixParameter extends CompoundParameter implements MatrixParameterInterface, VariableListener {

    private final MatrixParameterInterface matrix;
    private final Parameter mask;
    private ArrayList<Integer> rows = new ArrayList<>();


    public MaskedMatrixParameter(MatrixParameterInterface matrix, Parameter mask) {
        super(matrix.getParameterName() + ".mask");
        this.matrix = matrix;
        this.mask = mask;
        addParameter(matrix);
        addParameter(mask);
        this.rows = makeRowsFromMask();
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == mask) {
            ArrayList<Integer> oldRows = rows;
            this.rows = makeRowsFromMask();
            int ni = rows.size();
            int oi = oldRows.size();
            if (ni == oi) {
                type = ChangeType.ALL_VALUES_CHANGED;
            } else if (ni < oi) {
                type = ChangeType.REMOVED;
            } else {
                type = ChangeType.ADDED;
            }
            index = -1;
        }
        super.variableChangedEvent(variable, index, type);
    }

    private ArrayList<Integer> makeRowsFromMask() {
        ArrayList<Integer> newRows = new ArrayList<>();
        for (int i = 0; i < mask.getDimension(); i++) {
            if (mask.getParameterValue(i) == 1) {
                newRows.add(i);
            }
        }
        return newRows;
    }


    @Override
    public double getParameterValue(int row, int col) {
        return matrix.getParameterValue(rows.get(row), col);
    }


    @Override
    public void setParameterValue(int row, int col, double value) {
        matrix.setParameterValue(rows.get(row), col, value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        matrix.setParameterValueQuietly(rows.get(row), col, value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        matrix.setParameterValueNotifyChangedAll(rows.get(row), col, value);
    }

    @Override
    public double[] getColumnValues(int col) {
        double[] maskedValues = new double[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            maskedValues[i] = matrix.getParameterValue(rows.get(i), col);
        }
        return maskedValues;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        double[][] values = new double[matrix.getColumnDimension()][rows.size()];
        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            for (int j = 0; j < rows.size(); j++) {
                values[i][j] = getParameterValue(i, j);
            }
        }
        return values;
    }

    @Override
    public int getColumnDimension() {
        return matrix.getColumnDimension();
    }

    @Override
    public int getRowDimension() {
        return rows.size();
    }

    @Override
    public int getUniqueParameterCount() {
        return matrix.getUniqueParameterCount();
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return matrix.getUniqueParameter(index);
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
