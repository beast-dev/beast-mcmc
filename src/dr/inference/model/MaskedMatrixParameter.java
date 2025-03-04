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
    //    private ArrayList<Integer> rows = new ArrayList<>();
    private int[] rows;
    private int[] cols;


    public MaskedMatrixParameter(MatrixParameterInterface matrix, Parameter mask) {
        super(matrix.getParameterName() + ".mask");
        this.matrix = matrix;
        this.mask = mask;
        addParameter(matrix);
        addParameter(mask);
        RowsAndCols rowsAndCols = makeRowsFromMask();
        this.rows = rowsAndCols.rows;
        this.cols = rowsAndCols.cols;

    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == mask) {
//            int[] oldRows = rows;
//            int[] oldCols = cols;
            RowsAndCols rowsAndCols = makeRowsFromMask();
            this.rows = rowsAndCols.rows;
            this.cols = rowsAndCols.cols;

            type = ChangeType.ALL_VALUES_CHANGED;

//            int ni = rows.length;
//            int oi = oldRows.length;
//            if (ni == oi) {
//                type = ChangeType.ALL_VALUES_CHANGED;
//            } else if (ni < oi) {
//                type = ChangeType.REMOVED;
//            } else {
//                type = ChangeType.ADDED;
//            }
            index = -1;
        }
        super.variableChangedEvent(variable, index, type);
    }

    private RowsAndCols makeRowsFromMask() {
        ArrayList<ArrayList<Integer>> allRows = new ArrayList<>();
        ArrayList<Integer> newCols = new ArrayList<>();
        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            ArrayList<Integer> colRows = new ArrayList<>();
            for (int j = 0; j < matrix.getRowDimension(); j++) {
                if (mask.getParameterValue(i * matrix.getRowDimension() + j) == 1) {
                    colRows.add(j);
                }
            }
            if (colRows.size() > 0) {
                newCols.add(i);
//                int[] rowsArray = new int[colrows.length];
//                for (int j = 0; j < colrows.length; i++) {
//                    rowsArray[j] = colRows.get(j);
//                }
                allRows.add(colRows);
            }

        }

        if (newCols.size() > 1) {
            for (int i = 1; i < newCols.size(); i++) {
                if (allRows.get(0).size() != allRows.get(i).size()) {
                    throw new RuntimeException("Invalid mask structure");
                }
                for (int j = 0; j < allRows.get(0).size(); j++) {
                    if (allRows.get(0).get(j) != allRows.get(i).get(j)) {
                        throw new RuntimeException("Invalid mask structure");
                    }
                }
            }
        }

        int[] rows = new int[allRows.get(0).size()];
        for (int i = 0; i < rows.length; i++) {
            rows[i] = allRows.get(0).get(i);
        }

        int[] cols = new int[newCols.size()];
        for (int i = 0; i < cols.length; i++) {
            cols[i] = newCols.get(i);
        }

        return new RowsAndCols(rows, cols);
    }

    private class RowsAndCols {

        public final int[] rows;
        public final int[] cols;

        RowsAndCols(int[] rows, int[] cols) {
            this.rows = rows;
            this.cols = cols;
        }
    }


    @Override
    public double getParameterValue(int row, int col) {
        return matrix.getParameterValue(rows[row], cols[col]);
    }


    @Override
    public void setParameterValue(int row, int col, double value) {
        matrix.setParameterValue(rows[row], cols[col], value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        matrix.setParameterValueQuietly(rows[row], cols[col], value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        matrix.setParameterValueNotifyChangedAll(rows[row], cols[col], value);
    }

    @Override
    public double[] getColumnValues(int col) {
        double[] maskedValues = new double[rows.length];
        for (int i = 0; i < rows.length; i++) {
            maskedValues[i] = matrix.getParameterValue(rows[i], cols[col]);
        }
        return maskedValues;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        double[][] values = new double[matrix.getColumnDimension()][rows.length];
        for (int col : cols) {
            for (int row : rows) {
                values[col][row] = getParameterValue(col, row);
            }
        }
        return values;
    }

    @Override
    public int getColumnDimension() {
        return cols.length;
    }

    @Override
    public int getRowDimension() {
        return rows.length;
    }

    @Override
    public int getDimension() {
        return getColumnDimension() * getRowDimension();
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
