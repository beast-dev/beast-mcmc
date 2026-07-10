/*
 * InfinitesimalRatesLogger.java
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

package dr.evomodel.substmodel;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.util.Transform;
import dr.xml.Reportable;

public class InfinitesimalRatesLogger implements Loggable, Reportable {

    public InfinitesimalRatesLogger(SubstitutionModel substitutionModel, Transform transform,
                                    boolean diagonalElements, String order, Integer subset) {
        this.substitutionModel = substitutionModel;
        this.diagonalElements = diagonalElements;
        this.transform = transform;
        this.order = order;
        this.subset = subset;
        this.stateCount = substitutionModel.getDataType().getStateCount();
    }

    public InfinitesimalRatesLogger(SubstitutionModel substitutionModel,  Transform transform,
                                    boolean diagonalElements, String order) {
        this(substitutionModel, transform, diagonalElements, order, null);
    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = make();
        }
        return columns;
    }

    private LogColumn[] make() {
        int nOutputs;

        if (subset != null) {
            nOutputs = subset;
        } else {
            nOutputs = stateCount * stateCount;
            if (!diagonalElements) nOutputs -= stateCount;
        }
        LogColumn[] columns = new LogColumn[nOutputs];

        if (generator == null) {
            generator = new double[stateCount * stateCount];
        }

        int index = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                if (!diagonalElements && i == j) continue;
                final int row = i;
                final int col = j;
                int indexK = 0;
                if ("byrow".equals(order)) {
                    indexK = index++;
                } else if ("bycol".equals(order)) {
                    throw new IllegalArgumentException("Order " + order + " not implemented");
                } else if ("rowCol".equals(order)) {
                    indexK = rowColIndexCreator(i, j);
                } else {
                    throw new IllegalArgumentException("Invalid order: " + order);
                }
                if (subset != null && indexK >= subset) {continue;}
                final int k = indexK;
                columns[k] = new NumberColumn(substitutionModel.getId() + "." + (i + 1) + "." + (j + 1)) {
                    @Override
                    public double getDoubleValue() {
                        if (k == 0) { // Refresh at first-element read
                            substitutionModel.getInfinitesimalMatrix(generator);
                        }
                        if (transform != null) {
                            return transform.transform(generator[row * stateCount + col]);
                        } else {
                            return generator[row * stateCount + col];
                        }
                    }
                };
            }

        }

        return columns;
    }

    private int rowColIndexCreator(int i, int j) {
        int nUpperTri = stateCount * (stateCount + 1 ) / 2;
        if (!diagonalElements) nUpperTri -= stateCount + 1;
        int k = 0;
        if (i <= j) {
            k = (stateCount - 1) * i - i * (i - 1) / 2 + j - i - 1;
            if (diagonalElements) k += i + 1;
        } else {
            k = nUpperTri;
            k += (stateCount - 1) * j - j * (j - 1) / 2 + i - j;
        }
        return k;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        for (LogColumn column : getColumns()) {
            sb.append(column.getFormatted() + " ");
        }
        sb.append("\n");
        return sb.toString();
    }

    private LogColumn[] columns;
    private final int stateCount;
    private final Transform transform;
    private final SubstitutionModel substitutionModel;
    private final boolean diagonalElements;
    private final String order;
    private double[] generator;
    private final Integer subset;
}
