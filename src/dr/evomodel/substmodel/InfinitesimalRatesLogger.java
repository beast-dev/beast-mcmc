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

public class InfinitesimalRatesLogger implements Loggable {

    public InfinitesimalRatesLogger(SubstitutionModel substitutionModel) {
        this.substitutionModel = substitutionModel;
    }

    @Override
    public LogColumn[] getColumns() {
        int stateCount = substitutionModel.getDataType().getStateCount();

        if (generator == null) {
            generator = new double[stateCount * stateCount];
        }

        LogColumn[] columns = new LogColumn[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                final int k = i * stateCount + j;
                columns[k] = new NumberColumn(substitutionModel.getId() + "." + (i + 1) + "." + (j + 1)) {
                    @Override
                    public double getDoubleValue() {
                        if (k == 0) { // Refresh at first-element read
                            substitutionModel.getInfinitesimalMatrix(generator);
                        }
                        return generator[k];
                    }
                };
            }
        }

        return columns;
    }

    private final SubstitutionModel substitutionModel;
    private double[] generator;
}
