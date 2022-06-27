/*
 * PiecewiseLinearTimeDependentModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class PiecewiseLinearTimeDependentModelLogger implements Loggable {

    private final PiecewiseLinearTimeDependentModel model;
    private LogColumn[] columns;

    public PiecewiseLinearTimeDependentModelLogger(PiecewiseLinearTimeDependentModel model) {
        this.model = model;
    }



    @Override
    public LogColumn[] getColumns() {

        if (columns == null) {
            List<LogColumn> tmp = new ArrayList<>();

            final int dim = 3;

            for (int i = 0; i < dim; ++i) {
                final int index = i;
                tmp.add(new NumberColumn("intercept" + (i + 1)) {
                    @Override
                    public double getDoubleValue() {
                        PiecewiseLinearTimeDependentModel.SlopeInterceptPack pack = model.getSlopeInterceptPack();
                        return pack.intercepts[index];
                    }
                });
            }

            for (int i = 0; i < dim; ++i) {
                final int index = i;
                tmp.add(new NumberColumn("slope" + (i + 1)) {
                    @Override
                    public double getDoubleValue() {
                        PiecewiseLinearTimeDependentModel.SlopeInterceptPack pack = model.getSlopeInterceptPack();
                        return pack.slopes[index];
                    }
                });
            }

            for (int i = 0; i < dim; ++i) {
                final int index = i;
                tmp.add(new NumberColumn("break" + (i + 1)) {
                    @Override
                    public double getDoubleValue() {
                        PiecewiseLinearTimeDependentModel.SlopeInterceptPack pack = model.getSlopeInterceptPack();
                        return pack.breaks[index];
                    }
                });
            }

            columns = tmp.toArray(new LogColumn[] { });
        }

        return columns;
    }
}
