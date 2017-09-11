/*
 * ProductStatistic.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 */
public abstract class AbstractAlgebraStatistic extends Statistic.Abstract {

    private int firstDataDimension = 0;

    private int dimension = 0;
    private boolean elementwise;
    private final double[] constants;

    public AbstractAlgebraStatistic(String name, boolean elementwise, double[] constants) {
        super(name);
        this.elementwise = elementwise;
        this.constants = constants;
        if (constants != null) {
            firstDataDimension = constants.length;
        }
    }

    public void addStatistic(Statistic statistic) {
        if (dimension > 0) {
            throw new IllegalArgumentException("Cannot add more data after the statistic has been evaluated.");
        }

        if (!elementwise) {
            if (firstDataDimension == 0) {
                firstDataDimension = statistic.getDimension();
            } else  if (firstDataDimension != statistic.getDimension()) {
                throw new IllegalArgumentException("Data added doesn't match previous dimensions");
            }
        }

        statistics.add(statistic);
    }

    public int getDimension() {
        if (dimension == 0) {
            if ((constants != null || statistics.size() > 1) && !elementwise) {
                dimension = statistics.get(0).getDimension();
            } else {
                elementwise = true;
                dimension = 1;
            }
        }

        return dimension;
    }

    /**
     * @return product of contained statistics
     */
    public double getStatisticValue(int dim) {

        // Check we have set the dimension
        getDimension();

        double value = Double.NaN;

        for (Statistic statistic : statistics) {
            if (elementwise) {
                // do the operation across all the elements

                if (Double.isNaN(value)) {
                    value = statistic.getStatisticValue(0);
                }
                for (int j = 1; j < statistic.getDimension(); j++) {
                    value = doOperation(value, statistic.getStatisticValue(j));
                }
            } else {
                if (Double.isNaN(value)) {
                    value = statistic.getStatisticValue(dim);
                } else {
                    value = doOperation(value, statistic.getStatisticValue(dim));
                }
            }
        }
        if (constants != null && constants.length > 1) {
            if (constants.length == 1) {
                value = doOperation(value, constants[0]);
            } else {
                value = doOperation(value, constants[dim]);
            }
        }

        return value;
    }

    protected abstract double doOperation(double a, double b);

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private final List<Statistic> statistics = new ArrayList<Statistic>();
}
