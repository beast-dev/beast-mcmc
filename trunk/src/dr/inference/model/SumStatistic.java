/*
 * SumStatistic.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: SumStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 */
public class SumStatistic extends Statistic.Abstract {

    private int dimension = 0;
    private final boolean elementwise;
    private final boolean absolute;

    public SumStatistic(String name, boolean elementwise, boolean absolute) {
        super(name);
        this.elementwise = elementwise;
        this.absolute = absolute;
    }

    public SumStatistic(String name, boolean elementwise) {
        this(name, elementwise, false);
    }

    public void addStatistic(Statistic statistic) {
        if (!elementwise) {
            if (dimension == 0) {
                dimension = statistic.getDimension();
            } else if (dimension != statistic.getDimension()) {
                throw new IllegalArgumentException();
            }
        } else {
            dimension = 1;
        }
        statistics.add(statistic);
    }

    public int getDimension() {
        return elementwise ? 1 : dimension;
    }

    /**
     * @return mean of contained statistics
     */
    public double getStatisticValue(int dim) {

        double sum = 0.0;

        if (elementwise) {
            assert dim == 0;
        }

        for (Statistic statistic : statistics) {
            if (elementwise) {
                if (absolute) {
                    //  System.err.println("statistic.getDimension(): " + statistic.getDimension());
                    for (int j = 0; j < statistic.getDimension(); j++) {
                        sum += Math.abs(statistic.getStatisticValue(j));
                    }
                } else {
                    for (int j = 0; j < statistic.getDimension(); j++) {
                        sum += statistic.getStatisticValue(j);
                        //    if(statistic.getStatisticValue(j)<0) {
                        //      System.err.println("statisticValue: " + statistic.getStatisticValue(j));
                        //  }
                    }
                }
            } else {
                sum += statistic.getStatisticValue(dim);
            }
        }

        return sum;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private final List<Statistic> statistics = new ArrayList<Statistic>();
}
