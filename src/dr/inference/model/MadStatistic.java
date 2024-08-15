/*
 * MadStatistic.java
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

import dr.stats.DiscreteStatistics;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andy Magee
 */
public class MadStatistic extends Statistic.Abstract {

    public MadStatistic(String name, boolean centerIsMedian, boolean dispersionIsMedian) {
        super(name);
        this.centerIsMedian = centerIsMedian;
        this.dispersionIsMedian = dispersionIsMedian;
    }

    public void addStatistic(Statistic statistic) {
        statistics.add(statistic);
    }

    public int getDimension() { return 1; }

    /** @return mean of contained statistics */
    public double getStatisticValue(int dim) {
        List<Double> xList = new ArrayList<>();
        int dimensionCount = 0;
        int n;

        for (Statistic statistic : statistics) {
            n = statistic.getDimension();

            for (int j = 0; j < n; j++) {
                xList.add(statistic.getStatisticValue(j));
            }
            dimensionCount += n;
        }

        double[] x = new double[dimensionCount];
        for (int i = 0; i < dimensionCount; i++) {
            x[i] = xList.get(i);
        }

        double center;
        if ( centerIsMedian ) {
            center = DiscreteStatistics.median(x);
        } else {
            center = DiscreteStatistics.mean(x);
        }

        double[] absDiff = new double[dimensionCount];
        for (int i = 0; i < dimensionCount; i++) {
            absDiff[i] = Math.abs(x[i] - center);
        }

        double mad;
        if ( dispersionIsMedian ) {
            mad = DiscreteStatistics.median(absDiff);
        } else {
            mad = DiscreteStatistics.mean(absDiff);
        }
        return mad;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private List<Statistic> statistics = new ArrayList<Statistic>();
    private boolean centerIsMedian;
    private boolean dispersionIsMedian;

}
