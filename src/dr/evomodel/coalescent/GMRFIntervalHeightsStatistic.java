/*
 * GMRFIntervalHeightsStatistic.java
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

package dr.evomodel.coalescent;

import dr.inference.model.Statistic;


/**
 * A statistic that reports the heights of all intercoalescent intervals for a GMRF skyride likelihood plot 
 *
 * @author Marc A. Suchard
 */
public class GMRFIntervalHeightsStatistic extends Statistic.Abstract {


    public GMRFIntervalHeightsStatistic(String name, GMRFSkyrideLikelihood skyrideLikelihood) {
        super(name);
        this.skyrideLikelihood = skyrideLikelihood;
    }

    public int getDimension() {
        return skyrideLikelihood.getCorrectFieldLength();
    }

    /**
     * @return the end-points in time of each coalescent interval
     */
    public double getStatisticValue(int dim) {
        if (dim == 0) {
            // This assumes that each dimension will be called in turn, so
            // the call for dim 0 updates the array.
            heights = skyrideLikelihood.getCoalescentIntervalHeights();
        }

        return heights[dim];
    }

    final private GMRFSkyrideLikelihood skyrideLikelihood;
    private double[] heights = null;
}