/*
 * CoalescentTreeIntervalStatistic.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Guy Baele
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;


public class CoalescentTreeIntervalStatistic extends Statistic.Abstract {

    private Tree tree;
    private TreeIntervals intervals;

    public CoalescentTreeIntervalStatistic(Tree tree) {
        this.tree = tree;
        this.intervals = new TreeIntervals(tree);
    }

    public int getDimension() {
    	return intervals.getIntervalCount();
    }

    public double getStatisticValue(int i) {
    	//need more elegant way to trigger recalculation of intervals
    	if (i == 0) {
    		intervals.setIntervalsUnknown();
    	}
    	
        double interval = intervals.getInterval(i);

        return interval;
    }
}
