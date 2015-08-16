/*
 * CoalescentTreeIntervalStatistic.java
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

/**
 * @author Guy Baele
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;

/**
* @author Guy Baele
*/
public class CoalescentTreeIntervalStatistic extends Statistic.Abstract {

    private Tree tree;
    private TreeIntervals intervals;
    private int dimension;
    private double[] nonZeroIntervals;
    private int[] lineageCounts;

    public CoalescentTreeIntervalStatistic(Tree tree) {
        this.tree = tree;
        this.intervals = new TreeIntervals(tree);
        //detect here if there are zero-length intervals and set the dimension accordingly
        setDimension(intervals);
        this.lineageCounts = new int[this.dimension];
        this.nonZeroIntervals = new double[this.dimension];
        removeZeroLengthIntervals(intervals);
    }

    public int getDimension() {
    	return this.dimension;
    }

    public double getStatisticValue(int i) {
    	//need more elegant way to trigger recalculation of intervals
    	if (i == 0) {
    		intervals.setIntervalsUnknown();
    		removeZeroLengthIntervals(intervals);
    	}
    	
        //double interval = intervals.getInterval(i);
    	//return interval;
    	
    	return this.nonZeroIntervals[i];
    }
    
    public double getLineageCount(int i) {
    	//need more elegant way to trigger recalculation of intervals
    	if (i == 0) {
    		intervals.setIntervalsUnknown();
    		removeZeroLengthIntervals(intervals);
    	}
    	
    	return this.lineageCounts[i];
    }
    
    public String getStatisticName() {
    	return "coalescentTreeIntervalStatistic";
    }
    
    private void setDimension(TreeIntervals intervals) {
    	int dim = intervals.getIntervalCount();
    	for (int i = 0; i < intervals.getIntervalCount(); i++) {
    		if (intervals.getInterval(i) == 0.0) {
    			dim--;
    		}
    	}
    	this.dimension = dim;
    }
    
    private void removeZeroLengthIntervals(TreeIntervals intervals) {
    	int dim = 0;
    	for (int i = 0; i < intervals.getIntervalCount(); i++) {
    		double intervalLength = intervals.getInterval(i);
    		if (intervalLength != 0.0) {
    			this.nonZeroIntervals[dim] = intervalLength;
    			this.lineageCounts[dim] = intervals.getLineageCount(i);
    			dim++;
    		}
    	}
    }
}
