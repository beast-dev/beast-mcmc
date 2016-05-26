/*
 * EmpiricalPiecewiseConstant.java
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

package dr.evolution.coalescent;

/**
 *
 * @version $Id: EmpiricalPiecewiseConstant.java,v 1.4 2004/10/01 23:30:16 alexei Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class EmpiricalPiecewiseConstant extends DemographicFunction.Abstract {

	/**
	 * Creates a piecewise constant model with the given break points.
	 * @param intervals an array of intervals, Each interval represents time
	 * over which a different population size is assumed, the first interval starts from the
	 * present and proceeds back into the past.
	 * @param popSizes the effective population sizes of each interval. The last population size represents the
	 * time from the end of the last interval out to infinity. 
	 * @param lag a lag to align the break points with the most recent sample time (must be positive)
	 */
	public EmpiricalPiecewiseConstant(double[] intervals, double[] popSizes, double lag, Type units) {
		super(units);
		if (popSizes == null || intervals == null) { throw new IllegalArgumentException(); }
		if (popSizes.length != intervals.length + 1) { throw new IllegalArgumentException(); }
		if (lag < 0.0) throw new IllegalArgumentException("Lag must be greater than 1.");

		this.intervals = intervals;
		this.popSizes = popSizes;
		this.lag = lag;
	}
	
	public void setLag(double lag) {
		this.lag = lag;
	}
	
	public void setPopulationSizes(double[] popSizes) {
		this.popSizes = popSizes;
	}
	
	// **************************************************************
	// Implementation of abstract methods
	// **************************************************************
	
	public double getDemographic(double t) {
        int epoch = 0;
		double t1 = t+lag;
		while (t1 > getEpochDuration(epoch)) {
			t1 -= getEpochDuration(epoch);
			epoch += 1;
		}
		return getDemographic(epoch, t1);
	}

	public double getIntensity(double t) { 
	
		// find the first epoch that is involved
		double t2 = lag;
		int epoch = 0;
		while (t2 > getEpochDuration(epoch)) {
			t2 -= getEpochDuration(epoch);		
			epoch += 1;
		}
		
		// add last fraction of first epoch
		double intensity = getIntensity(epoch)-getIntensity(epoch, t2);
				
		double t1 = t-(getEpochDuration(epoch)-t2);
		epoch += 1;
		while (t1 > getEpochDuration(epoch)) {
			t1 -= getEpochDuration(epoch);		
			intensity += getIntensity(epoch);
			epoch += 1;
		}
		// add last fraction of intensity
        // when t1 may be negative (for example when t is in the first epoch) the intensity need
        // to be substracted
        intensity += t1 >= 0 ? getIntensity(epoch, t1) : getIntensity(epoch-1, t1);
	
		return intensity; 
	}
	
	public double getInverseIntensity(double x) { 
		throw new RuntimeException("Not implemented!");	
	}
	
	public double getUpperBound(int i) { return 1e9;}
	public double getLowerBound(int i) { return Double.MIN_VALUE;}
	
	public int getNumArguments() { return 1; }
	
	public String getArgumentName(int i) { 
		return "lag";
	}
	
	public double getArgument(int i) { 
		return lag;
	}
	
	public void setArgument(int i, double value) { 
		lag = value; 
	}
	
	public DemographicFunction getCopy() {
		EmpiricalPiecewiseConstant df = new EmpiricalPiecewiseConstant(new double[intervals.length], new double[popSizes.length], lag, getUnits());
		System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
		System.arraycopy(popSizes, 0, df.popSizes, 0, popSizes.length);
		
		return df;
	}

	
	/**
	 * @return the value of the demographic function for the given epoch and time relative to start of epoch.
	 */
	protected double getDemographic(int epoch, double t) {
		return getEpochDemographic(epoch);
	}
	
	/**
	 * @return the value of the intensity function for the given epoch.
	 */
	protected double getIntensity(int epoch) {
		return getEpochDuration(epoch) / getEpochDemographic(epoch);
	}
	
	/**
	 * @return the value of the intensity function for the given epoch and time relative to start of epoch.
	 */
	protected double getIntensity(int epoch, double relativeTime) {
		return relativeTime / getEpochDemographic(epoch);
	}
	
	/** 
	 * @return the duration of the specified epoch (in whatever units this demographic model is specified in).
	 */
    public double getEpochDuration(int epoch) {
        if (epoch < intervals.length) {
            return intervals[epoch];
        }
        return Double.POSITIVE_INFINITY;
    }

	/**
	 * @return the pop size of a given epoch.
	 */
	public double getEpochDemographic(int epoch) {
		if (epoch >= popSizes.length) { throw new IllegalArgumentException(); }
		return popSizes[epoch];
	}
	
	public String toString() {
		StringBuffer buffer = new StringBuffer();
		buffer.append(popSizes[0]);
		for (int i =1; i < popSizes.length; i++) {
            buffer.append("\t").append(popSizes[i]);
		}
		return buffer.toString();
	}
	
	double[] intervals;
	double[] popSizes;
	double lag;
}
