/*
 * PiecewiseLinearPopulation.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
 * @version $Id: PiecewiseLinearPopulation.java,v 1.7 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewiseLinearPopulation extends PiecewiseConstantPopulation {
		
	/**
	 * Construct demographic model with default settings
	 */
	public PiecewiseLinearPopulation(double[] intervals, double[] thetas, int units) {
	
		super(intervals, thetas, units);
	}
	
	// **************************************************************
	// Implementation of abstract methods
	// **************************************************************
	
	/**
	 * @return the value of the demographic function for the given epoch and time relative to start of epoch.
	 */
	protected final double getDemographic(int epoch, double t) {
		// if in last epoch then the population is flat.
		if (epoch == (thetas.length - 1)) {
			return getEpochDemographic(epoch);
		}
		
		double popSize1 = getEpochDemographic(epoch);
		double popSize2 = getEpochDemographic(epoch+1);
		
		double width = getEpochDuration(epoch);
		
		return popSize1 * (width-t) + (popSize2 * t) / width;
	}
	
	public DemographicFunction getCopy() {
		PiecewiseLinearPopulation df = new PiecewiseLinearPopulation(new double[intervals.length], new double[thetas.length], getUnits());
		System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
		System.arraycopy(thetas, 0, df.thetas, 0, thetas.length);
		
		return df;
	}
	
	/**
	 * @return the value of the intensity function for the given epoch.
	 */
	protected final double getIntensity(int epoch) {
		return 2.0 * getEpochDuration(epoch) / (getEpochDemographic(epoch) + getEpochDemographic(epoch+1));
	}
	
	/**
	 * @return the value of the intensity function for the given epoch and time relative to start of epoch.
	 */
	protected final double getIntensity(int epoch, double relativeTime) {
		return 2.0 * relativeTime / (getEpochDemographic(epoch) + getDemographic(epoch, relativeTime));
	}
}
