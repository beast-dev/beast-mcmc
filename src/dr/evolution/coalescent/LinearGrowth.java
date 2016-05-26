/*
 * LinearGrowth.java
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
 * Under the assumption that the effective population size at time zero was zero (and all coalescent events are
 * before zero in backwards time), this models linear population groth
 *
 * @author Matthew Hall
 */
public class LinearGrowth extends ConstantPopulation
{
	//
	// Public stuff
	//

	/**
	 * Construct demographic model with default settings
     * @param units
     */
	public LinearGrowth(Type units) {

		super(units);
	}

	/**
	 * @return initial population size.
	 */
	public double getN0() { return N0; }

	/**
	 * sets initial population size.
     * @param N0 new size
     */
	public void setN0(double N0) { this.N0 = N0; }

		
	// Implementation of abstract methods
	
	public double getDemographic(double t) {
		if(t>0){
			throw new RuntimeException("Negative times only!");
		}

		return -getN0()*t;
	}
	public double getIntensity(double t) {
		throw new RuntimeException("getIntensity is not implemented (and not finite); use getIntegral instead");
	}

	public double getInverseIntensity(double x) {
		throw new RuntimeException("Not implemented");
	}

	public double getIntegral(double start, double finish) {
		return 1/getN0() * Math.log((-start)/(-finish));
	}

	public double getInverseIntegral(double x, double start){

		return start*Math.exp(-(x*getN0()));

	}

	public int getNumArguments() {
		return 1;
	}
	
	public String getArgumentName(int n) {
		return "N0";
	}
	
	public double getArgument(int n) {
		return getN0();
	}
	
	public void setArgument(int n, double value) {
		setN0(value);
	}

	public double getLowerBound(int n) {
		return 0.0;
	}
	
	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	public DemographicFunction getCopy() {
		LinearGrowth df = new LinearGrowth(getUnits());
		df.N0 = N0;
		
		return df;
	}

	//
	// private stuff
	//
	
	private double N0;
}
