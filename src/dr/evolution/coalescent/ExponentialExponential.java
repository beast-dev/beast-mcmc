/*
 * ExponentialExponential.java
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
 * This class models a two phase exponential growth
 *
 * @author Andrew Rambaut
 */
public class ExponentialExponential extends ExponentialGrowth {

	/**
	 * Construct demographic model with default settings
	 */
	public ExponentialExponential(Type units) {
	
		super(units);
	}
	
	public double getTransitionTime() { return transitionTime; }
	public void setTransitionTime(double transitionTime) { this.transitionTime = transitionTime; }

    /**
     * @return growth rate.
     */
    public final double getAncestralGrowthRate() { return r1; }

    /**
     * sets growth rate to r1.
     * @param r1
     */
    public void setAncestralGrowthRate(double r1) { this.r1 = r1; }


    // Implementation of abstract methods

	public double getDemographic(double t) {
		
		double N0 = getN0();
		double r = getGrowthRate();
        double r1 = getAncestralGrowthRate();
        double changeTime = getTransitionTime();
		
		//return nOne + ((nZero - nOne) * Math.exp(-r*t));
	
        if(t < changeTime) {
            return N0*Math.exp(-r * t);
        }

        double N1 = N0 * Math.exp(-r * changeTime);

		return N1 * Math.exp(-r1* (t - changeTime));
	}

     /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    public double getIntegral(double start, double finish)
    {
        return getNumericalIntegral(start, finish);
    }

	public double getIntensity(double t) {

        throw new RuntimeException("Not implemented!");
  	}

	public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
	}
	
	public int getNumArguments() {
		return 4;
	}
	
	public String getArgumentName(int n) {
		switch (n) {
			case 0: return "N0";
			case 1: return "r0";
            case 2: return "r1";
			case 3: return "transitionTime";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
            case 2: return getAncestralGrowthRate();
			case 3: return getTransitionTime();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
            case 2: setAncestralGrowthRate(value); break;
			case 3: setTransitionTime(value); break;
			default: throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}

	public double getLowerBound(int n) {
        switch (n) {
            case 0: return 0;
            case 1: return Double.NEGATIVE_INFINITY;
            case 2: return Double.NEGATIVE_INFINITY;
            case 3: return Double.NEGATIVE_INFINITY;
            default: throw new IllegalArgumentException("Argument " + n + " does not exist");
        }
	}
	
	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	//
	// private stuff
	//

    private double transitionTime;
    private double r1;
}
