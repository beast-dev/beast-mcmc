/*
 * ExpConstant.java
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
 * This class models exponential growth up to a plateau
 *
 * For compatibility with the CaseToCase model, it is aware of the possibility that the change point is in negative
 * time
 * 
 * @author Matthew Hall
 */
public class ExpConstant extends ExponentialGrowth {

	/**
	 * Construct demographic model with default settings
	 */
	public ExpConstant(Type units) {
	
		super(units);
	}
	
	public double getTransitionTime() { return transitionTime; }
	public void setTransitionTime(double transitionTime) { this.transitionTime = transitionTime; }
	
	// Implementation of abstract methods

	public double getDemographic(double t) {
		
		double N0 = getN0();
		double r = getGrowthRate();
        double changeTime = getTransitionTime();
		
		//return nOne + ((nZero - nOne) * Math.exp(-r*t));
	
        if(t < changeTime){
            return N0*Math.exp(-r*changeTime);
        }
		
		return N0*Math.exp(-r*t);
	}

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
	public double getIntensity(double t) {

        double N0 = getN0();
        double r = getGrowthRate();
        double changeTime = getTransitionTime();

        double plateauLevel = N0*Math.exp(-r*changeTime);

        if (r == 0.0) return t/getN0();

        if(t==0){
            return 0;
        }
        if(changeTime <= 0){
            if(t >= changeTime){
                return super.getIntensity(t);
            } else {
                return (1/(N0*r)) * (Math.exp(r*changeTime) - 1) + (t-changeTime)/plateauLevel;
            }
        } else {
            if(t <= changeTime){
                return (t/plateauLevel);
            } else {
                return (changeTime/plateauLevel) + (1/(N0*r)) * (Math.exp(r*t) - Math.exp(r*changeTime));
            }
        }
  	}

	public double getInverseIntensity(double x) {
        double N0 = getN0();
        double r = getGrowthRate();
        double changeTime = getTransitionTime();

        double plateauLevel = N0*Math.exp(-r*changeTime);

        if (r == 0.0) {
            return getN0()*x;
        } else

        if(changeTime <= 0){
            if(x > (1/(N0*r)) * (Math.exp(r*changeTime) - 1)){
                return super.getInverseIntensity(x);
            } else {
                return Math.exp(-r*changeTime)*(N0*x - (Math.exp(r*changeTime)-1)/r) + changeTime;
            }
        } else {
            if(x < changeTime/plateauLevel){
                return plateauLevel*x;
            } else {
                return Math.log(r*(N0*x - changeTime*Math.exp(r*changeTime))+Math.exp(r*changeTime))/r;
            }
        }
	}
	
	public int getNumArguments() {
		return 3;
	}
	
	public String getArgumentName(int n) {
		switch (n) {
			case 0: return "N0";
			case 1: return "r";
			case 2: return "transitionTime";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getTransitionTime();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setTransitionTime(value); break;
			default: throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}

	public double getLowerBound(int n) {
        switch (n) {
            case 0: return 0;
            case 1: return 0;
            case 2: return Double.NEGATIVE_INFINITY;
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
}
