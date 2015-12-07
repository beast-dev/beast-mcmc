/*
 * CataclysmicDemographic.java
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
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 * 
 * @version $Id: CataclysmicDemographic.java,v 1.5 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class CataclysmicDemographic extends ExponentialGrowth {
	
	//
	// Public stuff
	//

	/**
	 * Construct demographic model with default settings
	 */
	public CataclysmicDemographic(Type units) {
	
		super(units);
	}

	/**
	 * returns the positive-valued decline rate
	 */
	public final double getDeclineRate() { return -d; }
	
	/**
	 * sets the decline rate.
	 */
	public void setDeclineRate(double d) { 
//		if (d <= 0) throw new IllegalArgumentException();
		this.d = d; 
	}
	
	public final double getCataclysmTime() { return catTime; }
	
	public final void setCataclysmTime(double t) { 
		if (t <= 0) throw new IllegalArgumentException();
		catTime = t; 
	}
			
	/**
	 * An alternative parameterization of this model. This
	 * function sets the decline rate using N0 & t which must
	 * already have been set.
	 */
	public final void setSpikeFactor(double f) {
		setDeclineRate( Math.log(f) / catTime );
	}
			
	// Implementation of abstract methods

	public double getDemographic(double t) {

		double d = getDeclineRate();

		if (t < catTime) {
			return getN0() * Math.exp(t * d);
		} else {
			double spikeHeight = getN0() * Math.exp(catTime * d);
			//System.out.println("Spike height = " + spikeHeight);
			t -= catTime;
		
			double r = getGrowthRate();
			if (r == 0) {
				return spikeHeight;
			} else {
				return spikeHeight * Math.exp(-t * r);
			}
		}
	}

	public double getIntensity(double t) {
	
		double d = getDeclineRate();
		double r = getGrowthRate();
		if (t < catTime) {
			return (Math.exp(t*-d)-1.0)/getN0()/-d;
		} else {
			
			double intensityUpToSpike = (Math.exp(catTime*-d)-1.0)/getN0()/-d;
			
			double spikeHeight = getN0() * Math.exp(catTime * d);
			t -= catTime;
			//System.out.println("Spike height = " + spikeHeight);
			
			if (r == 0) {
				return t/spikeHeight + intensityUpToSpike;
			} else {
				return (Math.exp(t*r)-1.0)/spikeHeight/r + intensityUpToSpike;
			}
		}
	}
	
	public double getInverseIntensity(double x) {
		double d = getDeclineRate();
		double r = getGrowthRate();
		double intensityUpToSpike = (Math.exp(catTime*-d)-1.0)/getN0()/-d;
		
		if(x < intensityUpToSpike){
			return -Math.log(1.0 - getN0() * d * x)/d;
		}
		
		
		double spikeHeight = getN0() * Math.exp(catTime * d);
		x -= intensityUpToSpike;
		
		if(r == 0){
			return spikeHeight*x + catTime;
		}
		
		return catTime + Math.log(1.0 + spikeHeight * x * r)/r;
		//throw new UnsupportedOperationException();
	}
	
	public int getNumArguments() {
		return 4;
	}
	
	public String getArgumentName(int n) {
		
		switch (n) {
			case 0: return "N0";
			case 1: return "r";
			case 2: return "d";
			case 3: return "t";
			default: throw new IllegalArgumentException();
		}
		
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getDeclineRate();
			case 3: return getCataclysmTime();
			default: throw new IllegalArgumentException();
		}
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setDeclineRate(value); break;
			case 3: setCataclysmTime(value); break;
			default: throw new IllegalArgumentException();
		}
	}

	public DemographicFunction getCopy() {
		CataclysmicDemographic df = new CataclysmicDemographic(getUnits());
		df.setN0(getN0());
		df.setGrowthRate(getGrowthRate());
		df.d = d;
		df.catTime = catTime;
		
		return df;
	}

	//
	// private stuff
	//

	private double d;
	private double catTime;	
}
