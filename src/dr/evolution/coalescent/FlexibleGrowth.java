/*
 * FlexibleGrowth.java
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
 * before zero in backwards time), this models population growth according to the function N(t)=N0*K(-t)^r/(1+K(-t)^(r-1)
 *
 * This takes a variety of shapes depending on the value of r. For r<0 it climbs to a peak and then declines. At r=0 it
 * approaches an asymptote; r=1 is linear.
 *
 * (if times are transformed to the log scale, this is logistic growth)
 *
 * @author Matthew Hall
 */
public class FlexibleGrowth extends PowerLawGrowth
{
	//
	// Public stuff
	//

	/**
	 * Construct demographic model with default settings
     * @param units
     */
	public FlexibleGrowth(Type units) {

		super(units);
	}

	/**
	 * @return initial population size.
	 */

	public double getK(){
		return K;
	}


	public void setK(double K) {

		this.K = K;
	}

    // Implementation of abstract methods
	
	public double getDemographic(double t) {
		if(t>0){
			throw new RuntimeException("Negative times only! t="+t);
		}

		return getN0()*K*Math.pow(-t,getR())/(1+K*Math.pow(-t, getR()-1));
	}
	public double getIntensity(double t) {
		throw new RuntimeException("getIntensity is not implemented (and not finite); use getIntegral instead");
	}

	public double getInverseIntensity(double x) {
		throw new RuntimeException("Not implemented");
	}

	public double getInverseIntegral(double x, double start) {
		throw new RuntimeException("Not implemented");
	}

	public double getIntegral(double start, double finish) {

		return 1/getN0() * (1/((getR()-1)*K) * (Math.pow(-finish, -getR() + 1) - Math.pow(-start, -getR() +1 )) + Math.log((-start)/(-finish)));
	}

	public int getNumArguments() {
		return 3;
	}

	public String getArgumentName(int n) {
		switch (n) {
			case 0:
				return "N0";
			case 1:
				return "r";
			case 2:
				return "K";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}

	public double getArgument(int n) {
		switch (n) {
			case 0:
				return getN0();
			case 1:
				return getR();
			case 2:
				return getK();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}


	public void setArgument(int n, double value) {
		switch (n) {
			case 0:
				setN0(value);
				break;
			case 1:
				setR(value);
				break;
			case 2:
				setK(value);
				break;
			default:
				throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}

	public double getLowerBound(int n) {
		switch (n) {
			case 0:
				return 0;
			case 1:
				return Double.NEGATIVE_INFINITY;
			case 2:
				return 0;
			default:
				throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}
	
	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	public DemographicFunction getCopy() {
		FlexibleGrowth df = new FlexibleGrowth(getUnits());
		df.setN0(getN0());
		df.setR(getR());
		df.K = K;
		
		return df;
	}

	//
	// private stuff
	//

	private double K;
}
