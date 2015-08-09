/*
 * PowerLawGrowth.java
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
 * before zero in backwards time), this models population growth according to a power law N(x)=N0*(-x)^r, where r>1
 * (Sadly this does not work for r<1)
 *
 * @author Matthew Hall
 */
public class PowerLawGrowth extends LinearGrowth
{
	//
	// Public stuff
	//

	/**
	 * Construct demographic model with default settings
     * @param units
     */
	public PowerLawGrowth(Type units) {

		super(units);
	}

	/**
	 * @return initial population size.
	 */

	public double getR(){
		return r;
	}


	public void setR(double r) {


		this.r = r;
	}

		
	// Implementation of abstract methods
	
	public double getDemographic(double t) {
		if(t>0){
			throw new RuntimeException("Negative times only!");
		}

		return getN0()*Math.pow(-t,r);
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

		return (Math.pow(-finish, -r+1) - Math.pow(-start, -r+1))/(getN0()*(r-1));
	}


	public int getNumArguments() {
		return 2;
	}

	public String getArgumentName(int n) {
		if (n == 0) {
			return "N0";
		} else {
			return "r";
		}
	}

	public double getArgument(int n) {
		if (n == 0) {
			return getN0();
		} else {
			return getR();
		}
	}

	public void setArgument(int n, double value) {
		if (n == 0) {
			setN0(value);
		} else {
			setR(value);
		}
	}

	public double getLowerBound(int n) {

		if (n == 0) {
			return 0;
		} else {
			return 1;
		}
	}


	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	public DemographicFunction getCopy() {
		PowerLawGrowth df = new PowerLawGrowth(getUnits());
		df.setN0(getN0());
		df.r=r;
		
		return df;
	}

	//
	// private stuff
	//

	private double r;
}
