/*
 * ConstExponential.java
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
 * This class models exponential growth from an initial population size which
 * then transitions back to a constant population size.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $ID$
 *
 */
public class ConstExpConst extends ConstExponential {

	/**
	 * Construct demographic model with default settings
	 */
	public ConstExpConst(Type units) {
		super(units);
	}

	public double getTime1() {
		return time1;
	}

	public void setTime1(double time1) {
		this.time1 = time1;
	}

	public double getTime2() {
		return time1 + (-Math.log(getN1()/getN0())/getGrowthRate());
	}

	public void setProportion(double p) {
		this.setN1(getN0() * p);
	}

	// Implementation of abstract methods

	public double getDemographic(double t) {

		double N0 = getN0();
		double N1 = getN1();
		double r = getGrowthRate();
		double t1 = getTime1();
				
		if (t < t1) {
			return N0;
		}

		double t2 = getTime2();
		if (t >= t2) {
			return N1;
		}

		return N0 * Math.exp(-r*(t - t1));
	}

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
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
			case 1: return "r";
			case 2: return "N1";
			case 3: return "time1";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}

	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getN1();
			case 3: return getTime1();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}

	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setN1(value); break;
			case 3: setTime1(value); break;
			default: throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}

	public double getLowerBound(int n) {
		return 0.0;
	}

	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	//
	// private stuff
	//

	private double time1 = 0.0;
}