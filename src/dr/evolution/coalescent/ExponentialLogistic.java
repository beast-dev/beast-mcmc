/*
 * ExponentialLogistic.java
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
 * Exponential growth followed by Logistic growth.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: ConstLogistic.java,v 1.4 2005/04/11 11:43:03 alexei Exp $
 *
 */
public class ExponentialLogistic extends LogisticGrowth {

	/**
	 * Construct demographic model with default settings
	 */
	public ExponentialLogistic(Type units) {

		super(units);
	}

	public double getTime() {
		return time;
	}
	public void setTime(double time) {
		this.time = time;
	}

	public double getR1() {
		return r1;
	}

	public void setR1(double r1) {
		this.r1 = r1;
	}

	// Implementation of abstract methods

	public double getDemographic(double t) {

		double transition_time = getTime();

		// size of the population under the logistic at transition_time
		if (t < transition_time) {
		    return super.getDemographic(t);
		} else {
            double r1 = getR1();
			double N1 = super.getDemographic(transition_time);
			return N1 * Math.exp(-r1*(t - transition_time));
		}
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

	public double getIntegral(double start, double finish) {
		//final double v1 = getIntensity(finish) - getIntensity(start);
		// Until the above getIntensity is implemented, numerically integrate
        final double numerical = getNumericalIntegral(start, finish);
        return numerical;
	}

	public int getNumArguments() {
		return 5;
	}

	public String getArgumentName(int n) {
		switch (n) {
			case 0: return "N0";
			case 1: return "r0";
			case 2: return "c";
			case 3: return "r1";
			case 4: return "t1";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}

	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getShape();
			case 3: return getR1();
			case 4: return getTime();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}

	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setShape(value); break;
			case 3: setR1(value); break;
			case 4: setTime(value); break;
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

	/**
	 * The transition time between exponential and logistic
	 */
	private double time = 0.0;

	/**
	 * The growth rate of the exponential phase
	 */
	private double r1 = 0.0;

}