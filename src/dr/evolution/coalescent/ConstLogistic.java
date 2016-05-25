/*
 * ConstLogistic.java
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
 * This class models logistic growth from an initial population size.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: ConstLogistic.java,v 1.4 2005/04/11 11:43:03 alexei Exp $
 *
 */
public class ConstLogistic extends LogisticGrowth {

	/**
	 * Construct demographic model with default settings
	 */
	public ConstLogistic(Type units) {
	
		super(units);
	}
	
	public double getN1() { return N1; }
	public void setN1(double N1) { this.N1 = N1; }
			
	// Implementation of abstract methods

	public double getDemographic(double t) {
		
		double nZero = getN0();
		double nOne = getN1();
		double r = getGrowthRate();
		double c = getShape();
		
		double common = Math.exp(-r*t);
		return nOne + ((nZero - nOne) * (1 + c) * common) / (c + common);
	}

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
	public double getIntensity(double t) {
        // Untested code
        double nZero = getN0();
		double nOne = getN1();
		double r = getGrowthRate();
		double c = getShape();

        double aa = nOne + (nZero-nOne)*(1+c);
        double bb = nOne*c;
        final double e = Math.exp(-r * t);
        double v1 = Math.log(bb + aa* e) / (aa*-r);
        double v2 = c * Math.log(aa + bb/e) / (bb*r);
        return v1 + v2;
        //throw new RuntimeException("Not implemented!");
	}

	public double getInverseIntensity(double x) {
		
		throw new RuntimeException("Not implemented!");
	}
	
	public double getIntegral(double start, double finish) {
		final double v1 = getIntensity(finish) - getIntensity(start);
		// Until the above getIntensity is implemented, numerically integrate
        final double numerical = getNumericalIntegral(start, finish);
        return numerical;
	}
		
	public int getNumArguments() {
		return 4;
	}
	
	public String getArgumentName(int n) {
		switch (n) {
			case 0: return "N0";
			case 1: return "r";
			case 2: return "c";
			case 3: return "N1";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getShape();
			case 3: return getN1();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setShape(value); break;
			case 3: setN1(value); break;
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
	
	private double N1 = 0.0;
}
