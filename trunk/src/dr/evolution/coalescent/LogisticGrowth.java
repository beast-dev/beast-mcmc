/*
 * LogisticGrowth.java
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
 * This class models logistic growth.
 * 
 * @version $Id: LogisticGrowth.java,v 1.15 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LogisticGrowth extends ExponentialGrowth {
	
	/**
	 * Construct demographic model with default settings
	 */
	public LogisticGrowth(Type units) {
	
		super(units);
	}
	
	public void setShape(double value) { c = value; }
	public double getShape() { return c; }
						
	/**
	 * An alternative parameterization of this model. This
	 * function sets the time at which there is a 0.5 proportion
	 * of N0.
	 */
	public void setTime50(double time50) { 
	
		c = 1.0 / (Math.exp(getGrowthRate() * time50) - 2.0);
		
        // The general form for any k where t50 is the time at which Nt = N0/k:
        //		c = (k - 1.0) / (Math.exp(getGrowthRate() * time50) - k);
	}
						
	// Implementation of abstract methods

    /**
     * Gets the value of the demographic function N(t) at time t.
     * @param t the time
     * @return the value of the demographic function N(t) at time t.
     */
	public double getDemographic(double t) {

		double nZero = getN0();
		double r = getGrowthRate();
		double c = getShape();

//		return nZero * (1 + c) / (1 + (c * Math.exp(r*t)));
//		AER rearranging this to use exp(-rt) may help
// 		with some overflow situations...

		double common = Math.exp(-r*t);
		return (nZero * (1 + c) * common) / (c + common);
	}

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
	public double getIntensity(double t) {
		throw new RuntimeException("Not implemented!");
	}

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
	public double getInverseIntensity(double x) {
		
		throw new RuntimeException("Not implemented!");
	}
	
	public double getIntegral(double start, double finish) {
		
		double intervalLength = finish - start;
		
		double nZero = getN0();
		double r = getGrowthRate();
		double c = getShape();
		double expOfMinusRT = Math.exp(-r*start);
		double expOfMinusRG = Math.exp(-r*intervalLength);

		double term1 = nZero*(1.0+c);
		if (term1==0.0) {
			throw new RuntimeException("Infinite integral!");
		}

		double term2 = c*(1.0 - expOfMinusRG);

		double term3 = (term1*expOfMinusRT) * r * expOfMinusRG;
		if (term3==0.0 && term2>0.0) {
			throw new RuntimeException("Infinite integral!");
		}

		double term4;
		if (term3!=0.0 && term2==0.0) {term4=0.0;}
		else if (term3==0.0 && term2==0.0) {
		    throw new RuntimeException("term3 and term2 are both zeros. N0=" + getN0() + " growthRate=" +  getGrowthRate() + "c=" + c);
		}    
		else {term4 = term2 / term3;}

		double term5 = intervalLength / term1;

		return term5 + term4;
	}
	
	public int getNumArguments() {
		return 3;
	}
	
	public String getArgumentName(int n) {
		switch (n) {
			case 0: return "N0";
			case 1: return "r";
			case 2: return "c";
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getShape();
		}
		throw new IllegalArgumentException("Argument " + n + " does not exist");
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setShape(value); break;
			default: throw new IllegalArgumentException("Argument " + n + " does not exist");

		}
	}

	public double getLowerBound(int n) {
		return 0.0;
	}
	
	public double getUpperBound(int n) {
		return Double.POSITIVE_INFINITY;
	}

	public DemographicFunction getCopy() {
		LogisticGrowth df = new LogisticGrowth(getUnits());
		df.setN0(getN0());
		df.setGrowthRate(getGrowthRate());
		df.c = c;
		
		return df;
	}

	//
	// private stuff
	//

	private double c;	
}
