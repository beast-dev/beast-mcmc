/*
 * LogisticGrowthN0.java
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

/*
 * LogisticGrowthN0.java
 *
 * Daniel Wilson 4th October 2011
 *
 */

package dr.evomodel.epidemiology;

import dr.evolution.coalescent.*;

/**
 * This class models logistic growth.
 *
 * @author Daniel Wilson
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Matthew Hall
 * @version $Id: LogisticGrowth.java,v 1.15 2008/03/21 20:25:56 rambaut Exp $
 */
public class LogisticGrowthN0 extends ExponentialGrowth {

    /**
     * Construct demographic model with default settings
     */
    public LogisticGrowthN0(Type units) {
        super(units);
    }

	public void setT50(double value) {
		t50 = value;
	}
	
	public double getT50() {
		return t50;
	}

    // Implementation of abstract methods

    /**
     * Gets the value of the demographic function N(t) at time t.
     *
     * @param t the time
     * @return the value of the demographic function N(t) at time t.
     */
    public double getDemographic(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();
		
		return N0 * (1 + Math.exp(-r * T50)) / (1 + Math.exp(-r * (T50-t)));
    }

    public double getLogDemographic(double t) {
		return Math.log(getDemographic(t));
	}

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();
		double exp_rT50 = Math.exp(-r*T50);
		
		return (t + exp_rT50 * (Math.exp(r * t) - 1)/r) / (N0 * (1 + exp_rT50));
    }

    /**
     * Returns the inverse function of getIntensity
     *
     * If exp(-qt) = a(t-k) then t = k + (1/q) * W(q*exp(-q*k)/a) where W(x) is the Lambert W function.
     *
     * for our purposes:
     *
     * q = -r
     * a = (q/exp(q*T50))
     * k = N0*(1+exp(q*T50))*x - (1/q)exp(q*T50)
     *
     * For large x, W0(x) is approximately equal to ln(x) - ln(ln(x)); if q*exp(-q*k)/a rounds to infinity, we log it
     * and use this instead
     */
    public double getInverseIntensity(double x) {

        double q = -getGrowthRate();
        double T50 = getT50();
        double N0 = getN0();
        double a = (q/Math.exp(q*T50));
        double k = N0*(1+Math.exp(q*T50))*x - (1/q)*Math.exp(q*T50);

        double lambertInput = q*Math.exp(-q*k)/a;

        double lambertResult;

        if(lambertInput==Double.POSITIVE_INFINITY){

            //use the asymptote; note q/a = exp(q*T50)

            double logInput = q*T50-q*k;
            lambertResult = logInput - Math.log(logInput);

        } else {
            lambertResult = LambertW.branch0(lambertInput);
        }

        return k + (1/q)*lambertResult;

    }

    public double getIntegral(double start, double finish) {
		return getIntensity(finish) - getIntensity(start);
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
                return "t50";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getT50();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0:
                setN0(value);
                break;
            case 1:
                setGrowthRate(value);
                break;
            case 2:
                setT50(value);
                break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public DemographicFunction getCopy() {
        LogisticGrowthN0 df = new LogisticGrowthN0(getUnits());
        df.setN0(getN0());
        df.setGrowthRate(getGrowthRate());
        df.setT50(getT50());
        return df;
    }

    //
    // private stuff
    //

    private double t50;
}
