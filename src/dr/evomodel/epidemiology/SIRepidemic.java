/*
 * SIRepidemic.java
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

package dr.evomodel.epidemiology;

import dr.evolution.coalescent.*;

/**
 * This class models an SIR epidemic.
 * Uses the new parameterization (23rd February 2012)
 *
 * @author Daniel Wilson
 */
public class SIRepidemic extends ODEDemographicFunction {

	public SIRepidemic(Type units) {
		super(units);
		// Number of variables in the Runge-Kutta scheme
		nvar = 2;
	}

	public void setN0(double value) { N0 = value; RKinit(); }
	public double getN0() { return N0; }

	public void setGrowthRate(double value) { r = value; RKinit(); }
	public double getGrowthRate() { return r; }
	public void setDoublingTime(double value) { setGrowthRate(Math.log(2)/value); }
	public double getDoublingTime() { return Math.log(2)/r; }
	
	public void setu0(double value) { u0 = value; RKinit(); }
	public double getu0() { return u0; }
	
	public void setR0(double value) { R0 = value; RKinit(); }
	public double getR0() { return R0; }
	
	public void setIntegrateIntensity() { integrate_intensity = true; }
	public void unsetIntegrateIntensity() { integrate_intensity = false; }
	
	public double getbeta() { return getGrowthRate()/(1.0-1.0/getR0()); }
	public double getgamma() { return getGrowthRate()/(getR0()-1.0); }

	public double gets0() { return s0; }
	public double geti0() { return i0; }
	
	/* Calculate minimum proportion infectious from R0 */
	public double getsmin(double _R0) {
		if(_R0<=1.0) return 0.0;
		return -W.branch0(-_R0*Math.exp(-_R0))/_R0;
	}
	
	/* Calculate proportion infectious from proportion susceptible */
	public double s_to_i(double susceptible) {
		double i = 1.0-susceptible+Math.log(susceptible)/R0;
		// Do not allow i to reach boundaries at 0 and 1.
		if(i<1e-12) i = 1e-12;
		if(i>1.0-1e-12) i = 1.0-1e-12;
		return i;
	}
	
	/* Convert proportion susceptible to transformed proportion susceptible, where smin is the minimum proportion of susceptibles */
	public double s_to_u(double susceptible, double smin) {
		return Math.log(susceptible-smin)-Math.log(1.0-susceptible);
	}

	/* Convert transformed proportion susceptible to proportion susceptible, where smin is the minimum proportion of susceptibles  */
	public double u_to_s(double u, double smin) {
		double s = smin+(1.0-smin)/(1.0+Math.exp(-u));
		// Do not allow s to reach boundaries at smin and 1
		if(s<smin+1e-12) s = smin+1e-12;
		if(s>1.0-1e-12) s = 1.0-1e-12;
		return s;
	}
	
	/* Convert transformed proportion susceptible to proportion infectious, where smin is the minimum proportion of susceptibles  */
	public double u_to_i(double u, double smin, double _R0) {
		double s = smin+(1.0-smin)/(1.0+Math.exp(-u));
		double i = 1.0-s+Math.log(s)/_R0;
		// Do not allow i to reach boundaries at 0 and 1.
		if(i<1e-12) i = 1e-12;
		if(i>1.0-1e-12) i = 1.0-1e-12;
		return i;
	}
	
	/**
	 * Variable 0 is the integrated intensity function (by definition)
	 * Variable 1 is the transformed proportion of susceptibles
	 */
	void derivs(double t, double[] y, double[] dydt) {
		double Lambda = y[0];
		if(!integrate_intensity) {
			// Switch off to avoid underflow/overflow problems when wishing to calculate prevalence only
			dydt[0] = 1.0;
		} else if(Lambda<LAMBDA_MAX) {
			double Ne = getDemographicFromPrevalence(y,t);
			dydt[0] = 1.0/Ne;
		} else {
			// When the coalescence intensity is already extremely large, set a ceiling. This will have the effect of inducing
			// a zero likelihood when the density of coalescence times is evaluated. The interpretation is that all coalescence
			// events must have already occurred by this point.
			dydt[0] = 0.0;
		}
		double s = u_to_s(y[1],smin);
		if(s==1.0) {
			dydt[1] = beta-gamma;
		} else {
			if(s<smin+1e-6) s = smin+1e-6;
			dydt[1] = (beta * s * (1.0-s) + gamma * s * Math.log(s)) * (1.0-smin) / (s-smin) / (1.0-s);
		}
		if(Double.isNaN(dydt[0]) | Double.isInfinite(dydt[0]) | Double.isNaN(dydt[1]) | Double.isInfinite(dydt[1])) {
			System.out.println("t = " + t + ", y = {" + y[0] + ", " + y[1] + "}, s = " + s + ", dydt = {" + dydt[0] + ", " + dydt[1] + "}");
			derivs(t,y,dydt);
		}
	}

	/**
	 * Calculate the effective population size at time t
	 */
	double getDemographicFromPrevalence(double[] y, double t) {
		double s = u_to_s(y[1],smin);
		double i = u_to_i(y[1],smin,getR0());
		double Ne = N0 * s0 / i0 * i / s;
		// Do not allow NaN or very small values
		if(Double.isNaN(Ne) | Ne<1e-12) Ne = 1e-12;
		return Ne;
	}

	/**
	 * 
	 */
	public double getDemographic(double t) {
    	if(!valid) return 0.0;
    	double ret = super.getDemographic(t);
		if(Ynow[0]>=LAMBDA_MAX) return 0.0;
		return ret;
	}
	
	/**
	 * 
	 */
	public double getIntensity(double t) {
    	if(!valid) return Math.log(0.0);
    	double ret = super.getIntensity(t);
		if(ret>=LAMBDA_MAX) return Math.log(0.0);
		return ret;
	}

	/**
	 * Obtain the proportion of susceptibles at time t (compare to super.getDemographic and super.getIntensity)
	 */
	public double getSusceptibles(double t) {
    	if(!valid) return Math.log(0.0);
		Evaluate(t);
		if(RKfail) return Math.log(0.0);
		if(Ynow[0]>=LAMBDA_MAX) return Math.log(0.0);
    	return u_to_s(Ynow[1],smin);
	}
	
	/**
	 * Obtain the transformed proportion of susceptibles at time t (compare to super.getDemographic and super.getIntensity)
	 */
	public double getTransformedSusceptibles(double t) {
    	if(!valid) return Math.log(0.0);
		Evaluate(t);
		if(RKfail) return Math.log(0.0);
    	return Ynow[1];
	}
	
	/**
	 * Overload RKinit to make use of "valid"
	 */
	void RKinit() {
		super.RKinit();
		/* Set internal variables */
		valid = true;
		smin = getsmin(getR0());
		s0 = u_to_s(u0,smin);
		i0 = u_to_i(u0,smin,getR0());
		beta = getbeta();
		gamma = getgamma();
	}

	/**
	 * Initially the integrated intensity function (Y[0]) is zero
	 * and the transformed proportion of susceptibles (Y[1]) is u0
	 */
	void setInit() {
		Y[0][0] = 0.0;
		Y[1][0] = u0;
	}

	public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getu0();
            case 3:
            	return getR0();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N0";
            case 1:
                return "r";
            case 2:
                return "u0";
            case 3:
            	return "R0";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public DemographicFunction getCopy() {
        SIRepidemic df = new SIRepidemic(getUnits());
        df.setN0(getN0());
        df.setGrowthRate(getGrowthRate());
        df.setu0(getu0());
        df.setR0(getR0());
        // Copy RK variables?
        return df;
    }

    /**
     * Returns the inverse function of getIntensity
     */
    public double getInverseIntensity(double x) {
		// Can be implemented using simple linear interpolation
        throw new RuntimeException("Not implemented!");
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public int getNumArguments() {
        return 4;
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
                setu0(value);
                break;
            case 3:
            	setR0(value);
            	break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

	// Parameters
    private boolean valid = false;
	private double N0 = 0;		// Initial effective population size (0<N0)
	private double r  = 0;		// Intrinsic growth rate (0<r)
	private double u0 = 0;		// Initial transformed progress of the epidemic (-Inf<u0<Inf)
	private double R0 = 0;		// Basic reproductive number (1<R0)
	private boolean integrate_intensity = true;	// Flag indicating whether the integrated intensity function is calculated
	// Pre-calculated from parameters
	private double smin = 0;	// Minimum proportion of susceptibles (depends on R0)
	private double s0 = 0;		// Initial progress of the epidemic (smin<s0<1)
	private double i0 = 0;		// Initial proportion infectious
	private double beta = 0;	// Initial proportion infectious
	private double gamma = 0;	// Initial proportion infectious
	LambertW W;					// Needed for calculating smin
	// Maximum possible value of Lambda
	private double LAMBDA_MAX = 1000.0;
}
