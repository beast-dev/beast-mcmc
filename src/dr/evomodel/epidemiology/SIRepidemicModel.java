/*
 * SIRepidemicModel.java
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
 * SIRepidemicModel.java
 *
 */

package dr.evomodel.epidemiology;

import dr.evomodel.coalescent.*;
import dr.evolution.coalescent.DemographicFunction;
import dr.inference.model.Parameter;

/**
 * SIR epidemic.
 * 
 * @author Daniel Wilson
 */
public class SIRepidemicModel extends DemographicModel
{
	//
	// Public stuff
	//
	
	/**
	 * Construct demographic model with default settings
	 */
	public SIRepidemicModel(Parameter N0Parameter, Parameter growthRateParameter, 
								Parameter tpeakParameter, Parameter gammaParameter, Type units, boolean usingGrowthRate, double minPrevalence) {
		this(SIRepidemicModelParser.SIREPI_MODEL, N0Parameter, growthRateParameter, tpeakParameter, gammaParameter, units, usingGrowthRate, minPrevalence);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public SIRepidemicModel(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter tpeakParameter, Parameter gammaParameter, Type units, boolean usingGrowthRate, double minPrevalence) {	
		super(name);
				
		sirepi = new SIRepidemic(units);
		
		this.N0Parameter = N0Parameter;
		addVariable(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addVariable(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.tpeakParameter = tpeakParameter;
		addVariable(tpeakParameter);
		tpeakParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

		this.gammaParameter = gammaParameter;
		addVariable(gammaParameter);
		gammaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		
		this.usingGrowthRate = usingGrowthRate;

		this.minPrevalence = minPrevalence;
		
		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		double r;
		if (usingGrowthRate) {
			r = growthRateParameter.getParameterValue(0);
		} else {
			double doublingTime = growthRateParameter.getParameterValue(0);
			r = Math.log(2)/doublingTime;
		}
		double tpeak = tpeakParameter.getParameterValue(0);
		double gamma = gammaParameter.getParameterValue(0);
		double R0 = 1.0 + (r/gamma);
		sirepi.setR0(R0);
		double smin = sirepi.getsmin(R0);
		double u0 = sirepi.s_to_u(1.0/R0,smin);
		// Calculate s0
		sirepi.setN0(1.0);
		sirepi.setu0(u0);
		sirepi.unsetIntegrateIntensity();
		if(tpeak>=0) {
			// this has the effect of running the process in reverse, so s0 is evaluated at -tpeak. (Negative times are not allowed in the ODE solver).
			sirepi.setGrowthRate(-r);
			u0 = sirepi.getTransformedSusceptibles(tpeak);
		} else {
			sirepi.setGrowthRate(r);
			u0 = sirepi.getTransformedSusceptibles(-tpeak);
		}
		// Did the RK fail when trying to calculate u0?
		boolean RKfailed = sirepi.RKfail;
		// If i0 is less than minPrevalence, force fail
		double i0 = sirepi.u_to_i(u0,smin,R0);
		if(i0<minPrevalence) RKfailed = true;
		// Set the remaining values to their intended values
		sirepi.setIntegrateIntensity();
		sirepi.setN0(N0Parameter.getParameterValue(0));
		sirepi.setGrowthRate(r);
		sirepi.setu0(u0);
		// If the RK failed when trying to calculate u0, impose RK failure now
		sirepi.RKfail = RKfailed;
		return sirepi;
	}
	
	//
	// protected stuff
	//

	Parameter N0Parameter = null;
	Parameter growthRateParameter = null;
	Parameter tpeakParameter = null;
	Parameter gammaParameter = null;
	SIRepidemic sirepi = null;
	boolean usingGrowthRate = true;
	double minPrevalence = 0.0;
}
