/*
 * LogisticGrowthN0Model.java
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
 * LogisticGrowthN0Model.java
 *
 * Daniel Wilson 4th October 2011
 *
 */

package dr.evomodel.epidemiology;

import dr.evomodel.coalescent.*;
import dr.evolution.coalescent.DemographicFunction;
import dr.inference.model.Parameter;

/**
 * Logistic growth.
 * 
 * @version $Id: LogisticGrowthN0Model.java,v 1.21 2008/03/21 20:25:57 rambaut Exp $
 *
 * @author Daniel Wilson
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LogisticGrowthN0Model extends DemographicModel
{
	//
	// Public stuff
	//
	
	/**
	 * Construct demographic model with default settings
	 */
	public LogisticGrowthN0Model(Parameter N0Parameter, Parameter growthRateParameter, 
								Parameter t50Parameter, Type units, boolean usingGrowthRate) {
		this(LogisticGrowthN0ModelParser.LOGISTIC_GROWTH_MODEL, N0Parameter, growthRateParameter, t50Parameter, units, usingGrowthRate);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public LogisticGrowthN0Model(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter t50Parameter, Type units, boolean usingGrowthRate) {	
		super(name);
		
		logisticGrowthN0 = new LogisticGrowthN0(units);
		
		this.N0Parameter = N0Parameter;
		addVariable(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addVariable(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

		this.t50Parameter = t50Parameter;
		addVariable(t50Parameter);
		t50Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

		this.usingGrowthRate = usingGrowthRate;

		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		logisticGrowthN0.setN0(N0Parameter.getParameterValue(0));
		if (usingGrowthRate) {
			double r = growthRateParameter.getParameterValue(0);
			logisticGrowthN0.setGrowthRate(r);
		} else {
			double doublingTime = growthRateParameter.getParameterValue(0);
			logisticGrowthN0.setDoublingTime(doublingTime);
		}
		logisticGrowthN0.setT50(t50Parameter.getParameterValue(0));
		
		return logisticGrowthN0;
	}
	
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter growthRateParameter = null;	
	Parameter t50Parameter = null;	
	LogisticGrowthN0 logisticGrowthN0 = null;
	boolean usingGrowthRate = true;
}
