/*
 * LogisticGrowthModel.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.LogisticGrowth;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Logistic growth.
 * 
 * @version $Id: LogisticGrowthModel.java,v 1.21 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LogisticGrowthModel extends DemographicModel {
	
	//
	// Public stuff
	//
	
	public static String POPULATION_SIZE = "populationSize";
	public static String LOGISTIC_GROWTH_MODEL = "logisticGrowth";

	public static String GROWTH_RATE = "growthRate";
	public static String DOUBLING_TIME = "doublingTime";
	public static String TIME_50 = "t50";
	public static String ALPHA = "alpha";
	

	/**
	 * Construct demographic model with default settings
	 */
	public LogisticGrowthModel(Parameter N0Parameter, Parameter growthRateParameter, 
								Parameter shapeParameter, double alpha, int units, 
								boolean usingGrowthRate) {
	
		this(LOGISTIC_GROWTH_MODEL, N0Parameter, growthRateParameter, shapeParameter, alpha, units, usingGrowthRate);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public LogisticGrowthModel(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, int units, boolean usingGrowthRate) {
	
		super(name);
		
		logisticGrowth = new LogisticGrowth(units);
		
		this.N0Parameter = N0Parameter;
		addParameter(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addParameter(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.shapeParameter = shapeParameter;
		addParameter(shapeParameter);
		shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.alpha = alpha;
		this.usingGrowthRate = usingGrowthRate;
		
		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
	
		logisticGrowth.setN0(N0Parameter.getParameterValue(0));
		
		if (usingGrowthRate) {
			double r = growthRateParameter.getParameterValue(0);
			logisticGrowth.setGrowthRate(r);
		} else {
			double doublingTime = growthRateParameter.getParameterValue(0);
			logisticGrowth.setDoublingTime(doublingTime);
		}
		
		logisticGrowth.setTime50(shapeParameter.getParameterValue(0));
		
		return logisticGrowth;
	}
	
	/**
	 * Parses an element from an XMLObject into LogisticGrowthModel. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return LOGISTIC_GROWTH_MODEL; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			int units = XMLParser.Utils.getUnitsAttr(xo);
				
			XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
			Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);		
					
			boolean usingGrowthRate = true;		
			Parameter rParam;		
			if (xo.getChild(GROWTH_RATE) != null) {			
				cxo = (XMLObject)xo.getChild(GROWTH_RATE);
				rParam = (Parameter)cxo.getChild(Parameter.class);
			} else {
				cxo = (XMLObject)xo.getChild(DOUBLING_TIME);
				rParam = (Parameter)cxo.getChild(Parameter.class);
				usingGrowthRate = false;
			}
			
			cxo = (XMLObject)xo.getChild(TIME_50);
			Parameter cParam = (Parameter)cxo.getChild(Parameter.class);
				
			return new LogisticGrowthModel(N0Param, rParam, cParam, 0.5, units, usingGrowthRate);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "Logistic growth demographic model.";
		}

		public Class getReturnType() { return LogisticGrowthModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			XMLUnits.SYNTAX_RULES[0],
			new ElementRule(POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
				"This parameter represents the carrying capacity (maximum population size). " + 
				"If the shape is very large then the current day population size will be very close to the carrying capacity."),
			new XORRule(
		
				new ElementRule(GROWTH_RATE, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"This parameter determines the rate of growth during the exponential phase. See exponentialGrowth for details."),
				new ElementRule(DOUBLING_TIME, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"This parameter determines the doubling time at peak growth rate.")
			),
			new ElementRule(TIME_50, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
				"This parameter represents the time in the past when the population had half of the carrying capacity (population size). " +
				"It is therefore a positive number with the same units as divergence times. " +
				"A scale operator is recommended with a starting value near zero. " + 
				"A lower bound of zero should be employed and an upper bound is required!" )
		};
	};
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter growthRateParameter = null;	
	Parameter shapeParameter = null;	
	double alpha = 0.5;
	LogisticGrowth logisticGrowth = null;
	boolean usingGrowthRate = true;
}
