/*
 * ExpansionModel.java
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
import dr.evolution.coalescent.Expansion;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Exponential growth from a constant ancestral population size.
 * 
 * @version $Id: ExpansionModel.java,v 1.5 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class ExpansionModel extends DemographicModel
{
	
	//
	// Public stuff
	//
	
	public static final String EXPANSION_MODEL = "expansion";
	public static final String POPULATION_SIZE = "populationSize";
	public static final String ANCESTRAL_POPULATION_PROPORTION = "ancestralPopulationProportion";

	public static final String GROWTH_RATE = "growthRate";
	public static final String DOUBLING_TIME = "doublingTime";

	/**
	 * Construct demographic model with default settings
	 */
	public ExpansionModel(Parameter N0Parameter, Parameter N1Parameter, 
									Parameter growthRateParameter, Type units, boolean usingGrowthRate) {
	
		this(EXPANSION_MODEL, N0Parameter, N1Parameter, growthRateParameter, units, usingGrowthRate);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public ExpansionModel(String name, Parameter N0Parameter, Parameter N1Parameter, 
									Parameter growthRateParameter, Type units, boolean usingGrowthRate) {
	
		super(name);
		
		expansion = new Expansion(units);
		
		this.N0Parameter = N0Parameter;
		addParameter(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.N1Parameter = N1Parameter;
		addParameter(N1Parameter);
		N1Parameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addParameter(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.usingGrowthRate = usingGrowthRate;

		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		
		double N0 = N0Parameter.getParameterValue(0);
		double N1 = N1Parameter.getParameterValue(0);
		double growthRate = growthRateParameter.getParameterValue(0);
		
		if (usingGrowthRate) {
			expansion.setGrowthRate(growthRate);
		} else {
			double doublingTime = growthRate;
			growthRate = Math.log(2)/doublingTime;
			expansion.setDoublingTime(doublingTime);
		}
		
		expansion.setN0(N0);
		expansion.setProportion(N1);
		
		return expansion;
	}
	
	/**
	 * Parses an element from an DOM document into a ConstantExponentialModel. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return EXPANSION_MODEL; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			Type units = XMLParser.Utils.getUnitsAttr(xo);
				
			XMLObject cxo = xo.getChild(POPULATION_SIZE);
			Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = xo.getChild(ANCESTRAL_POPULATION_PROPORTION);
			Parameter N1Param = (Parameter)cxo.getChild(Parameter.class);
			
			Parameter rParam;
			boolean usingGrowthRate = true;
			
			if (xo.getChild(GROWTH_RATE) != null) {			
				cxo = xo.getChild(GROWTH_RATE);
				rParam = (Parameter)cxo.getChild(Parameter.class);
			} else {
				cxo = xo.getChild(DOUBLING_TIME);
				rParam = (Parameter)cxo.getChild(Parameter.class);
				usingGrowthRate = false;
			}
							
			return new ExpansionModel(N0Param, N1Param, rParam, units, usingGrowthRate);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A demographic model of constant population size followed by exponential growth.";
		}

		public Class getReturnType() { return ConstantExponentialModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private final XMLSyntaxRule[] rules = {
			XMLUnits.SYNTAX_RULES[0],
			new ElementRule(POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(ANCESTRAL_POPULATION_PROPORTION, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new XORRule(
		
				new ElementRule(GROWTH_RATE, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"A value of zero represents a constant population size, negative values represent decline towards the present, " +
					"positive numbers represents exponential growth towards the present. " + 
					"A random walk operator is recommended for this parameter with a starting value of 0.0 and no upper or lower limits." ),
				new ElementRule(DOUBLING_TIME, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"This parameter determines the doubling time.")
			)
		};
	};
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter N1Parameter = null;	
	Parameter growthRateParameter = null;	
	Expansion expansion = null;
	boolean usingGrowthRate = true;
}
