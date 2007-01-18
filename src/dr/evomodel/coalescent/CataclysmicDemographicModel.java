/*
 * CataclysmicDemographicModel.java
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

import dr.evolution.coalescent.CataclysmicDemographic;
import dr.evolution.coalescent.DemographicFunction;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * This class models an exponentially growing model that suffers a
 * cataclysmic event and goes into exponential decline
 * This model is nested with the constant-population size model (r=0).
 * 
 * @version $Id: CataclysmicDemographicModel.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class CataclysmicDemographicModel extends DemographicModel
{
	
	//
	// Public stuff
	//
	
	public static String POPULATION_SIZE = "populationSize";
	public static String GROWTH_RATE = "growthRate";
	public static String SPIKE_SIZE = "spikeFactor";
	public static String TIME_OF_CATACLYSM = "timeOfCataclysm";
	
	public static String CATACLYSM_MODEL = "cataclysm";

	

	/**
	 * Construct demographic model with default settings
	 */
	public CataclysmicDemographicModel( Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter, Type units) {
	
		this(CATACLYSM_MODEL, N0Parameter, N1Parameter, growthRateParameter,  timeParameter, units);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public CataclysmicDemographicModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter, Type units) {
	
		super(name);
		
		cataclysm = new CataclysmicDemographic(units);
		
		this.N0Parameter = N0Parameter;
		addParameter(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		
		this.N1Parameter = N1Parameter;
		addParameter(N1Parameter);
		N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addParameter(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));

		this.timeParameter = timeParameter;
		addParameter(timeParameter);
		timeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));


		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		cataclysm.setN0(N0Parameter.getParameterValue(0));
		cataclysm.setGrowthRate(growthRateParameter.getParameterValue(0));
		cataclysm.setCataclysmTime(timeParameter.getParameterValue(0));
		
		// Doesn't this...
		/*
		double N0 = N0Parameter.getParameterValue(0);
		double N1 = N1Parameter.getParameterValue(0) * N0;
		double t = timeParameter.getParameterValue(0);
		double declineRate = Math.log(N1/N0)/t;
		*/ // ..collapse to...
		
		double t = timeParameter.getParameterValue(0);
		double declineRate = Math.log(N1Parameter.getParameterValue(0))/t;
		cataclysm.setDeclineRate(declineRate);
		
		
		return cataclysm;
	}
	
	/**
	 * Parses an element from an DOM document into a ExponentialGrowth. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return CATACLYSM_MODEL; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			Type units = XMLParser.Utils.getUnitsAttr(xo);
			
			XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
			Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(GROWTH_RATE);
			Parameter rParam = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(SPIKE_SIZE);
			Parameter N1Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(TIME_OF_CATACLYSM);
			Parameter tParam = (Parameter)cxo.getChild(Parameter.class);
			
			return new CataclysmicDemographicModel(N0Param, N1Param, rParam, tParam, units);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A demographic model of exponential growth.";
		}

		public Class getReturnType() { return CataclysmicDemographicModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(GROWTH_RATE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The rate of exponential growth before the cataclysmic event." ),
			new ElementRule(SPIKE_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The factor larger the population size was at its height." ),
			new ElementRule(TIME_OF_CATACLYSM, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The time of the cataclysmic event that lead to exponential decline." ),
			XMLUnits.SYNTAX_RULES[0]
		};	
	};
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter N1Parameter = null;	
	Parameter growthRateParameter = null;	
	Parameter timeParameter = null;	
	CataclysmicDemographic cataclysm = null;
}
