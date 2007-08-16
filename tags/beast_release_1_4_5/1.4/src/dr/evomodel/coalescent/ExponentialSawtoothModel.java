/*
 * ExponentialSawtoothModel.java
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
import dr.evolution.coalescent.ExponentialSawtooth;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * 
 * @version $Id: ExponentialSawtoothModel.java,v 1.4 2005/04/11 11:24:39 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class ExponentialSawtoothModel extends DemographicModel
{
	
	//
	// Public stuff
	//
	
	public static String POPULATION_SIZE = "populationSize";
	public static String GROWTH_RATE = "growthRate";
	public static String WAVELENGTH = "wavelength";
	public static String OFFSET = "offset";
	
	public static String EXPONENTIAL_SAWTOOTH = "exponentialSawtooth";

	

	/**
	 * Construct demographic model with default settings
	 */
	public ExponentialSawtoothModel( Parameter N0Parameter, Parameter growthRateParameter, Parameter wavelengthParameter, Parameter offsetParameter, int units) {
	
		this(EXPONENTIAL_SAWTOOTH, N0Parameter, growthRateParameter, wavelengthParameter, offsetParameter, units);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public ExponentialSawtoothModel(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter wavelengthParameter, Parameter offsetParameter, int units) {
	
		super(name);
		
		expSaw = new ExponentialSawtooth(units);
		
		this.N0Parameter = N0Parameter;
		addParameter(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));
		
		this.growthRateParameter = growthRateParameter;
		addParameter(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, -Double.MAX_VALUE, 1));

		this.wavelengthParameter = wavelengthParameter;
		addParameter(wavelengthParameter);
		wavelengthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));
		
		this.offsetParameter = offsetParameter;
		addParameter(offsetParameter);
		offsetParameter.addBounds(new Parameter.DefaultBounds(1.0, -1.0, 1));

		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		expSaw.setN0(N0Parameter.getParameterValue(0));
		expSaw.setGrowthRate(growthRateParameter.getParameterValue(0));
		expSaw.setWavelength(wavelengthParameter.getParameterValue(0));
		
		double offset = offsetParameter.getParameterValue(0);
		if (offset < 0.0) {
			offset += 1.0;
		}
		expSaw.setOffset(offset);
		
		return expSaw;
	}
	
	/**
	 * Parses an element from an DOM document into a ExponentialGrowth. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return EXPONENTIAL_SAWTOOTH; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			int units = XMLParser.Utils.getUnitsAttr(xo);
			
			XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
			Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(GROWTH_RATE);
			Parameter rParam = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(WAVELENGTH);
			Parameter wavelengthParam = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(OFFSET);
			Parameter tParam = (Parameter)cxo.getChild(Parameter.class);
			
			return new ExponentialSawtoothModel(N0Param, rParam, wavelengthParam, tParam, units);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A demographic model of succesive exponential growth and periodic population crashes.";
		}

		public Class getReturnType() { return TwoEpochDemographicModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(GROWTH_RATE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The rate of exponential growth during the growth phase." ),
			new ElementRule(WAVELENGTH, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The wavelength between successive population crashes." ),
			new ElementRule(OFFSET, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"The proportion of the last growth phase that is not achieved at the final sample time." ),
			XMLUnits.SYNTAX_RULES[0]
		};	
	};
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter growthRateParameter = null;	
	Parameter wavelengthParameter = null;	
	Parameter offsetParameter = null;	
	ExponentialSawtooth expSaw = null;
}
