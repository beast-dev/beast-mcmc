/*
 * ConstantLogisticModel.java
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

import dr.evolution.coalescent.ConstLogistic;
import dr.evolution.coalescent.DemographicFunction;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Logistic growth from a constant ancestral population size.
 * 
 * @version $Id: ConstantLogisticModel.java,v 1.7 2005/04/11 11:24:39 alexei Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class ConstantLogisticModel extends DemographicModel
{
	
	//
	// Public stuff
	//
	
	public static String CONSTANT_LOGISTIC_MODEL = "constantLogistic";
	public static String POPULATION_SIZE = "populationSize";
	public static String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";

	public static String GROWTH_RATE = "growthRate";
	public static String SHAPE = "shape";
	public static String ALPHA = "alpha";
	

	/**
	 * Construct demographic model with default settings
	 */
	public ConstantLogisticModel(Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, Type units) {
	
		this(CONSTANT_LOGISTIC_MODEL, N0Parameter, N1Parameter, growthRateParameter, shapeParameter, alpha, units);
	}

	/**
	 * Construct demographic model with default settings
	 */
	public ConstantLogisticModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, Type units) {
	
		super(name);
		
		constLogistic = new ConstLogistic(units);
		
		this.N0Parameter = N0Parameter;
		addParameter(N0Parameter);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.N1Parameter = N1Parameter;
		addParameter(N1Parameter);
		N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.growthRateParameter = growthRateParameter;
		addParameter(growthRateParameter);
		growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.shapeParameter = shapeParameter;
		addParameter(shapeParameter);
		shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.alpha = alpha;
		
		setUnits(units);
	}


	// general functions

	public DemographicFunction getDemographicFunction() {
		constLogistic.setN0(N0Parameter.getParameterValue(0));
		constLogistic.setN1(N1Parameter.getParameterValue(0));
		
		double r = growthRateParameter.getParameterValue(0);
		constLogistic.setGrowthRate(r);
		
		// AER 24/02/03 
		// logisticGrowth.setShape(Math.exp(shapeParameter.getParameterValue(0)));
		
		// New parameterization of logistic shape to be the time at which the 
		// population reached some proportion alpha:
		double C = ((1.0 - alpha) * Math.exp(- r * shapeParameter.getParameterValue(0))) / alpha;
		constLogistic.setShape(C);
		
		return constLogistic;
	}
	
	/**
	 * Parses an element from an DOM document into a ExponentialGrowth. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return CONSTANT_LOGISTIC_MODEL; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			Type units = XMLParser.Utils.getUnitsAttr(xo);
				
			XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
			Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(ANCESTRAL_POPULATION_SIZE);
			Parameter N1Param = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(GROWTH_RATE);
			Parameter rParam = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject)xo.getChild(SHAPE);
			Parameter cParam = (Parameter)cxo.getChild(Parameter.class);
			
			double alpha = xo.getDoubleAttribute(ALPHA);
				
			return new ConstantLogisticModel(N0Param, N1Param, rParam, cParam, alpha, units);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A demographic model of constant population size followed by logistic growth.";
		}

		public Class getReturnType() { return ConstantLogisticModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			XMLUnits.SYNTAX_RULES[0],
			new ElementRule(POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(ANCESTRAL_POPULATION_SIZE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(GROWTH_RATE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(SHAPE, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			AttributeRule.newDoubleRule(ALPHA)
		};
	};
	//
	// protected stuff
	//

	Parameter N0Parameter = null;	
	Parameter N1Parameter = null;	
	Parameter growthRateParameter = null;	
	Parameter shapeParameter = null;	
	double alpha = 0.5;
	ConstLogistic constLogistic = null;
}
