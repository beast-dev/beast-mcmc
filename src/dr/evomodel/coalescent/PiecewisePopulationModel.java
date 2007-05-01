/*
 * PiecewisePopulationModel.java
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
import dr.evolution.coalescent.PiecewiseConstantPopulation;
import dr.evolution.coalescent.PiecewiseLinearPopulation;
import dr.evoxml.XMLUnits;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 * @version $Id: PiecewisePopulationModel.java,v 1.13 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewisePopulationModel extends DemographicModel
{
	//
	// Public stuff
	//
	
	public static String PIECEWISE_POPULATION = "piecewisePopulation";
	public static String EPOCH_SIZES = "epochSizes";
	public static String EPOCH_WIDTHS = "epochWidths";

	/**
	 * Construct demographic model with default settings
	 */
	public PiecewisePopulationModel(Parameter N0Parameter, double[] epochLengths, boolean linear, int units) {
	
		this(PIECEWISE_POPULATION, N0Parameter, epochLengths, linear, units);
		
	}

	/**
	 * Construct demographic model with default settings
	 */
	public PiecewisePopulationModel(String name, Parameter N0Parameter,  double[] epochLengths, boolean linear, int units) {
	
		super(name);
		
		if (N0Parameter.getDimension() != (epochLengths.length + 1)) {
			throw new IllegalArgumentException(
				"epochSize parameter must have one less components than the number of epochs: (" + (epochLengths.length + 1) + 
				") but instead has " + N0Parameter.getDimension() + "!"
			);
		}
		
		this.N0Parameter = N0Parameter;
		this.epochLengths = epochLengths;
		addParameter(N0Parameter);
		//addParameter(epochLengths);
		N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, N0Parameter.getDimension()));
		//epochLengths.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, epochLengths.getDimension()));
		
		setUnits(units);
		
		if (linear) {
			piecewiseFunction = new PiecewiseLinearPopulation(epochLengths, new double[N0Parameter.getDimension()], units);
		} else {
			piecewiseFunction = new PiecewiseConstantPopulation(epochLengths, new double[N0Parameter.getDimension()], units);
		}
	}
	
	/**
	 * 
	 */
	public DemographicFunction getDemographicFunction() {
		for (int i = 0; i < N0Parameter.getDimension(); i++) {
			piecewiseFunction.setArgument(i, N0Parameter.getParameterValue(i));
		}
		return piecewiseFunction;
	}
	
	// **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************
	
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no intermediates need to be recalculated...
	}
	
	protected void handleParameterChangedEvent(Parameter parameter, int index) {

		if (parameter == N0Parameter) {
			//System.out.println("popSize parameter changed..");
		}
		
		// no intermediates need to be recalculated...
	}
	
	protected void storeState() {} // no additional state needs storing
	protected void restoreState() {} // no additional state needs restoring	
	protected void acceptState() {} // no additional state needs accepting	
	protected void adoptState(Model source) {} // no additional state needs adopting	
	
	/**
	 * Parses an element from an DOM document into a PiecewisePopulation. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return PIECEWISE_POPULATION; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			int units = XMLParser.Utils.getUnitsAttr(xo);
				
			XMLObject cxo = (XMLObject)xo.getChild(EPOCH_SIZES);
			Parameter epochSizes = (Parameter)cxo.getChild(Parameter.class);
	
			
			XMLObject obj = (XMLObject)xo.getChild(EPOCH_WIDTHS);
			double[] epochWidths = obj.getDoubleArrayAttribute("widths");
			
			boolean isLinear =  xo.getBooleanAttribute("linear");
					
			return new PiecewisePopulationModel(epochSizes, epochWidths, isLinear, units);
		}
	
	
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents a piecewise population model";
		}

		public Class getReturnType() { return PiecewisePopulationModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(EPOCH_SIZES, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(EPOCH_WIDTHS, 
				new XMLSyntaxRule[] { AttributeRule.newDoubleArrayRule("widths", false), }),
			XMLUnits.SYNTAX_RULES[0],
			AttributeRule.newBooleanRule("linear")
		};
	};


	//
	// protected stuff
	//
	
	Parameter N0Parameter;
	double[] epochLengths;
	DemographicFunction piecewiseFunction = null;
}
