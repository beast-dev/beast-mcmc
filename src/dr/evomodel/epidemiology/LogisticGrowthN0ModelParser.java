/*
 * LogisticGrowthN0ModelParser.java
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

import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into LogisticGrowthN0Model.
 * @author Daniel Wilson
 */
public class LogisticGrowthN0ModelParser extends AbstractXMLObjectParser
{
	
	public static String POPULATION_SIZE = "populationSize";
	public static String LOGISTIC_GROWTH_MODEL = "logisticGrowthN0";
	
	public static String GROWTH_RATE = "growthRate";
	public static String DOUBLING_TIME = "doublingTime";
	public static String TIME_50 = "t50";
	
	public String getParserName() { return LOGISTIC_GROWTH_MODEL; }
	
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
		
		XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
		Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);		
		
		boolean usingGrowthRate = true;		
		Parameter rParam = null;		
		if (xo.getChild(GROWTH_RATE) != null) {			
			cxo = (XMLObject)xo.getChild(GROWTH_RATE);
			rParam = (Parameter)cxo.getChild(Parameter.class);
		} else {
			cxo = (XMLObject)xo.getChild(DOUBLING_TIME);
			rParam = (Parameter)cxo.getChild(Parameter.class);
			usingGrowthRate = false;
		}
		
		cxo = (XMLObject)xo.getChild(TIME_50);
		Parameter t50Param = (Parameter)cxo.getChild(Parameter.class);
		
		return new LogisticGrowthN0Model(N0Param, rParam, t50Param, units, usingGrowthRate);
	}
	
	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************
	
	public String getParserDescription() {
		return "Logistic growth demographic model.";
	}
	
	public Class getReturnType() { return LogisticGrowthN0Model.class; }
	
	public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
	XMLUnits.SYNTAX_RULES[0],
	new ElementRule(POPULATION_SIZE, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"This parameter represents the present day population size."),
	new XORRule(
				new ElementRule(GROWTH_RATE, 
								new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
								"This parameter determines the rate of growth during the exponential phase. " +
                                        "See exponentialGrowth for details."),
				new ElementRule(DOUBLING_TIME, 
								new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
								"This parameter determines the doubling time at peak growth rate.")
				),
	new ElementRule(TIME_50, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"This parameter represents the time in the past (or future) when the population had (or will have) " + 
					"half of the carrying capacity (maximum population size). " +
					"It is a positive number if half the carrying capacity was attained in the past " +
					"or a negative number if it will be attained in the future.")
	};
}
