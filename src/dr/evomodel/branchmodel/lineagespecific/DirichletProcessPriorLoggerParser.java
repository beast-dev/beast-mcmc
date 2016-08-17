/*
 * DirichletProcessPriorLoggerParser.java
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

package dr.evomodel.branchmodel.lineagespecific;

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessPriorLoggerParser extends AbstractXMLObjectParser {

	public static final String DPP_LOGGER = "dppLogger";
	public static final String PRECISION = "precision";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//		ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) xo.getChild(ParametricMultivariateDistributionModel.class);
		Parameter precisionParameter = (Parameter)xo.getElementFirstChild(PRECISION); 
		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		Parameter categoriesParameter = (Parameter)xo.getElementFirstChild(DirichletProcessPriorParser.CATEGORIES); 
		
		return new DirichletProcessPriorLogger(precisionParameter, categoriesParameter, uniquelyRealizedParameters);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserName() {
		return DPP_LOGGER;
	}

	@Override
	public String getParserDescription() {
		return DPP_LOGGER;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessPriorLogger.class;
	}

}// END: class
