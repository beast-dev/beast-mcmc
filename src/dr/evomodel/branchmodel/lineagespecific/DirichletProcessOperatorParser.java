/*
 * DirichletProcessOperatorParser.java
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
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessOperatorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_OPERATOR = "dpOperator";
	public static final String DATA_LOG_LIKELIHOOD = "dataLogLikelihood";
	public static final String MH_STEPS = "mhSteps";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		DirichletProcessPrior dpp = (DirichletProcessPrior) xo.getChild(DirichletProcessPrior.class);
		Likelihood likelihood = (Likelihood) xo .getElementFirstChild(DATA_LOG_LIKELIHOOD);
		Parameter categoriesParameter = (Parameter) xo.getElementFirstChild(  DirichletProcessPriorParser.CATEGORIES);
		
		CountableRealizationsParameter allParameters = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);
		CompoundParameter uniquelyRealizedParameters = (CompoundParameter) xo.getChild(CompoundParameter.class);
		
		int M = xo.getIntegerAttribute(MH_STEPS);
		final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

		return new DirichletProcessOperator(dpp, // 
				categoriesParameter, //
				uniquelyRealizedParameters, //
				allParameters, //
				likelihood, //
				M, //
				weight //
				);
		
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {
		new ElementRule(DirichletProcessPrior.class, false),
		new ElementRule(CompoundParameter.class, false), //
		new ElementRule(CountableRealizationsParameter.class, false), //
		AttributeRule.newDoubleRule(MCMCOperator.WEIGHT) //
		};

	}// END: getSyntaxRules

	@Override
	public String getParserName() {
		return DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public String getParserDescription() {
		return DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessOperator.class;
	}

}// END: class
