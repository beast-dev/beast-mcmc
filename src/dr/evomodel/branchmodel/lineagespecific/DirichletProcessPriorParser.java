/*
 * DirichletProcessPriorParser.java
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

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessPriorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_PRIOR = "dirichletProcessPrior";
	public static final String BASE_MODEL = "baseModel";
	public static final String CONCENTRATION = "concentration";
	public static final String CATEGORIES = "categories";
	
	@Override
	public String getParserName() {
		return DIRICHLET_PROCESS_PRIOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) xo.getElementFirstChild(BASE_MODEL);
		Parameter gamma = (Parameter) xo.getElementFirstChild(CONCENTRATION);

		return new DirichletProcessPrior(categoriesParameter, //
				uniquelyRealizedParameters, //
				baseModel, //
				gamma);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				
				new ElementRule(CATEGORIES,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments
				
				new ElementRule(CompoundParameter.class, false), // realized parameters
				
		        new ElementRule(BASE_MODEL,
                        new XMLSyntaxRule[] {
                                new ElementRule(ParametricMultivariateDistributionModel.class, 1, Integer.MAX_VALUE),
                        }
		        ), // base models
		        
		        new ElementRule(CONCENTRATION,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),// gamma
				
		};
	}

	@Override
	public String getParserDescription() {
		return DIRICHLET_PROCESS_PRIOR;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessPrior.class;
	}
	
}//END: class
