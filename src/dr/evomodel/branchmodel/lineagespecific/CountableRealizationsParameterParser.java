/*
 * CountableRealizationsParameterParser.java
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
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 */
public class CountableRealizationsParameterParser extends
		AbstractXMLObjectParser {

	public static final String COUNTABLE_REALIZATIONS_PARAMETER = "countableRealizationsParameter";

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Parameter categoriesParameter = (Parameter) xo .getElementFirstChild(DirichletProcessPriorParser.CATEGORIES);
		CompoundParameter realizedParameters = (CompoundParameter) xo .getChild(CompoundParameter.class);

		return new CountableRealizationsParameter(categoriesParameter,
				realizedParameters);
	}

	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

				new ElementRule(DirichletProcessPriorParser.CATEGORIES,
						new XMLSyntaxRule[] { new ElementRule(Parameter.class,
								false) }), // categories assignments

				new ElementRule(CompoundParameter.class, false) // realized parameters

		};
	}

	public String getParserDescription() {
		return COUNTABLE_REALIZATIONS_PARAMETER;
	}

	public Class getReturnType() {
		return Parameter.class;
	}

	public String getParserName() {
		return COUNTABLE_REALIZATIONS_PARAMETER;
	}
}
