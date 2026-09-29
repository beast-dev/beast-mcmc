/*
 * FiniteMixtureParameterParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inferencexml.model;

import dr.inference.model.FiniteMixtureParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class FiniteMixtureParameterParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "finiteMixtureParameter";
    private static final String VALUES = "values";
    private static final String CATEGORIES = "categories";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter values = (Parameter) xo.getElementFirstChild(VALUES);
        Parameter categories = (Parameter) xo.getElementFirstChild(CATEGORIES);

        return new FiniteMixtureParameter(xo.getId(), values, categories);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(VALUES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
            new ElementRule(CATEGORIES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
    };

    public String getParserDescription() {
        return "A finite mixture parameter.";
    }

    public Class getReturnType() {
        return FiniteMixtureParameter.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
