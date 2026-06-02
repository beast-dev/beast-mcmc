/*
 * DefaultModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

import dr.inference.model.DefaultModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 *
 */
public class DefaultModelParser extends AbstractXMLObjectParser {

    public static final String DUMMY_MODEL = "dummyModel";

    public String getParserName() {
        return DUMMY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) {

        DefaultModel likelihood = new DefaultModel();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Parameter parameter = (Parameter) xo.getChild(i);
            likelihood.addVariable(parameter);
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A function wraps a component model that would otherwise not be registered with the MCMC. Always returns a log likelihood of zero.";
    }

    public Class getReturnType() {
        return DefaultModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
    };

}
