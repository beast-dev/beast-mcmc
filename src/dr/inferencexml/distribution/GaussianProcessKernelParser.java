/*
 * GaussianProcessKernelParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.distributions.gp.GaussianProcessKernel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class GaussianProcessKernelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "kernel";
    private static final String TYPE = "type";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        List<Parameter> parameters = new ArrayList<>();
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        if (parameter != null) {
            parameters.add(parameter);
        }

        final GaussianProcessKernel kernel;
        try {
            kernel = GaussianProcessKernel.factory(xo.getStringAttribute(TYPE), id, parameters);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }

        return kernel;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TYPE),
            new ElementRule(Parameter.class, true),
    };

    public String getParserDescription() { // TODO update
        return null;
    }

    public Class getReturnType() { return GaussianProcessKernel.class; }
}
