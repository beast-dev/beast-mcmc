/*
 * GaussianProcessKernelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.distributions.gp.GaussianProcessKernel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class GaussianProcessKernelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "kernel";
    private static final String TYPE = "type";
    private static final String SCALE = "scale";
    private static final String LENGTH = "length";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        List<Parameter> parameters = new ArrayList<>();

        Parameter scale = xo.hasChildNamed(SCALE) ?
                (Parameter) xo.getElementFirstChild(SCALE) :
                new Parameter.Default(1.0);

        Parameter length = xo.hasChildNamed(LENGTH) ?
                (Parameter) xo.getElementFirstChild(LENGTH) :
                new Parameter.Default(1.0);

        parameters.add(scale);
        parameters.add(length);

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
            new ElementRule(SCALE, Parameter.class, "",  true),
            new ElementRule(LENGTH, Parameter.class, "", true),
    };

    public String getParserDescription() { // TODO update
        return null;
    }

    public Class getReturnType() { return GaussianProcessKernel.class; }
}
