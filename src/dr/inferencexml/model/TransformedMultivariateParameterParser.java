/*
 * TransformedMultivariateParameter.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.model;

import dr.inference.model.Parameter;
import dr.inference.model.TransformedMultivariateParameter;
import dr.util.Transform;
import dr.xml.*;

public class TransformedMultivariateParameterParser extends AbstractXMLObjectParser {

    private static final String TRANSFORMED_MULTIVARIATE_PARAMETER = "transformedMultivariateParameter";
    public static final String INVERSE = "inverse";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Transform.MultivariableTransform transform = (Transform.MultivariableTransform)
                xo.getChild(Transform.MultivariableTransform.class);
        final boolean inverse = xo.getAttribute(INVERSE, false);

        return new TransformedMultivariateParameter(parameter, transform, inverse);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(Transform.MultivariableTransform.class),
            AttributeRule.newBooleanRule(INVERSE, true),

    };

    public String getParserDescription() {
        return "A transformed multivariate parameter.";
    }

    public Class getReturnType() {
        return TransformedMultivariateParameter.class;
    }

    public String getParserName() {
        return TRANSFORMED_MULTIVARIATE_PARAMETER;
    }
}
