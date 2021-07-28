/*
 * TransformedParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.MatrixTransformedParameter;
import dr.inference.model.Parameter;
import dr.inference.model.TransformedParameter;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TransformedParameterParser extends AbstractXMLObjectParser {

    public static final String TRANSFORMED_PARAMETER = "transformedParameter";
    public static final String INVERSE = "inverse";
    public static final String AS_MATRIX = "asMatrix";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final Transform.ParsedTransform parsedTransform = (Transform.ParsedTransform) xo.getChild(Transform.ParsedTransform.class);
        final boolean inverse = xo.getAttribute(INVERSE, false);

        final boolean asMatrix = xo.getAttribute(AS_MATRIX, false);
        if (asMatrix) {
            if (parameter instanceof MatrixParameterInterface) {
                return new MatrixTransformedParameter((MatrixParameterInterface) parameter, parsedTransform.transform, inverse);
            } else {
                throw new XMLParseException("'asMatrix' is 'true' but the supplied parameter is not a matrix. " +
                        "Not currently implemented.");
            }
        }

        TransformedParameter transformedParameter = new TransformedParameter(parameter, parsedTransform.transform, inverse);
        return transformedParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(Transform.ParsedTransform.class),
            AttributeRule.newBooleanRule(INVERSE, true),
            AttributeRule.newBooleanRule(AS_MATRIX, true)
    };

    public String getParserDescription() {
        return "A transformed parameter.";
    }

    public Class getReturnType() {
        return TransformedParameter.class;
    }

    public String getParserName() {
        return TRANSFORMED_PARAMETER;
    }
}
