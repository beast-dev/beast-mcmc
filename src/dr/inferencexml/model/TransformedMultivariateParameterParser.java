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

import dr.inference.model.*;
import dr.util.Transform;
import dr.xml.*;

public class TransformedMultivariateParameterParser extends AbstractXMLObjectParser {

    private static final String TRANSFORMED_MULTIVARIATE_PARAMETER = "transformedMultivariateParameter";
    public static final String INVERSE = "inverse";
    private static final String BOUNDS = "bounds";
    private static final String AS_MATRIX = "asMatrix";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Transform.MultivariableTransform transform = (Transform.MultivariableTransform)
                xo.getChild(Transform.MultivariableTransform.class);
        final boolean inverse = xo.getAttribute(INVERSE, false);

        final TransformedMultivariateParameter transformedParameter;
        final boolean asMatrix = xo.getAttribute(AS_MATRIX, false);
        if (asMatrix) {
            if (parameter instanceof MatrixParameterInterface) {
                transformedParameter = new TransformedMatrixParameter((MatrixParameterInterface) parameter, transform, inverse);
            } else {
                throw new XMLParseException("'asMatrix' is 'true' but the supplied parameter is not a matrix. " +
                        "Not currently implemented.");
            }
        } else {
            transformedParameter = new TransformedMultivariateParameter(parameter, transform, inverse);
        }

        if (xo.hasChildNamed(BOUNDS)) {
            Bounds<Double> bounds = ((Parameter) xo.getElementFirstChild(BOUNDS)).getBounds();
            transformedParameter.addBounds(bounds);
        } else {
            transformedParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,
                    Double.NEGATIVE_INFINITY, parameter.getDimension()));
        }
        return transformedParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(Transform.MultivariableTransform.class),
            AttributeRule.newBooleanRule(INVERSE, true),
            AttributeRule.newBooleanRule(AS_MATRIX, true),

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
