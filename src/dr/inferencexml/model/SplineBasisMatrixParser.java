/*
 * MatrixMatrixProductParser.java
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

import dr.inference.model.Parameter;
import dr.math.distributions.SplineBasisMatrix;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Pratyusa Datta
 */

public class SplineBasisMatrixParser extends AbstractXMLObjectParser {

    private final static String SPLINE_BASIS = "bSplineBasis";
    private final static String KNOTS = "knots";
    private final static String DEGREE = "degree";
    private final static String INTERCEPT = "includeIntercept";
    private final static String LOWER_BOUND = "lowerBound";
    private final static String UPPER_BOUND = "upperBound";
    private final static String ZERO_OUT_OF_BOUNDS = "zeroOutOfBounds";

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(KNOTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newIntegerRule(DEGREE, true),
            AttributeRule.newBooleanRule(INTERCEPT, true),
            AttributeRule.newDoubleRule(LOWER_BOUND, true),
            AttributeRule.newDoubleRule(UPPER_BOUND, true),
            AttributeRule.newBooleanRule(ZERO_OUT_OF_BOUNDS, true),
    };

    @Override
    public Object parseXMLObject (XMLObject xo) throws XMLParseException {

        Parameter evaluationPoints = (Parameter) xo.getChild(Parameter.class);

        Parameter knots = null;
        if (xo.hasChildNamed(KNOTS)) {
            knots = (Parameter) xo.getElementFirstChild(KNOTS);
        }

        int degree = xo.getAttribute(DEGREE, 3);
        boolean includeIntercept = xo.getAttribute(INTERCEPT, true);

        Double lowerBound = null;
        if (xo.hasAttribute(LOWER_BOUND)) {
            lowerBound = xo.getDoubleAttribute(LOWER_BOUND);
        }

        Double upperBound = null;
        if (xo.hasAttribute(UPPER_BOUND)) {
            upperBound = xo.getDoubleAttribute(UPPER_BOUND);
        }

        boolean zeroOfOutBounds = xo.getAttribute(ZERO_OUT_OF_BOUNDS, false);

        return new SplineBasisMatrix(xo.getId(), evaluationPoints, knots,
                degree, includeIntercept, lowerBound, upperBound, zeroOfOutBounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules () {
        return rules;
    }

    @Override
    public String getParserDescription () {
        return "Gets Latent Factor Model to return data with residuals computed";
    }

    @Override
    public Class getReturnType () {
        return SplineBasisMatrix.class;
    }

    @Override
    public String getParserName () {
        return SPLINE_BASIS;
    }
}
