/*
 * ScaleOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ScaleOperator;
import dr.xml.*;

/**
 */
public class ScaleOperatorParser extends AbstractXMLObjectParser {
    public static final String SCALE_OPERATOR = "scaleOperator";
    public static final String SCALE_ALL = "scaleAll";
    public static final String SCALE_ALL_IND = "scaleAllIndependently";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String DEGREES_OF_FREEDOM = "df";
    public static final String INDICATORS = "indicators";
    public static final String PICKONEPROB = "pickoneprob";
    public static final String IGNORE_BOUNDS = "ignoreBounds";

    public String getParserName() {
        return SCALE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final boolean scaleAll = xo.getAttribute(SCALE_ALL, false);
        final boolean scaleAllInd = xo.getAttribute(SCALE_ALL_IND, false);
        final int degreesOfFreedom = xo.getAttribute(DEGREES_OF_FREEDOM, 0);
        final boolean ignoreBounds = xo.getAttribute(IGNORE_BOUNDS, false);

        final AdaptationMode mode = AdaptationMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
        }

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Bounds<Double> bounds = parameter.getBounds();

        Boolean negativeSpace = null;

        for (int dim = 0; dim < parameter.getDimension(); dim++) {

            if (parameter.getParameterValue(dim) < 0.0) {
                // the parameter is in negative space

                if (negativeSpace != null && !negativeSpace) {
                    throw new XMLParseException("Scale operator can only be used on parameters where all elements are strictly positive or negative (" + parameter.getId() + ")");
                }

                negativeSpace = true;

                if (bounds.getUpperLimit(dim) > 0.0) {
                    throw new XMLParseException("Scale operator can only be used on parameters constrained to be strictly positive or negative (" + parameter.getId() + ")");
                }
                if (!ignoreBounds && !Double.isInfinite(bounds.getLowerLimit(dim))) {
                    throw new XMLParseException("Scale operator can only be used on parameters with an infinite upper or lower bound (use a RandomWalk) (" + parameter.getId() + ")");
                }
            } else if (parameter.getParameterValue(dim) > 0.0) {
                if (negativeSpace != null && negativeSpace) {
                    throw new XMLParseException("Scale operator can only be used on parameters where all elements are strictly positive or negative (" + parameter.getId() + ")");
                }

                negativeSpace = false;

                if (bounds.getLowerLimit(dim) < 0.0) {
                    throw new XMLParseException("Scale operator can only be used on parameters constrained to be strictly positive or negative (" + parameter.getId() + ")");
                }
                if (!ignoreBounds && !Double.isInfinite(bounds.getUpperLimit(dim))) {
                    throw new XMLParseException("Scale operator can only be used on parameters with an infinite upper or lower bound (use a RandomWalk) (" + parameter.getId() + ")");
                }
            } else {
                // disallow zero starting values
                throw new XMLParseException("Scale operator can only be used on parameters where all elements are strictly positive or negative (" + parameter.getId() + ")");
            }
        }

        Parameter indicator = null;
        double indicatorOnProb = 1.0;
        final XMLObject inds = xo.getChild(INDICATORS);

        if (inds != null) {
            indicator = (Parameter) inds.getChild(Parameter.class);
            if (inds.hasAttribute(PICKONEPROB)) {
                indicatorOnProb = inds.getDoubleAttribute(PICKONEPROB);
                if (!(0 <= indicatorOnProb && indicatorOnProb <= 1)) {
                    throw new XMLParseException("pickoneprob must be between 0.0 and 1.0");
                }
            }
        }
        ScaleOperator operator = new ScaleOperator(parameter, scaleAll,
                degreesOfFreedom, scaleFactor,
                mode, indicator, indicatorOnProb,
                scaleAllInd);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a scale operator on a given parameter.";
    }

    public Class getReturnType() {
        return ScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newBooleanRule(SCALE_ALL, true),
            AttributeRule.newBooleanRule(SCALE_ALL_IND, true),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newIntegerRule(DEGREES_OF_FREEDOM, true),
            AttributeRule.newBooleanRule(IGNORE_BOUNDS, true),

            new ElementRule(Parameter.class),
            new ElementRule(INDICATORS,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(PICKONEPROB, true),
                            new ElementRule(Parameter.class)}, true),
    };
}
