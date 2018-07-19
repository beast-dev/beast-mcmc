/*
 * RandomWalkOperatorParser.java
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
import dr.inference.model.TransformedMultivariateParameter;
import dr.inference.model.TransformedParameter;
import dr.inference.operators.*;
import dr.util.Transform;
import dr.xml.*;

import static dr.util.Transform.parseTransform;

/**
 */
public class RandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_OPERATOR = "randomWalkOperator";
    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";

    public static final String BOUNDARY_CONDITION = "boundaryCondition";

    public static final String INVERSE = "inverse";

        public String getParserName() {
            return RANDOM_WALK_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            if (xo.hasAttribute(LOWER) || xo.hasAttribute(UPPER)) {
                throw new XMLParseException("Do not provide lower/upper bounds on for a RandomWalkOperator; set these values are parameter bounds");
            }

            RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.valueOf(
                    xo.getAttribute(BOUNDARY_CONDITION, RandomWalkOperator.BoundaryCondition.reflecting.name()));


            final Bounds<Double> bounds = parameter.getBounds();
            final int dim = parameter.getDimension();

            boolean lowerBoundsSet = true;
            boolean upperBoundsSet = true;
            for (int i = 0; i < dim; ++i) {
                if (bounds.getLowerLimit(i) == null || Double.isInfinite(bounds.getLowerLimit(i))) {
                    lowerBoundsSet = false;
                }
                if (bounds.getUpperLimit(i) == null || Double.isInfinite(bounds.getUpperLimit(i))) {
                    upperBoundsSet = false;
                }
            }

            if (condition == RandomWalkOperator.BoundaryCondition.logit) {
                if (!lowerBoundsSet || !upperBoundsSet) {
                    throw new XMLParseException("The logit transformed RandomWalkOperator cannot be used on a parameter without bounds.");
                }
            }

            if (condition == RandomWalkOperator.BoundaryCondition.log) {
                if (!lowerBoundsSet) {
                    throw new XMLParseException("The log transformed RandomWalkOperator cannot be used on a parameter without lower bounds.");
                }
            }

            RandomWalkOperator randomWalk;

            if (xo.hasChildNamed(UPDATE_INDEX)) {
                XMLObject cxo = xo.getChild(UPDATE_INDEX);
                Parameter updateIndex = (Parameter) cxo.getChild(Parameter.class);
                if (updateIndex.getDimension() != parameter.getDimension())
                    throw new RuntimeException("Parameter to update and missing indices must have the same dimension");
                randomWalk = new RandomWalkOperator(parameter, updateIndex, windowSize, condition,
                        weight, mode);
            } else {
                randomWalk = new RandomWalkOperator(parameter, null, windowSize, condition, weight, mode);
            }

            final Transform transform = parseTransform(xo);

            if (transform == null) {
                return randomWalk;
            } else {
                final boolean inverse = xo.getAttribute(INVERSE, false);
                TransformedParameter transformedParameter;
                if (transform.isMultivariate()) {
                    transformedParameter
                            = new TransformedMultivariateParameter(parameter,
                            (Transform.MultivariableTransform) transform,
                            inverse);
                } else {
                    transformedParameter = new TransformedParameter(parameter, transform, inverse);
                }
                return new TransformedParameterRandomWalkOperator(transformedParameter, randomWalk);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
                new ElementRule(UPDATE_INDEX,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        },true),
                new StringAttributeRule(BOUNDARY_CONDITION, null, RandomWalkOperator.BoundaryCondition.values(), true),
                new ElementRule(Parameter.class)
        };
}
