/*
 * TransformedRandomWalkOperatorParser.java
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

import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.TransformedRandomWalkOperator;
import dr.util.Transform;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the TransformedRandomWalkOperator.
 *
 * @author Guy Baele
 */
public class TransformedRandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String TRANSFORMED_RANDOM_WALK_OPERATOR = "transformedRandomWalkOperator";
    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";

    public static final String BOUNDARY_CONDITION = "boundaryCondition";

        public String getParserName() {
            return TRANSFORMED_RANDOM_WALK_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AdaptationMode mode = AdaptationMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            int dim = parameter.getDimension();
            Transform[] transformations = new Transform[dim];
            for (int i = 0; i < dim; i++) {
                transformations[i] = Transform.NONE;
            }

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Transform.ParsedTransform) {
                    Transform.ParsedTransform thisObject = (Transform.ParsedTransform) child;

                    System.err.println("Transformations:");
                    for (int j = thisObject.start; j < thisObject.end; ++j) {
                        transformations[j] = thisObject.transform;
                        System.err.print(transformations[j].getTransformName() + " ");
                    }
                    System.err.println();
                }
            }
            
            Double lower = null;
            Double upper = null;

            if (xo.hasAttribute(LOWER)) {
                lower = xo.getDoubleAttribute(LOWER);
            }

            if (xo.hasAttribute(UPPER)) {
                upper = xo.getDoubleAttribute(UPPER);
            }

            TransformedRandomWalkOperator.BoundaryCondition condition = TransformedRandomWalkOperator.BoundaryCondition.valueOf(
                    xo.getAttribute(BOUNDARY_CONDITION, TransformedRandomWalkOperator.BoundaryCondition.reflecting.name()));

            if (xo.hasChildNamed(UPDATE_INDEX)) {
                XMLObject cxo = xo.getChild(UPDATE_INDEX);
                Parameter updateIndex = (Parameter) cxo.getChild(Parameter.class);
                if (updateIndex.getDimension() != parameter.getDimension())
                    throw new RuntimeException("Parameter to update and missing indices must have the same dimension");
                return new TransformedRandomWalkOperator(parameter, transformations, updateIndex, windowSize, condition,
                        weight, mode, lower, upper);
            }

            return new TransformedRandomWalkOperator(parameter, transformations, null, windowSize, condition, weight, mode, lower, upper);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a transformed random walk operator on a given parameter.";
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
                AttributeRule.newDoubleRule(LOWER, true),
                AttributeRule.newDoubleRule(UPPER, true),
                AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
                new ElementRule(UPDATE_INDEX,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        },true),
                new StringAttributeRule(BOUNDARY_CONDITION, null, TransformedRandomWalkOperator.BoundaryCondition.values(), true),
                new ElementRule(Parameter.class),
                new ElementRule(Transform.ParsedTransform.class, 0, Integer.MAX_VALUE)
        };
}
