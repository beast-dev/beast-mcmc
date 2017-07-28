/*
 * ApplyOperatorParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 */
public class ApplyOperatorParser extends AbstractXMLObjectParser {

    public static final String APPLY_OPERATOR = "jitter";
    public static final String TIMES = "times";
    public static final String SCALE = "scale";
;

        public String getParserName() {
            return APPLY_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double scale = xo.getAttribute(SCALE, 0.1);

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            NormalDistribution normal = new NormalDistribution(0, 1);

            Bounds bounds = parameter.getBounds();

            for (int i = 0; i < parameter.getDimension(); ++i) {

                double value = parameter.getParameterValue(i);
                value += scale * (Double) normal.nextRandom();

                if (bounds != null) {
                    final double lower = (Double) bounds.getLowerLimit(i);
                    final double upper = (Double) bounds.getUpperLimit(i);

                    value = RandomWalkOperator.reflectValue(value, lower, upper);
                }

                parameter.setParameterValue(i, value);
            }
            return null;
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
                AttributeRule.newIntegerRule(TIMES, true),
                AttributeRule.newDoubleRule(SCALE, true),
//                new ElementRule(SimpleMCMCOperator.class),
                new ElementRule(Parameter.class),
        };
}
