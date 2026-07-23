/*
 * RejectionOperator.java
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

package dr.inference.operators.rejection;

import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

public class RejectionOperator extends SimpleMCMCOperator {

    public interface RejectionProvider extends GibbsOperator {
        double[] getProposedUpdate();

        Parameter getParameter();
    }


    private final RejectionProvider gibbsOp;
    private final AcceptCondition condition;

    public RejectionOperator(RejectionProvider gibbsOp, AcceptCondition condition, double weight) {
        setWeight(weight);
        this.gibbsOp = gibbsOp;
        this.condition = condition;
    }

    @Override
    public String getOperatorName() {
        return REJECTION_OPERATOR;
    }

    @Override
    public double doOperation() {
        double[] values = gibbsOp.getProposedUpdate();
        if (condition.satisfiesCondition(values)) {

            Parameter parameter = gibbsOp.getParameter();
            for (int i = 0; i < parameter.getDimension(); i++) {
                parameter.setParameterValueQuietly(i, values[i]);
            }

            parameter.fireParameterChangedEvent();

            return Double.POSITIVE_INFINITY;
        }
        return Double.NEGATIVE_INFINITY;
    }


    private static final String REJECTION_OPERATOR = "rejectionOperator";
    private static final String CONDITION = "condition";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            RejectionProvider gibbsOp = (RejectionProvider) xo.getChild(RejectionProvider.class);
            double weight = xo.getAttribute(WEIGHT, gibbsOp.getWeight());

            AcceptCondition condition = null;
            if (xo.hasAttribute(CONDITION)) {
                String stringCondition = xo.getStringAttribute(CONDITION);

                for (AcceptCondition.SimpleAcceptCondition simpleCondition : AcceptCondition.SimpleAcceptCondition.values()) {
                    if (stringCondition.equalsIgnoreCase(simpleCondition.getName())) {
                        condition = simpleCondition;
                        break;
                    }
                }
                if (condition == null) {
                    throw new XMLParseException("Unrecognized condition type: " + stringCondition);
                }
            } else {
                condition = (AcceptCondition) xo.getChild(AcceptCondition.class);
            }
            return new RejectionOperator(gibbsOp, condition, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(RejectionProvider.class),
                    new XORRule(
                            AttributeRule.newStringRule(CONDITION),
                            new ElementRule(AcceptCondition.class)
                    ),
                    AttributeRule.newDoubleRule(WEIGHT, true)
            };
        }

        @Override
        public String getParserDescription() {
            return "A rejection sampler that always accepts a Gibbs proposal if a condition is met (and always rejects otherwise)";
        }

        @Override
        public Class getReturnType() {
            return RejectionOperator.class;
        }

        @Override
        public String getParserName() {
            return REJECTION_OPERATOR;
        }
    };


}
