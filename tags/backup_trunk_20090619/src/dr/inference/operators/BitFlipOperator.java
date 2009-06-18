/*
 * BitFlipOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * A generic operator that flips bits.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class BitFlipOperator extends SimpleMCMCOperator {

    public static final String BIT_FLIP_OPERATOR = "bitFlipOperator";
    public static final String BITS = "bits";
    public static final String USES_SUM_PRIOR = "usesPriorOnSum";

    public BitFlipOperator(Parameter parameter, double weight, boolean usesPriorOnSum) {
//       this(parameter,weight,1,usesPriorOnSum);
//    }
//
//    public BitFlipOperator(Parameter parameter, double weight, int bits, boolean usesPriorOnSum) {
        this.parameter = parameter;
//        this.bits = bits;
        this.usesPriorOnSum = usesPriorOnSum;
        setWeight(weight);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * Change the parameter and return the hastings ratio.
     * Flip (Switch a 0 to 1 or 1 to 0) for a random bit in a bit vector.
     * Return the hastings ratio which makes all subsets of vectors with the same number of 1 bits
     * equiprobable, unless usesPriorOnSum = false then all configurations are equiprobable
     */
    public final double doOperation() {
        final int dim = parameter.getDimension();
        double sum = 0.0;

        if(usesPriorOnSum) {
            for (int i = 0; i < dim; i++) {
                sum += parameter.getParameterValue(i);
            }
        }

        final int pos = MathUtils.nextInt(dim);

        final int value = (int) parameter.getParameterValue(pos);
        double logq = 0.0;
        if (value == 0) {
            parameter.setParameterValue(pos, 1.0);

            if(usesPriorOnSum)
                logq = -Math.log((dim - sum) / (sum + 1));

        } else if (value == 1) {
            parameter.setParameterValue(pos, 0.0);
            if(usesPriorOnSum)
                logq = -Math.log(sum / (dim - sum + 1));

        } else {
            throw new RuntimeException("expected 1 or 0");
        }

        // hastings ratio is designed to make move symmetric on sum of 1's
        return logq;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "bitFlip(" + parameter.getParameterName() + ")";
    }

    public final String getPerformanceSuggestion() {
        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIT_FLIP_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            boolean usesPriorOnSum = xo.getAttribute(USES_SUM_PRIOR,true);

            return new BitFlipOperator(parameter, weight, usesPriorOnSum);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a bit-flip operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
//                AttributeRule.newIntegerRule(BITS,true),
                AttributeRule.newBooleanRule(USES_SUM_PRIOR,true),
                new ElementRule(Parameter.class)
        };

    };
    // Private instance variables

    private Parameter parameter = null;
//    private int bits;
    private boolean usesPriorOnSum = true;
}
