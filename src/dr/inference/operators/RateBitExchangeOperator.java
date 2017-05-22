/*
 * RateBitExchangeOperator.java
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

package dr.inference.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 */
public class RateBitExchangeOperator extends SimpleMCMCOperator {

    public static final String OPERATOR_NAME = "rateBitExchangeOperator";
    public static final String BITS = "bits";
    public static final String RATES = "rates";
    public static final String MAX_TRIES = "maxTries";
    private static final int MAX_TRIES_MAGIC_NUMBER = 1000;

    RateBitExchangeOperator(Parameter rateParameter, Parameter bitParameter, double weight,
                            int maxTries) {
        this.rateParameter = rateParameter;
        this.bitParameter = bitParameter;
        setWeight(weight);
        this.maxTries = maxTries;
    }

    public double doOperation() {
        int dim = rateParameter.getDimension() / 2;
        int tries = 0;

        // Find a pair-set in which at least one bit is non-zero
        int index = -1;
        do {
            index = MathUtils.nextInt(dim);
            ++tries;
        } while( (bitParameter.getParameterValue(index) + bitParameter.getParameterValue(index+dim) < 1) &&
                (tries < maxTries) );

        if (tries >= maxTries) {
//            throw new RuntimeException("Too many attempts");
            return Double.NEGATIVE_INFINITY;
        }

        // Swap (bit,rate) values
        double tmpBit = bitParameter.getParameterValue(index);
        double tmpRate = rateParameter.getParameterValue(index);
        bitParameter.setParameterValue(index, bitParameter.getParameterValue(index+dim));
        rateParameter.setParameterValue(index, rateParameter.getParameterValue(index+dim));
        bitParameter.setParameterValue(index+dim,tmpBit);
        rateParameter.setParameterValue(index+dim,tmpRate);

        return 0;
    }

        // Interface MCMCOperator
    public final String getOperatorName() {
        return OPERATOR_NAME + "(" + bitParameter.getParameterName() +"," + rateParameter.getParameterName() + ")";
    }

    public final String getPerformanceSuggestion() {
        return "No performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            Parameter ratesParameter = (Parameter) ((XMLObject)xo.getChild(RATES)).getChild(Parameter.class);
            Parameter bitsParameter = (Parameter) ((XMLObject)xo.getChild(BITS)).getChild(Parameter.class);

            int maxTries = xo.getAttribute(MAX_TRIES, MAX_TRIES_MAGIC_NUMBER);

            return new RateBitExchangeOperator(ratesParameter,bitsParameter, weight, maxTries);
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
                new ElementRule(BITS, new XMLSyntaxRule[] {
                new ElementRule(Parameter.class)
                }),
                new ElementRule(RATES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                AttributeRule.newIntegerRule(MAX_TRIES, true),
        };

    };

    // Private instance variables
    private Parameter bitParameter = null;
    private Parameter rateParameter = null;
    final private int maxTries;
}
