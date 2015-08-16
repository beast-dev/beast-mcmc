/*
 * ExponentialMarkovModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class ExponentialMarkovModelParser extends AbstractXMLObjectParser {

    public static final String CHAIN_PARAMETER = "chainParameter";
    public static final String JEFFREYS = "jeffreys";
    public static final String REVERSE = "reverse";
    public static final String SHAPE = "shape";


    public String getParserName() {
        return ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL;
    }

    /**
     * Reads a gamma distribution model from a DOM Document element.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject object = xo.getChild(CHAIN_PARAMETER);

        Parameter chainParameter = (Parameter) object.getChild(0);

        boolean jeffreys = xo.getAttribute(JEFFREYS, false);
        boolean reverse = xo.getAttribute(REVERSE, false);

        double shape = xo.getAttribute(SHAPE, 1.0);
        if (shape < 1.0) {
            throw new XMLParseException("ExponentialMarkovModel: shape parameter must be >= 1.0");
        }

        if (shape == 1.0) {
            System.out.println("Exponential markov model on parameter " +
                    chainParameter.getParameterName() + " (jeffreys=" + jeffreys + ", reverse=" +
                    reverse + ")");
        } else {
            System.out.println("Gamma markov model on parameter " +
                    chainParameter.getParameterName() + " (jeffreys=" + jeffreys + ", reverse=" +
                    reverse + " shape=" + shape + ")");
        }

        return new ExponentialMarkovModel(chainParameter, jeffreys, reverse, shape);
    }

    public String getParserDescription() {
        return "A continuous state, discrete time markov chain in which each new state is an " +
                "exponentially distributed variable with a mean of the previous state.";
    }

    public Class getReturnType() {
        return ExponentialMarkovModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(JEFFREYS, true),
            AttributeRule.newBooleanRule(REVERSE, true),
            new ElementRule(CHAIN_PARAMETER, Parameter.class)
    };
}