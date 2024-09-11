/*
 * BitSwapOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.operators.BitSwapOperator;
import dr.inference.operators.MCMCOperator;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.xml.*;

/**
 *
 */
public class BitSwapOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_SWAP_OPERATOR = "bitSwapOperator";
    public static final String RADIUS = "radius";

    private static final String DATA = MixedDistributionLikelihoodParser.DATA;
    private static final String INDICATORS = MixedDistributionLikelihoodParser.INDICATORS;

    public String getParserName() {
        return BIT_SWAP_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter data = (Parameter) ((XMLObject) xo.getChild(DATA)).getChild(Parameter.class);
        Parameter indicators = (Parameter) ((XMLObject) xo.getChild(INDICATORS)).getChild(Parameter.class);
        int radius = -1;

        if (xo.hasAttribute(RADIUS)) {
            double rd = xo.getDoubleAttribute(RADIUS);

            if (rd > 0) {
                if (rd < 1) {
                    rd = Math.round(rd * indicators.getDimension());
                }
                radius = (int) Math.round(rd);
                if (!(radius >= 1 && radius < indicators.getDimension() - 1)) {
                    radius = -1;
                }
            }
            if (radius < 1) {
                throw new XMLParseException("invalid radius " + rd);
            }
        }

        return new BitSwapOperator(data, indicators, radius, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-swap operator on a given parameter and data.";
    }

    public Class getReturnType() {
        return BitSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(RADIUS, true),
            new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
            new ElementRule(INDICATORS, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
    };
}
