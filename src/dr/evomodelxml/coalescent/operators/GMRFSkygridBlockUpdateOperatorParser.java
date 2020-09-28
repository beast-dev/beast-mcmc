/*
 * GMRFSkyrideBlockUpdateOperatorParser.java
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

package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.GMRFSkygridLikelihood;
import dr.evomodel.coalescent.operators.GMRFSkygridBlockUpdateOperator;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.IOException;
import java.util.logging.*;

/**
 *
 */
public class GMRFSkygridBlockUpdateOperatorParser extends AbstractXMLObjectParser {

    public static final String BLOCK_UPDATE_OPERATOR = "gmrfSkygridBlockUpdateOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String MAX_ITERATIONS = "maxIterations";
    public static final String STOP_VALUE = "stopValue";

    public String getParserName() {
        return BLOCK_UPDATE_OPERATOR;
    }

    public String[] getParserNames() {
        return new String[]{
                getParserName(),
                BLOCK_UPDATE_OPERATOR,
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);
        if (mode == AdaptationMode.DEFAULT) mode = AdaptationMode.ADAPTATION_ON;

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
        if (scaleFactor == 1.0) {
            mode = AdaptationMode.ADAPTATION_OFF;
        }

        if (scaleFactor < 1.0) {
            throw new XMLParseException("scaleFactor must be greater than or equal to 1.0");
        }

        int maxIterations = xo.getAttribute(MAX_ITERATIONS, 200);

        double stopValue = xo.getAttribute(STOP_VALUE, 0.01);
        
        GMRFSkygridLikelihood gmrfSkygridLikelihood = (GMRFSkygridLikelihood) xo.getChild(GMRFSkygridLikelihood.class);
        return new GMRFSkygridBlockUpdateOperator(gmrfSkygridLikelihood, weight, mode, scaleFactor,
                maxIterations, stopValue);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a GMRF block-update operator for the joint distribution of the population sizes and precision parameter.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newDoubleRule(STOP_VALUE, true),
            AttributeRule.newIntegerRule(MAX_ITERATIONS, true),
            new ElementRule(GMRFSkygridLikelihood.class)
    };

}
