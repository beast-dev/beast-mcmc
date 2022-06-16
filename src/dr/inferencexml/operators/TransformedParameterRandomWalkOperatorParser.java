/*
 * TransformedParameterRandomWalkOperatorParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.BoundedSpace;
import dr.inference.model.TransformedParameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RandomWalkOperator;
import dr.inference.operators.TransformedParameterRandomWalkOperator;
import dr.xml.*;

public class TransformedParameterRandomWalkOperatorParser extends RandomWalkOperatorParser {

    public static final String TRANSFORMED_PARAMETER_RANDOM_WALK_OPERATOR = "transformedParameterRandomWalkOperator";

    public String getParserName() {
        return TRANSFORMED_PARAMETER_RANDOM_WALK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Object randomWalk;
        try {
            randomWalk = super.parseXMLObject(xo);
        } catch (XMLParseException e) {
            throw new XMLParseException("RandomWalkOperatorParser failled in TraansformedParameterRandomWalkOperator.");
        }
        BoundedSpace bounds = (BoundedSpace) xo.getChild(BoundedSpace.class);
        return new TransformedParameterRandomWalkOperator((RandomWalkOperator) randomWalk, bounds);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a random walk operator on a given transformed parameter.";
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
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(UPDATE_INDEX,
                    new XMLSyntaxRule[]{
                            new ElementRule(TransformedParameter.class),
                    }, true),
            new StringAttributeRule(BOUNDARY_CONDITION, null, RandomWalkOperator.BoundaryCondition.values(), true),
            new ElementRule(TransformedParameter.class),
            new ElementRule(BoundedSpace.class, true)
    };
}
