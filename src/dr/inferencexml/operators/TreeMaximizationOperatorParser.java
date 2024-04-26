/*
 * MaximizerWrtParameterOperatorParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.TreeMaximizationOperator;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class TreeMaximizationOperatorParser extends AbstractXMLObjectParser {

    public final static String NAME = "treeMaximizationOperator";
    private final String MAX_STEP = "maxSteps";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int maxSteps = (int) xo.getAttribute(MAX_STEP, 20);
        MaximizerWrtParameter maximizerWrtParameter = (MaximizerWrtParameter) xo.getChild(MaximizerWrtParameter.class);
        SimpleMCMCOperator treeOperator = (SimpleMCMCOperator) xo.getChild(MCMCOperator.class);
        return new TreeMaximizationOperator(treeOperator, maximizerWrtParameter, maxSteps, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(MaximizerWrtParameter.class),
            AttributeRule.newIntegerRule(MAX_STEP, true),
            new ElementRule(SimpleMCMCOperator.class, "Tree operator"),

    };
    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TreeMaximizationOperator.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
