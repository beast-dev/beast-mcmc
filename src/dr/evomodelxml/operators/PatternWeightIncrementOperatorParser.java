/*
 * PatternWeightIncrementOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.PatternWeightIncrementOperator;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class PatternWeightIncrementOperatorParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return PatternWeightIncrementOperator.PATTERN_WEIGHT_INCREMENT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

       BeagleTreeLikelihood treeLikelihood = (BeagleTreeLikelihood) xo.getChild(BeagleTreeLikelihood.class);
        final double weight = xo.getDoubleAttribute("weight");
        return new PatternWeightIncrementOperator(treeLikelihood, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element progressively adds sites into the pattern list.";
    }

    public Class getReturnType() {
        return PatternWeightIncrementOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule("weight"),
            new ElementRule(BeagleTreeLikelihood.class),
    };
}
