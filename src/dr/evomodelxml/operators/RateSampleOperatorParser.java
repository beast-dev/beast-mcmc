/*
 * RateSampleOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.oldevomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.operators.RateSampleOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateSampleOperatorParser extends AbstractXMLObjectParser {

    public static final String SAMPLE_OPERATOR = "rateSampleOperator";
    public static final String SAMPLE_ALL = "sampleAll";

    public String getParserName() {
        return SAMPLE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final boolean sampleAll = xo.getBooleanAttribute(SAMPLE_ALL);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        RateEvolutionLikelihood rateEvolution = (RateEvolutionLikelihood) xo.getChild(RateEvolutionLikelihood.class);

        RateSampleOperator operator = new RateSampleOperator(treeModel, sampleAll, rateEvolution);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a rateSample operator on a given parameter.";
    }

    public Class getReturnType() {
        return RateSampleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(SAMPLE_ALL, true),
            new ElementRule(TreeModel.class),
            new ElementRule(RateEvolutionLikelihood.class, true),
    };
}
