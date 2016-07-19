/*
 * UCLikelihoodParser.java
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

package dr.oldevomodelxml.clock;

import dr.oldevomodel.clock.RateEvolutionLikelihood;
import dr.oldevomodel.clock.UCLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class UCLikelihoodParser extends AbstractXMLObjectParser {

    public static final String UC_LIKELIHOOD = "UCLikelihood";

    public static final String VARIANCE = "variance";

    public String getParserName() {
        return UC_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.RATES);

        Parameter rootRate = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.ROOTRATE);

        Parameter variance = (Parameter) xo.getElementFirstChild(VARIANCE);


        boolean isLogSpace = xo.getAttribute(RateEvolutionLikelihood.LOGSPACE, false);

        return new UCLikelihood(tree, ratesParameter, variance, rootRate, isLogSpace);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an object that can calculate the likelihood " +
                        "of rates in a tree under the assumption of " +
                        "(log)normally distributed rates. ";
    }

    public Class getReturnType() {
        return UCLikelihood.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RateEvolutionLikelihood.RATES, Parameter.class, "The branch rates parameter", false),
            AttributeRule.newBooleanRule(RateEvolutionLikelihood.LOGSPACE, true, "true if model considers the log of the rates."),
            new ElementRule(RateEvolutionLikelihood.ROOTRATE, Parameter.class, "The root rate parameter"),
            new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
    };

}
