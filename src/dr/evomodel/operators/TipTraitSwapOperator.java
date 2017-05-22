/*
 * TipTraitSwapOperator.java
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

package dr.evomodel.operators;

import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TipTraitSwapOperator extends SimpleMCMCOperator {

    public static final String TIP_TRAIT_SWAP_OPERATOR = "tipTraitSwapOperator";

    public TipTraitSwapOperator(String traitName, IntegratedMultivariateTraitLikelihood traitLikelihood, double weight) {
        this.traitLikelihood = traitLikelihood;
        this.traitName = traitName;
        setWeight(weight);
    }

    private int index1;
    private int index2;

    public double doOperation() {

        int tipCount = traitLikelihood.getTreeModel().getExternalNodeCount();

        // Choose two tips to swap
        index1 = MathUtils.nextInt(tipCount);
        index2 = index1;
        while (index2 == index1)
            index2 = MathUtils.nextInt(tipCount);

        swap(index1, index2);

        traitLikelihood.makeDirty();

        return 0;
    }

    private void swap(int i, int j) {
        double[] trait1 = traitLikelihood.getTipDataValues(i);
        double[] trait2 = traitLikelihood.getTipDataValues(j);

        traitLikelihood.setTipDataValuesForNode(j, trait1);
        traitLikelihood.setTipDataValuesForNode(i, trait2);           
    }

    public void reject() {
        super.reject();
        // There is currently no restore functions for tip states, so manually adjust state
        swap(index1, index2);
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String getOperatorName() {
        return TIP_TRAIT_SWAP_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TIP_TRAIT_SWAP_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            IntegratedMultivariateTraitLikelihood traitLikelihood =
                    (IntegratedMultivariateTraitLikelihood) xo.getChild(IntegratedMultivariateTraitLikelihood.class);
            final String traitName = xo.getStringAttribute("trait");
            final double weight = xo.getDoubleAttribute("weight");
            return new TipTraitSwapOperator(traitName, traitLikelihood, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an operator to swap tip traits between two random tips.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("trait"),
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(IntegratedMultivariateTraitLikelihood.class),
        };

    };


    private final IntegratedMultivariateTraitLikelihood traitLikelihood;
    private final String traitName;
}
