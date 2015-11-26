/*
 * TipStateSwapOperator.java
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

import dr.evomodel.treelikelihood.AbstractLikelihoodCore;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TipStateSwapOperator extends SimpleMCMCOperator {

    public static final String TIP_STATE_OPERATOR = "tipStateSwapOperator";

    public TipStateSwapOperator(TreeLikelihood treeLikelihood, double weight) {
        this.treeLikelihood = treeLikelihood;
        setWeight(weight);
        int patternCount = treeLikelihood.getPatternCount();
        states1 = new int[patternCount];
        states2 = new int[patternCount];
    }

    private final int[] states1;
    private final int[] states2;
    private int index1;
    private int index2;

    public double doOperation() throws OperatorFailedException {

        int tipCount = treeLikelihood.getTreeModel().getExternalNodeCount();

        // Choose two tips to swap
        index1 = MathUtils.nextInt(tipCount);
        index2 = index1;
        while (index2 == index1)
            index2 = MathUtils.nextInt(tipCount);

        swap(index1, index2);

        treeLikelihood.makeDirty();

        return 0;
    }

    private void swap(int i, int j) {
        AbstractLikelihoodCore likelihoodCore = (AbstractLikelihoodCore) treeLikelihood.getLikelihoodCore();

        likelihoodCore.getNodeStates(i, states1);
        likelihoodCore.getNodeStates(j, states2);

        likelihoodCore.setNodeStates(j, states1);
        likelihoodCore.setNodeStates(i, states2);
    }

    public void reject() {
        super.reject();
        // There is currently no restore functions for tip states, so manually adjust state
        swap(index1, index2);
    }

    public String getPerformanceSuggestion() {
        if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String getOperatorName() {
        return TIP_STATE_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TIP_STATE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeLikelihood treeLikelihood =
                    (TreeLikelihood) xo.getChild(TreeLikelihood.class);
            final double weight = xo.getDoubleAttribute("weight");
            return new TipStateSwapOperator(treeLikelihood, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an operator to swap tip states between two random tips.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(TreeLikelihood.class),
        };

    };


    private final TreeLikelihood treeLikelihood;
}
