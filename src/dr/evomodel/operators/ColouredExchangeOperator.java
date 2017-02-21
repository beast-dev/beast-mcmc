/*
 * ColouredExchangeOperator.java
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Implements branch exchange operations.
 * There is a NARROW and WIDE variety.
 * The narrow exchange is very similar to a rooted-tree
 * nearest-neighbour interchange but with the restriction
 * that node height must remain consistent.
 * <p/>
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class ColouredExchangeOperator extends SimpleMCMCOperator {

    public static final String NARROW_EXCHANGE = "colouredNarrowExchange";
    public static final String WIDE_EXCHANGE = "colouredWideExchange";

    public static final int NARROW = 0;
    public static final int WIDE = 1;

    private static final int MAX_TRIES = 10000;

    private int mode = NARROW;
    private final TreeModel tree;

    private final ColourSamplerModel colouringModel;

    public ColouredExchangeOperator(int mode, TreeModel tree, ColourSamplerModel colouringModel, double weight) {
        this.mode = mode;
        this.tree = tree;

        this.colouringModel = colouringModel;

        setWeight(weight);
    }

    public double doOperation() {

        double logP = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

        int tipCount = tree.getExternalNodeCount();

        switch (mode) {
            case NARROW:
                narrow();
                break;
            case WIDE:
                wide();
                break;
        }

        if (tree.getExternalNodeCount() != tipCount) {
            throw new RuntimeException("Lost some tips in " + ((mode == NARROW) ? "NARROW mode." : "WIDE mode."));
        }

        colouringModel.resample();

        double logQ = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

        return logP - logQ;
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public void narrow() {

        NodeRef i = null, iP = null, j = null, jP = null;
        int tries = 0;

        //Echoose

        while (tries < MAX_TRIES) {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            while (tree.getRoot() == i || tree.getParent(i) == tree.getRoot()) {
                i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            }

            iP = tree.getParent(i);
            jP = tree.getParent(iP);
            j = tree.getChild(jP, 0);
            if (j == iP) {
                j = tree.getChild(jP, 1);
            }

            if ((tree.getNodeHeight(j) < tree.getNodeHeight(iP)) && (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
                break;
            }
            tries += 1;
        }
        //System.out.println("tries = " + tries);

        //Eupdate
        if (tries < MAX_TRIES) {
            eupdate(i, j, iP, jP);

            tree.pushTreeChangedEvent(iP);
            tree.pushTreeChangedEvent(jP);
        } else throw new RuntimeException("Couldn't find valid narrow move on this tree!!");
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public void wide() {

        NodeRef i = null, iP = null, j = null, jP = null;
        int tries = 0;

        //Echoose

        while (tries < MAX_TRIES) {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            while (tree.getRoot() == i) {
                i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            }

            j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            while (j == i || j == tree.getRoot()) {
                j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            }

            iP = tree.getParent(i);
            jP = tree.getParent(j);

            if ((iP != jP) && (i != jP) && (j != iP) &&
                    (tree.getNodeHeight(j) < tree.getNodeHeight(iP)) &&
                    (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
                break;
            }
            tries += 1;
        }
        //System.out.println("tries = " + tries);

        //Eupdate
        if (tries < MAX_TRIES) {
            eupdate(i, j, iP, jP);
        } else {
            throw new RuntimeException("Couldn't find valid wide move on this tree!");
        }
    }

    public int getMode() {
        return mode;
    }

    public String getOperatorName() {
        return ((mode == NARROW) ? "Narrow" : "Wide") + " Exchange";
    }

    public Element createOperatorElement(Document d) {
        throw new RuntimeException("not implemented");
    }

    private void eupdate(NodeRef i, NodeRef j, NodeRef iP, NodeRef jP)  {

        tree.beginTreeEdit();
        tree.removeChild(iP, i);
        tree.removeChild(jP, j);
        tree.addChild(jP, i);
        tree.addChild(iP, j);
        tree.endTreeEdit();
    }

    public double getMinimumAcceptanceLevel() {
        if (mode == NARROW) return 0.05;
        else return 0.01;
    }

    public double getMinimumGoodAcceptanceLevel() {
        if (mode == NARROW) return 0.05;
        else return 0.01;
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

    public static XMLObjectParser NARROW_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NARROW_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);
            double weight = xo.getDoubleAttribute("weight");

            return new ColouredExchangeOperator(NARROW, treeModel, colourSamplerModel, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a narrow exchange operator. " +
                    "This operator swaps a random subtree with its uncle.";
        }

        public Class getReturnType() {
            return ColouredExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(TreeModel.class),
                new ElementRule(ColourSamplerModel.class),
        };

    };

    public static XMLObjectParser WIDE_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return WIDE_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);
            double weight = xo.getDoubleAttribute("weight");

            return new ColouredExchangeOperator(WIDE, treeModel, colourSamplerModel, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a wide exchange operator. " +
                    "This operator swaps two random subtrees.";
        }

        public Class getReturnType() {
            return ColouredExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(ColourSamplerModel.class),
                new ElementRule(TreeModel.class)
        };

    };
}
