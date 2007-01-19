package dr.evomodel.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.ARGModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/*
* ARGExchangeOperator.java
*
* Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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


/**
 * Implements branch exchange operations.
 * There is a NARROW and WIDE variety.
 * The narrow exchange is very similar to a rooted-tree
 * nearest-neighbour interchange but with the restriction
 * that node height must remain consistent.
 * <p/>
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */
public class ARGExchangeOperator extends SimpleMCMCOperator {

    public static final String NARROW_EXCHANGE = "argNarrowExchange";
    public static final String WIDE_EXCHANGE = "argWideExchange";

    public static final int NARROW = 0;
    public static final int WIDE = 1;

    private static final int MAX_TRIES = 10000;

    private int mode = NARROW;
    private ARGModel tree;

    public ARGExchangeOperator(int mode, ARGModel tree, int weight) {
        this.mode = mode;
        this.tree = tree;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {

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

        return 0.0;
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public void narrow() throws OperatorFailedException {

        NodeRef i = null, iP = null, j = null, jP = null;
        int tries = 0;

        //Echoose

        while (tries < MAX_TRIES) {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            while (tree.getRoot() == i || tree.getParent(i, 0) == tree.getRoot() || tree.getParent(i, 1) == tree.getRoot()) {
                i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            }

            iP = tree.getParent(i, 0);
            if (tree.isReassortment(i) && MathUtils.nextBoolean())
                iP = tree.getParent(i, 1);

            jP = tree.getParent(iP, 0);
            if (tree.isReassortment(iP) && MathUtils.nextBoolean())
                jP = tree.getParent(iP, 1);

            j = tree.getChild(jP, 0);
            if (j == iP) {
                j = tree.getChild(jP, 1);
            }

            if (j != iP && i != j &&             // can still occur if i is child of doubly-linked reassortment
                    (tree.getNodeHeight(j) < tree.getNodeHeight(iP)) && (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
                break;
            }
            tries += 1;
        }
        //System.out.println("tries = " + tries);

        //Eupdate
        if (tries < MAX_TRIES) {
            eupdateARG(i, j, iP, jP);

            tree.pushTreeChangedEvent(iP);
            tree.pushTreeChangedEvent(jP);
        } else throw new OperatorFailedException("Couldn't find valid narrow move on this tree!!");
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public void wide() throws OperatorFailedException {

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
            eupdateARG(i, j, iP, jP);
        } else throw new OperatorFailedException("Couldn't find valid wide move on this tree!");
    }

    public int getMode() {
        return mode;
    }

    public String getOperatorName() {
        return ((mode == NARROW) ? "Narrow" : "Wide") + " Exchange";
    }

    private void eupdateARG(NodeRef i, NodeRef j, NodeRef iP, NodeRef jP) throws OperatorFailedException {

        // There are three different cases:
        // 1) neither i nor j are reassortments, 2) either i or j are reassortments, 3) both i and j are reassortments

        tree.beginTreeEdit();

        boolean iBifurcation = tree.isBifurcation(i);
        boolean jBifurcation = tree.isBifurcation(j);

        if (iBifurcation && jBifurcation) {
            tree.removeChild(iP, i);
            tree.removeChild(jP, j);
            tree.addChild(jP, i);
            tree.addChild(iP, j);
        } else if (!iBifurcation && !jBifurcation) {
            tree.singleRemoveChild(iP, i);
            tree.singleRemoveChild(jP, j);
            tree.singleAddChild(jP, i);
            tree.singleAddChild(iP, j);
        } else {
            if (jBifurcation) {
                // one reassortment; force i to be bifurcation and j to be reassortment
                NodeRef t = i;
                NodeRef tP = iP;
                i = j;
                iP = jP;
                j = t;
                jP = tP;
            }
//            System.err.println(tree.toGraphString());
//            ARGModel.Node iNode = (ARGModel.Node) i;
//            ARGModel.Node jNode = (ARGModel.Node) j;
//            System.err.println("i = "+iNode.number+" : j = "+jNode.number);
            tree.removeChild(iP, i);
            tree.singleRemoveChild(jP, j);
            tree.addChild(jP, i);
            tree.singleAddChild(iP, j);
        }

        try {
            tree.endTreeEdit();
        } catch (MutableTree.InvalidTreeException ite) {
            throw new OperatorFailedException(ite.toString());
        }
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
        if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
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

            ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
            int weight = xo.getIntegerAttribute("weight");

            return new ARGExchangeOperator(NARROW, treeModel, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a narrow exchange operator. " +
                    "This operator swaps a random subtree with its uncle.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule("weight"),
                new ElementRule(ARGModel.class)
        };

    };

    public static XMLObjectParser WIDE_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return WIDE_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
            int weight = xo.getIntegerAttribute("weight");

            return new ARGExchangeOperator(WIDE, treeModel, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a wide exchange operator. " +
                    "This operator swaps two random subtrees.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule("weight"),
                new ElementRule(ARGModel.class)
        };

    };
}