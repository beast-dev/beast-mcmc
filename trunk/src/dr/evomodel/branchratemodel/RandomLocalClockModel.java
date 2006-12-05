/*
 * RandomLocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.NodeAttributeProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class RandomLocalClockModel extends AbstractModel implements BranchRateModel, NodeAttributeProvider {

    public static final String LOCAL_BRANCH_RATES = "randomLocalClockModel";
    public static final String RATE_INDICATORS = "rateIndicator";
    public static final String RATES = "rates";

    // The rate categories of each branch

    // the index of the root node.
    private int rootNodeNumber;
    private int storedRootNodeNumber;

    double scaleFactor;
    TreeModel treeModel;

    public RandomLocalClockModel(TreeModel treeModel, Parameter rateIndicatorParameter, Parameter ratesParameter) {

        super(LOCAL_BRANCH_RATES);

        if (rateIndicatorParameter.getDimension() != treeModel.getNodeCount() - 1) {
            throw new IllegalArgumentException("The rate category parameter must be of length nodeCount-1");
        }

        for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
            rateIndicatorParameter.setParameterValue(i, 0.0);
            ratesParameter.setParameterValue(i, 1.0);
        }

        addModel(treeModel);
        this.treeModel = treeModel;

        addParameter(rateIndicatorParameter);
        addParameter(ratesParameter);

        rootNodeNumber = treeModel.getRoot().getNumber();
        storedRootNodeNumber = rootNodeNumber;

        recalculateScaleFactor();
    }

    private void recalculateScaleFactor() {

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchInTime =
                        treeModel.getNodeHeight(treeModel.getParent(node)) -
                                treeModel.getNodeHeight(node);

                double branchLength = branchInTime * getUnscaledBranchRate(treeModel, node);

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
        }

        scaleFactor = timeTotal / branchTotal;
    }


    public void handleModelChangedEvent(Model model, Object object, int index) {
        recalculateScaleFactor();
        fireModelChanged();
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        recalculateScaleFactor();
        fireModelChanged();
    }

    protected void storeState() {
        storedRootNodeNumber = rootNodeNumber;
    }

    protected void restoreState() {
        recalculateScaleFactor();
        rootNodeNumber = storedRootNodeNumber;
    }

    protected void acceptState() {
    }

    public double getBranchRate(Tree tree, NodeRef node) {

        return getUnscaledBranchRate(tree, node) * scaleFactor;
    }

    private double getUnscaledBranchRate(Tree tree, NodeRef node) {

        if (tree.isRoot(node)) {
            return 1.0;
        } else {

            double rate;
            if (isRateChangeOnBranchAbove(tree, node)) {
                rate = tree.getNodeRate(node);
            } else {
                rate = getUnscaledBranchRate(tree, tree.getParent(node));
            }
            return rate;
        }
    }

    public final boolean isRateChangeOnBranchAbove(Tree tree, NodeRef node) {
        return (int) Math.round(((TreeModel) tree).getNodeTrait(node)) == 1;
    }

    public String getNodeAttributeLabel() {
        return "changed";
    }

    public String getAttributeForNode(Tree tree, NodeRef node) {
        return (isRateChangeOnBranchAbove(tree, node) ? "1" : "0");
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOCAL_BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter rateIndicatorParameter = (Parameter) xo.getSocketChild(RATE_INDICATORS);
            Parameter ratesParameter = (Parameter) xo.getSocketChild(RATES);

            Logger.getLogger("dr.evomodel").info("Using random local clock (RLM) model.");

            return new RandomLocalClockModel(tree, rateIndicatorParameter, ratesParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an random local clock (RLM) model." +
                            "Each branch either has a new independent rate or " +
                            "inherits the rate of the branch above it depending on the indicator vector.";
        }

        public Class getReturnType() {
            return RandomLocalClockModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RATE_INDICATORS, Parameter.class, "The rate change indicators parameter", false),
                new ElementRule(RATES, Parameter.class, "The rate changes parameter", false)
        };
    };


}