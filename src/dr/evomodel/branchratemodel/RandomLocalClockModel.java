/*
 * RandomLocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * A model of rates evolving along a tree, such that at each node the rate may change or not depending on an
 * indicator which chooses whether the parent rate is inherited or a new rate begins on the branch above the node.
 * This model is implemented using stochastic variable selection.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class RandomLocalClockModel extends AbstractBranchRateModel
        implements RandomLocalTreeVariable, Citable {

    public RandomLocalClockModel(TreeModel treeModel,
                                 Parameter meanRateParameter,
                                 Parameter rateIndicatorParameter,
                                 Parameter ratesParameter,
                                 boolean ratesAreMultipliers,
                                 double threshold) {

        super(RandomLocalClockModelParser.LOCAL_BRANCH_RATES);

        this.ratesAreMultipliers = ratesAreMultipliers;

        indicators = new TreeParameterModel(treeModel, rateIndicatorParameter, false);
        rates = new TreeParameterModel(treeModel, ratesParameter, false);

        if (Double.isNaN(threshold)) { // NaN == binary values
            rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(1, 0, rateIndicatorParameter.getDimension()));
            this.threshold = 0.5;
            for (int i = 0; i < rateIndicatorParameter.getDimension(); ++i) {
                rateIndicatorParameter.setParameterValue(i, 0.0);
            }
        } else {
            rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, -Double.MAX_VALUE, rateIndicatorParameter.getDimension()));
            this.threshold = threshold;
        }
        ratesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0, ratesParameter.getDimension()));

        for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
            ratesParameter.setParameterValue(i, 1.0);
        }

        this.meanRateParameter = meanRateParameter;

        addModel(treeModel);
        this.treeModel = treeModel;

        addModel(indicators);
        addModel(rates);
        if (meanRateParameter != null) addVariable(meanRateParameter);

        unscaledBranchRates = new double[treeModel.getNodeCount()];

        Logger.getLogger("dr.evomodel").info("  indicator parameter name is '" + rateIndicatorParameter.getId() + "' with threshold = " + threshold);

        recalculateScaleFactor();
    }

    /**
     * @param tree the tree
     * @param node the node to retrieve the variable of
     * @return the raw real-valued variable at this node
     */
    public final double getVariable(Tree tree, NodeRef node) {
        return rates.getNodeValue(tree, node);
    }

    /**
     * @param tree the tree
     * @param node the node
     * @return true of the variable at this node is included in function, thus representing a change in the
     *         function looking down the tree.
     */
    public final boolean isVariableSelected(Tree tree, NodeRef node) {
        return indicators.getNodeValue(tree, node) > threshold;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        recalculationNeeded = true;
//        recalculateScaleFactor();  // Careful about change here; seems fine ... but who knows without unit-tests
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        if (recalculationNeeded) {
            recalculateScaleFactor();
            recalculationNeeded = false;
        }
        return unscaledBranchRates[node.getNumber()] * scaleFactor;
    }

    private void calculateUnscaledBranchRates(TreeModel tree) {
        cubr(tree, tree.getRoot(), 1.0);
    }

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch rates across the tree
     * taking into account the indicator variables.
     *
     * @param tree the tree
     * @param node the node
     * @param rate the rate of the parent node
     */
    private void cubr(TreeModel tree, NodeRef node, double rate) {

        int nodeNumber = node.getNumber();

        if (!tree.isRoot(node)) {
            if (isVariableSelected(tree, node)) {
                if (ratesAreMultipliers) {
                    rate *= getVariable(tree, node);
                } else {
                    rate = getVariable(tree, node);
                }
            }
        }
        unscaledBranchRates[nodeNumber] = rate;

        int childCount = tree.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            cubr(tree, tree.getChild(node, i), rate);
        }
    }

    private void recalculateScaleFactor() {

        calculateUnscaledBranchRates(treeModel);

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchInTime =
                        treeModel.getNodeHeight(treeModel.getParent(node)) -
                                treeModel.getNodeHeight(node);

                double branchLength = branchInTime * unscaledBranchRates[node.getNumber()];

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
        }

        scaleFactor = timeTotal / branchTotal;

        if (meanRateParameter != null)
            scaleFactor *= meanRateParameter.getParameterValue(0);
    }

    // AR - as TreeParameterModels are now loggable, the indicator parameter should be logged
    // directly.
//    private static String[] attributeLabel = {"changed"};
//
//    public String[] getNodeAttributeLabel() {
//        return attributeLabel;
//    }
//
//    public String[] getAttributeForNode(Tree tree, NodeRef node) {
//
//        if (tree.isRoot(node)) {
//            return new String[]{"false"};
//        }
//
//        return new String[]{(isVariableSelected((TreeModel) tree, node) ? "true" : "false")};
//    }

    // the scale factor necessary to maintain the mean rate
    private double scaleFactor;

    // the tree model
    private TreeModel treeModel;

    // true if the rate variables are treated as relative
    // to the parent rate rather than absolute rates
    private boolean ratesAreMultipliers = false;

    // the unscaled rates of each branch, taking into account the indicators
    private double[] unscaledBranchRates;

    // the mean rate across all the tree, if null then mean rate is scaled to 1.0
    private Parameter meanRateParameter;

    private TreeParameterModel indicators;
    private TreeParameterModel rates;

    boolean recalculationNeeded = true;

    private final double threshold;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        return "Local clock model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("AJ", "Drummond"),
                    new Author("MA", "Suchard")
            },
            "Bayesian random local clocks, or one rate to rule them all",
            2010,
            "BMC Biology",
            "8: 114",
            "10.1186/1741-7007-8-114"
    );
}