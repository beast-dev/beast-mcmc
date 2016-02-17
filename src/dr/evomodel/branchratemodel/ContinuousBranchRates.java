/*
 * ContinuousBranchRates.java
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
import dr.evomodelxml.branchratemodel.ContinuousBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.*;
import dr.math.MathUtils;

/**
 *
 * Uses the implementation of continuous branch rates described in:
 * Li, W.L.S. and A.J. Drummond, Model Averaging and Bayes Factor
 * Calculation of Relaxed Molecular Clocks in Bayesian Phylogenetics.
 * Molecular Biology and Evolution, 2012. 29(2): p. 751-761.
 *
 * @author Wai Lok Sibon Li
 * @version $Id: ContinuousBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class ContinuousBranchRates extends AbstractBranchRateModel {

    private final ParametricDistributionModel distributionModel;

    // The rate quantiles of each branch
    final TreeParameterModel rateCategoryQuantiles;

    private final double[] rates;
    private boolean normalize = false;
    private double normalizeBranchRateTo = Double.NaN;
    private double scaleFactor = 1.0;
    private TreeModel treeModel;

    private boolean updateScaleFactor = false;
    private boolean updateRates = true;

    public ContinuousBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel model) {
        this(tree, rateCategoryQuantilesParameter, model, false, Double.NaN);

    }

    public ContinuousBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel model,
            boolean normalize,
            double normalizeBranchRateTo) {

        super(ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES);

        this.rateCategoryQuantiles = new TreeParameterModel(tree, rateCategoryQuantilesParameter, false);

        rates = new double[tree.getNodeCount()];

        this.normalize = normalize;

        this.treeModel = tree;
        this.distributionModel = model;
        this.normalizeBranchRateTo = normalizeBranchRateTo;

        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(1.0, 0.0, rateCategoryQuantilesParameter.getDimension());
        rateCategoryQuantilesParameter.addBounds(bound);

        randomizeRates();

        addModel(model);
        addModel(rateCategoryQuantiles);

        if (normalize) {
            // if we want to normalize the rates then we need to listen for changes on the tree
            addModel(treeModel);
            updateScaleFactor = true;
        }

        updateRates = true;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distributionModel) {
            updateRates = true;
            fireModelChanged();
        } else if (model == rateCategoryQuantiles) {
            updateRates = true;
            fireModelChanged(null, index);
        } else if (model == treeModel && normalize) {
            updateScaleFactor = true;
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
        updateRates = true;
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

        if (updateRates) {
            computeRates();
        }

        if (updateScaleFactor) {
            computeFactor();
        }

        return rates[node.getNumber()] * scaleFactor;
    }

    /**
     * Calculates the actual rates corresponding to the quantiles.
     */
    private void randomizeRates() {

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if (!treeModel.isRoot(treeModel.getNode(i))) {
                double r = MathUtils.nextDouble();
                rateCategoryQuantiles.setNodeValue(treeModel, treeModel.getNode(i), r);
            }
        }

        updateRates = false;
    }

    /**
     * Calculates the actual rates corresponding to the quantiles.
     */
    private void computeRates() {

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if (!treeModel.isRoot(treeModel.getNode(i))) {
                rates[treeModel.getNode(i).getNumber()] = distributionModel.quantile(rateCategoryQuantiles.getNodeValue(treeModel, treeModel.getNode(i)));
            }
        }

        updateRates = false;
    }

    // compute scale factor
    private void computeFactor() {

        //scale mean rate to 1.0 or separate parameter

        double treeRate = 0.0;
        double treeTime = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
                treeTime += treeModel.getBranchLength(node);
            }
        }

        scaleFactor = normalizeBranchRateTo / (treeRate / treeTime);

        updateScaleFactor = false;
    }

}