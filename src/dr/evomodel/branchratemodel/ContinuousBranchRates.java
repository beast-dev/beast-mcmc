/*
 * ContinuousBranchRates.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodelxml.branchratemodel.ContinuousBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.*;

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

    // The rate categories of each branch
    //d final TreeParameterModel rateCategories;
    final TreeParameterModel rateCategoryQuantiles;

    //private final int categoryCount;
    //private final double step;
    private final double[] rates;
    private boolean normalize = false;
    private double normalizeBranchRateTo = Double.NaN;
    private double scaleFactor = 1.0;
    private TreeModel treeModel;

    private Tree tree;

    //overSampling control the number of effective categories

    public ContinuousBranchRates(
            TreeModel tree,
            /*Parameter rateCategoryParameter, */
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel model
            /*int overSampling*/) {
        this(tree, /* rateCategoryParameter, */rateCategoryQuantilesParameter, model, /*overSampling, */false, Double.NaN);

    }

    public ContinuousBranchRates(
            TreeModel tree,
            /* Parameter rateCategoryParameter, */
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel model,
            /*int overSampling,*/
            boolean normalize,
            double normalizeBranchRateTo) {

        super(ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES);

        //d this.rateCategories = new TreeParameterModel(tree, rateCategoryParameter, false);
        this.rateCategoryQuantiles = new TreeParameterModel(tree, rateCategoryQuantilesParameter, false);

        rates = new double[tree.getNodeCount()];

        this.normalize = normalize;

        this.treeModel = tree;
        this.distributionModel = model;
        this.normalizeBranchRateTo = normalizeBranchRateTo;

        this.tree = new SimpleTree(tree);

        //Force the boundaries of rateCategoryParameter to match the category count
        //d Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        //d rateCategoryParameter.addBounds(bound);
        //rateCategoryQuantilesParameter.;



        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(1.0, 0.0, rateCategoryQuantilesParameter.getDimension());
        rateCategoryQuantilesParameter.addBounds(bound);
        /*if(rateCategoryQuantilesParameter.getBounds()==null) {
            System.out.println("oh NO!!! " + rateCategoryQuantilesParameter.getBounds().getLowerLimit(0) + "\t"
             + rateCategoryQuantilesParameter.getBounds().getUpperLimit(1));
        }*/


        /* for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            int index = (int) Math.floor((i + 0.5) * overSampling);
            rateCategoryParameter.setParameterValue(i, index);
        } */

        addModel(model);
        // AR - commented out: changes to the tree are handled by model changed events fired by rateCategories
//        addModel(tree);
        //d addModel(rateCategories);

        addModel(rateCategoryQuantiles);

        //addModel(treeModel); // Maybe
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events fired by rateCategories
//        addVariable(rateCategoryParameter);

        if (normalize) {
            tree.addModelListener(new ModelListener() {

                public void modelChangedEvent(Model model, Object object, int index) {
                    computeFactor();
                }

                public void modelRestored(Model model) {
                    computeFactor();
                }
            });
        }

        setupRates();
    }

    // compute scale factor
    private void computeFactor() {

        //scale mean rate to 1.0 or separate parameter

        double treeRate = 0.0;
        double treeTime = 0.0;

        //normalizeBranchRateTo = 1.0;
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
//d                int rateCategory = (int) Math.round(rateCategories.getNodeValue(treeModel, node));
//d                 treeRate += rates[rateCategory] * treeModel.getBranchLength(node);
                treeTime += treeModel.getBranchLength(node);

// d              System.out.println("rates and time\t" + rates[rateCategory] + "\t" + treeModel.getBranchLength(node));
            }
        }
        //treeRate /= treeTime;

        scaleFactor = normalizeBranchRateTo / (treeRate / treeTime);
        System.out.println("scaleFactor\t\t\t\t\t" + scaleFactor);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distributionModel) {
            setupRates();
            fireModelChanged();
        } //else if (model == rateCategories) {
        // AR - commented out: if just the rate categories have changed the rates will be the same
//            setupRates();
        //  fireModelChanged(null, index);
        //}
        else if (model == rateCategoryQuantiles) {
            setupRates();   // Maybe
            //rateCategories.fireModelChanged();
            fireModelChanged(null, index);
        } /*else if (model == treeModel) {
            setupRates(); // Maybe
        }*/
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events
        //setupRates();   // Maybe
    }

    protected void storeState() {
        //setupRates();   // Maybe
    }

    protected void restoreState() {
        setupRates();
    }

    protected void acceptState() {
        //setupRates();   // Maybe
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

        //int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));
        //System.out.println("dslkjafsdf " + node.getNumber());
        return rates[node.getNumber()] * scaleFactor;
    }

    /**
     * Calculates the actual rates corresponding to the category indices.
     */
    protected void setupRates() {

        //rateCategoryQuantiles.

        //double z = step / 2.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            //rates[i] = distributionModel.quantile(rateCategoryQuantiles.getNodeValue(rateCategoryQuantiles.getTreeModel(), rateCategoryQuantiles.getTreeModel().getNode(i) ));
            if (!tree.isRoot(tree.getNode(i))) {
                //System.out.println(rateCategoryQuantiles.getNodeValue(tree, tree.getNode(i)));
                rates[tree.getNode(i).getNumber()] = distributionModel.quantile(rateCategoryQuantiles.getNodeValue(tree, tree.getNode(i)));
            }

//            System.out.println("road " + i + "\t" + rateCategoryQuantiles.getNodeValue(tree, tree.getNode(i) ) + "\t" + tree.toString() + "\t" + tree.toString());
            //System.out.println("road " + i + "\t" + rateCategoryQuantiles.getNodeValue(rateCategoryQuantiles.getTreeModel(), rateCategoryQuantiles.getTreeModel().getNode(i) ) + "\t" + rateCategoryQuantiles.getTreeModel().toString() + "\t" + treeModel.toString());

            //System.out.print(rates[i]+"\t");
            //z += step;
            //rates[i] = distributionModel.quantile(z);
            //System.out.println("ape " + tree);
        }
        //System.out.println("\n");
        //System.out.println();
        if (normalize) computeFactor();
    }
}