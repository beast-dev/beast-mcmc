/*
 * DiscretizedBranchRates.java
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
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodelxml.DiscretizedBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Michael Defoin Platel
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class DiscretizedBranchRates extends AbstractModel implements BranchRateModel {

    private final ParametricDistributionModel distributionModel;

    // The rate categories of each branch
    final TreeParameterModel rateCategories;

    private final int categoryCount;
    private final double step;
    private final double[] rates;

    //overSampling control the number of effective categories

    public DiscretizedBranchRates(
            TreeModel tree,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            int overSampling) {

        super(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

        this.rateCategories = new TreeParameterModel(tree, rateCategoryParameter, false);

        categoryCount = (tree.getNodeCount() - 1) * overSampling;
        step = 1.0 / (double) categoryCount;

        rates = new double[categoryCount];

        this.distributionModel = model;

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        rateCategoryParameter.addBounds(bound);

        for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            int index = (int) Math.floor((i + 0.5) * overSampling);
            rateCategoryParameter.setParameterValue(i, index);
        }

        addModel(model);
        // AR - commented out: changes to the tree are handled by model changed events fired by rateCategories
//        addModel(tree);
        addModel(rateCategories);
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events fired by rateCategories
//        addVariable(rateCategoryParameter);

        setupRates();
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distributionModel) {
            setupRates();
            fireModelChanged();
        } else if (model == rateCategories) {
            // AR - commented out: if just the rate categories have changed the rates will be the same
//            setupRates();
            fireModelChanged(null, index);
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events
//        setupRates();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        setupRates();
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

        int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));

        return rates[rateCategory];
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    /**
     * Calculates the actual rates corresponding to the category indices.
     */
    protected void setupRates() {

        double z = step / 2.0;
        for (int i = 0; i < categoryCount; i++) {
            rates[i] = distributionModel.quantile(z);
            z += step;
        }
    }
}
