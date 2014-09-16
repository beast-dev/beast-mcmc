/*
 * DiscretizedBranchRates.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Michael Defoin Platel
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class DiscretizedBranchRates extends AbstractBranchRateModel {
    // Turn on an off the caching on rates for categories -
    // if off then the rates will be flagged to update on
    // a restore.
    // Currently turned off as it is not working with multiple partitions for
    // some reason.
    private static final boolean CACHE_RATES = false;

    private final ParametricDistributionModel distributionModel;

    // The rate categories of each branch
    final TreeParameterModel rateCategories;

    private final int categoryCount;
    private final double step;
    private final double[][] rates;
    private final boolean normalize;
    private final double normalizeBranchRateTo;

    private final TreeModel treeModel;
    private final double logDensityNormalizationConstant;

    private double scaleFactor = 1.0;
    private double storedScaleFactor;

    private boolean updateRateCategories = true;
    private int currentRateArrayIndex = 0;
    private int storedRateArrayIndex;

    //overSampling control the number of effective categories

    public DiscretizedBranchRates(
            TreeModel tree,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            int overSampling) {
        this(tree, rateCategoryParameter, model, overSampling, false, Double.NaN, false, false);

    }

    public DiscretizedBranchRates(
            TreeModel tree,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            int overSampling,
            boolean normalize,
            double normalizeBranchRateTo,
            boolean randomizeRates,
            boolean keepRates) {

        super(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

        this.rateCategories = new TreeParameterModel(tree, rateCategoryParameter, false);

        categoryCount = (tree.getNodeCount() - 1) * overSampling;
        step = 1.0 / (double) categoryCount;

        rates = new double[2][categoryCount];

        this.normalize = normalize;

        this.treeModel = tree;
        this.distributionModel = model;
        this.normalizeBranchRateTo = normalizeBranchRateTo;

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        rateCategoryParameter.addBounds(bound);

        for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            if (!keepRates) {
                int index = (randomizeRates) ?
                        MathUtils.nextInt(rateCategoryParameter.getDimension() * overSampling) : // random rate
                        (int) Math.floor((i + 0.5) * overSampling); // default behavior
                rateCategoryParameter.setParameterValue(i, index);
            }
        }

        addModel(model);
        addModel(rateCategories);

        updateRateCategories = true;

        // Each parameter take any value in [1, \ldots, categoryCount]
        // NB But this depends on the transition kernel employed.  Using swap-only results in a different constant
        logDensityNormalizationConstant = -rateCategoryParameter.getDimension() * Math.log(categoryCount);
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
                int rateCategory = (int) Math.round(rateCategories.getNodeValue(treeModel, node));
                treeRate += rates[currentRateArrayIndex][rateCategory] * treeModel.getBranchLength(node);
                treeTime += treeModel.getBranchLength(node);

                //System.out.println("rates and time\t" + rates[rateCategory] + "\t" + treeModel.getBranchLength(node));
            }
        }
        //treeRate /= treeTime;

        scaleFactor = normalizeBranchRateTo / (treeRate / treeTime);
        //System.out.println("scaleFactor\t\t\t\t\t" + scaleFactor);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distributionModel) {
            updateRateCategories = true;
            fireModelChanged();
        } else if (model == rateCategories) {
            fireModelChanged(null, index);
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do here
   }

    protected void storeState() {
        if (CACHE_RATES) {
            storedRateArrayIndex = currentRateArrayIndex;
            storedScaleFactor = scaleFactor;
        }
    }

    protected void restoreState() {
        if (CACHE_RATES) {
            currentRateArrayIndex = storedRateArrayIndex;
            scaleFactor = storedScaleFactor;
        } else {
            updateRateCategories = true;
        }
    }

    protected void acceptState() {
    }

    public final double getBranchRate(final Tree tree, final NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

        if (updateRateCategories) {
            setupRates();
        }

        int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));

        //System.out.println(rates[rateCategory] + "\t"  + rateCategory);
        return rates[currentRateArrayIndex][rateCategory] * scaleFactor;
    }

    /**
     * Calculates the actual rates corresponding to the category indices.
     */
    private void setupRates() {

        if (CACHE_RATES) {
            // flip the current array index
            currentRateArrayIndex = 1 - currentRateArrayIndex;
        }

        double z = step / 2.0;
        for (int i = 0; i < categoryCount; i++) {
            rates[currentRateArrayIndex][i] = distributionModel.quantile(z);
            //System.out.print(rates[i]+"\t");
            z += step;
        }

        if (normalize) computeFactor();

        updateRateCategories = false;
    }

    public double getLogLikelihood() {
        return logDensityNormalizationConstant;
    }
}
