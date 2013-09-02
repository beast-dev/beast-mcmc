/*
 * CountableMixtureBranchRates.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodelxml.branchratemodel.CountableMixtureBranchRatesParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class CountableMixtureBranchRates extends AbstractBranchRateModel implements Loggable {

    private final Parameter ratesParameter;
    private final Parameter rateCategoryParameter;

    private final int categoryCount;

    public CountableMixtureBranchRates(TreeModel treeModel, Parameter ratesParameter, Parameter rateCategoryParameter) {

        super(CountableMixtureBranchRatesParser.COUNTABLE_CLOCK_BRANCH_RATES);

        this.rateCategoryParameter = rateCategoryParameter;
        rateCategories = new TreeParameterModel(treeModel, rateCategoryParameter, false);

        categoryCount = ratesParameter.getDimension();


        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        rateCategoryParameter.addBounds(bound);

        for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            int index = MathUtils.nextInt(categoryCount);
            rateCategoryParameter.setParameterValue(i, index);
        }

        addModel(rateCategories);

//        updateRateCategories = true;
//
//        // Each parameter take any value in [1, \ldots, categoryCount]
//        // NB But this depends on the transition kernel employed.  Using swap-only results in a different constant
//        logDensityNormalizationConstant = -rateCategoryParameter.getDimension() * Math.log(categoryCount);
//
//        treeParameter = new TreeParameterModel(treeModel, allocationParameter, false);

        this.ratesParameter = ratesParameter;
//        this.allocationParameter = allocationParameter;

        addVariable(ratesParameter);
//        addVariable(allocationParameter);
    }

    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[ratesParameter.getDimension()];
        for (int i = 0; i < ratesParameter.getDimension(); ++i) {
            columns[i] = new OccupancyColumn(i);
        }

        return columns;
    }

    class OccupancyColumn extends NumberColumn {
        private final int index;

        public OccupancyColumn(int index) {
            super("Occupancy");
            this.index = index;
        }

        public double getDoubleValue() {
            int occupancy = 0;
            for (int j = 0; j < rateCategoryParameter.getDimension(); ++j) {
                if (Math.round(rateCategoryParameter.getParameterValue(j)) == index) {
                    occupancy++;
                }
            }
            return occupancy;
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(null, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {


        assert !tree.isRoot(node) : "root node doesn't have a rate!";

//        if (updateRateCategories) {
//            setupRates();
//        }

        int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));

        //System.out.println(rates[rateCategory] + "\t"  + rateCategory);
//        return rates[currentRateArrayIndex][rateCategory] * scaleFactor;

//        final int i = node.getNumber();
//        final int k = (int) allocationParameter.getParameterValue(i);
        return ratesParameter.getParameterValue(rateCategory);
    }

    private final TreeParameterModel rateCategories;
}