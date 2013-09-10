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
import dr.evolution.util.TaxonList;
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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class CountableMixtureBranchRates extends AbstractBranchRateModel implements Loggable {

    private final Parameter ratesParameter;
    private final TreeModel treeModel;
    private final AbstractBranchRateModel randomEffectsModel;

    private final int categoryCount;

    public CountableMixtureBranchRates(TreeModel treeModel, Parameter ratesParameter, Parameter rateCategoryParameter) {
        this(treeModel, ratesParameter, rateCategoryParameter, null);
    }

    public CountableMixtureBranchRates(TreeModel treeModel, Parameter ratesParameter, Parameter rateCategoryParameter,
                                       AbstractBranchRateModel randomEffectsModel) {
        this(treeModel, ratesParameter, rateCategoryParameter, randomEffectsModel, false);
    }

    public CountableMixtureBranchRates(TreeModel treeModel, Parameter ratesParameter, Parameter rateCategoryParameter,
                                       AbstractBranchRateModel randomEffects, boolean inLogSpace) {
        super(CountableMixtureBranchRatesParser.COUNTABLE_CLOCK_BRANCH_RATES);

        this.treeModel = treeModel;
        rateCategories = new TreeParameterModel(treeModel, rateCategoryParameter, false);
        categoryCount = ratesParameter.getDimension();

        for (int i = 0; i < rateCategoryParameter.getDimension(); ++i) {
            rateCategoryParameter.setParameterValue(i, 0.0);
        }

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        rateCategoryParameter.addBounds(bound);

        addModel(rateCategories);
        addModel(treeModel);
        this.ratesParameter = ratesParameter;
        addVariable(ratesParameter);

        // Handle random effects
        this.randomEffectsModel = randomEffects;
        if (randomEffectsModel != null) {
            addModel(randomEffectsModel);
        }
        // TODO Check that randomEffectsModel mean is zero

        modelInLogSpace = inLogSpace;
    }

    public double getLogLikelihood() {
        return (randomEffectsModel != null) ? randomEffectsModel.getLogLikelihood() : 0.0;
    }

    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[ratesParameter.getDimension()];
        for (int i = 0; i < ratesParameter.getDimension(); ++i) {
            columns[i] = new OccupancyColumn(i);
        }

        return columns;
    }

    public void setClade(TaxonList taxonList, int rateCategory, boolean includeStem, boolean excludeClade) throws Tree.MissingTaxonException {
        if (!excludeClade) {
            throw new IllegalArgumentException("Excluding clades not yet implemented in countable branch rate mixture models.");
        }

        if (!includeStem) {
            throw new IllegalArgumentException("Excluding stems not yet implemented in countable branch rate mixture models.");
        }

        Set<String> leafSet = Tree.Utils.getLeavesForTaxa(treeModel, taxonList);
        if (leafSetList == null) {
            leafSetList = new ArrayList<CladeContainer>();
        }
        leafSetList.add(new CladeContainer(leafSet, rateCategory));
        cladesChanged = true;
    }

    private class CladeContainer {
        private Set<String> leafSet;
        private int rateCategory;

        public CladeContainer(Set<String> leafSet, int rateCategory) {
            this.leafSet = leafSet;
            this.rateCategory = rateCategory;
        }

        public Set<String> getLeafSet() {
            return leafSet;
        }

        public int getRateCategory() {
            return rateCategory;
        }
    }

    private void updateCladeRateCategories() {
        if (leafSetList != null) {
            for (NodeRef node : treeModel.getNodes()) {
                if (node != treeModel.getRoot()) {
                    rateCategories.setNodeValue(treeModel, node, 0.0);
                }
            }
            for (CladeContainer clade : leafSetList) {
                NodeRef node = Tree.Utils.getCommonAncestorNode(treeModel, clade.getLeafSet());
                if (node != treeModel.getRoot()) {
                    rateCategories.setNodeValue(treeModel, node, clade.getRateCategory());
                }
            }
        }
    }

    public void randomize() {
        for (NodeRef node : treeModel.getNodes()) {
            if (node != treeModel.getRoot()) {
                int index = MathUtils.nextInt(categoryCount);
                rateCategories.setNodeValue(treeModel, node, index);
            }
        }
    }

    private class OccupancyColumn extends NumberColumn {
        private final int index;

        public OccupancyColumn(int index) {
            super("Occupancy");
            this.index = index;
        }

        public double getDoubleValue() {
            int occupancy = 0;
            for (NodeRef node : treeModel.getNodes()) {
                if (node != treeModel.getRoot()) {
                    if (Math.round(rateCategories.getNodeValue(treeModel, node)) == index) {
                        occupancy++;
                    }
                }
            }
            return occupancy;
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == rateCategories) {
//            fireModelChanged();
        } else if (model == treeModel) {
            cladesChanged = true;
            fireModelChanged();
        } else if (model == randomEffectsModel) {
            if (object == randomEffectsModel) {
                fireModelChanged();
            } else if (object == null) {
                fireModelChanged(null, index);
            } else {
                throw new IllegalArgumentException("Unknown object component!");
            }
        } else {
            throw new IllegalArgumentException("Unknown model component!");
        }
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

        synchronized (this) {
            if (cladesChanged) {
                updateCladeRateCategories();
                cladesChanged = false;
            }
        }

        int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));
        double effect = ratesParameter.getParameterValue(rateCategory);
        if (randomEffectsModel != null) {
            if (modelInLogSpace) {
                effect += randomEffectsModel.getBranchRate(tree, node);
            } else {
                effect *= randomEffectsModel.getBranchRate(tree, node);
            }
        }
        if (modelInLogSpace) {
            effect = Math.exp(effect);
        }
        return effect;
    }

    private final TreeParameterModel rateCategories;
    private boolean cladesChanged = false;
    private List<CladeContainer> leafSetList = null;
    private final boolean modelInLogSpace;
}