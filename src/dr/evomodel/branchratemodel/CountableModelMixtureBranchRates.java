/*
 * CountableModelMixtureBranchRates.java
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
import dr.evomodelxml.branchratemodel.CountableMixtureBranchRatesParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

/**
 * A modification of CountableMixtureBranchRates so the fixed effects come from a list of branchRateModels
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class CountableModelMixtureBranchRates extends AbstractBranchRateModel implements Loggable {

    private final List<AbstractBranchRateModel> fixedEffectsModels;
    private final TreeModel treeModel;
    private final List<AbstractBranchRateModel> randomEffectsModels;
    private final int categoryCount;

    public CountableModelMixtureBranchRates(CountableBranchCategoryProvider rateCategories,
                                            TreeModel treeModel,
                                            List<AbstractBranchRateModel> fixedEffects,
                                            List<AbstractBranchRateModel> randomEffects, boolean inLogSpace) {
        super(CountableMixtureBranchRatesParser.COUNTABLE_CLOCK_BRANCH_RATES);

        this.treeModel = treeModel;
        categoryCount = fixedEffects.size();
        this.rateCategories = rateCategories;
        rateCategories.setCategoryCount(categoryCount);

        if (rateCategories instanceof Model) {
            addModel((Model)rateCategories);
        }
        this.fixedEffectsModels = fixedEffects;
        for (AbstractBranchRateModel model : fixedEffectsModels) {
            addModel(model);
        }

        // Handle random effects
        this.randomEffectsModels = randomEffects;
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                addModel(model);
            }
        }
        // TODO Check that randomEffectsModel means are zero

        modelInLogSpace = inLogSpace;
    }

    public double getLogLikelihood() {
        double logLike = 0.0;
        for (AbstractBranchRateModel model : fixedEffectsModels) {
            logLike += model.getLogLikelihood();
        }
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                logLike += model.getLogLikelihood();
            }
        }
        return logLike;
    }

    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[categoryCount];
        for (int i = 0; i < categoryCount; ++i) {
            columns[i] = new OccupancyColumn(i);
        }

        return columns;
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
                    if (rateCategories.getBranchCategory(treeModel, node) == index) {
                        occupancy++;
                    }
                }
            }
            return occupancy;
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == rateCategories) {
            fireModelChanged();
        } else {
            AbstractBranchRateModel foundModel = findModel(model);
            if (foundModel != null) {
                if (object == model) {
                    fireModelChanged();
                } else if (object == null) {
                    fireModelChanged(null, index);
                } else {
//                    throw new IllegalArgumentException("Unknown object component!");
                    fireModelChanged();
                }
            } else {
                throw new IllegalArgumentException("Unknown model component!");
            }
        }
    }

    private AbstractBranchRateModel findModel(Model model) {
        if (randomEffectsModels != null) {
            int index = randomEffectsModels.indexOf(model);
            if (index != -1) {
                return randomEffectsModels.get(index);
            }
        }
        int index = fixedEffectsModels.indexOf(model);
        if (index != -1) {
            return fixedEffectsModels.get(index);
        }
        return null;
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

        int rateCategory = rateCategories.getBranchCategory(tree, node);

        AbstractBranchRateModel fixedModel = fixedEffectsModels.get(rateCategory);
        double effect = fixedModel.getBranchRate(tree, node);

        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                if (modelInLogSpace) {
                    effect += model.getBranchRate(tree, node);
                } else {
                    effect *= model.getBranchRate(tree, node);
                }
            }
        }
        if (modelInLogSpace) {
            effect = Math.exp(effect);
        }
        return effect;
    }

    private final CountableBranchCategoryProvider rateCategories;
    private final boolean modelInLogSpace;
}