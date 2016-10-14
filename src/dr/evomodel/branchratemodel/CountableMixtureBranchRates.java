/*
 * CountableMixtureBranchRates.java
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
import dr.evolution.tree.TreeTrait;
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
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class CountableMixtureBranchRates extends AbstractBranchRateModel implements Loggable {

    private final Parameter ratesParameter;
    private final TreeModel treeModel;
    private final List<AbstractBranchRateModel> randomEffectsModels;
    private final int categoryCount;

    public CountableMixtureBranchRates(CountableBranchCategoryProvider rateCategories,
                                       TreeModel treeModel, Parameter ratesParameter,
                                       List<AbstractBranchRateModel> randomEffects, boolean inLogSpace) {
        super(CountableMixtureBranchRatesParser.COUNTABLE_CLOCK_BRANCH_RATES);

        this.treeModel = treeModel;
        categoryCount = ratesParameter.getDimension();
        this.rateCategories = rateCategories;
        rateCategories.setCategoryCount(categoryCount);

        if (rateCategories instanceof Model) {
            addModel((Model)rateCategories);
        }
        this.ratesParameter = ratesParameter;
        addVariable(ratesParameter);

        // Handle random effects
        this.randomEffectsModels = randomEffects;
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels)
            addModel(model);
        }
        // TODO Check that randomEffectsModel means are zero

        modelInLogSpace = inLogSpace;

        helper.addTrait(this);
        helper.addTrait(new TreeTrait.I() {

            @Override
            public String getTraitName() {
                return getCategoryTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Integer getTrait(Tree tree, NodeRef node) {
                return getBranchCategory(tree, node);
            }
        });
        helper.addTrait(new TreeTrait.D() {

            @Override
            public String getTraitName() {
                return getCategoryEffectTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return getBranchCategoryEffect(tree, node);
            }
        });
        helper.addTrait(new TreeTrait.D() {

            @Override
            public String getTraitName() {
                return getCategoryRateTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return getBranchCategoryRate(tree, node);
            }
        });
        helper.addTrait(new TreeTrait.D() {

            @Override
            public String getTraitName() {
                return getRandomEffectTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return getBranchRandomEffect(tree, node);
            }
        });
    }

    private String getCategoryTraitName() {
        return getTraitName() + ".category";
    }

    private String getCategoryEffectTraitName() {
        return getTraitName() + ".category.effect";
    }

    private String getCategoryRateTraitName() {
        return getTraitName() + ".category.rate";
    }

    private String getRandomEffectTraitName() {
        return getTraitName() + ".random.effect";
    }

    private int getBranchCategory(Tree tree, NodeRef node) {
        return rateCategories.getBranchCategory(tree, node);
    }

    private double getBranchCategoryRate(Tree tree, NodeRef node) {
        if (modelInLogSpace) {
            return ratesParameter.getParameterValue(getBranchCategory(tree, node));
        } else {
            return Math.exp(ratesParameter.getParameterValue(getBranchCategory(tree, node)));
        }
    }

    private double getBranchCategoryEffect(Tree tree, NodeRef node) {
        if (modelInLogSpace) {
            return (getBranchCategoryRate(tree, node) - ratesParameter.getParameterValue(0));
        }  else {
            return (getBranchCategoryRate(tree, node) / ratesParameter.getParameterValue(0));
        }
    }

    private double getBranchRandomEffect(Tree tree, NodeRef node) {
        double effect;
        if (modelInLogSpace) {
            effect = 0;
        } else {
            effect = 1;
        }
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                if (modelInLogSpace) {
                    effect += model.getBranchRate(tree, node);
                } else {
                    effect *= model.getBranchRate(tree, node);
                }
            }
        }
        return effect;
    }

    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }

    public double getLogLikelihood() {
        double logLike = 0.0;
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                logLike += model.getLogLikelihood();
            }
        }
        return logLike;
    }

    void test() {
        getTrait(null, null);
    }

    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[ratesParameter.getDimension()];
        for (int i = 0; i < ratesParameter.getDimension(); ++i) {
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
            AbstractBranchRateModel foundModel = findRandomEffectsModel(model);
            if (foundModel != null) {
                if (object == model) {
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
    }

    private AbstractBranchRateModel findRandomEffectsModel(Model model) {
        AbstractBranchRateModel found = null;
        int index = randomEffectsModels.indexOf(model);
        if (index != -1) {
            found = randomEffectsModels.get(index);
        }
        return found;
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
        double effect = ratesParameter.getParameterValue(rateCategory);
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

    private final Helper helper = new Helper();

    private final CountableBranchCategoryProvider rateCategories;
    private final boolean modelInLogSpace;
}