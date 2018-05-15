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
import dr.evomodelxml.branchratemodel.TimeDependentBranchRatesParser;
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
public class TimeDependentBranchRates extends AbstractBranchRateModel implements Loggable {

    private final Parameter intercept;
    private final Parameter timeCoefficient;
    private final TreeModel treeModel;
    private final List<AbstractBranchRateModel> randomEffectsModels;

    public TimeDependentBranchRates(TreeModel treeModel, Parameter intercept, Parameter timeCoefficient,
                                       List<AbstractBranchRateModel> randomEffects, boolean inLogSpace) {

        super(TimeDependentBranchRatesParser.TIME_DEPENDENT_BRANCH_RATES);

        this.treeModel = treeModel;
        this.intercept = intercept;
        addVariable(intercept);
        this.timeCoefficient = timeCoefficient;
        addVariable(timeCoefficient);

        // Handle random effects
        this.randomEffectsModels = randomEffects;
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels)
                addModel(model);
        }
        // TODO Check that randomEffectsModel means are zero (inherited from countableMixtureBranchRates)

        helper.addTrait(this);
        helper.addTrait(new TreeTrait.D() {

            @Override
            public String getTraitName() {
                return getBranchTimeEffectTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return getBranchTimeEffect(tree, node);
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

    private String getBranchTimeEffectTraitName() {
        return getTraitName() + ".time.effect";
    }

    private String getRandomEffectTraitName() {
        return getTraitName() + ".random.effect";
    }

    private double getMidpointHeight(Tree tree, NodeRef node, boolean log){
        double nodeHeight = tree.getNodeHeight(node);
        double parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
        double midpoint = nodeHeight+(parentNodeHeight-nodeHeight)/2;
        if(log){
            return Math.log(midpoint);
        } else {
            return midpoint;
        }
    }

    private double getBranchTimeEffect(Tree tree, NodeRef node) {
        double nodeHeight = tree.getNodeHeight(node);
        double parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
        return timeCoefficient.getParameterValue(0)*getMidpointHeight(tree, node, true);
    }

    private double getBranchRandomEffect(Tree tree, NodeRef node) {
        double effect = 0;
        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                effect += model.getBranchRate(tree, node);
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

//    public LogColumn[] getColumns() {
//        LogColumn[] columns = new LogColumn[ratesParameter.getDimension()];
//        for (int i = 0; i < ratesParameter.getDimension(); ++i) {
//            columns[i] = new OccupancyColumn(i);
//        }
//
//        return columns;
//    }

//    private class OccupancyColumn extends NumberColumn {
//        private final int index;
//
//        public OccupancyColumn(int index) {
//            super("Occupancy");
//            this.index = index;
//        }
//
//        public double getDoubleValue() {
//            int occupancy = 0;
//            for (NodeRef node : treeModel.getNodes()) {
//                if (node != treeModel.getRoot()) {
//                    if (rateCategories.getBranchCategory(treeModel, node) == index) {
//                        occupancy++;
//                    }
//                }
//            }
//            return occupancy;
//        }
//    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
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

        double rate = intercept.getParameterValue(0);
        rate = rate + timeCoefficient.getParameterValue(0)*getMidpointHeight(tree, node, true);

        if (randomEffectsModels != null) {
            for (AbstractBranchRateModel model : randomEffectsModels) {
                rate += model.getBranchRate(tree, node);
            }
        }
        return  Math.exp(rate);
    }

    private final Helper helper = new Helper();

}