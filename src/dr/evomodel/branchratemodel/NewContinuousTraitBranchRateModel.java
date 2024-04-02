/*
 * ContinuousTraitBranchRateModel.java
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
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.utilities.TreeTraitLogger;
import dr.evomodelxml.branchratemodel.NewContinuousTraitBranchRateModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

public class NewContinuousTraitBranchRateModel extends AbstractBranchRateModel {

    private final String trait;
    private final int dimension;
    private final Parameter rateParameter;
    private final Parameter ratioParameter;
    private TreeDataLikelihood traitLikelihood;

    public NewContinuousTraitBranchRateModel(String trait, TreeDataLikelihood traitLikelihood, int dimension) {
        super(NewContinuousTraitBranchRateModelParser.TRAIT_BRANCH_RATES);

        this.traitLikelihood = traitLikelihood;
        this.trait = trait;
        this.dimension = dimension;

        this.rateParameter = null;
        this.ratioParameter = null;
        addModel(traitLikelihood);
    }

    public NewContinuousTraitBranchRateModel(String trait, Parameter rateParameter, Parameter ratioParameter) {
        super(NewContinuousTraitBranchRateModelParser.TRAIT_BRANCH_RATES);

        this.trait = trait;
        dimension = 0;
        this.rateParameter = rateParameter;
        this.ratioParameter = ratioParameter;

        if (rateParameter != null) {
            addVariable(rateParameter);
        }

        if (ratioParameter != null) {
            addVariable(ratioParameter);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }


    public double getBranchRate(final Tree tree, final NodeRef node) {
        NodeRef parent = tree.getParent(node);
        if (parent == null) {
            throw new IllegalArgumentException("Root does not have a valid rate");
        }

        double rate = 1.0;
        TreeModel treeModel = (TreeModel) tree;

        if (rateParameter != null) {
            double scale = 1.0;
            double ratio = 1.0;

            if (rateParameter != null) {
                scale = rateParameter.getParameterValue(0);
            }

            if (ratioParameter != null) {
                ratio = ratioParameter.getParameterValue(0);
            }


            // get the log rate for the node and its parent
            double rate1 = ratio * treeModel.getMultivariateNodeTrait(node, trait)[0];
            double rate2 = ratio * treeModel.getMultivariateNodeTrait(parent, trait)[0];

            if (rate1 == rate2) {
                return scale * Math.exp(rate1);
            }

            rate = scale * (Math.exp(rate2) - Math.exp(rate1)) / (rate2 - rate1);
        } else {
            double rate1 =  treeModel.getMultivariateNodeTrait(node, trait)[dimension];
            double rate2 =  treeModel.getMultivariateNodeTrait(parent, trait)[dimension];

            if (rate1 == rate2) {
                return Math.exp(rate1);
            }

            rate = (Math.exp(rate2) - Math.exp(rate1)) / (rate2 - rate1); // TODO Should this not be averaged on the log-scale?
        }

        return rate;
    }

}