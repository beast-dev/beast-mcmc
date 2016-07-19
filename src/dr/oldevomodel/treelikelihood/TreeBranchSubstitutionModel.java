/*
 * TreeBranchSubstitutionModel.java
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * TreeBranchSubstitutionModel - provides a model for making substitution model epoches.
 *
 * @author Marc A. Suchard
 */

@Deprecated // Switching to BEAGLE
public class TreeBranchSubstitutionModel extends AbstractModel {

    public TreeBranchSubstitutionModel(String name) {
        super(name);
    }

    public TreeBranchSubstitutionModel(String name, SiteModel siteModel, SubstitutionModel substModel, BranchRateModel branchModel) {
        super(name);
        this.siteModel = siteModel;
        this.substModel = substModel;
        this.branchModel = branchModel;

        if (siteModel != null)
            addModel(siteModel);
        if (substModel != null)
            addModel(substModel);
        if (branchModel != null)
            addModel(branchModel);
    }

    public void getTransitionProbabilities(Tree tree, NodeRef node, int rateCategory, double[] probs) {

        NodeRef parent = tree.getParent(node);

        final double branchRate = branchModel.getBranchRate(tree, node);

        // Get the operational time of the branch
        final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

        if (branchTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchTime);
        }

        double branchLength = siteModel.getRateForCategory(rateCategory) * branchTime;
        substModel.getTransitionProbabilities(branchLength, probs);
    } // getTransitionProbabilities

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    protected SiteModel siteModel;
    protected SubstitutionModel substModel;
    protected BranchRateModel branchModel;

}
