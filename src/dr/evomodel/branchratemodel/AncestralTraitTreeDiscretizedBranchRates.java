/*
 * AncestralTraitTreeDiscretizedBranchRates.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.TransformableTree;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 * @author Paul Bastide
 */

public class AncestralTraitTreeDiscretizedBranchRates extends DiscretizedBranchRates {

    public AncestralTraitTreeDiscretizedBranchRates(
            AncestralTraitTreeModel ancestralTree,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            int overSampling) {
        super((TreeModel) ancestralTree.getOriginalTree(), rateCategoryParameter, model, overSampling);
    }

    public AncestralTraitTreeDiscretizedBranchRates(
            AncestralTraitTreeModel ancestralTree,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            int overSampling,
            boolean normalize,
            double normalizeBranchRateTo,
            boolean randomizeRates,
            boolean keepRates,
            boolean cacheRates) {
        super(ancestralTree.getOriginalTree(), rateCategoryParameter, model, overSampling, normalize, normalizeBranchRateTo, randomizeRates, keepRates, cacheRates);

        // adding the key word to the tree means the keyword will be logged in the
        // header of the tree file.
        ancestralTree.addKeyword("discretized_branch_rates");

        // adding the key word to the model means the keyword will be logged in the
        // header of the logfile.
        this.addKeyword("discretized_branch_rates");
    }

    @Override
    public final double getBranchRate(final Tree tree, final NodeRef node) {
        assert tree instanceof TransformableTree;
        TransformableTree transTree = (TransformableTree) tree;
        if (transTree.isInOriginalTree(node)) {
            return super.getBranchRate(transTree.getOriginalTree(), transTree.getOriginalNode(node));
        } else {
            return 1.0;
        }
    }

    @Override
    public int getBranchRateCategory(final Tree tree, final NodeRef node) {
        assert tree instanceof TransformableTree;
        TransformableTree transTree = (TransformableTree) tree;
        if (transTree.isInOriginalTree(node)) {
            return super.getBranchRateCategory(transTree.getOriginalTree(), transTree.getOriginalNode(node));
        } else {
            return -1;
        }
    }
}
