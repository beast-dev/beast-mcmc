/*
 * StarTreeLikelihood.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.treelikelihood;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Complete the likelihood of sequence data given a star-tree with height = treeModel.rootHeight
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class StarTreeLikelihood extends OldBeagleTreeLikelihood {

    public StarTreeLikelihood(PatternList patternList, TreeModel treeModel,
                              BranchSubstitutionModel branchSubstitutionModel, SiteRateModel siteRateModel,
                              BranchRateModel branchRateModel, boolean useAmbiguities,
                              PartialsRescalingScheme rescalingScheme,
                               Map<Set<String>, Parameter> partialsRestrictions) {
        super(patternList, treeModel, branchSubstitutionModel, siteRateModel, branchRateModel, null, useAmbiguities,
                rescalingScheme, partialsRestrictions);

        // Modify tree into star
        forceStarTree(treeModel);

        rootHeightParameter = treeModel.getRootHeightParameter();

        // Print info to screen
        StringBuilder sb = new StringBuilder();
        sb.append("Building a star-tree sequence likelihood model.  Please cite:");
        Logger.getLogger("dr.app.beagle.evomodel").info(sb.toString());
    }

    private void forceStarTree(TreeModel treeModel) {
        double rootHeight = treeModel.getNodeHeight(treeModel.getRoot());
        for (int i = 0; i < treeModel.getInternalNodeCount(); ++i) {
            NodeRef node = treeModel.getInternalNode(i);
            if (node != treeModel.getRoot()) {
                treeModel.setNodeHeight(node, rootHeight);
            }
        }
        fixedTree = true;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        boolean validUpdate = true;

        if (fixedTree && model == treeModel) {
            if (SAMPLE_ROOT) {
                if (object instanceof TreeModel.TreeChangedEvent) {
                    TreeModel.TreeChangedEvent event = (TreeModel.TreeChangedEvent) object;
                    if (event.getNode() != treeModel.getRoot()) {
                        validUpdate = false;
                    }
                } else if (object != rootHeightParameter) {
                    validUpdate = false;
                }
            } else {
                validUpdate = false;
            }

//            if (validUpdate) {
//                forceStarTree(treeModel);
//            }
        }
        if (validUpdate) {
            super.handleModelChangedEvent(model, object, index);
        } else {
            throw new IllegalArgumentException(
                    "Invalid operator; do not sample tree structure or internal node heights");
        }
    }

    private final Parameter rootHeightParameter;
    private boolean fixedTree = false;

    private static boolean SAMPLE_ROOT = false;
}
