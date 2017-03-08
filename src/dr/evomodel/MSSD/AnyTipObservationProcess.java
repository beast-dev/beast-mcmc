/*
 * AnyTipObservationProcess.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.MSSD;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

/**
 * Package: AnyTipObservationProcess
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 18, 2008
 * Time: 6:45:00 PM
 */
public class AnyTipObservationProcess extends AbstractObservationProcess {
    protected double[] u0;
    protected double[] p;

    public AnyTipObservationProcess(String modelName, TreeModel treeModel, PatternList patterns, SiteRateModel siteModel,
                                    BranchRateModel branchRateModel, Parameter mu, Parameter lam) {
        super(modelName, treeModel, patterns, siteModel, branchRateModel, mu, lam);
    }

    public double calculateLogTreeWeight() {
        int L = treeModel.getNodeCount();
        if (u0 == null || p == null) {
            u0 = new double[L];    // probability that the trait at node i survives to no leaf
            p = new double[L];     // probability of survival on the branch ancestral to i
        }
        int i, j, childNumber;
        NodeRef node;
        double logWeight = 0.0;

        double averageRate = getAverageRate();

        for (i = 0; i < L; ++i) {
            p[i] = 1.0 - getNodeSurvivalProbability(i, averageRate);
        }

        TreeUtils.postOrderTraversalList(treeModel, postOrderNodeList);

        for (int postOrderIndex = 0; postOrderIndex < nodeCount; postOrderIndex++) {

            i = postOrderNodeList[postOrderIndex];

            if (i < treeModel.getExternalNodeCount()) { // Is tip
                u0[i] = 0.0;
                logWeight += 1.0 - p[i];
            } else { // Is internal node or root
                u0[i] = 1.0;
                node = treeModel.getNode(i);
                for (j = 0; j < treeModel.getChildCount(node); ++j) {
                    childNumber = treeModel.getChild(node, j).getNumber();
                    u0[i] *= 1.0 - p[childNumber] * (1.0 - u0[childNumber]);
                }
                logWeight += (1.0 - u0[i]) * (1.0 - p[i]);
            }
        }

        return -logWeight * lam.getParameterValue(0) / (getAverageRate() * mu.getParameterValue(0));
    }


    private void setTipNodePatternInclusion() { // These values never change
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);

            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                extantInTipsBelow[i * patternCount + patternIndex] = 1;
                int taxonIndex = patterns.getTaxonIndex(treeModel.getNodeTaxon(node));
                int[] states = dataType.getStates(patterns.getPatternState(taxonIndex, patternIndex));
                for (int state : states) {
                    if (state == deathState) {
                        extantInTipsBelow[i * patternCount + patternIndex] = 0;
                    }
                }
                extantInTips[patternIndex] += extantInTipsBelow[i * patternCount + patternIndex];

            }
        }

        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                nodePatternInclusion[i * patternCount + patternIndex] =
                        (extantInTipsBelow[i * patternCount + patternIndex] >= extantInTips[patternIndex]);
            }
        }
    }

    void setNodePatternInclusion() {

        if (postOrderNodeList == null) {
            postOrderNodeList = new int[nodeCount];
        }

        if (nodePatternInclusion == null) {
            nodePatternInclusion = new boolean[nodeCount * patternCount];
            storedNodePatternInclusion = new boolean[nodeCount * patternCount];
        }

        if (extantInTips == null) {
            extantInTips = new int[patternCount];
            extantInTipsBelow = new int[nodeCount * patternCount];
            setTipNodePatternInclusion();
        }

        // Determine post-order traversal
        TreeUtils.postOrderTraversalList(treeModel, postOrderNodeList);

        // Do post-order traversal
        for (int postOrderIndex = 0; postOrderIndex < nodeCount; postOrderIndex++) {
            NodeRef node = treeModel.getNode(postOrderNodeList[postOrderIndex]);
            final int nChildren = treeModel.getChildCount(node);
            if (nChildren > 0) {
                final int nodeNumber = node.getNumber();
                for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                    extantInTipsBelow[nodeNumber * patternCount + patternIndex] = 0;
                    for (int j = 0; j < nChildren; j++) {
                        final int childIndex = treeModel.getChild(node, j).getNumber();
                        extantInTipsBelow[nodeNumber * patternCount + patternIndex] +=
                                extantInTipsBelow[childIndex * patternCount + patternIndex];
                    }
                }
            }
        }

        for (int i = treeModel.getExternalNodeCount(); i < treeModel.getNodeCount(); ++i) {
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                nodePatternInclusion[i * patternCount + patternIndex] =
                        (extantInTipsBelow[i * patternCount + patternIndex] >= extantInTips[patternIndex]);
            }
        }

        nodePatternInclusionKnown = true;
    }

    private int[] extantInTips;
    private int[] extantInTipsBelow; // Easier to store/restore (later) if 1D array

    private int[] postOrderNodeList;

}
