/*
 * TreeTraversal.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;

import java.util.*;


/**
 * @author Andrew Rambaut
 * @author Marc A Suchard
 */
public abstract class TreeTraversal {

    public enum TraversalType {
        POST_ORDER,          // likelihood
        REVERSE_LEVEL_ORDER, // likelihood
        PRE_ORDER            // simulation
    }

    protected TreeTraversal(final Tree treeModel,
                            final BranchRateModel branchRateModel,
                            final TraversalType traversalType) {
        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.traversalType = traversalType;

        updateNode = new boolean[treeModel.getNodeCount()];
        updateAllNodes();
    }

    public abstract void dispatchTreeTraversalCollectBranchAndNodeOperations();

    public final Tree getTree() {
        return treeModel;
    }

    public final void setAllNodesUpdated() {
        Arrays.fill(updateNode, false);
        updateAllNodes = false;
    }

    public final void updateAllNodes() {
        Arrays.fill(updateNode, true);
        updateAllNodes = true;
    }

    public final void updateNode(final NodeRef node) {
        updateNode[node.getNumber()] = true;
    }

    public final void updateNodeAndChildren(final NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            final NodeRef child = treeModel.getChild(node, i);
            updateNode[child.getNumber()] = true;
        }
    }

    public final void updateNodeAndDescendents(final NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            final NodeRef child = treeModel.getChild(node, i);
            updateNodeAndDescendents(child);
        }
    }

    public final void updateNodeAndAncestors(final NodeRef node) {
        updateNode[node.getNumber()] = true;

        if (!treeModel.isRoot(node)) {
            final NodeRef parent = treeModel.getParent(node);
            updateNodeAndAncestors(parent);
        }
    }

    protected final double computeBranchLength(final Tree tree, final NodeRef node) {
        final double branchRate;

        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }
        final double parentHeight = tree.getNodeHeight(tree.getParent(node));
        final double nodeHeight = tree.getNodeHeight(node);

        // Get the operational time of the branch
        final double branchLength = branchRate * (parentHeight - nodeHeight);

        assert branchLength >= 0.0 : "Negative branch length: " + branchLength + " for node " +
                node.getNumber() + (tree.isExternal(node) ?
                " (" + tree.getNodeTaxon(node).getId() + ")" : "");

        return branchLength;
    }

    protected final Tree treeModel;
    protected final BranchRateModel branchRateModel;
    protected final boolean[] updateNode;
    protected boolean updateAllNodes;

    protected final TraversalType traversalType;
}

