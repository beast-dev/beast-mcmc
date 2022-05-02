/*
 * DiscontinuousHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Parameter;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class NodeHeightBounds implements GraphicalParameterBound {

    final Parameter nodeHeight;
    final TreeModel treeModel;
    final private int externalNodeCount;

    public NodeHeightBounds(Parameter nodeHeight,
                            TreeModel treeModel) {

        this.nodeHeight = nodeHeight;
        this.treeModel = treeModel;
        this.externalNodeCount = treeModel.getExternalNodeCount();

    }

    @Override
    public Parameter getParameter() {
        return nodeHeight;
    }

    @Override
    public int[] getConnectedParameterIndices(int index) {
        NodeRef nodeRef = treeModel.getNode(index + externalNodeCount);
        int nodeCount = 0;
        for (int i = 0; i < treeModel.getChildCount(nodeRef); i++) {
            if (!treeModel.isExternal(treeModel.getChild(nodeRef, i))) {
                nodeCount++;
            }
        }
        if (!treeModel.isRoot(nodeRef)) {
            nodeCount++;
        }
        int[] connectedIndices = new int[nodeCount];
        for (int i = 0; i < treeModel.getChildCount(nodeRef); i++) {
            NodeRef childNode = treeModel.getChild(nodeRef, i);
            if (!treeModel.isExternal(childNode)) {
                connectedIndices[i] = childNode.getNumber() - externalNodeCount;
            }
        }
        if (!treeModel.isRoot(nodeRef)) {
            connectedIndices[nodeCount - 1] = treeModel.getParent(nodeRef).getNumber() - externalNodeCount;
        }
        return connectedIndices;
    }

    @Override
    public double getFixedLowerBound(int index) {
        NodeRef nodeRef = treeModel.getNode(index + externalNodeCount);
        double lowerBound = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < treeModel.getChildCount(nodeRef); i++) {
            NodeRef childNode = treeModel.getChild(nodeRef, i);
            if (treeModel.isExternal(childNode)) {
                if (treeModel.getNodeHeight(childNode) > lowerBound) {
                    lowerBound = treeModel.getNodeHeight(childNode);
                }
            }
        }
        return lowerBound;
    }

    @Override
    public double getFixedUpperBound(int index) {
        return Double.POSITIVE_INFINITY;
    }

}