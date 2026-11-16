/*
 * NodeHeightBounds.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
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
        int nodeIndex = 0;
        for (int i = 0; i < treeModel.getChildCount(nodeRef); i++) {
            NodeRef childNode = treeModel.getChild(nodeRef, i);
            if (!treeModel.isExternal(childNode)) {
                connectedIndices[nodeIndex] = childNode.getNumber() - externalNodeCount;
                nodeIndex++;
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

    @Override
    public boolean isBelow(int i, int j) {
        NodeRef nodeRef = treeModel.getNode(i + externalNodeCount);
        return !treeModel.isRoot(nodeRef) && treeModel.getParent(nodeRef).getNumber() - externalNodeCount == j;
    }

    @Override
    public void checkParameterMapping() {
        checkParameterMapping(nodeHeight, "the parameter inside graphicalParameterBounds");
    }

    @Override
    public void checkParameterMapping(Parameter parameter) {
        checkParameterMapping(parameter, "the operator's parameter");
    }

    private void checkParameterMapping(Parameter parameter, String role) {
        // Coordinate i is taken to be the height of node i + externalNodeCount (the internal nodes in node-number
        // order, root included). That holds for a nodeHeightProxyParameter with rootNode="true" on this tree only.
        final int internalNodeCount = treeModel.getInternalNodeCount();
        if (!(parameter instanceof NodeHeightProxyParameter)
                || ((NodeHeightProxyParameter) parameter).getTree() != treeModel
                || parameter.getDimension() != internalNodeCount) {
            final String rule = "graphicalParameterBounds on tree '" + treeModel.getId() +
                    "' take coordinate i to be the height of node i + " + externalNodeCount +
                    " (the internal nodes, root included), which holds only for a nodeHeightProxyParameter with" +
                    " rootNode=\"true\" on that tree";
            if (parameter == null) {
                throw new IllegalArgumentException(rule + ", but " + role + " is missing");
            }
            throw new IllegalArgumentException(rule + ", but " + role + ", " + parameter.getParameterName() +
                    (parameter instanceof NodeHeightProxyParameter ?
                            " on tree '" + ((NodeHeightProxyParameter) parameter).getTree().getId() + "'" : "") +
                    " (dimension " + parameter.getDimension() + "; '" + treeModel.getId() + "' has " +
                    internalNodeCount + " internal nodes), is not one: treeModel.allInternalNodeHeights and" +
                    " nodeHeightProxyParameter with rootNode=\"false\" map coordinates to the wrong nodes");
        }
    }

    @Override
    public boolean isWithinBounds(double[] values) {
        for (int i = 0; i < values.length; i++) {
            NodeRef nodeRef = treeModel.getNode(i + externalNodeCount);
            double upperBound = treeModel.isRoot(nodeRef) ? Double.POSITIVE_INFINITY :
                    values[treeModel.getParent(nodeRef).getNumber() - externalNodeCount];
            if (!(values[i] >= getFixedLowerBound(i) && values[i] <= upperBound)) {
                return false;
            }
        }
        return true;
    }

}