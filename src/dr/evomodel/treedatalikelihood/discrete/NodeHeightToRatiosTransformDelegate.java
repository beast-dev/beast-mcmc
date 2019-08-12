/*
 * NodeHeightToRatiosTransformDelegate.java
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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightToRatiosTransformDelegate extends AbstractNodeHeightTransformDelegate {

    private Parameter ratios;
    private final LikelihoodTreeTraversal postOrderTraversal;
    private final SimulationTreeTraversal nodeHeightUpdateTraversal;
    private Map<Integer, Epoch> nodeEpochMap = new HashMap<Integer, Epoch>();
    private List<Epoch> epochs = new ArrayList<Epoch>();

    public NodeHeightToRatiosTransformDelegate(TreeModel treeModel,
                  Parameter nodeHeights,
                  Parameter ratios,
                  BranchRateModel branchRateModel) {

        super(treeModel, nodeHeights);

        this.ratios = ratios;
        postOrderTraversal = new LikelihoodTreeTraversal(
                tree,
                branchRateModel,
                TreeTraversal.TraversalType.POST_ORDER);

        nodeHeightUpdateTraversal = new SimulationTreeTraversal(
                tree,
                branchRateModel,
                TreeTraversal.TraversalType.PRE_ORDER);


        addModel(treeModel);
        addVariable(ratios);
        constructEpochs();
    }

    private void constructEpochs() {
        nodeEpochMap.clear();
        epochs.clear();

        postOrderTraversal.updateAllNodes();
        postOrderTraversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

        final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = postOrderTraversal.getNodeOperations();

        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {
            final NodeRef node = tree.getNode(op.getNodeNumber());
            final NodeRef leftChild = tree.getNode(op.getLeftChild());
            final NodeRef rightChild = tree.getNode(op.getRightChild());

            final double leftAnchorHeight = getAnchorTipHeight(leftChild);
            final double rightAnchorHeight = getAnchorTipHeight(rightChild);

            if (tree.isRoot(node)) {
                if (!tree.isExternal(leftChild)) {
                    nodeEpochMap.get(leftChild.getNumber()).endEpoch(node, null);
                }
                if (!tree.isExternal(rightChild)) {
                    nodeEpochMap.get(rightChild.getNumber()).endEpoch(node, null);
                }
            } else {
                if (rightAnchorHeight > leftAnchorHeight) {
                    addToEpoch(node, rightChild, leftChild);
                } else {
                    addToEpoch(node, leftChild, rightChild);
                }
            }
        }
    }

    private void addToEpoch(NodeRef node, NodeRef anchorChild, NodeRef otherChild) {
        Epoch epoch = nodeEpochMap.get(anchorChild.getNumber());
        if (epoch == null) {
            if (!tree.isExternal(anchorChild)) {
                throw new RuntimeException("Internal node should be assigned to an epoch already.");
            }
            epoch = new Epoch(anchorChild);
        }
        epoch.addInternalNode(node);
        nodeEpochMap.put(node.getNumber(), epoch);

        Epoch endingEpoch = nodeEpochMap.get(otherChild.getNumber());
        if (endingEpoch != null) {
            endingEpoch.endEpoch(node, epoch);
        }
    }

    private double getAnchorTipHeight(NodeRef node) {
        double anchorTipHeight = tree.getNodeHeight(node);
        if (nodeEpochMap.containsKey(node.getNumber())) {
            anchorTipHeight = nodeEpochMap.get(node.getNumber()).getAnchorTipHeight();
        }
        return anchorTipHeight;
    }

    public double[] getRatios() {
        return ratios.getParameterValues();
    }

    public double[] getNodeHeights() {
        return nodeHeights.getParameterValues();
    }

    public void updateRatios() {
        for (Epoch epoch : epochs) {
            double previousNodeHeight = tree.getNodeHeight(epoch.getConnectingNode());
            final double anchorNodeHeight = epoch.getAnchorTipHeight();
            for (NodeRef node : epoch.getInternalNodes()) {
                final int ratioNum = getRatiosIndex(node);
                final double currentNodeHeight = tree.getNodeHeight(node);
                ratios.setParameterValueQuietly(ratioNum, (currentNodeHeight - anchorNodeHeight) / (previousNodeHeight - anchorNodeHeight));
                previousNodeHeight = currentNodeHeight;
            }
        }
    }

    public void setRatios(double[] ratios) {
        for (int i = 0; i < ratios.length; i++) {
            this.ratios.setParameterValueQuietly(i, ratios[i]);
        }
    }

    public void updateNodeHeights() {
        nodeHeightUpdateTraversal.updateAllNodes();
        nodeHeightUpdateTraversal.dispatchTreeTraversalCollectBranchAndNodeOperations();
        final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = nodeHeightUpdateTraversal.getNodeOperations();
        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {
            NodeRef node = tree.getNode(op.getLeftChild());
            if (!tree.isRoot(node) && !tree.isExternal(node)) {
                Epoch epoch = nodeEpochMap.get(node.getNumber());
                final NodeRef parentNode = tree.getParent(node);
                final double ratio = ratios.getParameterValue(getRatiosIndex(node));
                final double nodeHeight = ratio * (tree.getNodeHeight(parentNode) - epoch.getAnchorTipHeight()) + epoch.getAnchorTipHeight();
//                tree.setNodeHeight(node, nodeHeight);
                nodeHeights.setParameterValueQuietly(getRatiosIndex(node), nodeHeight);
            }
        }
        tree.pushTreeChangedEvent();
    }

    private int getRatiosIndex(NodeRef node) {
        return indexHelper.getParameterIndexFromNodeNumber(node.getNumber()) - tree.getExternalNodeCount();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            if (object instanceof TreeChangedEvent) {
                TreeModel.TreeChangedEvent changedEvent = (TreeModel.TreeChangedEvent) object;
                if (changedEvent.isTreeChanged()) {
                    constructEpochs();
                    updateRatios();
                }
            }
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratios) {
            updateNodeHeights();
        } else if (variable == nodeHeights) {
            updateRatios();
        }
    }

    @Override
    double[] transform(double[] values) {
        setNodeHeights(values);
        updateRatios();
        return getRatios();
    }

    @Override
    double[] inverse(double[] values) {
        setRatios(values);
        updateNodeHeights();
        return getNodeHeights();
    }

    @Override
    String getReport() {
        updateRatios();
        StringBuilder sb = new StringBuilder();
        sb.append("NodeHeight ratios: ").append(new dr.math.matrixAlgebra.Vector(getRatios()));
        sb.append("\n");

//            Parameter testRatios = new Parameter.Default(tree.getNodeCount() - tree.getExternalNodeCount() - 1, 0.99);
//            inverse(testRatios.getParameterValues(), 0, testRatios.getDimension());
//
//            sb.append("New NodeHeights: ").append(new dr.math.matrixAlgebra.Vector(getNodeHeights()));
//            sb.append("\n");
//            sb.append(tree.getNewick()).append("\n");
        return sb.toString();
    }

    @Override
    Parameter getParameter() {
        return ratios;
    }

    @Override
    double getLogJacobian(double[] values) {
        double logJacobian = 0.0;
        for (int i = tree.getExternalNodeCount(); i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                Epoch epoch = nodeEpochMap.get(i);
                logJacobian += Math.log(getNodePartial(node));
            }
        }
        return logJacobian;
    }

    private int getNodeHeightGradientIndex(NodeRef node) {
        return node.getNumber() - tree.getExternalNodeCount();
    }

    @Override
    double[] updateGradientLogDensity(double[] gradient, double[] value) {
        // gradient is wrt nodeHeight, value = ratios
        double[] ratiosGradientUnweightedLogDensity = new double[ratios.getDimension()];
        postOrderTraversal.updateAllNodes();
        postOrderTraversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

        final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = postOrderTraversal.getNodeOperations();

        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {
            final NodeRef node = tree.getNode(op.getNodeNumber());
            final NodeRef leftChild = tree.getNode(op.getLeftChild());
            final NodeRef rightChild = tree.getNode(op.getRightChild());

            final int nodeIndex = getRatiosIndex(node);

            if (!tree.isRoot(node)) {
                final double nodePartial = getNodePartial(node);
                ratiosGradientUnweightedLogDensity[nodeIndex] += nodePartial * gradient[getNodeHeightGradientIndex(node)];
                ratiosGradientUnweightedLogDensity[nodeIndex] += getEpochGradientAddition(node, leftChild, ratiosGradientUnweightedLogDensity);
                ratiosGradientUnweightedLogDensity[nodeIndex] += getEpochGradientAddition(node, rightChild, ratiosGradientUnweightedLogDensity);
            }
        }
        return ratiosGradientUnweightedLogDensity;
    }

    private double getNodePartial(NodeRef node) {
        final int nodeIndex = getRatiosIndex(node);
        return (tree.getNodeHeight(node) - nodeEpochMap.get(node.getNumber()).getAnchorTipHeight()) / ratios.getParameterValue(nodeIndex);
    }

    private double getEpochGradientAddition(NodeRef node, NodeRef child, double[] ratiosGradientUnweightedLogDensity) {
        final int childIndex = getRatiosIndex(child);
        final int nodeIndex = getRatiosIndex(node);
        if (childIndex < 0) {
            return 0.0;
        } else if (nodeEpochMap.get(child.getNumber()) == nodeEpochMap.get(node.getNumber())){
            return ratiosGradientUnweightedLogDensity[childIndex] *
                    ratios.getParameterValue(childIndex) / ratios.getParameterValue(nodeIndex);
        } else {
            return ratiosGradientUnweightedLogDensity[childIndex] * ratios.getParameterValue(childIndex)
                    / (tree.getNodeHeight(node) - nodeEpochMap.get(child.getNumber()).getAnchorTipHeight())
                    * getNodePartial(node);
        }
    }

    private class Epoch implements Comparable {
        private final int anchorTipNodeNumber;
        private List<NodeRef> internalNodes = new ArrayList<NodeRef>();
        private Epoch lastEpoch;
        private NodeRef connectingNode;

        private Epoch(NodeRef anchorTipNode) {
            this.anchorTipNodeNumber = anchorTipNode.getNumber();
            epochs.add(this);
        }

        public double getAnchorTipHeight() {
            return tree.getNodeHeight(tree.getNode(anchorTipNodeNumber));
        }

        public void endEpoch(NodeRef node, Epoch lastEpoch) {
            this.lastEpoch = lastEpoch;
            this.connectingNode = node;
        }

        public void addInternalNode(NodeRef node) {
            internalNodes.add(0, node);
        }

        public List<NodeRef> getInternalNodes() {
            return internalNodes;
        }

        public NodeRef getConnectingNode() {
            return connectingNode;
        }

        @Override
        public int compareTo(Object o) {
            return Double.compare(getAnchorTipHeight(), ((Epoch) o).getAnchorTipHeight());
        }
    }
}
