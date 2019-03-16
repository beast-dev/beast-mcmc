/*
 * NodeHeightTransform.java
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
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodelxml.continuous.hmc.NodeHeightTransformParser;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
abstract class NodeHeightTransformDelegate extends AbstractModel {
    protected TreeModel tree;
    protected Parameter nodeHeights;
    protected TreeParameterModel indexHelper;

    public NodeHeightTransformDelegate(TreeModel treeModel,
                                       Parameter nodeHeights) {
        super(NodeHeightTransformParser.NAME);
        this.tree = treeModel;
        this.nodeHeights = nodeHeights;
        indexHelper = new TreeParameterModel(treeModel, new Parameter.Default(tree.getNodeCount() - 1), false);
    }

    public void setNodeHeights(double[] nodeHeights) {
        if (nodeHeights.length != this.nodeHeights.getDimension()) {
            throw new RuntimeException("Dimension mismatch!");
        }
        for (int i = 0; i < nodeHeights.length; i++) {
            this.nodeHeights.setParameterValueQuietly(i, nodeHeights[i]);
        }
        tree.pushTreeChangedEvent();
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    abstract double[] transform(double[] values, int from, int to);

    abstract double[] inverse(double[] values, int from, int to);

    abstract String getReport();

    public static class CoalescentIntervals extends NodeHeightTransformDelegate {

        private GMRFSkyrideLikelihood skyrideLikelihood;
        private Parameter coalescentIntervals;
        private OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping;

        public CoalescentIntervals(TreeModel treeModel,
                                   Parameter nodeHeights,
                                   Parameter coalescentIntervals,  // TODO probably don't need
                                   GMRFSkyrideLikelihood skyrideLikelihood) {

            super(treeModel, nodeHeights);

            this.skyrideLikelihood = skyrideLikelihood;
            this.coalescentIntervals = coalescentIntervals;
            this.intervalNodeMapping = skyrideLikelihood.getIntervalNodeMapping();
            addVariable(coalescentIntervals);

            addVariable(nodeHeights);
            this.proxyValuesKnown = false;
        }

        @Override
        double[] transform(double[] values, int from, int to) {
            setNodeHeights(values);
            skyrideLikelihood.setupCoalescentIntervals();
            return coalescentIntervals.getParameterValues();
        }

        @Override
        double[] inverse(double[] values, int from, int to) {
            if (values.length != coalescentIntervals.getDimension()) {
                throw new RuntimeException("Dimension mismatch!");
            }

            double currentHeight = 0.0;

            for (int i = 0; i < values.length; i++) {
                int[] nodeNumbers = intervalNodeMapping.getNodeNumbersForInterval(i);
                currentHeight += values[i];
                TreeModel.Node node = (TreeModel.Node) tree.getNode(nodeNumbers[nodeNumbers.length - 1]);
                node.heightParameter.setParameterValueQuietly(0, currentHeight);
            }
            tree.pushTreeChangedEvent();
            return nodeHeights.getParameterValues();
        }

        @Override
        String getReport() {
            return null;
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {

        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == coalescentIntervals) {
                // TODO I believe this is unnecessary once we start using the proxy
                inverse(coalescentIntervals.getParameterValues(), 0, coalescentIntervals.getDimension());
            }
            proxyValuesKnown = false;
        }

        private boolean proxyValuesKnown;

        private Parameter createProxyForCoalescentIntervals() {
            return new Parameter.Proxy("coalescentIntervals",
                    skyrideLikelihood.getCoalescentIntervalDimension()) {

                private double[] proxy;

                {
                    proxy = new double[dim];
                }

                @Override
                public double getParameterValue(int dim) {
                    updateCoalescentIntervals();
                    return proxy[dim];
                }

                @Override
                public void setParameterValue(int dim, double value) { // This function is very expensive, avoid repeated calls
                    setParameterValueQuietly(dim, value);
                    updateAllNodeHeights();
                }

                @Override
                public void setParameterValueQuietly(int dim, double value) {
                    proxy[dim] = value;
                }

                @Override
                public void setParameterValueNotifyChangedAll(int dim, double value) {
                    setParameterValue(dim, value);
                }

                @Override
                public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {
                    updateAllNodeHeights();
                }

                private void updateCoalescentIntervals() {
                    if (!proxyValuesKnown) {
                         System.arraycopy(skyrideLikelihood.getCoalescentIntervals(), 0,
                                 proxy, 0, proxy.length);
                         proxyValuesKnown = true;
                     }
                }

                private void updateAllNodeHeights() {
                    updateCoalescentIntervals();
                    inverse(proxy, 0, proxy.length);
                }
            };
        }
    }


    public static class Ratios extends NodeHeightTransformDelegate {
        private Parameter ratios;
        private final LikelihoodTreeTraversal epochConstructionTraversal;
        private final SimulationTreeTraversal nodeHeightUpdateTraversal;
        private Map<Integer, Epoch> nodeEpochMap = new HashMap<Integer, Epoch>();
        private List<Epoch> epochs = new ArrayList<Epoch>();

        public Ratios(TreeModel treeModel,
                      Parameter nodeHeights,
                      Parameter ratios,
                      BranchRateModel branchRateModel) {

            super(treeModel, nodeHeights);

            this.ratios = ratios;
            epochConstructionTraversal = new LikelihoodTreeTraversal(
                    tree,
                    branchRateModel,
                    TreeTraversal.TraversalType.POST_ORDER);

            nodeHeightUpdateTraversal = new SimulationTreeTraversal(
                    tree,
                    branchRateModel,
                    TreeTraversal.TraversalType.PRE_ORDER);


            addModel(treeModel);
            addVariable(nodeHeights);
            addVariable(ratios);
            constructEpochs();
        }

        private void constructEpochs() {
            nodeEpochMap.clear();
            epochs.clear();

            epochConstructionTraversal.updateAllNodes();
            epochConstructionTraversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

            final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = epochConstructionTraversal.getNodeOperations();

            for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {
                final NodeRef node = tree.getNode(op.getNodeNumber());
                final NodeRef leftChild = tree.getNode(op.getLeftChild());
                final NodeRef rightChild = tree.getNode(op.getRightChild());

                final double leftAnchorHeight = getAnchorTipHeight(leftChild);
                final double rightAnchorHeight = getAnchorTipHeight(rightChild);

                if (tree.isRoot(node)) {
                    nodeEpochMap.get(leftChild.getNumber()).endEpoch(node, null);
                    nodeEpochMap.get(rightChild.getNumber()).endEpoch(node, null);
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
                    final int ratioNum = indexHelper.getParameterIndexFromNodeNumber(node.getNumber()) - tree.getExternalNodeCount();
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
                    final double ratio = ratios.getParameterValue(indexHelper.getParameterIndexFromNodeNumber(node.getNumber()) - tree.getExternalNodeCount());
                    final double nodeHeight = ratio * (tree.getNodeHeight(parentNode) - epoch.getAnchorTipHeight()) + epoch.getAnchorTipHeight();
//                tree.setNodeHeight(node, nodeHeight);
                    nodeHeights.setParameterValueQuietly(indexHelper.getParameterIndexFromNodeNumber(node.getNumber()) - tree.getExternalNodeCount(), nodeHeight);
                }
            }
            tree.pushTreeChangedEvent();
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
        double[] transform(double[] values, int from, int to) {
            setNodeHeights(values);
            updateRatios();
            return getRatios();
        }

        @Override
        double[] inverse(double[] values, int from, int to) {
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

        private class Epoch implements Comparable {
            private final NodeRef anchorTipNode;
            private List<NodeRef> internalNodes = new ArrayList<NodeRef>();
            private Epoch lastEpoch;
            private NodeRef connectingNode;

            private Epoch(NodeRef anchorTipNode) {
                this.anchorTipNode = anchorTipNode;
                epochs.add(this);
            }

            public double getAnchorTipHeight() {
                return tree.getNodeHeight(anchorTipNode);
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

}