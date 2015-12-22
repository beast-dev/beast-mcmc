/*
 * NewAbstractTreeLikelihood.java
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

package dr.app.beagle.evomodel.newtreelikelihood;

import beagle.Beagle;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.Reportable;

/**
 * AbstractTreeLikelihood - a base class for likelihood calculators of sites on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id: AbstractTreeLikelihood.java,v 1.16 2005/06/07 16:27:39 alexei Exp $
 */

public abstract class NewAbstractTreeLikelihood extends AbstractModelLikelihood implements Reportable {

    protected static final boolean COUNT_TOTAL_OPERATIONS = false;

    public NewAbstractTreeLikelihood(String name, TreeModel treeModel) {

        super(name);

        this.treeModel = treeModel;
        addModel(treeModel);

        nodeCount = treeModel.getNodeCount();

        this.tipCount = treeModel.getExternalNodeCount();

        internalNodeCount = nodeCount - tipCount;

        // one partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

        likelihoodKnown = false;
    }


    /**
     * Set update flag for a node and its children
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNode[child.getNumber()] = true;
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and all its descendents
     */
    protected void updateNodeAndDescendents(NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNodeAndDescendents(child);
        }

        likelihoodKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }
        likelihoodKnown = false;
    }

//    protected abstract boolean hasBranchSpecificEvolutionaryProcess();

//    protected boolean updateBranchSpecificEvolutionaryProcess(final Tree tree, final NodeRef parent,
//                                                              final NodeRef node, final boolean flip) {
//
//        final double branchRate;
//
//        synchronized (branchRateModel) {
//            branchRate = branchRateModel.getBranchRate(tree, node);
//        }
//        final double parentHeight = tree.getNodeHeight(parent);
//        final double nodeHeight = tree.getNodeHeight(node);
//
//        // Get the operational time of the branch
//        final double branchLength = branchRate * (parentHeight - nodeHeight);
//        if (branchLength < 0.0) {
//            throw new RuntimeException("Negative branch length: " + branchLength);
//        }
//
//        if (flip) {
//            substitutionModelDelegate.flipMatrixBuffer(nodeNum);
//        }
//        branchUpdateIndices[branchUpdateCount] = nodeNum;
//        branchLengths[branchUpdateCount] = branchLength;
//        branchUpdateCount++;
//
//        return true;
//    }

    protected boolean updateBranchSpecificEvolutionaryProcess(final Tree tree, final NodeRef parent,
                                                              final NodeRef node, final boolean flip) {
        return false;
    }


//    protected void handlePartialsScaling(final int nodeNum, final int x) {
//        if (useScaleFactors) {
//            // get the index of this scaling buffer
//            int n = nodeNum - tipCount;
//
//            if (recomputeScaleFactors) {
//                // flip the indicator: can take either n or (internalNodeCount + 1) - n
//                scaleBufferHelper.flipOffset(n);
//
//                // store the index
//                scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);
//
//                operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
//                operations[x + 2] = Beagle.NONE;
//
//            } else {
//                operations[x + 1] = Beagle.NONE;
//                operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
//            }
//
//        } else {
//
//            if (useAutoScaling) {
//                scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
//            }
//            operations[x + 1] = Beagle.NONE; // Not using scaleFactors
//            operations[x + 2] = Beagle.NONE;
//        }
//    }

    protected void handlePartialsScaling(final int nodeNum, final int x) {
        // Do nothing
    }

//    protected void handleEvolutionaryProcess(final NodeRef child1, final NodeRef child2, final int x) {
//        operations[x + 4] = substitutionModelDelegate.getMatrixIndex(child1.getNumber()); // source matrix 1
//        operations[x + 6] = substitutionModelDelegate.getMatrixIndex(child2.getNumber()); // source matrix 2
//    }

    protected void handleEvolutionaryProcess(final NodeRef child1, final NodeRef child2, final int x) {
        // Do nothing
    }


//    protected void handleRestrictedPartials(final int nodeNum) {
//        if (hasRestrictedPartials) {
//            // Test if this set of partials should be restricted
//            if (updateRestrictedNodePartials) {
//                // Recompute map
//                computeNodeToRestrictionMap();
//                updateRestrictedNodePartials = false;
//            }
//            if (partialsMap[nodeNum] != null) {
//
//            }
//        }
//    }

    protected void handleRestrictedPartials(final int nodeNum) {
        // Do nothing
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param tree           tree
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    protected boolean traverse(final Tree tree, final NodeRef node, final int[] operatorNumber, final boolean flip) {

        boolean update = false;

        final int nodeNum = node.getNumber();

        final NodeRef parent = tree.getParent(node);

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the evolutionary process for this branch
        if (parent != null && updateNode[nodeNum]) {
            update = updateBranchSpecificEvolutionaryProcess(tree, parent, node, flip);
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            final NodeRef child1 = tree.getChild(node, 0);
            final int[] op1 = {-1};
            final boolean update1 = traverse(tree, child1, op1, flip);

            final NodeRef child2 = tree.getChild(node, 1);
            final int[] op2 = {-1};
            final boolean update2 = traverse(tree, child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                final int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                handlePartialsScaling(nodeNum, x);

//                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNumber()); // source node 1
//                operations[x + 4] = substitutionModelDelegate.getMatrixIndex(child1.getNumber()); // source matrix 1
//                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNumber()); // source node 2
//                operations[x + 6] = substitutionModelDelegate.getMatrixIndex(child2.getNumber()); // source matrix 2

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNumber()); // source node 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNumber()); // source node 2

                handleEvolutionaryProcess(child1, child2, x);

                operationCount[operationListCount]++;

                update = true;

                handleRestrictedPartials(nodeNum);
            }
        }

        return update;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (COUNT_TOTAL_OPERATIONS)
            totalModelChangedCount++;
        likelihoodKnown = false;
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {

        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (COUNT_TOTAL_OPERATIONS)
            totalGetLogLikelihoodCount++;
        if (CompoundLikelihood.DEBUG_PARALLEL_EVALUATION) {
            System.err.println((likelihoodKnown ? "lazy" : "evaluate"));
        }
        if (!likelihoodKnown) {
            if (COUNT_TOTAL_OPERATIONS)
                totalCalculateLikelihoodCount++;
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        if (COUNT_TOTAL_OPERATIONS)
            totalMakeDirtyCount++;
        likelihoodKnown = false;
        updateAllNodes();
    }
    
    public boolean isLikelihoodKnown() {
    	return likelihoodKnown;
    }

    protected abstract double calculateLogLikelihood();

    public String getReport() {
        if (hasInitialized) {
            String rtnValue =  getClass().getName() + "(" + getLogLikelihood() + ")";
            if (COUNT_TOTAL_OPERATIONS)
             rtnValue += " total operations = " + totalOperationCount +
                         " matrix updates = " + totalMatrixUpdateCount + " model changes = " + totalModelChangedCount +
                         " make dirties = " + totalMakeDirtyCount +
                         " calculate likelihoods = " + totalCalculateLikelihoodCount +
                         " get likelihoods = " + totalGetLogLikelihoodCount +
                         " all rate updates = " + totalRateUpdateAllCount +
                         " partial rate updates = " + totalRateUpdateSingleCount;
            return rtnValue;
        } else {
            return getClass().getName() + "(uninitialized)";
        }
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the tree
     */
    protected TreeModel treeModel = null;

    /**
     * the number of nodes in the tree
     */
    protected final int nodeCount;

    protected final int tipCount;
    protected final int internalNodeCount;


    protected int[][] operations;
    protected int operationListCount;
    protected int[] operationCount;

    protected int[] branchUpdateIndices;
    protected double[] branchLengths;
    protected int branchUpdateCount;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    protected final BufferIndexHelper partialBufferHelper;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    protected boolean hasInitialized = false;

    protected int totalOperationCount = 0;
    protected int totalMatrixUpdateCount = 0;
    protected int totalGetLogLikelihoodCount = 0;
    protected int totalModelChangedCount = 0;
    protected int totalMakeDirtyCount = 0;
    protected int totalCalculateLikelihoodCount = 0;
    protected int totalRateUpdateAllCount = 0;
    protected int totalRateUpdateSingleCount = 0;

}