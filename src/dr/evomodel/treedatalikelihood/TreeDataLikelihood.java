/*
 * BeagleTreeLikelihood.java
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
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.*;
import java.util.logging.Logger;

/**
 * TreeDataLikelihood - uses plugin delegates to compute the likelihood of some data given a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

public final class TreeDataLikelihood extends AbstractModelLikelihood implements TreeTraitProvider, Reportable {

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    private static final long MAX_UNDERFLOWS_BEFORE_ERROR = 100;

    public enum TraversalType {
        POST_ORDER,
        REVERSE_LEVEL_ORDER
    };

    public TreeDataLikelihood(DataLikelihoodDelegate likelihoodDelegate,
                              TreeModel treeModel,
                              BranchRateModel branchRateModel) {
        this(likelihoodDelegate, treeModel, branchRateModel, null);
    }

    public TreeDataLikelihood(DataLikelihoodDelegate likelihoodDelegate,
                              TreeModel treeModel,
                              BranchRateModel branchRateModel,
                              DataSimulationDelegate simulationDelegate) {

        super("TreeDataLikelihood");  // change this to use a const once the parser exists

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing TreeDataLikelihood");

        this.likelihoodDelegate = likelihoodDelegate;
        addModel(likelihoodDelegate);
        likelihoodDelegate.setCallback(this);

        this.simulationDelegate = simulationDelegate;
        if (simulationDelegate != null) {
            addModel(simulationDelegate);
            simulationDelegate.setCallback(this);
            treeTraits.addTraits(simulationDelegate.getTreeTraits());
        }

        this.treeModel = treeModel;
        addModel(treeModel);

        updateNode = new boolean[treeModel.getNodeCount()];
        for (int i = 0; i < updateNode.length; i++) {
            updateNode[i] = true;
        }

        likelihoodKnown = false;

        if (branchRateModel != null) {
            this.branchRateModel = branchRateModel;
            logger.info("  Branch rate model used: " + branchRateModel.getModelName());
        } else {
            this.branchRateModel = new DefaultBranchRateModel();
        }
        addModel(this.branchRateModel);

        hasInitialized = true;
    }

    public DataLikelihoodDelegate getDataLikelihoodDelegate() {
        return likelihoodDelegate;
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    @Override
    public final Model getModel() {
        return this;
    }

    @Override
    public final double getLogLikelihood() {
        if (COUNT_TOTAL_OPERATIONS)
            totalGetLogLikelihoodCount++;

        if (!likelihoodKnown) {
            if (COUNT_TOTAL_OPERATIONS)
                totalCalculateLikelihoodCount++;

            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    @Override
    public final void makeDirty() {
        if (COUNT_TOTAL_OPERATIONS)
            totalMakeDirtyCount++;

        likelihoodKnown = false;
        likelihoodDelegate.makeDirty();
        updateAllNodes();
    }

    public final boolean isLikelihoodKnown() {
        return likelihoodKnown;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    @Override
    protected final void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {

                if (((TreeModel.TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeModel.TreeChangedEvent) object).getNode());

                } else if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }
        } else if (model == likelihoodDelegate) {

            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else if (model == branchRateModel) {

            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else if (model == simulationDelegate) {

            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else {

            assert false: "Unknown componentChangedEvent";
        }

        if (COUNT_TOTAL_OPERATIONS)
            totalModelChangedCount++;

        likelihoodKnown = false;

        fireModelChanged();
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    @Override
    protected final void storeState() {

        assert(likelihoodKnown) : "the likelihood should always be known at this point in the cycle";

        storedLogLikelihood = logLikelihood;

    }

    @Override
    protected final void restoreState() {

        // restore the likelihood and flag it as known
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
    }

    @Override
    protected void acceptState() {
    } // nothing to do


    /**
     * Simulate process along the current tree.
     */
    protected final void simulateProcess() {

    }


    /**
     * Calculate the log likelihood of the data for the current tree.
     *
     * @return the log likelihood.
     */
    private final double calculateLogLikelihood() {

        double logL = Double.NEGATIVE_INFINITY;
        boolean done = false;
        long underflowCount = 0;

        do {
            dispatchTreeTraversalCollectBranchAndNodeOperations(likelihoodDelegate);

            final NodeRef root = treeModel.getRoot();

            try {
                logL = likelihoodDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());

                done = true;
            } catch (DataLikelihoodDelegate.LikelihoodUnderflowException e) {

                // if there is an underflow, assume delegate will attempt to rescale
                // so flag all nodes to update and return to try again.
                updateAllNodes();
                underflowCount ++;
            }

        } while (!done && underflowCount < MAX_UNDERFLOWS_BEFORE_ERROR);

        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < updateNode.length; i++) {
            updateNode[i] = false;
        }

        return logL;
    }

    private void dispatchTreeTraversalCollectBranchAndNodeOperations(ProcessOnTreeDelegate delegate) {
        branchOperations.clear();
        nodeOperations.clear();

        switch (delegate.getOptimalTraversalType()) {

            case POST_ORDER:
                traversePostOrder(treeModel);
                break;
            case REVERSE_LEVEL_ORDER:
                traverseReverseLevelOrder(treeModel);
                break;
            default:
                assert false : "Unknown traversal type";
        }

        if (COUNT_TOTAL_OPERATIONS) {
            totalMatrixUpdateCount += branchOperations.size();
            totalOperationCount += nodeOperations.size();
        }
    }

    /**
     * Traverse the tree in post order.
     *
     * @param tree           tree
     * @return boolean
     */
    private void traversePostOrder(Tree tree) {
        traversePostOrder(tree, tree.getRoot());
    }

    /**
     * Traverse the tree in post order.
     *
     * @param tree           tree
     * @param node           node
     * @return boolean
     */
    private boolean traversePostOrder(Tree tree, NodeRef node) {

        boolean update = false;

        int nodeNum = node.getNumber();

        // First update the transition probability matrix(ices) for this branch
        if (tree.getParent(node) != null && updateNode[nodeNum]) {
            // @todo - at the moment a matrix is updated even if a branch length doesn't change

            addBranchUpdateOperation(tree, node);

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traversePostOrder(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traversePostOrder(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                nodeOperations.add(new DataLikelihoodDelegate.NodeOperation(nodeNum, child1.getNumber(), child2.getNumber()));

                update = true;

            }
        }

        return update;

    }

    /**
     * Traverse the tree in reverse level order.
     *
     * @param tree           tree
     */
    private void traverseReverseLevelOrder(final Tree tree) {

        // create a map of all the operations at each particular level
        Map<Integer, List<DataLikelihoodDelegate.NodeOperation>> operationMap =
                new HashMap<Integer, List<DataLikelihoodDelegate.NodeOperation>>();

        traverseLevelOrder(tree, tree.getRoot(), 0, operationMap);

        // get the levels as keys in reverse order (they are currently largest towards
        // the tips) and add the operations to the nodeOperation array.
        List<Integer> keyList = new ArrayList<Integer>(operationMap.keySet());
        Collections.sort(keyList, Collections.reverseOrder());

        for (Integer key : keyList) {
            List<DataLikelihoodDelegate.NodeOperation> opList = operationMap.get(key);
            for (DataLikelihoodDelegate.NodeOperation op : opList) {
                nodeOperations.add(op);
            }
        }
    }

    /**
     * Traverse the tree in level order.
     *
     * @param tree           tree
     * @param node          node
     * @return boolean
     */
    private boolean traverseLevelOrder(final Tree tree, final NodeRef node,
                                              final int level,
                                              Map<Integer, List<DataLikelihoodDelegate.NodeOperation>> operationMap) {
        boolean update = false;

        int nodeNum = node.getNumber();

        // First update the transition probability matrix(ices) for this branch
        if (tree.getParent(node) != null && updateNode[nodeNum]) {
            // @todo - at the moment a matrix is updated even if a branch length doesn't change

            addBranchUpdateOperation(tree, node);

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes incrementing the level (this will give
            // level order but we will reverse these later
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverseLevelOrder(tree, child1, level + 1, operationMap);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverseLevelOrder(tree, child2, level + 1, operationMap);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                List<DataLikelihoodDelegate.NodeOperation> ops = operationMap.get(level);
                if (ops == null) {
                    ops = new ArrayList<DataLikelihoodDelegate.NodeOperation>();
                    operationMap.put(level, ops);
                }
                ops.add(new DataLikelihoodDelegate.NodeOperation(nodeNum, child1.getNumber(), child2.getNumber()));

                update = true;

            }
        }

        return update;
    }

    /**
     * Add this node to the branchOperation list for updating of the transition probability matrix.
     *
     * @param tree           tree
     * @param node           node
     */
    private void addBranchUpdateOperation(Tree tree, NodeRef node) {
        final double branchRate;

        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }
        final double parentHeight = tree.getNodeHeight(tree.getParent(node));
        final double nodeHeight = tree.getNodeHeight(node);

        // Get the operational time of the branch
        final double branchLength = branchRate * (parentHeight - nodeHeight);

        assert branchLength > 0.0 : "Negative branch length: " + branchLength + " for node " +
                node.getNumber() + (tree.isExternal(node) ?
                " (" + tree.getNodeTaxon(node).getId() + ")": "");

        branchOperations.add(new DataLikelihoodDelegate.BranchOperation(node.getNumber(), branchLength));
    }

    /**
     * Set update flag for a node only
     */
    protected void updateNode(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

        updateNode[node.getNumber()] = true;
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            if (COUNT_TOTAL_OPERATIONS)
                totalRateUpdateSingleCount++;

            NodeRef child = treeModel.getChild(node, i);
            updateNode[child.getNumber()] = true;
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and all its descendents
     */
    protected void updateNodeAndDescendents(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

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
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateAllCount++;

        for (int i = 0; i < updateNode.length; i++) {
            updateNode[i] = true;
        }
        likelihoodKnown = false;
    }

    // **************************************************************
    // Reportable IMPLEMENTATION
    // **************************************************************

    @Override
    public String getReport() {
        if (hasInitialized) {
            String rtnValue =  getClass().getName() + "(" + getLogLikelihood() + ")";
            if (COUNT_TOTAL_OPERATIONS)
                rtnValue += "\n  total operations = " + totalOperationCount +
                        "\n  matrix updates = " + totalMatrixUpdateCount +
                        "\n  model changes = " + totalModelChangedCount +
                        "\n  make dirties = " + totalMakeDirtyCount +
                        "\n  calculate likelihoods = " + totalCalculateLikelihoodCount +
                        "\n  get likelihoods = " + totalGetLogLikelihoodCount +
                        "\n  all rate updates = " + totalRateUpdateAllCount +
                        "\n  partial rate updates = " + totalRateUpdateSingleCount;
            return rtnValue;
        } else {
            return getClass().getName() + "(uninitialized)";
        }
    }

    // **************************************************************
    // TreeTrait IMPLEMENTATION
    // **************************************************************

    /**
     * Returns an array of all the available traits
     *
     * @return the array
     */
    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    /**
     * Returns a trait that is stored using a specific key. This will often be the same
     * as the 'name' of the trait but may not be depending on the application.
     *
     * @param key a unique key
     * @return the trait
     */
    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * The data likelihood delegate
     */
    private final DataLikelihoodDelegate likelihoodDelegate;

    /**
     * the tree model
     */
    private final TreeModel treeModel;

    /**
     * the branch rate model
     */
    private final BranchRateModel branchRateModel;

    /**
     * data simulation delegates
     */
    private final DataSimulationDelegate simulationDelegate;

    /**
     * TreeTrait helper
     */
    private Helper treeTraits = new Helper();

    /**
     * Flags to specify which nodes are to be updated
     */
    private boolean[] updateNode;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;

    private boolean hasInitialized = false;

    private int totalOperationCount = 0;
    private int totalMatrixUpdateCount = 0;
    private int totalGetLogLikelihoodCount = 0;
    private int totalModelChangedCount = 0;
    private int totalMakeDirtyCount = 0;
    private int totalCalculateLikelihoodCount = 0;
    private int totalRateUpdateAllCount = 0;
    private int totalRateUpdateSingleCount = 0;

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private List<DataLikelihoodDelegate.BranchOperation> branchOperations = new ArrayList<DataLikelihoodDelegate.BranchOperation>();
    private List<DataLikelihoodDelegate.NodeOperation> nodeOperations = new ArrayList<DataLikelihoodDelegate.NodeOperation>();

}//END: class