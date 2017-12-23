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

import dr.evolution.tree.*;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.*;
import dr.xml.Reportable;

import java.util.List;
import java.util.logging.Logger;

/**
 * TreeDataLikelihood - uses plugin delegates to compute the likelihood of some data given a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

public final class TreeDataLikelihood extends AbstractModelLikelihood implements TreeTraitProvider, Reportable {

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
    private static final long MAX_UNDERFLOWS_BEFORE_ERROR = 100;

    public TreeDataLikelihood(DataLikelihoodDelegate likelihoodDelegate,
                              Tree treeModel,
                              BranchRateModel branchRateModel) {

        super("TreeDataLikelihood");  // change this to use a const once the parser exists

        assert likelihoodDelegate != null;
        assert treeModel != null;
        assert branchRateModel != null;

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing TreeDataLikelihood");

        this.likelihoodDelegate = likelihoodDelegate;
        addModel(likelihoodDelegate);
        likelihoodDelegate.setCallback(this);

        this.treeModel = treeModel;
        isTreeRandom = (treeModel instanceof AbstractModel) && ((AbstractModel) treeModel).isVariable();
        if (isTreeRandom) {
            addModel(((AbstractModel)treeModel));
        }

        likelihoodKnown = false;

        this.branchRateModel = branchRateModel;
        if (!(branchRateModel instanceof DefaultBranchRateModel)) {
            logger.info("  Branch rate model used: " + branchRateModel.getModelName());
        }
        addModel(this.branchRateModel);

        treeTraversalDelegate = new LikelihoodTreeTraversal(treeModel, branchRateModel,
                    likelihoodDelegate.getOptimalTraversalType()
        );

        rateRescalingScheme = likelihoodDelegate.getRateRescalingScheme();

        hasInitialized = true;
    }

    public final Tree getTree() {
        return treeModel;
    }

    public final BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    public final DataLikelihoodDelegate getDataLikelihoodDelegate() {
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
            if (COUNT_TOTAL_OPERATIONS) {
                totalCalculateLikelihoodCount++;
            }

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
            if (object instanceof TreeChangedEvent) {

                final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

                if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");

                if (treeChangedEvent.isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    if (rateRescalingScheme == RateRescalingScheme.NONE ||             // The usual behaviour
                            (rateRescalingScheme == RateRescalingScheme.TREE_HEIGHT && // Or, tree-height scaling
                                    !treeModel.isRoot(treeChangedEvent.getNode()))) {  // and not the root
                        updateNodeAndChildren(((TreeChangedEvent) object).getNode());
                    } else {
                        // Tree-length scaling or a change of the root
                        updateAllNodes();
                    }
                } else if (treeChangedEvent.isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
                    updateAllNodes();
                }
//                else {
//                    // Other event types are ignored (probably trait changes).
//                }
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
        } else {

            assert false : "Unknown componentChangedEvent";
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

        assert (likelihoodKnown) : "the likelihood should always be known at this point in the cycle";

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
     * Calculate the log likelihood of the data for the current tree.
     *
     * @return the log likelihood.
     */
    private double calculateLogLikelihood() {

        double logL = Double.NEGATIVE_INFINITY;
        boolean done = false;
        long underflowCount = 0;

        do {
            treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

            final List<DataLikelihoodDelegate.BranchOperation> branchOperations = treeTraversalDelegate.getBranchOperations();
            final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = treeTraversalDelegate.getNodeOperations();

            if (COUNT_TOTAL_OPERATIONS) {
                totalMatrixUpdateCount += branchOperations.size();
                totalOperationCount += nodeOperations.size();
            }

            final NodeRef root = treeModel.getRoot();

            try {
                logL = likelihoodDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());

                done = true;
            } catch (DataLikelihoodDelegate.LikelihoodException e) {

                // if there is an underflow, assume delegate will attempt to rescale
                // so flag all nodes to update and return to try again.
                updateAllNodes();
                underflowCount++;
            }

        } while (!done && underflowCount < MAX_UNDERFLOWS_BEFORE_ERROR);

        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        setAllNodesUpdated();

        return logL;
    }

    private void setAllNodesUpdated() {
        treeTraversalDelegate.setAllNodesUpdated();
    }

    /**
     * Set update flag for a node only
     */
    protected void updateNode(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

        treeTraversalDelegate.updateNode(node);
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount += 1 + treeModel.getChildCount(node);

        treeTraversalDelegate.updateNodeAndChildren(node);
        likelihoodKnown = false;
    }

//    /**
//     * Set update flag for a node and all its descendents
//     */
//    protected void updateNodeAndDescendents(NodeRef node) {
//        if (COUNT_TOTAL_OPERATIONS)
//            totalRateUpdateSingleCount++;
//
//        treeTraversalDelegate.updateNodeAndDescendents(node);
//        likelihoodKnown = false;
//    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateAllCount++;

        treeTraversalDelegate.updateAllNodes();
        likelihoodKnown = false;
    }

    // **************************************************************
    // Reportable IMPLEMENTATION
    // **************************************************************

    @Override
    public String getReport() {
        if (hasInitialized) {
            StringBuilder sb = new StringBuilder();

            String delegateString = likelihoodDelegate.getReport();
            if (delegateString != null) {
                sb.append(delegateString);
            }

            sb.append(getClass().getName()).append("(").append(getLogLikelihood()).append(")");

            if (COUNT_TOTAL_OPERATIONS)
                sb.append("\n  total operations = ").append(totalOperationCount).append(
                          "\n  matrix updates = ").append(totalMatrixUpdateCount).append(
                          "\n  model changes = ").append(totalModelChangedCount).append(
                          "\n  make dirties = ").append(totalMakeDirtyCount).append(
                          "\n  calculate likelihoods = ").append(totalCalculateLikelihoodCount).append(
                          "\n  get likelihoods = ").append(totalGetLogLikelihoodCount).append(
                          "\n  all rate updates = ").append(totalRateUpdateAllCount).append(
                          "\n  partial rate updates = ").append(totalRateUpdateSingleCount);

            return sb.toString();
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
    // Decorate with TreeTraitProviders
    // **************************************************************

    public void addTrait(TreeTrait trait) {
        treeTraits.addTrait(trait);
    }

    public void addTraits(TreeTrait[] traits) {
        treeTraits.addTraits(traits);
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
    private final Tree treeModel;

    /**
     * the branch rate model
     */
    private final BranchRateModel branchRateModel;

    /**
     * TreeTrait helper
     */
    private final Helper treeTraits = new Helper();

    private final LikelihoodTreeTraversal treeTraversalDelegate;

    private final RateRescalingScheme rateRescalingScheme;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;

    private boolean hasInitialized = false;

    private final boolean isTreeRandom;

    private int totalOperationCount = 0;
    private int totalMatrixUpdateCount = 0;
    private int totalGetLogLikelihoodCount = 0;
    private int totalModelChangedCount = 0;
    private int totalMakeDirtyCount = 0;
    private int totalCalculateLikelihoodCount = 0;
    private int totalRateUpdateAllCount = 0;
    private int totalRateUpdateSingleCount = 0;
}