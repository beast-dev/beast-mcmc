/*
 * PrefetchingTreeDataLikelihood.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
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

public final class PrefetchTreeDataLikelihood extends AbstractModelLikelihood implements TreeTraitProvider, Reportable, PrefetchableLikelihood {

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    private static final long MAX_UNDERFLOWS_BEFORE_ERROR = 100;
    private static final boolean PREFETCH_DEBUG = false;

    public PrefetchTreeDataLikelihood(DataLikelihoodDelegate likelihoodDelegate,
                                      TreeModel treeModel,
                                      BranchRateModel branchRateModel) {

        super("TreeDataLikelihood");  // change this to use a const once the parser exists

        assert likelihoodDelegate != null;
        assert treeModel != null;
        assert branchRateModel != null;

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing TreeDataLikelihood");

        if (!(likelihoodDelegate instanceof PrefetchMultiPartitionDataLikelihoodDelegate)) {
            throw new RuntimeException("DataLikelihoodDelegate in PrefetchTreeDataLikelihood " +
                    "should be PrefetchMultiPartitionDataLikelihoodDelegate");
        }
        this.likelihoodDelegate = (PrefetchMultiPartitionDataLikelihoodDelegate)likelihoodDelegate;
        int prefetchCount = this.likelihoodDelegate.getPrefetchCount();

        addModel(likelihoodDelegate);
        likelihoodDelegate.setCallback(null);

        this.treeModel = treeModel;
        isTreeRandom = treeModel.isTreeRandom();
        if (isTreeRandom) {
            addModel(treeModel);
        }

        likelihoodKnown = false;

        this.branchRateModel = branchRateModel;
        if (!(branchRateModel instanceof DefaultBranchRateModel)) {
            logger.info("  Branch rate model used: " + branchRateModel.getModelName());
        }
        addModel(this.branchRateModel);

        this.prefetchCount = prefetchCount;

        treeTraversalDelegates = new LikelihoodTreeTraversal[prefetchCount];
        for (int i = 0; i < prefetchCount; i++) {
            treeTraversalDelegates[i] = new LikelihoodTreeTraversal(treeModel, branchRateModel,
                    likelihoodDelegate.getOptimalTraversalType());
        }

        prefetchedLogLikelihoods = new double[prefetchCount];

        isPrefetching = false;
        currentPrefetch = -1;

        hasInitialized = true;
    }

    public Tree getTree() {
        return treeModel;
    }

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    // **************************************************************
    // PrefetchableLikelihood IMPLEMENTATION
    // **************************************************************

    @Override
    public int getPrefetchCount() {
        return prefetchCount;
    }

    @Override
    public void setCurrentPrefetch(int prefetch) {
        likelihoodDelegate.setCurrentPrefetch(prefetch);

        isPrefetching = !(currentPrefetch < 0);
        this.currentPrefetch = prefetch;
        likelihoodKnown = false;
    }

    public void collectOperations() {
        treeTraversalDelegates[currentPrefetch].dispatchTreeTraversalCollectBranchAndNodeOperations();
    }

    public void prefetchLogLikelihoods() {
        // todo this is currently calling sequentially. Needs to be concurrent for prefetching to be meaningful
        for (int i = 0; i < prefetchCount; i++) {
            setCurrentPrefetch(i);
            prefetchedLogLikelihoods[i] = calculateLogLikelihood(treeTraversalDelegates[i]);

            // restore the state so the next set of updates can be applied
            // If doing concurrently then this would just be done once for all the partitions
            likelihoodDelegate.restoreState();
        }
    }

    @Override
    public void acceptPrefetch(int prefetch) {
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

            if (!isPrefetching) {
                logLikelihood = calculateLogLikelihood(treeTraversalDelegates[0]);
                if (PREFETCH_DEBUG) {
                    System.err.println("PTDL calculated likelihood: " + logLikelihood);
                }
            } else {
                logLikelihood = prefetchedLogLikelihoods[currentPrefetch];
                if (PREFETCH_DEBUG) {
                    System.err.println("PTDL returning pre-calculated likelihood: " + prefetchedLogLikelihoods[currentPrefetch]);
                }
            }
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

                if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");

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
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                }
            }
        } else if (model == likelihoodDelegate) {

            assert !isPrefetching;

            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else if (model == branchRateModel) {

            assert !isPrefetching;

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
        currentPrefetch = -1;
        isPrefetching = false;
    } // nothing to do

    /**
     * Calculate the log likelihood of the data for the current tree.
     *
     * @return the log likelihood.
     * @param treeTraversalDelegate
     */
    private final double calculateLogLikelihood(LikelihoodTreeTraversal treeTraversalDelegate) {

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
        if (PREFETCH_DEBUG) {
            System.err.println("TDL setting all nodes to updated - " + (isPrefetching ? currentPrefetch + 1 : "no prefetch"));
        }
        treeTraversalDelegates[(isPrefetching ? currentPrefetch : 0)].setAllNodesUpdated();
    }

    /**
     * Set update flag for a node only
     */
    protected void updateNode(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

        if (PREFETCH_DEBUG) {
            System.err.println("TDL setting update for node " + node.getNumber() + " - " + (isPrefetching ? currentPrefetch + 1 : "no prefetch"));
        }

        treeTraversalDelegates[(isPrefetching ? currentPrefetch : 0)].updateNode(node);
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount += 1 + treeModel.getChildCount(node);

        if (PREFETCH_DEBUG) {
            System.err.println("TDL setting update for node " + node.getNumber() + " and children - " + (isPrefetching ? currentPrefetch + 1 : "no prefetch"));
        }

        treeTraversalDelegates[(isPrefetching ? currentPrefetch : 0)].updateNodeAndChildren(node);
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and all its descendents
     */
    protected void updateNodeAndDescendents(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateSingleCount++;

        treeTraversalDelegates[(isPrefetching ? currentPrefetch : 0)].updateNodeAndDescendents(node);
        likelihoodKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        if (COUNT_TOTAL_OPERATIONS)
            totalRateUpdateAllCount++;

        if (PREFETCH_DEBUG) {
            System.err.println("TDL setting updating all nodes - " + (isPrefetching ? currentPrefetch + 1 : "no prefetch"));
        }

        treeTraversalDelegates[(isPrefetching ? currentPrefetch : 0)].updateAllNodes();

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
                System.err.println(delegateString);
            }

            sb.append(getClass().getName() + "(" + getLogLikelihood() + ")");

            if (COUNT_TOTAL_OPERATIONS)
                sb.append("\n  total operations = " + totalOperationCount +
                        "\n  matrix updates = " + totalMatrixUpdateCount +
                        "\n  model changes = " + totalModelChangedCount +
                        "\n  make dirties = " + totalMakeDirtyCount +
                        "\n  calculate likelihoods = " + totalCalculateLikelihoodCount +
                        "\n  get likelihoods = " + totalGetLogLikelihoodCount +
                        "\n  all rate updates = " + totalRateUpdateAllCount +
                        "\n  partial rate updates = " + totalRateUpdateSingleCount);

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
    private final PrefetchMultiPartitionDataLikelihoodDelegate likelihoodDelegate;

    /**
     * the tree model
     */
    private final TreeModel treeModel;

    /**
     * the branch rate model
     */
    private final BranchRateModel branchRateModel;

    /**
     * TreeTrait helper
     */
    private final Helper treeTraits = new Helper();

    private final LikelihoodTreeTraversal[] treeTraversalDelegates;

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

    private int currentPrefetch = 0;
    private boolean isPrefetching = false;
    private final int prefetchCount;
    private double[] prefetchedLogLikelihoods;
}