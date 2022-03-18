/*
 * BastaLikelihood.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.basta;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Guy Baele
 * @author Marc A. Suchard
 * @version $Id$
 */

public final class BastaLikelihood extends AbstractModelLikelihood implements
        TreeTraitProvider, Citable, Profileable, Reportable {

    private static final boolean COUNT_TOTAL_OPERATIONS = true;

    private static final boolean TEST = false;

    public BastaLikelihood(String name,
                           BastaLikelihoodDelegate likelihoodDelegate,
                           Tree treeModel,
                           BranchRateModel branchRateModel,
                           int numberSubIntervals) {

        super(name);

        assert likelihoodDelegate != null;
        assert treeModel != null;
        assert branchRateModel != null;

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing BastaLikelihood");

        this.likelihoodDelegate = likelihoodDelegate;
        addModel(likelihoodDelegate);

        this.treeModel = treeModel;
        isTreeRandom = (treeModel instanceof AbstractModel) && ((AbstractModel) treeModel).isVariable();
        if (isTreeRandom) {
            addModel((AbstractModel)treeModel);
        }

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);

        treeIntervals = new BigFastTreeIntervals((TreeModel)treeModel);
        treeTraversalDelegate = new CoalescentIntervalTraversal(treeModel, treeIntervals, branchRateModel, numberSubIntervals);

        addModel(treeIntervals);

        likelihoodKnown = false;
        hasInitialized = true;
    }

    @Override
    public final Model getModel() {
        return this;
    }

    @Override @SuppressWarnings("Duplicates")
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

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // TODO
    }
    @Override @SuppressWarnings("Duplicates")
    protected final void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {

                final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

                if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");

                if (treeChangedEvent.isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNode(((TreeChangedEvent) object).getNode());
                } else if (treeChangedEvent.isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
                    updateAllNodes();
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
        } else {

            assert false : "Unknown componentChangedEvent";
        }

        if (COUNT_TOTAL_OPERATIONS)
            totalModelChangedCount++;

        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected final void storeState() {

        assert (likelihoodKnown) : "the likelihood should always be known at this point in the cycle";
        storedLogLikelihood = logLikelihood;

        if (TEST) treeIntervals.storeModelState();
    }

    @Override
    protected final void restoreState() {

        // restore the likelihood and flag it as known
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;

        if (TEST) treeIntervals.restoreModelState();
    }

    @Override
    protected void acceptState() {
        if (TEST) treeIntervals.acceptModelState();
    } // nothing to do

    /**
     * Calculate the log likelihood of the data for the current tree.
     *
     * @return the log likelihood.
     */
    private double calculateLogLikelihood() {

        double logL;

//        likelihoodDelegate.makeDirty();

        treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

        final List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<ProcessOnCoalescentIntervalDelegate.OtherOperation> nodeOperations =
                treeTraversalDelegate.getOtherOperations();

        if (COUNT_TOTAL_OPERATIONS) {
            totalMatrixUpdateCount += branchOperations.size();
            totalOperationCount += nodeOperations.size();
        }

        final NodeRef root = treeModel.getRoot();
        logL = likelihoodDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());

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
//    protected void updateNodeAndChildren(NodeRef node) {
//        if (COUNT_TOTAL_OPERATIONS)
//            totalRateUpdateSingleCount += 1 + treeModel.getChildCount(node);
//
//        treeTraversalDelegate.updateNodeAndChildren(node);
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
        return null;
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

    @Override
    public Citation.Category getCategory() { return Citation.Category.FRAMEWORK; }

    @Override
    public String getDescription() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getDescription();
        } else {
            return null;
        }
    }

    @Override
    public List<Citation> getCitations() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getCitations();
        } else {
            return new ArrayList<>();
        }
    }

    // **************************************************************
    // INSTANCE PROFILEABLE
    // **************************************************************

    @Override
    public long getTotalCalculationCount() {
        return likelihoodDelegate.getTotalCalculationCount();
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * The data likelihood delegate
     */
    private final
    BastaLikelihoodDelegate likelihoodDelegate;

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

    private final CoalescentIntervalTraversal treeTraversalDelegate;
    private final BigFastTreeIntervals treeIntervals;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown;

    private boolean hasInitialized;

    private final boolean isTreeRandom;

    private int totalOperationCount = 0;
    private int totalMatrixUpdateCount = 0;
    private int totalGetLogLikelihoodCount = 0;
    private int totalModelChangedCount = 0;
    private int totalMakeDirtyCount = 0;
    private int totalCalculateLikelihoodCount = 0;
    private int totalRateUpdateAllCount = 0;
    private int totalRateUpdateSingleCount = 0;
    private int totalPostOrderStatistics = 0;
    private int totalCalculatePostOrderStatistics = 0;
}
