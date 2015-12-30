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
import dr.app.beagle.evomodel.treelikelihood.EvolutionaryProcessDelegate;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.Reportable;

import java.util.Map;
import java.util.Set;

/**
 * AbstractTreeLikelihood - a base class for likelihood calculators of sites on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id: AbstractTreeLikelihood.java,v 1.16 2005/06/07 16:27:39 alexei Exp $
 */

public abstract class NewAbstractLikelihoodOnTree extends AbstractModelLikelihood implements Reportable {

    protected static final boolean COUNT_TOTAL_OPERATIONS = false;

    public NewAbstractLikelihoodOnTree(String name, TreeModel treeModel, Map<Set<String>, Parameter> partialsRestrictions) {

        super(name);

        this.treeModel = treeModel;
        addModel(treeModel);

        nodeCount = treeModel.getNodeCount();
        tipCount = treeModel.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        // one partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

        likelihoodKnown = false;

        // Partials restrictions
        this.partialsRestrictions = partialsRestrictions;
//            hasRestrictedPartials = (partialsRestrictions != null);
        if (hasRestrictedPartials) {
            numRestrictedPartials = partialsRestrictions.size();
            updateRestrictedNodePartials = true;
            partialsMap = new Parameter[treeModel.getNodeCount()];
//            partials = new double[stateCount * patternCount * categoryCount];
        } else {
            numRestrictedPartials = 0;
            updateRestrictedNodePartials = false;
        }
    }

    public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        /* No need to rescale partials */
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
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

    protected boolean updateBranchSpecificEvolutionaryProcess(final Tree tree, final int nodeNum, final NodeRef parent,
                                                              final NodeRef node, final boolean flip) {
        return false;
    }

    protected void handlePartialsScaling(final int[] operations, final int nodeNum, final int x) {
        // Do nothing
    }

//    protected void handleEvolutionaryProcess(final int[] operations, final NodeRef child1, final NodeRef child2, final int x) {
//        // Do nothing
//    }

    protected void handleRestrictedPartials(final int nodeNum) {
        // Do nothing
    }

    protected void setEvolutionaryProcessDelegate(final EvolutionaryProcessDelegate delegate) {
        this.evolutionaryProcessDelegate = delegate;
    }

    protected void prepareStorage() {
        if (branchUpdateIndices == null) {
            branchUpdateIndices = new int[nodeCount];
            branchLengths = new double[nodeCount];
//            scaleBufferIndices = new int[internalNodeCount];
//            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[numRestrictedPartials + 1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[numRestrictedPartials + 1];
        }
    }

    protected void prepareScaling() {
        // Do nothing
    }

    protected void prepareTips() {
        // Do nothing
    }

    protected void updateRootInformation() {
        // Do nothing
    }

    protected int accumulateScaleFactors() {
        return Beagle.NONE;
    }

    protected void updateSiteModelAction() {
        // Do nothing
    }

    protected double computedAscertainedLogLikelihood() {
        throw new RuntimeException("Not an ascertained model");
    }

    protected boolean doRescalingNow(boolean firstAttempt) {
        return false;
    }

    /**
      * Calculate the log likelihood of the current state.
      *
      * @return the log likelihood.
      */
    protected double calculateLogLikelihood() {

        prepareStorage();

        prepareScaling();

        prepareTips();

        branchUpdateCount = 0;
        operationListCount = 0;

        if (hasRestrictedPartials) {
            for (int i = 0; i <= numRestrictedPartials; i++) {
                operationCount[i] = 0;
            }
        } else {
            operationCount[0] = 0;
        }

        final NodeRef root = treeModel.getRoot();
        postOrderTraverse(treeModel, root, null, true);

        if (updateSubstitutionModel) { // TODO More efficient to update only the substitution model that changed, instead of all
            evolutionaryProcessDelegate.updateSubstitutionModels(beagle);

            // we are currently assuming a no-category model...
        }

        if (updateSiteModel) {
            updateSiteModelAction();
        }

        if (branchUpdateCount > 0) {
            evolutionaryProcessDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount);
        }

        if (COUNT_TOTAL_OPERATIONS) {
            totalMatrixUpdateCount += branchUpdateCount;

            for (int i = 0; i <= numRestrictedPartials; i++) {
                totalOperationCount += operationCount[i];
            }
        }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            if (hasRestrictedPartials) {
                for (int i = 0; i <= numRestrictedPartials; i++) {
                    beagle.updatePartials(operations[i], operationCount[i], Beagle.NONE);
                    if (i < numRestrictedPartials) {
                        //                        restrictNodePartials(restrictedIndices[i]);
                    }
                }
            } else {
                beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);
            }

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNumber());

            int cumulateScaleBufferIndex = accumulateScaleFactors();

            updateRootInformation();

            double[] sumLogLikelihoods = new double[1];

            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            logL = sumLogLikelihoods[0];

            if (ascertainedSitePatterns) {
                logL = computedAscertainedLogLikelihood();
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (doRescalingNow(firstRescaleAttempt)) {

                    branchUpdateCount = 0;

                    if (hasRestrictedPartials) {
                        for (int i = 0; i <= numRestrictedPartials; i++) {
                            operationCount[i] = 0;
                        }
                    } else {
                        operationCount[0] = 0;
                    }

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    postOrderTraverse(treeModel, root, null, false);

                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        // If these are needed...
        //beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
        //********************************************************************

        return logL;
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
    final protected boolean postOrderTraverse(final Tree tree, final NodeRef node, final int[] operatorNumber, final boolean flip) {

        boolean update = false;

        final int nodeNum = node.getNumber();

        final NodeRef parent = tree.getParent(node);

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the evolutionary process for this branch
        if (parent != null && updateNode[nodeNum]) {
            update = updateBranchSpecificEvolutionaryProcess(tree, nodeNum, parent, node, flip);
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            final NodeRef child1 = tree.getChild(node, 0);
            final int[] op1 = {-1};
            final boolean update1 = postOrderTraverse(tree, child1, op1, flip);

            final NodeRef child2 = tree.getChild(node, 1);
            final int[] op2 = {-1};
            final boolean update2 = postOrderTraverse(tree, child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                final int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;
                final int[] operations = this.operations[operationListCount];

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                handlePartialsScaling(operations, nodeNum, x);

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNumber()); // source node 1
                operations[x + 4] = evolutionaryProcessDelegate.getMatrixIndex(child1.getNumber()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNumber()); // source node 2
                operations[x + 6] = evolutionaryProcessDelegate.getMatrixIndex(child2.getNumber()); // source matrix 2

                handleRestrictedPartials(nodeNum);

                ++operationCount[operationListCount];
                update = true;
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

        fireModelChanged();

        if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {

                if (((TreeModel.TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeModel.TreeChangedEvent) object).getNode());
                    updateRestrictedNodePartials = true;

                } else if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
                    //                    System.err.println("Full tree update event - these events currently aren't used\n" +
                    //                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                    updateRestrictedNodePartials = true;
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }
        }

        if (COUNT_TOTAL_OPERATIONS) {
            totalModelChangedCount++;
        }

        likelihoodKnown = false;
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {
        partialBufferHelper.storeState();
        evolutionaryProcessDelegate.storeState();

        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore

        partialBufferHelper.restoreState();
        evolutionaryProcessDelegate.restoreState();

        updateRestrictedNodePartials = true;

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
        if (COUNT_TOTAL_OPERATIONS) {
            totalGetLogLikelihoodCount++;
        }

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
        if (COUNT_TOTAL_OPERATIONS) {
            totalMakeDirtyCount++;
        }

        likelihoodKnown = false;
        updateAllNodes();

        updateSiteModel = true;
        updateSubstitutionModel = true;
        updateRestrictedNodePartials = true;
    }
    
    public boolean isLikelihoodKnown() {
    	return likelihoodKnown;
    }

//    protected abstract double calculateLogLikelihood();

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

    protected static final boolean hasRestrictedPartials = false;
    protected final int numRestrictedPartials;
    protected final Map<Set<String>, Parameter> partialsRestrictions;
    protected Parameter[] partialsMap;
    protected double[] partials;
    protected boolean updateRestrictedNodePartials;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    protected final BufferIndexHelper partialBufferHelper;
    protected EvolutionaryProcessDelegate evolutionaryProcessDelegate;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    // Scale factor info
    protected boolean everUnderflowed = false;
//    protected boolean implementsScaling = false;
//    protected boolean recomputeScaleFactors = false;


    /**
     * the BEAGLE library instance
     */
    protected Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;

    /**
     * Flag to specify if site patterns are ascertained
     */
    protected boolean ascertainedSitePatterns = false;




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