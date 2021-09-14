package dr.evomodel.treelikelihood.thorneytreelikelihood;
/*
 * AbstractTreeLikelihood.java
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


import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * ApproximatePoissonTreeLikelihood - a tree likelihood which uses an ML tree as expected number
 * of substitutions and assumes that these are drawn from a gamma (with mean == variance).
 *
 * This is similar to work by Jeff Thorne and colleagues:
 *
 *
 * And more recently Didelot et al (2018) BioRxiv
 *
 * @author Andrew Rambaut
 */

public class ThorneyTreeLikelihood extends AbstractModelLikelihood implements Reportable {

    public ThorneyTreeLikelihood(String name, TreeModel treeModel, BranchLengthProvider branchLengthProvider, ThorneyBranchLengthLikelihoodDelegate thorneyBranchLengthLikelihoodDelegate) {

        super(name);

        this.treeModel = treeModel;
        addModel(treeModel);
        this.thorneyBranchLengthLikelihoodDelegate = thorneyBranchLengthLikelihoodDelegate;
        this.branchLengthProvider = branchLengthProvider;
         assert thorneyBranchLengthLikelihoodDelegate instanceof Model;
         addModel((Model) thorneyBranchLengthLikelihoodDelegate);



        updateNode = new boolean[treeModel.getNodeCount()];
        Arrays.fill(updateNode, true);

        branchLogL = new double[treeModel.getNodeCount()];
        storedBranchLogL = new double[treeModel.getNodeCount()];
        likelihoodKnown = false;

        cachedRoot = treeModel.getRoot().getNumber();
        cachedRootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        cachedRootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();
    }



    /**
     * Set update flag for node and remove it's old contribution to the likelihood.
     * Also handle the root and children so that the 1 branch between children is marked as updated.
     * @param node
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        NodeRef parent = treeModel.getParent(node);
        if (parent != null && !updateNode[parent.getNumber()]) {
            updateNode(parent);
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {

        updateNode(node);

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNode(child);
        }

//        likelihoodKnown = false;
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
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            updateNode[i] = true;
        }
        likelihoodKnown = false;
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

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        fireModelChanged();

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    NodeRef node= ((TreeChangedEvent) object).getNode();
                    updateNodeAndChildren(node);
                    if(treeModel.getRoot()== node){
                        // If the root is changed then the old root's children must be updated so that they are
                        // no longer counted as 1 branch in the likelihood
                        // This currently updates more often than needed.
                        NodeRef oldRootNode = treeModel.getNode(cachedRoot);
                        if(oldRootNode!=node){
                            updateNodeAndChildren(oldRootNode);
                        }
                    }

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
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

        } else if (model == thorneyBranchLengthLikelihoodDelegate) {
            if (index == -1) {
                updateAllNodes();
            }  else {
            updateNode(treeModel.getNode(index));
        }
        }
        else{
            throw new RuntimeException("Unknown componentChangedEvent");
        }
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
        storedCachedRoot = cachedRoot;
        storedCachedRootChild1 = cachedRootChild1;
        storedCachedRootChild2 = cachedRootChild2;

        System.arraycopy(branchLogL, 0, storedBranchLogL, 0, branchLogL.length);
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        cachedRoot = storedCachedRoot;
        cachedRootChild1 = storedCachedRootChild1;
        cachedRootChild2 = storedCachedRootChild2;
        System.arraycopy(storedBranchLogL, 0, branchLogL, 0, branchLogL.length);
    }

    protected void acceptState() {
        cachedRoot = treeModel.getRoot().getNumber();
        cachedRootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        cachedRootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************
    public double calculateLogLikelihood(NodeRef node, int root, int rootChild1, int rootChild2) {
        int nodeIndex = node.getNumber();

        if (updateNode[nodeIndex]) {
            double logL;
            if(nodeIndex==root){
                logL=0;
            }else{
                double time = treeModel.getBranchLength(node);
                double mutations = branchLengthProvider.getBranchLength(treeModel, node);

//                    if (nodeIndex == rootChild1) {
//                        // sum the branches on both sides of the root
//                        NodeRef node2 = treeModel.getNode(rootChild2);
//                        time += treeModel.getBranchLength(node2);
//                        mutations += branchLengthProvider.getBranchLength(treeModel, node2);
//                    }
//                gamma.setScale(1.0);
//                branchLogL[i] = gamma.logPdf(x);
                logL = thorneyBranchLengthLikelihoodDelegate.getLogLikelihood(mutations, treeModel,node); //SaddlePointExpansion.logBinomialProbability((int)x, sequenceLength, expected, 1.0D - expected);

            }
            for (int i = 0; i < treeModel.getChildCount(node); i++) {
                logL += this.calculateLogLikelihood(treeModel.getChild(node, i),root,rootChild1,rootChild2);
            }
            branchLogL[nodeIndex] = logL;
            updateNode[nodeIndex] = false;
        }
        return branchLogL[nodeIndex];
    }
    private double calculateLogLikelihoodLinear() {

//        makeDirty();
//Could make this faster by only adding the changed values
        int root = treeModel.getRoot().getNumber();
        int rootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        int rootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();

        //
        updateNode[rootChild1] = updateNode[rootChild1] || updateNode[rootChild2] ;

        branchLogL[root]=0.0;
        branchLogL[rootChild2]=0.0;

        double logL = 0.0;
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if (updateNode[i] && i != root) {
                NodeRef node = treeModel.getNode(i);
                // skip the root and the second child of the root (this is added to the first child)

                double time = treeModel.getBranchLength(node);
                double mutations = branchLengthProvider.getBranchLength(treeModel, node);

//                if (i == rootChild1) {
//                    // sum the branches on both sides of the root
//                    NodeRef node2 = treeModel.getNode(rootChild2);
//                    time += treeModel.getBranchLength(node2);
//                    mutations += branchLengthProvider.getBranchLength(treeModel, node2);
//                }
////                double mean = expected * sequenceLength;
//
////                gamma.setScale(1.0);
////                branchLogL[i] = gamma.logPdf(x);

                branchLogL[i] = thorneyBranchLengthLikelihoodDelegate.getLogLikelihood(mutations,treeModel,node);
            //SaddlePointExpansion.logPoissonProbability(mean,(int)x); //SaddlePointExpansion.logBinomialProbability((int)x, sequenceLength, expected, 1.0D - expected);
            }
            updateNode[i] = false;
            logL += branchLogL[i];
        }

        return logL;
    }

    private double calculateLogLikelihood() {

        int root = treeModel.getRoot().getNumber();
        int rootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        int rootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();
        updateNode[rootChild1] = updateNode[rootChild1] || updateNode[rootChild2] ;
        assert updateNode[root];

        return calculateLogLikelihood(treeModel.getRoot(), root, rootChild1, rootChild2);

    }
    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        likelihoodKnown = false;
        updateAllNodes();
    }

    public String getReport() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    private interface LikelihoodDelegate {
        double calculateLogLikelihood(NodeRef node, int root, int rootChild1, int rootChild2);
    }

    private class additiveRateLikelihood implements LikelihoodDelegate{

        public double calculateLogLikelihood(NodeRef node, int root, int rootChild1, int rootChild2) {
            int nodeIndex = node.getNumber();

            if (updateNode[nodeIndex]) {
                double logL;
                if(nodeIndex==root){
                    logL=0;
                }else{
                    double time = treeModel.getBranchLength(node);
                    double mutations = branchLengthProvider.getBranchLength(treeModel, node);

//                    if (nodeIndex == rootChild1) {
//                        // sum the branches on both sides of the root
//                        NodeRef node2 = treeModel.getNode(rootChild2);
//                        time += treeModel.getBranchLength(node2);
//                        mutations += branchLengthProvider.getBranchLength(treeModel, node2);
//                    }
//                gamma.setScale(1.0);
//                branchLogL[i] = gamma.logPdf(x);
                        logL = thorneyBranchLengthLikelihoodDelegate.getLogLikelihood(mutations, treeModel,node); //SaddlePointExpansion.logBinomialProbability((int)x, sequenceLength, expected, 1.0D - expected);

                }
                for (int i = 0; i < treeModel.getChildCount(node); i++) {
                    logL += this.calculateLogLikelihood(treeModel.getChild(node, i),root,rootChild1,rootChild2);
                }
                branchLogL[nodeIndex] = logL;
                updateNode[nodeIndex] = false;
            }
            return branchLogL[nodeIndex];
        }
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the tree
     */
    private final TreeModel treeModel;

    private final ThorneyBranchLengthLikelihoodDelegate thorneyBranchLengthLikelihoodDelegate;
    private final BranchLengthProvider branchLengthProvider;


    //private final double[][] distanceMatrix;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown = false;
    private int cachedRoot;
    private int storedCachedRoot;
    private int cachedRootChild1;
    private int cachedRootChild2;
    private int storedCachedRootChild1;
    private int storedCachedRootChild2;
    private LikelihoodDelegate likelihoodDelegate;
    private double[] branchLogL;
    private double[] storedBranchLogL;


}

