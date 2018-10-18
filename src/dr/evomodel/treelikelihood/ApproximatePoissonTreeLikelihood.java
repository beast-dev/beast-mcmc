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

package dr.evomodel.treelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GammaDistribution;
import dr.xml.Reportable;

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

public class ApproximatePoissonTreeLikelihood extends AbstractModelLikelihood implements Reportable {

    public ApproximatePoissonTreeLikelihood(String name, Tree dataTree, int sequenceLength, TreeModel treeModel, BranchRateModel branchRateModel) {

        super(name);

        this.sequenceLength = sequenceLength;

        this.treeModel = treeModel;
        addModel(treeModel);

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);

        updateNode = new boolean[treeModel.getNodeCount()];
        branchLengths = new double[treeModel.getNodeCount()];

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            updateNode[i] = true;

            if (!dataTree.isRoot(dataTree.getNode(i))) {
                // Adding a small epsilon to avoid zeros.
                branchLengths[i] = dataTree.getBranchLength(dataTree.getNode(i)) * sequenceLength + 1.0E-8;
            }
        }

        distanceMatrix = new double[treeModel.getExternalNodeCount()][];
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            distanceMatrix[i] = new double[treeModel.getExternalNodeCount()];
            for (int j = i + 1; j < treeModel.getExternalNodeCount(); j++) {
                distanceMatrix[i][j] = distanceMatrix[j][i] = TreeUtils.getPathLength(dataTree, dataTree.getNode(i), dataTree.getNode(j));
            }
        }

        branchLogL = new double[treeModel.getNodeCount()];
        storedBranchLogL = new double[treeModel.getNodeCount()];

        likelihoodKnown = false;

        gamma = new GammaDistribution(1.0, 1.0);
    }

    /**
     * Set update flag for a node only
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
                    updateNodeAndChildren(((TreeChangedEvent) object).getNode());

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

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        System.arraycopy(branchLogL, 0, storedBranchLogL, 0, branchLogL.length);
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        System.arraycopy(storedBranchLogL, 0, branchLogL, 0, branchLogL.length);
    }

    protected void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    private double calculateLogLikelihood() {

//        makeDirty();

        int root = treeModel.getRoot().getNumber();
        int rootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        int rootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();

        //
        updateNode[rootChild1] = updateNode[rootChild1] || updateNode[rootChild2];

        double logL = 0.0;
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if (updateNode[i] && i != root && i != rootChild2) {
                NodeRef node = treeModel.getNode(i);
                // skip the root and the second child of the root (this is added to the first child)

                double shape = treeModel.getBranchLength(node) * branchRateModel.getBranchRate(treeModel, node);
                double x = branchLengths[i];

                if (i == rootChild1) {
                    // sum the branches on both sides of the root
                    NodeRef node2 = treeModel.getNode(rootChild2);
                    shape += treeModel.getBranchLength(node2) * branchRateModel.getBranchRate(treeModel, node2);
                    x += branchLengths[rootChild2];
                }
                gamma.setShape(shape * sequenceLength);
                gamma.setScale(1.0);
                branchLogL[i] = gamma.logPdf(x);
            }
            updateNode[i] = false;
            logL += branchLogL[i];
        }

        return logL;
    }

    private void calculateLogLikelihood(NodeRef node) {
        NodeRef c1 = treeModel.getChild(node, 0);
        if (!treeModel.isExternal(c1)) {
            calculateLogLikelihood(c1);
        }
        NodeRef c2 = treeModel.getChild(node, 1);
        if (!treeModel.isExternal(c2)) {
            calculateLogLikelihood(c2);
        }


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

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the tree
     */
    private final TreeModel treeModel;

    private final BranchRateModel branchRateModel;

    private final int sequenceLength;

    private final double[] branchLengths;
    private final double[][] distanceMatrix;

    private final GammaDistribution gamma;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown = false;

    private double[] branchLogL;
    private double[] storedBranchLogL;


}