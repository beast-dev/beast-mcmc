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
import dr.evomodel.branchratemodel.BranchRateModel;
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
 * Thorne, Kishino, & Painter (1998) Estimating the Rate of Evolution of the Rate of Molecular Evolution. Mol. Biol. Evol. 15:1647–1657.
 *
 * And more recently Didelot et al (2018) BioRxiv
 *
 * @author Andrew Rambaut
 */

public abstract class ApproximatePoissonTreeLikelihood extends AbstractModelLikelihood implements Reportable {

    public ApproximatePoissonTreeLikelihood(String name, Tree dataTree, TreeModel treeModel, BranchRateModel branchRateModel) {

        super(name);

        this.dataTree = dataTree;

        this.treeModel = treeModel;
        addModel(treeModel);

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);

        updateNode = new boolean[treeModel.getNodeCount()];
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            updateNode[i] = true;
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
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
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
        double logL = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if (updateNode[i]) {
                NodeRef node = treeModel.getNode(i);
                double shape = treeModel.getBranchLength(node) * branchRateModel.getBranchRate(treeModel, node);
                gamma.setShape(shape);
                gamma.setScale(1.0);
                branchLogL[i] = gamma.logPdf(dataTree.getBranchLength(dataTree.getNode(i)));
            }
            logL += branchLogL[i];
        }

        return logL;
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

    public boolean isLikelihoodKnown() {
        return likelihoodKnown;
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

    /**
     * the data tree
     */
    private final Tree dataTree;

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