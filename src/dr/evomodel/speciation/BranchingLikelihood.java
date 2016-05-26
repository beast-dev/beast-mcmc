/*
 * BranchingLikelihood.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.speciation.BranchingLikelihoodParser;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;


/**
 * A likelihood function for branching processes. Takes a tree and a branching model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: BranchingLikelihood.java,v 1.4 2004/12/14 21:00:58 alexei Exp $
 */
public class BranchingLikelihood extends AbstractModelLikelihood {

    // PUBLIC STUFF

    public BranchingLikelihood(Tree tree, BranchingModel branchingModel) {
        this(BranchingLikelihoodParser.BRANCHING_LIKELIHOOD, tree, branchingModel);
    }

    public BranchingLikelihood(String name, Tree tree, BranchingModel branchingModel) {

        super(name);

        this.tree = tree;
        this.branchingModel = branchingModel;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
        if (branchingModel != null) {
            addModel(branchingModel);
        }
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************


    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        return calculateLogLikelihood();

    }

    public final void makeDirty() {
    }

    /**
     * Calculates the log likelihood of this set of tree nodes,
     * given a branching model.
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;
        for (int j = 0; j < tree.getInternalNodeCount(); j++) {
            logL += branchingModel.logNodeProbability(tree, tree.getInternalNode(j));
        }
        //System.err.println("logL=" + logL);
        return logL;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /**
     * The speciation model.
     */
    BranchingModel branchingModel = null;

    /**
     * The tree.
     */
    Tree tree = null;
}