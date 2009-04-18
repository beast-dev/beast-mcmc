/*
 * BranchingLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;


/**
 * A likelihood function for branching processes. Takes a tree and a branching model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: BranchingLikelihood.java,v 1.4 2004/12/14 21:00:58 alexei Exp $
 */
public class BranchingLikelihood extends AbstractModelLikelihood {

    // PUBLIC STUFF

    public static final String BRANCHING_LIKELIHOOD = "branchingLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "branchingTree";

    public BranchingLikelihood(Tree tree, BranchingModel branchingModel) {
        this(BRANCHING_LIKELIHOOD, tree, branchingModel);
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

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
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
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BRANCHING_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) {

            XMLObject cxo = (XMLObject) xo.getChild(MODEL);
            BranchingModel branchingModel = (BranchingModel) cxo.getChild(BranchingModel.class);

            cxo = (XMLObject) xo.getChild(TREE);
            TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

            return new BranchingLikelihood(treeModel, branchingModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the demographic function.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MODEL, new XMLSyntaxRule[]{
                        new ElementRule(BranchingModel.class)
                }),
                new ElementRule(TREE, new XMLSyntaxRule[]{
                        new ElementRule(TreeModel.class)
                }),
        };
    };

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