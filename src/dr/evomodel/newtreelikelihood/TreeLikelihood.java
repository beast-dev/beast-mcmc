/*
 * TreeLikelihood.java
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

package dr.evomodel.newtreelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * TreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLikelihood.java,v 1.31 2006/08/30 16:02:42 rambaut Exp $
 */

public class TreeLikelihood extends AbstractTreeLikelihood {

    public static final String TREE_LIKELIHOOD = "treeLikelihood";

    /**
     * Constructor.
     */
    public TreeLikelihood(PatternList patternList,
                          TreeModel treeModel,
                          SiteModel siteModel) {

        super(TREE_LIKELIHOOD, patternList, treeModel);

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.categoryCount = siteModel.getCategoryCount();

            final Logger logger = Logger.getLogger("dr.evomodel");

            likelihoodCore = LikelihoodCoreFactory.getLikelihoodCore(stateCount);

            logger.info("Crazy new TreeLikelihood");

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount);

            int extNodeCount = treeModel.getExternalNodeCount();

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                } else {
                    setPartials(likelihoodCore, patternList, categoryCount, index, i);
                }
            }

            updateSubstitutionMatrix = true;

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
    }

    public final LikelihoodCore getLikelihoodCore() {
        return likelihoodCore;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

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
                    // Currently this event type is not used.
                    System.err.println("Full tree update event - these events currently aren't used\n" +
                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == frequencyModel) {

            updateSubstitutionMatrix = true;
            updateAllNodes();

        } else if (model instanceof SiteModel) {

            updateSubstitutionMatrix = true;
            updateAllNodes();

        } else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }

        super.handleModelChangedEvent(model, object, index);
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {

        likelihoodCore.storeState();
        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodCore.restoreState();

        super.restoreState();

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (branchUpdateIndices == null) {
            branchUpdateIndices = new int[nodeCount];
            branchLengths = new double[nodeCount];
        }

        if (operations == null) {
            operations = new int[nodeCount];
        }

        branchUpdateCount = 0;
        operationCount = 0;

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root, null);

        if (updateSubstitutionMatrix) {
            likelihoodCore.updateSubstitutionModel(siteModel.getSubstitutionModel());
        }

        if (rates == null) {
            rates = new double[siteModel.getCategoryCount()];
        }

        for (int i = 0; i < categoryCount; i++) {
            rates[i] = siteModel.getRateForCategory(i);
        }

        likelihoodCore.updateMatrices(branchUpdateIndices, branchLengths, branchUpdateCount, rates);

        likelihoodCore.updatePartials(operations, operationCount);

        double[] frequencies = frequencyModel.getFrequencies();
        double[] proportions = siteModel.getCategoryProportions();

        likelihoodCore.calculateLogLikelihoods(root.getNumber(), frequencies, proportions, patternLogLikelihoods);

        double logL = 0.0;
        for (int i = 0; i < patternCount; i++) {
            logL += patternLogLikelihoods[i] * patternWeights[i];
        }

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }

        updateSubstitutionMatrix = false;
        //********************************************************************

        return logL;
    }

    private double[] rates;

    private int[] branchUpdateIndices;
    private double[] branchLengths;
    private int branchUpdateCount;

    private int[] operations;
    private int operationCount;

    /**
     * Traverse the tree calculating partial likelihoods.
     */
    private boolean traverse(Tree tree, NodeRef node, int[] operatorNumber) {

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate = 1.0;

            // Get the operational time of the branch
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            branchUpdateIndices[branchUpdateCount] = nodeNum;
            branchLengths[branchUpdateCount] = branchTime;

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final int[] op1 = new int[] { -1 };
            final boolean update1 = traverse(tree, child1, op1);

            NodeRef child2 = tree.getChild(node, 1);
            final int[] op2 = new int[] { -1 };
            final boolean update2 = traverse(tree, child2, op2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int x = operationCount * 5;
                operations[x] = -1; // dependent ancestor
                operations[x + 1] = 0; // isDependent?
                operations[x + 2] = child1.getNumber(); // source node 1
                operations[x + 3] = child2.getNumber(); // source node 2
                operations[x + 4] = nodeNum; // destination node

                // if one of the child nodes have an update then set the dependency
                // element to this operation.
                if (op1[0] != -1) {
                    operations[op1[0] * 4] = operationCount;
                    operations[x + 1] = 1; // isDependent?
                }
                if (op2[0] != -1) {
                    operations[op2[0] * 4] = operationCount;
                    operations[x + 1] = 1; // isDependent?
                }

                if (operatorNumber != null) {
                    operatorNumber[0] = operationCount;
                }

                operationCount ++;

                update = true;
            }
        }

        return update;

    }


    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

            return new TreeLikelihood(
                    patternList,
                    treeModel,
                    siteModel
            );
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of a patternlist on a tree given the site model.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
        };
    };

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the frequency model for these sites
     */
    protected final FrequencyModel frequencyModel;

    /**
     * the site model for these sites
     */
    protected final SiteModel siteModel;


    /**
     * the categories for each site
     */
    protected int[] siteCategories = null;


    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[] probabilities;


    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * the LikelihoodCore
     */
    protected LikelihoodCore likelihoodCore;

    /**
     * Flag to specify that the subsitution matrix has changed
     */
    protected boolean updateSubstitutionMatrix;
    
}