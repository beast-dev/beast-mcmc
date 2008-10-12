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

    public static final String TREE_LIKELIHOOD = "crazyTreeLikelihood";

    /**
     * Constructor.
     */
    public TreeLikelihood(PatternList patternList,
                          TreeModel treeModel,
                          SiteModel siteModel,
                          	                      BranchRateModel branchRateModel
    ) {

        super(TREE_LIKELIHOOD, patternList, treeModel);

        try {
            final Logger logger = Logger.getLogger("dr.evomodel");

            logger.info("Crazy new TreeLikelihood");

            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

			if (branchRateModel != null) {
				this.branchRateModel = branchRateModel;
				logger.info("Branch rate model used: " + branchRateModel.getModelName());
			} else {
				this.branchRateModel = new DefaultBranchRateModel();
			}
			addModel(this.branchRateModel);

            this.categoryCount = siteModel.getCategoryCount();


	        int[] configuration = new int[3];
	        configuration[0] = stateCount;
	        configuration[1] = patternCount;
	        configuration[2] = siteModel.getCategoryCount(); // matrixCount

            likelihoodCore = LikelihoodCoreFactory.loadLikelihoodCore(configuration, this);

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

            updateSubstitutionModel = true;
            updateSiteModel = true;

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
        hasInitialized = true;
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

		} else if (model == branchRateModel) {
			if (index == -1) {
				updateAllNodes();
			} else {
				updateNode(treeModel.getNode(index));
			}

        } else if (model == frequencyModel) {

            updateSubstitutionModel = true;
            updateAllNodes();

        } else if (model instanceof SiteModel) {

            updateSubstitutionModel = true;
            updateSiteModel = true;
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
            operations = new int[nodeCount * 3];
            dependencies = new int[nodeCount * 2];
        }

        branchUpdateCount = 0;
        operationCount = 0;

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root, null);

        if (updateSubstitutionModel) {
            likelihoodCore.updateSubstitutionModel(siteModel.getSubstitutionModel());
        }

        if (updateSiteModel) {
            likelihoodCore.updateSiteModel(siteModel);
        }

        likelihoodCore.updateMatrices(branchUpdateIndices, branchLengths, branchUpdateCount);

        likelihoodCore.updatePartials(operations, dependencies, operationCount);

        likelihoodCore.calculateLogLikelihoods(root.getNumber(), patternLogLikelihoods);

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

        updateSubstitutionModel = false;
        //********************************************************************

        return logL;
    }

    private double[] rates;

    private int[] branchUpdateIndices;
    private double[] branchLengths;
    private int branchUpdateCount;

    private int[] operations;
    private int[] dependencies;
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

			final double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            branchUpdateIndices[branchUpdateCount] = nodeNum;
            branchLengths[branchUpdateCount] = branchTime;
            branchUpdateCount++;
            
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

                int x = operationCount * 3;
                operations[x] = child1.getNumber(); // source node 1
                operations[x + 1] = child2.getNumber(); // source node 2
                operations[x + 2] = nodeNum; // destination node

                int y = operationCount * 3;
                dependencies[y] = -1; // dependent ancestor
                dependencies[y + 1] = 0; // isDependent?

                // if one of the child nodes have an update then set the dependency
                // element to this operation.
                if (op1[0] != -1) {
                    dependencies[op1[0] * 3] = operationCount;
                    dependencies[y + 1] = 1; // isDependent?
                }
                if (op2[0] != -1) {
                    dependencies[op2[0] * 3] = operationCount;
                    dependencies[y + 1] = 1; // isDependent?
                }

                if (operatorNumber != null) {
                    dependencies[y] = operationCount;
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
            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            return new TreeLikelihood(
                    patternList,
                    treeModel,
                    siteModel,
                    branchRateModel
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
                new ElementRule(BranchRateModel.class, true)
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
	 * the branch rate model
	 */
	protected final BranchRateModel branchRateModel;

    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * the LikelihoodCore
     */
    protected LikelihoodCore likelihoodCore;

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;

}