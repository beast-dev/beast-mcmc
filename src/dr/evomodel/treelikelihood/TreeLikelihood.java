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

package dr.evomodel.treelikelihood;

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
 *
 * @version $Id: TreeLikelihood.java,v 1.31 2006/08/30 16:02:42 rambaut Exp $
 */

public class TreeLikelihood extends AbstractTreeLikelihood {

    public static final String TREE_LIKELIHOOD = "treeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String USE_SCALING = "useScaling";

    /**
     * Constructor.
     */
    public TreeLikelihood(	PatternList patternList,
                              TreeModel treeModel,
                              SiteModel siteModel,
                              BranchRateModel branchRateModel,
                              boolean useAmbiguities,
                              boolean storePartials,
                              boolean useScaling)
    {

        super(TREE_LIKELIHOOD, patternList, treeModel);

        this.storePartials = storePartials;

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            integrateAcrossCategories = siteModel.integrateAcrossCategories();

            this.categoryCount = siteModel.getCategoryCount();

            if (integrateAcrossCategories)	{
                if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

	                if (NativeNucleotideLikelihoodCore.isAvailable()) {

                        Logger.getLogger("dr.evomodel").info("TreeLikelihood using native nucleotide likelihood core");
                        likelihoodCore = new NativeNucleotideLikelihoodCore();
                    } else {

                        Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java nucleotide likelihood core");
                        likelihoodCore = new NucleotideLikelihoodCore();
                    }

                } else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {
                    Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java amino acid likelihood core");
                    likelihoodCore = new AminoAcidLikelihoodCore();
                } else if (patternList.getDataType() instanceof dr.evolution.datatype.Codons) {
                    Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java codon likelihood core");
                    likelihoodCore = new CodonLikelihoodCore(patternList.getStateCount());
                    useAmbiguities = true;
                } else {
                    Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
                    likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
                }
            } else {
                Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            }
            Logger.getLogger("dr.evomodel").info( "  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            Logger.getLogger("dr.evomodel").info("  Partial likelihood scaling " + (useScaling ? "on." : "off."));

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                Logger.getLogger("dr.evomodel").info("Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultBranchRateModel();
            }
            addModel(this.branchRateModel);

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories, useScaling);

            int extNodeCount = treeModel.getExternalNodeCount();
            int intNodeCount = treeModel.getInternalNodeCount();

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                }

                if (useAmbiguities) {
                    setPartials(likelihoodCore, patternList, categoryCount, index, i);
                } else {
                    setStates(likelihoodCore, patternList, index, i);
                }
            }

            for (int i = 0; i < intNodeCount; i++) {
                likelihoodCore.createNodePartials(extNodeCount + i);
            }
        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
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

                if (((TreeModel.TreeChangedEvent)object).isNodeChanged()) {

                    updateNodeAndChildren(((TreeModel.TreeChangedEvent)object).getNode());

                } else {
                    updateAllNodes();

                }
            }

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else if (model == frequencyModel) {

           updateAllNodes();

        } else if (model instanceof SiteModel) {

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

        if (storePartials) {
            likelihoodCore.storeState();
        }
        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        if (storePartials) {
            likelihoodCore.restoreState();
        } else {
            updateAllNodes();
        }

        super.restoreState();

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (!integrateAcrossCategories) {
            if (siteCategories == null) {
                siteCategories = new int[patternCount];
            }
            for (int i = 0; i < patternCount; i++) {
                siteCategories[i] = siteModel.getCategoryOfSite(i);
            }
        }

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root);

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }
        //********************************************************************

        double logL = 0.0;

        for (int i = 0; i < patternCount; i++) {
            logL += patternLogLikelihoods[i] * patternWeights[i];
        }

        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     * @return whether the partials for this node were recalculated.
     */
    private boolean traverse(Tree tree, NodeRef node) {

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            final double branchTime = branchRate * ( tree.getNodeHeight(parent) - tree.getNodeHeight(node) );

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeNum);

            for (int i = 0; i < categoryCount; i++) {

                siteModel.getTransitionProbabilitiesForCategory(i, branchTime, probabilities);
                likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            int nodeCount = tree.getChildCount(node);
            if (nodeCount != 2)
                throw new RuntimeException("binary trees only!");

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int childNum1 = child1.getNumber();
                int childNum2 = child2.getNumber();

                likelihoodCore.setNodePartialsForUpdate(nodeNum);

                if (integrateAcrossCategories) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
                } else {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();

                    double[] partials = getRootPartials();

                    likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);
                }

                update = true;
            }
        }

        return update;

    }

	public final double[] getRootPartials() {
		if (rootPartials == null) {
		    rootPartials = new double[patternCount * stateCount];
		}

		int nodeNum = treeModel.getRoot().getNumber();
		if (integrateAcrossCategories) {

		    // moved this call to here, because non-integrating siteModels don't need to support it - AD
		    double[] proportions = siteModel.getCategoryProportions();
		    likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
		} else {
		    likelihoodCore.getPartials(nodeNum, rootPartials);
		}

		return rootPartials;
	}
	/** the root partial likelihoods (a temporary array that is used
	 * to fetch the partials - it should not be examined directly -
	 * use getRootPartials() instead).
	 */
	private double[] rootPartials = null;

    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return TREE_LIKELIHOOD; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = false;
            boolean storePartials = true;
            boolean useScaling = false;
            if (xo.hasAttribute(USE_AMBIGUITIES)) {
                useAmbiguities = xo.getBooleanAttribute(USE_AMBIGUITIES);
            }
            if (xo.hasAttribute(STORE_PARTIALS)) {
                storePartials = xo.getBooleanAttribute(STORE_PARTIALS);
            }
            if (xo.hasAttribute(USE_SCALING)) {
                useScaling = xo.getBooleanAttribute(USE_SCALING);
            }

            PatternList patternList = (PatternList)xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel)xo.getChild(SiteModel.class);

            BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);

            return new TreeLikelihood(patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials, useScaling);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of a patternlist on a tree given the site model.";
        }

        public Class getReturnType() { return Likelihood.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
                AttributeRule.newBooleanRule(STORE_PARTIALS, true),
                AttributeRule.newBooleanRule(USE_SCALING, true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true)
        };
    };

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /** the frequency model for these sites */
    protected FrequencyModel frequencyModel = null;

    /** the site model for these sites */
    protected SiteModel siteModel = null;

    /** the branch rate model  */
    protected BranchRateModel branchRateModel = null;

    private boolean storePartials = false;

    private boolean integrateAcrossCategories = false;

    /** the categories for each site */
    protected int[] siteCategories = null;


    /** the pattern likelihoods */
    protected double[] patternLogLikelihoods = null;

    /** the number of rate categories */
    protected int categoryCount;

    /** an array used to store transition probabilities */
    protected double[] probabilities;

    /** the LikelihoodCore */
    protected LikelihoodCore likelihoodCore;
}
