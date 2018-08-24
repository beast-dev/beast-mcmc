/*
 * PurifyingTreeLikelihood.java
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeChangedEvent;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * PurifyingTreeLikelihood - implements a Likelihood Function for sequences on a tree.
 *
 * @version $Id: PurifyingTreeLikelihood.java,v 1.6 2006/01/10 16:48:28 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

@Deprecated // Switching to BEAGLE
public class PurifyingTreeLikelihood extends AbstractTreeLikelihood {

    public static final String PURIFYING_TREE_LIKELIHOOD = "purifyingTreeLikelihood";
    public static final String HALF_LIFE = "halfLife";
    public static final String PROPORTION = "proportion";
    public static final String AVERAGE = "average";

    /**
     * Constructor.
     */
    public PurifyingTreeLikelihood(	PatternList patternList,
                                       TreeModel treeModel,
                                       SiteModel siteModel,
                                       Parameter proportionParameter,
                                       Parameter lambdaParameter,
                                       boolean useAmbiguities,
                                       boolean useAveraging ) throws TaxonList.MissingTaxonException
    {

        super(PURIFYING_TREE_LIKELIHOOD, patternList, treeModel);

        this.useAveraging = useAveraging;

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.proportionParameter = proportionParameter;
            addVariable(proportionParameter);

            this.lambdaParameter = lambdaParameter;
            addVariable(lambdaParameter);

            integrateAcrossCategories = siteModel.integrateAcrossCategories();

            this.categoryCount = siteModel.getCategoryCount();

            if (integrateAcrossCategories)	{
                if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

                    if (NativeNucleotideLikelihoodCore.isAvailable()) {

                        System.out.println("TreeLikelihood using native nucleotide likelihood core.");
                        likelihoodCore = new NativeNucleotideLikelihoodCore();
                    } else {

                        System.out.println("TreeLikelihood using Java nucleotide likelihood core.");
                        likelihoodCore = new NucleotideLikelihoodCore();
                    }

                } else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {
                    System.out.println("TreeLikelihood using Java amino acid likelihood core.");
                    likelihoodCore = new AminoAcidLikelihoodCore();
                } else {
                    System.out.println("TreeLikelihood using Java general likelihood core.");
                    likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
                }
            } else {
                System.out.println("TreeLikelihood using Java general likelihood core.");
                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            }

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

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
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent)object).isNodeChanged()) {

                    updateNodeAndChildren(((TreeChangedEvent)object).getNode());

                } else {
                    updateAllNodes();

                }
            }

            updateRates = true;

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

    public void handleParameterChangedEvent(Parameter parameter, int index) {
        // mu has changed
        if (parameter == proportionParameter || parameter == lambdaParameter) {

            updateRates = true;
            updateAllNodes();
        }
    }

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
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        NodeRef root = treeModel.getRoot();

        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

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

        double p = proportionParameter.getParameterValue(0);
        double lambda = Math.log(2)/lambdaParameter.getParameterValue(0);

        updateRates = true;
        updateAllNodes();
        if (updateRates) {
            if (nodeTimes == null) {
                nodeTimes = new double[treeModel.getNodeCount()];
            }

            calculateNodeRates(treeModel, root, 1.0, p, lambda);
        }

        traverse(treeModel, root);

        updateRates = false;

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

        if (Double.isNaN(logL)) {
            throw new RuntimeException("Likelihood NaN");
        }

        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     * @return whether the partials for this node were recalculated.
     */
    private double calculateNodeRates(TreeModel tree, NodeRef node, double mu, double p, double lambda) {

        NodeRef parent = tree.getParent(node);

        double time0 = 0.0;

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            double t1 = calculateNodeRates(tree, child1, mu, p, lambda);

            NodeRef child2 = tree.getChild(node, 1);
            double t2 = calculateNodeRates(tree, child2, mu, p, lambda);

            if (useAveraging) {
                time0 = (t1 + t2) / 2.0;
            } else {
                // pick larger of the two
                if (t1 > t2) {
                    time0 = t1;
                } else {
                    time0 = t2;
                }
            }
        }

        // don't bother if you are at the root because rate at root is ignored
        if (parent == null) return 0;

        double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
        double time1 = time0 + branchTime;

        double branchRate = rateIntegral(time1, mu, p, lambda);

        if (time0 > 0.0) {
            branchRate -= rateIntegral(time0, mu, p, lambda);
        }

        if (branchRate != tree.getNodeRate(node)) {
            updateNode(node);
            nodeTimes[node.getNumber()] = branchRate;
        }

        return time1;
    }

    private double rateIntegral(double time, double mu, double p, double lambda) {
        return mu * ( (p * time) - (((1.0 - p) / lambda) * (Math.exp(-lambda * time) - 1.0)));
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

            // Get the average rate over this branch

            // ***************************************************************
            // Rate at nodes model
            //double branchRate = (tree.getNodeRate(node) + tree.getNodeRate(parent)) / 2;
            // ***************************************************************

            // ***************************************************************
            // Rate at branches model
            double branchTime = nodeTimes[node.getNumber()];
            // ***************************************************************


            for (int i = 0; i < categoryCount; i++) {
                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
                likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int childNum1 = child1.getNumber();
                int childNum2 = child2.getNumber();

                if (integrateAcrossCategories) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
                } else {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();

                    if (integrateAcrossCategories) {

                        // moved this call to here, because non-integrating siteModels don't need to support it - AD
                        double[] proportions = siteModel.getCategoryProportions();
                        likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
                    } else {
                        likelihoodCore.getPartials(nodeNum, rootPartials);
                    }

                    likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
                }

                update = true;
            }
        }

        return update;

    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }


    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return PURIFYING_TREE_LIKELIHOOD; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject)xo.getChild(PROPORTION);
            Parameter proportionParam = (Parameter)cxo.getChild(Parameter.class);

            cxo = (XMLObject)xo.getChild(HALF_LIFE);
            Parameter lambdaParam = (Parameter)cxo.getChild(Parameter.class);

            PatternList patternList = (PatternList)xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel)xo.getChild(SiteModel.class);

            boolean useAveraging = xo.getBooleanAttribute(AVERAGE);


            try {
                return new PurifyingTreeLikelihood(patternList, treeModel, siteModel, proportionParam, lambdaParam, false, useAveraging);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of a patternlist on a tree given the site model.";
        }

        public Class getReturnType() { return Likelihood.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(AVERAGE, false),
                new ElementRule(PROPORTION, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(HALF_LIFE, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class)
        };
    };

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /** the frequency model for these sites */
    protected FrequencyModel frequencyModel = null;

    /** the site model for these sites */
    protected SiteModel siteModel = null;

    protected Parameter proportionParameter = null;
    protected Parameter lambdaParameter = null;
    private boolean updateRates = false;

    private double[] nodeTimes = null;

    // If true then the average of the two incoming branch lengths is used for rate function in internal branches (as opposed to longest)
    private boolean useAveraging = true;

    private boolean integrateAcrossCategories = false;

    /** the root partial likelihoods */
    protected double[] branchRates = null;

    /** the categories for each site */
    protected int[] siteCategories = null;

    /** the root partial likelihoods */
    protected double[] rootPartials = null;

    /** the pattern likelihoods */
    protected double[] patternLogLikelihoods = null;

    /** the number of rate categories */
    protected int categoryCount;

    /** an array used to store transition probabilities */
    protected double[] probabilities;

    /** the LikelihoodCore */
    protected LikelihoodCore likelihoodCore;
}