/*
 * PurifyingGammaTreeLikelihood.java
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
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * PurifyingGammaTreeLikelihood - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @version $Id: PurifyingGammaTreeLikelihood.java,v 1.6 2006/01/10 16:48:28 rambaut Exp $
 */

@Deprecated // Switching to BEAGLE
public class PurifyingGammaTreeLikelihood extends AbstractTreeLikelihood {

    public static final String PURIFYING_GAMMA_TREE_LIKELIHOOD = "purifyingGammaTreeLikelihood";
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String HALF_LIFE = "halfLife";
    public static final String SUBSTITUTION_RATE = "substitutionRate";

    /**
     * Constructor.
     */
    public PurifyingGammaTreeLikelihood(PatternList patternList,
                                        TreeModel treeModel,
                                        GammaSiteModel siteModel,
                                        Parameter substitutionRateParameter,
                                        Parameter halfLifeParameter,
                                        boolean useAmbiguities) throws TaxonList.MissingTaxonException {

        super(PURIFYING_GAMMA_TREE_LIKELIHOOD, patternList, treeModel);

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.substitutionRateParameter = substitutionRateParameter;
            addVariable(substitutionRateParameter);

            this.halfLifeParameter = halfLifeParameter;
            addVariable(halfLifeParameter);

            this.gammaCategoryCount = siteModel.getCategoryCount();
            this.categoryCount = gammaCategoryCount;

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

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, true);

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

                if (((TreeChangedEvent) object).isNodeChanged()) {

                    updateNodeAndChildren(((TreeChangedEvent) object).getNode());

                } else {
                    updateAllNodes();

                }
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

    public void handleParameterChangedEvent(Parameter parameter, int index) {
        if (parameter == substitutionRateParameter || parameter == halfLifeParameter) {

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
     *
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        final NodeRef root = treeModel.getRoot();

        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (gammaCategoryRates == null) {
            gammaCategoryRates = new double[gammaCategoryCount];
        }

        for (int i = 0; i < gammaCategoryCount; i++) {
            gammaCategoryRates[i] = siteModel.getRateForCategory(i);
        }

        final double mu = siteModel.getMutationRateParameter().getParameterValue(0);
        final double k = substitutionRateParameter.getParameterValue(0);
        final double lambda = Math.log(2) / halfLifeParameter.getParameterValue(0);

        if (nodeTimes == null) {
            nodeTimes = new double[treeModel.getNodeCount()];
        }

        calculateNodeTimes(treeModel, root);

        traverse(treeModel, root, mu, k, lambda);

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
     *
     * @return whether the partials for this node were recalculated.
     */
    private double calculateNodeTimes(TreeModel tree, NodeRef node) {

        NodeRef parent = tree.getParent(node);

        double time0 = 0.0;

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            double t1 = calculateNodeTimes(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            double t2 = calculateNodeTimes(tree, child2);

            time0 = (t1 + t2) / 2.0;
        }

        // don't bother if you are at the root because rate at root is ignored
        if (parent == null) return 0.0;

        double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
        double time1 = time0 + branchTime;

        nodeTimes[node.getNumber()] = time0;

        return time1;
    }

    private double rateIntegral(double time, double mu, double k, double lambda) {
        return (k * time) - (((mu - k) / lambda) * (Math.exp(-lambda * time) - 1.0));
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @return whether the partials for this node were recalculated.
     */
    private boolean traverse(Tree tree, NodeRef node, double mu, double k, double lambda) {

        boolean update = false;

        final int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            double time0 = nodeTimes[node.getNumber()];
            double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
            double time1 = time0 + branchTime;

            // ***************************************************************

            for (int i = 0; i < categoryCount; i++) {

                double branchLength = rateIntegral(time1, mu, k * gammaCategoryRates[i], lambda);

                if (time0 > 0.0) {
                    branchLength -= rateIntegral(time0, mu, k * gammaCategoryRates[i], lambda);
                }

                siteModel.getTransitionProbabilities(branchLength, probabilities);
                likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            boolean update1 = traverse(tree, child1, mu, k, lambda);

            NodeRef child2 = tree.getChild(node, 1);
            boolean update2 = traverse(tree, child2, mu, k, lambda);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int childNum1 = child1.getNumber();
                int childNum2 = child2.getNumber();

                likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();

                    // moved this call to here, because non-integrating siteModels don't need to support it - AD
                    double[] proportions = siteModel.getCategoryProportions();
                    likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);

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

        public String getParserName() {
            return PURIFYING_GAMMA_TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(SUBSTITUTION_RATE);
            Parameter substitutionRateParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(HALF_LIFE);
            Parameter halfLifeParam = (Parameter) cxo.getChild(Parameter.class);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            GammaSiteModel siteModel = (GammaSiteModel) xo.getChild(GammaSiteModel.class);

            try {
                return new PurifyingGammaTreeLikelihood(patternList, treeModel, siteModel, substitutionRateParam, halfLifeParam, false);
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

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(SUBSTITUTION_RATE, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(HALF_LIFE, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(GammaSiteModel.class)
        };
    };

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the frequency model for these sites
     */
    protected FrequencyModel frequencyModel = null;

    /**
     * the site model for these sites
     */
    protected GammaSiteModel siteModel = null;

    protected Parameter substitutionRateParameter = null;
    protected Parameter halfLifeParameter = null;

    private double[] nodeTimes = null;

    /**
     * the root partial likelihoods
     */
    protected double[] branchRates = null;

    /**
     * the categories for each site
     */
    protected int[] siteCategories = null;

    /**
     * the root partial likelihoods
     */
    protected double[] rootPartials = null;

    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the pattern likelihoods
     */
    protected double[] gammaCategoryRates = null;
    /**
     * the number of rate categories
     */
    protected int gammaCategoryCount;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to store transition probabilities
     */
    protected double[] probabilities;

    /**
     * the LikelihoodCore
    */
	protected LikelihoodCore likelihoodCore;
}