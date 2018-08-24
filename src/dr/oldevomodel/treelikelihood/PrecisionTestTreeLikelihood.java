/*
 * PrecisionTestTreeLikelihood.java
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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.*;

import java.math.BigDecimal;
import java.util.logging.Logger;

/**
 * TreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLikelihood.java,v 1.31 2006/08/30 16:02:42 rambaut Exp $
 */

@Deprecated // Switching to BEAGLE
public class PrecisionTestTreeLikelihood extends AbstractTreeLikelihood {

    public static final String TREE_LIKELIHOOD = "precisionTestTreeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String ALLOW_MISSING_TAXA = "allowMissingTaxa";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String USE_SCALING = "useScaling";
    public static final String FORCE_JAVA_CORE = "forceJavaCore";

    /**
     * Constructor.
     */
    public PrecisionTestTreeLikelihood(PatternList patternList,
                                      TreeModel treeModel,
                                      SiteModel siteModel,
                                      BranchRateModel branchRateModel,
                                      boolean useAmbiguities,
                                      boolean allowMissingTaxa,
                                      boolean storePartials) {

        super(TREE_LIKELIHOOD, patternList, treeModel);

        this.storePartials = storePartials;

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.categoryCount = siteModel.getCategoryCount();

            final Logger logger = Logger.getLogger("dr.evomodel");
            String coreName = "Java general";
            likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            precisionLikelihoodCore = new ArbitraryPrecisionLikelihoodCore(patternList.getStateCount(), 20);

            logger.info("PrecisionTestTreeLikelihood using " + coreName + " likelihood core");

            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                logger.info("Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultBranchRateModel();
            }
            addModel(this.branchRateModel);

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, true);
            precisionLikelihoodCore.initialize(nodeCount, patternCount, categoryCount, true);

            int extNodeCount = treeModel.getExternalNodeCount();
            int intNodeCount = treeModel.getInternalNodeCount();

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    if (!allowMissingTaxa) {
                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                                ", is not found in patternList, " + patternList.getId());
                    }
                    if (useAmbiguities) {
                        setMissingPartials(likelihoodCore, i);
                        setMissingPartials(precisionLikelihoodCore, i);
                    } else {
                        setMissingStates(likelihoodCore, i);
                        setMissingStates(precisionLikelihoodCore, i);
                    }
                } else {
                    if (useAmbiguities) {
                        setPartials(likelihoodCore, patternList, categoryCount, index, i);
                        setPartials(precisionLikelihoodCore, patternList, categoryCount, index, i);
                    } else {
                        setStates(likelihoodCore, patternList, index, i);
                        setStates(precisionLikelihoodCore, patternList, index, i);
                    }
                }
            }

            for (int i = 0; i < intNodeCount; i++) {
                likelihoodCore.createNodePartials(extNodeCount + i);
                precisionLikelihoodCore.createNodePartials(extNodeCount + i);
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
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeChangedEvent) object).getNode());

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
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
            precisionLikelihoodCore.storeState();
        }
        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        if (storePartials) {
            likelihoodCore.restoreState();
            precisionLikelihoodCore.restoreState();
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
     *
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (precisionPatternLogLikelihoods == null) {
            precisionPatternLogLikelihoods = new double[patternCount];
        }

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root);

        double logL = 0.0;

        for (int i = 0; i < patternCount; i++) {
            logL += patternLogLikelihoods[i] * patternWeights[i];
        }

        double precisionLogL = 0.0;

        for (int i = 0; i < patternCount; i++) {
            precisionLogL += precisionPatternLogLikelihoods[i] * patternWeights[i];
        }

        if (Math.abs(logL - precisionLogL) > 1.0E-5) {
            System.out.println("logL = " + logL + " precision logL = " + precisionLogL);
        }

        if (logL == Double.NEGATIVE_INFINITY) {
            // We probably had an underflow... turn on scaling
            likelihoodCore.setUseScaling(true);

            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverse(treeModel, root);

            logL = 0.0;
            for (int i = 0; i < patternCount; i++) {
                logL += patternLogLikelihoods[i] * patternWeights[i];
            }

            if (Math.abs(logL - precisionLogL) > 1.0E-5) {
                System.out.println("scaled logL = " + logL + " precision logL = " + precisionLogL);
            }

        }

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }
        //********************************************************************

        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
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
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeNum);
            precisionLikelihoodCore.setNodeMatrixForUpdate(nodeNum);

            for (int i = 0; i < categoryCount; i++) {

                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
                likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
                precisionLikelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                final int childNum1 = child1.getNumber();
                final int childNum2 = child2.getNumber();

                likelihoodCore.setNodePartialsForUpdate(nodeNum);
                precisionLikelihoodCore.setNodePartialsForUpdate(nodeNum);

                likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
                precisionLikelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();

                    double[] partials = getRootPartials();
                    likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);

                    BigDecimal[] precisionPartials = getPrecisionRootPartials();
                    precisionLikelihoodCore.calculateLogLikelihoods(precisionPartials, frequencies, precisionPatternLogLikelihoods);
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
        // moved this call to here, because non-integrating siteModels don't need to support it - AD
        double[] proportions = siteModel.getCategoryProportions();
        likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);

        return rootPartials;
    }

    private double[] rootPartials = null;

    public final BigDecimal[] getPrecisionRootPartials() {
        if (precisionRootPartials == null) {
            precisionRootPartials = new BigDecimal[patternCount * stateCount];
        }

        int nodeNum = treeModel.getRoot().getNumber();
        // moved this call to here, because non-integrating siteModels don't need to support it - AD
        double[] proportions = siteModel.getCategoryProportions();
        precisionLikelihoodCore.integratePartials(nodeNum, proportions, precisionRootPartials);

        return precisionRootPartials;
    }

    private BigDecimal[] precisionRootPartials = null;

    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
            boolean allowMissingTaxa = xo.getAttribute(ALLOW_MISSING_TAXA, false);
            boolean storePartials = xo.getAttribute(STORE_PARTIALS, true);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            return new PrecisionTestTreeLikelihood(
                    patternList,
                    treeModel,
                    siteModel,
                    branchRateModel,
                    useAmbiguities, allowMissingTaxa, storePartials);
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
                AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
                AttributeRule.newBooleanRule(ALLOW_MISSING_TAXA, true),
                AttributeRule.newBooleanRule(STORE_PARTIALS, true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true),
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

    private final boolean storePartials;

    /**
     * the categories for each site
     */
    protected int[] siteCategories = null;


    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;
    /**
     * the pattern likelihoods
     */
    protected double[] precisionPatternLogLikelihoods = null;

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
    protected ArbitraryPrecisionLikelihoodCore precisionLikelihoodCore;
}