/*
 * AdvancedTreeLikelihood.java
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
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * AdvancedTreeLikelihood - implements a Likelihood Function for sequences on a tree.
 * This one has some advanced models such as multiple site models for different clades.
 * This only makes sense if those clades are being constrained to remain monophyletic.
 *
 * @author Andrew Rambaut
 * @version $Id: AdvancedTreeLikelihood.java,v 1.11 2006/01/10 16:48:27 rambaut Exp $
 */

@Deprecated // Switching to BEAGLE
public class AdvancedTreeLikelihood extends AbstractTreeLikelihood {

    public static final String ADVANCED_TREE_LIKELIHOOD = "advancedTreeLikelihood";
    public static final String CLADE = "clade";
    public static final String INCLUDE_STEM = "includeStem";
    public static final String TIPS = "tips";
    public static final String DELTA = "delta";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String USE_SCALING = "useScaling";

    /**
     * Constructor.
     */
    public AdvancedTreeLikelihood(PatternList patternList,
                                  TreeModel treeModel,
                                  SiteModel siteModel,
                                  BranchRateModel branchRateModel,
                                  boolean useAmbiguities,
                                  boolean useScaling) {

        super(ADVANCED_TREE_LIKELIHOOD, patternList, treeModel);

        try {
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            if (!siteModel.integrateAcrossCategories()) {
                throw new RuntimeException("AdvancedTreeLikelihood can only use SiteModels that require integration across categories");
            }

            this.categoryCount = siteModel.getCategoryCount();

            if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

                if (NativeNucleotideLikelihoodCore.isAvailable()) {

                    Logger.getLogger("dr.evomodel").info("AdvancedTreeLikelihood using native nucleotide likelihood core.");
                    likelihoodCore = new NativeNucleotideLikelihoodCore();
                } else {

                    Logger.getLogger("dr.evomodel").info("AdvancedTreeLikelihood Java nucleotide likelihood core.");
                    likelihoodCore = new NucleotideLikelihoodCore();
                }

            } else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {
                Logger.getLogger("dr.evomodel").info("AdvancedTreeLikelihood Java amino acid likelihood core.");
                likelihoodCore = new AminoAcidLikelihoodCore();
            } else if (patternList.getDataType() instanceof dr.evolution.datatype.Codons) {
                Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
                useAmbiguities = true;
            } else {
                Logger.getLogger("dr.evomodel").info("AdvancedTreeLikelihood using Java general likelihood core");
                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            }
//            likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            Logger.getLogger("dr.evomodel").info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            Logger.getLogger("dr.evomodel").info("  Partial likelihood scaling " + (useScaling ? "on." : "off."));

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                Logger.getLogger("dr.evomodel").info("Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultBranchRateModel();
            }
            addModel(this.branchRateModel);

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

    /**
     * Add an additional siteModel for a clade in the tree.
     */
    public void addCladeSiteModel(SiteModel siteModel,
                                  TaxonList taxonList,
                                  boolean includeStem) throws TreeUtils.MissingTaxonException {
        Logger.getLogger("dr.evomodel").info("SiteModel added for clade.");
        cladeSiteModels.add(new Clade(siteModel, taxonList, includeStem));
        addModel(siteModel);
        commonAncestorsKnown = true;
    }

    /**
     * Add an additional siteModel for the tips of the tree.
     */
    public void addTipsSiteModel(SiteModel siteModel) {
        Logger.getLogger("dr.evomodel").info("SiteModel added for tips.");
        tipsSiteModel = siteModel;
        addModel(siteModel);
    }

    private void addDeltaParameter(Parameter deltaParameter, TaxonList taxa) {
        this.deltaParameter = deltaParameter;
        this.deltaTips = new HashSet<Integer>();

        if (taxa != null) {
            boolean first = true;
            StringBuffer sb = new StringBuffer("Delta parameter added for tips: {");
            for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                NodeRef node = treeModel.getExternalNode(i);
                Taxon taxon = treeModel.getNodeTaxon(node);
                if (taxa.getTaxonIndex(taxon) != -1) {
                    if (!first) {
                        sb.append(", ");
                    } else {
                        first = false;
                    }
                    sb.append(taxon.getId());
                    deltaTips.add(node.getNumber());
                }
            }
            sb.append("}");

            Logger.getLogger("dr.evomodel").info(sb.toString());
        } else {
            Logger.getLogger("dr.evomodel").info("Delta parameter added for all tips.");
        }

        addVariable(deltaParameter);
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // deltaParameter has changed...
        updateAllNodes();
        super.handleVariableChangedEvent(variable, index, type);
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

                    commonAncestorsKnown = false;

                }
            }

        } else if (model == branchRateModel) {
            updateAllNodes();

        } else if (model == frequencyModel) {

            updateAllNodes();

        } else if (model instanceof SiteModel) {

            if (model == siteModel) {

                updateAllNodes();

            } else if (model == tipsSiteModel) {

                updateAllNodes();

            } else {

                // find the siteModel in the additional siteModel list
                NodeRef node = null;
                for (int i = 0, n = cladeSiteModels.size(); i < n; i++) {
                    Clade clade = cladeSiteModels.get(i);

                    if (!commonAncestorsKnown) {
                        clade.findMRCA();
                    }

                    if (clade.getSiteModel() == model) {
                        node = treeModel.getNode(clade.getNode());
                    }
                }
                commonAncestorsKnown = true;

                updateNodeAndDescendents(node);
            }

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
     *
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

        if (!commonAncestorsKnown) {
            for (int i = 0, n = cladeSiteModels.size(); i < n; i++) {
                (cladeSiteModels.get(i)).findMRCA();
            }
            commonAncestorsKnown = true;
        }

        traverse(treeModel, root, siteModel);

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
     *
     * @return whether the partials for this node were recalculated.
     */
    private boolean traverse(Tree tree, NodeRef node, SiteModel currentSiteModel) {

        boolean update = false;

        int nodeNum = node.getNumber();

        SiteModel nextSiteModel = currentSiteModel;

        if (tipsSiteModel != null && tree.isExternal(node)) {
            currentSiteModel = tipsSiteModel;
        } else {
            for (int i = 0, n = cladeSiteModels.size(); i < n; i++) {
                Clade clade = cladeSiteModels.get(i);
                if (clade.getNode() == nodeNum) {
                    nextSiteModel = clade.getSiteModel();

                    if (clade.includeStem()) {
                        currentSiteModel = nextSiteModel;
                    }
                    break;
                }
            }
        }

        NodeRef parent = tree.getParent(node);

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeNum);

            if (tree.isExternal(node) && deltaParameter != null &&
                    (deltaTips.size() == 0 || deltaTips.contains(new Integer(node.getNumber())))) {
                branchTime += deltaParameter.getParameterValue(0);
            }

            for (int i = 0; i < categoryCount; i++) {
                double branchLength = currentSiteModel.getRateForCategory(i) * branchTime;
                currentSiteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
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
            boolean update1 = traverse(tree, child1, nextSiteModel);

            NodeRef child2 = tree.getChild(node, 1);
            boolean update2 = traverse(tree, child2, nextSiteModel);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int childNum1 = child1.getNumber();
                int childNum2 = child2.getNumber();

                likelihoodCore.setNodePartialsForUpdate(nodeNum);

                likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();
                    double[] proportions = currentSiteModel.getCategoryProportions();

                    likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
                    likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
                }

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
            return ADVANCED_TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
            boolean useScaling = xo.getAttribute(USE_SCALING, false);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            AdvancedTreeLikelihood treeLikelihood = new AdvancedTreeLikelihood(patternList, treeModel, siteModel,
                    branchRateModel, useAmbiguities,
                    useScaling);

            if (xo.hasChildNamed(TIPS)) {
                SiteModel siteModel2 = (SiteModel) xo.getElementFirstChild(TIPS);
                treeLikelihood.addTipsSiteModel(siteModel2);
            }

            XMLObject xoc = (XMLObject) xo.getChild(DELTA);
            if (xoc != null) {
                Parameter deltaParameter = (Parameter) xoc.getChild(Parameter.class);
                TaxonList taxa = (TaxonList) xoc.getChild(TaxonList.class);
                treeLikelihood.addDeltaParameter(deltaParameter, taxa);
            }

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof XMLObject) {

                    xoc = (XMLObject) xo.getChild(i);
                    if (xoc.getName().equals(CLADE)) {

                        SiteModel siteModel2 = (SiteModel) xoc.getChild(SiteModel.class);
                        TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);

                        boolean includeStem = false;

                        if (xoc.hasAttribute(INCLUDE_STEM)) {
                            includeStem = xoc.getBooleanAttribute(INCLUDE_STEM);

                            if (taxonList.getTaxonCount() == 1 && !includeStem) {
                                throw new XMLParseException("The site model is only applied to 1 taxon and therefore must include the stem branch");
                            }
                        } else if (taxonList.getTaxonCount() == 1) {
                            includeStem = true;
                        }

                        try {

                            treeLikelihood.addCladeSiteModel(siteModel2, taxonList, includeStem);

                        } catch (TreeUtils.MissingTaxonException mte) {
                            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                        }
                    }

                }
            }

            return treeLikelihood;
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

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
                // AttributeRule.newBooleanRule(STORE_PARTIALS, true),
                AttributeRule.newBooleanRule(USE_SCALING, true),
                new ElementRule(TIPS, SiteModel.class, "A siteModel that will be applied only to the tips.", 0, 1),
                new ElementRule(DELTA,
                        new XMLSyntaxRule[]{
                                new ElementRule(TaxonList.class, "A set of taxa to which to apply the delta model to", 0, 1),
                                new ElementRule(Parameter.class, "A parameter that specifies the amount of extra substitutions per site at each tip.", 0, 1),
                        }, true),
                new ElementRule(CLADE,
                        new XMLSyntaxRule[]{
                                AttributeRule.newBooleanRule(INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel."),
                                new ElementRule(TaxonList.class, "A set of taxa which defines a clade to apply a different site model to"),
                                new ElementRule(SiteModel.class, "A siteModel that will be applied only to this clade")
                        }, 0, Integer.MAX_VALUE),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true)
        };
    };

    private class Clade {
        Clade(SiteModel siteModel, TaxonList taxonList, boolean includeStem) throws TreeUtils.MissingTaxonException {
            this.siteModel = siteModel;
            this.leafSet = TreeUtils.getLeavesForTaxa(treeModel, taxonList);
            this.includeStem = includeStem;
            if (taxonList.getTaxonCount() == 1) {
                this.includeStem = true;
            }

            findMRCA();
        }

        void findMRCA() {
            node = TreeUtils.getCommonAncestorNode(treeModel, leafSet).getNumber();
        }

        int getNode() {
            return node;
        }

        boolean includeStem() {
            return includeStem;
        }

        SiteModel getSiteModel() {
            return this.siteModel;
        }

        SiteModel siteModel;
        Set<String> leafSet;
        int node;
        boolean includeStem;
    }

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
    protected SiteModel siteModel = null;

    /**
     * the branch rate model
     */
    protected BranchRateModel branchRateModel = null;

    private final boolean storePartials = false;


    /**
     * the site model for the tips
     */
    protected SiteModel tipsSiteModel = null;

    protected Parameter deltaParameter = null;
    protected Set<Integer> deltaTips = null;

    /**
     * the site models for specific clades
     */
    protected ArrayList<Clade> cladeSiteModels = new ArrayList<Clade>();

    private boolean commonAncestorsKnown = true;

    /**
     * the root partial likelihoods
     */
    protected double[] rootPartials = null;

    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

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
