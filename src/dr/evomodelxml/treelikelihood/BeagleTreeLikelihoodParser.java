/*
 * BeagleTreeLikelihoodParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class BeagleTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BEAGLE_THREAD_COUNT = "beagle.thread.count";
    public static final String THREAD_COUNT = "thread.count";

    public static final String TREE_LIKELIHOOD = "treeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String SCALING_SCHEME = "scalingScheme";
    public static final String DELAY_SCALING = "delayScaling";
    public static final String PARTIALS_RESTRICTION = "partialsRestriction";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(PatternList patternList, MutableTreeModel treeModel,
                                                        BranchModel branchModel,
                                                        GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        TipStatesModel tipStatesModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        boolean delayRescalingUntilUnderflow,
                                                        Map<Set<String>, Parameter> partialsRestrictions,
                                                        XMLObject xo) throws XMLParseException {
        return new BeagleTreeLikelihood(
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                delayRescalingUntilUnderflow,
                partialsRestrictions
        );
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
        if (branchModel == null) {
            SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
            if (substitutionModel == null) {
                substitutionModel = siteRateModel.getSubstitutionModel();
            }
            if (substitutionModel == null) {
                throw new XMLParseException("No substitution model available for TreeLikelihood: "+xo.getId());
            }
            branchModel = new HomogeneousBranchModel(substitutionModel, rootFreqModel);
        }

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        if (branchModel instanceof EpochBranchModel && rootFreqModel != null) {
            EpochBranchModel epochBranchModel = (EpochBranchModel) branchModel;
            epochBranchModel.setRootFrequencyModel(rootFreqModel);
        }

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);
//        if (xo.getChild(TipStatesModel.class) != null) {
//            throw new XMLParseException("Sequence Error Models are not supported under BEAGLE yet. Please use Native BEAST Likelihood.");
//        }

        PartialsRescalingScheme scalingScheme = PartialsRescalingScheme.DEFAULT;
        boolean delayScaling = true;
        if (xo.hasAttribute(SCALING_SCHEME)) {
            scalingScheme = PartialsRescalingScheme.parseFromString(xo.getStringAttribute(SCALING_SCHEME));
            if (scalingScheme == null)
                throw new XMLParseException("Unknown scaling scheme '"+xo.getStringAttribute(SCALING_SCHEME)+"' in "+
                        "OldBeagleTreeLikelihood object '"+xo.getId());

        }
        if (xo.hasAttribute(DELAY_SCALING)) {
            delayScaling = xo.getBooleanAttribute(DELAY_SCALING);
        }

        Map<Set<String>, Parameter> partialsRestrictions = null;

        if (xo.hasChildNamed(PARTIALS_RESTRICTION)) {
            XMLObject cxo = xo.getChild(PARTIALS_RESTRICTION);
            TaxonList taxonList = (TaxonList) cxo.getChild(TaxonList.class);
//            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);
            try {
                TreeUtils.getLeavesForTaxa(treeModel, taxonList);
            } catch (TreeUtils.MissingTaxonException e) {
                throw new XMLParseException("Unable to parse taxon list: " + e.getMessage());
            }
            throw new XMLParseException("Restricting internal nodes is not yet implemented.  Contact Marc");

        }

        int beagleThreadCount = -1;
        if (System.getProperty(BEAGLE_THREAD_COUNT) != null) {
            beagleThreadCount = Integer.parseInt(System.getProperty(BEAGLE_THREAD_COUNT));
        }

        if (beagleThreadCount == -1) {
            // no beagle_thread_count is given so use the number of available processors
            // (actually logical threads - so 2 x number of cores when hyperthreads are used).

            beagleThreadCount = Runtime.getRuntime().availableProcessors();
        }

        System.setProperty(BEAGLE_THREAD_COUNT, Integer.toString(beagleThreadCount));

        return createTreeLikelihood(
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                delayScaling,
                partialsRestrictions,
                xo
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

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchModel.class, true),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            AttributeRule.newStringRule(SCALING_SCHEME,true),
            new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[] {
                    new ElementRule(TaxonList.class),
                    new ElementRule(Parameter.class),
            }, true),
            new ElementRule(TipStatesModel.class, true),
            new ElementRule(FrequencyModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
