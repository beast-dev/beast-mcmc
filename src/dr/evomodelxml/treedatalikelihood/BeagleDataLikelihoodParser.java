/*
 * TreeDataLikelihoodParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class BeagleDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BEAGLE_INSTANCE_COUNT = "beagle.instance.count";

    public static final String TREE_DATA_LIKELIHOOD = "treeDataLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String INSTANCE_COUNT = "instanceCount";
    public static final String SCALING_SCHEME = "scalingScheme";
    public static final String DELAY_SCALING = "delayScaling";

    public String getParserName() {
        return TREE_DATA_LIKELIHOOD;
    }

    protected TreeDataLikelihood createTreeDataLikelihood(List<PatternList> patternLists,
                                                          TreeModel treeModel,
                                                          BranchModel branchModel,
                                                          SiteRateModel siteRateModel,
                                                          BranchRateModel branchRateModel,
                                                          TipStatesModel tipStatesModel,
                                                          boolean useAmbiguities,
                                                          PartialsRescalingScheme scalingScheme,
                                                          boolean delayRescalingUntilUnderflow,
                                                          XMLObject xo) throws XMLParseException {

        if (tipStatesModel != null) {
            throw new XMLParseException("Tip State Error models are not supported yet with TreeDataLikelihood");
        }

//        DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
//                treeModel,
//                patternLists.get(0),
//                branchModel,
//                siteRateModel,
//                useAmbiguities,
//                scalingScheme,
//                delayRescalingUntilUnderflow);

        DataLikelihoodDelegate dataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                Collections.singletonList(branchModel),
                Collections.singletonList(siteRateModel),
                useAmbiguities,
                scalingScheme,
                delayRescalingUntilUnderflow);

        return new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                branchRateModel);
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        int instanceCount = xo.getAttribute(INSTANCE_COUNT, 1);
        if (instanceCount < 1) {
            instanceCount = 1;
        }

        String ic = System.getProperty(BEAGLE_INSTANCE_COUNT);
        if (ic != null && ic.length() > 0) {
            instanceCount = Integer.parseInt(ic);
        }

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
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

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);

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

        if (instanceCount == 1 || patternList.getPatternCount() < instanceCount) {
            return createTreeDataLikelihood(
                    Collections.singletonList(patternList),
                    treeModel,
                    branchModel,
                    siteRateModel,
                    branchRateModel,
                    tipStatesModel,
                    useAmbiguities,
                    scalingScheme,
                    delayScaling,
                    xo
            );
        }

        // using multiple instances of BEAGLE...

        if (tipStatesModel != null) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with a TipStateModel (i.e., a sequence error model).");
        }

        List<PatternList> patternLists = new ArrayList<PatternList>();
        for (int i = 0; i < instanceCount; i++) {

            patternLists.add(new Patterns(patternList, i, instanceCount));
        }

        return createTreeDataLikelihood(
                    patternLists,
                    treeModel,
                    branchModel,
                    siteRateModel,
                    branchRateModel,
                    null,
                    useAmbiguities,
                    scalingScheme,
                    delayScaling,
                    xo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of data on a tree given the site model.";
    }

    public Class getReturnType() {
        return TreeDataLikelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteRateModel.class),
            new ElementRule(BranchModel.class, true),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            AttributeRule.newStringRule(SCALING_SCHEME,true),
            new ElementRule(TipStatesModel.class, true),
            new ElementRule(FrequencyModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
