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
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
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
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class TreeDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BEAGLE_INSTANCE_COUNT = "beagle.instance.count";

    public static final String TREE_DATA_LIKELIHOOD = "treeDataLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String INSTANCE_COUNT = "instanceCount";
    public static final String SCALING_SCHEME = "scalingScheme";
    public static final String DELAY_SCALING = "delayScaling";

    public static final String PARTITION = "partition";

    public String getParserName() {
        return TREE_DATA_LIKELIHOOD;
    }

    protected Likelihood createTreeDataLikelihood(List<PatternList> patternLists,
                                                  List<BranchModel> branchModels,
                                                  List<SiteRateModel> siteRateModels,
                                                  TreeModel treeModel,
                                                  BranchRateModel branchRateModel,
                                                  TipStatesModel tipStatesModel,
                                                  boolean useAmbiguities,
                                                  PartialsRescalingScheme scalingScheme,
                                                  boolean delayRescalingUntilUnderflow) throws XMLParseException {

        if (tipStatesModel != null) {
            throw new XMLParseException("Tip State Error models are not supported yet with TreeDataLikelihood");
        }

        List<Taxon> treeTaxa = treeModel.asList();
        List<Taxon> patternTaxa = patternLists.get(0).asList();

        if (!patternTaxa.containsAll(treeTaxa)) {
            throw new XMLParseException("TreeModel contains more taxa than the partition pattern list.");
        }

        if (!treeTaxa.containsAll(patternTaxa)) {
            throw new XMLParseException("TreeModel contains fewer taxa than the partition pattern list.");
        }

//        DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
//                treeModel,
//                patternLists.get(0),
//                branchModel,
//                siteRateModel,
//                useAmbiguities,
//                scalingScheme,
//                delayRescalingUntilUnderflow);

        boolean useBeagle3 = Boolean.parseBoolean(System.getProperty("USE_BEAGLE3", "true"));
        boolean useJava = Boolean.parseBoolean(System.getProperty("java.only", "false"));

        if ( useBeagle3 && MultiPartitionDataLikelihoodDelegate.IS_MULTI_PARTITION_COMPATIBLE() && !useJava) {
            DataLikelihoodDelegate dataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                    treeModel,
                    patternLists,
                    branchModels,
                    siteRateModels,
                    useAmbiguities,
                    scalingScheme,
                    delayRescalingUntilUnderflow);

            return new TreeDataLikelihood(
                    dataLikelihoodDelegate,
                    treeModel,
                    branchRateModel);
        } else {
            // The multipartition data likelihood isn't available so make a set of single partition data likelihoods
            List<Likelihood> treeDataLikelihoods = new ArrayList<Likelihood>();

            for (int i = 0; i < patternLists.size(); i++) {

                DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                        treeModel,
                        patternLists.get(i),
                        branchModels.get(i),
                        siteRateModels.get(i),
                        useAmbiguities,
                        scalingScheme,
                        delayRescalingUntilUnderflow);

                treeDataLikelihoods.add(
                        new TreeDataLikelihood(
                                dataLikelihoodDelegate,
                                treeModel,
                                branchRateModel));

            }

            if (treeDataLikelihoods.size() == 1) {
                return treeDataLikelihoods.get(0);
            }

            return new CompoundLikelihood(treeDataLikelihoods);
        }
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);

        // TreeDataLikelihood doesn't currently support Instances defined from the command line
//        int instanceCount = xo.getAttribute(INSTANCE_COUNT, 1);
//        if (instanceCount < 1) {
//            instanceCount = 1;
//        }
//
//        String ic = System.getProperty(BEAGLE_INSTANCE_COUNT);
//        if (ic != null && ic.length() > 0) {
//            instanceCount = Integer.parseInt(ic);
//        }

        List<PatternList> patternLists = new ArrayList<PatternList>();
        List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();
        List<BranchModel> branchModels = new ArrayList<BranchModel>();

        boolean hasSinglePartition = false;

        PatternList patternList = (PatternList)xo.getChild(PatternList.class);
        if (patternList != null) {
            hasSinglePartition = true;
            patternLists.add(patternList);

            GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
            siteRateModels.add(siteRateModel);

            FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

            BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
            if (branchModel == null) {
                SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
                if (substitutionModel == null) {
                    substitutionModel = siteRateModel.getSubstitutionModel();
                }
                if (substitutionModel == null) {
                    throw new XMLParseException("No substitution model available for partition in DataTreeLikelihood: "+xo.getId());
                }
                branchModel = new HomogeneousBranchModel(substitutionModel, rootFreqModel);
            }
            branchModels.add(branchModel);
        }

        int k = 0;
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChildName(i).equals(PARTITION)) {
                if (hasSinglePartition) {
                    throw new XMLParseException("Either a single set of patterns should be given or multiple 'partitions' elements within DataTreeLikelihood: "+xo.getId());
                }
                k += 1;

                XMLObject cxo = (XMLObject)xo.getChild(i);
                patternList = (PatternList) cxo.getChild(PatternList.class);
                patternLists.add(patternList);

                GammaSiteRateModel siteRateModel = (GammaSiteRateModel) cxo.getChild(GammaSiteRateModel.class);
                siteRateModels.add(siteRateModel);

                FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

                BranchModel branchModel = (BranchModel) cxo.getChild(BranchModel.class);
                if (branchModel == null) {
                    SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
                    if (substitutionModel == null) {
                        substitutionModel = siteRateModel.getSubstitutionModel();
                    }
                    if (substitutionModel == null) {
                        throw new XMLParseException("No substitution model available for partition " + k + " in DataTreeLikelihood: "+xo.getId());
                    }
                    branchModel = new HomogeneousBranchModel(substitutionModel, rootFreqModel);
                }
                branchModels.add(branchModel);

                BranchRateModel branchRateModel = (BranchRateModel) cxo.getChild(BranchRateModel.class);
                if (branchRateModel != null) {
                    throw new XMLParseException("Partitions are not currently allowed their own BranchRateModel in TreeDataLikelihood object '"+xo.getId());
                }
            }
        }

        if (patternLists.size() == 0) {
            throw new XMLParseException("Either a single set of patterns should be given or multiple 'partitions' elements within DataTreeLikelihood: "+xo.getId());
        }

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        if (branchRateModel == null) {
            throw new XMLParseException("BranchRateModel missing from TreeDataLikelihood object '"+xo.getId());
            // branchRateModel = new DefaultBranchRateModel();
        }

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);

        PartialsRescalingScheme scalingScheme = PartialsRescalingScheme.DEFAULT;
        boolean delayScaling = true;
        if (xo.hasAttribute(SCALING_SCHEME)) {
            scalingScheme = PartialsRescalingScheme.parseFromString(xo.getStringAttribute(SCALING_SCHEME));
            if (scalingScheme == null)
                throw new XMLParseException("Unknown scaling scheme '"+xo.getStringAttribute(SCALING_SCHEME)+"' in "+
                        "TreeDataLikelihood object '"+xo.getId());

        }
        if (xo.hasAttribute(DELAY_SCALING)) {
            delayScaling = xo.getBooleanAttribute(DELAY_SCALING);
        }

        if (tipStatesModel != null) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with a TipStateModel (i.e., a sequence error model).");
        }

        return createTreeDataLikelihood(
                patternLists,
                branchModels,
                siteRateModels,
                treeModel,
                branchRateModel,
                null,
                useAmbiguities,
                scalingScheme,
                delayScaling);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of data on a tree given the site model.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(SCALING_SCHEME,true),

            // really it should be this set of elements or the PARTITION elements
            new OrRule(new AndRule(new XMLSyntaxRule[]{
                    new ElementRule(PatternList.class, true),
                    new ElementRule(SiteRateModel.class, true),
                    new ElementRule(FrequencyModel.class, true),
                    new ElementRule(BranchModel.class, true)})
                    ,
                    new ElementRule(PARTITION, new XMLSyntaxRule[] {
                    new ElementRule(PatternList.class),
                    new ElementRule(SiteRateModel.class),
                    new ElementRule(FrequencyModel.class, true),
                    new ElementRule(BranchModel.class, true)
            }, 1, Integer.MAX_VALUE)),

            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TreeModel.class),
            new ElementRule(TipStatesModel.class, true)
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
