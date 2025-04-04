/*
 * TreeDataLikelihoodParser.java
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

package dr.evomodelxml.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 */
public class TreeDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BEAGLE_INSTANCE_COUNT = "beagle.instance.count";
    public static final String BEAGLE_THREAD_COUNT = "beagle.thread.count";
    public static final String THREAD_COUNT = "thread.count";

    public static final String TREE_DATA_LIKELIHOOD = "treeDataLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String INSTANCE_COUNT = "instanceCount";
    public static final String PREFER_GPU = "preferGPU";
    public static final String SCALING_SCHEME = "scalingScheme";
    public static final String DELAY_SCALING = "delayScaling";
    public static final String USE_PREORDER = "usePreOrder";
    public static final String BRANCHRATE_DERIVATIVE = "branchRateDerivative";
    public static final String BRANCHINFINITESIMAL_DERIVATIVE = "branchInfinitesimalDerivative";

    public static final String PARTITION = "partition";

    public String getParserName() {
        return TREE_DATA_LIKELIHOOD;
    }

    protected Likelihood createTreeDataLikelihood(String id,
                                                  int beagleInstanceCount,
                                                  List<PatternList> patternLists,
                                                  List<BranchModel> branchModels,
                                                  List<SiteRateModel> siteRateModels,
                                                  Tree treeModel,
                                                  BranchRateModel branchRateModel,
                                                  TipStatesModel tipStatesModel,
                                                  boolean useAmbiguities,
                                                  boolean preferGPU,
                                                  PartialsRescalingScheme scalingScheme,
                                                  boolean delayRescalingUntilUnderflow,
                                                  PreOrderSettings settings) throws XMLParseException {

        final Logger logger = Logger.getLogger("dr.evomodel");

        if (tipStatesModel != null) {
            throw new XMLParseException("Tip State Error models are not supported yet with TreeDataLikelihood");
        }

        List<Taxon> treeTaxa = treeModel.asList();
        List<Taxon> patternTaxa = patternLists.get(0).asList();

        if (!patternTaxa.containsAll(treeTaxa)) {
            throw new XMLParseException("TreeModel "+ treeModel.getId() + " contains more taxa (" + treeModel.getExternalNodeCount() + ") than the partition pattern list (" + patternTaxa.size() + ").");
        }

        if (!treeTaxa.containsAll(patternTaxa)) {
            throw new XMLParseException("TreeModel " + treeModel.getId() + " contains fewer taxa (" + treeModel.getExternalNodeCount() + ") than the partition pattern list (" + patternTaxa.size() +").");
        }

        boolean useBeagle3MultiPartition = false;

        if (patternLists.size() > 1) {
            // will currently recommend true if using GPU, CUDA or OpenCL.
            useBeagle3MultiPartition = MultiPartitionDataLikelihoodDelegate.IS_MULTI_PARTITION_RECOMMENDED();

            if (System.getProperty("USE_BEAGLE3_EXTENSIONS") != null) {
                useBeagle3MultiPartition = Boolean.parseBoolean(System.getProperty("USE_BEAGLE3_EXTENSIONS"));
            }

            if (System.getProperty("beagle.multipartition.extensions") != null &&
                    !System.getProperty("beagle.multipartition.extensions").equals("auto")) {
                useBeagle3MultiPartition = Boolean.parseBoolean(System.getProperty("beagle.multipartition.extensions"));
            }
        }

        boolean useJava = Boolean.parseBoolean(System.getProperty("java.only", "false"));
        if (useJava) {
            logger.warning("  Java-only computation is not available - ignoring this option.");
        }

        int beagleThreadCount = -1;
        if (System.getProperty(BEAGLE_THREAD_COUNT) != null) {
            // if beagle_thread_count is set then use that - this is a per-instance thread count
            beagleThreadCount = Integer.parseInt(System.getProperty(BEAGLE_THREAD_COUNT));
        }

        if (beagleThreadCount == -1) {
            // no beagle_thread_count is given so use the number of available processors
            // (actually logical threads - so 2 x number of cores when hyperthreads are used).

            beagleThreadCount = Runtime.getRuntime().availableProcessors();
            System.setProperty(BEAGLE_THREAD_COUNT, Integer.toString(beagleThreadCount));
            
            // 'threadCount' controls the top level number of Java threads holding the
            // likelihood/prior evaluations. Shouldn't be considered here - by default
            // this will use an autosizing thread pool so should probably be left alone.
            // if (System.getProperty(THREAD_COUNT) != null) {
            //  threadCount = Integer.parseInt(System.getProperty(THREAD_COUNT));
            // }
        }

        String plural = (patternLists.size() > 1 ? "s": "");
        if ( useBeagle3MultiPartition) {

            logger.info("\nCreating multi-partition tree data likelihood for " + patternLists.size() + " partition" + plural);

            if ( beagleInstanceCount > 1) {
                logger.warning("  -beagle_instances option is not compatible with BEAGLE multi-partition extensions");
            }

            try {
                DataLikelihoodDelegate dataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                        treeModel,
                        patternLists,
                        branchModels,
                        siteRateModels,
                        useAmbiguities,
                        scalingScheme,
                        delayRescalingUntilUnderflow
                );

                return new TreeDataLikelihood(
                        dataLikelihoodDelegate,
                        treeModel,
                        branchRateModel);
            } catch (DataLikelihoodDelegate.DelegateTypeException dte) {
                useBeagle3MultiPartition = false;
            }

        }

        // The multipartition data likelihood isn't available so make a set of single partition data likelihoods
        List<Likelihood> treeDataLikelihoods = new ArrayList<>();

        // Todo: allow for different number of threads per beagle instance according to pattern counts
//        if (beagleThreadCount == -1 && threadCount >= 0) {
//            System.setProperty(BEAGLE_THREAD_COUNT, Integer.toString(threadCount / patternLists.size()));
//        }

        logger.info("\nCreating tree data likelihood" + plural + " for " + patternLists.size() + " partition" + plural);

        if (beagleInstanceCount > 1) {
            logger.info("  dividing each partition between " + beagleInstanceCount + " BEAGLE instances:");
        }

        for (int i = 0; i < patternLists.size(); i++) {
            PatternList partitionPatterns = patternLists.get(i);
            // can't divide up a partition by more than the number of patterns...
            int bic = Math.min(partitionPatterns.getPatternCount(), beagleInstanceCount);

            for (int j = 0; j < bic; j++) {
                PatternList subPatterns = new Patterns(partitionPatterns, j, bic);
                DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                        treeModel,
                        subPatterns,
                        branchModels.get(i),
                        siteRateModels.get(i),
                        useAmbiguities,
                        preferGPU,
                        scalingScheme,
                        delayRescalingUntilUnderflow,
                        settings);

                TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(
                        dataLikelihoodDelegate,
                        treeModel,
                        branchRateModel);

                treeDataLikelihood.setId(id + "_" + (j + 1));

                treeDataLikelihoods.add(treeDataLikelihood);
            }
        }

        if (treeDataLikelihoods.size() == 1) {
            return treeDataLikelihoods.get(0);
        }

        return new CompoundLikelihood(treeDataLikelihoods);

    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        boolean usePreOrder = xo.getAttribute(USE_PREORDER, false);
        boolean branchRateDerivative = xo.getAttribute(BRANCHRATE_DERIVATIVE, usePreOrder);
        boolean branchInfinitesimalDerivative = xo.getAttribute(BRANCHINFINITESIMAL_DERIVATIVE, false);
        if (usePreOrder != (branchRateDerivative || branchInfinitesimalDerivative)) {
            throw new RuntimeException("Need to specify derivative types.");
        }
        PreOrderSettings settings = new PreOrderSettings(usePreOrder, branchRateDerivative, branchInfinitesimalDerivative, useAmbiguities);

        int beagleInstanceCount = xo.getAttribute(INSTANCE_COUNT, 1);
        String bic = System.getProperty(BEAGLE_INSTANCE_COUNT);
        if (bic != null && bic.length() > 0) {
            beagleInstanceCount = Math.max(1, Integer.parseInt(bic));
        }

        List<PatternList> patternLists = new ArrayList<PatternList>();
        List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();
        List<BranchModel> branchModels = new ArrayList<BranchModel>();

        boolean hasSinglePartition = false;

        PatternList patternList = (PatternList)xo.getChild(PatternList.class);
        if (patternList != null) {
            hasSinglePartition = true;
            patternLists.add(patternList);

            SiteRateModel siteRateModel = (SiteRateModel) xo.getChild(SiteRateModel.class);
            if (siteRateModel == null) {
                siteRateModel = new DiscretizedSiteRateModel("SiteRateModel");
            }
            siteRateModels.add(siteRateModel);

            FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

            BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
            if (branchModel == null) {
                SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
                if (substitutionModel == null && siteRateModel instanceof GammaSiteRateModel) {
                    // for backwards compatibility the old GammaSiteRateModelParser can provide the substitution model...
                    substitutionModel = ((GammaSiteRateModel)siteRateModel).getSubstitutionModel();
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

                SiteRateModel siteRateModel = (SiteRateModel) cxo.getChild(SiteRateModel.class);
                if (siteRateModel == null) {
                    siteRateModel = new DiscretizedSiteRateModel("SiteRateModel");
                }
                siteRateModels.add(siteRateModel);

                FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

                BranchModel branchModel = (BranchModel) cxo.getChild(BranchModel.class);
                if (branchModel == null) {
                    SubstitutionModel substitutionModel = (SubstitutionModel) cxo.getChild(SubstitutionModel.class);
                    if (substitutionModel == null && siteRateModel instanceof DiscretizedSiteRateModel) {
                        // for backwards compatibility the old GammaSiteRateModelParser can provide the substitution model...
                        substitutionModel = ((DiscretizedSiteRateModel)siteRateModel).getSubstitutionModel();
                    }
                    if (substitutionModel == null && siteRateModel instanceof GammaSiteRateModel) {
                        // for backwards compatibility the old GammaSiteRateModelParser can provide the substitution model...
                        substitutionModel = ((GammaSiteRateModel)siteRateModel).getSubstitutionModel();
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

        Tree treeModel = (Tree) xo.getChild(Tree.class);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        if (branchRateModel == null) {
            throw new XMLParseException("BranchRateModel missing from TreeDataLikelihood object '"+xo.getId());
            // branchRateModel = new DefaultBranchRateModel();
        }

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);

        final boolean preferGPU = xo.getAttribute(PREFER_GPU, false);

        PartialsRescalingScheme scalingScheme = PartialsRescalingScheme.DEFAULT;
        if (xo.hasAttribute(SCALING_SCHEME)) {
            scalingScheme = PartialsRescalingScheme.parseFromString(xo.getStringAttribute(SCALING_SCHEME));
            if (scalingScheme == null)
                throw new XMLParseException("Unknown scaling scheme '"+xo.getStringAttribute(SCALING_SCHEME)+"' in "+
                        "TreeDataLikelihood object '"+xo.getId());

        }

        final boolean delayScaling = xo.getAttribute(DELAY_SCALING, true);

        if (tipStatesModel != null) {
            throw new XMLParseException("TreeDataLikelihood is not currently compatible with TipStateModel (i.e., a sequence error model).");
        }

        return createTreeDataLikelihood(
                xo.getId(),
                beagleInstanceCount,
                patternLists,
                branchModels,
                siteRateModels,
                treeModel,
                branchRateModel,
                null,
                useAmbiguities,
                preferGPU,
                scalingScheme,
                delayScaling,
                settings);
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
            AttributeRule.newBooleanRule(PREFER_GPU, true),
            AttributeRule.newStringRule(SCALING_SCHEME,true),
            AttributeRule.newIntegerRule(INSTANCE_COUNT, true),

            // really it should be this set of elements or the PARTITION elements
            new OrRule(
                    new AndRule(new XMLSyntaxRule[]{
                            new ElementRule(PatternList.class, true),
                            new ElementRule(SubstitutionModel.class, true),
                            new ElementRule(SiteRateModel.class, true),
                            new ElementRule(FrequencyModel.class, true),
                            new ElementRule(BranchModel.class, true)}
                    ),
                    new ElementRule(PARTITION, new XMLSyntaxRule[] {
                            new ElementRule(PatternList.class),
                            new ElementRule(SubstitutionModel.class, true),
                            new ElementRule(SiteRateModel.class),
                            new ElementRule(FrequencyModel.class, true),
                            new ElementRule(BranchModel.class, true)
                    }, 1, Integer.MAX_VALUE)),

            new ElementRule(BranchRateModel.class, true),
            new ElementRule(Tree.class),
            new ElementRule(TipStatesModel.class, true)
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
