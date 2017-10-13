/*
 * ContinuousDataLikelihoodParser.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.continuous.*;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.preorder.ConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.MultivariateConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.TipRealizedValuesViaFullConditionalDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class ContinuousDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CONJUGATE_ROOT_PRIOR = AbstractMultivariateTraitLikelihood.CONJUGATE_ROOT_PRIOR;
    private static final String USE_TREE_LENGTH = AbstractMultivariateTraitLikelihood.USE_TREE_LENGTH;
    private static final String SCALE_BY_TIME = AbstractMultivariateTraitLikelihood.SCALE_BY_TIME;
    private static final String RECIPROCAL_RATES = AbstractMultivariateTraitLikelihood.RECIPROCAL_RATES;
    private static final String DRIFT_MODELS = AbstractMultivariateTraitLikelihood.DRIFT_MODELS;

    private static final String RECONSTRUCT_TRAITS = "reconstructTraits";
    private static final String FORCE_COMPLETELY_MISSING = "forceCompletelyMissing";
    private static final String ALLOW_SINGULAR = "allowSingular";
    private static final String FORCE_FULL_PRECISION = "forceFullPrecision";
    private static final String FORCE_DRIFT = "forceDrift";

    private static final String CONTINUOUS_DATA_LIKELIHOOD = "traitDataLikelihood";

    public String getParserName() {
        return CONTINUOUS_DATA_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree treeModel = (Tree) xo.getChild(Tree.class);
        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        
        boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);
        boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);
        boolean reciprocalRates = xo.getAttribute(RECIPROCAL_RATES, false);

        if (reciprocalRates) {
            throw new XMLParseException("Reciprocal rates are not yet implemented.");
        }

        if (rateModel == null) {
            rateModel = new DefaultBranchRateModel();
        }

        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, scaleByTime, useTreeLength);

        final int dim = diffusionModel.getPrecisionmatrix().length;
        ConjugateRootTraitPrior rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(xo, dim);

        String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
        List<Integer> missingIndices;
//        Parameter sampleMissingParameter = null;
        ContinuousTraitPartialsProvider dataModel;
        boolean useMissingIndices = true;

        if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            missingIndices = returnValue.missingIndices;
//            sampleMissingParameter = returnValue.sampleMissingParameter;
            traitName = returnValue.traitName;
            useMissingIndices = returnValue.useMissingIndices;

            PrecisionType precisionType = PrecisionType.SCALAR;

            if (xo.getAttribute(FORCE_FULL_PRECISION, false) ||
                    (useMissingIndices && !xo.getAttribute(FORCE_COMPLETELY_MISSING, false))) {
                precisionType = PrecisionType.FULL;
            }

//            System.err.println("Using precisionType == " + precisionType + " for data model.");

            dataModel = new ContinuousTraitDataModel(traitName,
                    traitParameter,
                    missingIndices, useMissingIndices,
                    dim, precisionType);
        } else {  // Has ContinuousTraitPartialsProvider
            dataModel = (ContinuousTraitPartialsProvider) xo.getChild(ContinuousTraitPartialsProvider.class);
        }

        final boolean allowSingular;
        if (dataModel instanceof IntegratedFactorAnalysisLikelihood) {
            if (xo.hasAttribute(ALLOW_SINGULAR)) {
                allowSingular = xo.getAttribute(ALLOW_SINGULAR, false);
            } else {
                allowSingular = true;
            }
        } else {
            allowSingular = xo.getAttribute(ALLOW_SINGULAR, false);
        }

        List<BranchRateModel> driftModels = AbstractMultivariateTraitLikelihood.parseDriftModels(xo, diffusionModel);

        DiffusionProcessDelegate diffusionProcessDelegate;
        if (driftModels != null || xo.getAttribute(FORCE_DRIFT, false)) {
            diffusionProcessDelegate = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
        } else {
            diffusionProcessDelegate = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        }

        ContinuousDataLikelihoodDelegate delegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, allowSingular);

        if (dataModel instanceof IntegratedFactorAnalysisLikelihood) {
            ((IntegratedFactorAnalysisLikelihood)dataModel).setLikelihoodDelegate(delegate);
        }

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(delegate, treeModel,
                rateModel);

        boolean reconstructTraits = xo.getAttribute(RECONSTRUCT_TRAITS, true);
        if (reconstructTraits) {

//            if (missingIndices != null && missingIndices.size() == 0) {
            if (!useMissingIndices) {

                ProcessSimulationDelegate simulationDelegate =
                        delegate.getPrecisionType()== PrecisionType.SCALAR ?
                                new ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, delegate) :
                                new MultivariateConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

                TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, simulationDelegate);

                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

            } else {

                ProcessSimulationDelegate simulationDelegate =
                        delegate.getPrecisionType()== PrecisionType.SCALAR ?
                                new ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, delegate) :
                                new MultivariateConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

                TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, simulationDelegate);

                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

                ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                        traitName, treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

                treeDataLikelihood.addTraits(new ProcessSimulation(treeDataLikelihood, fullConditionalDelegate).getTreeTraits());

//                String partialTraitName = getPartiallyMissingTraitName(traitName);
//
//                ProcessSimulationDelegate partialSimulationDelegate = new ProcessSimulationDelegate.ConditionalOnPartiallyMissingTipsDelegate(partialTraitName,
//                        treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate);
//
//                TreeTraitProvider partialTraitProvider = new ProcessSimulation(partialTraitName,
//                        treeDataLikelihood, partialSimulationDelegate);
//
//                treeDataLikelihood.addTraits(partialTraitProvider.getTreeTraits());
            }
        }

        return treeDataLikelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of trait data on a tree given a diffusion model.";
    }

    public Class getReturnType() {
        return TreeDataLikelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(MultivariateDiffusionModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(CONJUGATE_ROOT_PRIOR, ConjugateRootTraitPrior.rules),
            new XORRule(
                    new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[] {
                            new ElementRule(Parameter.class),
                    }),
                    new ElementRule(ContinuousTraitPartialsProvider.class)
            ),
            new ElementRule(DRIFT_MODELS, new XMLSyntaxRule[]{
                    new ElementRule(BranchRateModel.class, 1, Integer.MAX_VALUE),
            }, true),
            AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
            AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
            AttributeRule.newBooleanRule(RECIPROCAL_RATES, true),
            AttributeRule.newBooleanRule(RECONSTRUCT_TRAITS, true),
            AttributeRule.newBooleanRule(FORCE_COMPLETELY_MISSING, true),
            AttributeRule.newBooleanRule(ALLOW_SINGULAR, true),
            AttributeRule.newBooleanRule(FORCE_FULL_PRECISION, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
