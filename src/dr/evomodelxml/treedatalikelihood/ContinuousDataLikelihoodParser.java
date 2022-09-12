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
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.*;
import dr.evomodelxml.continuous.ContinuousTraitDataModelParser;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;

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
    private static final String OPTIMAL_TRAITS = AbstractMultivariateTraitLikelihood.OPTIMAL_TRAITS;

    private static final String RECONSTRUCT_TRAITS = "reconstructTraits";
    private static final String ALLOW_SINGULAR = "allowSingular";
    private static final String FORCE_DRIFT = "forceDrift";
    private static final String FORCE_OU = "forceOU";

    private static final String STRENGTH_OF_SELECTION_MATRIX = "strengthOfSelectionMatrix";


    public static final String CONTINUOUS_DATA_LIKELIHOOD = "traitDataLikelihood";

    public static final String FACTOR_NAME = "factors";

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

        final String traitName;

        ContinuousTraitPartialsProvider dataModel;

        if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
            dataModel = ContinuousTraitDataModelParser.parseContinuousTraitDataModel(xo);
        } else {  // Has ContinuousTraitPartialsProvider
            dataModel = (ContinuousTraitPartialsProvider) xo.getChild(ContinuousTraitPartialsProvider.class);
        }

        traitName = xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME);
        dataModel.setTipTraitName(getTipTraitName(traitName)); // TODO: not an ideal solution as the trait name could be set differently later

        ConjugateRootTraitPrior rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(xo, dataModel.getTraitDimension());

        final boolean allowSingular;

        if (xo.hasAttribute(ALLOW_SINGULAR)) {
            //TODO: check compatibility (there are cases where allowSingular=false guarantees the incorrect likelihood)
            allowSingular = xo.getBooleanAttribute(ALLOW_SINGULAR);
        } else {
            allowSingular = dataModel.getDefaultAllowSingular();
        }

        dataModel.addTreeAndRateModel(treeModel, rateTransformation);

        List<BranchRateModel> driftModels = AbstractMultivariateTraitLikelihood.parseDriftModels(xo, diffusionModel);
        List<BranchRateModel> optimalTraitsModels = AbstractMultivariateTraitLikelihood.parseOptimalValuesModels(xo, diffusionModel);

        MultivariateElasticModel elasticModel = null;
        if (xo.hasChildNamed(STRENGTH_OF_SELECTION_MATRIX)) {
            XMLObject cxo = xo.getChild(STRENGTH_OF_SELECTION_MATRIX);
            MatrixParameterInterface strengthOfSelectionMatrixParam;
            strengthOfSelectionMatrixParam = (MatrixParameterInterface) cxo.getChild(MatrixParameterInterface.class);
            if (strengthOfSelectionMatrixParam != null) {
                elasticModel = new MultivariateElasticModel(strengthOfSelectionMatrixParam);
            }
        }

        DiffusionProcessDelegate diffusionProcessDelegate;
        boolean integratedProcess = dataModel instanceof IntegratedProcessTraitDataModel; //TODO: can add to interface if that would be better
        if ((optimalTraitsModels != null && elasticModel != null) || xo.getAttribute(FORCE_OU, false)) {
            if (!integratedProcess) {
                diffusionProcessDelegate = new OUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, elasticModel);
            } else {
                diffusionProcessDelegate = new IntegratedOUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, elasticModel);
            }
        } else {
            if (driftModels != null || xo.getAttribute(FORCE_DRIFT, false)) {
                diffusionProcessDelegate = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
            } else {
                diffusionProcessDelegate = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
            }
        }
        ContinuousDataLikelihoodDelegate delegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, allowSingular);

        if (dataModel instanceof IntegratedFactorAnalysisLikelihood) {
            ((IntegratedFactorAnalysisLikelihood) dataModel).setLikelihoodDelegate(delegate);
        }

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(delegate, treeModel,
                rateModel);

        boolean reconstructTraits = xo.getAttribute(RECONSTRUCT_TRAITS, true);
        if (reconstructTraits) {

//            if (missingIndices != null && missingIndices.size() == 0) {

            ProcessSimulationDelegate simulationDelegate =
                    delegate.getPrecisionType() == PrecisionType.SCALAR ?
                            new ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                    diffusionModel, dataModel, rootPrior, rateTransformation, delegate) :
                            new MultivariateConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                    diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, simulationDelegate);

            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

            if (dataModel.usesMissingIndices()) {

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

            //TODO: remove below (should let ConditionalTraitSimulationHelper figure everything out)
            int[] partitionDimensions = dataModel.getPartitionDimensions();
            if (partitionDimensions.length > 1) {
                PartitionedTreeTraitProvider partitionedProvider =
                        new PartitionedTreeTraitProvider(treeDataLikelihood.getTreeTraits(), partitionDimensions);
                treeDataLikelihood.addTraits(partitionedProvider.getTreeTraits());
            }

        }

        delegate.setExtensionHelper();

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
            new XORRule(new XMLSyntaxRule[]{
                    new AndRule(ContinuousTraitDataModelParser.rules),
                    new ElementRule(ContinuousTraitPartialsProvider.class)
            }),
            new ElementRule(DRIFT_MODELS, new XMLSyntaxRule[]{
                    new ElementRule(BranchRateModel.class, 1, Integer.MAX_VALUE),
            }, true),
            new ElementRule(OPTIMAL_TRAITS, new XMLSyntaxRule[]{
                    new ElementRule(BranchRateModel.class, 1, Integer.MAX_VALUE),
            }, true),
            new ElementRule(STRENGTH_OF_SELECTION_MATRIX,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}, true),
            AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
            AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
            AttributeRule.newBooleanRule(RECIPROCAL_RATES, true),
            AttributeRule.newBooleanRule(RECONSTRUCT_TRAITS, true),
            AttributeRule.newBooleanRule(ALLOW_SINGULAR, true),
            AttributeRule.newBooleanRule(FORCE_DRIFT, true),
            AttributeRule.newBooleanRule(FORCE_OU, true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            TreeTraitParserUtilities.jitterRules(true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
