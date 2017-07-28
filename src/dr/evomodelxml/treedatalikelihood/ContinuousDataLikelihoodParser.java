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

import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.RestrictedPartials;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.ProcessSimulationDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;
import org.w3c.dom.Attr;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class ContinuousDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CONJUGATE_ROOT_PRIOR = AbstractMultivariateTraitLikelihood.CONJUGATE_ROOT_PRIOR;
    public static final String USE_TREE_LENGTH = AbstractMultivariateTraitLikelihood.USE_TREE_LENGTH;
    public static final String SCALE_BY_TIME = AbstractMultivariateTraitLikelihood.SCALE_BY_TIME;
    public static final String RECIPROCAL_RATES = AbstractMultivariateTraitLikelihood.RECIPROCAL_RATES;
    public static final String PRIOR_SAMPLE_SIZE = AbstractMultivariateTraitLikelihood.PRIOR_SAMPLE_SIZE;

    public static final String RECONSTRUCT_TRAITS = "reconstructTraits";
    public static final String FORCE_COMPLETELY_MISSING = "forceCompletelyMissing";
    public static final String ALLOW_SINGULAR = "allowSingular";
    public static final String FORCE_FULL_PRECISION = "forceFullPrecision";

    public static final String CONTINUOUS_DATA_LIKELIHOOD = "traitDataLikelihood";

    public String getParserName() {
        return CONTINUOUS_DATA_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);


        boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);
        boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);
//        boolean reciprocalRates = xo.getAttribute(RECIPROCAL_RATES, false); // TODO Still need to add

        if (rateModel == null) {
            rateModel = new DefaultBranchRateModel();
        }

        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, scaleByTime, useTreeLength);

        final int dim = diffusionModel.getPrecisionmatrix().length;
        ConjugateRootTraitPrior rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(xo, dim);

        String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
        List<Integer> missingIndices = null;
        Parameter sampleMissingParameter = null;
        ContinuousTraitPartialsProvider dataModel = null;
        boolean useMissingIndices = true;

        if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            missingIndices = returnValue.missingIndices;
            sampleMissingParameter = returnValue.sampleMissingParameter;
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

        List<RestrictedPartials> restrictedPartialsList =
                AbstractMultivariateTraitLikelihood.parseRestrictedPartials(xo, true);

        ContinuousDataLikelihoodDelegate delegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, allowSingular);

        if (dataModel instanceof IntegratedFactorAnalysisLikelihood) {
            ((IntegratedFactorAnalysisLikelihood)dataModel).setLikelihoodDelegate(delegate);
        }

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(delegate, treeModel, rateModel);

        boolean reconstructTraits = xo.getAttribute(RECONSTRUCT_TRAITS, true);
        if (reconstructTraits) {

//            if (missingIndices != null && missingIndices.size() == 0) {
            if (!useMissingIndices) {

                ProcessSimulationDelegate simulationDelegate = new ProcessSimulationDelegate.ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate);

                TreeTraitProvider traitProvider = new ProcessSimulation(traitName,
                        treeDataLikelihood, simulationDelegate);

                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

            } else {

                ProcessSimulationDelegate simulationDelegate =
                        delegate.getPrecisionType()== PrecisionType.SCALAR ?
                                new ProcessSimulationDelegate.ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate) :
                                new ProcessSimulationDelegate.MultivariateConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                        diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate);

                TreeTraitProvider traitProvider = new ProcessSimulation(traitName,
                        treeDataLikelihood, simulationDelegate);

                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

                ProcessSimulationDelegate fullConditionalDelegate = new ProcessSimulationDelegate.TipRealizedValuesViaFullConditionalDelegate(
                        traitName, treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate);

                treeDataLikelihood.addTraits(new ProcessSimulation(("fc." + traitName), treeDataLikelihood, fullConditionalDelegate).getTreeTraits());

//                String partialTraitName = getPartiallyMissingTraitName(traitName);
//
//                ProcessSimulationDelegate parialSimulationDelegate = new ProcessSimulationDelegate.ConditionalOnPartiallyMissingTipsDelegate(partialTraitName,
//                        treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, delegate);
//
//                TreeTraitProvider partialTraitProvider = new ProcessSimulation(partialTraitName,
//                        treeDataLikelihood, parialSimulationDelegate);
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
            new ElementRule(TreeModel.class),
            new ElementRule(MultivariateDiffusionModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(CONJUGATE_ROOT_PRIOR, ConjugateRootTraitPrior.rules),
            new XORRule(
                    new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[] {
                            new ElementRule(Parameter.class),
                    }),
                    new ElementRule(ContinuousTraitPartialsProvider.class)
            ),
            new ElementRule(RestrictedPartials.class, 0, Integer.MAX_VALUE),
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
