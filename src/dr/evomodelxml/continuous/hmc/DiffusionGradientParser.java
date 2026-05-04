/*
 * DiffusionGradientParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodel.treedatalikelihood.hmc.DiffusionParametersGradient;
import dr.evomodel.treedatalikelihood.hmc.PrecisionGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.CompoundDerivative;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;


/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */


public class DiffusionGradientParser extends AbstractXMLObjectParser {
    private final static String DIFFUSION_GRADIENT = "diffusionGradient";
    private final static String IMPLEMENTATION = "implementation";
    private final static String IMPLEMENTATION_LEGACY = "legacy";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return DIFFUSION_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        final String implementation = xo.getAttribute(IMPLEMENTATION, IMPLEMENTATION_LEGACY).toLowerCase();

        List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivationParametersList
                = new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>();

        CompoundParameter compoundParameter = new CompoundParameter(null);
        List<GradientWrtParameterProvider> derivativeList = new ArrayList<GradientWrtParameterProvider>();

        List<AbstractDiffusionGradient> diffGradients = xo.getAllChildren(AbstractDiffusionGradient.class);
        if (diffGradients != null) {
            if (diffGradients.size() == 1) {
                return diffGradients.get(0);
            }
            boolean canUseSharedBranchSpecificPath = true;
            for (AbstractDiffusionGradient grad : diffGradients) {
                try {
                    final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter derivation =
                            grad.getDerivationParameter();
                    derivationParametersList.add(derivation);
                } catch (ClassCastException | UnsupportedOperationException ex) {
                    // Some newer gradients (e.g., native orthogonal selection-strength)
                    // are already parameter-native and do not expose continuous-process
                    // derivation metadata required by the shared branch-specific combiner.
                    canUseSharedBranchSpecificPath = false;
                    break;
                }
                compoundParameter.addParameter(grad.getRawParameter());
                derivativeList.add(grad);
            }
            if (!canUseSharedBranchSpecificPath) {
                return new CompoundDerivative(new ArrayList<GradientWrtParameterProvider>(diffGradients));
            }
        }

        CompoundGradient parametersGradients = new CompoundDerivative(derivativeList);

//        testSameModel(precisionGradient, attenuationGradient);

        TreeDataLikelihood treeDataLikelihood = ((TreeDataLikelihood) diffGradients.get(0).getLikelihood());
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
        final ProcessGradientSpec gradientSpec = ProcessGradientSpec.fromList(derivationParametersList);
        final ContinuousTraitGradientForBranch traitGradient = gradientSpec.build(
                dim, tree, continuousData, implementation);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, compoundParameter);

        return new DiffusionParametersGradient(branchSpecificGradient, parametersGradients);

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AbstractDiffusionGradient.class, 1, Integer.MAX_VALUE),
            AttributeRule.newStringRule(IMPLEMENTATION, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PrecisionGradient.class;
    }
}
