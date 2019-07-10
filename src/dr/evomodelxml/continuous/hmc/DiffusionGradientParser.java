/*
 * DiffusionGradientParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient;
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
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return DIFFUSION_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivationParametersList
                = new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>();

        CompoundParameter compoundParameter = new CompoundParameter(null);
        List<GradientWrtParameterProvider> derivativeList = new ArrayList<GradientWrtParameterProvider>();

        AbstractPrecisionGradient precisionGradient = (AbstractPrecisionGradient) xo.getChild(AbstractPrecisionGradient.class);
        if (precisionGradient != null) {
            derivationParametersList.add(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_VARIANCE);
            compoundParameter.addParameter(precisionGradient.getRawParameter());
            derivativeList.add(precisionGradient);
        }

        List<AbstractDiffusionGradient.ParameterDiffusionGradient> diffGradients = xo.getAllChildren(AbstractDiffusionGradient.ParameterDiffusionGradient.class);
        if (diffGradients != null) {
            for (AbstractDiffusionGradient.ParameterDiffusionGradient grad : diffGradients) {
                derivationParametersList.add(grad.getDerivationParameter());
                compoundParameter.addParameter(grad.getRawParameter());
                derivativeList.add(grad);
            }
        }

        CompoundGradient parametersGradients = new CompoundDerivative(derivativeList);

//        testSameModel(precisionGradient, attenuationGradient);

        TreeDataLikelihood treeDataLikelihood = ((TreeDataLikelihood) precisionGradient.getLikelihood());
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, tree, continuousData,
                        derivationParametersList);
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

