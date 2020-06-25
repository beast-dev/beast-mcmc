/*
 * FullyConjugateTreeTipsPotentialDerivativeParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.branchratemodel.*;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.SumDerivative;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 */

public class BranchRateGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "branchRateGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String USE_HESSIAN = "useHessian";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        boolean useHessian = xo.getAttribute(USE_HESSIAN, false);

        final Object child = xo.getChild(TreeDataLikelihood.class);

        if (child != null) {
            return parseTreeDataLikelihood((TreeDataLikelihood) child, traitName, useHessian);
        } else {

            CompoundLikelihood compoundLikelihood = (CompoundLikelihood) xo.getChild(CompoundLikelihood.class);
            List<GradientWrtParameterProvider> providers = new ArrayList<>();

            for (Likelihood likelihood : compoundLikelihood.getLikelihoods()) {
                if (!(likelihood instanceof TreeDataLikelihood)) {
                    throw new XMLParseException("Unknown likelihood type");
                }

                GradientWrtParameterProvider provider = parseTreeDataLikelihood((TreeDataLikelihood) likelihood,
                        traitName, useHessian);

                providers.add(provider);
            }

            checkBranchRateModels(providers);

            return new SumDerivative(providers);
        }
    }

    static void checkBranchRateModels(List<GradientWrtParameterProvider> providers) throws XMLParseException {
        BranchRateModel rateModel = ((TreeDataLikelihood)providers.get(0).getLikelihood()).getBranchRateModel();
        for (GradientWrtParameterProvider provider : providers) {
            if (rateModel != ((TreeDataLikelihood)provider.getLikelihood()).getBranchRateModel()) {
                throw new XMLParseException("All TreeDataLikelihoods must use the same BranchRateModel");
            }
        }
    }

    private GradientWrtParameterProvider parseTreeDataLikelihood(TreeDataLikelihood treeDataLikelihood,
                                                                 String traitName,
                                                                 boolean useHessian) throws XMLParseException {



        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();

        if (branchRateModel instanceof DifferentiableBranchRates) {

            Parameter branchRates = ((DifferentiableBranchRates)branchRateModel).getRateParameter();

            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

            if (delegate instanceof ContinuousDataLikelihoodDelegate) {

                ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
                return new BranchRateGradient(traitName, treeDataLikelihood, continuousData, branchRates);

            } else if (delegate instanceof BeagleDataLikelihoodDelegate) {

                BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) delegate;
                return new BranchRateGradientForDiscreteTrait(traitName, treeDataLikelihood, beagleData, branchRates, useHessian);

            } else {
                throw new XMLParseException("Unknown likelihood delegate type");
            }

        } else {
            throw new XMLParseException("Only implemented for differentiable rates models");
        }
    }


    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),
            new XORRule(
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(CompoundLikelihood.class)
            ),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BranchRateGradient.class;
    }
}
