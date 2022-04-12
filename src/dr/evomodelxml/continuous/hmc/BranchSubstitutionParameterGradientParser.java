/*
 * BranchSubstitutionParameterGradientParser.java
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

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.BranchSpecificSubstitutionParameterBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchSubstitutionParameterGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionParameterGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BranchSubstitutionParameterGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "branchSubstitutionParameterGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String HOMOGENEOUS_PROCESS = "homogeneous";
    private static final String DIMENSION = "dim";
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    public static final String USE_HESSIAN = "useHessian";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        boolean useHessian = xo.getAttribute(USE_HESSIAN, false);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final Double tolerance = xo.hasAttribute(GRADIENT_CHECK_TOLERANCE) ? xo.getDoubleAttribute(GRADIENT_CHECK_TOLERANCE) : null;
        BranchSpecificSubstitutionParameterBranchModel branchModel = (BranchSpecificSubstitutionParameterBranchModel) xo.getChild(BranchModel.class);

        BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();

        Boolean homogeneous = xo.getAttribute(HOMOGENEOUS_PROCESS, false);

        int dim = xo.getAttribute(DIMENSION, 0);
        if (homogeneous) {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new HomogeneousSubstitutionParameterGradient(traitName, treeDataLikelihood, parameter, beagleData, dim);
        } else {
            DifferentiableBranchRates branchRateModel = (DifferentiableBranchRates) xo.getChild(DifferentiableBranchRates.class);
            CompoundParameter branchParameter = branchModel.getBranchSpecificParameters(branchRateModel);
            return new BranchSubstitutionParameterGradient(traitName, treeDataLikelihood, beagleData,
                    branchParameter, branchRateModel, tolerance, useHessian, dim);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),
            new ElementRule(TreeDataLikelihood.class),
            new XORRule(
                    new AndRule(
                            new ElementRule(BranchSpecificSubstitutionParameterBranchModel.class),
                            new ElementRule(DifferentiableBranchRates.class)),
                    new ElementRule(Parameter.class)),
            AttributeRule.newBooleanRule(HOMOGENEOUS_PROCESS, true),
            AttributeRule.newIntegerRule(DIMENSION, true),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BranchSubstitutionParameterGradient.class;
    }
}
