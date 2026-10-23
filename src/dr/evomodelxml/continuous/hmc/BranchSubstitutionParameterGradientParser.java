/*
 * BranchSubstitutionParameterGradientParser.java
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

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.BranchSpecificSubstitutionParameterBranchModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchSubstitutionParameterGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionParameterGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
import static dr.evomodel.substmodel.DifferentialMassProvider.Mode;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BranchSubstitutionParameterGradientParser extends AbstractXMLObjectParser {

    public static final String NAME = "branchSubstitutionParameterGradient";
    public static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    public static final String HOMOGENEOUS_PROCESS = "homogeneous";
    public static final String DIMENSION = "dim";
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    public static final String USE_HESSIAN = "useHessian";
    public static final String MODE = "mode";
    public static final String BACKEND = "backend";

    private static final double DEFAULT_EXACT_TOLERANCE = 1E-3;

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
        BranchSpecificSubstitutionParameterBranchModel branchModel =
                (BranchSpecificSubstitutionParameterBranchModel) xo.getChild(BranchModel.class);

        Boolean homogeneous = xo.getAttribute(HOMOGENEOUS_PROCESS, false);

        Mode mode = parseMode(xo);

        int dim = xo.getAttribute(DIMENSION, 0);
        if (homogeneous) {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            // firstOrder and affineCorrected are biased by design, so their reports only print
            final Double reportTolerance = mode == Mode.EXACT ?
                    (tolerance != null ? tolerance : DEFAULT_EXACT_TOLERANCE) : null;

            HomogeneousSubstitutionModelGradient.Backend backend = parseBackend(xo);
            getBeagleDelegate(treeDataLikelihood);
            // rate parameters keep the old path for the modes and Hessians only it provides
            final boolean rateOnAdjoint = HomogeneousSubstitutionModelGradient.isRateParameter(treeDataLikelihood, parameter)
                    && mode == Mode.EXACT && !useHessian;
            if (HomogeneousSubstitutionModelGradient.isFrequencyParameter(treeDataLikelihood, parameter)
                    || rateOnAdjoint
                    || backend == HomogeneousSubstitutionModelGradient.Backend.ADJOINT) {
                return parseHomogeneousSubstitutionModelGradient(xo, treeDataLikelihood, parameter, backend, mode,
                        useHessian, reportTolerance);
            }

            BeagleDataLikelihoodDelegate beagleData = getPreOrderBeagleDelegate(treeDataLikelihood);
            try {
                if (dim == -1) { // meaning all dimensions should be differentiated
                    return new HomogeneousSubstitutionParameterGradient(traitName, treeDataLikelihood, parameter, beagleData,
                            mode, reportTolerance);
                } else {
                    return new HomogeneousSubstitutionParameterGradient(traitName, treeDataLikelihood, parameter, beagleData,
                            dim, mode, reportTolerance);
                }
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(e.getMessage());
            }
        } else {
            BeagleDataLikelihoodDelegate beagleData = getPreOrderBeagleDelegate(treeDataLikelihood);
            DifferentiableBranchRates branchRateModel =
                    (DifferentiableBranchRates) xo.getChild(DifferentiableBranchRates.class);
            CompoundParameter branchParameter = branchModel.getBranchSpecificParameters(branchRateModel);
            return new BranchSubstitutionParameterGradient(traitName, treeDataLikelihood, beagleData,
                    branchParameter, branchRateModel, tolerance, useHessian, dim, mode);
        }
    }

    // All dimensions of one parameter, from one adjoint pass (or per dimension); frequencies add the root term
    private static Object parseHomogeneousSubstitutionModelGradient(XMLObject xo,
                                                                    TreeDataLikelihood treeDataLikelihood,
                                                                    Parameter parameter,
                                                                    HomogeneousSubstitutionModelGradient.Backend backend,
                                                                    Mode mode, boolean useHessian,
                                                                    Double reportTolerance)
            throws XMLParseException {

        // dim defaults to 0 when omitted, so route on whether it was given; dim="0" of a scalar is all of it
        if (xo.hasAttribute(DIMENSION)) {
            final int given = xo.getIntegerAttribute(DIMENSION);
            if (given != -1 && !(given == 0 && parameter.getDimension() == 1)) {
                throw new XMLParseException("The gradient wrt " + parameter.getId() + " covers all dimensions together;"
                        + " remove dim (or use dim=\"-1\"). Single frequency dimensions are not meaningful on the simplex");
            }
        }
        if (mode != Mode.EXACT) {
            throw new XMLParseException("Only mode=\"exact\" is supported for the gradient wrt " + parameter.getId());
        }
        if (useHessian) {
            throw new XMLParseException("useHessian is not supported for the gradient wrt " + parameter.getId());
        }
        if (xo.hasAttribute(TRAIT_NAME)) {
            Logger.getLogger("dr.evomodelxml").info("traitName is ignored for the gradient wrt " + parameter.getId());
        }

        try {
            return new HomogeneousSubstitutionModelGradient(treeDataLikelihood, parameter, backend, reportTolerance);
        } catch (IllegalArgumentException | UnsupportedOperationException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    private static HomogeneousSubstitutionModelGradient.Backend parseBackend(XMLObject xo) throws XMLParseException {
        try {
            return HomogeneousSubstitutionModelGradient.Backend.parse(
                    xo.getAttribute(BACKEND, HomogeneousSubstitutionModelGradient.Backend.AUTO.getLabel()));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    public static BeagleDataLikelihoodDelegate getBeagleDelegate(TreeDataLikelihood treeDataLikelihood)
            throws XMLParseException {
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new XMLParseException(NAME + " requires <treeDataLikelihood id=\"" + treeDataLikelihood.getId()
                    + "\"> to use a single BEAGLE delegate (no partition splits)");
        }
        return (BeagleDataLikelihoodDelegate) delegate;
    }

    // The per-dimension gradients write differential matrices into the branch-infinitesimal pre-order buffers
    public static BeagleDataLikelihoodDelegate getPreOrderBeagleDelegate(TreeDataLikelihood treeDataLikelihood)
            throws XMLParseException {
        BeagleDataLikelihoodDelegate beagleData = getBeagleDelegate(treeDataLikelihood);
        if (!beagleData.isUsePreOrder() || !beagleData.getPreOrderSettings().isBranchInfinitesimalDerivative()) {
            throw new XMLParseException(NAME + " requires usePreOrder=\"true\" and branchInfinitesimalDerivative=\"true\""
                    + " on <treeDataLikelihood id=\"" + treeDataLikelihood.getId() + "\">");
        }
        return beagleData;
    }

    private Mode parseMode(XMLObject xo) throws XMLParseException {
        String label = xo.getAttribute(MODE, Mode.EXACT.getLabel());
        try {
            return Mode.parse(label);
        } catch (Exception e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newBooleanRule(USE_HESSIAN, true),
            AttributeRule.newStringRule(MODE, true),
            AttributeRule.newStringRule(BACKEND, true),
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
        return GradientWrtParameterProvider.class;
    }
}
