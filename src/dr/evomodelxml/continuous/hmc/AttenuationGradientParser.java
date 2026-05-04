/*
 * AttenuationGradientParser.java
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
import dr.evomodel.treedatalikelihood.hmc.CanonicalSelectionParameterGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDiagonalAttenuationGradient;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createSelectionParameterGradient;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */


public class AttenuationGradientParser extends AbstractXMLObjectParser {
    private final static String PRECISION_GRADIENT = "attenuationGradient";
    private final static String PARAMETER = "parameter";
    private final static String IMPLEMENTATION = "implementation";
    private final static String IMPLEMENTATION_LEGACY = "legacy";
    private final static String IMPLEMENTATION_CANONICAL = "canonical";
    private final static String ATTENUATION_CORRELATION = "correlation";
    private final static String ATTENUATION_DIAGONAL = "diagonal";
    private final static String ATTENUATION_BOTH = "both";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return PRECISION_GRADIENT;
    }

    private ParameterMode parseParameterMode(XMLObject xo) throws XMLParseException {
        // Choose which parameter(s) to update.
        ParameterMode mode = ParameterMode.WRT_BOTH;
        String parameterString = xo.getAttribute(PARAMETER, ATTENUATION_BOTH).toLowerCase();
        if (parameterString.compareTo(ATTENUATION_CORRELATION) == 0) {
            mode = ParameterMode.WRT_CORRELATION;
        } else if (parameterString.compareTo(ATTENUATION_DIAGONAL) == 0) {
            mode = ParameterMode.WRT_DIAGONAL;
        }
        return mode;
    }

    enum ParameterMode {
        WRT_BOTH {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter,
                                  Parameter requestedParameter,
                                  AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter) {
                if (nativeBlockParameter != null) {
                    // The orthogonal block path is the supported tree-HMC route.
                    // Non-orthogonal block charts are still wired here, but should
                    // be treated as experimental until their live HMC gradient
                    // behavior is hardened to the same standard.
                    return createSelectionParameterGradient(
                            branchSpecificGradient,
                            treeDataLikelihood,
                            requestedParameter,
                            nativeBlockParameter);
                }
                throw new RuntimeException("Gradient wrt full attenuation not yet implemented.");
            }
        },
        WRT_CORRELATION {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter,
                                  Parameter requestedParameter,
                                  AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter) {
                throw new RuntimeException("Gradient wrt correlation of attenuation not yet implemented.");
            }
        },
        WRT_DIAGONAL {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter,
                                  Parameter requestedParameter,
                                  AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter) {
                return createDiagonalAttenuationGradient(branchSpecificGradient, treeDataLikelihood, parameter);
            }
        };

        abstract Object factory(BranchSpecificGradient branchSpecificGradient,
                                TreeDataLikelihood treeDataLikelihood,
                                MatrixParameterInterface parameter,
                                Parameter requestedParameter,
                                AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter);
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
        final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter =
                parameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter
                        ? (AbstractBlockDiagonalTwoByTwoMatrixParameter) parameter
                        : null;
        final Parameter requestedParameter =
                nativeBlockParameter != null ? nativeBlockParameter.getParameter() : null;
        final String implementation = xo.getAttribute(IMPLEMENTATION, IMPLEMENTATION_LEGACY).toLowerCase();

        if (nativeBlockParameter != null && IMPLEMENTATION_CANONICAL.equals(implementation)) {
            if (!continuousData.usesCanonicalOULikelihood()) {
                throw new XMLParseException(
                        "attenuationGradient implementation=\"canonical\" requires "
                                + "traitDataLikelihood implementation=\"canonical\".");
            }
            requestedParameter.setId(parameter.getId());
            return new CanonicalSelectionParameterGradient(
                    treeDataLikelihood,
                    continuousData,
                    requestedParameter,
                    nativeBlockParameter);
        }

        final ContinuousTraitGradientForBranch traitGradient;
        if (nativeBlockParameter != null) {
            requestedParameter.setId(parameter.getId());
            traitGradient = new ContinuousTraitGradientForBranch.SelectionParameterGradient(
                    dim, tree, continuousData, requestedParameter, nativeBlockParameter);
        } else {
            final ProcessGradientSpec gradientSpec = ProcessGradientSpec.single(
                    ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH);
            traitGradient = gradientSpec.build(dim, tree, continuousData, implementation);
        }
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient(
                        traitName,
                        treeDataLikelihood,
                        continuousData,
                        traitGradient,
                        nativeBlockParameter != null ? requestedParameter : parameter);

        ParameterMode parameterMode = parseParameterMode(xo);
        return parameterMode.factory(
                branchSpecificGradient,
                treeDataLikelihood,
                parameter,
                nativeBlockParameter != null ? requestedParameter : null,
                nativeBlockParameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class),
            new ElementRule(MatrixParameterInterface.class),
            AttributeRule.newStringRule(IMPLEMENTATION, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AbstractDiffusionGradient.ParameterDiffusionGradient.class;
    }
}
