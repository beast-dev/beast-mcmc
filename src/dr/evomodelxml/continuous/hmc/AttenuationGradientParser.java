/*
 * AttenuationGradientParser.java
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
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.hmc.DiagonalAttenuationGradient;
import dr.evomodel.treedatalikelihood.hmc.PrecisionGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */


public class AttenuationGradientParser extends AbstractXMLObjectParser {
    private final static String PRECISION_GRADIENT = "attenuationGradient";
    private final static String PARAMETER = "parameter";
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
                                  MatrixParameterInterface parameter) {
                throw new RuntimeException("Gradient wrt full attenuation not yet implemented.");
            }
        },
        WRT_CORRELATION {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter) {
                throw new RuntimeException("Gradient wrt correlation of attenuation not yet implemented.");
            }
        },
        WRT_DIAGONAL {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter) {
                return new DiagonalAttenuationGradient(branchSpecificGradient, treeDataLikelihood, parameter);
            }
        };

        abstract Object factory(BranchSpecificGradient branchSpecificGradient,
                                TreeDataLikelihood treeDataLikelihood,
                                MatrixParameterInterface parameter);
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
        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, tree, continuousData,
                        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, parameter);

        ParameterMode parameterMode = parseParameterMode(xo);
        return parameterMode.factory(branchSpecificGradient, treeDataLikelihood, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(WishartStatisticsWrapper.class, true),
            new ElementRule(BranchRateGradient.class, true),
            new ElementRule(Likelihood.class),
            new ElementRule(MatrixParameterInterface.class),
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
