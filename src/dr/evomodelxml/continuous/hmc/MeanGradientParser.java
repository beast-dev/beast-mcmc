/*
 * MeanGradientParser.java
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
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDriftGradient;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;


/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */


public class MeanGradientParser extends AbstractXMLObjectParser {
    private final static String MEAN_GRADIENT = "meanGradient";
    private final static String PARAMETER = "parameter";
    private final static String ROOT = "root";
    private final static String DRIFT = "drift";
    private final static String OPT = "opt";
    private final static String BOTH = "both";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return MEAN_GRADIENT;
    }

    private ParameterMode parseParameterMode(XMLObject xo,
                                             ContinuousDataLikelihoodDelegate continuousData,
                                             Parameter parameter) throws XMLParseException {
        // Choose which parameter(s) to update.
        ParameterMode mode = ParameterMode.WRT_DRIFT;
        String parameterString = xo.getAttribute(PARAMETER, DRIFT).toLowerCase();
        if (parameterString.compareTo(ROOT) == 0) {
            mode = ParameterMode.WRT_ROOT;
        } else if (parameterString.compareTo(BOTH) == 0) {
            mode = ParameterMode.WRT_BOTH;
        } else {
            DiffusionProcessDelegate diffusionDelegate = continuousData.getDiffusionProcessDelegate();
            assert diffusionDelegate instanceof AbstractDriftDiffusionModelDelegate : "Model does not have drift.";
            assert ((AbstractDriftDiffusionModelDelegate) diffusionDelegate).isConstantDrift() : "Model does not have constant drift.";
            if (continuousData.getRootPrior().getMeanParameter() == parameter) mode = ParameterMode.WRT_BOTH;
        }
        return mode;
    }

    enum ParameterMode {
        WRT_BOTH {
            @Override
            public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
                return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_CONSTANT_DRIFT_AND_ROOT_MEAN;
            }
        },
        WRT_DRIFT {
            @Override
            public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
                return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_CONSTANT_DRIFT;
            }
        },
        WRT_ROOT {
            @Override
            public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
                return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_ROOT_MEAN;
            }
        };

        abstract ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter();
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

        ParameterMode parameterMode = parseParameterMode(xo, continuousData, parameter);

        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, tree, continuousData,
                        new ArrayList<>(
                                Arrays.asList(parameterMode.getDerivationParameter())
                        ));
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, parameter);

        return createDriftGradient(branchSpecificGradient, treeDataLikelihood, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class),
            new ElementRule(Parameter.class),
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
