/*
 * MeanGradientParser.java
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
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodel.treedatalikelihood.hmc.CanonicalMeanParameterGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

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
    private final static String IMPLEMENTATION = "implementation";
    private final static String IMPLEMENTATION_LEGACY = "legacy";
    private final static String IMPLEMENTATION_CANONICAL = "canonical";
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
        } else if (parameterString.compareTo(OPT) == 0) {
            if (sharesRootMeanParameter(continuousData.getRootPrior().getMeanParameter(), parameter)) {
                mode = ParameterMode.WRT_BOTH;
            } else {
                mode = ParameterMode.WRT_DRIFT;
            }
        } else if (parameterString.compareTo(DRIFT) == 0) {
            mode = ParameterMode.WRT_DRIFT;
        } else {
            DiffusionProcessDelegate diffusionDelegate = continuousData.getDiffusionProcessDelegate();
            assert diffusionDelegate instanceof AbstractDriftDiffusionModelDelegate : "Model does not have drift.";
            assert ((AbstractDriftDiffusionModelDelegate) diffusionDelegate).isConstantDrift() : "Model does not have constant drift.";
            if (sharesRootMeanParameter(continuousData.getRootPrior().getMeanParameter(), parameter)) {
                mode = ParameterMode.WRT_BOTH;
            }
        }
        return mode;
    }

    private boolean sharesRootMeanParameter(final Parameter rootMeanParameter, final Parameter candidateParameter) {
        if (rootMeanParameter == candidateParameter) {
            return true;
        }
        if (rootMeanParameter == null || candidateParameter == null) {
            return false;
        }

        final String rootId = rootMeanParameter.getId();
        final String candidateId = candidateParameter.getId();
        if (rootId != null && rootId.equals(candidateId)) {
            return true;
        }

        if (rootMeanParameter instanceof CompoundParameter && candidateParameter instanceof CompoundParameter) {
            final CompoundParameter rootCompound = (CompoundParameter) rootMeanParameter;
            final CompoundParameter candidateCompound = (CompoundParameter) candidateParameter;
            if (rootCompound.getParameterCount() != candidateCompound.getParameterCount()) {
                return false;
            }
            for (int i = 0; i < rootCompound.getParameterCount(); ++i) {
                final Parameter rootPart = rootCompound.getParameter(i);
                final Parameter candidatePart = candidateCompound.getParameter(i);
                if (rootPart == candidatePart) {
                    continue;
                }
                final String rootPartId = rootPart.getId();
                final String candidatePartId = candidatePart.getId();
                if (rootPartId == null || !rootPartId.equals(candidatePartId)) {
                    return false;
                }
            }
            return true;
        }

        return false;
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

        final String implementation = xo.getAttribute(IMPLEMENTATION, IMPLEMENTATION_LEGACY).toLowerCase();
        final boolean useCanonicalImplementation =
                IMPLEMENTATION_CANONICAL.equals(implementation)
                        || (IMPLEMENTATION_LEGACY.equals(implementation) && continuousData.usesCanonicalOULikelihood());

        if (useCanonicalImplementation) {
            if (!continuousData.usesCanonicalOULikelihood()) {
                throw new XMLParseException(
                        "meanGradient implementation=\"canonical\" requires "
                                + "traitDataLikelihood implementation=\"canonical\".");
            }
            return new CanonicalMeanParameterGradient(
                    treeDataLikelihood,
                    continuousData,
                    parameter,
                    parameterMode != ParameterMode.WRT_ROOT,
                    parameterMode != ParameterMode.WRT_DRIFT);
        }

        final ProcessGradientSpec gradientSpec = ProcessGradientSpec.single(parameterMode.getDerivationParameter());
        final ContinuousTraitGradientForBranch traitGradient = gradientSpec.build(
                dim, tree, continuousData, implementation);
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
