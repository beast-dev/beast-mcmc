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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.ArbitrarySubstitutionParameterBranchModel;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchSubstitutionParameterGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.BranchParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BranchSubstitutionParameterGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "branchSubstitutionParameterGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
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
        ArbitrarySubstitutionParameterBranchModel branchModel = (ArbitrarySubstitutionParameterBranchModel) xo.getChild(BranchModel.class);

        BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();

        BranchParameter branchParameter = (BranchParameter) xo.getChild(Parameter.class);

        return new DiscreteTraitBranchSubstitutionParameterGradient(traitName, treeDataLikelihood, beagleData,
                branchParameter, useHessian);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(ArbitrarySubstitutionParameterBranchModel.class),
            new ElementRule(Parameter.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DiscreteTraitBranchSubstitutionParameterGradient.class;
    }
}
