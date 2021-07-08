/*
 * BranchSpecificOptimaGradientParser.java
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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificOptimaGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Alexander Fisher
 */

public class BranchSpecificOptimaGradientParser extends AbstractXMLObjectParser {
    private final static String OPTIMA_GRADIENT = "optimaLikelihoodGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return OPTIMA_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        TreeDataLikelihood treeDataLikelihood = ((TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class));

        List<ArbitraryBranchRates> optimaBranchRates = xo.getAllChildren(ArbitraryBranchRates.class);

        CompoundParameter compoundParameter = new CompoundParameter(null);

        if (optimaBranchRates != null) {
            for (ArbitraryBranchRates brm : optimaBranchRates) {
                compoundParameter.addParameter(brm.getRateParameter());
            }
        }

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, tree, continuousData,
                        new ArrayList<>(
                                Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_BRANCH_SPECIFIC_DRIFT)
                        ));

        return new BranchSpecificOptimaGradient(traitName, treeDataLikelihood, continuousData, traitGradient, compoundParameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(ArbitraryBranchRates.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BranchSpecificOptimaGradient.class;
    }

}
