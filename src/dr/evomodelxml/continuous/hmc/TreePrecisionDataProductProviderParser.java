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

import dr.evomodel.continuous.hmc.CubicOrderTreePrecisionTraitProductProvider;
import dr.evomodel.continuous.hmc.OldLinearOrderTreePrecisionTraitProductProvider;
import dr.evomodel.continuous.hmc.LinearOrderTreePrecisionTraitProductProvider;
import dr.evomodel.continuous.hmc.TreePrecisionTraitProductProvider;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.inferencexml.model.MaskedParameterParser;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class TreePrecisionDataProductProviderParser extends AbstractXMLObjectParser {

    private static final String PRODUCT_PROVIDER = "precisionTraitProductOnTree";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String MASKING = MaskedParameterParser.MASKING;
    private static final String MODE = "mode";

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

//        Parameter mask = null;
//        if (xo.hasChildNamed(MASKING)) {
//            mask = (Parameter) xo.getElementFirstChild(MASKING);
//        }

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

        if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new XMLParseException("May not provide a sequence data likelihood to compute tip trait gradient");
        }

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

        return parseComputeMode(xo, treeDataLikelihood, continuousData, traitName);

    }

    private TreePrecisionTraitProductProvider parseComputeMode(XMLObject xo,
                                                                TreeDataLikelihood treeDataLikelihood,
                                                                ContinuousDataLikelihoodDelegate continuousData,
                                                                String traitName) throws XMLParseException {
        String mode = xo.getAttribute(MODE, "linear");
        if (mode.toLowerCase().compareTo("cubic") == 0) {
            return new CubicOrderTreePrecisionTraitProductProvider(treeDataLikelihood, continuousData);
        } else if (mode.toLowerCase().compareTo("old") == 0) {
            return new OldLinearOrderTreePrecisionTraitProductProvider(treeDataLikelihood, continuousData, traitName);
        } else {
            return new LinearOrderTreePrecisionTraitProductProvider(treeDataLikelihood, continuousData, traitName);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),
            AttributeRule.newStringRule(MODE, true),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TreePrecisionTraitProductProvider.class;
    }
}
