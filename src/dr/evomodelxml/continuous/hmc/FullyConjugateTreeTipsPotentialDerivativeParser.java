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

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.hmc.FullyConjugateTreeTipsPotentialDerivative;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.TreeTipGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.inferencexml.model.MaskedParameterParser;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class FullyConjugateTreeTipsPotentialDerivativeParser extends AbstractXMLObjectParser {

    private final static String FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE = "fullyConjugateTreeTipsPotentialDerivative";
    private final static String FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE2 = "traitGradientOnTree";
    public static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String MASKING = MaskedParameterParser.MASKING;

    @Override
    public String getParserName() {
        return FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE, FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        String name = xo.hasId() ? xo.getId() : FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE2;
        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
//        Object co = xo.getChild(0);

        final FullyConjugateMultivariateTraitLikelihood fcTreeLikelihood = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        Parameter mask = null;
        if (xo.hasChildNamed(MASKING)) {
            mask = (Parameter) xo.getElementFirstChild(MASKING);
        }

        if (fcTreeLikelihood != null) {

            return new FullyConjugateTreeTipsPotentialDerivative(fcTreeLikelihood, mask);

        } else if (treeDataLikelihood != null){

            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

            if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("May not provide a sequence data likelihood to compute tip trait gradient");
            }

            final ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

            return new TreeTipGradient(traitName, treeDataLikelihood, continuousData, mask);
        } else {
            throw new XMLParseException("Must provide a tree likelihood");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
                    new ElementRule(TreeDataLikelihood.class)
            ),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
            new ElementRule(Parameter.class, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FullyConjugateTreeTipsPotentialDerivative.class;
    }
}
