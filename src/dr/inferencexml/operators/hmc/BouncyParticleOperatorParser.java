/*
 * HamiltonianMonteCarloOperatorParser.java
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

package dr.inferencexml.operators.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.NewBouncyParticleOperator;
import dr.inference.operators.hmc.NoUTurnOperator;
import dr.inferencexml.model.MaskedParameterParser;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperatorParser extends AbstractXMLObjectParser {
    public final static String BPO_OPERATOR = "bouncyParticleOperator";

    public final static String DRAW_VARIANCE = "drawVariance";
    public static final String MODE = "mode";

    @Override
    public String getParserName() {
        return BPO_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        System.err.println("HERE?");
        System.exit(-1);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double drawVariance = xo.getDoubleAttribute(DRAW_VARIANCE);

        CoercionMode coercionMode = CoercionMode.parseMode(xo);

        GradientWrtParameterProvider derivative =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        ContinuousDataLikelihoodDelegate likelihoodDelegate =
                (ContinuousDataLikelihoodDelegate) xo.getChild(ContinuousDataLikelihoodDelegate.class);

        String traitName = (String) xo.getChild(String.class);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        Transform transform = (Transform.MultivariableTransformWithParameter)
                xo.getChild(Transform.MultivariableTransformWithParameter.class);

        if (derivative.getDimension() != parameter.getDimension()) {
            throw new XMLParseException("Gradient (" + derivative.getDimension() +
                    ") must be the same dimensions as the parameter (" + parameter.getDimension() + ")");
        }


        return new NewBouncyParticleOperator(coercionMode, weight, treeDataLikelihood, likelihoodDelegate, traitName,
                    parameter, drawVariance);

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(DRAW_VARIANCE),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newStringRule(MODE, true),
            new ElementRule(Parameter.class),
            new ElementRule(Transform.MultivariableTransformWithParameter.class, true),
            new ElementRule(GradientWrtParameterProvider.class),
    };

    @Override
    public String getParserDescription() {
        return "Returns a Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return HamiltonianMonteCarloOperator.class;
    }
}