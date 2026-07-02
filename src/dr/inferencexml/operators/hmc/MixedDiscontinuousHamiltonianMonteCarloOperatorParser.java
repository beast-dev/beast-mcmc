/*
 * MixedDiscontinuousHamiltonianMonteCarloOperatorParser.java
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperator;
import dr.util.Transform;
import dr.xml.*;

import static dr.util.Transform.Util.parseTransform;

/**
 * @author Filippo Monti (powered by OpenAI)
 */
public class MixedDiscontinuousHamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "mixedDiscontinuousHamiltonianMonteCarloOperator";
    public static final String STEP_SIZE = "stepSize";
    public static final String RANDOM_STEP_SIZE_FRACTION =
            DiscontinuousHamiltonianMonteCarloOperatorParser.RANDOM_STEP_SIZE_FRACTION;
    public static final String N_STEPS = "nSteps";
    public static final String MASK = "mask";
    public static final String CONTINUOUS_MASSES = "continuousMasses";
    public static final String DISCONTINUOUS_SCALES = "discontinuousScales";
    public static final String TARGET_ACCEPTANCE_PROBABILITY = "targetAcceptanceProbability";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        final double randomStepSizeFraction = xo.getAttribute(RANDOM_STEP_SIZE_FRACTION, 0.0);
        if (randomStepSizeFraction < 0.0 || randomStepSizeFraction > 1.0) {
            throw new XMLParseException("Random step size fraction must be in [0, 1]");
        }
        final int nSteps = xo.getIntegerAttribute(N_STEPS);
        final AdaptationMode adaptationMode = AdaptationMode.parseMode(xo);
        final double targetAcceptanceProbability = xo.getAttribute(TARGET_ACCEPTANCE_PROBABILITY,
                MixedDiscontinuousHamiltonianMonteCarloOperator.DEFAULT_TARGET_ACCEPTANCE_PROBABILITY);

        final GradientWrtParameterProvider gradientProvider =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        final DiscontinuousPotentialProvider discontinuousProvider =
                (DiscontinuousPotentialProvider) xo.getChild(DiscontinuousPotentialProvider.class);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        if (parameter == null) {
            parameter = discontinuousProvider.getParameter();
        }
        final Transform transform = parseTransform(xo);

        final int dimension = parameter.getDimension();
        if (gradientProvider.getDimension() != dimension || discontinuousProvider.getDimension() != dimension) {
            throw new XMLParseException("Gradient, discontinuous provider, and parameter dimensions must match");
        }
        if (transform instanceof Transform.MultivariableTransform &&
                ((Transform.MultivariableTransform) transform).getDimension() != dimension) {
            throw new XMLParseException("Transform dimension must match the parameter dimension");
        }

        if (!xo.hasChildNamed(MASK)) {
            throw new XMLParseException("A discontinuous mask is required");
        }
        final Parameter maskParameter = (Parameter) xo.getElementFirstChild(MASK);
        if (maskParameter.getDimension() != dimension) {
            throw new XMLParseException("Mask dimension must match the parameter dimension");
        }

        final boolean[] discontinuousMask = new boolean[dimension];
        for (int i = 0; i < dimension; i++) {
            discontinuousMask[i] = maskParameter.getParameterValue(i) != 0.0;
        }

        final double[] continuousMasses =
                DiscontinuousHamiltonianMonteCarloOperatorParser.parseVector(xo, CONTINUOUS_MASSES, dimension, 1.0);
        final double[] discontinuousScales =
                DiscontinuousHamiltonianMonteCarloOperatorParser.parseVector(xo, DISCONTINUOUS_SCALES, dimension, 1.0);

        final int gradientCheckCount =
                xo.getAttribute(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_COUNT, 0);
        final double gradientCheckTolerance =
                xo.getAttribute(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE, 1E-3);

        return new MixedDiscontinuousHamiltonianMonteCarloOperator(
                adaptationMode,
                gradientProvider, discontinuousProvider, discontinuousMask,
                continuousMasses, discontinuousScales, stepSize, randomStepSizeFraction, nSteps, weight,
                gradientCheckCount, gradientCheckTolerance, parameter, transform,
                targetAcceptanceProbability
        );
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(STEP_SIZE),
                AttributeRule.newDoubleRule(RANDOM_STEP_SIZE_FRACTION, true,
                        "Uniformly jitters each proposal's step size around stepSize. " +
                                "A value of 1.0 draws from [0.5 * stepSize, 1.5 * stepSize]."),
                AttributeRule.newIntegerRule(N_STEPS),
                AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
                AttributeRule.newDoubleRule(TARGET_ACCEPTANCE_PROBABILITY, true),
                AttributeRule.newIntegerRule(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_COUNT, true),
                AttributeRule.newDoubleRule(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE, true),
                new ElementRule(Parameter.class, true),
                new ElementRule(Transform.class, true),
                new ElementRule(Transform.ParsedTransform.class, true),
                new ElementRule(GradientWrtParameterProvider.class),
                new ElementRule(DiscontinuousPotentialProvider.class),
                new ElementRule(MASK, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }),
                new ElementRule(CONTINUOUS_MASSES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }, true),
                new ElementRule(DISCONTINUOUS_SCALES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns a mixed smooth/discontinuous HMC operator.";
    }

    @Override
    public Class getReturnType() {
        return MixedDiscontinuousHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
