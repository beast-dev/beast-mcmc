/*
 * EllipticalSliceOperatorParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MaskedParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.EllipticalSliceOperator;
import dr.inference.operators.MCMCOperator;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;
import org.ejml.ops.CommonOps;

/**
 */
public class EllipticalSliceOperatorParser extends AbstractXMLObjectParser {

    public static final String ELLIPTICAL_SLICE_SAMPLER = "ellipticalSliceSampler";
    public static final String SIGNAL_CONSTITUENT_PARAMETERS = "signalConstituentParameters";
    public static final String BRACKET_ANGLE = "bracketAngle";
    public static final String TRANSLATION_INVARIANT = "translationInvariant";
    public static final String ROTATION_INVARIANT = "rotationInvariant";

    public static final String DRAW_BY_ROW = "drawByRow";  // TODO What is this?

    public String getParserName() {
        return ELLIPTICAL_SLICE_SAMPLER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final Parameter variable = (Parameter) xo.getChild(Parameter.class);
        boolean drawByRowTemp = false;
        if (xo.hasAttribute(DRAW_BY_ROW)) {
            drawByRowTemp = xo.getBooleanAttribute(DRAW_BY_ROW);
        }
        final boolean drawByRow = drawByRowTemp;

        boolean signal = xo.getAttribute(SIGNAL_CONSTITUENT_PARAMETERS, true);
        if (!signal) {
            Parameter possiblyCompound = variable;
            if (variable instanceof MaskedParameter) {
                possiblyCompound = ((MaskedParameter) variable).getUnmaskedParameter();
            }
            if (!(possiblyCompound instanceof CompoundParameter)) {
                signal = true;
            }
        }

        double bracketAngle = xo.getAttribute(BRACKET_ANGLE, 0.0);

        boolean translationInvariant = xo.getAttribute(TRANSLATION_INVARIANT, false);
        boolean rotationInvariant = xo.getAttribute(ROTATION_INVARIANT, false);

        GaussianProcessRandomGenerator gaussianProcess = (GaussianProcessRandomGenerator)
                xo.getChild(GaussianProcessRandomGenerator.class);

        if (gaussianProcess == null) {

            final MultivariateDistributionLikelihood likelihood =
                    (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

            if (!(likelihood.getDistribution() instanceof GaussianProcessRandomGenerator)) {
                throw new XMLParseException("Elliptical slice sampling only works for multivariate normally distributed random variables");
            }


            if (likelihood.getDistribution() instanceof MultivariateNormalDistribution)
                gaussianProcess = (MultivariateNormalDistribution) likelihood.getDistribution();

            if (likelihood.getDistribution() instanceof MultivariateNormalDistributionModel)
                gaussianProcess = (MultivariateNormalDistributionModel) likelihood.getDistribution();

        }

        EllipticalSliceOperator operator = new EllipticalSliceOperator(variable, gaussianProcess,
                drawByRow, signal, bracketAngle,
                translationInvariant, rotationInvariant);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An elliptical slice sampler for parameters with Gaussian priors.";
    }

    public Class getReturnType() {
        return EllipticalSliceOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(SIGNAL_CONSTITUENT_PARAMETERS, true),
            AttributeRule.newDoubleRule(BRACKET_ANGLE, true),

            AttributeRule.newBooleanRule(TRANSLATION_INVARIANT, true),
            AttributeRule.newBooleanRule(ROTATION_INVARIANT, true),

            new ElementRule(Parameter.class),
            new XORRule(
                    new ElementRule(GaussianProcessRandomGenerator.class),
                    new ElementRule(MultivariateDistributionLikelihood.class)
            ),
            AttributeRule.newBooleanRule(DRAW_BY_ROW, true),
//            new ElementRule(MultivariateNormalDistribution.class),
    };
}
