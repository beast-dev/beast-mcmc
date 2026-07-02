/*
 * DiscontinuousHamiltonianMonteCarloOperatorParser.java
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
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.DiscontinuousHamiltonianMonteCarloOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filippo Monti (powered by OpenAI)
 */
public class DiscontinuousHamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "discontinuousHamiltonianMonteCarloOperator";
    public static final String STEP_SIZE = "stepSize";
    public static final String RANDOM_STEP_SIZE_FRACTION = "randomStepSizeFraction";
    public static final String N_STEPS = "nSteps";
    public static final String MOMENTUM_SCALES = "momentumScales";
    public static final String GRADIENT_CHECK_COUNT =
            HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_COUNT;
    public static final String GRADIENT_CHECK_TOLERANCE =
            HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        final double randomStepSizeFraction = xo.getAttribute(RANDOM_STEP_SIZE_FRACTION, 0.0);
        if (randomStepSizeFraction < 0.0 || randomStepSizeFraction > 1.0) {
            throw new XMLParseException("Random step size fraction must be in [0, 1]");
        }
        final int nSteps = xo.getIntegerAttribute(N_STEPS);

        final DiscontinuousPotentialProvider provider =
                (DiscontinuousPotentialProvider) xo.getChild(DiscontinuousPotentialProvider.class);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        if (parameter == null) {
            parameter = provider.getParameter();
        }

        if (parameter.getDimension() != provider.getDimension()) {
            throw new XMLParseException("Provider and parameter dimensions must match");
        }

        final double[] scales = parseVector(xo, MOMENTUM_SCALES, provider.getDimension(), 1.0);

        return new DiscontinuousHamiltonianMonteCarloOperator(
                provider, scales, stepSize, randomStepSizeFraction, nSteps, weight);
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
                AttributeRule.newIntegerRule(GRADIENT_CHECK_COUNT, true),
                AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
                new ElementRule(Parameter.class, true),
                new ElementRule(DiscontinuousPotentialProvider.class),
                new ElementRule(MOMENTUM_SCALES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns an all-discontinuous HMC operator with Laplace momentum.";
    }

    @Override
    public Class getReturnType() {
        return DiscontinuousHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    static double[] parseVector(XMLObject xo, String childName, int dimension, double defaultValue)
            throws XMLParseException {
        final double[] values = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            values[i] = defaultValue;
        }

        if (xo.hasChildNamed(childName)) {
            final Parameter parameter = (Parameter) xo.getElementFirstChild(childName);
            if (parameter.getDimension() != dimension) {
                throw new XMLParseException("Dimension mismatch in " + childName);
            }
            System.arraycopy(parameter.getParameterValues(), 0, values, 0, dimension);
        }

        return values;
    }
}
