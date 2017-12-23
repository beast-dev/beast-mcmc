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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.NoUTurnOperator;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser {

    private final static String HMC_OPERATOR = "hamiltonianMonteCarloOperator";
    private final static String N_STEPS = "nSteps";
    private final static String STEP_SIZE = "stepSize";
    private final static String DRAW_VARIANCE = "drawVariance";
    private final static String MODE = "mode";
    private final static String NUTS = "nuts";
    private final static String VANILLA = "vanilla";

    @Override
    public String getParserName() {
        return HMC_OPERATOR;
    }

    private int parseRunMode(XMLObject xo) throws XMLParseException {
        int mode = 0;
        if (xo.getAttribute(MODE, VANILLA).toLowerCase().compareTo(NUTS) == 0) {
            mode = 1;
        }
        return mode;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int nSteps = xo.getIntegerAttribute(N_STEPS);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double drawVariance = xo.getDoubleAttribute(DRAW_VARIANCE);
        int runMode = parseRunMode(xo);

        CoercionMode coercionMode = CoercionMode.parseMode(xo);

        GradientWrtParameterProvider derivative =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        Transform transform = (Transform.MultivariableTransformWithParameter)
                xo.getChild(Transform.MultivariableTransformWithParameter.class);

        if (derivative.getDimension() != parameter.getDimension()) {
            throw new XMLParseException("Gradient (" + derivative.getDimension() +
                    ") must be the same dimensions as the parameter (" + parameter.getDimension() + ")");
        }


        if (runMode == 0) {
            return new HamiltonianMonteCarloOperator(coercionMode, weight, derivative, parameter, transform,
                    stepSize, nSteps, drawVariance);
        } else {
            return new NoUTurnOperator(coercionMode, weight, derivative, parameter,transform,
                    stepSize, nSteps, drawVariance);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS),
            AttributeRule.newDoubleRule(STEP_SIZE),
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
