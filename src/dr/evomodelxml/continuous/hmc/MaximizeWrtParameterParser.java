/*
 * MaximumLikelihoodEstimatorParser.java
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

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class MaximizeWrtParameterParser extends AbstractXMLObjectParser {

    private static final String NAME = "maximizeWrtParameter";
    private static final String DENSITY = "densityWrtParameter";
    private static final String N_ITERATIONS = "nIterations";
    private static final String INITIAL_GUESS = "startAtCurrentState";
    private static final String PRINT_SCREEN = "printToScreen";
    private static final String INCLUDE_JACOBIAN = "includeJacobian";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider gradient =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        Parameter parameter;
        Likelihood likelihood;

        int nIterations = Math.abs(xo.getAttribute(N_ITERATIONS, 0));
        boolean initialGuess = xo.getAttribute(INITIAL_GUESS, true);
        boolean printScreen = xo.getAttribute(PRINT_SCREEN, false);
        boolean includeJacobian = xo.getAttribute(INCLUDE_JACOBIAN, false);

        if (gradient != null) {
            parameter = gradient.getParameter();
            likelihood = gradient.getLikelihood();
        } else {
            XMLObject cxo = xo.getChild(DENSITY);
            parameter = (Parameter) cxo.getChild(Parameter.class);
            likelihood = (Likelihood) cxo.getChild(Likelihood.class);
        }

        Transform transform = (Transform) xo.getChild(Transform.class);

        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(likelihood, parameter, gradient, transform,
                new MaximizerWrtParameter.Settings(nIterations, initialGuess, printScreen, includeJacobian));
        maximizer.maximize();

        return maximizer;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return MaximizerWrtParameter.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

    private static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(GradientWrtParameterProvider.class),
                    new ElementRule(DENSITY,
                            new XMLSyntaxRule[]{
                                    new ElementRule(Likelihood.class),
                                    new ElementRule(Parameter.class),
                            })
            ),
            new ElementRule(Transform.MultivariableTransformWithParameter.class, true),
            AttributeRule.newDoubleRule(N_ITERATIONS, true),
            AttributeRule.newBooleanRule(INITIAL_GUESS, true),
            AttributeRule.newBooleanRule(PRINT_SCREEN, true),
            AttributeRule.newBooleanRule(INCLUDE_JACOBIAN, true),
    };
}
