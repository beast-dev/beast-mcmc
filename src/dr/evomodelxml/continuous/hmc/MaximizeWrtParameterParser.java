/*
 * MaximizeWrtParameterParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.hmc.NumericalGradient;
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
    private static final String MAX_SECONDS = "maxSeconds";
    private static final String MAX_EVALUATIONS = "maxEvaluations";
    private static final String REL_TOL = "relTol";
    private static final String STALL_WINDOW = "stallWindow";
    private static final String CHECK_GRADIENT_EVERY = "checkGradientEvery";
    private static final String EPSILON = "epsilon";
    private static final String PAST = "past";
    private static final String DELTA = "delta";

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

        if (gradient != null && !(gradient instanceof NumericalGradient)
                && TreeModel.containsSampledTipHeight(parameter)) {
            throw new XMLParseException(TreeModel.TIP_DATE_GRADIENT_NOT_IMPLEMENTED);
        }

        Transform transform = (Transform) xo.getChild(Transform.class);

        final double numericGradientStepSize = NumericGradientStepSizeProvider.parseStepSizeRatio(xo);

        double maxSeconds = xo.getAttribute(MAX_SECONDS, 0.0);
        long maxEvaluations = xo.getLongIntegerAttribute(MAX_EVALUATIONS, 0);
        double relTol = xo.getAttribute(REL_TOL, 0.0);
        int stallWindow = xo.getAttribute(STALL_WINDOW, 0);
        int checkGradientEvery = xo.getAttribute(CHECK_GRADIENT_EVERY, 0);
        if (maxSeconds < 0 || maxEvaluations < 0 || relTol < 0 || stallWindow < 0 || checkGradientEvery < 0) {
            throw new XMLParseException(MAX_SECONDS + ", " + MAX_EVALUATIONS + ", " + REL_TOL + ", " + STALL_WINDOW +
                    " and " + CHECK_GRADIENT_EVERY + " must not be negative");
        }

        Double epsilon = xo.hasAttribute(EPSILON) ? xo.getDoubleAttribute(EPSILON) : null;
        Integer past = xo.hasAttribute(PAST) ? xo.getIntegerAttribute(PAST) : null;
        Double delta = xo.hasAttribute(DELTA) ? xo.getDoubleAttribute(DELTA) : null;
        // lbfgs4j's own checks (LBFGSERR_INVALID_EPSILON, _TESTPERIOD, _DELTA), made at parse time
        if ((epsilon != null && epsilon < 0) || (past != null && past < 0) || (delta != null && delta < 0)) {
            throw new XMLParseException(EPSILON + ", " + PAST + " and " + DELTA + " must not be negative");
        }

        MaximizerWrtParameter.Settings settings =
                new MaximizerWrtParameter.Settings(nIterations, numericGradientStepSize, initialGuess, printScreen, includeJacobian)
                        .setBudgets(maxSeconds, maxEvaluations, relTol, stallWindow)
                        .setCheckGradientEvery(checkGradientEvery)
                        .setConvergence(epsilon, past, delta);

        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(likelihood, parameter, gradient, transform,
                settings);
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
            new StringAttributeRule(NumericGradientStepSizeProvider.NUMERIC_STEP_SIZE,
                    "a step-size ratio, or tiny, small, medium or large", true),
            AttributeRule.newDoubleRule(MAX_SECONDS, true, "stop after this many seconds (default 0: no limit)"),
            AttributeRule.newLongIntegerRule(MAX_EVALUATIONS, true, "stop after this many objective evaluations (default 0: no limit)"),
            AttributeRule.newDoubleRule(REL_TOL, true, "the relative improvement below which stallWindow iterations count as a stall"),
            AttributeRule.newIntegerRule(STALL_WINDOW, true, "stop when the objective has improved by at most relTol over this many iterations (default 0: off)"),
            AttributeRule.newIntegerRule(CHECK_GRADIENT_EVERY, true, "check the gradient against finite differences every this many evaluations (default 0: never)"),
            AttributeRule.newDoubleRule(EPSILON, true, "lbfgs4j's gradient-norm convergence tolerance"),
            AttributeRule.newIntegerRule(PAST, true, "lbfgs4j's delta-test period in iterations (0: off)"),
            AttributeRule.newDoubleRule(DELTA, true, "lbfgs4j's relative-decrease convergence tolerance"),
    };
}
