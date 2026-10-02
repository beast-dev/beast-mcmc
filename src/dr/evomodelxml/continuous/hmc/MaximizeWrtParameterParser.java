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

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
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
    private static final String LBFGS_EPSILON = "lbfgsEpsilon";
    private static final String LBFGS_PAST = "lbfgsPast";
    private static final String LBFGS_DELTA = "lbfgsDelta";
    private static final String LBFGS_FTOL = "lbfgsFtol";
    private static final String LBFGS_GTOL = "lbfgsGtol";
    private static final String LBFGS_XTOL = "lbfgsXtol";
    private static final String LBFGS_MAX_LINE_SEARCH = "lbfgsMaxLineSearch";
    private static final String FISHER_GRADIENT = "fisherGradient";
    private static final String FISHER_PRIOR_GRADIENT = "fisherPriorGradient";
    private static final String FISHER_REFRESH_ON_FAILURE = "fisherRefreshOnFailure";
    private static final String FISHER_MAXIMUM_REFRESHES = "fisherMaximumRefreshes";
    private static final String FISHER_GRADIENT_THRESHOLD = "fisherGradientThreshold";
    private static final String FISHER_EIGEN_FLOOR_FRACTION = "fisherEigenFloorFraction";
    private static final String FISHER_NULL_PRECISION = "fisherNullPrecision";
    private static final String FISHER_TRUST_RADIUS = "fisherTrustRadius";

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
        Double lbfgsEpsilon = xo.hasAttribute(LBFGS_EPSILON) ? xo.getDoubleAttribute(LBFGS_EPSILON) : null;
        Integer lbfgsPast = xo.hasAttribute(LBFGS_PAST) ? xo.getIntegerAttribute(LBFGS_PAST) : null;
        Double lbfgsDelta = xo.hasAttribute(LBFGS_DELTA) ? xo.getDoubleAttribute(LBFGS_DELTA) : null;
        Double lbfgsFtol = xo.hasAttribute(LBFGS_FTOL) ? xo.getDoubleAttribute(LBFGS_FTOL) : null;
        Double lbfgsGtol = xo.hasAttribute(LBFGS_GTOL) ? xo.getDoubleAttribute(LBFGS_GTOL) : null;
        Double lbfgsXtol = xo.hasAttribute(LBFGS_XTOL) ? xo.getDoubleAttribute(LBFGS_XTOL) : null;
        Integer lbfgsMaxLineSearch = xo.hasAttribute(LBFGS_MAX_LINE_SEARCH) ?
                xo.getIntegerAttribute(LBFGS_MAX_LINE_SEARCH) : null;
        boolean fisherRefreshOnFailure = xo.getAttribute(FISHER_REFRESH_ON_FAILURE, false);
        int fisherMaximumRefreshes = xo.getAttribute(FISHER_MAXIMUM_REFRESHES, 2);
        double fisherGradientThreshold = xo.getAttribute(FISHER_GRADIENT_THRESHOLD, 1.0e-2);
        double fisherEigenFloorFraction = xo.getAttribute(FISHER_EIGEN_FLOOR_FRACTION, 10.0);
        double fisherNullPrecision = xo.getAttribute(FISHER_NULL_PRECISION, 1.0);
        double fisherTrustRadius = xo.getAttribute(FISHER_TRUST_RADIUS, Double.POSITIVE_INFINITY);

        if (gradient != null) {
            parameter = gradient.getParameter();
            likelihood = gradient.getLikelihood();
        } else {
            XMLObject cxo = xo.getChild(DENSITY);
            parameter = (Parameter) cxo.getChild(Parameter.class);
            likelihood = (Likelihood) cxo.getChild(Likelihood.class);
        }

        Transform transform = (Transform) xo.getChild(Transform.class);

        final double numericGradientStepSize = NumericGradientStepSizeProvider.parseStepSizeRatio(xo);

        final GradientWrtParameterProvider fisherGradient = xo.hasChildNamed(FISHER_GRADIENT)
                ? (GradientWrtParameterProvider) xo.getElementFirstChild(FISHER_GRADIENT)
                : null;
        final GradientWrtParameterProvider fisherPriorGradient = xo.hasChildNamed(FISHER_PRIOR_GRADIENT)
                ? (GradientWrtParameterProvider) xo.getElementFirstChild(FISHER_PRIOR_GRADIENT)
                : null;

        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(likelihood, parameter, gradient, transform,
                new MaximizerWrtParameter.Settings(nIterations, numericGradientStepSize, initialGuess, printScreen,
                        includeJacobian, lbfgsEpsilon, lbfgsPast, lbfgsDelta, lbfgsFtol, lbfgsGtol, lbfgsXtol,
                        lbfgsMaxLineSearch, fisherRefreshOnFailure, fisherMaximumRefreshes,
                        fisherGradientThreshold, fisherEigenFloorFraction, fisherNullPrecision,
                        fisherTrustRadius), fisherGradient, fisherPriorGradient);
        maximizer.maximize();

        return maximizer;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Maximizes a density with L-BFGS, optionally rebuilding a rank-reduced " +
                "Fisher whitening and restarting with fresh curvature history after selected line-search failures.";
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
            new ElementRule(FISHER_GRADIENT,
                    new XMLSyntaxRule[]{new ElementRule(GradientWrtParameterProvider.class)}, true),
            new ElementRule(FISHER_PRIOR_GRADIENT,
                    new XMLSyntaxRule[]{new ElementRule(GradientWrtParameterProvider.class)}, true),
            AttributeRule.newDoubleRule(N_ITERATIONS, true),
            AttributeRule.newBooleanRule(INITIAL_GUESS, true),
            AttributeRule.newBooleanRule(PRINT_SCREEN, true),
            AttributeRule.newBooleanRule(INCLUDE_JACOBIAN, true),
            AttributeRule.newDoubleRule(LBFGS_EPSILON, true),
            AttributeRule.newIntegerRule(LBFGS_PAST, true),
            AttributeRule.newDoubleRule(LBFGS_DELTA, true),
            AttributeRule.newDoubleRule(LBFGS_FTOL, true),
            AttributeRule.newDoubleRule(LBFGS_GTOL, true),
            AttributeRule.newDoubleRule(LBFGS_XTOL, true),
            AttributeRule.newIntegerRule(LBFGS_MAX_LINE_SEARCH, true),
            AttributeRule.newBooleanRule(FISHER_REFRESH_ON_FAILURE, true),
            AttributeRule.newIntegerRule(FISHER_MAXIMUM_REFRESHES, true),
            AttributeRule.newDoubleRule(FISHER_GRADIENT_THRESHOLD, true),
            AttributeRule.newDoubleRule(FISHER_EIGEN_FLOOR_FRACTION, true),
            AttributeRule.newDoubleRule(FISHER_NULL_PRECISION, true),
            AttributeRule.newDoubleRule(FISHER_TRUST_RADIUS, true),
    };
}
