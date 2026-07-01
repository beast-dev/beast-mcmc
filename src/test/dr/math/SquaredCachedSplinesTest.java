/*
 * IntegratedTransformedSplinesTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.math;

import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import junit.framework.TestCase;

public class SquaredCachedSplinesTest extends TestCase {

    private static final double[] KNOTS = new double[]{0.25, 0.5, 0.75};
    private static final double LOWER = 0.0;
    private static final double UPPER = 1.0;
    private static final int DEGREE = 3;
    private static final double TOLERANCE = 5.0e-7;

    public void testSquaredCoefficientGradientMatchesFiniteDifferenceInsideBoundary() {
        checkCoefficientGradient(0.1, 0.9);
    }

    public void testSquaredCoefficientGradientMatchesFiniteDifferencePastBoundary() {
        checkCoefficientGradient(0.35, 1.35);
    }

    public void testSquaredInterceptGradientMatchesFiniteDifference() {
        Parameter coefficients = coefficients();
        Parameter intercept = new Parameter.Default(0.21);
        SquaredCachedSplines splines = squaredSplines(coefficients, intercept);

        double start = 0.2;
        double end = 1.25;
        double analytic = splines.getGradientWrtIntercept(start, end);
        double numerical = finiteDifference(intercept, 0, splines, start, end);

        assertEquals(numerical, analytic, TOLERANCE);
    }

    private void checkCoefficientGradient(double start, double end) {
        Parameter coefficients = coefficients();
        Parameter intercept = new Parameter.Default(0.21);
        SquaredCachedSplines splines = squaredSplines(coefficients, intercept);

        for (int i = 0; i < coefficients.getDimension(); i++) {
            double analytic = splines.getGradientWrtCoefficient(start, end, i);
            double numerical = finiteDifference(coefficients, i, splines, start, end);
            assertEquals(numerical, analytic, TOLERANCE);
        }
    }

    private static SquaredCachedSplines squaredSplines(Parameter coefficients, Parameter intercept) {
        return new SquaredCachedSplines(coefficients, intercept, KNOTS, LOWER, UPPER, DEGREE);
    }

    private static Parameter coefficients() {
        return new Parameter.Default(new double[]{0.13, -0.27, 0.41, -0.11, 0.08, 0.29});
    }

    private static double finiteDifference(Parameter parameter,
                                           int index,
                                           SquaredCachedSplines splines,
                                           double start,
                                           double end) {
        double h = 1.0e-6;
        double original = parameter.getParameterValue(index);

        parameter.setParameterValue(index, original + h);
        double plus = splines.getIntegral(start, end);

        parameter.setParameterValue(index, original - h);
        double minus = splines.getIntegral(start, end);

        parameter.setParameterValue(index, original);

        return (plus - minus) / (2.0 * h);
    }
}
