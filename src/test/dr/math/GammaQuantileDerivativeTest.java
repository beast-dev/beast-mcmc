/*
 * GammaQuantileDerivativeTest.java
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
package test.dr.math;

import dr.math.distributions.GammaQuantileDerivative;
import junit.framework.TestCase;
import org.apache.commons.math.special.Gamma;

/**
 * Shape derivatives of the regularized incomplete gamma and of the mean-one gamma quantile, against central
 * differences of a high-precision reference (commons-math P at 1e-16, quantiles by bisection).
 * BEAST's own GammaFunction.incompleteGammaP (AS32, accuracy 1e-8) is too coarse for this.
 */
public class GammaQuantileDerivativeTest extends TestCase {

    private static final double[] SHAPES = {0.05, 0.3, 1.0, 2.0, 5.0, 50.0};

    public GammaQuantileDerivativeTest(String name) {
        super(name);
    }

    public void testLowerRegularizedShapeDerivative() throws Exception {
        for (double a : SHAPES) {
            for (double u : new double[]{0.0625, 0.375, 0.625, 0.9375}) {
                final double x = preciseStandardQuantile(u, a);
                final double h = 1E-5 * a;
                final double numeric = (p(a + h, x) - p(a - h, x)) / (2 * h);
                final double analytic = GammaQuantileDerivative.lowerRegularizedShapeDerivative(a, x);
                assertEquals("a=" + a + " u=" + u, numeric, analytic, 1E-7 * Math.max(1.0, Math.abs(numeric)));
            }
        }
    }

    public void testMeanOneQuantileDerivative() throws Exception {
        for (int categories : new int[]{4, 8}) {
            for (double a : SHAPES) {
                for (int i = 0; i < categories; i++) {
                    final double u = (2.0 * i + 1.0) / (2.0 * categories);
                    final double h = 1E-5 * a;
                    final double numeric = (preciseStandardQuantile(u, a + h) / (a + h)
                            - preciseStandardQuantile(u, a - h) / (a - h)) / (2 * h);
                    final double analytic = GammaQuantileDerivative.meanOneQuantileDerivative(u, a);
                    // the derivative is a difference of two O(0.1) terms, so allow 1e-8 absolute near zero
                    assertEquals("K=" + categories + " a=" + a + " u=" + u, numeric, analytic,
                            1E-6 * Math.max(1E-2, Math.abs(numeric)));
                }
            }
        }
    }

    private static double p(double a, double x) throws Exception {
        return Gamma.regularizedGammaP(a, x, 1E-16, 1000000);
    }

    // bisection on log x, so tiny quantiles at small shapes stay accurate
    private static double preciseStandardQuantile(double u, double a) throws Exception {
        double lower = -800.0;
        double upper = Math.log(a + 60.0 * Math.sqrt(a) + 60.0);
        for (int i = 0; i < 200; i++) {
            final double middle = 0.5 * (lower + upper);
            if (p(a, Math.exp(middle)) < u) {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        return Math.exp(0.5 * (lower + upper));
    }
}
