/*
 * NumericalGradientTest.java
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

package test.dr.inference.hmc;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.NumericalGradient;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import junit.framework.TestCase;

/**
 * numericalGradient's central differences near a bound: the step stays small relative to the room to the bound, so a
 * positive parameter close to 0 (a clock rate of 7E-8) gets its slope, not a difference over [x / 2, 3x / 2].
 */
public class NumericalGradientTest extends TestCase {

    public NumericalGradientTest(String name) {
        super(name);
    }

    // log density of Gamma(shape 3, scale 0.5): 2 log x - 2x + c, with slope 2 / x - 2
    public void testCloseToTheLowerBound() {
        for (double value : new double[]{7E-8, 1E-3, 0.4, 2.2}) {
            final double exact = 2.0 / value - 2.0;
            assertEquals("x = " + value, exact, gradientAt(value), 1E-6 * Math.abs(exact));
        }
    }

    // at the bound itself, a one-sided difference
    public void testAtTheLowerBound() {
        final double gradient = gradientAt(0.0);
        assertFalse(Double.isNaN(gradient));
    }

    private static double gradientAt(double value) {
        Parameter x = new Parameter.Default("x", new double[]{value});
        x.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        DistributionLikelihood likelihood = new DistributionLikelihood(new GammaDistribution(3.0, 0.5));
        likelihood.addData(x);
        final double gradient = new NumericalGradient(likelihood, x).getGradientLogDensity()[0];
        assertEquals("the parameter is restored", value, x.getParameterValue(0), 0.0);
        return gradient;
    }
}
