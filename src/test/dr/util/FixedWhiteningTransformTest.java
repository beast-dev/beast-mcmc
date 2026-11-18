/*
 * FixedWhiteningTransformTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 */

package test.dr.util;

import dr.inference.model.Parameter;
import dr.util.FixedWhiteningTransform;
import dr.util.Transform;
import junit.framework.TestCase;

/** Focused round-trip and gradient-pullback tests for fixed whitening. */
public final class FixedWhiteningTransformTest extends TestCase {

    public void testRoundTripAndFiniteDifferenceGradient() {
        final Parameter parameter = new Parameter.Default(new double[]{1.2, 0.8, -0.3});
        final Transform.MultivariableTransform nativeTransform =
                new Transform.Array(Transform.NONE, 3, parameter);
        final double[] reference = new double[]{1.0, 1.0, 0.0};
        final double[][] inverseMap = new double[][]{
                {2.0, 0.2, -0.1},
                {0.0, 1.5, 0.3},
                {0.1, -0.2, 0.8},
        };
        final FixedWhiteningTransform transform = new FixedWhiteningTransform(
                parameter, nativeTransform, reference, inverseMap);

        final double[] x = parameter.getParameterValues();
        final double[] z = transform.transform(x, 0, x.length);
        final double[] recovered = transform.inverse(z, 0, z.length);
        for (int i = 0; i < x.length; ++i) {
            assertEquals(x[i], recovered[i], 1.0e-12);
        }

        // f(x) = 0.5 ||x||^2, so df/dx = x.  Validate the transform's
        // analytic pullback against finite differences in whitened z.
        final double[] analytic = transform.updateGradientUnWeightedLogDensity(
                x, x, 0, x.length);
        final double step = 1.0e-7;
        for (int j = 0; j < z.length; ++j) {
            final double[] plus = z.clone();
            final double[] minus = z.clone();
            plus[j] += step;
            minus[j] -= step;
            final double numeric = (quadratic(transform.inverse(plus, 0, plus.length)) -
                    quadratic(transform.inverse(minus, 0, minus.length))) / (2.0 * step);
            assertEquals(numeric, analytic[j], 1.0e-7);
        }
    }

    private static double quadratic(double[] values) {
        double result = 0.0;
        for (double value : values) {
            result += 0.5 * value * value;
        }
        return result;
    }
}
