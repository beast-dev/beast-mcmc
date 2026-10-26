/*
 * GammaQuantileDerivative.java
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

package dr.math.distributions;

import dr.math.GammaFunction;
import org.apache.commons.math.special.Gamma;

/**
 * Derivatives of gamma quantiles with respect to the shape, for the discretized gamma site-rate model.
 * <p>
 * For the standard gamma (scale one) the u-quantile x(a) solves P(a, x) = u, where P is the regularized lower
 * incomplete gamma function, so dx/da = -(dP/da) / f(x) with f the gamma density. The shape derivative of P uses the
 * series P(a, x) = e^{-x} sum_n x^{a+n} / Gamma(a+n+1), differentiated term by term:
 * dP/da = e^{-x} sum_n x^{a+n} / Gamma(a+n+1) [log x - psi(a+n+1)].
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class GammaQuantileDerivative {

    private static final int MAX_TERMS = 10000000;
    private static final double RELATIVE_TOLERANCE = 1E-17;

    /**
     * @return d q / d shape for the u-quantile q of the mean-one gamma, Gamma(shape, scale = 1 / shape)
     */
    public static double meanOneQuantileDerivative(double u, double shape) {
        final double x = shape * GammaDistribution.quantile(u, shape, 1.0 / shape); // the standard gamma u-quantile
        return standardQuantileDerivative(x, shape) / shape - x / (shape * shape);
    }

    /**
     * @return dx/da with P(a, x) held fixed, where x is a quantile of the standard gamma with shape a
     */
    public static double standardQuantileDerivative(double x, double a) {
        final double logDensity = (a - 1.0) * Math.log(x) - x - GammaFunction.lnGamma(a);
        return -lowerRegularizedShapeDerivative(a, x) / Math.exp(logDensity);
    }

    /**
     * @return d P(a, x) / d a for the regularized lower incomplete gamma function P
     */
    public static double lowerRegularizedShapeDerivative(double a, double x) {
        final double logX = Math.log(x);
        double logTerm = a * logX - x - GammaFunction.lnGamma(a + 1.0); // log of e^{-x} x^{a+n} / Gamma(a+n+1)
        double psi = Gamma.digamma(a + 1.0);                               // psi(a+n+1)

        double sum = 0.0;
        double mass = 0.0;
        for (int n = 0; n < MAX_TERMS; ++n) {
            final double term = Math.exp(logTerm);
            sum += term * (logX - psi);
            mass += term;
            // the terms peak near n = x - a; stop once past the peak and negligible
            if (n > x - a && term <= RELATIVE_TOLERANCE * mass) {
                return sum;
            }
            logTerm += logX - Math.log(a + n + 1.0);
            psi += 1.0 / (a + n + 1.0);
        }
        throw new ArithmeticException("Series for dP(a, x)/da did not converge at a = " + a + ", x = " + x);
    }
}
