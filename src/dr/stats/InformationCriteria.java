/*
 * InformationCriteria.java
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

package dr.stats;

/**
 * Information criteria of a maximum-likelihood fit, on the -2 log-likelihood scale (lower is better), as IQ-TREE
 * computes them (computeInformationScores).
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class InformationCriteria {

    private InformationCriteria() {
    }

    /**
     * @return -2 logL + 2k
     */
    public static double aic(double logLikelihood, int k) {
        return -2.0 * logLikelihood + 2.0 * k;
    }

    /**
     * @return AIC + 2k(k + 1) / (n - k - 1), or NaN when n <= k + 1, where IQ-TREE instead clamps the denominator at 1
     */
    public static double aicc(double logLikelihood, int k, double n) {
        if (!(n - k - 1.0 > 0.0)) {
            return Double.NaN;
        }
        return aic(logLikelihood, k) + 2.0 * k * (k + 1.0) / (n - k - 1.0);
    }

    /**
     * @return -2 logL + k log n
     */
    public static double bic(double logLikelihood, int k, double n) {
        return -2.0 * logLikelihood + k * Math.log(n);
    }
}
