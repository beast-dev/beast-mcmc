/*
 * SimplexTangent.java
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

package dr.inference.hmc;

/**
 * Gradients of functions of x on the simplex {x > 0, sum x = S}. Only the tangent part is defined there:
 * gradients that differ by c * 1 describe the same function on the simplex.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SimplexTangent {

    private SimplexTangent() {
    }

    /**
     * @return t_k = g_k - (x.g) / S, the directional derivative along e_k - x / S, which keeps sum x = S
     */
    public static double[] project(double[] gradient, double[] x, double simplexSum) {
        double inner = 0.0;
        for (int i = 0; i < x.length; ++i) {
            inner += x[i] * gradient[i];
        }
        double[] tangent = new double[gradient.length];
        for (int i = 0; i < gradient.length; ++i) {
            tangent[i] = gradient[i] - inner / simplexSum;
        }
        return tangent;
    }
}
