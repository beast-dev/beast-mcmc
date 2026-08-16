/*
 * ConjugateGradientSolver.java
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

package dr.math;

/**
 * Solves the symmetric positive-definite linear system A x = b via the
 * conjugate-gradient method.  A is supplied as a flat, row-major double[]
 * of length n*n; b is a double[] of length n.
 *
 * @author Marc A. Suchard
 */
public class ConjugateGradientSolver {

    private final int maxIterations;
    private final double tolerance;

    public ConjugateGradientSolver(int maxIterations, double tolerance) {
        this.maxIterations = maxIterations;
        this.tolerance = tolerance;
    }

    /**
     * @param A   symmetric positive-definite matrix, row-major, length n*n
     * @param b   right-hand side vector, length n
     * @return    solution x such that A x ≈ b
     */
    public double[] solve(double[] A, double[] b) {

        final int n = b.length;
        final double[] x = new double[n];
        final double[] r = new double[n];
        final double[] p = new double[n];
        final double[] Ap = new double[n];

        // r = b - A*x0  (x0 = 0)
        System.arraycopy(b, 0, r, 0, n);
        System.arraycopy(r, 0, p, 0, n);

        double rsold = dot(r, r);

        for (int iter = 0; iter < maxIterations; ++iter) {

            if (Math.sqrt(rsold) < tolerance) {
                break;
            }

            multiply(A, p, Ap, n);

            final double alpha = rsold / dot(p, Ap);

            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            final double rsnew = dot(r, r);
            final double beta = rsnew / rsold;

            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * p[i];
            }

            rsold = rsnew;
        }

        return x;
    }

    private static void multiply(double[] A, double[] v, double[] out, int n) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            final int row = i * n;
            for (int j = 0; j < n; ++j) {
                sum += A[row + j] * v[j];
            }
            out[i] = sum;
        }
    }

    private static double dot(double[] u, double[] v) {
        double sum = 0.0;
        for (int i = 0; i < u.length; ++i) {
            sum += u[i] * v[i];
        }
        return sum;
    }
}
