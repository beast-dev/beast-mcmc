/*
 * KroneckerOperation.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math;

import java.util.Arrays;

/**
 * A class that implements the Kronecker product and Kronecker sum matrix operations
 *
 * @author Marc A. Suchard
 */

public class KroneckerOperation {

    // Computes the Kronecker sum C of A (m-by-1) and B (1-by-n)
    // C = A %x% I_n + I_m %x% B
    public static double[] sum(double[] A, double[] B) {
        final int dim = A.length * B.length;

        double[] out = new double[dim];
        double[] tmpA = new double[dim];
        double[] tmpB = new double[dim];
        double[] Im = makeIdentityVector(A.length);
        double[] In = makeIdentityVector(B.length);

        sum(A, B, Im, In, tmpA, tmpB, out);
        return out;
    }

    public static void sum(double[] A, double[] B,
                           double[] Im, double[] In,
                           double[] tmpA, double[] tmpB,
                           double[] C) {
        final int m = A.length;
        final int n = B.length;
        final int dim = m * n;

        if (C.length != dim
                || Im.length != m || In.length != n
                || tmpA.length != dim || tmpB.length != dim) {
            throw new RuntimeException("Wrong dimensions in Kronecker sum");
        }

        product(A, m, 1, In, n, 1, tmpA);
        product(Im, m, 1, B, n, 1, tmpB);

        for (int i = 0; i < dim; i++) {
            C[i] = tmpA[i] + tmpB[i];
        }
    }

    // Computes the Kronecker sum C of A (m-by-m) and B (n-by-n)
    // C = A %x% I_n + I_m %x% B
    public static double[] sum(double[] A, int m, double[] B, int n) {
        final int dim = m * n;

        double[] out = new double[dim * dim];
        double[] tmpA = new double[dim * dim];
        double[] tmpB = new double[dim * dim];
        double[] Im = makeIdentityMatrix(m);
        double[] In = makeIdentityMatrix(n);

        sum(A, m, B, n, Im, In, tmpA, tmpB, out);
        return out;
    }

    public static void sum(double[] A, int m, double[] B, int n,
                           double[] Im, double[] In,
                           double[] tmpA, double[] tmpB,
                           double[] C) {

        final int dim = m * n;

        if (C.length != dim * dim || A.length != m * m || B.length != n * n
                || Im.length != m * m || In.length != n * n
                || tmpA.length != dim * dim || tmpB.length != dim * dim) {
            throw new RuntimeException("Wrong dimensions in Kronecker sum");
        }

        product(A, m, m, In, n, n, tmpA);
        product(Im, m, m, B, n, n, tmpB);

        for (int i = 0; i < dim * dim; i++) {
            C[i] = tmpA[i] + tmpB[i];
        }
    }

    private static double[] makeIdentityMatrix(int dim) {
        double[] out = new double[dim * dim];
        for (int i = 0; i < dim; i++) {
            out[i * dim + i] = 1.0;
        }
        return out;
    }

    private static double[] makeIdentityVector(int dim) {
        double[] out = new double[dim];
        Arrays.fill(out, 1.0);
        return out;
    }

    // Computes the Kronecker product of A (m-by-n) and B (p-by-q)
    public static double[] product(double[] A, int m, int n, double[] B, int p, int q) {
        double[] out = new double[m * p * n * q];
        product(A, m, n, B, p, q, out);
        return out;
    }

    public static void product(double[] A, int m, int n, double[] B, int p, int q, double[] out) {

        final int dimi = m * p;
        final int dimj = n * q;

        if (out.length != dimi * dimj || A.length != m * n || B.length != p * q) {
            throw new RuntimeException("Wrong dimensions in Kronecker product");
        }

        for (int i = 0; i < m; i++) { // 1,\ldots,m

            final int iOffset = i * p;

            for (int j = 0; j < n; j++) { // 1,\ldots,n

                final int jOffset = j * q;

                final double aij = A[i * n + j];

                for (int k = 0; k < p; k++) { // 1,\ldots,p
                    for (int l = 0; l < q; l++) { // 1,\ldots,q
                        out[(iOffset + k) * dimj + jOffset + l] = aij * B[k * q + l];
                    }
                }
            }
        }
    }

    public static double[][] product(double[][] A, double[][] B) {
        final int m = A.length;
        final int n = A[0].length;
        final int p = B.length;
        final int q = B[0].length;

        double[][] out = new double[m * p][n * q];
        product(A, B, out);
        return out;
    }

    public static void product(double[][] A, double[][] B, double[][] out) {
        final int m = A.length;
        final int n = A[0].length;
        final int p = B.length;
        final int q = B[0].length;

        if (out == null || out.length != m * p || out[0].length != n * q) {
            throw new RuntimeException("Wrong dimensions in Kronecker product");
        }

        for (int i = 0; i < m; i++) {
            final int iOffset = i * p;
            for (int j = 0; j < n; j++) {
                final int jOffset = j * q;
                final double aij = A[i][j];

                for (int k = 0; k < p; k++) {
                    for (int l = 0; l < q; l++) {
                        out[iOffset + k][jOffset + l] = aij * B[k][l];
                    }
                }

            }
        }        
    }

    public static double[] vectorize(double[][] A) {
        double[] out = new double[A.length * A[0].length];
        vectorize(A, out);
        return out;
    }

    public static void vectorize(double[][] A, double[] out) {
        final int m = A.length;
        final int n = A[0].length;

        if (out == null || out.length != m * n) {
            throw new RuntimeException("Wrong dimensions in vectorize");
        }

        int offset = 0;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                out[offset] = A[i][j];
                ++offset;
            }
        }
    }
}
