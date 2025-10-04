/*
 * NormalDistributionTest.java
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

package test.dr.distibutions;

import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.matrix.SparseSquareUpperTriangular;
import dr.stats.DiscreteStatistics;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 * 
 */
public class MultivariateNormalSamplingTest extends MathTestCase {

    public void setUp() { }

    public void testSparseBackSolve() {
        // first
        double[][] u = new double[][] { {1, 0, 0}, {0, 1, 0}, {0, 0, 2} };
        SparseSquareUpperTriangular U = new SparseSquareUpperTriangular(u);
        assertEquals(U.getEntryCount(), 3);

        double[] x = new double[] { 1, 2, 4 };
        double[] b = new double[3];

        scramble(b);
        U.multiplyInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(b, new double[] { 1, 2, 8 }, 0.0);

        scramble(x);
        U.backSolveInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(x, new double[] { 1, 2, 4}, 0.0);

        // second
        u = new double[][] { {1, 0, 0}, {0, 1, 1}, {0, 0, 2} };
        U = new SparseSquareUpperTriangular(u);
        assertEquals(U.getEntryCount(), 4);

        scramble(b);
        U.multiplyInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(b, new double[] { 1, 6, 8 }, 0.0);

        scramble(x);
        U.backSolveInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(x, new double[] { 1, 2, 4}, 0.0);

        // third
        u = new double[][] { {1, 0, -3}, {0, 1, 1}, {0, 0, 2} };
        U = new SparseSquareUpperTriangular(u);
        assertEquals(U.getEntryCount(), 5);

        scramble(b);
        U.multiplyInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(b, new double[] { -11, 6, 8 }, 0.0);

        scramble(x);
        U.backSolveInPlaceMatrixVector(x, 0, b, 0);
        assertEquals(x, new double[] { 1, 2, 4}, 0.0);
    }

    private void scramble(double[] x) {
        for (int i = 0; i < x.length; ++i) {
            x[i] = MathUtils.nextDouble();
        }
    }

    private double[][] getL(double[][] matrix) {
        CholeskyDecomposition chol = null;
        try {
            chol = new CholeskyDecomposition(matrix);
        } catch (IllegalDimension e) {
            fail();
        }

        return chol.getL();
    }

    private SparseSquareUpperTriangular getU(double[][] matrix) {
        double[][] L = getL(matrix);

        double[][] tmp = new double[L.length][L.length];
        for (int i = 0; i < L.length; ++i) {
            for (int j = 0; j < L.length; ++j) {
                tmp[i][j] = L[j][i];
            }
        }

        return new SparseSquareUpperTriangular(tmp);
    }

    public void testSampling() {

        double[][] precision = new double[][] { {1.0, -0.5, 0.0}, {-0.5, 1.0, -0.5}, {0.0, -0.5, 2} };
        SparseSquareUpperTriangular U = getU(precision);

        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
        double[][] L = getL(variance);

        double[] mean = new double[] { 1, 2, 3 };
        final int length = 1000000;

        double[] x = new double[3];
        double[] epsilon = new double[3];

        double scale = 0.5;

        // via variance
        MathUtils.setSeed(666);
        double[][] draws = new double[3][length];
        for (int i = 0; i < length; ++i) {
            MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, 0,
                    L, Math.sqrt(scale), x, 0, epsilon);
            for (int j = 0; j < 3; ++j) {
                draws[j][i] = x[j];
            }
        }

        assertEquals(DiscreteStatistics.mean(draws[0]), 1, 1E-2);
        assertEquals(DiscreteStatistics.mean(draws[1]), 2, 1E-2);
        assertEquals(DiscreteStatistics.mean(draws[2]), 3, 1E-2);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double v = DiscreteStatistics.covariance(draws[i], draws[j],
                        DiscreteStatistics.mean(draws[i]), DiscreteStatistics.mean(draws[j]),
                        1, 1);
                assertEquals(v, scale * variance[i][j], 1E-2);
            }
        }

        // via precision
        MathUtils.setSeed(666);
        draws = new double[3][length];
        for (int i = 0; i < length; ++i) {
            MultivariateNormalDistribution.nextMultivariateNormalViaBackSolvePrecision(mean, 0,
                    U, Math.sqrt(scale), x, 0, epsilon);
            for (int j = 0; j < 3; ++j) {
                draws[j][i] = x[j];
            }
        }

        assertEquals(DiscreteStatistics.mean(draws[0]), 1, 1E-2);
        assertEquals(DiscreteStatistics.mean(draws[1]), 2, 1E-2);
        assertEquals(DiscreteStatistics.mean(draws[2]), 3, 1E-2);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double v = DiscreteStatistics.covariance(draws[i], draws[j],
                        DiscreteStatistics.mean(draws[i]), DiscreteStatistics.mean(draws[j]),
                        1, 1);
                assertEquals(v, scale * variance[i][j], 1E-2);
            }
        }
    }
}
