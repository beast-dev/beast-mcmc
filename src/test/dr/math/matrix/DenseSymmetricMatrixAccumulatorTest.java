/*
 * DenseSymmetricMatrixAccumulatorTest.java
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

package test.dr.math.matrix;

import dr.math.matrix.DenseSymmetricMatrixAccumulator;
import junit.framework.TestCase;

/**
 * Tests for DenseSymmetricMatrixAccumulator.
 */
public class DenseSymmetricMatrixAccumulatorTest extends TestCase {

    private static final double TOL = 1e-14;

    public void testAddAndExtract() {
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(3);

        G.add(0, 0, 4.0);
        G.add(1, 0, 2.0);  // lower triangle entry [1][0]
        G.add(0, 1, 0.5);  // same entry [1][0] via upper-to-lower normalisation: adds 0.5
        G.add(1, 1, 3.0);
        G.add(2, 0, 1.0);
        G.add(2, 1, 0.0);
        G.add(2, 2, 5.0);

        // Dense extraction via addScaledToDense
        double[][] dense = new double[3][3];
        G.addScaledToDense(1.0, dense);

        assertEquals(4.0,  dense[0][0], TOL);
        assertEquals(2.5,  dense[1][0], TOL); // 2.0 + 0.5
        assertEquals(2.5,  dense[0][1], TOL); // symmetric
        assertEquals(3.0,  dense[1][1], TOL);
        assertEquals(1.0,  dense[2][0], TOL);
        assertEquals(1.0,  dense[0][2], TOL);
        assertEquals(0.0,  dense[2][1], TOL);
        assertEquals(5.0,  dense[2][2], TOL);
    }

    public void testQuadraticForm() {
        // G = [[2, 1], [1, 3]], x = [1, 2]
        // x'Gx = 2*1 + 2*1*2 + 3*4 = 2 + 4 + 12 = 18
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(2);
        G.add(0, 0, 2.0);
        G.add(1, 0, 1.0);
        G.add(1, 1, 3.0);

        double[] x = {1.0, 2.0};
        assertEquals(18.0, G.quadraticForm(x), TOL);
    }

    public void testAddMatVec() {
        // G = [[2, 1], [1, 3]], x = [1, 2], scale = 1.0
        // Gx = [2+2, 1+6] = [4, 7]
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(2);
        G.add(0, 0, 2.0);
        G.add(1, 0, 1.0);
        G.add(1, 1, 3.0);

        double[] x   = {1.0, 2.0};
        double[] y   = new double[2];
        G.addMatVec(x, 1.0, y);
        assertEquals(4.0, y[0], TOL);
        assertEquals(7.0, y[1], TOL);
    }

    public void testAddMatVecScale() {
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(2);
        G.add(0, 0, 2.0);
        G.add(1, 0, 1.0);
        G.add(1, 1, 3.0);

        double[] x = {1.0, 2.0};
        double[] y = new double[2];
        G.addMatVec(x, -2.0, y);
        assertEquals(-8.0, y[0], TOL);
        assertEquals(-14.0, y[1], TOL);
    }

    public void testZeroMatrix() {
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(3);
        double[] x = {1.0, 2.0, 3.0};
        assertEquals(0.0, G.quadraticForm(x), TOL);
        double[] y = new double[3];
        G.addMatVec(x, 1.0, y);
        for (double v : y) assertEquals(0.0, v, TOL);
    }

    public void testDimension() {
        assertEquals(5, new DenseSymmetricMatrixAccumulator(5).dimension());
    }

    public void testAddScaledToDenseScale() {
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(2);
        G.add(0, 0, 1.0);
        G.add(1, 0, 2.0);
        G.add(1, 1, 3.0);

        double[][] dense = new double[2][2];
        G.addScaledToDense(3.0, dense);
        assertEquals(3.0, dense[0][0], TOL);
        assertEquals(6.0, dense[1][0], TOL);
        assertEquals(6.0, dense[0][1], TOL);
        assertEquals(9.0, dense[1][1], TOL);
    }

    public void testQuadraticFormMatchesBruteForce() {
        // Random-ish 4×4 symmetric matrix
        DenseSymmetricMatrixAccumulator G = new DenseSymmetricMatrixAccumulator(4);
        double[][] M = {
            { 5.0,  1.5, -0.3,  2.0},
            { 1.5,  4.0,  0.7, -1.0},
            {-0.3,  0.7,  3.0,  0.5},
            { 2.0, -1.0,  0.5,  6.0}
        };
        // Add lower triangle to accumulator
        for (int i = 0; i < 4; i++)
            for (int j = 0; j <= i; j++)
                G.add(i, j, M[i][j]);

        double[] x = {1.0, -0.5, 2.0, 0.3};

        // Brute force
        double expected = 0.0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                expected += x[i] * M[i][j] * x[j];

        assertEquals(expected, G.quadraticForm(x), 1e-12);
    }
}
