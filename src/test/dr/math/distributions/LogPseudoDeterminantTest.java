/*
 * KroneckerOperationTest.java
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

package test.dr.math.distributions;

import dr.math.distributions.GaussianMarkovRandomField;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class LogPseudoDeterminantTest extends MathTestCase {

    public void testPseudoDeterminant() {

        GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q =
                new GaussianMarkovRandomField.SymmetricTriDiagonalMatrix(
                        new double[] { 2.0, 2.0, 2.0, 2.0},
                        new double[] { -0.75, -0.75, -0.75} );

        double logDet1 = GaussianMarkovRandomField.getLogPseudoDeterminantViaDenseEigendecomposition(q);
        double logDet2 = GaussianMarkovRandomField.getLogDeterminantViaRecursion(q);

        assertEquals(logDet1, logDet2, 1E-10);

        // mat <- matrix(c(2, -0.75, 0, 0, -0.75, 2, -0.75, 0, 0, -0.75, 2, -0.75, 0, 0, -0.75, 2), nrow = 4)
        // log(det(mat))
        assertEquals(logDet1, 2.258258, 1E-5);


        q = new GaussianMarkovRandomField.SymmetricTriDiagonalMatrix(
                        new double[] { 1.0, 2.0, 2.0, 1.0},
                        new double[] { -1.0, -1.0, -1.0} );

        logDet1 = GaussianMarkovRandomField.getLogPseudoDeterminantViaDenseEigendecomposition(q);
        logDet2 = GaussianMarkovRandomField.getLogPseudoDeterminantViaTriangularEigenDecomposition(q);

        assertEquals(logDet1, logDet2, 1E-10);

        // mat <- matrix(c(1, -1, 0, 0, -1, 2, -0.75, 0, 0, -0.75, 2, -1, 0, 0, -1, 1), nrow = 4, byrow = TRUE)
        // log(det(mat))
        assertEquals(logDet1,  1.386294, 1E-5);

    }
}
