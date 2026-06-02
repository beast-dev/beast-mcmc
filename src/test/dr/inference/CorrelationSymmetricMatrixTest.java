/*
 * CorrelationSymmetricMatrixTest.java
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

package test.dr.inference;

import dr.inference.model.CorrelationSymmetricMatrix;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import test.dr.math.MathTestCase;

public class CorrelationSymmetricMatrixTest extends MathTestCase {

    public void testGetAttributeValue() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_CORRELATION);
        double[] values = matrix.getAttributeValue();

        assertEquals(values.length, matrix.getRowDimension() * matrix.getColumnDimension());
    }

    public void testGetParameterValue() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_CORRELATION);
        double value = matrix.getParameterValue(0, 0);
        assertEquals(value, 1.0);

        double top = matrix.getParameterValue(2, 1);
        double bottom = matrix.getParameterValue(1, 2);
        assertEquals(top, bottom);

        double test = 4.0 * Math.sqrt(2 * 3);
        assertEquals(top, test);
    }

    public void testGetParameterAsMatrix() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_IS);
        double[][] M = matrix.getParameterAsMatrix();
        System.out.println(new Matrix(M));

        double sum = 0.0;
        for (double x : M[0]) {
            sum += x;
        }

        assertEquals(sum, 7.0);

    }

    private CorrelationSymmetricMatrix getMatrix(CorrelationSymmetricMatrix.Type type) {
            Parameter diagonals = new Parameter.Default(new double[] { 1.0, 2.0, 3.0, 4.0 });
            Parameter offDiagonals = new Parameter.Default(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

            return new CorrelationSymmetricMatrix(diagonals, offDiagonals, type);
    }
}