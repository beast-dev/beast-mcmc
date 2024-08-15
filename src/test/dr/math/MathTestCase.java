/*
 * MathTestCase.java
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

package test.dr.math;

import dr.math.matrixAlgebra.Vector;
import junit.framework.TestCase;

/**
 * @author Marc A. Suchard
 */
public abstract class MathTestCase extends TestCase {

    public MathTestCase() {        
    }

    public MathTestCase(String name) {
        super(name);
    }

    protected void assertEquals(double[][] a, double[][] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; ++i) {
            assertEquals(a[i], b[i], accuracy);
        }
    }

    protected void assertEquals(double[] a, double[] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], accuracy);
        }
    }

    protected void printSquareMatrix(double[] A, int dim) {
        double[] row = new double[dim];
        for (int i = 0; i < dim; i++) {
            System.arraycopy(A, i * dim, row, 0, dim);
            System.out.println(new Vector(row));
        }
    }

    protected static double accumulate(double[] vector) {
        double total = 0;
        for (double x : vector) {
            total += x;
        }
        return total;
    }
}
