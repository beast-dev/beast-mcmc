/*
 * MDSTest.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.multidimensionalscaling;

import dr.inference.operators.EllipticalSliceOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class MDSRotationTest extends MathTestCase {


    private static final double tolerance = 1E-6;

    public void testRotation() {
        final int dim = 3;
        final int length = 2;

        MathUtils.setSeed(666);

        double[] x = new double[dim * length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = 10.0 * MathUtils.nextGaussian();
        }

        System.err.println("Starting vector: " + new Vector(x));

        double[] initialLength = new double[length];
        for (int i = 0; i < length; ++i) {
            initialLength[i] = norm(x, i * dim, dim);
            System.err.println("\tlength = " + initialLength[i]);
        }

        EllipticalSliceOperator.transformPoint(x, false, true, dim);

        System.err.println("Ending vector : " + new Vector(x));

        double[] finalLength = new double[length];
        for (int i = 0; i < length; ++i) {
            finalLength[i] = norm(x, i * dim, dim);
            System.err.println("\tlength = " + finalLength[i]);
        }

        assertEquals(initialLength, finalLength, tolerance);

    }

    private double norm(double[] x, int offset, int dim) {
        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += x[offset + i] * x[offset + i];
        }
        return Math.sqrt(sum);
    }

//    public void testInitialization() {
//        if (mds != null) {
//            int i = mds.initialize(2, 100, 0);
//            assertEquals(i, 0);
//            i = mds.initialize(2, 100, 0);
//            assertEquals(i, 1);
//        } else {
//            System.out.println("testInitialization skipped");
//        }
//
//    }
//
//    public void testMakeDirty() {
//        if (mds != null) {
//            mds.makeDirty(0);
//        } else {
//            System.out.println("testMakeDirty skipped");
//        }
//
//    }


}
