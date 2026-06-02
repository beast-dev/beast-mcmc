/*
 * VectorSliceParameterTest.java
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

package test.dr.inference.model;

import dr.inference.model.Parameter;
import dr.inference.model.VectorSliceParameter;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class VectorSliceParameterTest extends MathTestCase {
    public void testVectorSlice() {

        Parameter pA = new Parameter.Default(new double[] { 0.0, 1.0, 2.0 });
        Parameter pB = new Parameter.Default(new double[] { 3.0, 4.0, 5.0 });

        VectorSliceParameter slice = new VectorSliceParameter("slice",1);
        slice.addParameter(pA);
        slice.addParameter(pB);

        assertEquals(2, slice.getDimension());
        assertEquals(1.0, slice.getParameterValue(0), tolerance);
        assertEquals(4.0, slice.getParameterValue(1), tolerance);

        pA.setParameterValue(1, 41.0);
        assertEquals(41.0, slice.getParameterValue(0), tolerance);

        slice.setParameterValue(1, 44.0);
        assertEquals(44.0, pB.getParameterValue(1), tolerance);

    }
    private static final double tolerance = 1E-6;
}
