/*
 * LKJTransformTest.java
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

package test.dr.util;

import dr.util.*;
import test.dr.math.MathTestCase;

/**
 * @author Marc Suchard
 */

public class UnitSimplexTest extends MathTestCase {

    public void testTransformation() {

        System.out.println("\nTest Unit-Simplex transform.");

        double[] valuesOnReals = new double[]{ 0.0, 0.0, 0.0 };
        double[] valuesOnSimplex = new double[]{ 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

        UnitSimplexToRealsTransform transform = new UnitSimplexToRealsTransform(valuesOnReals.length);
        double[] result1 = transform.inverse(valuesOnReals, 0, valuesOnReals.length);

        assertEquals(result1, valuesOnSimplex, 1E-10);

        double[] result2 = transform.transform(valuesOnSimplex, 0, valuesOnSimplex.length);

        System.out.println("Success");
    }

    public void testReducedDimension() {

        System.out.println("\nTest reduced dimension");

        double[] valuesOnReals = new double[]{ 0.0, 0.0, 10.0 };
        double[] valuesOnSimplex = new double[]{ 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

        UnitSimplexToRealsTransform transform = new UnitSimplexToRealsTransform(valuesOnReals.length);
        double[] result1 = transform.inverse(valuesOnReals, 0, valuesOnReals.length);

        assertEquals(result1, valuesOnSimplex, 1E-10);
    }



}
