/*
 * SmoothSkygridTest.java
 *
 * Copyright (c) 2002-2021 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.coalescent;

import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import junit.framework.TestCase;

/**
 * @author Marc A. Suchard
 */
public class SmoothSkygridTest extends TestCase {

    public void testSmoothSkyrideIntegrals() throws Exception {
        double x;
        double y;

        // Linear model case
        x = SmoothSkygridLikelihood.getNumericIntensityInInterval(1.1, 1.6,
                1, 2,
                Math.log(5), Math.log(10), 1.0);
        y = SmoothSkygridLikelihood.getAnalyticIntensityForLinearModel(1.1, 1.6,
                1,2,
                Math.log(5), Math.log(10));
        System.err.println(x + " " + y);
        assertEquals(x, y, tolerance);

        // Approximated step function
        x = SmoothSkygridLikelihood.getPopSizeInInterval(1.6,
                1, 2,
                Math.log(5), Math.log(10), 100.0);
        System.err.println(x);
        y = 10;
        assertEquals(x, y, tolerance);

        x = SmoothSkygridLikelihood.getIntensityInInterval(1, 2,
                1, 2,
                Math.log(5), Math.log(10), 100.0);
        System.err.println(x);
        y = 0.5 * 1 / 5 + 0.5 * 1 / 10;
        assertEquals(x, y, tolerance);
    }

    private final static double tolerance = 1E-4;
}
