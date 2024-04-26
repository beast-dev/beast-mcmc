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

import dr.evomodel.coalescent.smooth.OldSmoothSkygridLikelihood;
import junit.framework.TestCase;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

/**
 * @author Marc A. Suchard
 */
public class SmoothSkygridTest extends TestCase {

//    public void testStartUp() throws Exception {
//        long startTime = System.nanoTime();
//        UnivariateRealIntegrator integrator = new RombergIntegrator();
//        integrator.integrate(v -> 2.0 * v + 1.0, 1, 2);
//        long endTime = System.nanoTime();
//        System.err.println("Total startup time: " + (endTime-startTime) + "ns");
//    }

    public void testSmoothSkyrideIntegrals() throws Exception {
        double x;
        double y;
        
        double[] xx;
        double[] yy;

        long startTime;
        long endTime;

        startTime = System.nanoTime();
        UnivariateRealIntegrator integrator = new RombergIntegrator();
        integrator.integrate(v -> 2.0 * v + 1.0, 1, 2);
        endTime = System.nanoTime();
        System.err.println("Total startup time: " + (endTime-startTime) + "ns");

        // Test gradients
        startTime = System.nanoTime();
        xx = OldSmoothSkygridLikelihood.getGradientWrtLogPopSizesInIntervalViaCentralDifference(
                1.1, 1.9,
                1, 2,
                Math.log(5), Math.log(10), 2.0);
        endTime = System.nanoTime();
        System.err.println("Total cDiff execution time: " + (endTime-startTime) + "ns");
        
        startTime = System.nanoTime();
        yy = OldSmoothSkygridLikelihood.getGradientWrtLogPopSizesInInterval(
                1.1, 1.9,
                1, 2,
                Math.log(5), Math.log(10), 2.0);
        endTime = System.nanoTime();
        System.err.println("Total grad execution time: " + (endTime-startTime) + "ns");

        System.err.println(xx[0] + " " + xx[1]);
        System.err.println(yy[0] + " " + yy[1]);
        assertEquals(xx[0], yy[0], tolerance);
        assertEquals(xx[1], yy[1], tolerance);
        
        // New time-dependent refactoring
//        x = SmoothSkygridLikelihood.getNewIntensityInInterval(1.1, 1.6,
//                1,2,
//                Math.log(5), Math.log(10), 2.0);
//        y = SmoothSkygridLikelihood.getIntensityInInterval(1.1, 1.6,
//                1,2,
//                Math.log(5), Math.log(10), 2.0);
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);

//        // Constant model case
//        x = SmoothSkygridLikelihood.getIntensityInInterval(1.1, 1.6,
//                1,2,
//                Math.log(5), Math.log(10), 0.0);
//        y = (1.6 - 1.1) / Math.exp((Math.log(5) + Math.log(10)) / 2);
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);
//
//        // Linear model case
//        x = SmoothSkygridLikelihood.getNumericIntensityInInterval(1.1, 1.6,
//                1, 2,
//                Math.log(5), Math.log(10), 1.0);
//        y = SmoothSkygridLikelihood.getIntensityInInterval(1.1, 1.6,
//                1,2,
//                Math.log(5), Math.log(10), 1.0);
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);
//
//        // True step function
//        x = SmoothSkygridLikelihood.getIntensityInInterval(1.1, 1.4,
//                1, 2,
//                Math.log(5), Math.log(10), Double.POSITIVE_INFINITY);
//        y = (1.4 - 1.1) / 5;
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);
//
//        x = SmoothSkygridLikelihood.getIntensityInInterval(1.1, 1.6,
//                1, 2,
//                Math.log(5), Math.log(10), Double.POSITIVE_INFINITY);
//        y = (1.5 - 1.1) / 5 + (1.6 - 1.5) / 10;
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);
//
//        x = SmoothSkygridLikelihood.getIntensityInInterval(1.6, 1.8,
//                1, 2,
//                Math.log(5), Math.log(10), Double.POSITIVE_INFINITY);
//        y = (1.8 - 1.6) / 10;
//        System.err.println(x + " " + y);
//        assertEquals(x, y, tolerance);
//
        // Approximated step function
        startTime = System.nanoTime();
        x = OldSmoothSkygridLikelihood.getPopSizeInInterval(1.6,
                1, 2,
                Math.log(5), Math.log(10), 100.0);
        System.err.println(x);
        y = 10;
        assertEquals(x, y, tolerance);

        x = OldSmoothSkygridLikelihood.getIntensityInInterval(1, 2,
                1, 2,
                Math.log(5), Math.log(10), 100.0);
        System.err.println(x);
        y = 0.5 * 1 / 5 + 0.5 * 1 / 10;
        assertEquals(x, y, tolerance);
        endTime = System.nanoTime();
        System.err.println("Total execution time: " + (endTime-startTime) + "ns");
    }

    private final static double tolerance = 1E-4;
}
