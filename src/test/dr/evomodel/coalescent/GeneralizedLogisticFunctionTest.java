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
import dr.math.GeneralizedLogisticFunction;
import junit.framework.TestCase;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

/**
 * @author Marc A. Suchard
 */
public class GeneralizedLogisticFunctionTest extends TestCase {

//    public void testStartUp() throws Exception {
//        long startTime = System.nanoTime();
//        UnivariateRealIntegrator integrator = new RombergIntegrator();
//        integrator.integrate(v -> 2.0 * v + 1.0, 1, 2);
//        long endTime = System.nanoTime();
//        System.err.println("Total startup time: " + (endTime-startTime) + "ns");
//    }

    public void testGeneralizedLogisticFunction() throws Exception {

        double y;
        y = GeneralizedLogisticFunction.evaluate(1,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(2,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.1,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.9,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.4,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.6,
                1,2,
                10,20,
                10,1.5, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.4,
                1,2,
                10,20,
                10,1.6, 1);
        System.err.println(y);

        y = GeneralizedLogisticFunction.evaluate(1.6,
                1,2,
                10,20,
                10,1.6, 1);
        System.err.println(y);

//        assertEquals(xx[0], yy[0], tolerance);
    }

    private final static double tolerance = 1E-4;
}
