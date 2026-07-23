/*
 * EmpiricalBayesPoissonSmootherTest.java
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

import dr.math.EmpiricalBayesPoissonSmoother;
import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;

/**
 * @author Marc A. Suchard
 */
public class EmpiricalBayesPoissonSmootherTest extends MathTestCase {

    public void testSmootherMean() {

        MathUtils.setSeed(666);
        int length = 100;

        double[] in = new double[length];
        for (int i = 0; i < length; ++i) {
            in[i] = MathUtils.nextDouble() * 100.0;
        }

        doSmoothing(in);
    }

    public void testUnstableMean() {
        double[] in = new double[] {0, 1};
        doSmoothing(in);
    }

    private void doSmoothing(double[] in) {

        double meanX = DiscreteStatistics.mean(in);
        double varX = DiscreteStatistics.variance(in);
        System.err.println("Original mean: " + meanX);
        System.err.println("Original var : " + varX + "\n");

        double[] out = EmpiricalBayesPoissonSmoother.smooth(in);

        double meanY = DiscreteStatistics.mean(out);
        double varY = DiscreteStatistics.variance(out);
        System.err.println("Smoothed mean: " + meanY);
        System.err.println("Smoothed var : " + varY);
        System.err.println("");

        assertEquals(meanX, meanY, tolerance);
        assertTrue(varY <= varX);
    }

    private static final double tolerance = 10E-6;


}
