/*
 * BinomialLikelihood.java
 *
 * Copyright (C) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc A. Suchard
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.distribution;

import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 *
 * @author Marc A. Suchard
 */

public class ExponentialTiltedStabledDistributionTest extends TestCase {

    public ExponentialTiltedStabledDistributionTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(666);
    }

    private static final int sampleLength = (int) 1E6;

    double[] exponents = new double[] {0.5, 0.1, 0.5, 0.5, 0.9};
    double[] tilts = new double[] {2.0, 2.0, 0.75, 17.0, 17.0};
    double[] expectations = new double[] {0.3535, 0.0535, 0.5776, 0.1212, 0.6779}; // R:  copula::retstable()
    double[] variances = new double[] {0.0884,0.0241,0.3848, 0.0036, 0.0040};

    public void testSampler() {

        for (int setting = 0; setting < exponents.length; ++setting) {
            doTestSample(exponents[setting], tilts[setting],
                    expectations[setting], variances[setting]);
        }
    }

    private void doTestSample(double exponent, double tilt,
                              double expectation, double variance) {

        double sum = 0.0;
        double sumSquared = 0.0;
        for (int i = 0; i < sampleLength; ++i) {
            double draw = ExponentialTiltedStableDistribution.nextTiltedStable(exponent, tilt);
            sum += draw;
            sumSquared += draw * draw;
        }

        double mean = sum /sampleLength;
        double var = sumSquared / sampleLength - mean * mean;

        System.out.println(exponent + " " + tilt + " " + mean + " " + var);

        assertEquals(expectation, mean, 1E-3);
        assertEquals(variance, var,1E-2);

    }

    public static Test suite() {
        return new TestSuite(ExponentialTiltedStabledDistributionTest.class);
    }
}