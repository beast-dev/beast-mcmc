/*
 * TransformedMultivariateDistributionTest.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.math.distributions;

import dr.math.distributions.LKJCholeskyCorrelationDistribution;
import dr.math.distributions.LKJCorrelationDistribution;
import dr.math.distributions.TransformedMultivariateDistribution;
import dr.util.CorrelationToCholesky;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * @author Paul Bastide
 */

public class TransformedMultivariateDistributionTest extends TraceCorrelationAssert {

    private int dimTrait;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public TransformedMultivariateDistributionTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        dimTrait = 5;

        format.setMaximumFractionDigits(5);
    }

    public void test1() {
        System.out.println("\nTest Transformed Multivariate Distribution shape=1.0:");

        double shape = 1.0;

        LKJCholeskyCorrelationDistribution LKJCholDistribution = new LKJCholeskyCorrelationDistribution(dimTrait, shape);
        CorrelationToCholesky transform = new CorrelationToCholesky(dimTrait);
        TransformedMultivariateDistribution transformedDistribution = new TransformedMultivariateDistribution(LKJCholDistribution, transform);

        LKJCorrelationDistribution LKJCorDistribution = new LKJCorrelationDistribution(dimTrait, shape);

        double[] x = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};

        double value = LKJCorDistribution.logPdf(x);
        double valuebis = transformedDistribution.logPdf(x);

        assertEquals("value", format.format(value), format.format(valuebis));

        double[] gradient = LKJCorDistribution.getGradientLogDensity(x);
        double[] gradientbis = transformedDistribution.getGradientLogDensity(x);

        assertEquals("length", gradient.length, gradientbis.length);
        for (int i = 0; i < gradient.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient[i]),
                    format.format(gradientbis[i]));
        }

        double[] x2 = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        double value2 = LKJCorDistribution.logPdf(x);
        double value2bis = transformedDistribution.logPdf(x);

        assertEquals("value", format.format(value2), format.format(value2bis));

        double[] gradient2 = LKJCorDistribution.getGradientLogDensity(x);
        double[] gradient2bis = transformedDistribution.getGradientLogDensity(x);

        assertEquals("length", gradient2.length, gradient2bis.length);
        for (int i = 0; i < gradient2.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient2[i]),
                    format.format(gradient2bis[i]));
        }
    }

    public void test2() {
        System.out.println("\nTest Transformed Multivariate Distribution shape=2.5:");

        double shape = 2.5;

        LKJCholeskyCorrelationDistribution LKJCholDistribution = new LKJCholeskyCorrelationDistribution(dimTrait, shape);
        CorrelationToCholesky transform = new CorrelationToCholesky(dimTrait);
        TransformedMultivariateDistribution transformedDistribution = new TransformedMultivariateDistribution(LKJCholDistribution, transform);

        LKJCorrelationDistribution LKJCorDistribution = new LKJCorrelationDistribution(dimTrait, shape);

        double[] x = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};

        double value = LKJCorDistribution.logPdf(x);
        double valuebis = transformedDistribution.logPdf(x);

        assertEquals("value", format.format(value), format.format(valuebis));

        double[] gradient = LKJCorDistribution.getGradientLogDensity(x);
        double[] gradientbis = transformedDistribution.getGradientLogDensity(x);

        assertEquals("length", gradient.length, gradientbis.length);
        for (int i = 0; i < gradient.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient[i]),
                    format.format(gradientbis[i]));
        }

        double[] x2 = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        double value2 = LKJCorDistribution.logPdf(x);
        double value2bis = transformedDistribution.logPdf(x);

        assertEquals("value", format.format(value2), format.format(value2bis));

        double[] gradient2 = LKJCorDistribution.getGradientLogDensity(x);
        double[] gradient2bis = transformedDistribution.getGradientLogDensity(x);

        assertEquals("length", gradient2.length, gradient2bis.length);
        for (int i = 0; i < gradient2.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient2[i]),
                    format.format(gradient2bis[i]));
        }
    }
}
