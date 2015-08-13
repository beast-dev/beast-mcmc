/*
 * GeneralizedIntegerGammaTest.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.math;

import dr.evomodel.branchratemodel.LatentStateBranchRateModel;
import dr.inference.markovjumps.MarkovReward;
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import dr.math.IntegrableUnivariateFunction;
import dr.math.distributions.GeneralizedIntegerGammaDistribution;
import dr.math.matrixAlgebra.Vector;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;

/**
 * @author Marc A. Suchard
 */

public class GeneralizedIntegerGammaTest extends MathTestCase {

    private static final double tolerance = 10E-6;

    private double sum(double[] v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x;
        }
        return sum;
    }

    int[] shape1s = new int[] { 10, 6, 2 };
    int[] shape2s = new int[] { 1, 3, 4 };
    double[] rate1s = new double[] { 0.1, 2, 3 };
    double[] rate2s = new double[] { 2, 0.1, 0.2 };
    double[] xs = new double[] { 0.1, 0.5, 1.2 };

    public void testGeneratingFunction() {

        for (int i = 0; i < shape1s.length; ++i) {
            GeneralizedIntegerGammaDistribution gig = new GeneralizedIntegerGammaDistribution(
                    shape1s[i], shape2s[i], rate1s[i], rate2s[i]
            );

            double a = gig.generatingFunction(xs[i]);
            double b = gig.generatingFunctionPartialFraction(xs[i]);
            assertEquals(a, b, tolerance);
        }
    }

    public void testPdf() {

        for (int i = 0; i < shape1s.length; ++i) {
            final GeneralizedIntegerGammaDistribution gig = new GeneralizedIntegerGammaDistribution(
                    shape1s[i], shape2s[i], rate1s[i], rate2s[i]
            );

            TrapezoidIntegrator integrator = new TrapezoidIntegrator();
            double m0 = 0.0;
            double m1 = 0.0;
            double m2 = 0.0;
            try {
                m0 = integrator.integrate(
                        new UnivariateRealFunction() {
                            @Override
                            public double value(double x) throws FunctionEvaluationException {
                                final double pdf = gig.pdf(x);
                                return pdf;
                            }
                        }, 0.0, 1000.0);

                m1 = integrator.integrate(
                        new UnivariateRealFunction() {
                            @Override
                            public double value(double x) throws FunctionEvaluationException {
                                final double pdf = gig.pdf(x);
                                return x * pdf;
                            }
                        }, 0.0, 1000.0);


                m2 = integrator.integrate(
                        new UnivariateRealFunction() {
                            @Override
                            public double value(double x) throws FunctionEvaluationException {
                                final double pdf = gig.pdf(x);
                                return x * x * pdf;
                            }
                        }, 0.0, 1000.0);

            } catch (MaxIterationsExceededException e) {
                e.printStackTrace();
            } catch (FunctionEvaluationException e) {
                e.printStackTrace();
            }

            // Check normalization
            assertEquals(1.0, m0, tolerance);

            // Check mean
            double mean = shape1s[i] / rate1s[i] + shape2s[i] / rate2s[i];
            assertEquals(mean, m1, tolerance);

            // Check variance
            m2 -= m1 * m1;
            double variance = shape1s[i] / rate1s[i] / rate1s[i] + shape2s[i] / rate2s[i] / rate2s[i];
            assertEquals(variance, m2, tolerance * 10); // Harder to approximate
        }
    }
}
