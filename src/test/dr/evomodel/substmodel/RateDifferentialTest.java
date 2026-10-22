/*
 * RateDifferentialTest.java
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

package test.dr.evomodel.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.WrtRateParameter;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import junit.framework.TestCase;

/**
 * Checks the generic derivative of the normalized generator Q~ = Q / beta with respect to relative rates
 * (TN93 kappa1/kappa2, GTR exchangeabilities in both parameterizations).
 */
public class RateDifferentialTest extends TestCase {

    private static final int K = 4;
    private static final double[] PI = {0.1, 0.2, 0.3, 0.4};
    private static final double[] RATES = {1.0, 2.0, 0.5, 1.5, 3.0, 1.0};

    public RateDifferentialTest(String name) {
        super(name);
    }

    public void testTn93MatchesFiniteDifferenceOfGenerator() {
        Parameter kappa1 = new Parameter.Default(2.0);
        Parameter kappa2 = new Parameter.Default(4.0);
        TN93 tn93 = new TN93(kappa1, kappa2, frequencyModel());
        checkAgainstFiniteDifference(tn93, kappa1);
        checkAgainstFiniteDifference(tn93, kappa2);
        assertEquals(1, ((WrtRateParameter) tn93.factory(kappa1, 0)).getRateIndex());
        assertEquals(4, ((WrtRateParameter) tn93.factory(kappa2, 0)).getRateIndex());
    }

    public void testGtrMatchesFiniteDifferenceOfGenerator() {
        Parameter rates = new Parameter.Default(RATES.clone());
        checkAgainstFiniteDifference(new GTR(rates, frequencyModel()), rates);

        Parameter[] variables = new Parameter[RATES.length];
        for (int m = 0; m < RATES.length; m++) {
            variables[m] = new Parameter.Default(RATES[m]);
        }
        GTR gtr = new GTR(variables[0], variables[1], variables[2], variables[3], variables[4], variables[5],
                frequencyModel());
        for (int m = 0; m < RATES.length; m++) {
            checkAgainstFiniteDifference(gtr, variables[m]);
            assertEquals(m, ((WrtRateParameter) gtr.factory(variables[m], 0)).getRateIndex());
        }
    }

    public void testGtrWithHkyRatesMatchesHkyKappa() {
        // kappa scales both transitions, so dQ~/dkappa = dQ~/ds_AG + dQ~/ds_CT; HKY uses its own enum differential
        final double kappa = 2.7;
        Parameter rates = new Parameter.Default(new double[]{1.0, kappa, 1.0, 1.0, kappa, 1.0});
        GTR gtr = new GTR(rates, frequencyModel());
        Parameter kappaParameter = new Parameter.Default(kappa);
        HKY hky = new HKY(kappaParameter, frequencyModel());

        double[][] ag = toArray(gtr.getInfinitesimalDifferentialMatrix(gtr.factory(rates, 1)));
        double[][] ct = toArray(gtr.getInfinitesimalDifferentialMatrix(gtr.factory(rates, 4)));
        double[][] viaHky = toArray(hky.getInfinitesimalDifferentialMatrix(hky.factory(kappaParameter, 0)));
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                assertEquals("(" + i + "," + j + ")", viaHky[i][j], ag[i][j] + ct[i][j], 1E-14);
            }
        }
    }

    public void testRateWeightedSumVanishes() {
        // Q is homogeneous of degree 1 and beta of degree 1 in the rates, so sum_m s_m dQ~/ds_m = 0
        Parameter rates = new Parameter.Default(RATES.clone());
        GTR gtr = new GTR(rates, frequencyModel());
        double[][] sum = new double[K][K];
        for (int m = 0; m < RATES.length; m++) {
            double[][] d = toArray(gtr.getInfinitesimalDifferentialMatrix(gtr.factory(rates, m)));
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    sum[i][j] += RATES[m] * d[i][j];
                }
            }
        }
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                assertEquals(0.0, sum[i][j], 1E-13);
            }
        }
    }

    public void testRateFactoryRange() {
        Parameter rates = new Parameter.Default(RATES.clone());
        GTR gtr = new GTR(rates, frequencyModel());
        try {
            gtr.factory(rates, RATES.length);
            fail("dimension " + RATES.length + " must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    // dQ~/dtheta against central differences of the model's own normalized generator
    private static void checkAgainstFiniteDifference(BaseSubstitutionModel model, Parameter parameter) {
        final double h = 1E-6;
        for (int d = 0; d < parameter.getDimension(); d++) {
            WrappedMatrix analytic = ((DifferentiableSubstitutionModel) model)
                    .getInfinitesimalDifferentialMatrix(((DifferentiableSubstitutionModel) model).factory(parameter, d));

            final double value = parameter.getParameterValue(d);
            double[] plus = new double[K * K];
            double[] minus = new double[K * K];
            parameter.setParameterValue(d, value + h);
            model.getInfinitesimalMatrix(plus);
            parameter.setParameterValue(d, value - h);
            model.getInfinitesimalMatrix(minus);
            parameter.setParameterValue(d, value);

            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    double numeric = (plus[i * K + j] - minus[i * K + j]) / (2 * h);
                    assertEquals(model.getModelName() + " d=" + d + " (" + i + "," + j + ")",
                            numeric, analytic.get(i, j), 1E-7 * Math.max(1.0, Math.abs(numeric)));
                }
            }
        }
    }

    private static FrequencyModel frequencyModel() {
        return new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(PI.clone()));
    }

    private static double[][] toArray(WrappedMatrix matrix) {
        double[][] out = new double[K][K];
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                out[i][j] = matrix.get(i, j);
            }
        }
        return out;
    }
}
