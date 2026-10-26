/*
 * SiteRateDifferentialTest.java
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
package test.dr.evomodel.siteratemodel;

import dr.evomodel.siteratemodel.DifferentiableSiteRateModel;
import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateDelegate;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;
import org.apache.commons.math.special.Gamma;

/**
 * Derivatives of the discretized site-rate categories (rates and proportions) wrt the gamma shape, the proportion
 * of invariant sites and the relative rate, for GammaSiteRateModel and DiscretizedSiteRateModel(GammaSiteRateDelegate).
 */
public class SiteRateDifferentialTest extends TestCase {

    private static final int K = 4;

    public SiteRateDifferentialTest(String name) {
        super(name);
    }

    public void testShapeAgainstPreciseRates() throws Exception {
        // the models' own quantiles carry ~1e-9 noise from AS32, so differentiate the same arithmetic with precise ones
        for (boolean invariant : new boolean[]{false, true}) {
            for (double alpha : new double[]{0.2, 0.5, 2.0, 20.0}) {
                final double mu = 1.7;
                for (DifferentiableSiteRateModel model : models(alpha, invariant ? 0.3 : Double.NaN, mu)) {
                    final Parameter shape = shapeOf(model);
                    final int n = model.getCategoryCount();
                    double[] dRates = new double[n];
                    double[] dProportions = new double[n];
                    model.getCategoryDifferentials(shape, dRates, dProportions);

                    final double h = 1E-5 * alpha;
                    double[] plus = preciseRates(alpha + h, invariant, mu);
                    double[] minus = preciseRates(alpha - h, invariant, mu);
                    for (int c = 0; c < n; c++) {
                        double numeric = (plus[c] - minus[c]) / (2 * h);
                        assertEquals("alpha=" + alpha + " c=" + c, numeric, dRates[c], 1E-6 * Math.max(1E-2, Math.abs(numeric)));
                        assertEquals(0.0, dProportions[c], 0.0);
                    }
                }
            }
        }
    }

    public void testInvariantAndRelativeRateAgainstFiniteDifferences() {
        for (double alpha : new double[]{0.5, Double.NaN}) {
            for (DifferentiableSiteRateModel model : models(alpha, 0.3, 1.7)) {
                checkFiniteDifference(model, invariantOf(model));
                checkFiniteDifference(model, nuOf(model));
            }
        }
    }

    public void testParameterDetection() {
        for (DifferentiableSiteRateModel model : models(0.5, 0.3, 1.7)) {
            assertTrue(model.isDifferentiableWrtParameter(shapeOf(model)));
            assertTrue(model.isDifferentiableWrtParameter(invariantOf(model)));
            assertTrue(model.isDifferentiableWrtParameter(nuOf(model)));
            assertFalse(model.isDifferentiableWrtParameter(new Parameter.Default(1.0)));
        }
    }

    // the invariant proportion and the relative rate enter without quantiles, so plain differences are exact enough
    private static void checkFiniteDifference(DifferentiableSiteRateModel model, Parameter parameter) {
        final int n = model.getCategoryCount();
        double[] dRates = new double[n];
        double[] dProportions = new double[n];
        model.getCategoryDifferentials(parameter, dRates, dProportions);

        final double value = parameter.getParameterValue(0);
        final double h = 1E-6;
        parameter.setParameterValue(0, value + h);
        double[] ratesPlus = model.getCategoryRates().clone();
        double[] proportionsPlus = model.getCategoryProportions().clone();
        parameter.setParameterValue(0, value - h);
        double[] ratesMinus = model.getCategoryRates().clone();
        double[] proportionsMinus = model.getCategoryProportions().clone();
        parameter.setParameterValue(0, value);

        for (int c = 0; c < n; c++) {
            assertEquals(parameter.getId() + " rate " + c, (ratesPlus[c] - ratesMinus[c]) / (2 * h), dRates[c], 1E-7);
            assertEquals(parameter.getId() + " proportion " + c,
                    (proportionsPlus[c] - proportionsMinus[c]) / (2 * h), dProportions[c], 1E-7);
        }
    }

    private static DifferentiableSiteRateModel[] models(double alpha, double pInv, double mu) {
        Parameter shape = Double.isNaN(alpha) ? null : new Parameter.Default("alpha", alpha);
        Parameter invariant = Double.isNaN(pInv) ? null : new Parameter.Default("pInv", pInv);
        Parameter nu = new Parameter.Default("nu", mu);
        GammaSiteRateModel direct = new GammaSiteRateModel("siteModel", nu, 1.0, shape,
                shape == null ? 0 : K, invariant);

        Parameter shape2 = Double.isNaN(alpha) ? null : new Parameter.Default("alpha", alpha);
        Parameter invariant2 = Double.isNaN(pInv) ? null : new Parameter.Default("pInv", pInv);
        Parameter nu2 = new Parameter.Default("nu", mu);
        DiscretizedSiteRateModel discretized = new DiscretizedSiteRateModel("siteModel", nu2, 1.0,
                new GammaSiteRateDelegate("delegate", shape2, K, invariant2));
        return new DifferentiableSiteRateModel[]{direct, discretized};
    }

    private static Parameter shapeOf(DifferentiableSiteRateModel model) {
        return find(model, "alpha");
    }

    private static Parameter invariantOf(DifferentiableSiteRateModel model) {
        return find(model, "pInv");
    }

    private static Parameter nuOf(DifferentiableSiteRateModel model) {
        return find(model, "nu");
    }

    private static Parameter find(DifferentiableSiteRateModel model, String id) {
        for (int i = 0; i < model.getVariableCount(); i++) {
            if (id.equals(model.getVariable(i).getId())) {
                return (Parameter) model.getVariable(i);
            }
        }
        for (int i = 0; i < model.getModelCount(); i++) {
            dr.inference.model.Model inner = model.getModel(i);
            for (int j = 0; j < inner.getVariableCount(); j++) {
                if (id.equals(inner.getVariable(j).getId())) {
                    return (Parameter) inner.getVariable(j);
                }
            }
        }
        throw new IllegalArgumentException(id);
    }

    // GammaSiteRateDelegate.getCategories followed by mu, with quantiles from commons-math at 1e-16
    private static double[] preciseRates(double alpha, boolean invariant, double mu) throws Exception {
        final int offset = invariant ? 1 : 0;
        final int n = K + offset;
        double[] rates = new double[n];
        double sum = 0.0;
        for (int i = 0; i < K; i++) {
            rates[i + offset] = preciseStandardQuantile((2.0 * i + 1.0) / (2.0 * K), alpha) / alpha;
            sum += rates[i + offset];
        }
        for (int c = 0; c < n; c++) {
            rates[c] *= mu * n / sum;
        }
        return rates;
    }

    private static double preciseStandardQuantile(double u, double a) throws Exception {
        double lower = -800.0;
        double upper = Math.log(a + 60.0 * Math.sqrt(a) + 60.0);
        for (int i = 0; i < 200; i++) {
            final double middle = 0.5 * (lower + upper);
            if (Gamma.regularizedGammaP(a, Math.exp(middle), 1E-16, 1000000) < u) {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        return Math.exp(0.5 * (lower + upper));
    }
}
