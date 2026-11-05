/*
 * GammaQuantileExtremeShapeTest.java
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

import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.math.distributions.GammaDistribution;
import junit.framework.TestCase;
import org.apache.commons.math.special.Gamma;

/**
 * Plan 7 H3: a maximizer's trial point under the scaled logit of (0, 1000) can map the gamma shape to ~1E-300, and
 * {@code GammaDistribution.quantile} then looped forever inside {@code GammaFunction.incompleteGamma} (four office
 * runs spun for 15-25 h in one likelihood evaluation). Every shape now returns, NaN where there is no quantile, so the
 * likelihood is non-finite and the line search backs off; shapes up to 1E14 keep their quantiles.
 */
public class GammaQuantileExtremeShapeTest extends TestCase {

    public GammaQuantileExtremeShapeTest(String name) {
        super(name);
    }

    public void testExtremeShapesReturn() throws Exception {
        final double[] shapes = {Double.NaN, Double.POSITIVE_INFINITY, 1E-300, 1E-100, 1E-20, 1E16, 1E20, 1E50, 1E300};
        for (double shape : shapes) {
            final long start = System.nanoTime();
            for (int i = 0; i < 4; i++) {
                final double q = GammaDistribution.quantile((2.0 * i + 1.0) / 8.0, shape, 1.0 / shape);
                assertTrue("shape=" + shape + " i=" + i + ": " + q, Double.isNaN(q) || q >= 0);
            }
            assertTrue("shape=" + shape, (System.nanoTime() - start) / 1E9 < 30.0);
        }
    }

    // the discrete gamma's category rates at such a shape are NaN, not a hang
    public void testCategoryRatesAtATinyShape() throws Exception {
        final double[] rates = new double[4];
        final double[] proportions = new double[4];
        GammaSiteRateModel.setEqualRates(rates, proportions, 1E-300, 4, 0);
        for (double rate : rates) {
            assertTrue(Double.isNaN(rate));
        }
    }

    // ordinary and large shapes: a valid, increasing set of quantiles around the mean 1 (the iteration caps never bind)
    public void testOrdinaryShapesUnchanged() throws Exception {
        for (double shape : new double[]{1E-3, 0.05, 0.35, 1.0, 10.0, 1000.0, 1E6, 1E12, 1E14}) {
            double previous = -1.0;
            for (int i = 0; i < 4; i++) {
                final double q = GammaDistribution.quantile((2.0 * i + 1.0) / 8.0, shape, 1.0 / shape);
                assertFalse("shape=" + shape, Double.isNaN(q));
                assertTrue("shape=" + shape, q >= previous);
                previous = q;
            }
        }
        // the quantiles invert the CDF (commons-math's regularized gamma) at alpha = 0.5
        for (double u : new double[]{0.125, 0.375, 0.625, 0.875}) {
            final double q = GammaDistribution.quantile(u, 0.5, 2.0);
            assertEquals(u, Gamma.regularizedGammaP(0.5, q / 2.0), 1E-6);
        }
    }
}
