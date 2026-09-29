/*
 * CauchyDistributionTest.java
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

package test.dr.inference.distribution;

import dr.inference.distribution.CauchyDistributionModel;
import dr.inference.model.Parameter;
import dr.math.NumericalDerivative;
import dr.math.UnivariateFunction;
import junit.framework.TestCase;
import org.apache.commons.math.distribution.CauchyDistributionImpl;

/**
 * @author Marc A Suchard
 */

public class CauchyDistributionTest extends TestCase {

    public CauchyDistributionTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
    }

    public void testCauchyDistribution() {
        CauchyDistributionModel d1 = new CauchyDistributionModel(
                new Parameter.Default(2.0),
                new Parameter.Default(3.0)
        );

        CauchyDistributionImpl d2 = new CauchyDistributionImpl(2, 3);

        assertEquals(d1.pdf(-1), d2.density(-1), 1E-10);
        assertEquals(d1.logPdf(-1), Math.log(d2.density(-1)), 1E-10);

        double dx = NumericalDerivative.firstDerivative(new UnivariateFunction() {
            @Override
            public double evaluate(double argument) {
                return d1.logPdf(argument);
            }

            @Override
            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        }, 3.0);

        assertEquals(d1.getGradientLogDensity(new double[] { 3.0 } )[0], dx, 1E-4);
    }
}