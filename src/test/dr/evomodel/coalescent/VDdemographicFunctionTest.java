/*
 * VDdemographicFunctionTest.java
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

package test.dr.evomodel.coalescent;

import junit.framework.TestCase;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evolution.util.Units;

/**
 * @author Joseph Heled
 *         Date: 23/06/2009
 */
public class VDdemographicFunctionTest extends TestCase {
    public void testExp() {
        // test that numerical and exact integration match (up to a point, numerical is not that good for those 
        // exponential gone to constant transitions.
        {
            double[] times = {1, 3};
            double[] logPops = {0, 2, 3};
            VariableDemographicModel.Type type = VariableDemographicModel.Type.EXPONENTIAL;
            Units.Type units = Units.Type.SUBSTITUTIONS;
            final VDdemographicFunction f = new VDdemographicFunction(times, logPops, units, type);

            double[][] vals = {{0, 2, 1e-9}, {1, 2, 1e-9}, {0, 5, 1e-6}, {2, 6, 1e-6}};
            for(double[] c : vals) {
                double v1 = f.getIntegral(c[0], c[1]);
                double v2 = f.getNumericalIntegral(c[0], c[1]);
                assertEquals(Math.abs(1 - v1 / v2), 0, c[2]);
            }
        }

        {
            double[] times = {1, 3};
            // try a const interval
            double[] logPops = {0, 0, 3};
            VariableDemographicModel.Type type = VariableDemographicModel.Type.EXPONENTIAL;
            Units.Type units = Units.Type.SUBSTITUTIONS;
            final VDdemographicFunction f = new VDdemographicFunction(times, logPops, units, type);

            double[][] vals = {{0, .7, 1e-9}, {1, 2, 1e-9}, {0, 5, 1e-6}, {2, 6, 1e-6}};
            for(double[] c : vals) {
                double v1 = f.getIntegral(c[0], c[1]);
                double v2 = f.getNumericalIntegral(c[0], c[1]);
                assertEquals(Math.abs(1 - v1 / v2), 0, c[2]);
            }
        }
    }
}
