/*
 * ModifiedBesselFirstKindTest.java
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

import dr.math.MathUtils;
import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

import static dr.math.ModifiedBesselFirstKind.*;

public class ModifiedBesselFirstKindTest extends TestCase {

    private static final double TOL = 1e-2;


    public void testModifiedBesselReal() {

        double[] xs = new double[]{1e-2, 1e-1, 1, 1e1, 1e2};
        int[] orders = new int[]{1, 2, 5, 10};

        double[] ys = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            ys[i] = 1 / xs[i];
        }

        for (int i = 0; i < xs.length; i++) {
            for (int j = 0; j < orders.length; j++) {

                double bInt = bessi(xs[i], orders[j]);
                double bReal = bessi(xs[i], (double) orders[j]);

                if (!MathUtils.isRelativelyClose(bInt, bReal, TOL)) {
                    throw new AssertionFailedError("failed");
                }

                double byReal = bessi(ys[i], (double) orders[j]);
                double ratio = bessIRatio(xs[i], ys[i], orders[j]);

                if (!MathUtils.isRelativelyClose(bReal / byReal, ratio, TOL)) {
                    throw new AssertionFailedError("failed");
                }
            }
        }


    }


}
