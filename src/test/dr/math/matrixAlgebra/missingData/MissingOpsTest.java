/*
 * MissingOpsTest.java
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

package test.dr.math.matrixAlgebra.missingData;

import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Paul Bastide
 */

public class MissingOpsTest {

    interface Instance {

        int getRank();

        double getLogDet();

        DenseMatrix64F getMatrix();

        abstract class Basic implements Instance {

        }
    }

    Instance test0 = new Instance.Basic() {
        public int getRank() {
            return 2;
        }

        public double getLogDet() {
            return -42.160726233513714;
        }

        public DenseMatrix64F getMatrix() {
            return new DenseMatrix64F(2, 2, true,
                    0.004255873897787696, 0.010962329615067505,
                    0.010962329615067505, 0.02823689645782389);
        }
    };

    Instance test1 = new Instance.Basic() {
        public int getRank() {
            return 1;
        }

        public double getLogDet() {
            return -3.4267272543138203;
        }

        public DenseMatrix64F getMatrix() {
            return new DenseMatrix64F(2, 2, true,
                    0.004255918204391931, 0.010962443740752292,
                    0.010962443740752292, 0.028237190424652364);
        }
    };


    Instance[] all = {test0, test1};

    @Test
    public void safeDeterminant() throws Exception {
        for (Instance test : all) {

            int rank_test = test.getRank();
            double logDet_test = test.getLogDet();

            DenseMatrix64F P = test.getMatrix();

            InversionResult c = MissingOps.safeDeterminant(P, false);

            assertEquals(c.getLogDeterminant(), logDet_test, 1e-6);
            assertEquals(c.getEffectiveDimension(), rank_test, 1e-6);
        }
    }
}