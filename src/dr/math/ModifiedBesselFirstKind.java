/*
 * ModifiedBesselFirstKind.java
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

package dr.math;


import cern.jet.math.Bessel;


public class ModifiedBesselFirstKind {

    //Adapted from Numerical Recipes for C

    public static final int ACC = 40;
    public static final double BIGNO = 1.0e10;
    public static final double BIGNI = 1.0e-10;

    /**
     * @param x argument
     * @param n order
     * @return the modified Bessel function of the first kind and nth order
     */
    public static double bessi(double x, int n) {

        if (n == 0) return Bessel.i0(x);
        if (n == 1) return Bessel.i1(x);

        int j;
        double bi, bim, bip, tox, ans;

        if (x == 0.0)
            return 0.0;
        else {
            tox = 2.0 / Math.abs(x);
            bip = ans = 0.0;
            bi = 1.0;
            for (j = 2 * (n + (int)Math.sqrt(ACC * n)); j > 0; j--) {
                bim = bip + j * tox * bi;
                bip = bi;
                bi = bim;
                if (Math.abs(bi) > BIGNO) {
                    ans *= BIGNI;
                    bi *= BIGNI;
                    bip *= BIGNI;
                }
                if (j == n) ans = bip;
            }
            ans *= Bessel.i0(x) / bi;
            return (x < 0.0 && ((n & 1) != 0)) ? -ans : ans;
        }
    }
}
