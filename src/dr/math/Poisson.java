/*
 * Poisson.java
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


/**
 * @author Alexei Drummond
 *
 * @version $Id: Poisson.java,v 1.4 2005/04/28 16:51:44 rambaut Exp $
 */
public class Poisson {

	private Poisson() {
	}

	// oldm is a flag for whether xm has been changed since last call.
    private static double oldm = -1.0;
    private static double sq, alxm, g;

    // used by gammaln
    private static final double[] cof = {   76.18009172947146, -86.50532032941677,
				                            24.01409824083091,  -1.231739572450155,
				                            0.1208650973866179e-2, -0.5395239384953e-5};

    /**
     * Returns an integer value that is a
     * random deviate drawn from a Poisson distribution of mean xm.
     */
    public static int nextPoisson(double xm) {

        double em; // expected mean
        double t, y;

        if (xm < 12.0) { // use direct method
			if (xm != oldm) {
			    oldm = xm;
			    g = Math.exp(-xm);  // if xm is new compute the exponential
			}
			em = -1.0;
			t = 1.0;
			do {
			    ++em;
			    t *= MathUtils.nextDouble();
			} while (t > g);
		} else {
			if (xm != oldm) {
			    oldm = xm;
			    sq = Math.sqrt(2.0 * xm);
			    alxm = Math.log(xm);
			    g = xm * alxm - gammln(xm + 1.0);
			    // The function gammln is the natural log of the gamma function
			}
			do {
			    do {
			        y = Math.tan(Math.PI * MathUtils.nextDouble());
			        em = sq * y + xm; // em is y shifted and scaled
			    } while (em < 0.0); // reject if in realm of zero probability
			    em = Math.floor(em);
			    t = 0.9 *(1.0 + y*y) * Math.exp(em*alxm-gammln(em + 1.0)-g);
			    // The ratio of the desired distribution to the comparison function;
			    // we accept or reject by comparing it to another uniform deviate.
			    // The factor 0.9 is chosen so that t never exceeds 1.
			} while (MathUtils.nextDouble() > t);
        }
        return (int)em;
    }

    /**
     * Returns the value of log gamma(xx) for xx > 0.
     */
    public static double gammln(double xx) {
		double x, y, tmp, ser;

		int j;

		y = x = xx;
		tmp = x + 5.5;
		tmp -= (x + 0.5) * Math.log(tmp);
		ser = 1.000000000190015;
		for (j = 0; j <= 5; j++) {
			ser += cof[j] / ++y;
		}
		return -tmp + Math.log(2.5066282746310005 * ser / x);
	}
}
