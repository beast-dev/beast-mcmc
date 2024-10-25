/*
 * LambertW.java
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

/* Lambert W function, code by various authors */
/* Ported to Java by Daniel Wilson			*/

package dr.evomodel.epidemiology;

/* Author:  G. Jungman */

/* Started with code donated by K. Briggs; added
 * error estimates, GSL foo, and minor tweaks.
 * Some Lambert-ology from
 *  [Corless, Gonnet, Hare, and Jeffrey, "On Lambert's W Function".]
 */

public class LambertW { 
	
	static public class gsl_sf_result {
		public double val, err;
		public gsl_sf_result() {
			val = err = 0;
		};
	}
	
	public enum GSL_RETURN {GSL_SUCCESS, GSL_EMAXITER, GSL_EDOM} ;

	static double M_E = 2.71828182845904523536028747135266250;
	static double GSL_DBL_EPSILON = 2.2204460492503131e-16;
	
	/* Halley iteration (eqn. 5.12, Corless et al) */
	public static GSL_RETURN
	halley_iteration(
			double x,
			double w_initial,
			int max_iters,
			gsl_sf_result result
	)
	{
		double w = w_initial;
		int i;

		for(i=0; i<max_iters; i++) {
			double tol;
			final double e = Math.exp(w);
			final double p = w + 1.0;
			double t = w*e - x;
			/* printf("FOO: %20.16g  %20.16g\n", w, t); */

			if (w > 0) {
				t = (t/p)/e;  /* Newton iteration */
			} else {
				t /= e*p - 0.5*(p + 1.0)*t/p;  /* Halley iteration */
			};

			w -= t;

			tol = 10 * GSL_DBL_EPSILON * Math.max(Math.abs(w), 1.0/(Math.abs(p)*e));

			if(Math.abs(t) < tol)
			{
				result.val = w;
				result.err = 2.0*tol;
				return GSL_RETURN.GSL_SUCCESS;
			}
		}

		/* should never get here */
		result.val = w;
		result.err = Math.abs(w);
		return GSL_RETURN.GSL_EMAXITER;
	}


	static final double[] c = {
			-1.0,
			2.331643981597124203363536062168,
			-1.812187885639363490240191647568,
			1.936631114492359755363277457668,
			-2.353551201881614516821543561516,
			3.066858901050631912893148922704,
			-4.175335600258177138854984177460,
			5.858023729874774148815053846119,
			-8.401032217523977370984161688514,
			12.250753501314460424,
			-18.100697012472442755,
			27.029044799010561650
		};
	/* series which appears for q near zero;
	 * only the argument is different for the different branches
	 */
	public static double
	series_eval(double r)
	{
		final double t_8 = c[8] + r*(c[9] + r*(c[10] + r*c[11]));
		final double t_5 = c[5] + r*(c[6] + r*(c[7]  + r*t_8));
		final double t_1 = c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
		return c[0] + r*t_1;
	}


	/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
	static GSL_RETURN
	gsl_sf_lambert_W0_e(double x, gsl_sf_result result)
	{
		final double one_over_E = 1.0/M_E;
		final double q = x + one_over_E;

		if(x == 0.0) {
			result.val = 0.0;
			result.err = 0.0;
			return GSL_RETURN.GSL_SUCCESS;
		}
		else if(q < 0.0) {
			/* Strictly speaking this is an error. But because of the
			 * arithmetic operation connecting x and q, I am a little
			 * lenient in case of some epsilon overshoot. The following
			 * answer is quite accurate in that case. Anyway, we have
			 * to return GSL_EDOM.
			 */
			result.val = -1.0;
			result.err =  Math.sqrt(-q);
			return GSL_RETURN.GSL_EDOM;
		}
		else if(q == 0.0) {
			result.val = -1.0;
			result.err =  GSL_DBL_EPSILON; /* cannot error is zero, maybe q == 0 by "accident" */
			return GSL_RETURN.GSL_SUCCESS;
		}
		else if(q < 1.0e-03) {
			/* series near -1/E in sqrt(q) */
			final double r = Math.sqrt(q);
			result.val = series_eval(r);
			result.err = 2.0 * GSL_DBL_EPSILON * Math.abs(result.val);
			return GSL_RETURN.GSL_SUCCESS;
		}
		else {
			final int MAX_ITERS = 10;
			double w;

			if (x < 1.0) {
				/* obtain initial approximation from series near x=0;
				 * no need for extra care, since the Halley iteration
				 * converges nicely on this branch
				 */
				final double p = Math.sqrt(2.0 * M_E * q);
				w = -1.0 + p*(1.0 + p*(-1.0/3.0 + p*11.0/72.0)); 
			}
			else {
				/* obtain initial approximation from rough asymptotic */
				w = Math.log(x);
				if(x > 3.0) w -= Math.log(w);
			}

			return halley_iteration(x, w, MAX_ITERS, result);
		}
	}


	static GSL_RETURN
	gsl_sf_lambert_Wm1_e(double x, gsl_sf_result result)
	{
		if(x > 0.0) {
			return gsl_sf_lambert_W0_e(x, result);
		}
		else if(x == 0.0) {
			result.val = 0.0;
			result.err = 0.0;
			return GSL_RETURN.GSL_SUCCESS;
		}
		else {
			final int MAX_ITERS = 32;
			final double one_over_E = 1.0/M_E;
			final double q = x + one_over_E;
			double w;

			if (q < 0.0) {
				/* As in the W0 branch above, return some reasonable answer anyway. */
				result.val = -1.0; 
				result.err =  Math.sqrt(-q);
				return GSL_RETURN.GSL_EDOM;
			}

			if(x < -1.0e-6) {
				/* Obtain initial approximation from series about q = 0,
				 * as long as we're not very close to x = 0.
				 * Use full series and try to bail out if q is too small,
				 * since the Halley iteration has bad convergence properties
				 * in finite arithmetic for q very small, because the
				 * increment alternates and p is near zero.
				 */
				final double r = -Math.sqrt(q);
				w = series_eval(r);
				if(q < 3.0e-3) {
					/* this approximation is good enough */
					result.val = w;
					result.err = 5.0 * GSL_DBL_EPSILON * Math.abs(w);
					return GSL_RETURN.GSL_SUCCESS;
				}
			}
			else {
				/* Obtain initial approximation from asymptotic near zero. */
				final double L_1 = Math.log(-x);
				final double L_2 = Math.log(-L_1);
				w = L_1 - L_2 + L_2/L_1;
			}

			return halley_iteration(x, w, MAX_ITERS, result);
		}
	}


	/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

	static public double branch0(double x)
	{
		gsl_sf_result result = new gsl_sf_result();
		GSL_RETURN res = gsl_sf_lambert_W0_e(x, result); 
		if(res==GSL_RETURN.GSL_EMAXITER) {
			throw new RuntimeException("Too many iterations");
		}
		return result.val;
	}

	static public double branchNeg1(double x)
	{
		gsl_sf_result result = new gsl_sf_result();
		GSL_RETURN res = gsl_sf_lambert_Wm1_e(x, result); 
		if(res==GSL_RETURN.GSL_EMAXITER) {
			throw new RuntimeException("Too many iterations");
		}
		return result.val;
	}

}
