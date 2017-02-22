/*
 * NegativeBinomialDistribution.java
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

package dr.math.distributions;
import dr.math.ErrorFunction;

import org.apache.commons.math.MathException;
import org.apache.commons.math.special.Beta;

import dr.math.*;

/**
 * @author Trevor Bedford
 * @version $Id$
 */
public class NegativeBinomialDistribution implements Distribution {

    double mean;
    double alpha;

    public NegativeBinomialDistribution(double mean, double alpha) {
        this.mean = mean;
        this.alpha = alpha;
    }

    public double pdf(double x) {
        return pdf(x, mean, alpha);
    }

    public double logPdf(double x) {
        return logPdf(x, mean, alpha);
    }

    public double cdf(double x) {
        return cdf(x, mean, alpha);
    }

    public double quantile(double y) {
    	// TB - I'm having trouble implementing this
    	// LM - A first stab using simple minimisation to invert the function (under absolute loss)
    	// Implementation based on the qnbinom.c function used in R
    	 final double stdev = Math.sqrt(mean + (mean * mean * alpha));
    	 final double r = -1 * (mean*mean) / (mean - stdev*stdev);
         final double p = mean / (stdev*stdev);
         final double prob = y;
         
         final double Q = 1.0 / p;
         final double P = (1.0 - p) * Q;
         final double gamma = (Q + P)/stdev;
         final double z = Math.sqrt(2.0) * ErrorFunction.inverseErf(2.0 * y - 1.0);
         final double crudeY = mean + stdev * (z + gamma * (z*z - 1) / 6);
         
    	UnivariateFunction f = new UnivariateFunction() {
    		double tent = Double.NaN;
    		public double evaluate(final double argument) {
    			try {
    	            tent = Beta.regularizedBeta(p, r, argument+1);
    	        } catch (MathException e) {
    	            return Double.NaN;
    	        }
    			double score = Math.abs(tent-prob);
    			return score;
    		}
    		public int getNumArguments() {
    			return 1;
    		}
    		public double getLowerBound() { // 20% window should cut it. Probably too large even...
    			return Math.min(crudeY - .2*crudeY, 0);
    		}
    		public double getUpperBound() {
    			return crudeY + .2*crudeY;
    		}
    	};
    	UnivariateMinimum minimum = new UnivariateMinimum();
    	double q = minimum.findMinimum(f);
     	return Math.ceil(q);
    }

    public double mean() {
        return mean;
    }

    public double variance() {
        return mean + (mean * mean * alpha);
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException();
    }


    public static double pdf(double x, double mean, double alpha) {
        if (x < 0)  return 0;
        return Math.exp(logPdf(x, mean, alpha));
    }

    public static double logPdf(double x, double mean, double alpha) {
        if (x < 0)  return Double.NEGATIVE_INFINITY;
//        double r = -1 * (mean*mean) / (mean - stdev*stdev);
//        double p = mean / (stdev*stdev);
//        return Math.log(Math.pow(1-p,x)) + Math.log(Math.pow(p, r)) + GammaFunction.lnGamma(r+x) - GammaFunction.lnGamma(r) - GammaFunction.lnGamma(x+1);
        double theta = 1.0 / alpha;

        double p = theta / (theta + mean);
        return Math.log(1 - p) * x + Math.log(p) * theta + GammaFunction.lnGamma(theta + x) - GammaFunction.lnGamma(theta) - GammaFunction.lnGamma(x+1);
    }

    public static double cdf(double x, double mean, double alpha) {
        double theta = 1.0 / alpha;
        double p = theta / (theta + mean);
        try {
            return Beta.regularizedBeta(p, theta, x+1);
        } catch (MathException e) {
            // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
            // exceptions is bad. Better to return NaN and let the calling code deal with it.
            return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                "Couldn't calculate beta cdf for alpha = " + alpha + ", beta = " + beta + ": " +e.getMessage());
        }
    }


    public static void main(String[] args) {
        System.out.println("Test negative binomial");
        double mean = 5;
        double stdev = 5;
//         double r = -1 * (mean*mean) / (mean - stdev*stdev);
        double alpha = (stdev * stdev - mean) / (mean * mean);

        NegativeBinomialDistribution dist = new NegativeBinomialDistribution(5, alpha);
        System.out.println("Mean 5, sd 5, x 5, pdf 0.074487, logPdf -2.59713, median 4");
        System.out.println("pdf = " + dist.pdf(5));
        System.out.println("quantile(0.5) aka median = " + dist.quantile(0.5));
        System.out.println("logPdf = " + dist.logPdf(5));
    }

}
