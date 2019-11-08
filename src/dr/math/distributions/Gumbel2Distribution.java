/*
 * Gumbel2Distribution.java
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
import dr.math.MathUtils;
import dr.math.GammaFunction;
import dr.math.UnivariateFunction;

/**
 * Type 2 Gumbel distribution.
 * https://en.wikipedia.org/wiki/Type-2_Gumbel_distribution
 * <p/>
 * (Parameters : shape (a),  scale (b);
 *  mean: pow(scale, 1/shape) * Gamma(1-1/shape);
 * 	variance: pow(scale, 2/shape)* Gamma(1-1/shape) - pow(Gamma(1-1/shape), 2) ; 
 * <p/>
 * @author Luiz Max Carvalho
 */
public class Gumbel2Distribution implements Distribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     *
     * @param a, b the parameters of the Gumbel type II distribution
     */
    public Gumbel2Distribution(double shape, double scale) {
        this.shape = shape;
        this.scale = scale;
    }

    public double getShape() {
        return shape;
    }

    public void setShape(double value) {
    	shape = value;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double value) {
    	scale = value;
    }
    
    public double pdf(double x) {
        return pdf(x, shape, scale);
    }

    public double logPdf(double x) {
        return logPdf(x, shape, scale);
    }

    public double cdf(double x) {
        return cdf(x, shape, scale);
    }

    public double quantile(double y) {
        return quantile(y, shape, scale);
    }

    public double mean() {
        return mean(shape, scale);
    }

    public double variance() {
        return variance(shape, scale);
    }

    public double nextGumbel() {
        return nextGumbel(shape, scale);
    }
    
    public UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * probability density function of the Gumbel II distribution
     *
     * @param x      argument
     * @param shape, scale parameters of Gumbel II distribution
     * @return pdf value
     */
    public static double pdf(double x, double shape, double scale) {
    	if (x < 0) return 0;
        return Math.exp(logPdf(x, shape, scale));
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @param shape, scale parameters of Gumbel II distribution
     * @return log pdf value
     */
    public static double logPdf(double x, double shape, double scale) {
    	if (x < 0) return Double.NEGATIVE_INFINITY;
    	double ans =  Math.log(shape) + Math.log(scale) -( shape + 1)*x -scale*Math.pow(x, -shape);
        return ans;
    }

    /**
     * cumulative density function of the Gumbel II distribution
     *
     * @param x      argument
     * @param shape,scale parameters of distribution
     * @return cdf value
     */
    public static double cdf(double x, double shape, double scale) {
    	if(x < 0) return(0);
        return Math.exp(-scale*Math.pow(x, -shape));
    }


    /**
     * quantile (inverse cumulative density function) of the Gumbel type II distribution
     *
     * @param p  argument
     * @param shape, scale parameters  of Gumbel type II distribution
     * @return icdf value
     */
    public static double quantile(double p, double shape, double scale) {
        return Math.pow(-scale/Math.log(p), 1/shape);
    }

    /**
     * mean of the Gumbel type II distribution
     *
     * @param shape, scale parameters of Gumbel type  II distribution
     * @return mean
     */
    public static double mean(double shape, double scale) {
        return Math.pow(scale, 1/shape)* Math.exp(GammaFunction.lnGamma(1- 1/shape));
    }

    /**
     * variance of the Gumbel type II distribution
     *
     * @param a, scale parameters of  Gumbel type II  distribution
     * @return variance
     */
    public static double variance(double shape, double scale) {
    	final double term = Math.exp(GammaFunction.lnGamma(1-1/shape));
        return Math.pow(scale, 2/shape)*(term - Math.pow(term, 2));
    }
    
    public static double nextGumbel(double shape, double scale) { // Inverse CDF generating method
    	final double p = MathUtils.nextDouble();
    	return quantile(p, shape, scale);
    }
    
    public static void main(String[] args) {
        System.out.println("Test Gumbel type II");
        double shape = 0.5;
        double scale = 2.302585;
        Gumbel2Distribution dist = new Gumbel2Distribution(shape, scale);
        System.out.println("x = 2, pdf 0.01125111, logPdf -4.487288, median 11.03521");
        System.out.println("pdf = " + dist.pdf(2));
        System.out.println("logPdf = " + dist.logPdf(2));
        System.out.println("quantile(0.5) aka median = " + dist.quantile(0.5));
    }
    
    // parameters of the Gumbel type II
    double shape;
    double scale;
}
