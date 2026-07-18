/*
 * GumbelDistribution.java
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

package dr.math.distributions;

import dr.math.MathUtils;
import dr.inference.model.GradientProvider;
import dr.math.UnivariateFunction;

/**
 * Type 1 (location-scale) Gumbel distribution, minimum orientation.
 * https://en.wikipedia.org/wiki/Gumbel_distribution
 * <p/>
 * This is the "reversed" (minimum) orientation of the standard Type 1 Gumbel:
 * with z = (x - location) / scale,
 * <pre>
 *   pdf(x) = exp(z - exp(z)) / scale
 *   cdf(x) = 1 - exp(-exp(z))
 * </pre>
 * (the textbook/maximum orientation instead has pdf(x) = exp(-z - exp(-z)) / scale).
 * <p/>
 * This orientation is exactly the distribution of log(R) when R is Exponential
 * distributed: if R ~ Exponential(mean = theta), then log(R) ~ GumbelDistribution
 * (location = log(theta), scale = 1). Its right tail (large x) decays
 * double-exponentially, matching the fast decay of the Exponential's right tail in R;
 * its left tail (x to -infinity) decays only single-exponentially, matching the
 * Exponential density remaining flat/nonzero as R approaches 0.
 * <p/>
 * (Parameters: location (mu), scale (beta > 0);
 *  mean: location - beta * gamma (gamma = Euler-Mascheroni constant);
 *  variance: beta^2 * pi^2 / 6)
 */
public class GumbelDistribution implements Distribution, GradientProvider {

    // Euler-Mascheroni constant
    public static final double EULER_MASCHERONI = 0.5772156649015328606;

    /**
     * Constructor
     *
     * @param location, scale the parameters of the Gumbel type I distribution
     */
    public GumbelDistribution(double location, double scale) {
        this.location = location;
        this.scale = scale;
    }

    public double getLocation() {
        return location;
    }

    public void setLocation(double value) {
        location = value;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double value) {
        scale = value;
    }

    public double pdf(double x) {
        return pdf(x, location, scale);
    }

    public double logPdf(double x) {
        return logPdf(x, location, scale);
    }

    public double cdf(double x) {
        return cdf(x, location, scale);
    }

    public double quantile(double y) {
        return quantile(y, location, scale);
    }

    public double mean() {
        return mean(location, scale);
    }

    public double variance() {
        return variance(location, scale);
    }

    public double nextGumbel() {
        return nextGumbel(location, scale);
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * probability density function of the (minimum-oriented) Gumbel type I distribution
     *
     * @param x                 argument
     * @param location, scale   parameters of Gumbel type I distribution
     * @return pdf value
     */
    public static double pdf(double x, double location, double scale) {
        return Math.exp(logPdf(x, location, scale));
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x                 argument
     * @param location, scale   parameters of Gumbel type I distribution
     * @return log pdf value
     */
    public static double logPdf(double x, double location, double scale) {
        final double z = (x - location) / scale;
        return z - Math.exp(z) - Math.log(scale);
    }

    public static double gradLogPdf(double x, double location, double scale) {
        final double z = (x - location) / scale;
        return (1.0 - Math.exp(z)) / scale;
    }

    public static double hessianLogPdf(double x, double location, double scale) {
        final double z = (x - location) / scale;
        return -Math.exp(z) / (scale * scale);
    }

    /**
     * cumulative density function of the (minimum-oriented) Gumbel type I distribution
     *
     * @param x                 argument
     * @param location, scale   parameters of distribution
     * @return cdf value
     */
    public static double cdf(double x, double location, double scale) {
        final double z = (x - location) / scale;
        return 1.0 - Math.exp(-Math.exp(z));
    }

    /**
     * quantile (inverse cumulative density function) of the (minimum-oriented) Gumbel type I distribution
     *
     * @param p                 argument
     * @param location, scale   parameters of Gumbel type I distribution
     * @return icdf value
     */
    public static double quantile(double p, double location, double scale) {
        return location + scale * Math.log(-Math.log(1.0 - p));
    }

    /**
     * mean of the (minimum-oriented) Gumbel type I distribution
     *
     * @param location, scale   parameters of Gumbel type I distribution
     * @return mean
     */
    public static double mean(double location, double scale) {
        return location - scale * EULER_MASCHERONI;
    }

    /**
     * variance of the Gumbel type I distribution
     *
     * @param location, scale   parameters of Gumbel type I distribution
     * @return variance
     */
    public static double variance(double location, double scale) {
        return scale * scale * Math.PI * Math.PI / 6.0;
    }

    public static double nextGumbel(double location, double scale) { // Inverse CDF generating method
        final double p = MathUtils.nextDouble();
        return quantile(p, location, scale);
    }

    public static void main(String[] args) {
        System.out.println("Test Gumbel type I (minimum orientation)");
        double location = -1.9635100260214235;
        double scale = 1.0;
        GumbelDistribution dist = new GumbelDistribution(location, scale);
        System.out.println("location = " + location + ", scale = " + scale);
        System.out.println("mean = " + dist.mean() + " (location - scale*gamma)");
        System.out.println("variance = " + dist.variance() + " (scale^2 * pi^2/6)");
        System.out.println("logPdf(location) = " + dist.logPdf(location));
        System.out.println("quantile(0.5) aka median = " + dist.quantile(0.5));
    }

    // parameters of the Gumbel type I distribution
    double location;
    double scale;

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = gradLogPdf(x[i], location, scale);
        }
        return result;
    }
}
