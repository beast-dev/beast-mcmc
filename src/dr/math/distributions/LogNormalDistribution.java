/*
 * LogNormalDistribution.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.math.UnivariateFunction;
import dr.math.interfaces.OneVariableFunction;
import dr.math.iterations.BisectionZeroFinder;

/**
 * normal distribution (pdf, cdf, quantile)
 *
 * @author Korbinian Strimmer
 * @version $Id: LogNormalDistribution.java,v 1.3 2005/06/21 16:25:15 beth Exp $
 */
public class LogNormalDistribution implements Distribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     */
    public LogNormalDistribution(double M, double S) {
        this.M = M;
        this.S = S;
    }

    public final double getM() {
        return M;
    }

    public final void setM(double M) {
        this.M = M;
    }

    public final double getS() {
        return S;
    }

    public final void setS(double S) {
        this.S = S;
    }

    public double pdf(double x) {
        return pdf(x, M, S);
    }

    public double logPdf(double x) {
        return logPdf(x, M, S);
    }

    public double cdf(double x) {
        return cdf(x, M, S);
    }

    public double quantile(double y) {
        return quantile(y, M, S);
    }

    public double mean() {
        return mean(M, S);
    }

    public double variance() {
        return variance(M, S);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
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
     * probability density function
     *
     * @param x argument
     * @param M log mean
     * @param S log standard deviation
     * @return pdf at x
     */
    public static double pdf(double x, double M, double S) {
        return NormalDistribution.pdf(Math.log(x), M, S) / x;
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @param M log mean
     * @param S log standard deviation
     * @return log pdf at x
     */
    public static double logPdf(double x, double M, double S) {
        return NormalDistribution.logPdf(Math.log(x), M, S) - Math.log(x);
    }

    /**
     * cumulative density function
     *
     * @param x argument
     * @param M log mean
     * @param S log standard deviation
     * @return cdf at x
     */
    public static double cdf(double x, double M, double S) {
        return NormalDistribution.cdf(Math.log(x), M, S);
    }

    /**
     * quantiles (=inverse cumulative density function)
     *
     * @param z argument
     * @param M log mean
     * @param S log standard deviation
     * @return icdf at z
     */
    public static double quantile(double z, double M, double S) {
        return Math.exp(NormalDistribution.quantile(z, M, S));
    }

    /**
     * mean
     *
     * @param M log mean
     * @param S log standard deviation
     * @return mean
     */
    public static double mean(double M, double S) {
        return Math.exp(M + (S * S / 2));
    }

    /**
     * variance
     *
     * @param M log mean
     * @param S log standard deviation
     * @return variance
     */
    public static double variance(double M, double S) {
        double S2 = S * S;

        return Math.exp(S2 + 2 * M) * (Math.exp(S2) - 1);
    }

    // Private

    protected double M, S;


    public static void main(String[] args) {

        final LogNormalDistribution f = new LogNormalDistribution(1, 1);
        for (double i = 0.01; i < 1; i += 0.01) {
            final double y = i;
            BisectionZeroFinder zeroFinder = new BisectionZeroFinder(new OneVariableFunction() {
                public double value(double x) {
                    return f.cdf(x) - y;
                }
            }, 0, 100);
            zeroFinder.evaluate();
            System.err.println("" + f.quantile(i) + "\t" + zeroFinder.getResult());
        }
    }
}
