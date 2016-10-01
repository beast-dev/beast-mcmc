/*
 * PoissonDistribution.java
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

import dr.math.Poisson;
import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

/**
 * @author Alexei Drummond
 * @version $Id$
 */
public class PoissonDistribution implements Distribution {

    org.apache.commons.math.distribution.PoissonDistribution distribution;

    public PoissonDistribution(double mean) {
        distribution = new org.apache.commons.math.distribution.PoissonDistributionImpl(mean);
    }

    public double pdf(double x) {
        return distribution.probability(x);
    }

    public double logPdf(double x) {

        double pdf = distribution.probability(x);
        if (pdf == 0 || Double.isNaN(pdf)) { // bad estimate
            final double mean = mean();
            return x * Math.log(mean) - Poisson.gammln(x + 1) - mean;
        }
        return Math.log(pdf);

    }

    public double cdf(double x) {
        try {
            return distribution.cumulativeProbability(x);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double quantile(double y) {
        try {
            return distribution.inverseCumulativeProbability(y);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double mean() {
        return distribution.getMean();
    }

    public double variance() {
        return distribution.getMean();
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException();
    }

    public double truncatedMean(int max) {

        double CDF = 0;
        double mean = 0;
        for(int i=0; i<=max; i++) {
            double p = distribution.probability(i);
            mean += i*p;
            CDF += p;
        }
        return mean / CDF;        
    }

    public static double pdf(double x, double mean) {
        PoissonDistributionImpl dist = new PoissonDistributionImpl(mean);
        return dist.probability(x);
    }

    public static double logPdf(double x, double mean) {
        PoissonDistributionImpl dist = new PoissonDistributionImpl(mean);
        double pdf = dist.probability(x);
        if (pdf == 0 || Double.isNaN(pdf)) { // bad estimate
            return x * Math.log(mean) - Poisson.gammln(x + 1) - mean;
        }
        return Math.log(pdf);

    }

    public static double cdf(double x, double mean) {
        try {
            PoissonDistributionImpl dist = new PoissonDistributionImpl(mean);
            return dist.cumulativeProbability(x);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public static double quantile(double y, double mean) {
        try {
            PoissonDistributionImpl dist = new PoissonDistributionImpl(mean);
            return dist.inverseCumulativeProbability(y);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }
}
