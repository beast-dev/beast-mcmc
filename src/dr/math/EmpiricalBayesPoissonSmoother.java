/*
 * EmpiricalBayesPoissonSmoother.java
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

import dr.math.distributions.GammaDistribution;
import dr.stats.DiscreteStatistics;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 * @author Philippe Lemey
 */

public class EmpiricalBayesPoissonSmoother {

    /**
     * Provides an empirical Bayes estimate of counts under a Poisson sampling density with a Gamma prior
     * on the unknown intensity.  The marginal distribution of the data is then a Negative-Binomial.
     *
     * @param in vector of Poisson counts
     * @return smoothed count estimates
     */
    public static double[] smooth(double[] in) {
        final int length = in.length;
        double[] out = new double[length];
        double[] gammaStats = getNegBin(in);
        double alpha = gammaStats[0];
        double beta = gammaStats[1]; // As defined on wiki page (scale)
        double mean = gammaStats[2];

        if (beta == 0) {
            for (int i = 0; i < length; i++) {
                out[i] = mean;
            }
        } else {
            for (int i = 0; i < length; i++) {
                out[i] = (in[i] + alpha) / (1 + 1 / beta);
            }
        }
        return out;
    }

    public static double[] smoothWithSample(double[] in) {
        final int length = in.length;
        double[] out = new double[length];
        double[] gammaStats = getNegBin(in);
        double alpha = gammaStats[0];
        double beta = gammaStats[1]; // As defined on wiki page (scale)
        double mean = gammaStats[2];

        if (beta == 0) {
            for (int i = 0; i < length; i++) {
                out[i] = mean;
            }
        } else {
            for (int i = 0; i < length; i++) {
                double shape = in[i] + alpha;
                double scale = 1 / (1 + 1 / beta);
                out[i] = GammaDistribution.nextGamma(shape, scale);
            }
        }
        return out;
    }

    public static double[] smoothOld(double[] in) {
        final int length = in.length;
        double[] out = new double[length];
        double[] gammaStats = getNegBin(in);
        for (int i = 0; i < length; i++) {
            out[i] = (in[i] + gammaStats[0]) / (1 + 1 / gammaStats[1]);
        }
        return out;
    }

    // Method of moments estimators following Martiz 1969

    private static double[] getNegBin(double[] array) {
        double mean = DiscreteStatistics.mean(array);
        double variance = DiscreteStatistics.variance(array, mean);
        double returnArray0 = (1 - (mean / variance));
        double returnArray1 = (mean * ((1 - returnArray0) / returnArray0));

        double shape = returnArray1;
        double scale = (returnArray0 / (1 - returnArray0));

        if (variance <= mean) {
            shape = 0.0;
            scale = 0.0;
        }

//        // Check against Martiz 1969 (beta = shape, alpha = rate in the 1969 paper)
//        double matrizBeta = mean * mean / (variance - mean);
//        double matrizAlphaInv = mean / matrizBeta; // scale
//        System.err.println("mb = " + matrizBeta + " shape = " + shape);
//        System.err.println("ma = " + matrizAlphaInv + " scale = " + scale);

        return new double[]{shape, scale, mean};
    }
}
