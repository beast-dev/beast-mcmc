/*
 * BayesianBridgeRNG.java
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

package dr.inference.distribution.shrinkage;

import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

public class BayesianBridgeRNG {

    private static double drawUnregularizedCoefficient(double globalScale, double exponent) {
        // This is how rgnorm works in the R package gnorm
        double lambda = Math.pow((1.0 / globalScale),exponent);
        double u = MathUtils.nextDouble();
        double unsigned = Math.pow(GammaDistribution.quantile(u, 1.0 / exponent, 1.0 / lambda),1.0 / exponent);
        return MathUtils.nextBoolean() ? unsigned : -unsigned;
    }

    private static double drawRegularizedCoefficientBridgeProposal(double globalScale, double exponent, double slabWidth) {
        double x = Double.NaN;
        boolean done = false;
        double twoSlabSquared = 2.0 * slabWidth * slabWidth;
        while (!done) {
            double prop = drawUnregularizedCoefficient(globalScale, exponent);
            double logAcceptProb = -(prop * prop)/twoSlabSquared;
            if ( Math.log(MathUtils.nextDouble()) <= logAcceptProb) {
                x = prop;
                done = true;
            }
        }
        return x;
    }

    public static double drawRegularizedCoefficientNormalProposal(double globalScale, double exponent, double slabWidth) {
        double x = Double.NaN;
        boolean done = false;
        double logSlab2SqrtPi = Math.log(slabWidth * Math.sqrt(2.0 * Math.PI));
        while (!done) {
            double prop = MathUtils.nextGaussian() * slabWidth;
            double logAcceptProb = -Math.pow(Math.abs(prop / globalScale),exponent) - logSlab2SqrtPi;;
            if ( Math.log(MathUtils.nextDouble()) <= logAcceptProb) {
                x = prop;
                done = true;
            }
        }
        return x;
    }

    public static double drawRegularizedCoefficient(double globalScale, double exponent, double slabWidth) {
        double x;
        if ( globalScale > slabWidth ) {
            x = drawRegularizedCoefficientNormalProposal(globalScale, exponent, slabWidth);
        } else {
            x = drawRegularizedCoefficientBridgeProposal(globalScale, exponent, slabWidth);
        }
        return x;
    }

    // For Bridge without slab
    public static double[] nextRandom(double globalScale, double exponent, int dim) {
        double[] draws = new double[dim];
        for (int i = 0; i < dim; i++) {
            draws[i] = drawUnregularizedCoefficient(globalScale, exponent);
        }
        return draws;
    }

    // For Bridge with slab
    public static double[] nextRandom(double globalScale, double exponent, double slabWidth, int dim) {
        double[] draws = new double[dim];
        for (int i = 0; i < dim; i++) {
            draws[i] = drawRegularizedCoefficient(globalScale, exponent, slabWidth);
        }
        return draws;
    }
}
