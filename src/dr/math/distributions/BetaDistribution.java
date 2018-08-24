/*
 * BetaDistribution.java
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

import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;
import org.apache.commons.math.special.Beta;
import org.apache.commons.math.special.Gamma;

/**
 * User: dkuh004
 * Date: Mar 25, 2011
 * Time: 11:32:25 AM
 */
public class BetaDistribution extends AbstractContinuousDistribution implements Distribution {

    // Default inverse cumulative probability accurac
    public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

    // first shape parameter
    private double alpha;

    // second shape parameter
    private double beta;

    // Normalizing factor used in density computations. updated whenever alpha or beta are changed.
    private double z;

    // Inverse cumulative probability accuracy
    private final double solverAbsoluteAccuracy;


    /**
     * This general constructor creates a new beta distribution with a
     * specified mean and scale
     *
     * @param alpha   shape parameter
     * @param beta    shape parameter
     */
    public BetaDistribution(double alpha, double beta) {
        this.alpha = alpha;
        this.beta = beta;
        z = Double.NaN;
        solverAbsoluteAccuracy = DEFAULT_INVERSE_ABSOLUTE_ACCURACY;
    }


    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }

    // Recompute the normalization factor.
    private void recomputeZ() {
        if (Double.isNaN(z)) {
            z = Gamma.logGamma(alpha) + Gamma.logGamma(beta) - Gamma.logGamma(alpha + beta);
        }
    }

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        recomputeZ();
        if (x < 0 || x > 1) {
            return 0;
        } else if (x == 0) {
            if (alpha < 1) {
                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                        "Cannot compute beta density at 0 when alpha = {0,number}", alpha);
            }
            return 0;
        } else if (x == 1) {
            if (beta < 1) {
                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                        "Cannot compute beta density at 1 when beta = %.3g", beta);
            }
            return 0;
        } else {
            double logX = Math.log(x);
            double log1mX = Math.log1p(-x);
            return Math.exp((alpha - 1) * logX + (beta - 1) * log1mX - z);
        }
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x){
        recomputeZ();
        if (x < 0 || x > 1) {
            return Double.NEGATIVE_INFINITY;
        } else if (x == 0) {
            if (alpha < 1) {
                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                        "Cannot compute beta density at 0 when alpha = {0,number}", alpha);
            }
            if (alpha == 1) {
                return 0;
            }
            return Double.NEGATIVE_INFINITY;
        } else if (x == 1) {
            if (beta < 1) {
                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                        "Cannot compute beta density at 1 when beta = %.3g", beta);
            }
            if (beta == 1) {
                return 0;
            }
            return Double.NEGATIVE_INFINITY;
        } else {
            double logX = Math.log(x);
            double log1mX = Math.log1p(-x);
            return (alpha - 1) * logX + (beta - 1) * log1mX - z;
        }
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y){
        if (y == 0) {
            return 0;
        } else if (y == 1) {
            return 1;
        } else {
            try{
                return super.inverseCumulativeProbability(y);
            } catch (MathException e) {
//                throw MathRuntimeException.createIllegalArgumentException(                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;

//                    "Couldn't calculate beta quantile for alpha = " + alpha + ", beta = " + beta + ": " +e.getMessage());
            }
        }
    }

    @Override
    protected double getInitialDomain(double p) {
        return p;
    }

    @Override
    protected double getDomainLowerBound(double p) {
        return 0;
    }

    @Override
    protected double getDomainUpperBound(double p) {
        return 1;
    }

    public double cdf(double x)  {
        if (x <= 0) {
            return 0;
        } else if (x >= 1) {
            return 1;
        } else {
            try {
                return Beta.regularizedBeta(x, alpha, beta);
            } catch (MathException e) {
                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
                // exceptions is bad. Better to return NaN and let the calling code deal with it.
                return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                "Couldn't calculate beta cdf for alpha = " + alpha + ", beta = " + beta + ": " +e.getMessage());
            }
        }

    }

    public double cumulativeProbability(double x) throws MathException {
        if (x <= 0) {
            return 0;
        } else if (x >= 1) {
            return 1;
        } else {
            return Beta.regularizedBeta(x, alpha, beta);
        }
    }

    @Override
    public double cumulativeProbability(double x0, double x1) throws MathException {
        return cumulativeProbability(x1) - cumulativeProbability(x0);
    }

    //Return the absolute accuracy setting of the solver used to estimate inverse cumulative probabilities.
    protected double getSolverAbsoluteAccuracy() {
        return solverAbsoluteAccuracy;
    }


    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
//    public double cdf(double x){
//        throw new UnsupportedOperationException();
//    }


    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean(){
        return (alpha / (alpha + beta));
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance(){
        return (alpha * beta) / ((alpha + beta)* (alpha + beta) * ( alpha + beta + 1) );
    }

    /**
     * @return a probability density function representing this distribution
     */
    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
        }

        public final double getUpperBound() {
            return 1.0;
        }
    };


}


