/*
 * KernelDensityEstimatorDistribution.java
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
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.DifferentiableUnivariateRealFunction;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math.analysis.solvers.NewtonSolver;

/**
 * @author Marc Suchard
 */
public abstract class KernelDensityEstimatorDistribution implements Distribution {

    KernelDensityEstimatorDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {

        this.sample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.sample[i] = sample[i];
        }
        this.N = sample.length;
        processBounds(lowerBound, upperBound);
        setBandWidth(bandWidth);
    }

    abstract protected double evaluateKernel(double x);

    abstract protected void processBounds(Double lowerBound, Double upperBound);

    abstract protected void setBandWidth(Double bandWidth);

    abstract public double getFromPoint();

    abstract public double getToPoint();

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return evaluateKernel(x);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {

        double cdf;
        try {
            cdf = cdfFunction.value(x);
        } catch (FunctionEvaluationException e) {
            throw new RuntimeException(e.getMessage());
        }

        return cdf;
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return cdf value
     */
    public double quantile(final double y) {

        final DifferentiableUnivariateRealFunction root =
                new DifferentiableUnivariateRealFunction() {

            @Override
            public UnivariateRealFunction derivative() {
                return pdfFunction;
            }

            @Override
            public double value(double x) throws FunctionEvaluationException {

                return cdfFunction.value(x) - y;
            }
        };

        NewtonSolver solver = new NewtonSolver(root);

        double q;
        try {
            q = solver.solve(getFromPoint(), getToPoint());
        } catch (MaxIterationsExceededException e) {
            throw new RuntimeException(e.getMessage());
        } catch (FunctionEvaluationException e) {
            throw new RuntimeException(e.getMessage());
        }

        return q;
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * @return a probability density function representing this distribution
     */
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not Implemented.");
    }

    public double getBandWidth() {
        return bandWidth;
    }

    public enum Type {
        GAUSSIAN("Gaussian"),
        GAMMA("Gamma"),
        LOG_TRANSFORMED_GAUSSIAN("LogTransformedGaussian"),
        BETA("Beta");

        Type(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        public static Type parseFromString(String text) {
            for (Type format : Type.values()) {
                if (format.getText().compareToIgnoreCase(text) == 0)
                    return format;
            }
            return null;
        }

        private final String text;
    }

    protected int N;
    protected double lowerBound;
    protected double upperBound;
    double bandWidth;
    protected double[] sample;

    private final SimpsonIntegrator integrator = new SimpsonIntegrator();

    private final UnivariateRealFunction pdfFunction = new UnivariateRealFunction() {

        @Override
        public double value(double x) throws FunctionEvaluationException {
            return pdf(x);
        }
    };

    private final UnivariateRealFunction cdfFunction = new UnivariateRealFunction() {

        @Override
        public double value(double x) throws FunctionEvaluationException {
            double rangeMin = getFromPoint();

            double cdf;
            try {
                cdf = integrator.integrate(pdfFunction, rangeMin, x);
            } catch (MaxIterationsExceededException e) {
                throw new RuntimeException(e.getMessage());
            }

            return cdf;
        }
    };
}
