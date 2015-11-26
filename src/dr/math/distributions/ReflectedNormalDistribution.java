/*
 * ReflectedNormalDistribution.java
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
import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

/**
 * reflected normal distribution (pdf, cdf, quantile)
 *
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class ReflectedNormalDistribution implements Distribution {
    //
    // Public stuff
    //

    double lower;
    double upper;
    double precision;
    double m;
    double sd;

    /**
     * Constructor
     */
    public ReflectedNormalDistribution(double mean, double sd, double lower, double upper, double precision) {
        this.m = mean;
        this.sd = sd;
        this.lower = lower;
        this.upper = upper;
        this.precision = precision;
    }

    public double getMean() {
        return m;
    }

    public void setMean(double value) {
        m = value;
    }

    public double getSD() {
        return sd;
    }

    public void setSD(double value) {
        sd = value;
    }

    public double pdf(double x) {

        if (x < lower) return 0;
        if (x > upper) return 0;

        double pdf = NormalDistribution.pdf(x, m, sd);
        double newPDF = pdf;

        System.out.println("N(" + x + ",m,sd)=" + pdf);

        int i = 1;
        do {
            pdf = newPDF;

            int A = 2 * ((i + 1) / 2);       // 2  2  4  4  6  6
            int B = 2 * (i / 2);           // 0  2  2  4  4  6
            int C = i % 2 == 0 ? -1 : 1; // 1 -1  1 -1  1 -1

            double leftPDF = NormalDistribution.pdf(A * lower - C * x - B * upper, m, sd);
            double rightPDF = NormalDistribution.pdf(A * upper - C * x - B * lower, m, sd);
            newPDF = leftPDF + rightPDF + pdf;

            i += 1;

            System.out.println("newPDF=" + newPDF + " A=" + A + " B=" + B + " C=" + C);


        } while (newPDF - pdf > precision);

        return newPDF;
    }

    public double logPdf(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public double cdf(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public double quantile(double y) {
        throw new RuntimeException("Not implemented!");
    }

    public double mean() {
        throw new RuntimeException("Not implemented!");
    }

    public double variance() {
        throw new RuntimeException("Not implemented!");
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return lower;
        }

        public final double getUpperBound() {
            return upper;
        }
    };

    public static void main(String[] args) {

//        final ReflectedNormalDistribution rnd = new ReflectedNormalDistribution(2, 2, 0.5, 2, 1e-6);
        final ReflectedNormalDistribution rnd = new ReflectedNormalDistribution(2, 2, 1, 2, 1e-6);

        rnd.pdf(1);

        UnivariateRealFunction f = new UnivariateRealFunction() {
            public double value(double v) throws FunctionEvaluationException {
                return rnd.pdf(v);
            }
        };

        final UnivariateRealIntegrator integrator = new RombergIntegrator();
        integrator.setAbsoluteAccuracy(1e-14);
        integrator.setMaximalIterationCount(16);  // fail if it takes too much time

        double x;
        try {
            x = integrator.integrate(f, rnd.lower, rnd.upper);
//                      ptotErr += cdf != 0.0 ? Math.abs(x-cdf)/cdf : x;
//                      np += 1;
            //assertTrue("" + shape + "," + scale + "," + value + " " + Math.abs(x-cdf)/x + "> 1e-6", Math.abs(1-cdf/x) < 1e-6);
            System.out.println("Integrated pdf = " + x);

            //System.out.println(shape + ","  + scale + " " + value);
        } catch (ConvergenceException e) {
            // can't integrate , skip test
            //  System.out.println(shape + ","  + scale + " skipped");
        } catch (FunctionEvaluationException e) {

            throw new RuntimeException("I have no idea why I am here.");

        }

    }

}