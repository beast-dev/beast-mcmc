/*
 * RiemannApproximation.java
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

import dr.math.distributions.NormalDistribution;

/**
 * Approximates the integral of a given function using Riemann integration
 *
 * @author Alexei Drummond
 * @version $Id: RiemannApproximation.java,v 1.3 2005/05/24 20:26:01 rambaut Exp $
 */
public class RiemannApproximation implements Integral {

    public enum variant {
        UPPER,
        LOWER,
        MIDPOINT,
        TRAPEZOID
    }

    public RiemannApproximation(int sampleSize, variant mode) {
        this.mode = mode;
        this.sampleSize = sampleSize;
    }

    public RiemannApproximation(int sampleSize) {
        this(sampleSize, variant.UPPER);
    }

    /**
     * @param f   the function whose integral is of interest
     * @param min the minimum value of the function
     * @param max the  upper limit of the function
     * @return the approximate integral of the given function
     *         within the given range.
     */
    public double integrate(UnivariateFunction f, double min, double max) {

        double integral = 0.0;

        double gridpoint = min;
        double step = (max - min) / sampleSize;

        switch (mode){
            case UPPER:
                for (int i = 1; i <= sampleSize; i++) {
                    gridpoint += step;
                    integral += f.evaluate(gridpoint);
                }
                break;
            case LOWER:
                for (int i = 1; i <= sampleSize; i++) {
                    integral += f.evaluate(gridpoint);
                    gridpoint += step;
                }
                break;
            case MIDPOINT:
                for (int i = 1; i <= sampleSize; i++) {
                    integral += f.evaluate(gridpoint + step/2);
                    gridpoint += step;
                }
                break;
            case TRAPEZOID:
                for (int i = 1; i <= sampleSize; i++) {
                    integral += (f.evaluate(gridpoint) + f.evaluate(gridpoint+step))/2;
                    gridpoint += step;
                }
                break;
        }
        integral *= (max - min) / (double) sampleSize;
        return integral;
    }

    // this does the addition on a log scale, in case of underflow

    public double logIntegrate(UnivariateFunction f, double min, double max){
        double logIntegral = Double.NEGATIVE_INFINITY;

        double gridpoint = min;
        double step = (max - min) / sampleSize;

        switch (mode){
            case UPPER:
                for (int i = 1; i <= sampleSize; i++) {
                    gridpoint += step;
                    logIntegral = LogTricks.logSum(logIntegral, Math.log(f.evaluate(gridpoint)));
                }
                break;
            case LOWER:
                for (int i = 1; i <= sampleSize; i++) {
                    logIntegral = LogTricks.logSum(logIntegral, Math.log(f.evaluate(gridpoint)));
                    gridpoint += step;
                }
                break;
            case MIDPOINT:
                for (int i = 1; i <= sampleSize; i++) {
                    logIntegral = LogTricks.logSum(logIntegral, Math.log(f.evaluate(gridpoint+step/2)));
                    gridpoint += step;
                }
                break;
            case TRAPEZOID:
                for (int i = 1; i <= sampleSize; i++) {
                    logIntegral = LogTricks.logSum(logIntegral,
                            Math.log((f.evaluate(gridpoint) + f.evaluate(gridpoint+step))/2));
                    gridpoint += step;
                }
                break;
        }
        logIntegral += Math.log((max - min) / (double) sampleSize);
        return logIntegral;
    }

    public variant getMode() {
        return mode;
    }


    public static void main(String[] args) {

        UnivariateFunction normalPDF = new NormalDistribution(0.0, 1.0).getProbabilityDensityFunction();
        UnivariateFunction normalPDF2 = new NormalDistribution(0.0, 1.0).getProbabilityDensityFunction();
        UnivariateFunction normalPDF3 = new NormalDistribution(0.0, 1.0).getProbabilityDensityFunction();

        double Z = 1.0;
        //double Z = Math.sqrt(2*Math.PI);
        CompoundFunction threeNormals = new CompoundFunction(new UnivariateFunction[]{normalPDF, normalPDF2, normalPDF3}, Z);

        System.out.println("Riemann approximation to the integral of a three normal distribution:");
        RiemannApproximation integrator = new RiemannApproximation(100000);
        System.out.println("integral(N(0.0, 1.0))=" + integrator.integrate(normalPDF, -4.0, 4.0));
        System.out.println("integral(N(1.0, 2.0))=" + integrator.integrate(normalPDF2, -8.0, 8.0));
        System.out.println("integral(N(2.0, 3.0))=" + integrator.integrate(normalPDF3, -16.0, 16.0));


        double integral = integrator.integrate(threeNormals, -16.0, 16.0);
        System.out.println("Riemann approximation to the integral of the compound of three normal distribution:");
        System.out.println("integral(N(0.0, 1.0)*N(1.0, 2.0)*N(2.0, 3.0))=" + integral);
        System.out.println("Estimate normalizing constant is " + (1.0 / integral));

        /*System.out.println("Ten monte carlo approximations to the integral of a normal distribution:");
          MonteCarloIntegral integrator2 = new MonteCarloIntegral(10000);
          for (int i = 0; i < 10; i++) {
              System.out.println(integrator2.integrate(normalPDF, -4.0, 4.0));
          }*/
    }

    public int getSampleSize(){
        return sampleSize;
    }

    private int sampleSize;

    private final variant mode;
}
