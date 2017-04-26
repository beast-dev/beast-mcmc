/*
 * DirichletDistribution.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.GammaFunction;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public class DirichletDistribution implements MultivariateDistribution {

    public static final String TYPE = "dirichletDistribution";
    public static final boolean DEBUG = false;

    //4.0 != 3.9999999999999996
    //Other BEAST classes uses 1E-8 or 1E-6
    public static final double ACCURACY_THRESHOLD = 1E-12;

    private double[] counts;
    private double countSum = 0.0;
    private double countParameterSum;
    private int dim;

    private boolean sumToNumberOfElements;

    private double logNormalizingConstant;

    public DirichletDistribution(double[] counts, boolean sumToNumberOfElements) {
        this.counts = counts;
        this.sumToNumberOfElements = sumToNumberOfElements;
        if (this.sumToNumberOfElements) {
            countParameterSum = (double)counts.length;
        } else {
            countParameterSum = 1.0;
        }
        dim = counts.length;
        for (int i = 0; i < dim; i++) {
            countSum += counts[i];
        }

        computeNormalizingConstant();
    }

    public DirichletDistribution(double[] counts, double countParameterSum) {
        this.counts = counts;
        this.countParameterSum = countParameterSum;

        dim = counts.length;
        for (int i = 0; i < dim; i++) {
            countSum += counts[i];
        }

        computeNormalizingConstant();
    }

    private void computeNormalizingConstant() {
        logNormalizingConstant = GammaFunction.lnGamma(countSum);
        for (int i = 0; i < dim; i++) {
            logNormalizingConstant -= GammaFunction.lnGamma(counts[i]);
        }
        logNormalizingConstant -= dim * Math.log(countParameterSum);
    }


    public double logPdf(double[] x) {

        if (x.length != dim) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }

        double logPDF = logNormalizingConstant;
        double parameterSum = 0.0;
        for (int i = 0; i < dim; i++) {
            logPDF += (counts[i] - 1) * (Math.log(x[i]) - Math.log(countParameterSum));
            parameterSum += x[i];
//            if ((!sumToNumberOfElements && x[i] >= 1.0) || x[i] <= 0.0) {
//                if (DEBUG) {
//                    System.out.println("Invalid parameter value");
//                }
//                logPDF = Double.NEGATIVE_INFINITY;
//                break;
//            }
        }
        if (Math.abs(parameterSum - countParameterSum) > ACCURACY_THRESHOLD) {
            if (DEBUG) {
                System.out.println("Parameters do not sum to " + countParameterSum);
                for (int i = 0; i < dim; i++) {
                    System.out.println("x[" + i + "] = " + x[i]);
                }
                System.out.println("Current sum = " + parameterSum);
            }
            logPDF = Double.NEGATIVE_INFINITY;
        }

        return logPDF;
    }

    public double[][] getScaleMatrix() {
        return null;
    }

    public double[] getMean() {
        double[] mean = new double[dim];
        for (int i = 0; i < dim; i++)
            mean[i] = counts[i] / countSum;
        return mean;
    }

    public String getType() {
        return TYPE;
    }

    public static void main(String[] args) {

        //Test Dirichlet distribution for the standard n-simplex
        System.out.println("Test Dirichlet distribution for the standard n-simplex");
        //R: log(ddirichlet(c(0.5,0.2,0.3),c(1,2,3))) = 0.07696104
        double[] counts = new double[3];
        counts[0] = 1.0;
        counts[1] = 2.0;
        counts[2] = 3.0;
        DirichletDistribution dd = new DirichletDistribution(counts, false);
        double[] parameterValues = new double[3];
        parameterValues[0] = 0.5;
        parameterValues[1] = 0.2;
        parameterValues[2] = 0.3;
        System.out.println(dd.logPdf(parameterValues));

        //Test Scaled Dirichlet distribution
        System.out.println("Test Scaled Dirichlet distribution");
        //R: log(ddirichlet(c(1.5,0.6,0.9)/3,c(1,2,3))/(3^3)) = -3.218876
        dd = new DirichletDistribution(counts, true);
        parameterValues[0] = 1.5;
        parameterValues[1] = 0.6;
        parameterValues[2] = 0.9;
        System.out.println(dd.logPdf(parameterValues));

        parameterValues[0] = 1.0;
        parameterValues[1] = 1.0;
        parameterValues[2] = 1.0;
        System.out.println(dd.logPdf(parameterValues));

        counts = new double[4];
        counts[0] = 1.0;
        counts[1] = 1.0;
        counts[2] = 1.0;
        counts[3] = 1.0;
        dd = new DirichletDistribution(counts, true);
        parameterValues = new double[4];
        parameterValues[0] = 0.5;
        parameterValues[1] = 1.2;
        parameterValues[2] = 1.3;
        parameterValues[3] = 1.0;
        System.out.println(dd.logPdf(parameterValues));
        parameterValues[0] = 1.0;
        parameterValues[1] = 1.0;
        parameterValues[2] = 1.0;
        parameterValues[3] = 1.0;
        System.out.println(dd.logPdf(parameterValues));

    }

}
