/*
 * DirichletDistribution.java
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

import dr.math.GammaFunction;

/**
 * @author Marc A. Suchard
 */
public class DirichletDistribution implements MultivariateDistribution {

    public static final String TYPE = "dirichletDistribution";

    private double[] counts;
    private double countSum = 0.0;
    private int dim;

    private double logNormalizingConstant;

    public DirichletDistribution(double[] counts) {
        this.counts = counts;
        dim = counts.length;
        for (int i = 0; i < dim; i++)
            countSum += counts[i];

        computeNormalizingConstant();
    }

    private void computeNormalizingConstant() {
        logNormalizingConstant = GammaFunction.lnGamma(countSum);
        for (int i = 0; i < dim; i++)
            logNormalizingConstant -= GammaFunction.lnGamma(counts[i]);
    }


    public double logPdf(double[] x) {

        if (x.length != dim) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }

        double logPDF = logNormalizingConstant;
        for (int i = 0; i < dim; i++) {
            logPDF += (counts[i] - 1) * Math.log(x[i]);
            if (x[i] <= 0.0 || x[i] >= 1.0) {
                logPDF = Double.NEGATIVE_INFINITY;
                break;
            }
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
}
