/*
 * LKJCholeskyCorrelationDistribution.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.matrixAlgebra.WrappedMatrix;

import static dr.math.matrixAlgebra.WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal;


/**
 * @author Paul Bastide
 */

public class LKJCholeskyCorrelationDistribution extends AbstractLKJDistribution {

    public static final String TYPE = "LKJCholeskyCorrelation";

    public LKJCholeskyCorrelationDistribution(int dim, double shape) {
        super(dim, shape);
    }

    public LKJCholeskyCorrelationDistribution(int dim) { // returns a non-informative (uniform) density
        super(dim);
    }

    public double logPdf(double[] x) { // x must be of length dim*(dim-1)/2 [upper triangular]
        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(x, dim);
        return logPdf(L);
    }

    private double logPdf(WrappedMatrix L) {
        // See Stan manual, p. 558
        // See also: http://mc-stan.org/math/d7/d74/lkj__corr__cholesky__lpdf_8hpp_source.html
        double logDensity = 0.0;
        if (shape == 1) {
            for (int i = 1; i < dim; i++) {
                logDensity += (dim - i - 1) * Math.log(L.get(i, i));
            }
        } else {
            for (int i = 1; i < dim; i++) {
                logDensity += (dim - i - 1 + 2 * shape - 2) * Math.log(L.get(i, i));
            }
        }
        logDensity += logNormalizationConstant;
        return logDensity;
    }

    public double[] gradLogPdf(double[] x) { // x must be of length dim*(dim+1)/2 [upper triangular]
        WrappedMatrix.WrappedUpperTriangularMatrix L = new WrappedMatrix.WrappedUpperTriangularMatrix(x, dim);
        return gradLogPdf(L, shape);
    }

    public static double[] gradLogPdf(WrappedMatrix L, double shape) {
        int dim = L.getMajorDim();
        WrappedMatrix.WrappedUpperTriangularMatrix gradient = new WrappedMatrix.WrappedUpperTriangularMatrix(dim);
        if (shape == 1) {
            for (int i = 0; i < dim; ++i) {
                gradient.set(i, i, (dim - i - 1) / L.get(i, i));
            }
            return gradient.getBuffer();
        } else {
            for (int i = 0; i < dim; ++i) {
                gradient.set(i, i, (dim - i - 1 + 2 * shape - 2) / L.get(i, i));
            }
            return gradient.getBuffer();
        }
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }

    public String getType() {
        return TYPE;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

}

