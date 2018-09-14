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
import dr.util.CorrelationToCholesky;

import static dr.math.matrixAlgebra.WrappedMatrix.WrappedUpperTriangularMatrix.fillDiagonal;


/**
 * @author Paul Bastide
 */

public class LKJCholeskyCorrelationDistribution extends AbstractLKJDistribution {

    public static final String TYPE = "LKJCholeskyCorrelation";
    private CorrelationToCholesky correlationToCholeskyTransform;

    public LKJCholeskyCorrelationDistribution(int dim, double shape) {
        super(dim, shape);
        this.correlationToCholeskyTransform = new CorrelationToCholesky(dim);
    }

    public LKJCholeskyCorrelationDistribution(int dim) { // returns a non-informative (uniform) density
        super(dim);
        this.correlationToCholeskyTransform = new CorrelationToCholesky(dim);
    }

    public double logPdf(double[] x) { //x must be vechu(chol), so the V matrix.

        assert (x.length == upperTriangularSize(dim));
        double[] vechuCholesky = correlationToCholeskyTransform.transform(x, 0, x.length);
        WrappedMatrix.WrappedUpperTriangularMatrix L = fillDiagonal(vechuCholesky, dim);
        return logPdf(L);
    }

    private double logPdf(WrappedMatrix L) {
        // See Stan manual, p. 558
        // See also: http://mc-stan.org/math/d7/d74/lkj__corr__cholesky__lpdf_8hpp_source.html
        double shapeConst = 2 * shape - 2;
        double logDensity = 0.0;
        for (int i = 1; i < dim; i++) {
            logDensity += (dim - i - 1 + shapeConst) * Math.log(L.get(i, i));
        }
        logDensity += logNormalizationConstant;
        return logDensity;
    }

    private double[] gradLogPdf(double[] vechuCholesky) {

        assert (vechuCholesky.length == upperTriangularSize(dim));

        WrappedMatrix.WrappedUpperTriangularMatrix vechCholesky = fillDiagonal(vechuCholesky, dim);
        return gradLogPdf(vechCholesky, shape);
    }

    private static double[] gradLogPdf(WrappedMatrix L, double shape) {
        int dim = L.getMajorDim();
        double shapeConst = 2 * shape - 2;
        double[] gradient = new double[dim * (dim - 1) / 2];
        int k = 0;
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                gradient[k] = -(dim - j - 1 + shapeConst) * L.get(i, j) / Math.pow(L.get(j, j), 2);
                k++;
            }
        }
        return gradient;
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
    public double[] getGradientLogDensity(Object x) { //here x needs to be cholesky. well when the method was called,

        double[] cholesky = correlationToCholeskyTransform.transform((double[]) x, 0, ((double[]) x).length);
        return gradLogPdf(cholesky);
    }
}

