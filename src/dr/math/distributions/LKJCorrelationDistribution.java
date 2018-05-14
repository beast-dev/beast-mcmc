/*
 * LKJCorrelationDistribution.java
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

import dr.inference.model.GradientProvider;
import dr.math.GammaFunction;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.extractUpperTriangular;

/**
 * @author Paul Bastide
 */
public class LKJCorrelationDistribution implements MultivariateDistribution, GradientProvider {

    public static final String TYPE = "LKJCorrelation";

    private double shape;
    private int dim;
    private double logNormalizationConstant;

    public LKJCorrelationDistribution(int dim, double shape) {

        assert (shape < 0);

        this.shape = shape;
        this.dim = dim;
        this.logNormalizationConstant = computelogNormalizationConstant();

    }

    public LKJCorrelationDistribution(int dim) { // returns a non-informative (uniform) density

        this.shape = 1.0;
        this.dim = dim;
        this.logNormalizationConstant = computelogNormalizationConstant();

    }

    private double computelogNormalizationConstant() {
        // Lewandowski, Kurowicka, and Joe (2009)
        // See also Stan: http://mc-stan.org/math/db/d4f/lkj__corr__lpdf_8hpp_source.html
        // And: http://discourse.mc-stan.org/t/question-about-lkj-normalizing-constant/2001/11 (for the sign)
        double res = 0.0;
        if (shape == 1.0) {
            // Lewandowski et al. (2009) Theorem 5
            for (int k = 1; k <= (dim - 1) / 2; k++) {
                res -= GammaFunction.lnGamma(2.0 * k);
            }
            if ((dim % 2) == 1) {
                res -= 0.25 * (dim * dim - 1) * Math.log(Math.PI)
                        - 0.25 * (dim - 1) * (dim - 1) * Math.log(2.0)
                        - (dim - 1) * GammaFunction.lnGamma(0.5 * (dim + 1));
            } else {
                res -= 0.25 * dim * (dim - 2) * Math.log(Math.PI)
                        + 0.25 * (3 * dim * dim - 4 * dim) * Math.log(2.0)
                        + dim * GammaFunction.lnGamma(0.5 * dim)
                        - (dim - 1) * GammaFunction.lnGamma(dim);
            }
        } else {
            // Lewandowski et al. (2009), expression in proof of eq. (17)
            res = (dim - 1) * GammaFunction.lnGamma(shape + 0.5 * (dim - 1));
            for (int k = 1; k <= (dim - 1); k++) {
                res -= 0.5 * k * Math.log(Math.PI)
                        + GammaFunction.lnGamma(shape + 0.5 * (dim - 1 - k));
            }
        }
        return res;
    }

    public double logPdf(double[] x) { // x must be of length (2 choose dim) [upper triangular]
        if (shape == 1.0) { // Uniform
            return logNormalizationConstant;
        } else {
            SymmetricMatrix R = compoundCorrelationSymmetricMatrix(x, dim);
//            Matrix R = new Matrix(x, dim, dim);
            return logPdf(R);
        }
    }

    private double logPdf(Matrix R) {
        // See WishartDistribution.logPdf
        double logDensity = 0;

        try {
            logDensity = R.logDeterminant(); // Returns NaN is R is not positive-definite.

            if (Double.isInfinite(logDensity) || Double.isNaN(logDensity)) {
                return Double.NEGATIVE_INFINITY;
            }

            logDensity *= shape - 1;


        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        logDensity += logNormalizationConstant;
        return logDensity;
    }

    public double[] gradLogPdf(double[] x) { // x must be of length (2 choose dim) [upper triangular]
        if (shape == 1.0) { // Uniform
            return new double[x.length];
        } else {
            SymmetricMatrix R = compoundCorrelationSymmetricMatrix(x, dim);
            return gradLogPdf(R, shape);
        }
    }

    public static double[] gradLogPdf(SymmetricMatrix R, double shape) {

        double[] gradient = extractUpperTriangular((SymmetricMatrix) R.inverse());
        for (int i = 0; i < gradient.length; ++i) {
            gradient[i] = (shape - 1) * gradient[i];
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
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

}
