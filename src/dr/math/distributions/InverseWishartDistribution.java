/*
 * InverseWishartDistribution.java
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
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc Suchard
 */
public class InverseWishartDistribution implements MultivariateDistribution, WishartStatistics {

    public static final String TYPE = "InverseWishart";

    private double df;
    private int dim;
    private double[][] scaleMatrix;
    private Matrix S;
    private double logNormalizationConstant;

    /**
     * An Inverser Wishart distribution class for \nu degrees of freedom and scale matrix S with dim k
     * Expectation = (\nu - k - 1)^{-1} * S
     *
     * @param df
     * @param scaleMatrix
     */

    public InverseWishartDistribution(double df, double[][] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
        this.dim = scaleMatrix.length;

        S = new Matrix(scaleMatrix);
        computeNormalizationConstant();
    }

    private void computeNormalizationConstant() {
        logNormalizationConstant = 0;
        try {
            logNormalizationConstant = df / 2.0 * Math.log(new Matrix(scaleMatrix).determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        logNormalizationConstant -= df * dim / 2.0 * Math.log(2);
        logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
        for (int i = 1; i <= dim; i++) {
            logNormalizationConstant -= GammaFunction.lnGamma((df + 1 - i) / 2.0);
        }

    }


    public String getType() {
        return TYPE;
    }

    public double[][] getScaleMatrix() {
        return scaleMatrix;
    }

    public double[] getMean() {
        return null;
    }

    public double getDF() {
        return df;
    }

    public double[][] scaleMatrix() {
        return scaleMatrix;
    }


    public double logPdf(double[] x) {
        Matrix W = new Matrix(x, dim, dim);
        double logDensity = 0;

//	    System.err.println("here");
//	    double det = 0;
//	    try {
//	        det = W.determinant();
//	    }   catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();
//        }
//	    if( det < 0 ) {
//		    System.err.println("not positive definite");
//		    return Double.NEGATIVE_INFINITY;
//	    }


        try {
            logDensity = Math.log(W.determinant());

            logDensity *= -0.5;
            logDensity *= df + dim + 1;

            Matrix product = S.product(W.inverse());

            for (int i = 0; i < dim; i++)
                logDensity -= 0.5 * product.component(i, i);

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        logDensity += logNormalizationConstant;
        return logDensity;
    }

}
