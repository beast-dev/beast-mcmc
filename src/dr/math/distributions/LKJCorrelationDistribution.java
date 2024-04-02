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

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.extractUpperTriangular;

/**
 * @author Paul Bastide
 */
public class LKJCorrelationDistribution extends AbstractLKJDistribution {

    public static final String TYPE = "LKJCorrelation";

    public LKJCorrelationDistribution(int dim, double shape) {
        super(dim, shape);
    }

    public LKJCorrelationDistribution(int dim) { // returns a non-informative (uniform) density
        super(dim);
    }

    public double logPdf(double[] x) {

        assert (x.length == upperTriangularSize(dim));

//        if (shape == 1.0) { // Uniform //even when it's uniform, you still want to return -inf if it's not pos. def.
//            return logNormalizationConstant;
//        } else {
        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(x, dim);
        return logPdf(R);
//        }
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

    private double[] gradLogPdf(double[] x) {

        assert (x.length == upperTriangularSize(dim));

        if (shape == 1.0) { // Uniform
            return new double[x.length];
        } else {
            SymmetricMatrix R = compoundCorrelationSymmetricMatrix(x, dim);
            return gradLogPdf(R, shape);
        }
    }

    private static double[] gradLogPdf(SymmetricMatrix R, double shape) {

        double[] gradient = extractUpperTriangular((SymmetricMatrix) R.inverse());
        for (int i = 0; i < gradient.length; ++i) {
            gradient[i] = 2 * (shape - 1) * gradient[i];
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
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

}
