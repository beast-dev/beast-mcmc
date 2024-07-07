/*
 * NormalSufficientStatistics.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

/**
 * @author Marc A. Suchard
 */
public class NormalSufficientStatistics {

    private final DenseMatrix64F mean;
    private final DenseMatrix64F precision;

    private DenseMatrix64F variance = null;

    NormalSufficientStatistics(double[] buffer,
                                      int index,
                                      int dim,
                                      DenseMatrix64F Pd,
                                      PrecisionType precisionType) {

        int partialOffset = (precisionType.getPartialsDimension(dim)) * index;
        this.mean = MissingOps.wrap(buffer, partialOffset, dim, 1);
        this.precision = DenseMatrix64F.wrap(dim, dim,
                precisionType.getScaledPrecision(buffer, partialOffset, Pd.data, dim));

    }

    @SuppressWarnings("unused")
    NormalSufficientStatistics(double[] mean,
                                      double[] precision,
                                      int index,
                                      int dim,
                                      DenseMatrix64F Pd,
                                      PrecisionType precisionType) {

        int meanOffset = dim * index;
        this.mean = MissingOps.wrap(mean, meanOffset, dim, 1);

        int precisionOffset = (dim * dim) * index;
//        this.precision = new DenseMatrix64F(dim, dim);
        this.precision = MissingOps.wrap(precision, precisionOffset, dim, dim);
//                DenseMatrix64F.wrap(dim, dim,
//                        precisionType.getScaledPrecision(precision, precisionOffset, Pd.data, dim));

    }

    @SuppressWarnings("unused")
    public NormalSufficientStatistics(DenseMatrix64F mean,
                                      DenseMatrix64F precision) {
        this.mean = mean;
        this.precision = precision;
    }

    public NormalSufficientStatistics(DenseMatrix64F mean, DenseMatrix64F precision, DenseMatrix64F variance) {
        this.mean = mean;
        this.precision = precision;
        this.variance = variance;
    }

    public double getMean(int row) {
        return mean.get(row);
    }

    public double getPrecision(int row, int col) {
        return precision.unsafe_get(row, col);
    }

    public double getVariance(int row, int col) {
        computeVariance();
        return variance.unsafe_get(row, col);
    }

    @Deprecated
    public DenseMatrix64F getRawPrecision() { return precision; }

    @Deprecated
    public DenseMatrix64F getRawMean() { return mean; }

    @Deprecated
    public DenseMatrix64F getRawVariance() {
        computeVariance();
        return variance;
    }

    public DenseMatrix64F getRawPrecisionCopy() {
        return precision.copy();
    }

    public DenseMatrix64F getRawMeanCopy() { return mean.copy(); }

    public DenseMatrix64F getRawVarianceCopy() {
        computeVariance();
        return variance.copy();
    }

    private void computeVariance() {
        if (variance == null) {
            variance = new DenseMatrix64F(precision.numRows, precision.numCols);
            safeInvert2(precision, variance, false);
        }
    }
    public String toString() {
        return mean + " " + precision;
    }

    String toVectorizedString() {
        StringBuilder sb = new StringBuilder();
        sb. append(toVectorizedString(mean.getData())).append(" ").append(toVectorizedString(precision.getData()));
        if (variance != null) {
            sb.append(" ").append(toVectorizedString(variance.getData()));
        }
        return sb.toString();
    }

    public static String toVectorizedString(DenseMatrix64F matrix) {
        return toVectorizedString(matrix.getData());
    }

    private static String toVectorizedString(double[] vector) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < vector.length - 1; ++i) {
            sb.append(vector[i]).append(" ");
        }
        sb.append(vector[vector.length - 1]);
        return sb.toString();
    }

}
