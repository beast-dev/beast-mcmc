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
import dr.matrix.SparseSquareUpperTriangular;
import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

/**
 * @author Marc A. Suchard
 */
public class NormalSufficientStatistics {

    public interface Precision {

        double getPrecision(int row, int col);

        double getVariance(int row, int col);

        DenseMatrix64F getRawPrecision();

        DenseMatrix64F getRawVariance();

        DenseMatrix64F getRawPrecisionCopy();

        DenseMatrix64F getRawVarianceCopy();

        String toVectorizedString();

        class Dense implements Precision {

            final DenseMatrix64F precision;
            DenseMatrix64F variance = null;

            public Dense(double[] buffer,
                         int index,
                         int dim,
                         Precision Pd,
                         PrecisionType precisionType) {
                this(DenseMatrix64F.wrap(dim, dim,
                        precisionType.getScaledPrecision(buffer,
                                (precisionType.getPartialsDimension(dim)) * index,
                                Pd.getRawPrecision().data, dim)));
            }

            public Dense(double[] buffer,
                         int index,
                         int dim,
                         DenseMatrix64F Pd,
                         PrecisionType precisionType) {
                this(DenseMatrix64F.wrap(dim, dim,
                        precisionType.getScaledPrecision(buffer,
                                (precisionType.getPartialsDimension(dim)) * index,
                                Pd.data, dim)));
            }

            public Dense(double[] matrices,
                         int index,
                         int dim) {
                this(MissingOps.wrap(matrices,  (dim * dim) * index, dim, dim));
            }

            public Dense(DenseMatrix64F precision) {
                this.precision = precision;
            }

            public Dense(DenseMatrix64F precision, DenseMatrix64F variance) {
                this.precision = precision;
                this.variance = variance;
            }

            @Override
            public String toVectorizedString() {
                String string = NormalSufficientStatistics.toVectorizedString(precision);
                if (variance != null) {
                    string += " " + NormalSufficientStatistics.toVectorizedString(variance);
                }
                return string;
            }

            @Override
            public double getPrecision(int row, int col) {
                return precision.unsafe_get(row, col);
            }

            @Override
            public double getVariance(int row, int col) {
                computeVariance();
                return variance.unsafe_get(row, col);
            }

            @Override
            public DenseMatrix64F getRawPrecision() { return precision; }

            @Override
            public DenseMatrix64F getRawVariance() {
                computeVariance();
                return variance;
            }

            @Override
            public DenseMatrix64F getRawPrecisionCopy() {
                return precision.copy();
            }

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
        }

        class Sparse implements Precision {

            @Override
            public double getPrecision(int row, int col) {
                throw new RuntimeException("Not implemented");
            }

            @Override
            public double getVariance(int row, int col) {
                throw new RuntimeException("Not implemented");
            }

            @Override
            public DenseMatrix64F getRawPrecision() {
                return null;
            }

            @Override
            public DenseMatrix64F getRawVariance() {
                return null;
            }

            @Override
            public DenseMatrix64F getRawPrecisionCopy() {
                return null;
            }

            @Override
            public DenseMatrix64F getRawVarianceCopy() {
                return null;
            }

            @Override
            public String toVectorizedString() {
                return null;
            }
        }
    }

    private final DenseMatrix64F mean;
    private final Precision precision;

    NormalSufficientStatistics(double[] buffer,
                               int index,
                               int dim,
                               Precision precision,
                               PrecisionType precisionType) {
        int partialOffset = (precisionType.getPartialsDimension(dim)) * index;
        this.mean = MissingOps.wrap(buffer, partialOffset, dim, 1);
        this.precision = new Precision.Dense(buffer, index, dim, precision, precisionType);
    }

//    NormalSufficientStatistics(double[] buffer,
//                               int index,
//                               int dim,
//                               DenseMatrix64F Pd,
//                               SparseSquareUpperTriangular cPd,
//                               PrecisionType precisionType) {
//
//        int partialOffset = (precisionType.getPartialsDimension(dim)) * index;
//        this.mean = MissingOps.wrap(buffer, partialOffset, dim, 1);
//        if (Pd != null) {
//            this.precision = new Precision.Dense(buffer, index, dim, Pd, precisionType);
//        } else if (cPd != null) {
//            this.precision = new Precision.Sparse();
//        } else {
//            throw new RuntimeException("Not yet implemented");
//        }
//    }

    NormalSufficientStatistics(double[] mean,
                               double[] precision,
                               int index,
                               int dim,
                               DenseMatrix64F Pd,
                               SparseSquareUpperTriangular cPd,
                               PrecisionType precisionType) {

        int meanOffset = dim * index;
        this.mean = MissingOps.wrap(mean, meanOffset, dim, 1);
        this.precision = new Precision.Dense(precision, index, dim);
    }

    public NormalSufficientStatistics(DenseMatrix64F mean,
                                      DenseMatrix64F precision,
                                      SparseSquareUpperTriangular cPd) {
        this.mean = mean;
        this.precision = new Precision.Dense(precision);
    }

    public NormalSufficientStatistics(DenseMatrix64F mean, DenseMatrix64F precision, DenseMatrix64F variance,
                                      SparseSquareUpperTriangular cPd) {
        this.mean = mean;
        this.precision = new Precision.Dense(precision, variance);
//        this.precision = precision;
//        this.variance = variance;
//        this.cPd = null;
    }

    public double getMean(int row) {
        return mean.get(row);
    }

//    public double getPrecision(int row, int col) {
//        return precision.unsafe_get(row, col);
//    }

//    public double getVariance(int row, int col) {
//        computeVariance();
//        return variance.unsafe_get(row, col);
//    }

    public SparseSquareUpperTriangular getSparseCholeskyPrecision() { return null; } //cPd.getMatrix(); }

    public double getSparseCholeskyPrecisionScalar() { return 0.0; } //cPd.getScalar(); }

    @Deprecated
    public DenseMatrix64F getRawPrecision() { return precision.getRawPrecision(); }

    @Deprecated
    public DenseMatrix64F getRawMean() { return mean; }

    @Deprecated
    public DenseMatrix64F getRawVariance() { return precision.getRawVariance(); }

    public DenseMatrix64F getRawPrecisionCopy() { return precision.getRawPrecisionCopy(); }

    public DenseMatrix64F getRawMeanCopy() { return mean.copy(); }

    public DenseMatrix64F getRawVarianceCopy() { return precision.getRawVarianceCopy(); }

    public String toString() {
        return mean + " " + precision.toString();
    }

    String toVectorizedString() {
        return toVectorizedString(mean.getData()) + " " + precision.toVectorizedString();
    }

//    String toVectorizedString() {
//        StringBuilder sb = new StringBuilder();
//        sb. append(toVectorizedString(mean.getData())).append(" ").append(toVectorizedString(precision.getData()));
//        if (variance != null) {
//            sb.append(" ").append(toVectorizedString(variance.getData()));
//        }
//        return sb.toString();
//    }

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
