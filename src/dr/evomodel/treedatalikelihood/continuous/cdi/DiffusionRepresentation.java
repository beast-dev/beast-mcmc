/*
 * PrecisionType.java
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.evomodel.treedatalikelihood.preorder.PreOrderPrecision;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.matrix.SparseCompressedMatrix;
import dr.matrix.SparseSquareUpperTriangular;

/**
 * @author Marc A. Suchard
 */
public interface DiffusionRepresentation {

    enum Type {
        DENSE,
        SPARSE
    }

    Type getType();

    int getDimension();

    double[][] getVarianceCholeskyDecomposition(int index);

    SparseSquareUpperTriangular getPrecisionCholeskyDecomposition(int index);

    double[] getVariance(int index);

    PreOrderPrecision getPrecision(int index);

    void fill(double[] src, int denseOffset, int denseLength,
              SparseCompressedMatrix mat, int sparseIndex,
              double scalar);

    class Dense implements DiffusionRepresentation {

        public Dense(double[] matrix, int count, int dim) {
            this.matrix = matrix;
            this.count = count;
            this.dim = dim;
        }

        public final double[] matrix;
        public final int dim;
        public final int count;

        private double[][] variances;
        private double[][][] choleskyDecompositions;

        @Override
        public double[][] getVarianceCholeskyDecomposition(int index) {
            if (choleskyDecompositions == null) {
                choleskyDecompositions = new double[count][][];
            }

            if (choleskyDecompositions[index] == null) {
                choleskyDecompositions[index] = getCholeskyOfVariance(getVariance(index), dim);
            }

            return choleskyDecompositions[index];
        }

        @Override
        public SparseSquareUpperTriangular getPrecisionCholeskyDecomposition(int index) {
            return null;
        }

        @Override
        public double[] getVariance(int index) {
            if (variances == null) {
                variances = new double[count][];
            }

            if (variances[index] == null) {
                variances[index] = computeVarianceFromPrecision(matrix, index, dim);
            }

            return variances[index];
        }

        @Override
        public PreOrderPrecision getPrecision(int index) {
            return new PreOrderPrecision.Dense(
                    matrix, index, dim);
        }

        private static double[] computeVarianceFromPrecision(double[] matrix, int index, int dim) {
            final double[] precision;
            if (matrix.length != dim * dim) {
                double[] buffer = new double[dim * dim];
                System.arraycopy(matrix, dim * dim * index, buffer, 0, dim * dim);
                precision = buffer;
            } else {
                precision = matrix;
            }
            return new SymmetricMatrix(precision, dim).inverse().toArrayComponents();
        }

        public static double[][] getCholeskyOfVariance(double[] variance, final int dim) {
            return CholeskyDecomposition.execute(variance, 0, dim);
        }

        @Override
        public Type getType() {
            return Type.DENSE;
        }

        @Override
        public int getDimension() {
            return matrix.length;
        }

        @Override
        public void fill(double[] src, int denseOffset, int denseLength,
                         SparseCompressedMatrix mat, int sparseIndex,
                         double scalar) {
            for (int i = 0; i < denseLength; ++i) {
                matrix[i] = scalar * src[denseOffset + i];
            }
        }
    }

    class Sparse implements DiffusionRepresentation {

        private final PreOrderPrecision.Sparse[] precisions;

        public Sparse(int count) {
            this.precisions = new PreOrderPrecision.Sparse[count];
        }

        public Sparse(SparseCompressedMatrix sparsePrecisionMatrix,
                      SparseSquareUpperTriangular precisionCholeskyDecomposition,
                      double scalar) {
            this.precisions = new PreOrderPrecision.Sparse[] {
                    new PreOrderPrecision.Sparse(sparsePrecisionMatrix, precisionCholeskyDecomposition, scalar)
            };
        }

        @Override
        public Type getType() {
            return Type.SPARSE;
        }

        @Override
        public int getDimension() {
            return -1;
        }

        @Override
        public double[][] getVarianceCholeskyDecomposition(int index) {
            return null;
        }

        public SparseSquareUpperTriangular getPrecisionCholeskyDecomposition(int index) {
            return precisions[index].getCholeskyDecomposition();
        }

        @Override
        public double[] getVariance(int index) {
            return null;
        }

        @Override
        public PreOrderPrecision getPrecision(int index) {
            return precisions[index];
        }

        @Override
        public void fill(double[] src, int denseOffset, int denseLength,
                         SparseCompressedMatrix mat, int sparseIndex,
                         double scalar) {
            precisions[sparseIndex] = new PreOrderPrecision.Sparse(
                    mat, null, scalar);
        }
    }
}
