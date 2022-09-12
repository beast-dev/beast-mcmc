/*
 * WrappedMatrix.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.matrixAlgebra;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Variable;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

import static dr.math.matrixAlgebra.WrappedMatrix.Utils.makeString;

/**
 * @author Marc A. Suchard
 */

public interface WrappedMatrix extends ReadableMatrix, WritableVector, WritableMatrix {

    void set(final int i, final int j, final double x);

    double[] getBuffer();

    int getOffset();

    double[] data = null;

    String toString();

    abstract class Base implements WrappedMatrix {

        final public String toString() {
            return makeString(this);
        }
    }

    abstract class Abstract extends Base {
        final double[] buffer;
        final int offset;
        final int dimMajor;
        final int dimMinor;

        public Abstract(final double[] buffer, final int offset, final int dimMajor, final int dimMinor) {
            this.buffer = buffer;
            this.offset = offset;
            this.dimMajor = dimMajor;
            this.dimMinor = dimMinor;
        }

        final public double[] getBuffer() {
            return buffer;
        }

        final public int getOffset() {
            return offset;
        }

        final public int getMajorDim() {
            return dimMajor;
        }

        final public int getMinorDim() {
            return dimMinor;
        }

        public int getDim() {
            return getMajorDim() * getMinorDim();
        }
    }

    final class WrappedDenseMatrix extends Base {

        private final DenseMatrix64F matrix;

        public WrappedDenseMatrix(DenseMatrix64F matrix) {
            this.matrix = matrix;
        }

        @Override
        public void set(int i, double x) {
            matrix.data[i] = x;
        }

        @Override
        public double get(int i) {
            return matrix.data[i];
        }

        @Override
        public int getDim() {
            return matrix.getNumElements();
        }

        @Override
        public double get(int i, int j) {
            return matrix.unsafe_get(i, j);
        }

        @Override
        public int getMajorDim() {
            return matrix.getNumRows();
        }

        @Override
        public int getMinorDim() {
            return matrix.getNumCols();
        }

        @Override
        public void set(int i, int j, double x) {
            matrix.unsafe_set(i, j, x);
        }

        @Override
        public double[] getBuffer() {
            return matrix.data;
        }

        @Override
        public int getOffset() {
            return 0;
        }
    }

    final class ArrayOfArray extends Base {

        private final double[][] arrays;
        private final int dimMajor;
        private final int dimMinor;

        public ArrayOfArray(final double[][] arrays) {
            this.arrays = arrays;
            this.dimMajor = arrays.length;
            this.dimMinor = arrays[0].length;
        }

        @Override
        public double get(int i) {
            return arrays[i / dimMajor][i % dimMajor];
        }

        @Override
        public void set(int i, double x) {
            arrays[i / dimMajor][i % dimMajor] = x;
        }

        @Override
        public double get(int i, int j) {
            return arrays[i][j];
        }

        @Override
        public void set(int i, int j, double x) {
            arrays[i][j] = x;
        }

        @Override
        public int getMajorDim() {
            return dimMajor;
        }

        @Override
        public int getMinorDim() {
            return dimMinor;
        }

        @Override
        public int getDim() {
            return getMajorDim() * getMinorDim();
        }

        @Override
        public double[] getBuffer() {
            throw new RuntimeException("Not yet implemented");
        }

        public double[][] getArrays() {
            return arrays;
        }

        @Override
        public int getOffset() {
            return 0;
        }
    }

    final class Raw extends Abstract {

        public Raw(double[] buffer, int offset, int dimMajor, int dimMinor) {
            super(buffer, offset, dimMajor, dimMinor);
        }

        @Override
        final public double get(final int i) {
            return buffer[offset + i];
        }

        @Override
        final public void set(final int i, final double x) {
            buffer[offset + i] = x;
        }

        @Override
        final public double get(final int i, final int j) {
            return buffer[offset + i * dimMajor + j];
        }

        @Override
        final public void set(final int i, final int j, final double x) {
            buffer[offset + i * dimMajor + j] = x;
        }
    }

    final class Parameter extends Abstract {

        private final Variable<Double> variable;

        public Parameter(Variable<Double> variable, int offset, int dimMajor, int dimMinor) {
            super(null, offset, dimMajor, dimMinor);

            assert (variable.getSize() == dimMajor * dimMinor);

            this.variable = variable;
        }

        @Override
        final public double get(final int i, final int j) {
            return variable.getValue(offset + i * dimMajor + j);
        }

        @Override
        final public void set(final int i, final double x) {
            variable.setValue(offset + i, x);
        }

        @Override
        public void set(int i, int j, double x) {
            variable.setValue(offset + i * dimMajor + j, x);
        }

        @Override
        public double get(int i) {
            return variable.getValue(i);
        }
    }

    final class MatrixParameter extends Abstract {

        private final MatrixParameterInterface matrix;

        public MatrixParameter(MatrixParameterInterface matrix) {
            super(null, 0, matrix.getRowDimension(), matrix.getColumnDimension());

            this.matrix = matrix;
        }

        @Override
        public void set(int i, int j, double x) {
            matrix.setParameterValue(i, j, x);
        }

        @Override
        public double get(int i, int j) {
            return matrix.getParameterValue(i, j);
        }

        @Override
        public double get(int i) {
            return matrix.getParameterValue(i);
        }

        @Override
        public void set(int i, double x) {
            matrix.setParameterValue(i, x);
        }
    }

    final class Indexed extends Abstract {

        final private int[] indicesMajor;
        final private int[] indicesMinor;
        final int dimMajorFull;
        final int dimMinorFull;

        public Indexed(double[] buffer, int offset, int[] indicesMajor, int[] indicesMinor, int dimMajorFull, int dimMinorFull) {
            super(buffer, offset, indicesMajor.length, indicesMinor.length);
            this.indicesMajor = indicesMajor;
            this.indicesMinor = indicesMinor;
            this.dimMajorFull = dimMajorFull;
            this.dimMinorFull = dimMinorFull;
        }

        @Override
        final public double get(int i, int j) {
            return buffer[getIndex(i, j)];
        }

        @Override
        final public void set(int i, int j, double x) {
            buffer[getIndex(i, j)] = x;
        }

        @Override
        final public double get(int i) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        final public void set(int i, double x) {
            throw new RuntimeException("Not yet implemented");
        }

        private int getIndex(final int i, final int j) {
            return offset + indicesMajor[i] * dimMajorFull + indicesMinor[j];
        }
    }

    final class Utils {

        public static String makeString(ReadableMatrix matrix) {
            StringBuilder sb = new StringBuilder("[ ");
            if (matrix.getDim() > 0) {
                sb.append(matrix.get(0));
            }
            for (int i = 1; i < matrix.getDim(); ++i) {
                sb.append(", ").append(matrix.get(i));
            }
            sb.append(" ]");
            return sb.toString();
        }

        public static void gatherRowsAndColumns(final WrappedMatrix source, final WrappedMatrix destination,
                                                final int[] rowIndices, final int[] colIndices) {

            final int rowLength = rowIndices.length;
            final int colLength = colIndices.length;

            for (int i = 0; i < rowLength; ++i) {
                final int rowIndex = rowIndices[i];
                for (int j = 0; j < colLength; ++j) {
                    destination.set(i, j, source.get(rowIndex, colIndices[j]));
                }
            }
        }

        public static void scatterRowsAndColumns(final WrappedMatrix source, final WrappedMatrix destination,
                                                 final int[] rowIndices, final int[] colIndices) {
            scatterRowsAndColumns(source, destination, rowIndices, colIndices, true);
        }

        public static void scatterRowsAndColumns(final WrappedMatrix source, final WrappedMatrix destination,
                                                 final int[] rowIndices, final int[] colIndices, final boolean clear) {
            if (clear) {
                double[] data = destination.getBuffer();
                final int start = destination.getOffset();
                final int end = start + destination.getMajorDim() * destination.getMinorDim();
                Arrays.fill(data, start, end, 0.0);
            }

            final int rowLength = rowIndices.length;
            final int colLength = colIndices.length;

            for (int i = 0; i < rowLength; ++i) {
                final int rowIndex = rowIndices[i];
                for (int j = 0; j < colLength; ++j) {
                    destination.set(rowIndex, colIndices[i], source.get(i, j));
                }
            }
        }

        public static WrappedMatrix.Indexed wrapBlockDiagonalMatrix(double[] buffer, int offset, int startDim, int dimMat) {
            int[] inds = new int[dimMat];
            int dim = startDim;

            for (int i = 0; i < dimMat; i++) {
                inds[i] = dim;
                dim++;
            }

            return new WrappedMatrix.Indexed(buffer, offset, inds, inds, dimMat, dimMat);
        }

        public static void transferSymmetricBlockDiagonal(WrappedMatrix srcMat, WrappedMatrix destMat, int destOffset) {
            int srcDim = srcMat.getMajorDim();
            int destDim = destMat.getMajorDim();
            if (srcMat.getMinorDim() != srcDim || destMat.getMinorDim() != destDim) {
                throw new RuntimeException("Matrices must be square.");
            }

            int destI = destOffset;

            for (int i = 0; i < srcDim; i++) {

                destMat.set(destI, destI, srcMat.get(i, i));
                int destJ = destI + 1;

                for (int j = i + 1; j < srcDim; j++) {

                    double val = srcMat.get(i, j);
                    destMat.set(destI, destJ, val);
                    destMat.set(destJ, destI, val);
                    destJ++;
                }
                destI++;
            }
        }
    }

    final class WrappedUpperTriangularMatrix extends Abstract {
        // Class of square upper triangular matrices
        // Represented by an array of size dim*(dim+1)/2

        public WrappedUpperTriangularMatrix(final double[] buffer, int dim) {
            super(buffer, 0, dim, dim);
            assert buffer.length == dim * (dim + 1) / 2;
        }

        public WrappedUpperTriangularMatrix(int dim) {
            super(new double[dim * (dim + 1) / 2], 0, dim, dim);
        }

        @Override
        final public double get(final int i) {
            return buffer[offset + i];
        }

        @Override
        final public void set(final int i, final double x) {
            buffer[offset + i] = x;
        }

        @Override
        final public double get(final int i, final int j) {
            assert i <= j;
            return buffer[offset + pos(i, j)];
        }

        @Override
        final public void set(final int i, final int j, final double x) {
            assert i <= j;
            buffer[offset + pos(i, j)] = x;
        }

        private int pos(int i, int j) {
            return i * (2 * dimMajor - i + 1) / 2 + (j - i);
        }

        final public SymmetricMatrix transposedProduct() {
            SymmetricMatrix result = new SymmetricMatrix(dimMajor);
            for (int i = 0; i < dimMajor; i++) {
                for (int j = i; j < dimMajor; j++) {
                    double coef = 0;
                    for (int k = 0; k < i + 1; k++) {
                        coef += this.get(k, i) * this.get(k, j);
                    }
                    result.setSymmetric(i, j, coef);
                }
            }
            return result;
        }

        // Create the Cholesky of a correlation matrix (squared columns sum to 1)
        public static WrappedMatrix.WrappedUpperTriangularMatrix fillDiagonal(double[] x, int dim) {
            WrappedMatrix.WrappedStrictlyUpperTriangularMatrix V
                    = new WrappedMatrix.WrappedStrictlyUpperTriangularMatrix(x, dim);
            WrappedMatrix.WrappedUpperTriangularMatrix W
                    = new WrappedMatrix.WrappedUpperTriangularMatrix(dim);

            for (int j = 0; j < dim; j++) {
                double sum = 0;
                for (int i = 0; i < j; i++) {
                    double temp = V.get(i, j);
                    W.set(i, j, temp);
                    sum += temp * temp;
                }
                if (sum > 1.0) {
                    if (Math.abs(sum - 1.0) > 1E-6) {
                        throw new RuntimeException("Values are not consistent with the cholesky decomposition of " +
                                "a correlation matrix. Sum of squared values must be less than 1 (got " + sum + ")");
                    }
                    sum = 1.0;
                }
                W.set(j, j, Math.sqrt(1 - sum));
            }

            return W;
        }

        @Override
        final public int getDim() {
            return buffer.length;
        }
    }

    final class WrappedStrictlyUpperTriangularMatrix extends Abstract {
        // Class of square strictly upper triangular matrices
        // Represented by an array of size dim*(dim-1)/2

        private final double diagonalValue; // value on the diagonal (default to 0)

        public WrappedStrictlyUpperTriangularMatrix(final double[] buffer, int dim, double value) {
            super(buffer, 0, dim, dim);
            assert buffer.length == dim * (dim - 1) / 2;
            diagonalValue = value;
        }

        public WrappedStrictlyUpperTriangularMatrix(final double[] buffer, int dim) {
            super(buffer, 0, dim, dim);
            assert buffer.length == dim * (dim - 1) / 2;
            diagonalValue = 0.0;
        }

        public WrappedStrictlyUpperTriangularMatrix(int dim, double value) {
            super(new double[dim * (dim - 1) / 2], 0, dim, dim);
            diagonalValue = value;
        }

        public WrappedStrictlyUpperTriangularMatrix(int dim) {
            super(new double[dim * (dim - 1) / 2], 0, dim, dim);
            diagonalValue = 0.0;
        }

        @Override
        final public double get(final int i) {
            return buffer[offset + i];
        }

        @Override
        final public void set(final int i, final double x) {
            buffer[offset + i] = x;
        }

        @Override
        final public double get(final int i, final int j) {
            assert i <= j;
            if (i == j) return diagonalValue;
            return buffer[offset + pos(i, j)];
        }

        @Override
        final public void set(final int i, final int j, final double x) {
            assert i < j;
            buffer[offset + pos(i, j)] = x;
        }

        @Override
        final public int getDim() {
            return buffer.length;
        }

        private int pos(int i, int j) {
            return i * (2 * dimMajor - i - 1) / 2 + (j - i - 1);
        }

    }


}
