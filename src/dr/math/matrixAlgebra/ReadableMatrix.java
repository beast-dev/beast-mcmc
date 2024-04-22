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

import static dr.math.matrixAlgebra.WrappedMatrix.Utils.makeString;

/**
 * @author Marc A. Suchard
 */

public interface ReadableMatrix extends ReadableVector {

    double get(final int i, final int j);

    int getMajorDim();

    int getMinorDim();

    class Utils {

        public static double[][] toMatrixArray(ReadableMatrix Y) { // defeats the purpose of wrapping but useful helper function sometimes
            double[][] X = new double[Y.getMajorDim()][Y.getMinorDim()];
            for (int i = 0; i < Y.getMajorDim(); i++) {
                for (int j = 0; j < Y.getMinorDim(); j++) {
                    X[i][j] = Y.get(i, j);
                }
            }
            return X;
        }

        public static double[] toArray(ReadableMatrix matrix) {
            double[] array = new double[matrix.getDim()];
            int offset = 0;
            for (int i = 0; i < matrix.getMajorDim(); ++i) {
                for (int j = 0; j < matrix.getMinorDim(); ++j) {
                    array[offset] = matrix.get(i, j);
                    ++offset;
                }
            }
            return array;
        }

        public static WrappedVector product(ReadableMatrix matrix, ReadableVector vector) {

            final int majorDim = matrix.getMajorDim();
            final int minorDim = matrix.getMinorDim();

            assert (vector.getDim() == minorDim);

            final double[] result = new double[majorDim];

            for (int row = 0; row < majorDim; ++row) {

                double sum = 0.0;
                for (int col = 0; col < minorDim; ++col) {
                    sum += matrix.get(row, col) * vector.get(col);
                }
                result[row] = sum;
            }

            return new WrappedVector.Raw(result);
        }

        public static ReadableMatrix transposeProxy(final ReadableMatrix matrix) {
            return new ReadableMatrix() {
                @Override
                public double get(int i, int j) {
                    return matrix.get(j, i);
                }

                @Override
                public int getMajorDim() {
                    return matrix.getMinorDim();
                }

                @Override
                public int getMinorDim() {
                    return matrix.getMajorDim();
                }

                @Override
                public double get(int i) {
                    throw new RuntimeException("Not yet implemented");
                }

                @Override
                public int getDim() {
                    return matrix.getDim();
                }
            };
        }

        public static ReadableMatrix productProxy(final ReadableMatrix lhs, final ReadableMatrix rhs) {
            final int majorDim = lhs.getMajorDim();
            final int innerDim = lhs.getMinorDim();
            final int minorDim = rhs.getMinorDim();

            assert (innerDim == rhs.getMajorDim());

            return new ReadableMatrix() {

                @Override
                public double get(int i, int j) {

                    double sum = 0.0;
                    for (int k = 0; k < innerDim; ++k) {
                        sum += lhs.get(i, k) * rhs.get(k,j);
                    }

                    return sum;
                }

                @Override
                public int getMajorDim() {
                    return majorDim;
                }

                @Override
                public int getMinorDim() {
                    return minorDim;
                }

                @Override
                public double get(int i) {
                    return get(i / minorDim, i % minorDim);
                }

                @Override
                public int getDim() {
                    return majorDim * minorDim;
                }

                @Override
                public String toString() {
                    return makeString(this);
                }
            };
        }

        public static WrappedVector product(ReadableMatrix lhs, ReadableMatrix rhs) {

            final int majorDim = lhs.getMajorDim();
            final int innerDim = lhs.getMinorDim();
            final int minorDim = rhs.getMinorDim();

            assert (innerDim == rhs.getMajorDim());

            final double[] result = new double[majorDim * minorDim];

            for (int row = 0; row < majorDim; ++row) {
                for (int col = 0; col < minorDim; ++col) {
                    double sum = 0.0;
                    for (int inner = 0; inner < innerDim; ++inner) {
                        sum += lhs.get(row, inner) * rhs.get(inner, col);
                    }
                    result[row * minorDim + col] = sum;
                }
            }

            return new WrappedVector.Raw(result);
        }
    }
}
