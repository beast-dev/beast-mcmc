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

import java.util.Arrays;

/**
 * Created by msuchard on 1/27/17.
 */

public interface WrappedMatrix {

    double get(final int i);

    void set(final int i, final double x);

    double get(final int i, final int j);

    void set(final int i, final int j, final double x);

    int getMajorDim();

    int getMinorDim();

    double[] getBuffer();

    int getOffset();

    double[] data = null;

    abstract class Abstract implements WrappedMatrix {
        final protected double[] buffer;
        final protected int offset;
        final protected int dimMajor;
        final protected int dimMinor;

        public Abstract(final double[] buffer, final int offset, final int dimMajor, final int dimMinor) {
            this.buffer = buffer;
            this.offset = offset;
            this.dimMajor = dimMajor;
            this.dimMinor = dimMinor;
        }

        final public double[] getBuffer() {
            return buffer;
        }

        final public int getOffset() { return offset; }

        final public int getMajorDim() {
            return dimMajor;
        }

        final public int getMinorDim() {
            return dimMinor;
        }

        final public String toString() {
            StringBuilder sb = new StringBuilder("[ ");
            // TODO
//            if (dim > 0) {
//                sb.append(get(0));
//            }
//            for (int i = 1; i < dim; ++i) {
//                sb.append(", ").append(get(i));
//            }
//            sb.append(" ]");

            return sb.toString();
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

    final class Indexed extends Abstract {

        final private int[] indicesMajor;
        final private int[] indicesMinor;

        public Indexed(double[] buffer, int offset, int[] indicesMajor, int[] indicesMinor, int dimMajor, int dimMinor) {
            super(buffer, offset, dimMajor, dimMinor);
            this.indicesMajor = indicesMajor;
            this.indicesMinor = indicesMinor;
        }

        @Override
        final public double get(int i, int j) { return buffer[getIndex(i, j)]; }

        @Override
        final public void set(int i, int j, double x) { buffer[getIndex(i, j)] = x; }

        @Override
        final public double get(int i) { throw new RuntimeException("Not yet implemented"); }

        @Override
        final public void set(int i, double x) { throw new RuntimeException("Not yet implemented"); }

        private final int getIndex(final int i, final int j) {
            return offset + indicesMajor[i] * dimMajor + indicesMinor[j];
        }
    }

    final class Utils {

        public static void gatherRowsAndColumns(final WrappedMatrix source, final WrappedMatrix destination,
                                                final int[] rowIndices, final int[] colIndices) {

            final int rowLength = rowIndices.length;
            final int colLength = colIndices.length;

            for (int i = 0; i < rowLength; ++i) {
                final int rowIndex = rowIndices[i];
                for (int j = 0; j < colLength; ++j) {
                    destination.set(i, j, source.get(rowIndex, colIndices[j]) );
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
    }
}
