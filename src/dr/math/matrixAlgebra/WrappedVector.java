/*
 * WrappedVector.java
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

import dr.inference.model.Parameter;

/**
 * Created by msuchard on 1/27/17.
 */

public interface WrappedVector {

    double get(final int i);

    void set(final int i, final double x);

    int getDim();

    double[] getBuffer();

    int getOffset();

    abstract class Abstract implements  WrappedVector {
        final protected double[] buffer;
        final protected int offset;
        final protected int dim;

        public Abstract(final double[] buffer, final int offset, final int dim) {
            this.buffer = buffer;
            this.offset = offset;
            this.dim = dim;
        }

        final public double[] getBuffer() {
            return buffer;
        }

        final public int getOffset() { return offset; }

        final public int getDim() {
            return dim;
        }

        final public String toString() {
            StringBuilder sb = new StringBuilder("[ ");
            if (dim > 0) {
                sb.append(get(0));
            }
            for (int i = 1; i < dim; ++i) {
                sb.append(", ").append(get(i));
            }
            sb.append(" ]");

            return sb.toString();
        }
    }

    final class Raw extends Abstract {

        public Raw(double[] buffer, int offset, int dim) {
            super(buffer, offset, dim);
        }

        @Override
        final public double get(final int i) {
            return buffer[offset + i];
        }

        @Override
        final public void set(final int i, final double x) {
            buffer[offset + i] = x;
        }
    }

    final class Indexed extends Abstract {

        final private int[] indices;

        public Indexed(double[] buffer, int offset, int[] indices, int dim) {
            super(buffer, offset, dim);
            this.indices = indices;
        }

        @Override
        final public double get(int i) { return buffer[offset + indices[i]]; }

        @Override
        final public void set(int i, double x) { buffer[offset + indices[i]] = x; }
    }
}
