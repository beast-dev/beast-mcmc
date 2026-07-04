/*
 * DenseSymmetricMatrixAccumulator.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.matrix;

/**
 * Dense lower-triangle implementation of {@link SymmetricMatrixAccumulator}.
 *
 * Storage: {@code values[i][j]} for {@code i >= j} holds G[i][j] = G[j][i].
 * Off-diagonal entries appear only once; {@link #quadraticForm} and
 * {@link #addMatVec} multiply them by 2 to account for symmetry.
 *
 * @author Filippo Monti
 */
public final class DenseSymmetricMatrixAccumulator implements SymmetricMatrixAccumulator {

    private final double[][] values; // lower-triangle: values[i][j] for i >= j

    public DenseSymmetricMatrixAccumulator(int dim) {
        this.values = new double[dim][dim];
    }

    @Override
    public int dimension() {
        return values.length;
    }

    @Override
    public void add(int i, int j, double value) {
        if (i >= j) values[i][j] += value;
        else        values[j][i] += value;
    }

    @Override
    public double quadraticForm(double[] x) {
        int n = values.length;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += values[i][i] * x[i] * x[i];
            for (int j = 0; j < i; j++) {
                sum += 2.0 * values[i][j] * x[i] * x[j];
            }
        }
        return sum;
    }

    @Override
    public void addMatVec(double[] x, double scale, double[] y) {
        int n = values.length;
        for (int i = 0; i < n; i++) {
            y[i] += scale * values[i][i] * x[i];
            for (int j = 0; j < i; j++) {
                y[i] += scale * values[i][j] * x[j];
                y[j] += scale * values[i][j] * x[i];
            }
        }
    }

    @Override
    public void addScaledToDense(double scale, double[][] target) {
        int n = values.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double v = scale * values[i][j];
                target[i][j] += v;
                if (i != j) target[j][i] += v;
            }
        }
    }
}
