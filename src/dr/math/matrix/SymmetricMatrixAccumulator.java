/*
 * SymmetricMatrixAccumulator.java
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
 * In-place accumulator for a symmetric matrix G.
 *
 * Callers write entries via {@link #add}, which may be called for any (i,j) pair.
 * Implementations maintain only one triangle of storage; symmetry is implicit.
 *
 * @author Filippo Monti
 */
public interface SymmetricMatrixAccumulator {

    int dimension();

    /**
     * Add {@code value} to G[i][j] (and G[j][i]).
     * Indices may be given in any order; the implementation normalises to the stored triangle.
     */
    void add(int i, int j, double value);

    /** Returns x' G x. */
    double quadraticForm(double[] x);

    /** y += scale * G * x. */
    void addMatVec(double[] x, double scale, double[] y);

    /** target[i][j] += scale * G[i][j] for all i, j (both triangles filled in target). */
    void addScaledToDense(double scale, double[][] target);
}
