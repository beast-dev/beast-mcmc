/*
 * PermutationIndices.java
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

package dr.math.matrixAlgebra.missingData;

import org.ejml.data.DenseMatrix64F;

public class PermutationIndices {

    private final DenseMatrix64F matrix;
    private final int dim;

    private int zeroCount;
    private int nonZeroFiniteCount;
    private int infiniteCount;

    private int[] nonZeroFiniteIndices;
    private int[] zeroIndices;
    private int[] infiniteIndices;

    public PermutationIndices(DenseMatrix64F matrix) {

        this.matrix = matrix;
        dim = matrix.getNumCols();
        assert (dim == matrix.getNumRows());

        for (int i = 0; i < dim; ++i) {
            double diagonal = matrix.get(i, i);
            if (Double.isInfinite(diagonal)) {
                ++infiniteCount;
            } else if (diagonal == 0.0) {
                ++zeroCount;
            } else {
                ++nonZeroFiniteCount;
            }
        }
    }

    @SuppressWarnings("unused")
    public int getNumberOfZeroDiagonals() {
        return zeroCount;
    }

    public int getNumberOfNonZeroFiniteDiagonals() {
        return nonZeroFiniteCount;
    }

    @SuppressWarnings("unused")
    public int getNumberOfInfiniteDiagonals() {
        return infiniteCount;
    }

    public int[] getNonZeroFiniteIndices() {
        if (nonZeroFiniteIndices == null) {
            makeIndices();
        }
        return nonZeroFiniteIndices;
    }

    public int[] getZeroIndices() {
        if (zeroIndices == null) {
            makeIndices();
        }
        return zeroIndices;
    }

    @SuppressWarnings("unused")
    public int[] getInfiniteIndices() {
        if (infiniteIndices == null) {
            makeIndices();
        }
        return infiniteIndices;
    }

    private void makeIndices() {
        nonZeroFiniteIndices = new int[nonZeroFiniteCount];
        zeroIndices = new int[zeroCount];
        infiniteIndices = new int[infiniteCount];

        int zeroIndex = 0;
        int nonZeroFiniteIndex = 0;
        int infiniteIndex = 0;

        for (int i = 0; i < dim; ++i) {
            double diagonal = matrix.get(i, i);
            if (Double.isInfinite(diagonal)) {
                infiniteIndices[infiniteIndex] = i;
                ++infiniteIndex;
            } else if (diagonal == 0.0) {
                zeroIndices[zeroIndex] = i;
                ++zeroIndex;
            } else {
                nonZeroFiniteIndices[nonZeroFiniteIndex] = i;
                ++nonZeroFiniteIndex;
            }
        }
    }
}
