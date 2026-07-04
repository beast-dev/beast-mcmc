/*
 * OffsetBlockSymmetricAccumulator.java
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
 * A view into a parent {@link SymmetricMatrixAccumulator} that translates block-local indices
 * {@code (i, j)} to absolute parent indices {@code (rowOffset + i, colOffset + j)}.
 *
 * Used by {@code CovariateAugmentedBasisExpansion} to route interaction-covariate Gram blocks
 * (§9.2 of the refactor plan) without sub-allocating a separate accumulator.
 *
 * Only {@link #add} is supported; all other operations throw {@link UnsupportedOperationException}.
 * Pre-allocate instances in the outer class constructor to avoid hot-path allocation.
 *
 * @author Filippo Monti
 */
public final class OffsetBlockSymmetricAccumulator implements SymmetricMatrixAccumulator {

    private final SymmetricMatrixAccumulator parent;
    private final int rowOffset;
    private final int colOffset;

    public OffsetBlockSymmetricAccumulator(SymmetricMatrixAccumulator parent,
                                           int rowOffset,
                                           int colOffset) {
        this.parent    = parent;
        this.rowOffset = rowOffset;
        this.colOffset = colOffset;
    }

    /**
     * Translates (i, j) to (rowOffset+i, colOffset+j) and delegates to the parent.
     *
     * One call is sufficient for both diagonal and off-diagonal blocks:
     * {@link DenseSymmetricMatrixAccumulator#add} normalises any (r,c) pair to lower-triangle
     * storage via max(r,c)/min(r,c), so it correctly represents both G[r][c] and G[c][r]
     * regardless of whether r > c or r < c.
     */
    @Override
    public void add(int i, int j, double value) {
        parent.add(rowOffset + i, colOffset + j, value);
    }

    @Override
    public int dimension() {
        throw new UnsupportedOperationException("OffsetBlockSymmetricAccumulator: dimension() not supported");
    }

    @Override
    public double quadraticForm(double[] x) {
        throw new UnsupportedOperationException("OffsetBlockSymmetricAccumulator: quadraticForm() not supported");
    }

    @Override
    public void addMatVec(double[] x, double scale, double[] y) {
        throw new UnsupportedOperationException("OffsetBlockSymmetricAccumulator: addMatVec() not supported");
    }

    @Override
    public void addScaledToDense(double scale, double[][] target) {
        throw new UnsupportedOperationException("OffsetBlockSymmetricAccumulator: addScaledToDense() not supported");
    }
}
