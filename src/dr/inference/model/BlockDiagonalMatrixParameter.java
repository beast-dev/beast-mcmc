/*
 * BlockDiagonalMatrixParameter.java
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

package dr.inference.model;

import dr.math.distributions.GaussianMarkovRandomField;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class BlockDiagonalMatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    // feed with a DiagonalCorrelationMatrix

    private final MatrixParameterInterface block;
    private final GaussianMarkovRandomField field;

    private final int replicates;
    private final int blockDim;
    private final int matrixDim;

    public BlockDiagonalMatrixParameter(String name,
                                        MatrixParameterInterface block,
                                        GaussianMarkovRandomField field) {
        this(name, block, field.getDimension() / block.getColumnDimension(), field);
    }

    public BlockDiagonalMatrixParameter(String name,
                                        MatrixParameterInterface block, int replicates) {
        this(name, block, replicates, null);
    }

    public BlockDiagonalMatrixParameter(String name,
                                        MatrixParameterInterface block, int replicates,
                                        GaussianMarkovRandomField field) {
        super(name, new Parameter[] { block });

        this.block = block;
        this.field = field;
        this.replicates = replicates;
        this.blockDim = block.getRowDimension();
        this.matrixDim = replicates * blockDim;

        if (blockDim != block.getColumnDimension()) {
            throw new IllegalArgumentException("Must specify a square block");
        }

        if (matrixDim % blockDim != 0) {
            throw new IllegalArgumentException("Matrix dimension is not divisible by block dimension");
        }
    }

    public double getParameterValue(int row, int col) {

        int blockRow = row / blockDim;
        int innerRow = row % blockDim;

        int blockCol = col / blockDim;
        int innerCol = col % blockDim;

        if (blockRow != blockCol) {
            return 0.0;
        } else {
            return block.getParameterValue(innerRow, innerCol);
        }
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {

        double[][] matrix = new double[getRowDimension()][getColumnDimension()];
        for (int r = 0; r < replicates; ++r) {
            int offset = r * blockDim;
            for (int i = 0; i < blockDim; ++i) {
                for (int j = 0; j < blockDim; ++j) {
                    matrix[offset + i][offset + j] = block.getParameterValue(i, j);
                }
            }
        }
        return matrix;
    }

    @Override
    public int getUniqueParameterCount() {
        return uniqueParameters.size();
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return uniqueParameters.get(index);
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public boolean isConstrainedSymmetric() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getRowDimension() {
        return matrixDim;
    }

    @Override
    public int getColumnDimension() {
        return matrixDim;
    }
}
