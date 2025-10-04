/*
 * MultivariateDiffusionModel.java
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

package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.matrix.SparseCompressedMatrix;
import dr.matrix.SparseSquareUpperTriangular;
import no.uib.cipr.matrix.*;

import static dr.evomodel.continuous.MultivariateDiffusionModel.PreferredSimulationSpace.PRECISION;

/**
 * @author Marc Suchard
 */

public class SparseBandedMultivariateDiffusionModel extends AbstractBandedMultivariateDiffusionModel {

    private SparseCompressedMatrix sparseMatrix;
    private SparseCompressedMatrix savedSparseMatrix;

    public SparseBandedMultivariateDiffusionModel(MatrixParameterInterface block,
                                                  GaussianMarkovRandomField field) {
        this(block, field.getDimension(), field);
    }

    public SparseBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates) {
        this(block, replicates, null);
    }

    private SparseBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates,
                                                   GaussianMarkovRandomField field) {
        super(block, replicates, field);
        sparseMatrix = makeSparseMatrix();
    }

    private SparseCompressedMatrix makeSparseMatrix() {

        if (field == null) {
            throw new RuntimeException("Not yet implemented");
        }

        int nonZeroEntryCount = field.getNonZeroEntryCount() * blockDim * blockDim;

        int[] columnIndices = new int[nonZeroEntryCount];
        double[] values = new double[nonZeroEntryCount];

        int[] rowStarts = new int[precisionDim + 1];
        for (int i = 0; i < replicates; ++i) { // outer-loop over rows of field
            for (int k = 0; k < blockDim; ++k) { // inner-loop over rows of block

                final int finalK = k;
                final int startsIndex = i * blockDim + k + 1;
                rowStarts[startsIndex] = rowStarts[startsIndex - 1];

                field.mapOverNonZeroEntriesInRow((fieldI, fieldJ, fieldValue) -> { // outer-loop over columns of field
                    for (int l = 0; l < blockDim; ++l) { // inner-loop over columns of block
                        final int index = rowStarts[startsIndex];
                        columnIndices[index] = fieldJ * blockDim + l;
                        values[index] = fieldValue * block.getParameterValue(finalK, l);
                        ++rowStarts[startsIndex];
                    }
                }, i);
            }
        }

        return new SparseCompressedMatrix(rowStarts, columnIndices, values, precisionDim, precisionDim);
    }

    public PreferredSimulationSpace getPreferredSimulationSpace() {
        return PRECISION;
    }

    private SparseSquareUpperTriangular makeSparseUpperTriangularMatrix(
            UpperTriangBandMatrix fieldCholesky,
            UpperTriangDenseMatrix blockCholesky) {

        if (field == null) {
            throw new RuntimeException("Not yet implemented");
        }

        int bandWidth = fieldCholesky.numSuperDiagonals();

        int nonZeroEntryCount = (field.getNonZeroEntryCount() * blockDim * blockDim + precisionDim) / 2;

        int[] columnIndices = new int[nonZeroEntryCount];
        double[] values = new double[nonZeroEntryCount];
        int[] rowStarts = new int[precisionDim + 1];

        for (int i = 0; i < replicates; ++i) { // outer-loop over rows of field
            for (int k = 0; k < blockDim; ++k) { // inner-loop over rows of block

                final int startsIndex = i * blockDim + k + 1;
                rowStarts[startsIndex] = rowStarts[startsIndex - 1];

                for (int j = i; j < Math.min(i + bandWidth + 1, replicates); ++j) { // outer-loop over columns of field
                    for (int l = k; l < blockDim; ++l) { // inner-loop over columns of block
                        int index = rowStarts[startsIndex];
                        columnIndices[index] = j * blockDim + l;
                        values[index] = fieldCholesky.get(i, j) * blockCholesky.get(k, l);
                        ++rowStarts[startsIndex];
                    }
                }
            }
        }

        return new SparseSquareUpperTriangular(rowStarts, columnIndices, values, precisionDim);
    }

    @Override
    public double[][] getPrecisionMatrix() {
        checkVariableChanged();
        precisionMatrix = new double[precisionDim][precisionDim];

        // TODO delegate this work
        if (field != null) {
            field.mapOverAllNonZeroEntries((fieldI, fieldJ, fieldValue) ->
                    fillSquareBlock(
                            precisionMatrix, fieldI * blockDim, fieldJ * blockDim,
                            block, fieldValue));
        } else {
            for (int r = 0; r < replicates; ++r) {
                int offset = r * blockDim;
                fillSquareBlock(precisionMatrix, offset, offset, block, 1.0);
            }
        }
        return precisionMatrix;
    }

    private void fillSquareBlock(double[][] destination, int rowOffset, int columnOffset,
                                 MatrixParameterInterface source,
                                 double scalar) {
        assert source.getRowDimension() == source.getColumnDimension();

        final int dim = source.getRowDimension();
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                destination[rowOffset + i][columnOffset + j] = scalar * source.getParameterValue(i, j);
            }
        }
    }

    public SparseCompressedMatrix getSparsePrecisionMatrix() {
        checkVariableChanged();
//        getPrecisionCholeskyDecomposition();
        return makeSparseMatrix(); // TODO cache!!!
    }

    @SuppressWarnings("unused")
    public SparseSquareUpperTriangular getPrecisionCholeskyDecomposition() {

        checkVariableChanged();
        UpperTriangBandMatrix fieldCholeskyU = field.getCholeskyDecomposition();

        // TODO returns different decomposition when matrix is not PD
//        UpperSymmDenseMatrix prec = new UpperSymmDenseMatrix(blockDim);
//        for (int i = 0; i < blockDim; ++i) {
//            for (int j = i; j < blockDim; ++j) {
//                prec.set(i, j, blockMatrix[i][j]);
//            }
//        }
//        DenseCholesky chol = DenseCholesky.factorize(prec);
//        UpperTriangDenseMatrix blockCholeskyU = chol.getU();

        double[][] L;
        try {
            CholeskyDecomposition cd = new CholeskyDecomposition(blockMatrix);
            L = cd.getL();
        } catch (IllegalDimension e) {
            L = new double[0][];
        }

        boolean isPD = true;
        for (int i = 0; i < blockDim; ++i) {
            if (L[i][i] == 0.0) {
                isPD = false;
                break;
            }
        }

        UpperTriangDenseMatrix blockCholeskyU = new UpperTriangDenseMatrix(blockDim);
        for (int i = 0; i < blockDim; ++i) {
            for (int j = i; j < blockDim; ++j) {
                blockCholeskyU.set(i, j, L[j][i]);
            }
        }

        if (TEST_CHOLESKY && isPD) {
            double[][] test = GaussianMarkovRandomField.testCholeskyUpper(blockCholeskyU, blockDim);
            for (int i = 0; i < blockDim; ++i) {
                for (int j = 0; j < blockDim; ++j) {
                    if (Math.abs(blockMatrix[i][j] - test[i][j]) > 1E-8) {
                        throw new RuntimeException("Bad Cholesky decomposition in block matrix");
                    }
                }
            }
        }

        // Need to Kronecker product
        SparseSquareUpperTriangular kronecker = makeSparseUpperTriangularMatrix(fieldCholeskyU, blockCholeskyU);

        if (TEST_CHOLESKY && isPD) {

            double[][] dense = kronecker.makeDense();
            double[][] UtU = new double[precisionDim][precisionDim];

            for (int i = 0; i < precisionDim; ++i) {
                for (int j = 0; j < precisionDim; ++j) {
                    for (int k = 0; k < precisionDim; ++k) {
                        double Ut = dense[k][i];
                        double U = dense[k][j];
                        UtU[i][j] += Ut * U;
                    }
                }
            }

            double[][] check = getPrecisionMatrix();

            for (int i = 0; i < precisionDim; ++i) {
                for (int j = 0; j < precisionDim; ++j) {
                    if (Math.abs(check[i][j] - UtU[i][j]) > 1E-8) {
                        throw new RuntimeException("Bad Cholesky decomposition in whole precision");
                    }
                }
            }
        }

        return kronecker;
    }

    private static final boolean TEST_CHOLESKY = false;

    protected void storeState() {
        super.storeState();
        if (sparseMatrix != null) {
            if (savedSparseMatrix == null) {
                savedSparseMatrix = sparseMatrix.makeCopy();
            } else {
                sparseMatrix.copyTo(savedSparseMatrix);
            }
        }
    }

    protected void restoreState() {
        super.restoreState();
        SparseCompressedMatrix swap = sparseMatrix;
        sparseMatrix = savedSparseMatrix;
        savedSparseMatrix = swap;
    }

    private double[][] precisionMatrix;
}

