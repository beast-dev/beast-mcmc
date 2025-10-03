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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.evolution.tree.Tree;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.matrix.SparseCompressedMatrix;
import no.uib.cipr.matrix.UpperTriangBandMatrix;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;

/**
 * @author Marc Suchard
 */

public abstract class AbstractBandedMultivariateDiffusionModel extends MultivariateDiffusionModel {

    final GaussianMarkovRandomField field;
    final int precisionDim;
//    private final int fieldDim;
    private final int blockDim;

    private SparseCompressedMatrix sparseMatrix;
    private SparseCompressedMatrix savedSparseMatrix;

    public AbstractBandedMultivariateDiffusionModel(MatrixParameterInterface block,
                                                    GaussianMarkovRandomField field) {
        this(block, field.getDimension(), field);
    }

    public AbstractBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates) {
        this(block, replicates, null);
    }

    private AbstractBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates,
                                                     GaussianMarkovRandomField field) {
        super();

        this.block = block;
        this.field = field;
        this.replicates = replicates;
        this.blockDim = block.getRowDimension();
        this.precisionDim = blockDim * replicates;

        variableChanged = true;

        addVariable(block);

        sparseMatrix = makeSparseMatrix();

        if (field != null) {
            addModel(field);
        }
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

    private SparseCompressedMatrix makeSparseUpperTriangularMatrix(UpperTriangBandMatrix fieldCholesky,
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

        SparseCompressedMatrix matrix = new SparseCompressedMatrix(
                rowStarts, columnIndices, values, precisionDim, precisionDim);

        return matrix;
    }

    public MatrixParameterInterface getPrecisionParameter() {
        // TODO Delegate
        if (field != null) {
            return block;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public int getDimension() {
        return precisionDim;
    }

    public double[][] getPrecisionmatrix() {
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

//    public SparseCompressedMatrix

    public double[] getPrecisionmatrixAsVector() {
        return(flatten(getPrecisionmatrix()));
    }

    private static final boolean CHECK_DETERMINANT = false;

    @Override
    public double getDeterminantPrecisionMatrix() {
        return Math.exp(getLogDeterminantPrecisionMatrix());
    }

    private SparseCompressedMatrix getPrecisionCholeskyDecomposition() {

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
        for (int i = 0; i < blockDim && isPD; ++i) {
            if (L[i][i] == 0.0) {
                isPD = false;
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
        SparseCompressedMatrix kronecker = makeSparseUpperTriangularMatrix(fieldCholeskyU, blockCholeskyU);

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

            double[][] check = getPrecisionmatrix();

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

    @Override
    public double getLogDeterminantPrecisionMatrix() {
        checkVariableChanged();

        double logDet;

        // TODO delegate
        if (field != null) {
            int effectiveDim = field.isImproper() ? replicates - 1 : replicates;
            logDet = blockDim * field.getLogDeterminant() + effectiveDim * logDeterminateBlockMatrix;
        } else {
            logDet = replicates * logDeterminateBlockMatrix;
        }

        if (CHECK_DETERMINANT) {

            double[][] precision = getPrecisionmatrix();
            RobustEigenDecomposition ed = new RobustEigenDecomposition(new DenseDoubleMatrix2D(precision));
            DoubleMatrix1D values = ed.getRealEigenvalues();
            double sum = 0.0;
            for (int i = 0; i < values.size(); ++i) {
                double v = values.get(i);
                if (Math.abs(v) > 1E-6) {
                    sum += Math.log(v);
                }
            }

            if (Math.abs(sum - logDet) > 1E-6) {
                throw new RuntimeException("Incorrect (pseudo-) determinant");
            }
        }

        return logDet;
    }

    protected void checkVariableChanged() {
        if (variableChanged) {
            calculatePrecisionInfo();
            variableChanged = false;
        }
    }

    protected void calculatePrecisionInfo() {
        blockMatrix = block.getParameterAsMatrix();
        logDeterminateBlockMatrix = MultivariateNormalDistribution.calculatePrecisionMatrixLogDeterminate(
                        blockMatrix);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        variableChanged = true;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        variableChanged = true;
    }

    protected void storeState() {
        logSavedDeterminateBlockMatrix = logDeterminateBlockMatrix;
        savedBlockMatrix = blockMatrix; // TODO This is WRONG! need to make copy
        storedVariableChanged = variableChanged;

        if (sparseMatrix != null) {
            if (savedSparseMatrix == null) {
                savedSparseMatrix = sparseMatrix.makeCopy();
            } else {
                sparseMatrix.copyTo(savedSparseMatrix);
            }
        }
    }

    protected void restoreState() {
        logDeterminateBlockMatrix = logSavedDeterminateBlockMatrix;
        blockMatrix = savedBlockMatrix; // TODO This is WRONG! need to sway pointers
        variableChanged = storedVariableChanged;

        SparseCompressedMatrix swap = sparseMatrix;
        sparseMatrix = savedSparseMatrix;
        savedSparseMatrix = swap;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    public String[] getTreeAttributeLabel() {
        return new String[]{PRECISION_TREE_ATTRIBUTE};
    }

    public String[] getAttributeForTree(Tree tree) {
        if (block != null) {
            return new String[]{block.toSymmetricString()};
        } else {
            return new String[]{"null"};
        }
    }

    final protected MatrixParameterInterface block;
    final int replicates;

    private double logDeterminateBlockMatrix;
    private double logSavedDeterminateBlockMatrix;

    private double[][] blockMatrix;
    private double[][] savedBlockMatrix;

    private double[][] precisionMatrix;
    private double[][] savedPrecisionMatrix;

    private boolean variableChanged;
    private boolean storedVariableChanged;

}

