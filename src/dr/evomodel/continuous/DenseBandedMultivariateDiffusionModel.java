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
import dr.inference.model.*;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.RobustEigenDecomposition;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;

/**
 * @author Marc Suchard
 */

public class DenseBandedMultivariateDiffusionModel extends MultivariateDiffusionModel {
//        AbstractModel implements TreeAttributeProvider {

    // TODO create interface for `MultivariateDIffusion`

    private final GaussianMarkovRandomField field;
    private final int precisionDim;
    private final int blockDim;

    public DenseBandedMultivariateDiffusionModel(MatrixParameterInterface block,
                                                 GaussianMarkovRandomField field) {
        this(block, field.getDimension(), field);
    }

    public DenseBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates) {
        this(block, replicates, null);
    }

    private DenseBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates,
                                                  GaussianMarkovRandomField field) {
        super();

        this.block = block;
        this.field = field;
        this.replicates = replicates;
        this.blockDim = block.getRowDimension();
        this.precisionDim = blockDim * replicates;

        variableChanged = true;

        addVariable(block);

        if (field != null) {
            addModel(field);
        }
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
            double[][] gmrf = field.getScaleMatrix();

            assert replicates == gmrf.length;

            // TODO do O(replicates) loops instead of O(replicates^2) here
            for (int k = 0; k < replicates; ++k) {
                for (int l = 0; l < replicates; ++l) {
                    double scalar = gmrf[k][l];
                    if (scalar != 0.0) {
                        fillSquareBlock(
                                precisionMatrix, k * blockDim, l * blockDim,
                                block, scalar);
                    }
                }
            }
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

    public double[] getPrecisionmatrixAsVector() {
        return(flatten(getPrecisionmatrix()));
    }

    private static final boolean CHECK_DETERMINANT = true;

    public double getDeterminantPrecisionMatrix() {
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

        return Math.exp(logDet);
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
        savedBlockMatrix = blockMatrix;
        storedVariableChanged = variableChanged;
    }

    protected void restoreState() {
        logDeterminateBlockMatrix = logSavedDeterminateBlockMatrix;
        blockMatrix = savedBlockMatrix;
        variableChanged = storedVariableChanged;
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

