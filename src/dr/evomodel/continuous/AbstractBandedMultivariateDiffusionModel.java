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
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.matrix.SparseCompressedMatrix;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;

/**
 * @author Marc Suchard
 */

public abstract class AbstractBandedMultivariateDiffusionModel extends MultivariateDiffusionModel {

    final GaussianMarkovRandomField field;
    final int precisionDim;
//    private final int fieldDim;
    final int blockDim;

    AbstractBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates,
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

    public abstract double[][] getPrecisionmatrix();

    public abstract SparseCompressedMatrix getSparsePrecisionMatrix();

    public double[] getPrecisionmatrixAsVector() {
        return(flatten(getPrecisionmatrix()));
    }

    private static final boolean CHECK_DETERMINANT = false;

    @Override
    public double getDeterminantPrecisionMatrix() {
        return Math.exp(getLogDeterminantPrecisionMatrix());
    }

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
    }

    protected void restoreState() {
        logDeterminateBlockMatrix = logSavedDeterminateBlockMatrix;
        blockMatrix = savedBlockMatrix; // TODO This is WRONG! need to sway pointers
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

    double[][] blockMatrix;
    private double[][] savedBlockMatrix;

    private boolean variableChanged;
    private boolean storedVariableChanged;

}

