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
import dr.matrix.SparseCompressedMatrix;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;

/**
 * @author Marc Suchard
 */

public class DenseBandedMultivariateDiffusionModel extends AbstractBandedMultivariateDiffusionModel {

    public DenseBandedMultivariateDiffusionModel(MatrixParameterInterface block,
                                                 GaussianMarkovRandomField field) {
        super(block, field.getDimension(), field);
    }

    public DenseBandedMultivariateDiffusionModel(MatrixParameterInterface block, int replicates) {
        super(block, replicates, null);
    }

    public SparseCompressedMatrix getSparsePrecisionMatrix() {
        throw new RuntimeException("Unsupported");
    }

    public double[][] getPrecisionMatrix() {
        checkVariableChanged();
        double[][] precisionMatrix = new double[precisionDim][precisionDim];

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

    public double[] getPrecisionMatrixAsVector() {
        return(flatten(getPrecisionMatrix()));
    }
}

