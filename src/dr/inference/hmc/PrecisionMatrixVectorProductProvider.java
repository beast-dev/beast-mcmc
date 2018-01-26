/*
 * GradientWrtParameterProvider.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public interface PrecisionMatrixVectorProductProvider {

    double[] getProduct(Parameter vector);

    double[] getMassVector();

    double getTimeScale();

    class Generic implements PrecisionMatrixVectorProductProvider {

        private final MatrixParameterInterface matrix;

        public Generic(MatrixParameterInterface matrix) {
            this.matrix = matrix;
        }

        @Override
        public double[] getProduct(Parameter vector) {

            final int nRows = matrix.getRowDimension();
            final int nCols = matrix.getColumnDimension();

            assert (vector.getDimension() == nCols);

            double[] result = new double[nRows];

            for (int row = 0; row < nRows; ++row) {
                double sum = 0.0;
                for (int col = 0; col < nCols; ++col) {
                    sum += matrix.getParameterValue(row, col) * vector.getParameterValue(col);
                }
                result[row] = sum;
            }

            return result;
        }

        @Override
        public double[] getMassVector() {
            final int dim = Math.min(matrix.getRowDimension(), matrix.getColumnDimension());

            double[] mass = new double[dim];
            for (int i = 0; i < dim; ++i) {
                mass[i] = matrix.getParameterValue(i, i);
            }

            return mass;
        }

        @Override
        public double getTimeScale() {
            final int dim = Math.min(matrix.getRowDimension(), matrix.getColumnDimension());

            double max = Double.MIN_VALUE;
            for (int i = 0; i < dim; ++i) {
                max = Math.max(max, matrix.getParameterValue(i,i));
            }

            return Math.sqrt(max);
        }
    }
}