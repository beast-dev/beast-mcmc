/*
 * CompoundSymmetricMatrix.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.model;

import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.CorrelationToCholesky;
import dr.util.Transform;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;


/**
 * @author Marc Suchard
 * @author Paul Bastide
 */
public class CompoundSymmetricMatrix extends AbstractTransformedCompoundMatrix {

    private final boolean asCorrelation;
    private final boolean isCholesky;

    public CompoundSymmetricMatrix(Parameter diagonals, Parameter offDiagonal, boolean asCorrelation, boolean isCholesky) {
        super(diagonals, offDiagonal, getTransformation(diagonals.getDimension(), isCholesky), true);
        assert asCorrelation || !isCholesky; // cholesky only allowed when used as correlation.
        this.asCorrelation = asCorrelation;
        this.isCholesky = isCholesky;
    }

    private static Transform.MultivariableTransform getTransformation(int dim, Boolean isCholesky) {
        return isCholesky ? new CorrelationToCholesky(dim) : null;
    }

    @Override
    public String toString() {
        return toStringCompoundParameter(getVechDimension(dim));
    }

    private static int getVechuDimension(int dim) {
        return dim * (dim - 1) / 2;
    }

    private static int getVechDimension(int dim) {
        return dim * (dim + 1) / 2;
    }

    @Override
    public double getParameterValue(int row, int col) {
        if (row != col) {
            if (asCorrelation) {
                return offDiagonalParameter.getParameterValue(getUpperTriangularIndex(row, col)) *
                        Math.sqrt(diagonalParameter.getParameterValue(row) * diagonalParameter.getParameterValue(col));
            }
            return offDiagonalParameter.getParameterValue(getUpperTriangularIndex(row, col));
        } else if (isStrictlyUpperTriangular) {
            return diagonalParameter.getParameterValue(row);
        }
        return diagonalParameter.getParameterValue(row) *
                offDiagonalParameter.getParameterValue(getUpperTriangularIndex(row, row));
    }

    @Override
    public double[][] getParameterAsMatrix() {
        final int I = dim;
        double[][] parameterAsMatrix = new double[I][I];
        for (int i = 0; i < I; i++) {
            parameterAsMatrix[i][i] = getParameterValue(i, i);
            for (int j = i + 1; j < I; j++) {
                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
            }
        }
        return parameterAsMatrix;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return true;
    }

    public boolean isCholesky() {
        return isCholesky;
    }

    public boolean asCorrelation() {
        return asCorrelation;
    }

    private double[][] getCorrelationMatrix() {
        SymmetricMatrix correlation
                = SymmetricMatrix.compoundCorrelationSymmetricMatrix(offDiagonalParameter.getParameterValues(), dim);
        if (!asCorrelation) {
            for (int i = 0; i < dim; i++) {
                for (int j = i + 1; j < dim; j++) {
                    correlation.setSymmetric(i, j,
                            correlation.component(i, j) / Math.sqrt(diagonalParameter.getParameterValue(i) * diagonalParameter.getParameterValue(j)));
                }
            }
        }
        return correlation.toComponents();
    }

    public double[] updateGradientOffDiagonal(double[] vecX) {

        assert vecX.length == dim * dim;

        double[] diagQ = diagonalParameter.getParameterValues();

        double[] vechuGradient = new double[getVechuDimension(dim)];

        int k = 0;
        for (int i = 0; i < dim - 1; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                vechuGradient[k] = 2.0 * vecX[i * dim + j] * Math.sqrt(diagQ[i] * diagQ[j]);
                ++k;
            }
        }

        return updateGradientCorrelation(vechuGradient);
    }

    public double[] updateGradientFullOffDiagonal(double[] gradient) {
        assert gradient.length == dim * dim;

        double[] diagQ = diagonalParameter.getParameterValues();

        double[] offDiagGradient = new double[gradient.length];

        int k = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                offDiagGradient[k] = gradient[i * dim + j] * Math.sqrt(diagQ[i] * diagQ[j]);
                ++k;
            }
        }

        return offDiagGradient;
    }

    public double[] updateGradientCorrelation(double[] gradient) {
        if (!isCholesky) {
            return gradient;
        } else {
            CorrelationToCholesky transform = new CorrelationToCholesky(dim);
            return transform.updateGradientInverseUnWeightedLogDensity(gradient,
                    ((TransformedMultivariateParameter) offDiagonalParameter).getParameterUntransformedValues(),
                    0, gradient.length);
        }
    }

    public double[] updateGradientDiagonal(double[] vecX) {

        assert vecX.length == dim * dim;

        double[] diagQ = diagonalParameter.getParameterValues();

        double[] vecC = flatten(getCorrelationMatrix());

        double[] diagGradient = new double[dim];

        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;
            for (int k = 0; k < dim; ++k) {
                sum += vecX[i * dim + k] * Math.sqrt(diagQ[k] / diagQ[i]) * vecC[i * dim + k];
            }

            diagGradient[i] = sum;
        }

        return diagGradient;
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }
}
