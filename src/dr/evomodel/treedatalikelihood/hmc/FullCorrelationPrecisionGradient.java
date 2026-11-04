/*
 * FullCorrelationPrecisionGradient.java
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.model.*;
import dr.util.MatrixInnerProductTransform;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class FullCorrelationPrecisionGradient extends CorrelationPrecisionGradient {

    private final MatrixParameterInterface decomposedMatrix;
    private static final RuntimeException PARAMETER_EXCEPTION = new RuntimeException("off-diagonal parameter must be a mask of a inner product transform.");

    public FullCorrelationPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider, Likelihood likelihood, MatrixParameterInterface parameter) {
        super(gradientWrtPrecisionProvider, likelihood, parameter);

        // TODO: this is super messy but get's the job done (for now). Maybe create a class just to wrap these different transforms/masks
        Parameter correlationParameter = compoundSymmetricMatrix.getOffDiagonalParameter();
        if (correlationParameter instanceof MaskedParameter) {
            MaskedParameter maskedCorrelation = (MaskedParameter) correlationParameter;
            if (maskedCorrelation.getUnmaskedParameter() instanceof TransformedMultivariateParameter) {
                TransformedMultivariateParameter transformedParameter =
                        (TransformedMultivariateParameter) maskedCorrelation.getUnmaskedParameter();
                if (transformedParameter.getTransform() instanceof MatrixInnerProductTransform) { //TODO: chain rul in transforms
                    decomposedMatrix = (MatrixParameterInterface) transformedParameter.getUntransformedParameter();
                } else {
                    throw PARAMETER_EXCEPTION;
                }
            } else {
                throw PARAMETER_EXCEPTION;
            }
        } else {
            throw PARAMETER_EXCEPTION;
        }

    }


    @Override
    public int getDimension() {
        return getDimensionFull();
    }

    @Override
    public Parameter getParameter() {
        return decomposedMatrix;
    }

    @Override
    double[] getGradientParameter(double[] gradient) {
        int dim = decomposedMatrix.getRowDimension();

        double[] correlationGradient = compoundSymmetricMatrix.updateGradientFullOffDiagonal(gradient);
        double[] decomposedGradient = new double[gradient.length];

        DenseMatrix64F corGradMat = DenseMatrix64F.wrap(dim, dim, correlationGradient);
        DenseMatrix64F decompMat = DenseMatrix64F.wrap(dim, dim, decomposedMatrix.getParameterValues());
        DenseMatrix64F decompGradMat = DenseMatrix64F.wrap(dim, dim, decomposedGradient);

        CommonOps.multTransA(corGradMat, decompMat, decompGradMat);

        CommonOps.scale(2.0, decompGradMat);

        return decomposedGradient;
    }


}
