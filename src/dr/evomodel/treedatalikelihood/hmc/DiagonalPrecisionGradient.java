/*
 * DiagonalPrecisionGradient.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class DiagonalPrecisionGradient extends AbstractPrecisionGradient implements HessianWrtParameterProvider {

    public DiagonalPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider,
                                     Likelihood likelihood,
                                     MatrixParameterInterface parameter) {

        super(gradientWrtPrecisionProvider, likelihood, parameter, 0.0, Double.POSITIVE_INFINITY);
    }

    @Override
    public int getDimension() {
        return getDimensionDiagonal();
    }

    @Override
    public Parameter getParameter() {
        return compoundSymmetricMatrix.getDiagonalParameter();
    }

    @Override
    double[] getGradientParameter(double[] gradient) {
        return getGradientDiagonal(gradient);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        NumericalHessianFromGradient hessianFromGradient = new NumericalHessianFromGradient(this);
        return hessianFromGradient.getDiagonalHessianLogDensity();
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getReport() {
        return "diagonalPrecisionGradient." + compoundSymmetricMatrix.getParameterName() + "\n" +
                super.getReport();
    }
}
