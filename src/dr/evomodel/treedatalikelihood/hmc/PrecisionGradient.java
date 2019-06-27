/*
 * PrecisionGradient.java
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */
public class PrecisionGradient extends AbstractPrecisionGradient implements HessianWrtParameterProvider {

    public PrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider,
                             Likelihood likelihood,
                             MatrixParameterInterface parameter) {

        super(gradientWrtPrecisionProvider, likelihood, parameter);
    }

    NumericalHessianFromGradient hessianFromGradient;

    @Override
    double[] getGradientParameter(double[] gradient) {

        double[] gradientCorrelation = getGradientCorrelation(gradient);

        double[] gradientDiagonal = getGradientDiagonal(gradient);

        return mergeGradients(gradientDiagonal, gradientCorrelation);
    }

    @Override
    public Parameter getParameter() {
        return compoundSymmetricMatrix.getUntransformedCompoundParameter();
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    private double[] mergeGradients(double[] gradientDiagonal, double[] gradientCorrelation) {
        double[] gradient = new double[gradientDiagonal.length + gradientCorrelation.length];
        System.arraycopy(gradientDiagonal, 0, gradient, 0, gradientDiagonal.length);
        System.arraycopy(gradientCorrelation, 0, gradient, gradientDiagonal.length, gradientCorrelation.length);
        return gradient;
    }


    MultivariateFunction getNumeric() {
        return null;
    }

    @Override
    String checkNumeric(double[] analytic) {
        return "";
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        return new double[0];
    }

    @Override
    public double[][] getHessianLogDensity() {
        return new double[0][];
    }
}
