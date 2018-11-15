/*
 * CorrelationPrecisionGradient.java
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

import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class CorrelationPrecisionGradient extends AbstractPrecisionGradient {

    public CorrelationPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider,
                                        Likelihood likelihood,
                                        MatrixParameterInterface parameter) {

        super(gradientWrtPrecisionProvider, likelihood, parameter);
    }

    @Override
    public int getDimension() {
        return getDimensionCorrelation();
    }

    @Override
    public Parameter getParameter() {
        return compoundSymmetricMatrix.getUntransformedOffDiagonalParameter();
    }

    MultivariateFunction getNumeric() {

        return new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    compoundSymmetricMatrix.getOffDiagonalParameter().setParameterValue(i, argument[i]);
                }

                likelihood.makeDirty();
                System.err.println("likelihood in numeric:" + likelihood.getLogLikelihood());
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return compoundSymmetricMatrix.getOffDiagonalParameter().getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return -1.0;
            }

            @Override
            public double getUpperBound(int n) {
                return 1.0;
            }
        };
    }

    @Override
    String checkNumeric(double[] analytic) {

        System.err.println("Numeric at: \n" + new Vector(compoundSymmetricMatrix.getOffDiagonalParameter().getParameterValues()));

        double[] storedValues = compoundSymmetricMatrix.getOffDiagonalParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(getNumeric(), storedValues);
        double[] testGradientTrans = compoundSymmetricMatrix.updateGradientCorrelation(testGradient);
        for (int i = 0; i < storedValues.length; ++i) {
            compoundSymmetricMatrix.getOffDiagonalParameter().setParameterValue(i, storedValues[i]);
        }

        return getReportString(analytic, testGradient, testGradientTrans);
    }

    @Override
    double[] getGradientParameter(double[] gradient) {
        return getGradientCorrelation(gradient);
    }

}
