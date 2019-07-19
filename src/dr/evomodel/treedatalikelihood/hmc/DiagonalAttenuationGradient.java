/*
 * DiagonalAttenuationGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class DiagonalAttenuationGradient implements GradientWrtParameterProvider, Reportable {

    private final Likelihood likelihood;
    private final int dim;
    private final BranchSpecificGradient branchSpecificGradient;

    private final DiagonalMatrix attenuation;

    public DiagonalAttenuationGradient(BranchSpecificGradient branchSpecificGradient,
                                       Likelihood likelihood,
                                       MatrixParameterInterface parameter) {

        assert (parameter instanceof DiagonalMatrix)
                : "DiagonalAttenuationGradient can only be applied to a DiagonalMatrix.";

        this.attenuation = (DiagonalMatrix) parameter;

        this.branchSpecificGradient = branchSpecificGradient;
        this.likelihood = likelihood;
        this.dim = parameter.getColumnDimension();

    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return attenuation.getDiagonalParameter();
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradient = branchSpecificGradient.getGradientLogDensity();
        return extractDiagonalGradient(gradient);
    }

    private double[] extractDiagonalGradient(double[] gradient) {
        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = gradient[i];
        }
        return result;
    }

    String getReportString(double[] analytic, double[] numeric) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric: " + new Vector(numeric) +
                "\n";
    }

    @Override
    public String getReport() {
        return checkNumeric(getGradientLogDensity());
    }

    MultivariateFunction getNumeric() {

        return new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    attenuation.setParameterValue(i, argument[i]);
                }

                likelihood.makeDirty();
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return attenuation.getColumnDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return 0.0;
            }

            @Override
            public double getUpperBound(int n) {
                return Double.POSITIVE_INFINITY;
            }
        };
    }

    String checkNumeric(double[] analytic) {

        System.err.println("Numeric at: \n" + new Vector(attenuation.getParameterValues()));

        double[] storedValues = attenuation.getDiagonalParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(getNumeric(), storedValues);
        for (int i = 0; i < storedValues.length; ++i) {
            attenuation.setParameterValue(i, storedValues[i]);
        }

        return getReportString(analytic, testGradient);
    }
}