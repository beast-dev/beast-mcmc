/*
 * MaximumLikelihoodEstimator.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.MultivariateMinimum;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;
import com.github.lbfgs4j.LbfgsMinimizer;
import com.github.lbfgs4j.liblbfgs.Function;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class MaximizerWrtParameter implements Reportable {

    private final GradientWrtParameterProvider gradient;
    private final Parameter parameter;
    private final Likelihood likelihood;
    private final Transform transform;

    private Function function = null;
    private long time = 0;
    private double minimumValue = Double.NaN;
    private double[] minimumPoint = null;

    public MaximizerWrtParameter(Likelihood likelihood,
                                 Parameter parameter,
                                 GradientWrtParameterProvider gradient,
                                 Transform transform) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.gradient = gradient;
        this.transform = transform;

        // If gradient == null, then construct a numerical version using NumericalDerivative.gradient()
    }

    public void maximize() {

        function = new Function() {
            public int getDimension() {
                return gradient.getDimension();
            }

            public double valueAt(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, Math.exp(argument[i])); // TODO Handle transform
                }

                return -gradient.getLikelihood().getLogLikelihood();
            }

            public double[] gradientAt(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, Math.exp(argument[i])); // TODO Handle transform
                }

                double lnL = valueAt(argument);

                double[] gradient = MaximizerWrtParameter.this.gradient.getGradientLogDensity();

                for (int i = 0; i < argument.length; ++i) {
                    gradient[i] *= -Math.exp(argument[i]); // TODO Handle transform
                }

                return gradient;
            }
        };

        boolean printScreen = false;
        LbfgsMinimizer minimizer = new LbfgsMinimizer(printScreen);

        minimumPoint = minimizer.minimize(function);
        minimumValue = function.valueAt(minimumPoint);

        setParameter(new WrappedVector.Raw(minimumPoint), parameter);
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        if (function == null) {
            sb.append("Not yet executed.");
        } else {

            sb.append("lnX: ").append(new dr.math.matrixAlgebra.Vector(minimumPoint)).append("\n");
            sb.append("Gradient: ").append(new dr.math.matrixAlgebra.Vector(function.gradientAt(minimumPoint))).append("\n");
            sb.append("fx: ").append(String.valueOf(minimumValue)).append("\n");

        }

        return sb.toString();
    }

}
