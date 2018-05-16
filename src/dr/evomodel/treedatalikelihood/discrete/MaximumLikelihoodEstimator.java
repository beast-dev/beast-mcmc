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
import dr.inference.model.Parameter;
import dr.xml.Reportable;
import com.github.lbfgs4j.LbfgsMinimizer;
import com.github.lbfgs4j.liblbfgs.Function;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class MaximumLikelihoodEstimator implements Reportable {
    private final GradientWrtParameterProvider gradientWrtParameterProvider;
    private final Parameter parameter;

    public MaximumLikelihoodEstimator(GradientWrtParameterProvider gradientWrtParameterProvider) {
        this.gradientWrtParameterProvider = gradientWrtParameterProvider;
        this.parameter = gradientWrtParameterProvider.getParameter();
    }


    @Override
    public String getReport() {

        Function f = new Function() {
            public int getDimension() {
                return gradientWrtParameterProvider.getDimension();
            }

            public double valueAt(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, Math.exp(argument[i]));
                }

                return -gradientWrtParameterProvider.getLikelihood().getLogLikelihood();
            }

            public double[] gradientAt(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, Math.exp(argument[i]));
                }

                double lnL = valueAt(argument);

                double[] gradient = gradientWrtParameterProvider.getGradientLogDensity();

                for (int i = 0; i < argument.length; ++i) {
                    gradient[i] *= -Math.exp(argument[i]);
                }

                return gradient;
            }
        };
        boolean printScreen = false;
        LbfgsMinimizer minimizer = new LbfgsMinimizer(printScreen);
        double[] x = minimizer.minimize(f);
        double min = f.valueAt(x);

        StringBuilder sb = new StringBuilder();
        sb.append("lnX: ").append(new dr.math.matrixAlgebra.Vector(x)).append("\n");
        sb.append("Gradient: ").append(new dr.math.matrixAlgebra.Vector(f.gradientAt(x))).append("\n");
        sb.append("fx: ").append(String.valueOf(min)).append("\n");

        return sb.toString();
    }

}
