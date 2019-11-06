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
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Timer;
import dr.util.Transform;
import dr.xml.Reportable;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.LbfgsMinimizer;
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class MaximizerWrtParameter implements Reportable {

    private final GradientWrtParameterProvider gradient;
    private final GradientType gradientType;
    private final Parameter parameter;
    private final Likelihood likelihood;
    private final Transform transform;
    private final Function function;
    private final Settings settings;

    private double time = 0.0;
    private long count = 0;
    private double minimumValue = Double.NaN;
    private double[] minimumPoint = null;

    public static class Settings {
        int numberIterations;
        boolean startAtCurrentState;
        boolean printToScreen;

        public Settings(int numberIterations, boolean startAtCurrentState, boolean printToScreen) {
            this.numberIterations = numberIterations;
            this.startAtCurrentState = startAtCurrentState;
            this.printToScreen = printToScreen;
        }
    }

    public MaximizerWrtParameter(Likelihood likelihood,
                                 Parameter parameter,
                                 GradientWrtParameterProvider gradient,
                                 Transform transform,
                                 Settings settings) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.transform = transform;

        if (gradient == null) {
            this.gradient = constructGradient();
            this.gradientType = GradientType.NUMERICAL;
        } else {
            this.gradient = gradient;
            this.gradientType = GradientType.ANALYTIC;
        }

        this.function = constructFunction();
        this.settings = settings;
    }

    public Likelihood getLikelihood() {
        return likelihood;
    }

    public void maximize() {

        LBFGS_Param paramsBFGS = Lbfgs.defaultParams();

        if (settings.numberIterations > 0) {
            paramsBFGS.max_iterations = settings.numberIterations;
        }

        LbfgsMinimizer minimizer = new LbfgsMinimizer(paramsBFGS, settings.printToScreen);
        double[] x0 = null;

        if (settings.startAtCurrentState) {
            x0 = parameter.getParameterValues();

            if (transform != null) {
                x0 = transform.transform(x0, 0, x0.length);
            }

        }

        Timer timer = new Timer();
        timer.start();

        minimumPoint = minimizer.minimize(function, x0);

        timer.stop();

        time += timer.toSeconds();
        minimumValue = function.valueAt(minimumPoint);

        if (transform != null) {
            setParameter(new WrappedVector.Raw(transform.inverse(minimumPoint, 0, minimumPoint.length)), parameter);
        } else {
            setParameter(new WrappedVector.Raw(minimumPoint), parameter);
        }
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        if (function == null) {
            sb.append("Not yet executed.");
        } else {

            if (transform != null) {
                sb.append("Gradient is taken with respect to the transformed paramter values.\n");
                sb.append("Untransformed X: ").append(new dr.math.matrixAlgebra.Vector(transform.inverse(minimumPoint, 0, minimumPoint.length))).append("\n");
            }
            sb.append("X: ").append(new dr.math.matrixAlgebra.Vector(minimumPoint)).append("\n");
            sb.append("Gradient: ").append(new dr.math.matrixAlgebra.Vector(function.gradientAt(minimumPoint))).append("\n");
            sb.append("Gradient type: ").append(gradientType).append("\n");
            sb.append("Fx: ").append(minimumValue).append("\n");
            sb.append("Time: ").append(time).append("s\n");
            sb.append("Count: ").append(count).append("\n");

        }

        return sb.toString();
    }

    private double evaluateLogLikelihood() {
        ++count;
        return likelihood.getLogLikelihood();
    }

    private Function constructFunction() {

        return new Function() {

            @Override
            public int getDimension() {
                return gradient.getDimension();
            }

            @Override
            public double valueAt(double[] argument) {

                if (transform != null) {
                    argument = transform.inverse(argument, 0, argument.length);
                }

                setParameter(new WrappedVector.Raw(argument), parameter);
                return -evaluateLogLikelihood();
            }

            @Override
            public double[] gradientAt(double[] argument) {

                if (transform != null) {
                    argument = transform.inverse(argument, 0, argument.length);
                }

                setParameter(new WrappedVector.Raw(argument), parameter);

                double[] result = gradient.getGradientLogDensity();

                if (transform != null) {

                    result = transform.updateGradientUnWeightedLogDensity(result, argument, 0, argument.length);

                }
                for (int i = 0; i < result.length; ++i) {
                    result[i] = -result[i];
                }

                return result;
            }
        };
    }

    private GradientWrtParameterProvider constructGradient() {

        final MultivariateFunction function = new MultivariateFunction() {
            @Override
            public double evaluate(double[] argument) {

                setParameter(new WrappedVector.Raw(argument), parameter);
                return evaluateLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound(int n) {
                return Double.POSITIVE_INFINITY;
            }
        };

        return new GradientWrtParameterProvider() {

            @Override
            public Likelihood getLikelihood() {
                return likelihood;
            }

            @Override
            public Parameter getParameter() {
                return parameter;
            }

            @Override
            public int getDimension() {
                return parameter.getDimension();
            }

            @Override
            public double[] getGradientLogDensity() {
                return NumericalDerivative.gradient(function, parameter.getParameterValues());
            }
        };
    }

    private enum GradientType {
        ANALYTIC("analytic"),
        NUMERICAL("numerical");

        private String type;

        GradientType(String type) {
            this.type = type;
        }

        public String toString() { return type; }
    }
}
