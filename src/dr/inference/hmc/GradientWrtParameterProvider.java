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

import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.logging.Logger;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public interface GradientWrtParameterProvider {

    Likelihood getLikelihood();

    Parameter getParameter();

    int getDimension();

    double[] getGradientLogDensity();

    class ParameterWrapper implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

        final GradientProvider provider;
        final Parameter parameter;
        final Likelihood likelihood;

        public ParameterWrapper(GradientProvider provider, Parameter parameter, Likelihood likelihood) {
            this.provider = provider;
            this.parameter = parameter;
            this.likelihood = likelihood;
        }

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
            return provider.getGradientLogDensity(parameter.getParameterValues());
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
            return GradientWrtParameterProvider.getReportAndCheckForError(this,
                    parameter.getBounds().getLowerLimit(0),
                    parameter.getBounds().getUpperLimit(0), null);
        }
    }

    class MismatchException extends Exception {
    }

    class CheckGradientNumerically {

        private final GradientWrtParameterProvider provider;
        private final Parameter parameter;
        private final double lowerBound;
        private final double upperBound;

        private final boolean checkValues;
        private final double tolerance;
        private final double smallThreshold;

        CheckGradientNumerically(GradientWrtParameterProvider provider,
                                 double lowerBound, double upperBound,
                                 Double nullableTolerance,
                                 Double nullableSmallNumberThreshold) {
            this.provider = provider;
            this.parameter = provider.getParameter();
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;

            this.checkValues = nullableTolerance != null;
            this.tolerance = checkValues ? nullableTolerance : 0.0;

            this.smallThreshold = nullableSmallNumberThreshold != null ? nullableSmallNumberThreshold : 0.0;
        }


        private MultivariateFunction numeric = new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                setParameter(argument);
                return provider.getLikelihood().getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return lowerBound;
            }

            @Override
            public double getUpperBound(int n) {
                return upperBound;
            }
        };

        private void setParameter(double[] values) {

            for (int i = 0; i < values.length; ++i) {
                parameter.setParameterValueQuietly(i, values[i]);
            }

            parameter.fireParameterChangedEvent();
        }

        public double[] getNumericalGradient() {

            double[] savedValues = parameter.getParameterValues();
            double[] testGradient = NumericalDerivative.gradient(numeric, parameter.getParameterValues());

            setParameter(savedValues);
            return testGradient;
        }

        public String getReport() throws MismatchException {

            double[] analytic = provider.getGradientLogDensity();
            double[] numeric = getNumericalGradient();

            return makeReport("Gradient\n", analytic, numeric, checkValues, tolerance, smallThreshold);
        }
    }

    static String makeReport(String header,
                             double[] analytic,
                             double[] numeric,
                             boolean checkValues,
                             double tolerance,
                             double smallNumberThreshold) throws MismatchException {

        StringBuilder sb = new StringBuilder(header);
        sb.append("analytic: ").append(new dr.math.matrixAlgebra.Vector(analytic));
        sb.append("\n");
        sb.append("numeric : ").append(new dr.math.matrixAlgebra.Vector(numeric));
        sb.append("\n");

        if (checkValues) {
            for (int i = 0; i < analytic.length; ++i) {
                double relativeDifference = 2 * (analytic[i] - numeric[i]) / (analytic[i] + numeric[i]);
                boolean testFailed = Math.abs(relativeDifference) > tolerance &&
                        Math.abs(analytic[i]) > smallNumberThreshold && Math.abs(numeric[i]) > smallNumberThreshold ||
                        ((analytic[i] == 0.0 || numeric[i] == 0.0) && Math.abs(analytic[i] + numeric[i]) > tolerance);

                if (testFailed) {
                    sb.append("\nDifference @ ").append(i + 1).append(": ")
                            .append(analytic[i]).append(" ").append(numeric[i])
                            .append(" ").append(relativeDifference).append("\n");
                    Logger.getLogger("dr.inference.hmc").info(sb.toString());
                    throw new MismatchException();
                }
            }
        }

        return sb.toString();
    }

    static String getReportAndCheckForError(GradientWrtParameterProvider provider,
                                            double lowerBound, double upperBound,
                                            Double nullableTolerance) {
        return getReportAndCheckForError(provider, lowerBound, upperBound, nullableTolerance, null);
    }

    static String getReportAndCheckForError(GradientWrtParameterProvider provider,
                                            double lowerBound, double upperBound,
                                            Double nullableTolerance,
                                            Double nullableSmallNumberThreshold) {
        String report;
        try {
            report = new CheckGradientNumerically(provider,
                    lowerBound, upperBound,
                    nullableTolerance, nullableSmallNumberThreshold
            ).getReport();
        } catch (MismatchException e) {
            String message = e.getMessage();
            if (message == null) {
                message = provider.getParameter().getParameterName();
            }
            if (message == null) {
                message = "Gradient check failure";
            }
            throw new RuntimeException(message);
        }

        return report;
    }

    Double TOLERANCE = 1E-1;
}
