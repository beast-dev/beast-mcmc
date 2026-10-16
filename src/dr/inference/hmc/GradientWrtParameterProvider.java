/*
 * GradientWrtParameterProvider.java
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

package dr.inference.hmc;

import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.util.Attribute;
import dr.xml.Reportable;

import java.util.List;
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

    class Negative implements GradientWrtParameterProvider {

        private final GradientWrtParameterProvider provider;

        public Negative(GradientWrtParameterProvider provider) {
            this.provider = provider;
        }

        @Override
        public Likelihood getLikelihood() {
            return provider.getLikelihood();
        }

        @Override
        public Parameter getParameter() {
            return provider.getParameter();
        }

        @Override
        public int getDimension() {
            return provider.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {

            double[] gradient = provider.getGradientLogDensity();
            for (int i = 0; i < gradient.length; ++i) {
                gradient[i] =-gradient[i];
            }

            return gradient;
        }
    }

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

    class GaussianMeanWrapper implements GradientWrtParameterProvider, Reportable {

        final GradientProvider provider;
        final Parameter mean;
        final List<Attribute<double[]>> dataList;
        final Likelihood likelihood;

        public GaussianMeanWrapper(GradientProvider provider, Parameter mean,
                                   List<Attribute<double[]>> dataList, Likelihood likelihood) {
            this.provider = provider;
            this.mean = mean;
            this.dataList = dataList;
            this.likelihood = likelihood;
        }

        @Override
        public Likelihood getLikelihood() {
            return likelihood;
        }

        @Override
        public Parameter getParameter() {
            return mean;
        }

        @Override
        public int getDimension() {
            return mean.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {

            // d/dx log N(x | mean, P) = P (mean - x) is linear in x, so
            // d/d mean sum_r log N(x_r | mean, P) = P sum_r (x_r - mean) = -R * d/dx log N(xBar | mean, P)
            final int dim = mean.getDimension();
            final int R = dataList.size();

            final double[] xBar = new double[dim];
            for (Attribute<double[]> data : dataList) {
                final double[] x = data.getAttributeValue();
                for (int i = 0; i < dim; ++i) {
                    xBar[i] += x[i];
                }
            }
            for (int i = 0; i < dim; ++i) {
                xBar[i] /= R;
            }

            final double[] gradient = provider.getGradientLogDensity(xBar);
            for (int i = 0; i < dim; ++i) {
                gradient[i] *= -R;
            }

            return gradient;
        }

        @Override
        public String getReport() {
            return GradientWrtParameterProvider.getReportAndCheckForError(this,
                    Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null);
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


        private final MultivariateFunction numeric = new MultivariateFunction() {

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
