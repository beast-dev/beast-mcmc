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

import dr.inference.model.DerivativeOrder;
import dr.inference.model.DerivativeProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.List;
import java.util.logging.Logger;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public interface DerivativeWrtParameterProvider {

    Likelihood getLikelihood();

    Parameter getParameter();

    int getDimension(DerivativeOrder order);

    double[] getDerivativeLogDensity(DerivativeOrder order);

    DerivativeOrder getHighestOrder();

    static DerivativeOrder getHighestOrder(List<DerivativeWrtParameterProvider> providers) {

        if (providers.size() == 0) { return DerivativeOrder.ZEROTH; }

        DerivativeOrder highest = DerivativeOrder.FULL_HESSIAN;
        for (DerivativeWrtParameterProvider provider : providers) {
            if (provider.getHighestOrder().getValue() < highest.getValue()) {
                highest = provider.getHighestOrder();
            }
        }
        return highest;
    }

    class ParameterWrapper implements DerivativeWrtParameterProvider, Reportable {

        final DerivativeProvider provider;
        final Parameter parameter;
        final Likelihood likelihood;

        public ParameterWrapper(DerivativeProvider provider, Parameter parameter, Likelihood likelihood) {
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
        public int getDimension(DerivativeOrder order) {
            return order.getDerivativeDimension(parameter.getDimension());
        }

        @Override
        public double[] getDerivativeLogDensity(DerivativeOrder type) {

            assert (provider.getHighestOrder().getValue() >= type.getValue());

            return provider.getDerivativeLogDensity(parameter.getParameterValues(), type);
        }

        @Override
        public DerivativeOrder getHighestOrder() {
            return provider.getHighestOrder();
        }

        @Override
        public String getReport() {
//            return getReportAndCheckForError(this, parameter.getBounds().getLowerLimit(0),
//                    parameter.getBounds().getUpperLimit(0), null);
            return null;
        }
    }

    class MismatchException extends Exception { }

    class CheckDerivativeNumerically {

        private final DerivativeWrtParameterProvider provider;
        private final DerivativeOrder type;
        private final Parameter parameter;
        private final double lowerBound;
        private final double upperBound;

        private final boolean checkValues;
        private final double tolerance;

        CheckDerivativeNumerically(DerivativeWrtParameterProvider provider,
                                   DerivativeOrder type,
                                   double lowerBound, double upperBound,
                                   Double nullableTolerance) {
            this.provider = provider;
            this.type = type;
            this.parameter = provider.getParameter();
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;

            this.checkValues = nullableTolerance != null;
            this.tolerance = checkValues ? nullableTolerance : 0.0;
        }

        private MultivariateFunction numeric = new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                setParameter(argument);

                if (type == DerivativeOrder.GRADIENT) {
                    return provider.getLikelihood().getLogLikelihood();
                } else {
                    throw new RuntimeException("Not yet implemented");
                }
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

        private double[] getNumericalGradient() {

            double[] savedValues = parameter.getParameterValues();
            double[] testGradient = NumericalDerivative.gradient(numeric, parameter.getParameterValues());

            setParameter(savedValues);
            return testGradient;
        }

        public String getReport() throws MismatchException {

            double[] analytic = provider.getDerivativeLogDensity(type);
            double[] numeric = getNumericalGradient();

            return makeReport("Gradient\n", analytic, numeric, checkValues, tolerance);
        }
    }

    static String makeReport(String header,
                             double[] analytic,
                             double[] numeric,
                             boolean checkValues,
                             double tolerance) throws MismatchException {

        StringBuilder sb = new StringBuilder(header);
        sb.append("analytic: ").append(new dr.math.matrixAlgebra.Vector(analytic));
        sb.append("\n");
        sb.append("numeric : ").append(new dr.math.matrixAlgebra.Vector(numeric));
        sb.append("\n");

        if (checkValues) {
            for (int i = 0; i < analytic.length; ++i) {
                double relativeDifference = 2 * (analytic[i] - numeric[i]) / (analytic[i] + numeric[i]);
                if (Math.abs(relativeDifference) > tolerance) {
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

    static String getReportAndCheckForError(DerivativeWrtParameterProvider provider,
                                            DerivativeOrder type,
                                            double lowerBound, double upperBound,
                                            Double nullableTolerance) {
        String report;
        try {
            report = new CheckDerivativeNumerically(provider, type,
                    lowerBound, upperBound,
                    nullableTolerance
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
