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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.distributions.DirichletDistribution;
import dr.util.Attribute;
import dr.xml.Reportable;

import java.util.logging.Logger;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public interface GradientWrtParameterProvider extends NumericGradientStepSizeProvider {

    Likelihood getLikelihood();

    Parameter getParameter();

    int getDimension();

    double[] getGradientLogDensity();

    /**
     * @return S if the parameter is constrained to the simplex {x > 0, sum x = S}, otherwise NaN.
     * Gradient reports of simplex providers compare only the tangent part g - (x.g) 1 / S.
     */
    default double getSimplexSum() {
        return Double.NaN;
    }

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
        public double getSimplexSum() {
            return provider instanceof DirichletDistribution ?
                    ((DirichletDistribution) provider).getSumsTo() : Double.NaN;
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

        public MismatchException() {
            super();
        }

        public MismatchException(String message) {
            super(message);
        }
    }

    class CheckGradientNumerically {

        private final GradientWrtParameterProvider provider;
        private final Parameter parameter;
        private final double lowerBound;
        private final double upperBound;

        private final boolean checkValues;
        private final double tolerance;
        private final double smallThreshold;
        private final double numericStepSize;

        public CheckGradientNumerically(GradientWrtParameterProvider provider,
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
            this.numericStepSize = provider.getNumericGradientStepSize();
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
            double[] testGradient = new double[parameter.getDimension()];
            NumericalDerivative.gradient(numeric, parameter.getParameterValues(), testGradient, numericStepSize);

            setParameter(savedValues);
            return testGradient;
        }

        /**
         * Central differences along e_k - x / S, which keep sum x = S, so a Dirichlet on x stays finite.
         * The parameter is restored with one event afterwards.
         */
        public double[] getTangentNumerical(double simplexSum) throws MismatchException {

            final double[] x = parameter.getParameterValues();
            final int dim = x.length;

            double sum = 0.0;
            for (double v : x) {
                sum += v;
            }
            if (Math.abs(sum - simplexSum) > SIMPLEX_PRECONDITION * simplexSum) {
                throw new MismatchException("off-simplex: parameter " + parameter.getParameterName()
                        + " sums to " + sum + ", not " + simplexSum);
            }

            final double[] tangent = new double[dim];
            final double[] point = new double[dim];
            for (int k = 0; k < dim; ++k) {
                for (int i = 0; i < dim; ++i) {
                    point[i] = x[i] + TANGENT_STEP * ((i == k ? 1.0 : 0.0) - x[i] / simplexSum);
                }
                final double plus = numeric.evaluate(point);
                for (int i = 0; i < dim; ++i) {
                    point[i] = x[i] - TANGENT_STEP * ((i == k ? 1.0 : 0.0) - x[i] / simplexSum);
                }
                final double minus = numeric.evaluate(point);
                tangent[k] = (plus - minus) / (2 * TANGENT_STEP);
            }

            setParameter(x);
            return tangent;
        }

        public String getReport() throws MismatchException {

            double[] analytic = provider.getGradientLogDensity();

            final double simplexSum = provider.getSimplexSum();
            if (!Double.isNaN(simplexSum)) {
                return getSimplexReport(analytic, simplexSum);
            }

            double[] numeric = getNumericalGradient();

            return makeReport("Gradient\n", analytic, numeric, checkValues, tolerance, smallThreshold);
        }

        private String getSimplexReport(double[] analytic, double simplexSum) throws MismatchException {

            final double[] numericTangent = getTangentNumerical(simplexSum);
            final double[] analyticTangent = SimplexTangent.project(analytic, parameter.getParameterValues(), simplexSum);

            StringBuilder sb = new StringBuilder("Gradient (simplex, sum = " + simplexSum + ")\n");
            sb.append("analytic (tangent): ").append(new dr.math.matrixAlgebra.Vector(analyticTangent)).append("\n");
            sb.append("numeric (tangent): ").append(new dr.math.matrixAlgebra.Vector(numericTangent)).append("\n");
            sb.append("analytic (raw, extension-dependent, unchecked): ")
                    .append(new dr.math.matrixAlgebra.Vector(analytic)).append("\n");
            if (hasDirichletOn(provider.getLikelihood(), parameter)) {
                sb.append("numeric (raw): skipped (Dirichlet off-simplex)\n");
            } else {
                sb.append("numeric (raw, extension-dependent, unchecked): ")
                        .append(new dr.math.matrixAlgebra.Vector(getNumericalGradient())).append("\n");
            }

            if (checkValues) {
                String failure = checkLengthAndFinite(analyticTangent, numericTangent);
                for (int i = 0; failure == null && i < analyticTangent.length; ++i) {
                    final double a = analyticTangent[i];
                    final double n = numericTangent[i];
                    if (Math.abs(a - n) > tolerance * Math.max(1.0, Math.max(Math.abs(a), Math.abs(n)))) {
                        failure = "Tangent difference @ " + (i + 1) + ": " + a + " " + n;
                    }
                }
                if (failure != null) {
                    sb.append("\n").append(failure).append("\n");
                    Logger.getLogger("dr.inference.hmc").info(sb.toString());
                    throw new MismatchException(failure);
                }
            }

            return sb.toString();
        }

        private static final double TANGENT_STEP = 1E-6;
        private static final double SIMPLEX_PRECONDITION = 1E-10;

        private static boolean hasDirichletOn(Likelihood likelihood, Parameter parameter) {
            if (likelihood instanceof MultivariateDistributionLikelihood) {
                MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) likelihood;
                if (mdl.getDistribution() instanceof DirichletDistribution) {
                    for (Attribute<double[]> data : mdl.getDataList()) {
                        if (data == parameter) {
                            return true;
                        }
                    }
                }
            } else if (likelihood instanceof CompoundLikelihood) {
                for (Likelihood child : ((CompoundLikelihood) likelihood).getLikelihoods()) {
                    if (hasDirichletOn(child, parameter)) {
                        return true;
                    }
                }
            }
            return false;
        }

        private static String checkLengthAndFinite(double[] analytic, double[] numeric) {
            if (analytic.length != numeric.length) {
                return "Length mismatch: analytic " + analytic.length + ", numeric " + numeric.length;
            }
            for (int i = 0; i < analytic.length; ++i) {
                if (!Double.isFinite(analytic[i]) || !Double.isFinite(numeric[i])) {
                    return "Non-finite entry @ " + (i + 1) + ": " + analytic[i] + " " + numeric[i];
                }
            }
            return null;
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
            String failure = CheckGradientNumerically.checkLengthAndFinite(analytic, numeric);
            if (failure != null) {
                sb.append("\n").append(failure).append("\n");
                Logger.getLogger("dr.inference.hmc").info(sb.toString());
                throw new MismatchException(failure);
            }
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
