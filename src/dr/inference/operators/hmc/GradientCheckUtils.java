/*
 * GradientCheckUtils.java
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

package dr.inference.operators.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

import java.util.ArrayList;

/**
 * Numeric-vs-analytic gradient check, following the same convention as
 * {@code HamiltonianMonteCarloOperator.checkGradient} (compare an analytic
 * gradient against {@link NumericalDerivative#gradient}, throw with detailed
 * mismatch information on failure), factored out so discontinuous HMC
 * operators don't each reimplement it.
 *
 * Unlike the plain HMC operator, discontinuous operators must restrict the
 * check to a subset of coordinates: {@code joint.getLogLikelihood()} is not
 * differentiable across a discontinuous coordinate, so finite-differencing
 * it there would compare a real gradient against a meaningless one.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
class GradientCheckUtils {

    private GradientCheckUtils() {
    }

    /**
     * Overrides {@link NumericalDerivative}'s default centered-difference step ratio
     * ({@code EPSILON^(1/3)}, the standard truncation/round-off balance point) when set.
     * Diagnostic-only: exists to test whether a numeric-vs-analytic gradient mismatch is
     * finite-difference truncation error (should shrink with a smaller step) versus a real
     * analytic-gradient bug (should not). Unset by default so existing validated checks keep
     * their established step size.
     */
    private static final String STEP_SIZE_RATIO_PROPERTY =
            "dr.inference.operators.hmc.GradientCheckUtils.stepSizeRatio";

    static void checkGradient(final Likelihood joint, final Parameter parameter,
                              double[] analyticalGradient, boolean[] checkMask, double tolerance) {
        checkGradient(joint, parameter, analyticalGradient, checkMask, tolerance, null);
    }

    static void checkGradient(final Likelihood joint, final Parameter parameter,
                              double[] analyticalGradient, boolean[] checkMask,
                              double tolerance, final Transform transform) {

        final int dimension = parameter.getDimension();
        if (analyticalGradient.length != dimension || checkMask.length != dimension) {
            throw new RuntimeException("Unequal dimensions");
        }

        final double[] currentValues = parameter.getParameterValues();
        final double[] currentArgument = transform == null ? currentValues :
                transform.transform(currentValues, 0, dimension);
        final int[] activeIndices = activeIndices(checkMask);

        MultivariateFunction numeric = new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {
                double[] fullArgument = currentArgument.clone();
                for (int k = 0; k < activeIndices.length; k++) {
                    fullArgument[activeIndices[k]] = argument[k];
                }
                double[] fullValues = transform == null ? fullArgument :
                        transform.inverse(fullArgument, 0, fullArgument.length);
                ReadableVector.Utils.setParameter(fullValues, parameter);
                return transform == null ? joint.getLogLikelihood() :
                        joint.getLogLikelihood() - transform.logJacobian(fullValues, 0, fullValues.length);
            }

            @Override
            public int getNumArguments() {
                return activeIndices.length;
            }

            @Override
            public double getLowerBound(int n) {
                if (transform != null) {
                    return Double.NEGATIVE_INFINITY;
                }
                return parameter.getBounds().getLowerLimit(activeIndices[n]);
            }

            @Override
            public double getUpperBound(int n) {
                if (transform != null) {
                    return Double.POSITIVE_INFINITY;
                }
                return parameter.getBounds().getUpperLimit(activeIndices[n]);
            }
        };

        double[] activeArgument = new double[activeIndices.length];
        double[] analyticalActive = new double[activeIndices.length];
        double[] analyticalArgumentGradient = transform == null ? analyticalGradient :
                transform.updateGradientLogDensity(analyticalGradient, currentValues, 0, dimension);
        for (int k = 0; k < activeIndices.length; k++) {
            activeArgument[k] = currentArgument[activeIndices[k]];
            analyticalActive[k] = analyticalArgumentGradient[activeIndices[k]];
        }

        final String stepSizeRatioProperty = System.getProperty(STEP_SIZE_RATIO_PROPERTY);
        double[] numericActive = new double[activeIndices.length];
        if (stepSizeRatioProperty == null) {
            numericActive = NumericalDerivative.gradient(numeric, activeArgument);
        } else {
            NumericalDerivative.gradient(numeric, activeArgument, numericActive,
                    Double.parseDouble(stepSizeRatioProperty));
        }
        boolean isClose = MathUtils.isClose(analyticalActive, numericActive, tolerance);

        ReadableVector.Utils.setParameter(currentValues, parameter);

        if (!isClose) {
            String message = "Gradients do not match on the checked (non-discontinuous) coordinates:\n" +
                    "\tAnalytic: " + new WrappedVector.Raw(analyticalActive) + "\n" +
                    "\tNumeric : " + new WrappedVector.Raw(numericActive) + "\n" +
                    mismatchInformation(activeIndices, analyticalActive, numericActive, tolerance);
            throw new RuntimeException(message);
        }
    }

    private static int[] activeIndices(boolean[] mask) {
        int count = 0;
        for (boolean active : mask) {
            if (active) {
                count++;
            }
        }
        int[] indices = new int[count];
        int k = 0;
        for (int i = 0; i < mask.length; i++) {
            if (mask[i]) {
                indices[k++] = i;
            }
        }
        return indices;
    }

    private static String mismatchInformation(int[] activeIndices, double[] analyticGradient,
                                              double[] numericGradient, double tolerance) {
        int n = analyticGradient.length;
        double maxDiff = 0;
        int maxPosition = -1;
        double meanDiff = 0;
        ArrayList<Integer> overPositions = new ArrayList<>();
        double[] absDiffs = new double[n];

        for (int i = 0; i < n; i++) {
            double absDiff = Math.abs(analyticGradient[i] - numericGradient[i]);
            absDiffs[i] = absDiff;
            meanDiff += absDiff;
            if (absDiff > tolerance) {
                overPositions.add(i);
            }
            if (absDiff > maxDiff) {
                maxDiff = absDiff;
                maxPosition = i;
            }
        }
        meanDiff /= n;

        StringBuilder sb = new StringBuilder();
        sb.append("\tMaximum absolute difference: ").append(maxDiff)
                .append(" (at parameter index ").append(maxPosition >= 0 ? activeIndices[maxPosition] : -1).append(")\n");
        sb.append("\tAverage absolute difference: ").append(meanDiff).append("\n");
        sb.append("\tList of all checked coordinates exceeding the tolerance:\n");
        sb.append("\t\tparameterIndex    analytic    numeric    absoluteDifference\n");

        for (int position : overPositions) {
            sb.append("\t\t").append(activeIndices[position]).append("    ")
                    .append(analyticGradient[position]).append("    ")
                    .append(numericGradient[position]).append("    ")
                    .append(absDiffs[position]).append("\n");
        }

        return sb.toString();
    }
}
