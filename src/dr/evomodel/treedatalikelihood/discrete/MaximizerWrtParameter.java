/*
 * MaximizerWrtParameter.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ScoreContributionsWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Timer;
import dr.util.FixedWhiteningTransform;
import dr.util.Transform;
import dr.xml.Reportable;
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Evaluate;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Progress;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue;
import com.github.lbfgs4j.liblbfgs.MutableDouble;

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
    private final Transform nativeTransform;
    private Transform transform;
    private Function function;
    private final Settings settings;
    private final GradientWrtParameterProvider fisherGradient;
    private final GradientWrtParameterProvider fisherPriorGradient;

    private double time = 0.0;
    private long count = 0;
    private double minimumValue = Double.NaN;
    private double[] minimumPoint = null;
    private ReturnValue returnValue = null;
    private int fisherRefreshCount = 0;
    private int fisherRank = -1;
    private int fisherNullity = -1;
    private double fisherMinimumPositiveEigenvalue = Double.NaN;
    private double fisherMaximumEigenvalue = Double.NaN;
    private boolean fisherUsedEmpiricalScores = false;

    public static class Settings {
        int numberIterations;
        double stepSize;
        boolean startAtCurrentState;
        boolean printToScreen;
        boolean includeJacobian;
        Double lbfgsEpsilon;
        Integer lbfgsPast;
        Double lbfgsDelta;
        Double lbfgsFtol;
        Double lbfgsGtol;
        Double lbfgsXtol;
        Integer lbfgsMaxLineSearch;
        boolean fisherRefreshOnFailure;
        int fisherMaximumRefreshes;
        double fisherGradientThreshold;
        double fisherEigenFloorFraction;
        double fisherNullPrecision;
        double fisherTrustRadius;

        public Settings(int numberIterations, double stepSize, boolean startAtCurrentState, boolean printToScreen, boolean includeJacobian) {
            this(numberIterations, stepSize, startAtCurrentState, printToScreen, includeJacobian,
                    null, null, null, null, null, null, null,
                    false, 0, 1.0e-2, 10.0, 1.0, Double.POSITIVE_INFINITY);
        }

        public Settings(int numberIterations,
                        double stepSize,
                        boolean startAtCurrentState,
                        boolean printToScreen,
                        boolean includeJacobian,
                        Double lbfgsEpsilon,
                        Integer lbfgsPast,
                        Double lbfgsDelta,
                        Double lbfgsFtol,
                        Double lbfgsGtol,
                        Double lbfgsXtol,
                        Integer lbfgsMaxLineSearch) {
            this(numberIterations, stepSize, startAtCurrentState, printToScreen, includeJacobian,
                    lbfgsEpsilon, lbfgsPast, lbfgsDelta, lbfgsFtol, lbfgsGtol, lbfgsXtol,
                    lbfgsMaxLineSearch, false, 0, 1.0e-2, 10.0, 1.0,
                    Double.POSITIVE_INFINITY);
        }

        public Settings(int numberIterations,
                        double stepSize,
                        boolean startAtCurrentState,
                        boolean printToScreen,
                        boolean includeJacobian,
                        Double lbfgsEpsilon,
                        Integer lbfgsPast,
                        Double lbfgsDelta,
                        Double lbfgsFtol,
                        Double lbfgsGtol,
                        Double lbfgsXtol,
                        Integer lbfgsMaxLineSearch,
                        boolean fisherRefreshOnFailure,
                        int fisherMaximumRefreshes,
                        double fisherGradientThreshold,
                        double fisherEigenFloorFraction,
                        double fisherNullPrecision,
                        double fisherTrustRadius) {
            this.numberIterations = numberIterations;
            this.stepSize = stepSize;
            this.startAtCurrentState = startAtCurrentState;
            this.printToScreen = printToScreen;
            this.includeJacobian = includeJacobian;
            this.lbfgsEpsilon = lbfgsEpsilon;
            this.lbfgsPast = lbfgsPast;
            this.lbfgsDelta = lbfgsDelta;
            this.lbfgsFtol = lbfgsFtol;
            this.lbfgsGtol = lbfgsGtol;
            this.lbfgsXtol = lbfgsXtol;
            this.lbfgsMaxLineSearch = lbfgsMaxLineSearch;
            this.fisherRefreshOnFailure = fisherRefreshOnFailure;
            this.fisherMaximumRefreshes = fisherMaximumRefreshes;
            this.fisherGradientThreshold = fisherGradientThreshold;
            this.fisherEigenFloorFraction = fisherEigenFloorFraction;
            this.fisherNullPrecision = fisherNullPrecision;
            this.fisherTrustRadius = fisherTrustRadius;
        }
    }

    public MaximizerWrtParameter(Likelihood likelihood,
                                 Parameter parameter,
                                 GradientWrtParameterProvider gradient,
                                 Transform transform,
                                 Settings settings) {
        this(likelihood, parameter, gradient, transform, settings, null, null);
    }

    public MaximizerWrtParameter(Likelihood likelihood,
                                 Parameter parameter,
                                 GradientWrtParameterProvider gradient,
                                 Transform transform,
                                 Settings settings,
                                 GradientWrtParameterProvider fisherGradient,
                                 GradientWrtParameterProvider fisherPriorGradient) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.nativeTransform = transform;
        this.transform = transform;

        if (gradient == null) {
            this.gradient = constructGradient();
            this.gradientType = GradientType.NUMERICAL;
        } else {
            this.gradient = gradient;
            this.gradientType = GradientType.ANALYTIC;
        }
        this.fisherGradient = fisherGradient == null ? this.gradient : fisherGradient;
        this.fisherPriorGradient = fisherPriorGradient;

        this.settings = settings;
        validateFisherSettings();
        this.function = constructFunction();
    }

    public Likelihood getLikelihood() {
        return likelihood;
    }

    public double[] getMinimumPoint(boolean inParameterSpace) {
        if (inParameterSpace && (transform != null)) {
            return transform.inverse(minimumPoint, 0, minimumPoint.length);
        } else {
            return minimumPoint;
        }
    }

    public GradientWrtParameterProvider getGradient() {
        return gradient;
    }

    public boolean wasExecuted() {
        return function != null;
    }

    public Transform getTransform() {
        return transform;
    }

    public Function getFunction() {
        return function;
    }

    public ReturnValue getReturnValue() {
        return returnValue;
    }

    public int getFisherRefreshCount() {
        return fisherRefreshCount;
    }

    public int getFisherRank() {
        return fisherRank;
    }

    public int getFisherNullity() {
        return fisherNullity;
    }

    public boolean fisherUsedEmpiricalScores() {
        return fisherUsedEmpiricalScores;
    }

    public void maximize() {

        LBFGS_Param paramsBFGS = Lbfgs.defaultParams();

        if (settings.numberIterations > 0) {
            paramsBFGS.max_iterations = settings.numberIterations;
        }
        if (settings.lbfgsEpsilon != null) {
            paramsBFGS.epsilon = settings.lbfgsEpsilon;
        }
        if (settings.lbfgsPast != null) {
            paramsBFGS.past = settings.lbfgsPast;
        }
        if (settings.lbfgsDelta != null) {
            paramsBFGS.delta = settings.lbfgsDelta;
        }
        if (settings.lbfgsFtol != null) {
            paramsBFGS.ftol = settings.lbfgsFtol;
        }
        if (settings.lbfgsGtol != null) {
            paramsBFGS.gtol = settings.lbfgsGtol;
        }
        if (settings.lbfgsXtol != null) {
            paramsBFGS.xtol = settings.lbfgsXtol;
        }
        if (settings.lbfgsMaxLineSearch != null) {
            paramsBFGS.max_linesearch = settings.lbfgsMaxLineSearch;
        }

        Timer timer = new Timer();
        timer.start();

        double[] bestPointInParameterSpace = parameter.getParameterValues();
        double bestValue = objectiveAtParameter(bestPointInParameterSpace);
        boolean firstStage = true;
        while (true) {
            double[] x0 = null;
            if (settings.startAtCurrentState || !firstStage) {
                x0 = parameter.getParameterValues();
                if (transform != null) {
                    x0 = transform.transform(x0, 0, x0.length);
                }
            }

            minimumPoint = minimize(paramsBFGS, x0);
            setParameterFromOptimizationPoint(minimumPoint);
            final double stageValue = function.valueAt(minimumPoint);
            if (Double.isFinite(stageValue) && stageValue < bestValue) {
                bestValue = stageValue;
                bestPointInParameterSpace = parameter.getParameterValues();
            }
            final double[] stageGradient = function.gradientAt(minimumPoint);
            final double gradientNorm = euclideanNorm(stageGradient);
            printLargestGradientComponents(stageGradient);
            if (!shouldRefresh(returnValue, gradientNorm)) {
                break;
            }

            setParameter(new WrappedVector.Raw(bestPointInParameterSpace), parameter);
            if (!refreshFisherWhitening()) {
                break;
            }
            ++fisherRefreshCount;
            firstStage = false;
            if (settings.printToScreen) {
                System.out.printf("Fisher refresh %d/%d: rank=%d, nullity=%d, eigenvalue range=[%.6g, %.6g]%n",
                        fisherRefreshCount, settings.fisherMaximumRefreshes, fisherRank, fisherNullity,
                        fisherMinimumPositiveEigenvalue, fisherMaximumEigenvalue);
            }
        }

        timer.stop();

        time += timer.toSeconds();
        setParameter(new WrappedVector.Raw(bestPointInParameterSpace), parameter);
        minimumPoint = parameter.getParameterValues();
        if (transform != null) {
            minimumPoint = transform.transform(minimumPoint, 0, minimumPoint.length);
        }
        minimumValue = bestValue;
    }

    private void printLargestGradientComponents(double[] stageGradient) {
        if (!settings.printToScreen || !Boolean.getBoolean("beast.maximizer.printGradientComponents")) {
            return;
        }
        final Integer[] order = new Integer[stageGradient.length];
        for (int i = 0; i < order.length; ++i) {
            order[i] = i;
        }
        java.util.Arrays.sort(order, (left, right) ->
                Double.compare(Math.abs(stageGradient[right]), Math.abs(stageGradient[left])));
        final int requested = Integer.getInteger("beast.maximizer.gradientComponentCount", 20);
        final int count = Math.min(Math.max(requested, 0), order.length);
        System.out.printf("Largest transformed-gradient components (top %d of %d):%n", count, order.length);
        for (int rank = 0; rank < count; ++rank) {
            final int index = order[rank];
            System.out.printf("  %d\t%s\t%.12g%n", index, parameter.getDimensionName(index), stageGradient[index]);
        }
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        if (function == null) {
            sb.append("Not yet executed.");
        } else {

            if (transform != null) {
                sb.append("Gradient is taken with respect to the transformed parameter values.\n");
                sb.append("Untransformed X: ").append(new dr.math.matrixAlgebra.Vector(transform.inverse(minimumPoint, 0, minimumPoint.length))).append("\n");
            }
            sb.append("X: ").append(new dr.math.matrixAlgebra.Vector(minimumPoint)).append("\n");
            sb.append("Gradient: ").append(new dr.math.matrixAlgebra.Vector(function.gradientAt(minimumPoint))).append("\n");
            sb.append("Gradient type: ").append(gradientType).append("\n");
            sb.append("Fx: ").append(minimumValue).append("\n");
            sb.append("LBFGS return: ").append(returnValue).append("\n");
            sb.append("Fisher refreshes: ").append(fisherRefreshCount).append("\n");
            if (fisherRefreshCount > 0) {
                sb.append("Fisher metric source: ")
                        .append(fisherUsedEmpiricalScores ? "empirical branch scores" : "central derivative of analytical score")
                        .append("\n");
                sb.append("Fisher rank/nullity: ").append(fisherRank).append("/").append(fisherNullity).append("\n");
                sb.append("Fisher positive eigenvalue range: [")
                        .append(fisherMinimumPositiveEigenvalue).append(", ")
                        .append(fisherMaximumEigenvalue).append("]\n");
            }
            sb.append("Time: ").append(time).append("s\n");
            sb.append("Count: ").append(count).append("\n");

        }

        return sb.toString();
    }

    private double[] minimize(LBFGS_Param paramsBFGS, double[] x0) {
        final int dimension = function.getDimension();
        final double[] point = x0 == null ? new double[dimension] : x0;
        final MutableDouble objective = new MutableDouble();
        final double[] fallbackPoint = point.clone();
        final double[][] bestFinitePoint = new double[1][];
        final double[] bestFiniteObjective = new double[] { Double.POSITIVE_INFINITY };
        final boolean[] sawInvalidEvaluation = new boolean[] { false };

        LBFGS_Evaluate evaluate = new LBFGS_Evaluate() {
            @Override
            public double eval(Object instance, double[] x, double[] g, int n, double step) {
                try {
                    final double value = function.valueAt(x);
                    if (!Double.isFinite(value)) {
                        throw new IllegalArgumentException("Non-finite objective in LBFGS evaluation.");
                    }

                    final double[] gradient = function.gradientAt(x);
                    if (!allFinite(gradient)) {
                        throw new IllegalArgumentException("Non-finite gradient in LBFGS evaluation.");
                    }

                    System.arraycopy(gradient, 0, g, 0, gradient.length);
                    if (value < bestFiniteObjective[0]) {
                        bestFiniteObjective[0] = value;
                        bestFinitePoint[0] = x.clone();
                    }
                    return value;
                } catch (IllegalArgumentException e) {
                    sawInvalidEvaluation[0] = true;
                    fillPenaltyGradient(x, g, bestFinitePoint[0] == null ? fallbackPoint : bestFinitePoint[0]);
                    return invalidObjective(x, bestFinitePoint[0] == null ? fallbackPoint : bestFinitePoint[0]);
                }
            }
        };

        LBFGS_Progress progress = settings.printToScreen ? new LBFGS_Progress() {
            @Override
            public ReturnValue eval(Object instance,
                                    double[] x,
                                    double[] g,
                                    double fx,
                                    double xnorm,
                                    double gnorm,
                                    double step,
                                    int n,
                                    int k,
                                    int ls) {
                System.out.printf("Iteration %d:\n", k);
                System.out.printf("\tfx = %f, xnorm = %f, gnorm = %f, step = %f\n", fx, xnorm, gnorm, step);
                System.out.printf("\tn = %d, k = %d, ls = %d\n\n", n, k, ls);
                return ReturnValue.LBFGS_SUCCESS;
            }
        } : null;

        returnValue = Lbfgs.lbfgs(dimension, point, objective, evaluate, progress, null, paramsBFGS);
        if (settings.printToScreen) {
            System.out.println("LBFGS return: " + returnValue);
        }
        return sawInvalidEvaluation[0] && bestFinitePoint[0] != null ? bestFinitePoint[0] : point;
    }

    private static boolean allFinite(final double[] values) {
        if (values == null) {
            return false;
        }
        for (final double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private boolean shouldRefresh(final ReturnValue value, final double gradientNorm) {
        if (!settings.fisherRefreshOnFailure || fisherRefreshCount >= settings.fisherMaximumRefreshes
                || !Double.isFinite(gradientNorm) || gradientNorm <= settings.fisherGradientThreshold) {
            return false;
        }
        return value == ReturnValue.LBFGSERR_ROUNDING_ERROR
                || value == ReturnValue.LBFGSERR_MINIMUMSTEP
                || value == ReturnValue.LBFGSERR_MAXIMUMLINESEARCH;
    }

    private boolean refreshFisherWhitening() {
        try {
            final Transform baseTransform = nativeTransform;
            final double[] parameterValues = parameter.getParameterValues();
            final double[] reference = baseTransform == null
                    ? parameterValues.clone()
                    : baseTransform.transform(parameterValues, 0, parameterValues.length);
            final double[][] fisher = numericalInformation(fisherGradient, baseTransform, reference);
            final double[][] prior = fisherPriorGradient == null
                    ? null
                    : numericalInformation(fisherPriorGradient, baseTransform, reference);
            final FisherWhiteningBuilder.Result result = FisherWhiteningBuilder.build(
                    fisher, prior, settings.fisherEigenFloorFraction,
                    settings.fisherNullPrecision, settings.fisherTrustRadius);
            final Transform.MultivariableTransform identity = new Transform.NoTransformMultivariable() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }
            };
            final Transform.MultivariableTransform inner = baseTransform == null
                    ? identity
                    : (Transform.MultivariableTransform) baseTransform;
            transform = new FixedWhiteningTransform(parameter, inner, reference, result.map);
            function = constructFunction();
            fisherRank = result.rank;
            fisherNullity = result.nullity;
            fisherMinimumPositiveEigenvalue = result.minimumPositiveEigenvalue;
            fisherMaximumEigenvalue = result.maximumEigenvalue;
            return true;
        } catch (RuntimeException exception) {
            if (settings.printToScreen) {
                System.err.println("Fisher refresh skipped: " + exception.getMessage());
            }
            return false;
        }
    }

    private double[][] numericalInformation(final GradientWrtParameterProvider provider,
                                            final Transform baseTransform,
                                            final double[] reference) {
        if (provider == null) {
            throw new IllegalArgumentException("Fisher refresh requires an analytical gradient provider.");
        }
        if (provider.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Fisher gradient dimension differs from optimizer dimension.");
        }
        if (provider instanceof ScoreContributionsWrtParameterProvider) {
            final double[][] information = empiricalInformation(
                    (ScoreContributionsWrtParameterProvider) provider,
                    baseTransform,
                    reference);
            if (information != null) {
                fisherUsedEmpiricalScores = true;
                return information;
            }
        }
        final int dimension = reference.length;
        final double[][] plus = new double[dimension][];
        final double[][] minus = new double[dimension][];
        final double[] endpoint = parameter.getParameterValues();
        for (int j = 0; j < dimension; ++j) {
            final double h = Math.sqrt(Math.sqrt(Math.ulp(1.0))) * (Math.abs(reference[j]) + 1.0);
            final double[] point = reference.clone();
            point[j] += h;
            plus[j] = transformedGradient(provider, baseTransform, point);
            point[j] -= 2.0 * h;
            minus[j] = transformedGradient(provider, baseTransform, point);
        }
        setParameter(new WrappedVector.Raw(endpoint), parameter);
        final double[][] information = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            final double hi = Math.sqrt(Math.sqrt(Math.ulp(1.0))) * (Math.abs(reference[i]) + 1.0);
            for (int j = 0; j < dimension; ++j) {
                final double hj = Math.sqrt(Math.sqrt(Math.ulp(1.0))) * (Math.abs(reference[j]) + 1.0);
                final double derivativeIJ = (plus[j][i] - minus[j][i]) / (2.0 * hj);
                final double derivativeJI = (plus[i][j] - minus[i][j]) / (2.0 * hi);
                information[i][j] = -0.5 * (derivativeIJ + derivativeJI);
            }
        }
        return information;
    }

    private double[][] empiricalInformation(final ScoreContributionsWrtParameterProvider provider,
                                            final Transform baseTransform,
                                            final double[] reference) {
        final double[] endpoint = parameter.getParameterValues();
        final double[][] rawScores = provider.getScoreContributionsLogDensity();
        if (rawScores == null) {
            return null;
        }
        if (rawScores.length == 0) {
            throw new IllegalArgumentException("Empirical Fisher provider returned an empty score array.");
        }
        final int dimension = reference.length;
        final double[][] information = new double[dimension][dimension];
        for (final double[] rawScore : rawScores) {
            if (rawScore == null || rawScore.length != dimension || !allFinite(rawScore)) {
                throw new IllegalArgumentException("Empirical Fisher score contribution is invalid.");
            }
            final double[] score = baseTransform == null
                    ? rawScore
                    : baseTransform.updateGradientUnWeightedLogDensity(
                    rawScore, endpoint, 0, endpoint.length);
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j <= i; ++j) {
                    information[i][j] += score[i] * score[j];
                }
            }
        }
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < i; ++j) {
                information[j][i] = information[i][j];
            }
        }
        return information;
    }

    private double[] transformedGradient(final GradientWrtParameterProvider provider,
                                         final Transform baseTransform,
                                         final double[] transformedPoint) {
        final double[] rawPoint = baseTransform == null
                ? transformedPoint
                : baseTransform.inverse(transformedPoint, 0, transformedPoint.length);
        setParameter(new WrappedVector.Raw(rawPoint), parameter);
        double[] result = provider.getGradientLogDensity();
        if (baseTransform != null) {
            result = baseTransform.updateGradientUnWeightedLogDensity(
                    result, rawPoint, 0, rawPoint.length);
        }
        if (!allFinite(result)) {
            throw new IllegalArgumentException("Non-finite gradient while building Fisher metric.");
        }
        return result;
    }

    private void setParameterFromOptimizationPoint(final double[] point) {
        final double[] rawPoint = transform == null
                ? point
                : transform.inverse(point, 0, point.length);
        setParameter(new WrappedVector.Raw(rawPoint), parameter);
    }

    private double objectiveAtParameter(final double[] point) {
        setParameter(new WrappedVector.Raw(point), parameter);
        return -evaluateLogLikelihood();
    }

    private static double euclideanNorm(final double[] values) {
        double sum = 0.0;
        for (final double value : values) {
            sum += value * value;
        }
        return Math.sqrt(sum);
    }

    private void validateFisherSettings() {
        if (!settings.fisherRefreshOnFailure) {
            return;
        }
        if (!(nativeTransform == null || nativeTransform instanceof Transform.MultivariableTransform)) {
            throw new IllegalArgumentException("Fisher refresh requires a multivariable transform.");
        }
        if (settings.includeJacobian) {
            throw new IllegalArgumentException(
                    "Fisher refresh cannot be combined with includeJacobian because recentering changes its constant term.");
        }
        if (fisherGradient == null || fisherGradient.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Fisher gradient must match the optimizer parameter dimension.");
        }
        if (fisherPriorGradient != null
                && fisherPriorGradient.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Fisher prior gradient must match the optimizer parameter dimension.");
        }
        if (settings.fisherMaximumRefreshes < 1) {
            throw new IllegalArgumentException("Fisher maximum refreshes must be at least one.");
        }
        if (!(settings.fisherGradientThreshold >= 0.0) || !Double.isFinite(settings.fisherGradientThreshold)) {
            throw new IllegalArgumentException("Fisher gradient threshold must be finite and non-negative.");
        }
    }

    private static double invalidObjective(final double[] point, final double[] referencePoint) {
        double squaredDistance = 0.0;
        for (int i = 0; i < point.length; ++i) {
            final double diff = finiteDifference(point[i], referencePoint[i]);
            squaredDistance += diff * diff;
        }
        return 1.0e100 + 0.5 * squaredDistance;
    }

    private static void fillPenaltyGradient(final double[] point,
                                            final double[] gradient,
                                            final double[] referencePoint) {
        for (int i = 0; i < gradient.length; ++i) {
            gradient[i] = finiteDifference(point[i], referencePoint[i]);
        }
    }

    private static double finiteDifference(final double value, final double reference) {
        if (!Double.isFinite(value) || !Double.isFinite(reference)) {
            return 0.0;
        }
        final double difference = value - reference;
        if (difference > 1.0e6) {
            return 1.0e6;
        }
        if (difference < -1.0e6) {
            return -1.0e6;
        }
        return difference;
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

                if (settings.includeJacobian) {
                    assert transform != null;
                    return -evaluateLogLikelihood() - transform.logJacobian(argument, 0, argument.length);
                }
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

                    if (settings.includeJacobian) {
                        result = transform.updateGradientLogDensity(result, argument, 0, argument.length);
                    } else {
                        result = transform.updateGradientUnWeightedLogDensity(result, argument, 0, argument.length);
                    }


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
                double[] gradient = new double[parameter.getDimension()];
                NumericalDerivative.gradient(function, parameter.getParameterValues(), gradient, settings.stepSize);
                return gradient;
            }
        };
    }

    private enum GradientType {
        ANALYTIC("analytic"),
        NUMERICAL("numerical");

        private final String type;

        GradientType(String type) {
            this.type = type;
        }

        public String toString() { return type; }
    }
}
