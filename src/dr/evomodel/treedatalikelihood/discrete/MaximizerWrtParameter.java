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
import dr.inference.model.Bounds;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Profileable;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Evaluate;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Progress;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue;
import com.github.lbfgs4j.liblbfgs.MutableDouble;

import java.util.Arrays;
import java.util.logging.Logger;

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
    private double[] minimumGradient = null;
    private boolean executed = false;

    // the last (y, f, g) evaluated during maximize(), where L-BFGS asks for the value and the gradient at the same
    // point; outside maximize() other models may change between calls, so every call recomputes
    private boolean useCache = false;
    private double[] cachedPoint = null;
    private double cachedValue = Double.NaN;
    private double[] cachedGradient = null;
    private double[] parameterPoint = null;

    // the latest run of maximize()
    private ReturnValue status = null;
    private Reason reason = null;
    private int iterations = 0;
    private long valueCalls = 0;
    private long gradientCalls = 0;
    private long nonFiniteEvaluations = 0;
    private long likelihoodRecomputes = 0;
    private double seconds = 0.0;
    private Reason stopRequest = null;
    private long startNanos = 0;
    private double[] stallHistory = null;
    private double worstGradientCheckError = Double.NaN;

    private double lastFiniteValue = 0.0;
    private boolean lastEvaluationNonFinite = false;
    private static final double NON_FINITE_PENALTY = 1E10;

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel");

    public static class Settings {
        int numberIterations;
        double stepSize;
        boolean startAtCurrentState;
        boolean printToScreen;
        boolean includeJacobian;

        // optional budgets and lbfgs4j settings; the defaults keep the behaviour of earlier versions
        double maxSeconds = 0.0;
        long maxEvaluations = 0;
        double relTol = 0.0;
        int stallWindow = 0;
        int checkGradientEvery = 0;
        Double epsilon = null;
        Integer past = null;
        Double delta = null;

        /**
         * @param includeJacobian false (the default, used for optimization in TreeLime): maximize the unweighted log
         *                        density, whose optimum does not depend on the transform; true: maximize the density
         *                        of the transformed variable, log L - logJacobian, as HMC samples it
         */
        public Settings(int numberIterations, double stepSize, boolean startAtCurrentState, boolean printToScreen,
                        boolean includeJacobian) {
            this.numberIterations = numberIterations;
            this.stepSize = stepSize;
            this.startAtCurrentState = startAtCurrentState;
            this.printToScreen = printToScreen;
            this.includeJacobian = includeJacobian;
        }

        /**
         * Stop after maxSeconds of wall-clock time or maxEvaluations objective evaluations (0: no limit), or when the
         * objective has not improved by more than relTol (relative) over the last stallWindow iterations (0: off).
         */
        public Settings setBudgets(double maxSeconds, long maxEvaluations, double relTol, int stallWindow) {
            this.maxSeconds = maxSeconds;
            this.maxEvaluations = maxEvaluations;
            this.relTol = relTol;
            this.stallWindow = stallWindow;
            return this;
        }

        /**
         * Every checkGradientEvery gradient evaluations (0: never), compare the gradient along a random direction with
         * central finite differences of the objective and log the relative error.
         */
        public Settings setCheckGradientEvery(int checkGradientEvery) {
            this.checkGradientEvery = checkGradientEvery;
            return this;
        }

        /**
         * lbfgs4j's convergence tests: epsilon (gradient norm relative to max(1, |x|)), and past/delta (relative
         * decrease of the objective over the last past iterations); null keeps lbfgs4j's default.
         */
        public Settings setConvergence(Double epsilon, Integer past, Double delta) {
            this.epsilon = epsilon;
            this.past = past;
            this.delta = delta;
            return this;
        }
    }

    /**
     * Why the latest run of maximize() stopped.
     */
    public enum Reason {
        CONVERGED("converged"),
        MAX_ITERATIONS("maxIterations"),
        MAX_SECONDS("maxSeconds"),
        MAX_EVALUATIONS("maxEvaluations"),
        STALL("stall"),
        LINE_SEARCH_FAILURE("lineSearchFailure"),
        NON_FINITE("nonFinite"),
        ERROR("error");

        private final String name;

        Reason(String name) {
            this.name = name;
        }

        public String toString() { return name; }
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
        return executed;
    }

    public Transform getTransform() {
        return transform;
    }

    public Function getFunction() {
        return function;
    }

    public ReturnValue getStatus() {
        return status;
    }

    public Reason getReason() {
        return reason;
    }

    public int getIterations() {
        return iterations;
    }

    public long getValueCalls() {
        return valueCalls;
    }

    public long getGradientCalls() {
        return gradientCalls;
    }

    public long getNonFiniteEvaluations() {
        return nonFiniteEvaluations;
    }

    /**
     * Warnings about the set-up that would otherwise surface only as a run that never stops or misbehaves.
     */
    public String[] getWarnings() {
        java.util.List<String> warnings = new java.util.ArrayList<>();
        if (settings.numberIterations <= 0 && settings.maxSeconds <= 0 && settings.maxEvaluations <= 0
                && settings.stallWindow <= 0 && settings.past == null) {
            warnings.add("nIterations is 0 or absent and no budget is set: the maximizer runs until the gradient " +
                    "test passes, which may never happen");
        }
        if (transform == null) {
            Bounds<Double> bounds;
            try {
                bounds = parameter.getBounds();
            } catch (NullPointerException e) { // Parameter.Default without bounds
                bounds = null;
            }
            for (int i = 0; i < parameter.getDimension(); ++i) {
                if (bounds != null && (bounds.getLowerLimit(i) > Double.NEGATIVE_INFINITY
                        || bounds.getUpperLimit(i) < Double.POSITIVE_INFINITY)) {
                    warnings.add("parameter " + parameterName() + " is bounded but has no transform: " +
                            "steps can leave its domain");
                    break;
                }
            }
        }
        if (!settings.startAtCurrentState) {
            warnings.add("startAtCurrentState is false: the search starts at 0 in " +
                    (transform == null ? "parameter" : "transformed") + " space");
        }
        return warnings.toArray(new String[0]);
    }

    public void maximize() {

        LBFGS_Param paramsBFGS = Lbfgs.defaultParams();

        if (settings.numberIterations > 0) {
            paramsBFGS.max_iterations = settings.numberIterations;
        }
        if (settings.epsilon != null) {
            paramsBFGS.epsilon = settings.epsilon;
        }
        if (settings.past != null) {
            paramsBFGS.past = settings.past;
        }
        if (settings.delta != null) {
            paramsBFGS.delta = settings.delta;
        }

        final int dim = function.getDimension();
        double[] x;

        if (settings.startAtCurrentState) {
            x = parameter.getParameterValues();

            if (transform != null) {
                x = transform.transform(x, 0, x.length);
            }
        } else {
            x = new double[dim];
        }

        // operators call maximize() repeatedly: warn on the first run only
        final boolean firstRun = !executed;
        if (settings.printToScreen) {
            LOGGER.info(getBanner());
        }
        if (firstRun) {
            for (String warning : getWarnings()) {
                LOGGER.info("WARNING: maximizeWrtParameter: " + warning);
            }
        }

        status = null;
        reason = null;
        iterations = 0;
        valueCalls = 0;
        gradientCalls = 0;
        nonFiniteEvaluations = 0;
        stopRequest = null;
        stallHistory = settings.stallWindow > 0 ? new double[settings.stallWindow + 1] : null;
        worstGradientCheckError = Double.NaN;
        final long recomputesAtStart = getLikelihoodCalculationCount();
        cachedPoint = null;
        useCache = true;
        lastFiniteValue = 0.0;
        lastEvaluationNonFinite = false;

        startNanos = System.nanoTime();

        final LBFGS_Evaluate evaluate = (instance, point, g, n, step) -> {
            final double value = evaluate(point);
            System.arraycopy(cachedGradient, 0, g, 0, n);
            return value;
        };

        final LBFGS_Progress progress = (instance, point, g, fx, xnorm, gnorm, step, n, k, ls) -> {
            iterations = k;
            if (settings.printToScreen) {
                System.out.printf("Iteration %d:\n", k);
                System.out.printf("\tfx = %f, xnorm = %f, gnorm = %f, step = %f\n", fx, xnorm, gnorm, step);
                System.out.printf("\tn = %d, k = %d, ls = %d\n", n, k, ls);
                System.out.printf("\tevaluations = %d, seconds = %.3f\n\n", valueCalls, elapsedSeconds());
            }
            stopRequest = checkBudgets(fx, k);
            return stopRequest == null ? ReturnValue.LBFGS_SUCCESS : ReturnValue.LBFGS_STOP;
        };

        final MutableDouble fx = new MutableDouble();
        status = Lbfgs.lbfgs(dim, x, fx, evaluate, progress, null, paramsBFGS);

        seconds = elapsedSeconds();
        time += seconds;
        reason = stopRequest != null ? stopRequest : reasonFor(status);

        minimumPoint = x;
        minimumValue = function.valueAt(minimumPoint);
        minimumGradient = function.gradientAt(minimumPoint).clone();
        likelihoodRecomputes = getLikelihoodCalculationCount() - recomputesAtStart;
        useCache = false;
        cachedPoint = null;
        executed = true;

        if (transform != null) {
            setParameter(new WrappedVector.Raw(transform.inverse(minimumPoint, 0, minimumPoint.length)), parameter);
        } else {
            setParameter(new WrappedVector.Raw(minimumPoint), parameter);
        }

        final boolean failed = reason == Reason.LINE_SEARCH_FAILURE || reason == Reason.NON_FINITE
                || reason == Reason.ERROR;
        if (settings.printToScreen || (failed && firstRun)) {
            LOGGER.info("maximizeWrtParameter " + parameterName() + " stopped: " + reason + " (" + status +
                    ") after " + iterations + " iterations, " + valueCalls + " evaluations, " +
                    String.format("%.3f", seconds) + "s; -log density = " + minimumValue);
        }
    }

    private String getBanner() {
        StringBuilder sb = new StringBuilder("Maximizing wrt ").append(parameterName())
                .append(" (dim ").append(function.getDimension())
                .append(", ").append(gradientType).append(" gradient")
                .append(", nIterations = ").append(settings.numberIterations > 0 ? settings.numberIterations : "unlimited");
        if (settings.maxSeconds > 0) sb.append(", maxSeconds = ").append(settings.maxSeconds);
        if (settings.maxEvaluations > 0) sb.append(", maxEvaluations = ").append(settings.maxEvaluations);
        if (settings.stallWindow > 0) {
            sb.append(", stall: relTol = ").append(settings.relTol).append(" over ").append(settings.stallWindow);
        }
        return sb.append(")").toString();
    }

    private String parameterName() {
        return parameter.getId() != null ? parameter.getId() : parameter.getParameterName();
    }

    private double elapsedSeconds() {
        return (System.nanoTime() - startNanos) / 1E9;
    }

    private Reason checkBudgets(double fx, int k) {
        if (settings.maxSeconds > 0 && elapsedSeconds() >= settings.maxSeconds) {
            return Reason.MAX_SECONDS;
        }
        if (settings.maxEvaluations > 0 && valueCalls >= settings.maxEvaluations) {
            return Reason.MAX_EVALUATIONS;
        }
        if (stallHistory != null) {
            final int window = settings.stallWindow;
            stallHistory[k % (window + 1)] = fx;
            if (k > window) {
                final double before = stallHistory[(k - window) % (window + 1)];
                if (before - fx <= settings.relTol * Math.max(1.0, Math.abs(fx))) {
                    return Reason.STALL;
                }
            }
        }
        return null;
    }

    private Reason reasonFor(ReturnValue value) {
        switch (value) {
            case LBFGS_SUCCESS:
            case LBFGS_CONVERGENCE:
            case LBFGS_STOP: // the past/delta test
            case LBFGS_ALREADY_MINIMIZED:
                return Reason.CONVERGED;
            case LBFGSERR_MAXIMUMITERATION:
                return Reason.MAX_ITERATIONS;
            case LBFGSERR_ROUNDING_ERROR:
            case LBFGSERR_MINIMUMSTEP:
            case LBFGSERR_MAXIMUMSTEP:
            case LBFGSERR_MAXIMUMLINESEARCH:
            case LBFGSERR_WIDTHTOOSMALL:
            case LBFGSERR_INCREASEGRADIENT:
            case LBFGSERR_OUTOFINTERVAL:
            case LBFGSERR_INCORRECT_TMINMAX:
            case LBFGSERR_INVALIDPARAMETERS:
                return lastEvaluationNonFinite ? Reason.NON_FINITE : Reason.LINE_SEARCH_FAILURE;
            default:
                return Reason.ERROR;
        }
    }

    private long getLikelihoodCalculationCount() {
        return likelihood instanceof Profileable ? ((Profileable) likelihood).getTotalCalculationCount() : 0;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        if (!executed) {
            sb.append("Not yet executed.");
        } else {

            if (transform != null) {
                sb.append("Gradient is taken with respect to the transformed parameter values.\n");
                sb.append("Untransformed X: ").append(new dr.math.matrixAlgebra.Vector(transform.inverse(minimumPoint, 0, minimumPoint.length))).append("\n");
            }
            sb.append("X: ").append(new dr.math.matrixAlgebra.Vector(minimumPoint)).append("\n");
            sb.append("Gradient: ").append(new dr.math.matrixAlgebra.Vector(minimumGradient)).append("\n");
            sb.append("Gradient type: ").append(gradientType).append("\n");
            sb.append("Fx: ").append(minimumValue).append("\n");
            sb.append("Time: ").append(time).append("s\n");
            // calls through evaluateLogLikelihood, over every run
            sb.append("Count: ").append(count).append("\n");
            // the latest run
            sb.append("Status: ").append(status).append("\n");
            sb.append("Reason: ").append(reason).append("\n");
            sb.append("Iterations: ").append(iterations).append("\n");
            sb.append("ValueCalls: ").append(valueCalls).append("\n");
            sb.append("GradientCalls: ").append(gradientCalls).append("\n");
            sb.append("LikelihoodRecomputes: ").append(likelihoodRecomputes).append("\n");
            sb.append("NonFiniteEvaluations: ").append(nonFiniteEvaluations).append("\n");
            sb.append("Seconds: ").append(seconds).append("\n");
            final double gradientNorm = norm(minimumGradient);
            sb.append("GradientNorm: ").append(gradientNorm)
                    .append(" (relative ").append(gradientNorm / Math.max(1.0, norm(minimumPoint))).append(")\n");
            if (settings.checkGradientEvery > 0) {
                sb.append("WorstGradientCheckError: ").append(worstGradientCheckError).append("\n");
            }
        }

        return sb.toString();
    }

    private static double norm(double[] x) {
        double sum = 0.0;
        for (double v : x) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    private double evaluateLogLikelihood() {
        ++count;
        return likelihood.getLogLikelihood();
    }

    private void setPoint(double[] point) {
        parameterPoint = point.clone();
        double[] argument = transform != null ? transform.inverse(point, 0, point.length) : point;
        setParameter(new WrappedVector.Raw(argument), parameter);
    }

    /**
     * The objective (-log density) and its gradient at the transformed point y, computed together and cached. A
     * non-finite log density returns a large penalty at once, without computing its gradient, so the line search
     * backs off.
     */
    private double evaluate(double[] point) {

        if (useCache && cachedPoint != null && Arrays.equals(point, cachedPoint)) {
            if (!Arrays.equals(point, parameterPoint)) {
                setPoint(point);
            }
            return cachedValue;
        }

        setPoint(point);
        final double[] argument = parameter.getParameterValues();

        ++valueCalls;
        final double logLikelihood = evaluateLogLikelihood();

        double value;
        double[] result;

        lastEvaluationNonFinite = !Double.isFinite(logLikelihood);
        if (lastEvaluationNonFinite) {
            ++nonFiniteEvaluations;
            // lbfgs4j's default (More-Thuente) line search interpolates the value, and +Inf turns its next step into
            // NaN; a large finite value with a zero gradient makes it step back towards the last finite point instead
            value = NON_FINITE_PENALTY * (1.0 + Math.abs(lastFiniteValue));
            result = new double[point.length];
        } else {
            value = -logLikelihood;
            if (settings.includeJacobian) {
                assert transform != null;
                // the transformed density, log L - logJacobian (HamiltonianMonteCarloOperator's convention), whose
                // gradient is updateGradientLogDensity below
                value += transform.logJacobian(argument, 0, argument.length);
            }
            // the unweighted log density otherwise: the transform only reparameterizes, so the optimum is unchanged

            ++gradientCalls;
            result = gradient.getGradientLogDensity();

            if (transform != null) {
                if (settings.includeJacobian) {
                    result = transform.updateGradientLogDensity(result, argument, 0, argument.length);
                } else {
                    // chain rule only, matching the unweighted log density
                    result = transform.updateGradientUnWeightedLogDensity(result, argument, 0, argument.length);
                }
            }
            result = result.clone();
            for (int i = 0; i < result.length; ++i) {
                result[i] = -result[i];
            }
            lastFiniteValue = value;
        }

        cachedPoint = point.clone();
        cachedValue = value;
        cachedGradient = result;

        if (settings.checkGradientEvery > 0 && Double.isFinite(value)
                && gradientCalls % settings.checkGradientEvery == 0) {
            checkGradient(point.clone(), value, result.clone());
        }

        return value;
    }

    /**
     * Compare g.d with (f(y + hd) - f(y - hd)) / 2h along a random unit direction d, log the relative error and
     * restore the cache and the parameter at y.
     */
    private void checkGradient(double[] point, double value, double[] grad) {

        final double[] direction = new double[point.length];
        for (int i = 0; i < direction.length; ++i) {
            direction[i] = MathUtils.nextGaussian();
        }
        final double length = norm(direction);
        double directional = 0.0;
        for (int i = 0; i < direction.length; ++i) {
            direction[i] /= length;
            directional += grad[i] * direction[i];
        }

        final double h = 1E-5 * Math.max(1.0, norm(point));
        final double[] plus = point.clone();
        final double[] minus = point.clone();
        for (int i = 0; i < point.length; ++i) {
            plus[i] += h * direction[i];
            minus[i] -= h * direction[i];
        }
        final double numeric = (objectiveOnly(plus) - objectiveOnly(minus)) / (2 * h);
        final double error = Math.abs(numeric - directional) / Math.max(1.0, Math.abs(numeric));

        if (Double.isFinite(error) && !(error <= worstGradientCheckError)) {
            worstGradientCheckError = error;
        }
        LOGGER.info("maximizeWrtParameter gradient check: directional derivative " + directional +
                ", finite difference " + numeric + ", relative error " + error);

        setPoint(point);
        cachedPoint = point;
        cachedValue = value;
        cachedGradient = grad;
    }

    private double objectiveOnly(double[] point) {
        setPoint(point);
        final double[] argument = parameter.getParameterValues();
        double value = -evaluateLogLikelihood();
        if (settings.includeJacobian) {
            value += transform.logJacobian(argument, 0, argument.length);
        }
        return value;
    }

    private Function constructFunction() {

        return new Function() {

            @Override
            public int getDimension() {
                return gradient.getDimension();
            }

            @Override
            public double valueAt(double[] argument) {
                // outside maximize() nobody asks for the gradient at the same point: skip it
                return useCache ? evaluate(argument) : objectiveOnly(argument);
            }

            @Override
            public double[] gradientAt(double[] argument) {
                evaluate(argument);
                return cachedGradient.clone();
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
                final double[] x = parameter.getParameterValues();
                NumericalDerivative.gradient(function, x, gradient, settings.stepSize);
                // NumericalDerivative leaves the last coordinate at x - h
                setParameter(new WrappedVector.Raw(x), parameter);
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
