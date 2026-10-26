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
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Profileable;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.util.UnitSimplexToRealsTransform;
import dr.xml.Reportable;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Evaluate;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Progress;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue;
import com.github.lbfgs4j.liblbfgs.MutableDouble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

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
    private boolean cachedPointApplied = true; // false when the point's untransformed values were not finite

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
    private int restarts = 0;
    private boolean convergedToPrecision = false; // CONVERGED after a line search failure no restart could improve on
    private double valueAtLastRestart = Double.POSITIVE_INFINITY;
    private final List<Double> gradientCheckErrors = new ArrayList<Double>();

    // coordinates held at a lower bound (holdAtLowerBound) in the latest run, by index; null when none is set
    private boolean[] heldCoordinates = null;
    private double[] heldValues = null; // the bound of each held element, set exactly
    private double[] cachedUnmaskedGradient = null;
    private double[] holdPoint = null;
    private int holds = 0;

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
        double maxRunTime = 0.0; // seconds
        long maxEvaluations = 0;
        double relTol = 0.0;
        int stallWindow = 0;
        int checkGradientEvery = 0;
        int lineSearchRestarts = 0;
        Double epsilon = null;
        Integer past = null;
        Double delta = null;
        final List<LowerBoundHold> lowerBoundHolds = new ArrayList<LowerBoundHold>();

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
         * Stop after maxRunTime seconds of wall-clock time or maxEvaluations objective evaluations (0: no limit), or when the
         * objective has not improved by more than relTol (relative) over the last stallWindow iterations (0: off).
         */
        public Settings setBudgets(double maxRunTime, long maxEvaluations, double relTol, int stallWindow) {
            this.maxRunTime = maxRunTime;
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
         * After a failed line search, restart L-BFGS from the best point so far with a cleared memory, up to this
         * many times, while each restart makes progress (0: stop at the failure). With delta set (setConvergence), a
         * restart that finds no decrease at all ends the run as converged, to numerical precision.
         */
        public Settings setLineSearchRestarts(int lineSearchRestarts) {
            this.lineSearchRestarts = lineSearchRestarts;
            return this;
        }

        /**
         * Once an accepted iterate takes an element of the parameter (part of the maximized parameter) below the
         * bound, set that element to the bound and keep it there for the rest of the run: its gradient coordinate is
         * zeroed and L-BFGS restarts with a cleared memory. For a parameter whose optimum lies on its boundary, such
         * as a gamma shape that runs to 0, where L-BFGS would otherwise chase a plateau until its line search fails.
         * The parameter's transform must act element by element (log, logit, none).
         */
        public Settings addLowerBoundHold(Parameter parameter, double bound) {
            lowerBoundHolds.add(new LowerBoundHold(parameter, bound));
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

    private static final class LowerBoundHold {
        final Parameter parameter;
        final double bound;
        int offset = -1; // of the parameter's first element in the maximized parameter

        LowerBoundHold(Parameter parameter, double bound) {
            this.parameter = parameter;
            this.bound = bound;
        }

        String name(int element) {
            String name = parameter.getId() != null ? parameter.getId() : parameter.getParameterName();
            return parameter.getDimension() > 1 ? name + "[" + element + "]" : name;
        }
    }

    /**
     * An element held at its lower bound at the end of the latest run, with the derivative of the log density along
     * its transformed coordinate there: negative when the density still increases towards the boundary, which
     * confirms an optimum on the boundary; positive when it points back inside, where the bound may be too high.
     */
    public static final class Held {
        public final String name;
        public final double bound;
        public final double gradient;

        Held(String name, double bound, double gradient) {
            this.name = name;
            this.bound = bound;
            this.gradient = gradient;
        }

        public String toString() {
            return name + " = " + bound + " (dlogDensity/dy = " + gradient + (gradient <= 0 ?
                    ": the optimum lies on the boundary)" : ": points back inside, the bound may be too high)");
        }
    }

    /**
     * Why the latest run of maximize() stopped.
     */
    public enum Reason {
        CONVERGED("converged"),
        MAX_ITERATIONS("maxIterations"),
        MAX_RUN_TIME("maxRunTime"),
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

        this.gradient = gradient == null ? constructGradient() : gradient;

        this.function = constructFunction();
        this.settings = settings;

        for (LowerBoundHold hold : settings.lowerBoundHolds) {
            hold.offset = locate(this.parameter, hold.parameter);
            if (hold.offset < 0) {
                throw new IllegalArgumentException("holdAtLowerBound: " + hold.name(0) + " is not part of " +
                        parameterName());
            }
        }
    }

    // the offset of target's first element in root, or -1
    private static int locate(Parameter root, Parameter target) {
        if (root == target) {
            return 0;
        }
        if (root instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) root;
            int offset = 0;
            for (int i = 0; i < compound.getParameterCount(); ++i) {
                Parameter child = compound.getParameter(i);
                int at = locate(child, target);
                if (at >= 0) {
                    return offset + at;
                }
                offset += child.getDimension();
            }
        }
        return -1;
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

    public int getRestarts() {
        return restarts;
    }

    /**
     * @return the elements held at their lower bound at the end of the latest run
     */
    public List<Held> getHeld() {
        List<Held> held = new ArrayList<Held>();
        if (heldCoordinates == null || minimumGradient == null) {
            return held;
        }
        for (LowerBoundHold hold : settings.lowerBoundHolds) {
            for (int j = 0; j < hold.parameter.getDimension(); ++j) {
                int index = hold.offset + j;
                if (heldCoordinates[index]) {
                    held.add(new Held(hold.name(j), hold.bound, -heldGradient[index]));
                }
            }
        }
        return held;
    }

    private double[] heldGradient = null; // the objective's unmasked gradient at the optimum

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
     * @return the relative errors of the latest run's gradient checks (checkGradientEvery), in order; finite
     * differences are unreliable where the objective barely changes, so single large values need not mean a wrong
     * gradient
     */
    public List<Double> getGradientCheckErrors() {
        return new ArrayList<Double>(gradientCheckErrors);
    }

    public Parameter getParameter() {
        return parameter;
    }

    /**
     * @return the objective at the optimum of the latest run: -log density (plus the log-Jacobian with includeJacobian)
     */
    public double getMinimumValue() {
        return minimumValue;
    }

    public boolean includesJacobian() {
        return settings.includeJacobian;
    }

    /**
     * @return the wall-clock seconds of the latest run
     */
    public double getSeconds() {
        return seconds;
    }

    public double getWorstGradientCheckError() {
        return worstGradientCheckError;
    }

    /**
     * @return the settings that decide when a run stops, as name = value pairs
     */
    public String describeSettings() {
        StringBuilder sb = new StringBuilder();
        sb.append("nIterations = ").append(settings.numberIterations);
        if (settings.maxRunTime > 0) {
            sb.append(", maxRunTime = ").append(settings.maxRunTime).append("s");
        }
        if (settings.maxEvaluations > 0) {
            sb.append(", maxEvaluations = ").append(settings.maxEvaluations);
        }
        if (settings.stallWindow > 0) {
            sb.append(", relTol = ").append(settings.relTol).append(", stallWindow = ").append(settings.stallWindow);
        }
        if (settings.epsilon != null) {
            sb.append(", epsilon = ").append(settings.epsilon);
        }
        if (settings.past != null) {
            sb.append(", past = ").append(settings.past);
        }
        if (settings.delta != null) {
            sb.append(", delta = ").append(settings.delta);
        }
        sb.append(", lineSearchRestarts = ").append(settings.lineSearchRestarts);
        if (settings.checkGradientEvery > 0) {
            sb.append(", checkGradientEvery = ").append(settings.checkGradientEvery);
        }
        sb.append(", includeJacobian = ").append(settings.includeJacobian);
        return sb.toString();
    }

    /**
     * One coordinate of the optimum. The transformed value is the coordinate L-BFGS moves; with a multivariate
     * transform (node-height ratios, the stick-breaking simplex) it belongs to the block, not to the element of the
     * same name.
     */
    public static final class Coordinate {
        public final String name;
        public final double value; // natural scale
        public final double transformed;
        public final String transform;
        public final boolean free; // false for the stick-breaking simplex's last coordinate, which is always 0
        public final boolean held; // at its lower bound (holdAtLowerBound)
        public final double gradient; // dlogDensity/d(transformed), before held coordinates are masked

        Coordinate(String name, double value, double transformed, String transform, boolean free, boolean held,
                   double gradient) {
            this.name = name;
            this.value = value;
            this.transformed = transformed;
            this.transform = transform;
            this.free = free;
            this.held = held;
            this.gradient = gradient;
        }

        /**
         * @return optimized, held_at_bound or dependent
         */
        public String getStatus() {
            return held ? "held_at_bound" : (free ? "optimized" : "dependent");
        }
    }

    /**
     * @return the coordinates of the latest run's optimum, named by the maximized parameter's dimension names, or an
     * empty list before a run
     */
    public List<Coordinate> getOptimum() {
        List<Coordinate> coordinates = new ArrayList<Coordinate>();
        if (!executed) {
            return coordinates;
        }
        final double[] values = modelValues(minimumPoint);
        final String[] transforms = new String[minimumPoint.length];
        final boolean[] free = new boolean[minimumPoint.length];
        Arrays.fill(free, true);
        describeTransform(transform, 0, minimumPoint.length, transforms, free);
        final double[] grad = heldGradient != null ? heldGradient : minimumGradient;
        for (int i = 0; i < minimumPoint.length; ++i) {
            coordinates.add(new Coordinate(parameter.getDimensionName(i), values[i], minimumPoint[i], transforms[i],
                    free[i], heldCoordinates != null && heldCoordinates[i], -grad[i]));
        }
        return coordinates;
    }

    /**
     * @return the number of free coordinates (k), held ones included
     */
    public int getFreeCoordinateCount() {
        int k = 0;
        for (Coordinate coordinate : getOptimum()) {
            if (coordinate.free) {
                ++k;
            }
        }
        return k;
    }

    // the name of the transform acting on each coordinate in [from, to)
    private static void describeTransform(Transform transform, int from, int to, String[] names, boolean[] free) {
        if (transform == null) {
            Arrays.fill(names, from, to, "none");
        } else if (transform instanceof Transform.MultivariateArray) {
            int start = from;
            for (Transform.MultivariableTransform part : ((Transform.MultivariateArray) transform).getTransforms()) {
                int end = Math.min(to, start + part.getDimension());
                describeTransform(part, start, end, names, free);
                start = end;
            }
        } else if (transform instanceof Transform.Collection) {
            for (Transform.ParsedTransform segment : ((Transform.Collection) transform).getSegments()) {
                describeTransform(segment.transform, from + segment.start, Math.min(to, from + segment.end), names,
                        free);
            }
        } else if (transform instanceof Transform.Array) {
            for (int i = from; i < to; ++i) {
                describeTransform(((Transform.Array) transform).getTransform(i - from), i, i + 1, names, free);
            }
        } else {
            Arrays.fill(names, from, to, transform.getTransformName());
            if (transform instanceof UnitSimplexToRealsTransform && to > from) {
                free[to - 1] = false;
            }
        }
    }

    /**
     * Warnings about the set-up that would otherwise surface only as a run that never stops or misbehaves.
     */
    public String[] getWarnings() {
        java.util.List<String> warnings = new java.util.ArrayList<>();
        if (settings.numberIterations <= 0 && settings.maxRunTime <= 0 && settings.maxEvaluations <= 0
                && settings.stallWindow <= 0 && settings.past == null) {
            warnings.add("nIterations is 0 or absent and no budget is set: the maximizer runs until the gradient " +
                    "test passes, which may never happen");
        }
        if (transform == null && isBounded(parameter)) {
            warnings.add("parameter " + parameterName() + " is bounded but has no transform: " +
                    "steps can leave its domain");
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
        gradientCheckErrors.clear();
        final long recomputesAtStart = getLikelihoodCalculationCount();
        cachedPoint = null;
        useCache = true;
        lastFiniteValue = 0.0;
        lastEvaluationNonFinite = false;
        heldCoordinates = settings.lowerBoundHolds.isEmpty() ? null : new boolean[dim];
        heldValues = settings.lowerBoundHolds.isEmpty() ? null : new double[dim];
        holdPoint = null;
        holds = 0;
        if (heldCoordinates != null) {
            // a start below a bound is held at once
            double[] held = holdAtBounds(x);
            if (held != null) {
                x = held;
            }
        }

        startNanos = System.nanoTime();

        final LBFGS_Evaluate evaluate = (instance, point, g, n, step) -> {
            final double value = evaluate(point);
            System.arraycopy(cachedGradient, 0, g, 0, n);
            return value;
        };

        // iterations before the current L-BFGS run, and the best accepted point, for restarts
        final int[] previousIterations = {0};
        final double[] bestValue = {Double.POSITIVE_INFINITY};
        final double[][] bestPoint = {x.clone()};

        final LBFGS_Progress progress = (instance, point, g, fx, xnorm, gnorm, step, n, k, ls) -> {
            iterations = previousIterations[0] + k;
            if (fx < bestValue[0]) {
                bestValue[0] = fx;
                bestPoint[0] = point.clone();
            }
            if (settings.printToScreen) {
                System.out.printf("Iteration %d:\n", iterations);
                System.out.printf("\tfx = %f, xnorm = %f, gnorm = %f, step = %f\n", fx, xnorm, gnorm, step);
                System.out.printf("\tn = %d, k = %d, ls = %d\n", n, k, ls);
                System.out.printf("\tevaluations = %d, seconds = %.3f\n\n", valueCalls, elapsedSeconds());
            }
            stopRequest = checkBudgets(fx, iterations);
            if (stopRequest == null && heldCoordinates != null && belowBound(point)) {
                holdPoint = point.clone();
            }
            return stopRequest == null && holdPoint == null ? ReturnValue.LBFGS_SUCCESS : ReturnValue.LBFGS_STOP;
        };

        final MutableDouble fx = new MutableDouble();
        restarts = 0;
        while (true) {
            status = Lbfgs.lbfgs(dim, x, fx, evaluate, progress, null, paramsBFGS);
            reason = stopRequest != null ? stopRequest : reasonFor(status);

            if (holdPoint != null) {
                // an element went below its bound: hold it there and restart with a cleared memory
                x = holdAtBounds(holdPoint);
                holdPoint = null;
                holds++;
                previousIterations[0] = iterations;
                bestValue[0] = Double.POSITIVE_INFINITY;
                bestPoint[0] = x.clone();
                valueAtLastRestart = Double.POSITIVE_INFINITY;
                if (settings.numberIterations > 0) {
                    if (settings.numberIterations - iterations <= 0) {
                        reason = Reason.MAX_ITERATIONS;
                        break;
                    }
                    paramsBFGS.max_iterations = settings.numberIterations - iterations;
                }
                if (settings.printToScreen) {
                    LOGGER.info("maximizeWrtParameter: holding " + describeHeld() + " at the lower bound; restart");
                }
                continue;
            }

            // a failed line search ends the run; restarting from the best point with a cleared memory (a steepest-
            // descent first step) often gets past it, as long as each restart has made progress
            final int remaining = settings.numberIterations > 0 ? settings.numberIterations - iterations : 1;
            if (reason == Reason.LINE_SEARCH_FAILURE && restarts < settings.lineSearchRestarts && remaining > 0 &&
                    bestValue[0] < valueAtLastRestart) {
                restarts++;
                valueAtLastRestart = bestValue[0];
                previousIterations[0] = iterations;
                x = bestPoint[0].clone();
                if (settings.numberIterations > 0) {
                    paramsBFGS.max_iterations = remaining;
                }
                if (settings.printToScreen) {
                    LOGGER.info("maximizeWrtParameter: line search failed (" + status + "); restart " + restarts +
                            " from the best point");
                }
                continue;
            }
            convergedToPrecision = false;
            if (reason == Reason.LINE_SEARCH_FAILURE && restarts > 0 && !(bestValue[0] < valueAtLastRestart) &&
                    settings.delta != null) {
                // a restart from the best point with a cleared memory found no decrease at all: the objective has
                // stopped improving, the past/delta stopping rule's condition, so the best point is a minimum to
                // numerical precision
                reason = Reason.CONVERGED;
                convergedToPrecision = true;
                x = bestPoint[0].clone();
            }
            break;
        }
        valueAtLastRestart = Double.POSITIVE_INFINITY;

        seconds = elapsedSeconds();
        time += seconds;

        minimumPoint = x;
        minimumValue = function.valueAt(minimumPoint);
        minimumGradient = function.gradientAt(minimumPoint).clone();
        heldGradient = cachedUnmaskedGradient != null ? cachedUnmaskedGradient.clone() : minimumGradient;
        likelihoodRecomputes = getLikelihoodCalculationCount() - recomputesAtStart;
        useCache = false;
        cachedPoint = null;
        executed = true;

        setParameter(new WrappedVector.Raw(modelValues(minimumPoint)), parameter);

        final boolean failed = reason == Reason.LINE_SEARCH_FAILURE || reason == Reason.NON_FINITE
                || reason == Reason.ERROR;
        if (settings.printToScreen || (failed && firstRun)) {
            LOGGER.info("maximizeWrtParameter " + parameterName() + " stopped: " + reason + " (" + status +
                    ") after " + iterations + " iterations, " + valueCalls + " evaluations, " +
                    String.format("%.3f", seconds) + "s; -log density = " + minimumValue);
        }
        for (Held held : getHeld()) {
            if (settings.printToScreen || held.gradient > 0) {
                LOGGER.info((held.gradient > 0 ? "WARNING: " : "") + "maximizeWrtParameter: held at the lower bound: "
                        + held);
            }
        }
    }

    private String describeHeld() {
        StringBuilder sb = new StringBuilder();
        for (LowerBoundHold hold : settings.lowerBoundHolds) {
            for (int j = 0; j < hold.parameter.getDimension(); ++j) {
                if (heldCoordinates[hold.offset + j]) {
                    sb.append(sb.length() > 0 ? ", " : "").append(hold.name(j)).append(" = ").append(hold.bound);
                }
            }
        }
        return sb.toString();
    }

    private double[] untransformed(double[] point) {
        return transform != null ? transform.inverse(point, 0, point.length) : point.clone();
    }

    // whether an element not yet held lies below its bound at the transformed point
    private boolean belowBound(double[] point) {
        final double[] values = untransformed(point);
        for (LowerBoundHold hold : settings.lowerBoundHolds) {
            for (int j = 0; j < hold.parameter.getDimension(); ++j) {
                int index = hold.offset + j;
                if (!heldCoordinates[index] && values[index] < hold.bound) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Set the elements below their bounds to the bounds and mark them held.
     *
     * @return the new transformed point, or null when no element is below its bound
     */
    private double[] holdAtBounds(double[] point) {
        final double[] values = untransformed(point);
        final List<Integer> newlyHeld = new ArrayList<Integer>();
        for (LowerBoundHold hold : settings.lowerBoundHolds) {
            for (int j = 0; j < hold.parameter.getDimension(); ++j) {
                int index = hold.offset + j;
                if (!heldCoordinates[index] && values[index] < hold.bound) {
                    values[index] = hold.bound;
                    newlyHeld.add(index);
                }
            }
        }
        if (newlyHeld.isEmpty()) {
            return null;
        }
        final double[] transformed = transform != null ? transform.transform(values, 0, values.length) : values;
        // change only the held coordinates, and check that each maps to its own element alone
        final double[] result = point.clone();
        for (int index : newlyHeld) {
            result[index] = transformed[index];
        }
        final double[] check = untransformed(result);
        for (int i = 0; i < check.length; ++i) {
            if (!(Math.abs(check[i] - values[i]) <= 1E-8 * Math.max(1.0, Math.abs(values[i])))) {
                throw new IllegalStateException("holdAtLowerBound needs a transform that acts element by element " +
                        "on the held parameter (log, logit or none)");
            }
        }
        for (int index : newlyHeld) {
            heldCoordinates[index] = true;
            heldValues[index] = values[index];
        }
        return result;
    }

    private String getBanner() {
        StringBuilder sb = new StringBuilder("Maximizing wrt ").append(parameterName())
                .append(" (dim ").append(function.getDimension())
                .append(", ").append(gradient.getGradientType()).append(" gradient")
                .append(", nIterations = ").append(settings.numberIterations > 0 ? settings.numberIterations : "unlimited");
        if (settings.maxRunTime > 0) sb.append(", maxRunTime = ").append(settings.maxRunTime).append("s");
        if (settings.maxEvaluations > 0) sb.append(", maxEvaluations = ").append(settings.maxEvaluations);
        if (settings.stallWindow > 0) {
            sb.append(", stall: relTol = ").append(settings.relTol).append(" over ").append(settings.stallWindow);
        }
        return sb.append(")").toString();
    }

    private static boolean isBounded(Parameter parameter) {
        try {
            Bounds<Double> bounds = parameter.getBounds();
            for (int i = 0; i < parameter.getDimension(); ++i) {
                if (bounds != null && (bounds.getLowerLimit(i) > Double.NEGATIVE_INFINITY
                        || bounds.getUpperLimit(i) < Double.POSITIVE_INFINITY)) {
                    return true;
                }
            }
        } catch (NullPointerException e) { // Parameter.Default, alone or inside a CompoundParameter, without bounds
            return false;
        }
        return false;
    }

    private String parameterName() {
        return parameter.getId() != null ? parameter.getId() : parameter.getParameterName();
    }

    private double elapsedSeconds() {
        return (System.nanoTime() - startNanos) / 1E9;
    }

    private Reason checkBudgets(double fx, int k) {
        if (settings.maxRunTime > 0 && elapsedSeconds() >= settings.maxRunTime) {
            return Reason.MAX_RUN_TIME;
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
            sb.append("Gradient type: ").append(gradient.getGradientType()).append("\n");
            sb.append("Fx: ").append(minimumValue).append("\n");
            sb.append("Time: ").append(time).append("s\n");
            // calls through evaluateLogLikelihood, over every run
            sb.append("Count: ").append(count).append("\n");
            // the latest run
            sb.append("Status: ").append(status).append("\n");
            sb.append("Reason: ").append(reason);
            if (convergedToPrecision) {
                sb.append(" (to numerical precision: a restart from the best point after a failed line search found " +
                        "no decrease)");
            }
            sb.append("\n");
            sb.append("Iterations: ").append(iterations).append("\n");
            sb.append("Restarts: ").append(restarts).append("\n");
            if (!settings.lowerBoundHolds.isEmpty()) {
                sb.append("Holds: ").append(holds).append("\n");
                for (Held held : getHeld()) {
                    sb.append("HeldAtLowerBound: ").append(held).append("\n");
                }
            }
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
            // the same rows as timeTreeExport's params.tsv
            sb.append("Optimum:\n");
            for (Coordinate coordinate : getOptimum()) {
                sb.append("  ").append(coordinate.name).append(" = ").append(coordinate.value);
                if (!coordinate.getStatus().equals("optimized")) {
                    sb.append(" (").append(coordinate.getStatus()).append(")");
                }
                sb.append("\n");
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
        setParameter(new WrappedVector.Raw(modelValues(point)), parameter);
    }

    // the parameter's values at the transformed point, with held elements exactly at their bounds: the transform's
    // round trip can be one ulp off (e.g. exp(log(0.05)) on some JVMs)
    private double[] modelValues(double[] point) {
        final double[] values = transform != null ? transform.inverse(point, 0, point.length) : point.clone();
        if (heldCoordinates != null) {
            for (int i = 0; i < values.length; ++i) {
                if (heldCoordinates[i]) {
                    values[i] = heldValues[i];
                }
            }
        }
        return values;
    }

    /**
     * The objective (-log density) and its gradient at the transformed point y, computed together and cached. A
     * non-finite log density returns a large penalty at once, without computing its gradient, so the line search
     * backs off.
     */
    private double evaluate(double[] point) {

        if (useCache && cachedPoint != null && Arrays.equals(point, cachedPoint)) {
            if (cachedPointApplied && !Arrays.equals(point, parameterPoint)) {
                setPoint(point);
            }
            return cachedValue;
        }

        ++valueCalls;
        // a step can overflow a transform (e.g. a log-scale root height to +Inf); leave the model where it was
        final double[] proposed = transform != null ? transform.inverse(point, 0, point.length) : point;
        final boolean finitePoint = allFinite(proposed);
        double[] argument = proposed;
        double logLikelihood = Double.NaN;
        if (finitePoint) {
            setPoint(point);
            argument = parameter.getParameterValues();
            logLikelihood = evaluateLogLikelihood();
        }

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

        if (heldCoordinates != null) {
            // held coordinates do not move: L-BFGS's steps are zero along a zero gradient after a cleared memory
            cachedUnmaskedGradient = result.clone();
            for (int i = 0; i < result.length; ++i) {
                if (heldCoordinates[i]) {
                    result[i] = 0.0;
                }
            }
        }

        cachedPoint = point.clone();
        cachedPointApplied = finitePoint;
        cachedValue = value;
        cachedGradient = result;

        if (settings.checkGradientEvery > 0 && Double.isFinite(value)
                && gradientCalls % settings.checkGradientEvery == 0) {
            checkGradient(point.clone(), value, result.clone());
        }

        return value;
    }

    private static boolean allFinite(double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Compare g.d with (f(y + hd) - f(y - hd)) / 2h along a random unit direction d, log the relative error and
     * restore the cache and the parameter at y.
     */
    private void checkGradient(double[] point, double value, double[] grad) {

        final double[] direction = new double[point.length];
        for (int i = 0; i < direction.length; ++i) {
            // held coordinates have a zeroed gradient: check the others
            direction[i] = heldCoordinates != null && heldCoordinates[i] ? 0.0 : MathUtils.nextGaussian();
        }
        final double length = norm(direction);
        if (length == 0.0) {
            return;
        }
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
        if (!Double.isFinite(numeric)) {
            // no finite difference here: skip the check
            setPoint(point);
            cachedPoint = point;
            cachedValue = value;
            cachedGradient = grad;
            return;
        }
        final double error = Math.abs(numeric - directional) / Math.max(1.0, Math.abs(numeric));

        gradientCheckErrors.add(error);
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
        if (transform != null && !allFinite(transform.inverse(point, 0, point.length))) {
            return Double.NaN; // the transform overflows here; leave the model where it is
        }
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
            public GradientType getGradientType() {
                return GradientType.NUMERICAL;
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

}
