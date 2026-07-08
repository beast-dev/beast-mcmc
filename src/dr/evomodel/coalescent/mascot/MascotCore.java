/*
 * MascotCore.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */

package dr.evomodel.coalescent.mascot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Allocation-friendly reference engine for a MASCOT-style marginal structured
 * coalescent likelihood and the reverse-mode derivative of the same RK4
 * discretization.
 *
 * The flat parameter layout is per epoch:
 *
 * <pre>
 * log m[0,1], log m[0,2], ..., log m[K-1,K-2], log N[0], ..., log N[K-1]
 * </pre>
 *
 * where rows are source states and columns are destination states.
 */
public final class MascotCore {

    private static final double TIME_TOLERANCE = 1.0e-14;
    private static final double EPS = 1.0e-300;

    private final int stateCount;
    private final double[] boundaries;
    private final double maxStep;
    private final int migrationParametersPerEpoch;
    private final int parametersPerEpoch;
    private final int parameterCount;

    public MascotCore(int stateCount, double[] boundaries, double maxStep) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if (boundaries == null || boundaries.length < 2) {
            throw new IllegalArgumentException("boundaries must include at least 0.0 and infinity");
        }
        if (Math.abs(boundaries[0]) > TIME_TOLERANCE) {
            throw new IllegalArgumentException("the first boundary must be 0.0");
        }
        for (int i = 1; i < boundaries.length; i++) {
            if (!(boundaries[i] > boundaries[i - 1])) {
                throw new IllegalArgumentException("boundaries must be strictly increasing");
            }
        }
        if (!Double.isInfinite(boundaries[boundaries.length - 1])) {
            throw new IllegalArgumentException("the final boundary must be positive infinity");
        }
        if (maxStep <= 0.0 || !Double.isFinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }

        this.stateCount = stateCount;
        this.boundaries = boundaries.clone();
        this.maxStep = maxStep;
        this.migrationParametersPerEpoch = stateCount * (stateCount - 1);
        this.parametersPerEpoch = migrationParametersPerEpoch + stateCount;
        this.parameterCount = getEpochCount() * parametersPerEpoch;
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getEpochCount() {
        return boundaries.length - 1;
    }

    public int getParameterCount() {
        return parameterCount;
    }

    public double logLikelihood(Event[] events, double[] theta) {
        return evaluate(events, theta, false, false).logLikelihood;
    }

    public Result likelihoodAndGradient(Event[] events, double[] theta) {
        return evaluate(events, theta, true, false);
    }

    public Result evaluate(Event[] events, double[] theta, boolean computeGradient, boolean checkProbabilities) {
        if (events == null || events.length == 0) {
            throw new IllegalArgumentException("at least one tree event is required");
        }
        if (theta == null || theta.length != parameterCount) {
            throw new IllegalArgumentException("theta dimension " + (theta == null ? -1 : theta.length) +
                    " does not match expected dimension " + parameterCount);
        }

        Event[] sortedEvents = events.clone();
        Arrays.sort(sortedEvents, EVENT_COMPARATOR);

        ActiveState state = new ActiveState(events.length, stateCount);
        double currentTime = 0.0;
        List<OperationTape> operations = computeGradient ? new ArrayList<OperationTape>() : null;

        for (Event event : sortedEvents) {
            if (event.time < currentTime - TIME_TOLERANCE) {
                throw new IllegalArgumentException("events must be sorted by nondecreasing time");
            }

            while (event.time > currentTime + TIME_TOLERANCE) {
                if (state.activeCount == 0) {
                    currentTime = event.time;
                    break;
                }
                double segmentEnd = nextBoundaryAfter(currentTime, event.time);
                int epoch = epochAt(currentTime + TIME_TOLERANCE);
                integrateSegment(state, theta, currentTime, segmentEnd, epoch, operations);
                currentTime = segmentEnd;
                if (checkProbabilities) {
                    checkProbabilities(state);
                }
            }

            currentTime = event.time;
            if (event.type == EventType.SAMPLE) {
                applySampleEvent(state, event, operations);
            } else {
                applyCoalescentEvent(state, event, theta, operations);
            }

            if (checkProbabilities) {
                checkProbabilities(state);
            }
        }

        double[] gradient = null;
        if (computeGradient) {
            gradient = reverse(operations, theta, state.activeCount);
        }

        return new Result(state.logLikelihood, gradient, state.copyProbabilities(), state.copyActiveIds());
    }

    private void integrateSegment(ActiveState state, double[] theta, double start, double end, int epoch,
                                  List<OperationTape> operations) {
        if (!(end > start)) {
            throw new IllegalArgumentException("empty integration segment");
        }

        int steps = Math.max(1, (int) Math.ceil((end - start) / maxStep));
        double h = (end - start) / steps;
        double t = start;
        double[] y = packState(state.probabilities, state.activeCount, state.logLikelihood);
        IntervalTape intervalTape = operations == null ? null : new IntervalTape();

        for (int i = 0; i < steps; i++) {
            if (intervalTape == null) {
                y = rk4Step(y, state.activeCount, theta, epoch, h);
            } else {
                StepTape stepTape = rk4StepWithTape(y, state.activeIds, state.activeCount, theta, epoch, h, t);
                y = stepTape.yOut;
                intervalTape.steps.add(stepTape);
            }
            t += h;
        }

        unpackState(y, state);
        if (intervalTape != null) {
            operations.add(intervalTape);
        }
    }

    private double[] rk4Step(double[] y, int activeCount, double[] theta, int epoch, double h) {
        double[] k1 = rhs(y, activeCount, theta, epoch);
        double[] y2 = addScaled(y, k1, 0.5 * h);
        double[] k2 = rhs(y2, activeCount, theta, epoch);
        double[] y3 = addScaled(y, k2, 0.5 * h);
        double[] k3 = rhs(y3, activeCount, theta, epoch);
        double[] y4 = addScaled(y, k3, h);
        double[] k4 = rhs(y4, activeCount, theta, epoch);

        double[] out = y.clone();
        for (int i = 0; i < out.length; i++) {
            out[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        return out;
    }

    private StepTape rk4StepWithTape(double[] y, int[] activeIds, int activeCount, double[] theta, int epoch,
                                     double h, double t0) {
        double[] y0 = y.clone();
        double[] k1 = rhs(y0, activeCount, theta, epoch);
        double[] y2 = addScaled(y0, k1, 0.5 * h);
        double[] k2 = rhs(y2, activeCount, theta, epoch);
        double[] y3 = addScaled(y0, k2, 0.5 * h);
        double[] k3 = rhs(y3, activeCount, theta, epoch);
        double[] y4 = addScaled(y0, k3, h);
        double[] k4 = rhs(y4, activeCount, theta, epoch);

        double[] yOut = y0.clone();
        for (int i = 0; i < yOut.length; i++) {
            yOut[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        int[] activeIdsCopy = new int[activeCount];
        System.arraycopy(activeIds, 0, activeIdsCopy, 0, activeCount);
        return new StepTape(h, y0, y2, y3, y4, yOut, activeCount, activeIdsCopy, epoch, t0);
    }

    private double[] rhs(double[] y, int activeCount, double[] theta, int epoch) {
        int stateSize = activeCount * stateCount;
        double[] out = new double[stateSize + 1];
        double[] migration = new double[stateCount * stateCount];
        double[] q = new double[stateCount];
        fillRates(theta, epoch, migration, q);

        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int source = 0; source < stateCount; source++) {
                double p = y[offset + source];
                int row = source * stateCount;
                for (int sink = 0; sink < stateCount; sink++) {
                    out[offset + sink] += p * migration[row + sink];
                }
            }
        }

        double[] sums = new double[stateCount];
        double[] sumsSquares = new double[stateCount];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                double p = y[offset + state];
                sums[state] += p;
                sumsSquares[state] += p * p;
            }
        }

        double[] hValues = new double[stateSize];
        double hazard = 0.0;
        for (int state = 0; state < stateCount; state++) {
            hazard += 0.5 * q[state] * (sums[state] * sums[state] - sumsSquares[state]);
        }

        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            double r = 0.0;
            for (int state = 0; state < stateCount; state++) {
                double h = (sums[state] - y[offset + state]) * q[state];
                hValues[offset + state] = h;
                r += y[offset + state] * h;
            }
            for (int state = 0; state < stateCount; state++) {
                out[offset + state] += -y[offset + state] * (hValues[offset + state] - r);
            }
        }

        out[stateSize] = -hazard;
        return out;
    }

    private VjpResult rhsVjp(double[] y, int activeCount, double[] theta, int epoch, double[] adjointRhs) {
        int stateSize = activeCount * stateCount;
        double[] migration = new double[stateCount * stateCount];
        double[] q = new double[stateCount];
        fillRates(theta, epoch, migration, q);

        double[] gradY = new double[stateSize + 1];
        double[] gradTheta = new double[parameterCount];

        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int source = 0; source < stateCount; source++) {
                double v = 0.0;
                int row = source * stateCount;
                for (int sink = 0; sink < stateCount; sink++) {
                    v += adjointRhs[offset + sink] * migration[row + sink];
                }
                gradY[offset + source] += v;
            }
        }

        int thetaOffset = epoch * parametersPerEpoch;
        int rateIndex = 0;
        for (int source = 0; source < stateCount; source++) {
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = Math.exp(theta[thetaOffset + rateIndex]);
                double contribution = 0.0;
                for (int lineage = 0; lineage < activeCount; lineage++) {
                    int offset = lineage * stateCount;
                    contribution += y[offset + source] *
                            (adjointRhs[offset + sink] - adjointRhs[offset + source]);
                }
                gradTheta[thetaOffset + rateIndex] += rate * contribution;
                rateIndex++;
            }
        }

        double[] sums = new double[stateCount];
        double[] sumsSquares = new double[stateCount];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                double p = y[offset + state];
                sums[state] += p;
                sumsSquares[state] += p * p;
            }
        }

        double[] hValues = new double[stateSize];
        double[] rValues = new double[activeCount];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            double r = 0.0;
            for (int state = 0; state < stateCount; state++) {
                double h = (sums[state] - y[offset + state]) * q[state];
                hValues[offset + state] = h;
                r += y[offset + state] * h;
            }
            rValues[lineage] = r;
        }

        double[] b = new double[activeCount];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                b[lineage] += adjointRhs[offset + state] * y[offset + state];
            }
        }

        double[] cSums = new double[stateCount];
        double[] cValues = new double[stateSize];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                double upstream = adjointRhs[offset + state];
                gradY[offset + state] += upstream * (rValues[lineage] - hValues[offset + state]) +
                        b[lineage] * hValues[offset + state];
                double c = y[offset + state] * (b[lineage] - upstream);
                cValues[offset + state] = c;
                cSums[state] += c;
            }
        }

        double[] gradQ = new double[stateCount];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                double c = cValues[offset + state];
                gradY[offset + state] += q[state] * (cSums[state] - c);
                gradQ[state] += c * (sums[state] - y[offset + state]);
            }
        }

        double ellAdjoint = adjointRhs[stateSize];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * stateCount;
            for (int state = 0; state < stateCount; state++) {
                gradY[offset + state] += -ellAdjoint * hValues[offset + state];
            }
        }
        for (int state = 0; state < stateCount; state++) {
            double pairSum = 0.5 * (sums[state] * sums[state] - sumsSquares[state]);
            gradQ[state] += -ellAdjoint * pairSum;
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < stateCount; state++) {
            gradTheta[etaOffset + state] += -q[state] * gradQ[state];
        }

        return new VjpResult(gradY, gradTheta);
    }

    private void applySampleEvent(ActiveState state, Event event, List<OperationTape> operations) {
        if (event.lineage < 0) {
            throw new IllegalArgumentException("sample event has no lineage id");
        }
        if (findLineage(state.activeIds, state.activeCount, event.lineage) >= 0) {
            throw new IllegalArgumentException("lineage " + event.lineage + " is already active");
        }

        double[] probabilities = normalizedSampleProbabilities(event);
        int index = state.activeCount;
        state.ensureCapacity(index + 1);
        state.activeIds[index] = event.lineage;
        System.arraycopy(probabilities, 0, state.probabilities, index * stateCount, stateCount);
        state.activeCount++;

        if (operations != null) {
            operations.add(new SampleTape(index));
        }
    }

    private void applyCoalescentEvent(ActiveState state, Event event, double[] theta, List<OperationTape> operations) {
        int first = findLineage(state.activeIds, state.activeCount, event.child1);
        int second = findLineage(state.activeIds, state.activeCount, event.child2);
        if (first < 0 || second < 0) {
            throw new IllegalArgumentException("coalescent children are not both active: " +
                    event.child1 + ", " + event.child2);
        }
        if (first == second) {
            throw new IllegalArgumentException("coalescent children must be distinct");
        }
        if (findLineage(state.activeIds, state.activeCount, event.parent) >= 0) {
            throw new IllegalArgumentException("coalescent parent is already active: " + event.parent);
        }

        int epoch = epochAt(event.time);
        double[] q = new double[stateCount];
        fillInversePopulation(theta, epoch, q);

        double[] p1 = new double[stateCount];
        double[] p2 = new double[stateCount];
        System.arraycopy(state.probabilities, first * stateCount, p1, 0, stateCount);
        System.arraycopy(state.probabilities, second * stateCount, p2, 0, stateCount);

        double[] parentProbabilities = new double[stateCount];
        double lambda = 0.0;
        for (int s = 0; s < stateCount; s++) {
            parentProbabilities[s] = p1[s] * p2[s] * q[s];
            lambda += parentProbabilities[s];
        }
        if (!(lambda > 0.0) || !Double.isFinite(lambda)) {
            throw new IllegalArgumentException("invalid coalescent rate: " + lambda);
        }
        for (int s = 0; s < stateCount; s++) {
            parentProbabilities[s] /= lambda;
        }

        int beforeCount = state.activeCount;
        int afterCount = beforeCount - 1;
        int[] newActiveIds = new int[Math.max(state.activeIds.length, afterCount)];
        double[] newProbabilities = new double[Math.max(state.probabilities.length, afterCount * stateCount)];
        int[] keepBeforeIndices = new int[afterCount - 1];
        int out = 0;
        for (int i = 0; i < beforeCount; i++) {
            if (i == first || i == second) {
                continue;
            }
            keepBeforeIndices[out] = i;
            newActiveIds[out] = state.activeIds[i];
            System.arraycopy(state.probabilities, i * stateCount, newProbabilities, out * stateCount, stateCount);
            out++;
        }
        int parentIndexAfter = out;
        newActiveIds[parentIndexAfter] = event.parent;
        System.arraycopy(parentProbabilities, 0, newProbabilities, parentIndexAfter * stateCount, stateCount);

        state.activeIds = newActiveIds;
        state.probabilities = newProbabilities;
        state.activeCount = afterCount;
        state.logLikelihood += Math.log(lambda);

        if (operations != null) {
            operations.add(new CoalescentTape(epoch, first, second, parentIndexAfter, keepBeforeIndices,
                    p1, p2, parentProbabilities, lambda));
        }
    }

    private double[] normalizedSampleProbabilities(Event event) {
        double[] probabilities = new double[stateCount];
        if (event.stateProbabilities == null) {
            if (event.state < 0 || event.state >= stateCount) {
                throw new IllegalArgumentException("sample state out of range: " + event.state);
            }
            probabilities[event.state] = 1.0;
            return probabilities;
        }
        if (event.stateProbabilities.length != stateCount) {
            throw new IllegalArgumentException("sample probability dimension mismatch");
        }
        double sum = 0.0;
        for (int s = 0; s < stateCount; s++) {
            if (event.stateProbabilities[s] < 0.0) {
                throw new IllegalArgumentException("sample probabilities must be nonnegative");
            }
            sum += event.stateProbabilities[s];
        }
        if (!(sum > 0.0)) {
            throw new IllegalArgumentException("sample probabilities must have positive sum");
        }
        for (int s = 0; s < stateCount; s++) {
            probabilities[s] = event.stateProbabilities[s] / sum;
        }
        return probabilities;
    }

    private double[] reverse(List<OperationTape> operations, double[] theta, int finalActiveCount) {
        double[] gradTheta = new double[parameterCount];
        double[] adjointY = new double[finalActiveCount * stateCount + 1];
        adjointY[adjointY.length - 1] = 1.0;

        for (int opIndex = operations.size() - 1; opIndex >= 0; opIndex--) {
            OperationTape operation = operations.get(opIndex);
            if (operation instanceof IntervalTape) {
                IntervalTape interval = (IntervalTape) operation;
                for (int i = interval.steps.size() - 1; i >= 0; i--) {
                    ReverseResult result = reverseStep(interval.steps.get(i), theta, adjointY);
                    adjointY = result.adjointY;
                    addInPlace(gradTheta, result.gradient);
                }
            } else if (operation instanceof CoalescentTape) {
                ReverseResult result = reverseCoalescent((CoalescentTape) operation, adjointY);
                adjointY = result.adjointY;
                addInPlace(gradTheta, result.gradient);
            } else if (operation instanceof SampleTape) {
                adjointY = reverseSample((SampleTape) operation, adjointY);
            } else {
                throw new IllegalArgumentException("unknown tape operation: " + operation.getClass());
            }
        }

        return gradTheta;
    }

    private ReverseResult reverseStep(StepTape step, double[] theta, double[] adjointYOut) {
        double h = step.h;
        double[] gradient = new double[parameterCount];
        double[] adjointY0 = adjointYOut.clone();
        double[] adjointK1 = scale(adjointYOut, h / 6.0);
        double[] adjointK2 = scale(adjointYOut, h / 3.0);
        double[] adjointK3 = scale(adjointYOut, h / 3.0);
        double[] adjointK4 = scale(adjointYOut, h / 6.0);

        VjpResult vjp = rhsVjp(step.y4, step.activeCount, theta, step.epoch, adjointK4);
        addInPlace(gradient, vjp.gradient);
        addInPlace(adjointY0, vjp.adjointY);
        addScaledInPlace(adjointK3, vjp.adjointY, h);

        vjp = rhsVjp(step.y3, step.activeCount, theta, step.epoch, adjointK3);
        addInPlace(gradient, vjp.gradient);
        addInPlace(adjointY0, vjp.adjointY);
        addScaledInPlace(adjointK2, vjp.adjointY, 0.5 * h);

        vjp = rhsVjp(step.y2, step.activeCount, theta, step.epoch, adjointK2);
        addInPlace(gradient, vjp.gradient);
        addInPlace(adjointY0, vjp.adjointY);
        addScaledInPlace(adjointK1, vjp.adjointY, 0.5 * h);

        vjp = rhsVjp(step.y0, step.activeCount, theta, step.epoch, adjointK1);
        addInPlace(gradient, vjp.gradient);
        addInPlace(adjointY0, vjp.adjointY);

        return new ReverseResult(adjointY0, gradient);
    }

    private double[] reverseSample(SampleTape tape, double[] adjointAfter) {
        int afterCount = (adjointAfter.length - 1) / stateCount;
        int beforeCount = afterCount - 1;
        double[] adjointBefore = new double[beforeCount * stateCount + 1];

        int out = 0;
        for (int i = 0; i < afterCount; i++) {
            if (i == tape.sampleIndexAfter) {
                continue;
            }
            System.arraycopy(adjointAfter, i * stateCount, adjointBefore, out * stateCount, stateCount);
            out++;
        }
        adjointBefore[adjointBefore.length - 1] = adjointAfter[adjointAfter.length - 1];
        return adjointBefore;
    }

    private ReverseResult reverseCoalescent(CoalescentTape tape, double[] adjointAfter) {
        int afterCount = (adjointAfter.length - 1) / stateCount;
        int beforeCount = afterCount + 1;
        double[] adjointBefore = new double[beforeCount * stateCount + 1];
        double[] gradient = new double[parameterCount];

        for (int afterIndex = 0; afterIndex < tape.keepBeforeIndices.length; afterIndex++) {
            int beforeIndex = tape.keepBeforeIndices[afterIndex];
            System.arraycopy(adjointAfter, afterIndex * stateCount,
                    adjointBefore, beforeIndex * stateCount, stateCount);
        }

        int parentOffset = tape.parentIndexAfter * stateCount;
        double dot = 0.0;
        for (int s = 0; s < stateCount; s++) {
            dot += adjointAfter[parentOffset + s] * tape.parentProbabilities[s];
        }

        int thetaOffset = tape.epoch * parametersPerEpoch + migrationParametersPerEpoch;
        double ellAdjoint = adjointAfter[adjointAfter.length - 1];
        for (int s = 0; s < stateCount; s++) {
            double p1 = Math.max(tape.p1[s], EPS);
            double p2 = Math.max(tape.p2[s], EPS);
            double centered = adjointAfter[parentOffset + s] - dot;
            double bar = tape.parentProbabilities[s] * centered;

            adjointBefore[tape.child1Index * stateCount + s] += bar / p1;
            adjointBefore[tape.child2Index * stateCount + s] += bar / p2;
            gradient[thetaOffset + s] += -bar;

            adjointBefore[tape.child1Index * stateCount + s] += ellAdjoint * tape.parentProbabilities[s] / p1;
            adjointBefore[tape.child2Index * stateCount + s] += ellAdjoint * tape.parentProbabilities[s] / p2;
            gradient[thetaOffset + s] += -ellAdjoint * tape.parentProbabilities[s];
        }

        adjointBefore[adjointBefore.length - 1] = ellAdjoint;
        return new ReverseResult(adjointBefore, gradient);
    }

    private void fillRates(double[] theta, int epoch, double[] migration, double[] q) {
        int thetaOffset = epoch * parametersPerEpoch;
        int index = 0;
        for (int source = 0; source < stateCount; source++) {
            double rowSum = 0.0;
            int row = source * stateCount;
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = Math.exp(theta[thetaOffset + index]);
                migration[row + sink] = rate;
                rowSum += rate;
                index++;
            }
            migration[row + source] = -rowSum;
        }
        for (int state = 0; state < stateCount; state++) {
            q[state] = Math.exp(-theta[thetaOffset + migrationParametersPerEpoch + state]);
        }
    }

    private void fillInversePopulation(double[] theta, int epoch, double[] q) {
        int offset = epoch * parametersPerEpoch + migrationParametersPerEpoch;
        for (int state = 0; state < stateCount; state++) {
            q[state] = Math.exp(-theta[offset + state]);
        }
    }

    private int epochAt(double t) {
        if (t < -TIME_TOLERANCE) {
            throw new IllegalArgumentException("time is before zero: " + t);
        }
        int epoch = 0;
        while (epoch + 1 < boundaries.length && t >= boundaries[epoch + 1] - TIME_TOLERANCE) {
            epoch++;
        }
        if (epoch >= getEpochCount()) {
            return getEpochCount() - 1;
        }
        return epoch;
    }

    private double nextBoundaryAfter(double start, double stop) {
        for (int i = 1; i < boundaries.length; i++) {
            double boundary = boundaries[i];
            if (boundary > start + TIME_TOLERANCE) {
                if (boundary <= stop + TIME_TOLERANCE) {
                    return Math.min(boundary, stop);
                }
                return stop;
            }
        }
        return stop;
    }

    private double[] packState(double[] probabilities, int activeCount, double logLikelihood) {
        double[] y = new double[activeCount * stateCount + 1];
        System.arraycopy(probabilities, 0, y, 0, activeCount * stateCount);
        y[y.length - 1] = logLikelihood;
        return y;
    }

    private void unpackState(double[] y, ActiveState state) {
        int stateSize = state.activeCount * stateCount;
        state.ensureCapacity(state.activeCount);
        System.arraycopy(y, 0, state.probabilities, 0, stateSize);
        state.logLikelihood = y[stateSize];
    }

    private double[] addScaled(double[] x, double[] dx, double scale) {
        double[] out = x.clone();
        addScaledInPlace(out, dx, scale);
        return out;
    }

    private static double[] scale(double[] x, double scale) {
        double[] out = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            out[i] = scale * x[i];
        }
        return out;
    }

    private static void addInPlace(double[] destination, double[] source) {
        for (int i = 0; i < destination.length; i++) {
            destination[i] += source[i];
        }
    }

    private static void addScaledInPlace(double[] destination, double[] source, double scale) {
        for (int i = 0; i < destination.length; i++) {
            destination[i] += scale * source[i];
        }
    }

    private static int findLineage(int[] activeIds, int activeCount, int lineage) {
        for (int i = 0; i < activeCount; i++) {
            if (activeIds[i] == lineage) {
                return i;
            }
        }
        return -1;
    }

    private void checkProbabilities(ActiveState state) {
        for (int lineage = 0; lineage < state.activeCount; lineage++) {
            double sum = 0.0;
            int offset = lineage * stateCount;
            for (int s = 0; s < stateCount; s++) {
                double p = state.probabilities[offset + s];
                if (p < -1.0e-8 || !Double.isFinite(p)) {
                    throw new IllegalStateException("invalid probability for lineage " +
                            state.activeIds[lineage] + ": " + p);
                }
                sum += p;
            }
            if (Math.abs(sum - 1.0) > 1.0e-6) {
                throw new IllegalStateException("lineage probabilities do not sum to one for lineage " +
                        state.activeIds[lineage] + ": " + sum);
            }
        }
    }

    private static final Comparator<Event> EVENT_COMPARATOR = new Comparator<Event>() {
        @Override
        public int compare(Event first, Event second) {
            int time = Double.compare(first.time, second.time);
            if (time != 0) {
                return time;
            }
            int type = first.type.compareRank - second.type.compareRank;
            if (type != 0) {
                return type;
            }
            return first.sortId() - second.sortId();
        }
    };

    public enum EventType {
        SAMPLE(0),
        COALESCENCE(1);

        private final int compareRank;

        EventType(int compareRank) {
            this.compareRank = compareRank;
        }
    }

    public static final class Event {
        private final double time;
        private final EventType type;
        private final int lineage;
        private final int state;
        private final double[] stateProbabilities;
        private final int child1;
        private final int child2;
        private final int parent;

        private Event(double time, EventType type, int lineage, int state, double[] stateProbabilities,
                      int child1, int child2, int parent) {
            if (time < -TIME_TOLERANCE || !Double.isFinite(time)) {
                throw new IllegalArgumentException("invalid event time: " + time);
            }
            this.time = time;
            this.type = type;
            this.lineage = lineage;
            this.state = state;
            this.stateProbabilities = stateProbabilities == null ? null : stateProbabilities.clone();
            this.child1 = child1;
            this.child2 = child2;
            this.parent = parent;
        }

        public static Event sample(double time, int lineage, int state) {
            return new Event(time, EventType.SAMPLE, lineage, state, null, -1, -1, -1);
        }

        public static Event sample(double time, int lineage, double[] stateProbabilities) {
            return new Event(time, EventType.SAMPLE, lineage, -1, stateProbabilities, -1, -1, -1);
        }

        public static Event coalescence(double time, int child1, int child2, int parent) {
            return new Event(time, EventType.COALESCENCE, -1, -1, null, child1, child2, parent);
        }

        public double getTime() {
            return time;
        }

        public EventType getType() {
            return type;
        }

        public int getLineage() {
            return lineage;
        }

        public int getState() {
            return state;
        }

        public int getChild1() {
            return child1;
        }

        public int getChild2() {
            return child2;
        }

        public int getParent() {
            return parent;
        }

        private int sortId() {
            if (type == EventType.SAMPLE) {
                return lineage;
            }
            return parent;
        }
    }

    public static final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;

        private Result(double logLikelihood, double[] gradient, double[] rootProbabilities, int[] activeLineages) {
            this.logLikelihood = logLikelihood;
            this.gradient = gradient;
            this.rootProbabilities = rootProbabilities;
            this.activeLineages = activeLineages;
        }
    }

    private static final class ActiveState {
        private int[] activeIds;
        private double[] probabilities;
        private int activeCount;
        private double logLikelihood;
        private final int stateCount;

        private ActiveState(int capacity, int stateCount) {
            this.stateCount = stateCount;
            this.activeIds = new int[Math.max(1, capacity)];
            this.probabilities = new double[Math.max(1, capacity * stateCount)];
            this.activeCount = 0;
            this.logLikelihood = 0.0;
        }

        private void ensureCapacity(int capacity) {
            if (activeIds.length < capacity) {
                activeIds = Arrays.copyOf(activeIds, capacity);
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = Arrays.copyOf(probabilities, probabilityCapacity);
            }
        }

        private double[] copyProbabilities() {
            return Arrays.copyOf(probabilities, activeCount * stateCount);
        }

        private int[] copyActiveIds() {
            return Arrays.copyOf(activeIds, activeCount);
        }
    }

    private interface OperationTape {
    }

    private static final class IntervalTape implements OperationTape {
        private final List<StepTape> steps = new ArrayList<StepTape>();
    }

    private static final class StepTape {
        private final double h;
        private final double[] y0;
        private final double[] y2;
        private final double[] y3;
        private final double[] y4;
        private final double[] yOut;
        private final int activeCount;
        @SuppressWarnings("unused")
        private final int[] activeIds;
        private final int epoch;
        @SuppressWarnings("unused")
        private final double t0;

        private StepTape(double h, double[] y0, double[] y2, double[] y3, double[] y4, double[] yOut,
                         int activeCount, int[] activeIds, int epoch, double t0) {
            this.h = h;
            this.y0 = y0;
            this.y2 = y2;
            this.y3 = y3;
            this.y4 = y4;
            this.yOut = yOut;
            this.activeCount = activeCount;
            this.activeIds = activeIds;
            this.epoch = epoch;
            this.t0 = t0;
        }
    }

    private static final class SampleTape implements OperationTape {
        private final int sampleIndexAfter;

        private SampleTape(int sampleIndexAfter) {
            this.sampleIndexAfter = sampleIndexAfter;
        }
    }

    private static final class CoalescentTape implements OperationTape {
        private final int epoch;
        private final int child1Index;
        private final int child2Index;
        private final int parentIndexAfter;
        private final int[] keepBeforeIndices;
        private final double[] p1;
        private final double[] p2;
        private final double[] parentProbabilities;
        @SuppressWarnings("unused")
        private final double lambda;

        private CoalescentTape(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                               int[] keepBeforeIndices, double[] p1, double[] p2,
                               double[] parentProbabilities, double lambda) {
            this.epoch = epoch;
            this.child1Index = child1Index;
            this.child2Index = child2Index;
            this.parentIndexAfter = parentIndexAfter;
            this.keepBeforeIndices = keepBeforeIndices;
            this.p1 = p1;
            this.p2 = p2;
            this.parentProbabilities = parentProbabilities;
            this.lambda = lambda;
        }
    }

    private static final class VjpResult {
        private final double[] adjointY;
        private final double[] gradient;

        private VjpResult(double[] adjointY, double[] gradient) {
            this.adjointY = adjointY;
            this.gradient = gradient;
        }
    }

    private static final class ReverseResult {
        private final double[] adjointY;
        private final double[] gradient;

        private ReverseResult(double[] adjointY, double[] gradient) {
            this.adjointY = adjointY;
            this.gradient = gradient;
        }
    }
}
