/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;

import dr.util.TaskPool;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.function.BinaryOperator;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class ReversibleZigZagOperator extends AbstractZigZagOperator implements Reportable {

    public ReversibleZigZagOperator(GradientWrtParameterProvider gradientProvider,
                                    PrecisionMatrixVectorProductProvider multiplicationProvider,
                                    PrecisionColumnProvider columnProvider,
                                    double weight, Options runtimeOptions, Parameter mask,
                                    int threadCount) {

        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, mask, threadCount);
    }

    @Override
    public String getOperatorName() {
        return "Zig-zag particle operator";
    }

    @Override
    double integrateTrajectory(WrappedVector position) {

        String signString;
        if (DEBUG_SIGN) {
            signString = printSign(position);
            System.err.println(signString);
        }

        if (TIMING) {
            timer.startTimer("warmUp");
        }

        WrappedVector momentum = drawInitialMomentum();
        WrappedVector velocity = drawInitialVelocity(momentum);
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);

        BounceState bounceState = new BounceState(drawTotalTravelTime());

        if (TIMING) {
            timer.stopTimer("warmUp");
        }

        int count = 0;

        double[] p, v, a, g, m;
        if (TEST_NATIVE_OPERATOR) {

            p = position.getBuffer().clone();
            v = velocity.getBuffer().clone();
            a = action.getBuffer().clone();
            g = gradient.getBuffer().clone();
            m = momentum.getBuffer().clone();

            nativeZigZag.operate(columnProvider, p, v, a, g, m,
                    bounceState.remainingTime);
        }

        if (TEST_CRITICAL_REGION) {
            nativeZigZag.enterCriticalRegion(position.getBuffer(), velocity.getBuffer(),
                    action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
        }

        if (TIMING) {
            timer.startTimer("integrateTrajectory");
        }

        while (bounceState.isTimeRemaining()) {

            if (DEBUG) {
                debugBefore(position, count);
            }

            final MinimumTravelInformation firstBounce;

            if (taskPool != null) {

                if (FUSE) {

                    if (TIMING) {
                        timer.startTimer("getNext");
                    }

                    firstBounce = getNextBounceParallel(position,
                            velocity, action, gradient, momentum);

                    if (TIMING) {
                        timer.stopTimer("getNext");
                    }

                    if (TEST_NATIVE_BOUNCE) {
                        testNative(firstBounce, position, velocity, action, gradient, momentum);
                    }

                } else {

                    MinimumTravelInformation boundaryBounce = getNextBoundaryBounce(
                            position, velocity);
                    MinimumTravelInformation gradientBounce = getNextGradientBounceParallel(
                            action, gradient, momentum);

                    firstBounce = (boundaryBounce.time < gradientBounce.time) ?
                             new MinimumTravelInformation(boundaryBounce.time, boundaryBounce.index, Type.BOUNDARY) :
                             new MinimumTravelInformation(gradientBounce.time, gradientBounce.index, Type.GRADIENT);
                }

            } else {

                if (FUSE) {

                    if (TIMING) {
                        timer.startTimer("getNext");
                    }

                    firstBounce = getNextBounce(position,
                            velocity, action, gradient, momentum);

                    if (TIMING) {
                        timer.stopTimer("getNext");
                    }

                    if (TEST_NATIVE_BOUNCE) {
                        testNative(firstBounce, position, velocity, action, gradient, momentum);
                    }

                } else {

                    if (TIMING) {
                        timer.startTimer("getNextBoundary");
                    }

                    MinimumTravelInformation boundaryBounce = getNextBoundaryBounce(
                            position, velocity);

                    if (TIMING) {
                        timer.stopTimer("getNextBoundary");
                        timer.startTimer("getNextGradient");
                    }
                    MinimumTravelInformation gradientBounce = getNextGradientBounce(action, gradient, momentum);

                    if (TIMING) {
                        timer.stopTimer("getNextGradient");
                    }

                    firstBounce = (boundaryBounce.time < gradientBounce.time) ?
                            new MinimumTravelInformation(boundaryBounce.time, boundaryBounce.index, Type.BOUNDARY) :
                            new MinimumTravelInformation(gradientBounce.time, gradientBounce.index, Type.GRADIENT);

                }
            }

            bounceState = doBounce(bounceState, firstBounce, position, velocity, action, gradient, momentum);

            if (DEBUG) {
                debugAfter(bounceState, position);
                String newSignString = printSign(position);
                System.err.println(newSignString);
                if (bounceState.type != Type.BOUNDARY && signString.compareTo(newSignString) != 0) {
                    System.err.println("Sign error");
                }
            }

            ++count;
        }

        if (TIMING) {
            timer.stopTimer("integrateTrajectory");
        }

        if (TEST_CRITICAL_REGION) {
            nativeZigZag.exitCriticalRegion();
        }

        if (TEST_NATIVE_OPERATOR) {

            if (!close(p, position.getBuffer())) {

                System.err.println("c: " + new WrappedVector.Raw(p, 0, 10));
                System.err.println("c: " + new WrappedVector.Raw(position.getBuffer(), 0, 10));
            } else {
                System.err.println("close");
            }
        }

        if (DEBUG_SIGN) {
            printSign(position);
        }

        return 0.0;
    }

    private boolean close(double[] lhs, double[] rhs) {
        for (int i = 0; i < lhs.length; ++i) {
            if (Math.abs((lhs[i] - rhs[i]) / (lhs[i] + rhs[i])) > 0.00001) {
                return false;
            }
        }
        return true;
    }

    private void testNative(MinimumTravelInformation firstBounce,
                            WrappedVector position,
                            WrappedVector velocity,
                            WrappedVector action,
                            WrappedVector gradient,
                            WrappedVector momentum) {

        if (TIMING) {
            timer.startTimer("getNextC++");
        }

        final MinimumTravelInformation mti;
        if (TEST_CRITICAL_REGION) {

            if (!nativeZigZag.inCriticalRegion()) {
                nativeZigZag.enterCriticalRegion(position.getBuffer(), velocity.getBuffer(),
                        action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
            }

            mti = nativeZigZag.getNextEventInCriticalRegion();
        } else {
            mti = nativeZigZag.getNextEvent(position.getBuffer(), velocity.getBuffer(),
                    action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
        }

        if (TIMING) {
            timer.stopTimer("getNextC++");
        }

        if (!firstBounce.equals(mti)) {
            System.err.println(mti + " ?= " + firstBounce + "\n");
            System.exit(-1);
        }

    }

    private String printSign(ReadableVector position) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < position.getDim(); ++i) {
            double p = position.get(i);
            if (p < 0.0) sb.append("- ");
            else if (p > 0.0) sb.append("+ ");
            else sb.append("0 ");
        }
        return sb.toString();
    }

    private void debugAfter(BounceState bounceState, ReadableVector position) {
        System.err.println("post position: " + position);
        System.err.println(bounceState);
        System.err.println();
    }


    private void debugBefore(ReadableVector position, int count) {
        System.err.println("before number: " + count);
        System.err.println("init position: " + position);
    }

    private MinimumTravelInformation getNextBounce(WrappedVector position,
                                                   WrappedVector velocity,
                                                   WrappedVector action,
                                                   WrappedVector gradient,
                                                   WrappedVector momentum) {

        return getNextBounce(0, position.getDim(),
                position.getBuffer(), velocity.getBuffer(),
                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());

    }

    private MinimumTravelInformation getNextBounce(final int begin, final int end,
                                                   final double[] position,
                                                   final double[] velocity,
                                                   final double[] action,
                                                   final double[] gradient,
                                                   final double[] momentum) {

        double minimumTime = Double.POSITIVE_INFINITY;
        int index = -1;
        Type type = Type.NONE;

        for (int i = begin; i < end; ++i) {

            double boundaryTime = findBoundaryTime(i, position[i], velocity[i]);

            if (boundaryTime < minimumTime) {
                minimumTime = boundaryTime;
                index = i;
                type = Type.BOUNDARY;
            }

            double gradientTime = findGradientRoot(action[i], gradient[i], momentum[i]);

            if (gradientTime < minimumTime) {
                minimumTime = gradientTime;
                index = i;
                type = Type.GRADIENT;
            }
        }

        return new MinimumTravelInformation(minimumTime, index, type);
    }

    private MinimumTravelInformation getNextGradientBounce(WrappedVector action,
                                                           WrappedVector gradient,
                                                           WrappedVector momentum) {

        return getNextGradientBounce(0, action.getDim(),
                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
    }

    private MinimumTravelInformation getNextGradientBounce(final int begin, final int end,
                                                           final double[] action,
                                                           final double[] gradient,
                                                           final double[] momentum) {

        double minimumRoot = Double.POSITIVE_INFINITY;
        int index = -1;

        for (int i = begin; i < end; ++i) {

            double root = findGradientRoot(action[i], gradient[i], momentum[i]);

            if (root < minimumRoot) {
                minimumRoot = root;
                index = i;
            }
        }

        return new MinimumTravelInformation(minimumRoot, index);
    }

    private MinimumTravelInformation getNextGradientBounceParallel(WrappedVector inAction,
                                                                   WrappedVector inGradient,
                                                                   WrappedVector inMomentum) {

        final double[] action = inAction.getBuffer();
        final double[] gradient = inGradient.getBuffer();
        final double[] momentum = inMomentum.getBuffer();

        TaskPool.RangeCallable<MinimumTravelInformation> map =
                (start, end, thread) -> getNextGradientBounce(start, end, action, gradient, momentum);

        BinaryOperator<MinimumTravelInformation> reduce =
                (lhs, rhs) -> (lhs.time < rhs.time) ? lhs : rhs;

        return taskPool.mapReduce(map, reduce);
    }

    private MinimumTravelInformation getNextBounceParallel(WrappedVector inPosition,
                                                           WrappedVector inVelocity,
                                                           WrappedVector inAction,
                                                           WrappedVector inGradient,
                                                           WrappedVector inMomentum) {

        final double[] position = inPosition.getBuffer();
        final double[] velocity = inVelocity.getBuffer();
        final double[] action = inAction.getBuffer();
        final double[] gradient = inGradient.getBuffer();
        final double[] momentum = inMomentum.getBuffer();

        TaskPool.RangeCallable<MinimumTravelInformation> map =
                (start, end, thread) -> getNextBounce(start, end,
                        position, velocity, action, gradient, momentum);

        BinaryOperator<MinimumTravelInformation> reduce =
                (lhs, rhs) -> (lhs.time < rhs.time) ? lhs : rhs;

        return taskPool.mapReduce(map, reduce);
    }

    private static double findGradientRoot(double action,
                                           double gradient,
                                           double momentum) {

        if (gradient == 0.0) { //for fixed dimension, gradient = 0
            return Double.POSITIVE_INFINITY;
        } else {
            return minimumPositiveRoot(-0.5 * action, gradient, momentum);
        }
    }

    private double findBoundaryTime(int index, double position,
                                    double velocity) {

        double time = Double.POSITIVE_INFINITY;

        if (headingTowardsBoundary(position, velocity, index)) { // Also ensures x != 0.0
            time = Math.abs(position / velocity);
        }

        return time;
    }

    private MinimumTravelInformation getNextBoundaryBounce(WrappedVector inPosition,
                                                           WrappedVector inVelocity) {

        @SuppressWarnings("duplicate")
        final double[] position = inPosition.getBuffer();
        final double[] velocity = inVelocity.getBuffer();

        double minimumTime = Double.POSITIVE_INFINITY;
        int index = -1;

        for (int i = 0, len = position.length; i < len; ++i) {

            double time = findBoundaryTime(i, position[i], velocity[i]);

            if (time < minimumTime) {
                minimumTime = time;
                index = i;
            }
        }

        return new MinimumTravelInformation(minimumTime, index);
    }

    private static double minimumPositiveRoot(double a,
                                              double b,
                                              double c) {
        double signA = sign(a);
        b = b * signA;
        c = c * signA;
        a = a * signA;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        double sqrtDiscriminant = Math.sqrt(discriminant);
        double root = (-b - sqrtDiscriminant) / (2 * a);
        if (root <= 0.0) {
            root = (-b + sqrtDiscriminant) / (2 * a);
        }
        if (root <= 0.0) {
            root = Double.POSITIVE_INFINITY;
        }

        return root;
    }

//    private static double minimumPositiveRoot(double a,
//                                              double b,
//                                              double c) {
//
//        double discriminant = b * b - 4 * a * c;
//        if (discriminant < 0.0) {
//            return Double.POSITIVE_INFINITY;
//        }
//
//        double sqrtDiscriminant = Math.sqrt(discriminant);
//        double root1 = (-b - sqrtDiscriminant) / (2 * a);
//        double root2 = (-b + sqrtDiscriminant) / (2 * a);
//
//        root1 = (root1 > 0.0) ? root1 : Double.POSITIVE_INFINITY;
//        root2 = (root2 > 0.0) ? root2 : Double.POSITIVE_INFINITY;
//
//        return (root1 < root2) ? root1 : root2;
//    }

    private WrappedVector drawInitialMomentum() {

        ReadableVector mass = preconditioning.mass;
        double[] momentum = new double[mass.getDim()];

        for (int i = 0, len = momentum.length; i < len; i++) {
            int sign = (MathUtils.nextDouble() > 0.5) ? 1 : -1;
            momentum[i] = sign * MathUtils.nextExponential(1) * Math.sqrt(mass.get(i));
        }

        if (mask != null) {
            applyMask(momentum);
        }

        return new WrappedVector.Raw(momentum);
    }

    private static int sign(double x) {
        int sign = 0;
        if (x > 0.0) {
            sign = 1;
        } else if (x < 0.0) {
            sign = -1;
        }
        return sign;
    }

    private WrappedVector drawInitialVelocity(WrappedVector momentum) {

        ReadableVector mass = preconditioning.mass;
        double[] velocity = new double[momentum.getDim()];

        for (int i = 0, len = momentum.getDim(); i < len; ++i) {
            velocity[i] = sign(momentum.get(i)) / Math.sqrt(mass.get(i));
        }

        return new WrappedVector.Raw(velocity);
    }

    private BounceState doBounce(BounceState initialBounceState,
                                 MinimumTravelInformation firstBounce,
                                 WrappedVector position, WrappedVector velocity,
                                 WrappedVector action, WrappedVector gradient, WrappedVector momentum) {

        if (TIMING) {
            timer.startTimer("doBounce");
        }

        double remainingTime = initialBounceState.remainingTime;
        double eventTime = firstBounce.time;

        final BounceState finalBounceState;
        if (remainingTime < eventTime) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            finalBounceState = new BounceState(Type.NONE, -1, 0.0);

        } else {

            if (TIMING) {
                timer.startTimer("notUpdateAction");
            }

            final Type eventType = firstBounce.type;
            final int eventIndex = firstBounce.index;

            if (TEST_FUSED_DYNAMICS) {

                WrappedVector column = getPrecisionColumn(eventIndex);

                if (!TEST_NATIVE_INNER_BOUNCE) {

                    updateDynamics(position.getBuffer(), velocity.getBuffer(),
                            action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
                            column.getBuffer(), eventTime, eventIndex);

                } else {

                    nativeZigZag.updateDynamics(position.getBuffer(), velocity.getBuffer(),
                            action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
                            column.getBuffer(), eventTime, eventIndex, eventType.ordinal());
                }

                if (firstBounce.type == Type.BOUNDARY) { // Reflect against boundary

                    reflectMomentum(momentum, position, eventIndex);

                } else { // Bounce caused by the gradient

                    setZeroMomentum(momentum, eventIndex);

                }

                reflectVelocity(velocity, eventIndex);

            } else {

                if (!TEST_NATIVE_INNER_BOUNCE) {
                    updatePosition(position, velocity, eventTime);
                    updateMomentum(momentum, gradient, action, eventTime);

                    if (firstBounce.type == Type.BOUNDARY) { // Reflect against boundary

                        reflectMomentum(momentum, position, eventIndex);

                    } else { // Bounce caused by the gradient

                        setZeroMomentum(momentum, eventIndex);

                    }

                    reflectVelocity(velocity, eventIndex);
                    updateGradient(gradient, eventTime, action);

                } else {

                    if (TEST_CRITICAL_REGION) {
                        nativeZigZag.innerBounceCriticalRegion(eventIndex, eventIndex, eventType.ordinal());
                    } else {
                        nativeZigZag.innerBounce(position.getBuffer(), velocity.getBuffer(),
                                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
                                eventTime, eventIndex, eventType.ordinal());
                    }
                }

                if (TIMING) {
                    timer.stopTimer("notUpdateAction");
                }

                updateAction(action, velocity, eventIndex);
            }

            finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
        }

        if (TIMING) {
            timer.stopTimer("doBounce");
        }

        return finalBounceState;
    }

    private void updateDynamics(double[] p,
                                double[] v,
                                double[] a,
                                double[] g,
                                double[] m,
                                double[] c,
                                double time,
                                int index) {

        final double halfTimeSquared = time * time / 2;
        final double twoV = 2 * v[index];

        for (int i = 0, len = p.length; i < len; ++i) {
            final double gi = g[i]; final double ai = a[i];

            p[i] = p[i] + time * v[i];
            m[i] = m[i] + time * gi - halfTimeSquared * ai;
            g[i] = gi - time * ai;
            a[i] = ai - twoV * c[i];
        }

//        if (mask != null) { // TODO Appears unnecessary
//            applyMask(m);
//        }
    }

    private static void reflectMomentum(WrappedVector momentum,
                                        WrappedVector position,
                                        int eventIndex) {

        momentum.set(eventIndex, -momentum.get(eventIndex));
        position.set(eventIndex, 0.0); // Exactly on boundary to avoid potential round-off error
    }

    private static void setZeroMomentum(WrappedVector momentum,
                                        int gradientEventIndex) {

        momentum.set(gradientEventIndex, 0.0); // Exactly zero on gradient event to avoid potential round-off error
    }

    private static void reflectVelocity(WrappedVector velocity,
                                        int eventIndex) {

        velocity.set(eventIndex, -velocity.get(eventIndex));
    }

    private void updateMomentum(WrappedVector momentum,
                                WrappedVector gradient,
                                WrappedVector action,
                                double eventTime) {

        final double[] m = momentum.getBuffer();
        final double[] g = gradient.getBuffer();
        final double[] a = action.getBuffer();

        final double halfEventTimeSquared = eventTime * eventTime / 2;

        for (int i = 0, len = m.length; i < len; ++i) {
            m[i] += eventTime * g[i] - halfEventTimeSquared * a[i];
        }

        if (mask != null) {
            applyMask(m);
        }
    }

    private final static boolean DEBUG = false;
    private final static boolean DEBUG_SIGN = false;
    private final static boolean FUSE = true;
}
