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

import dr.evomodel.continuous.hmc.TaxonTaskPool;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.function.BinaryOperator;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class ZigZagOperator extends AbstractParticleOperator {

    public ZigZagOperator(GradientWrtParameterProvider gradientProvider,
                          PrecisionMatrixVectorProductProvider multiplicationProvider,
                          PrecisionColumnProvider columnProvider,
                          double weight, Options runtimeOptions, Parameter mask) {

        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, mask);

        if (PARALLEL) {
            int numberOfThreads = 4;
//            customThreadPool = new ForkJoinPool(numberOfThreads);
            taskPool = new TaxonTaskPool(gradientProvider.getDimension(), numberOfThreads);
        } else {
//            customThreadPool = null;
            taskPool = null;
        }
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

        WrappedVector momentum = drawInitialMomentum();
        WrappedVector velocity = drawInitialVelocity(momentum);
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);

        BounceState bounceState = new BounceState(drawTotalTravelTime());

        int count = 0;

        while (bounceState.isTimeRemaining()) {

            if (DEBUG) {
                debugBefore(position, count);
            }

            MinimumTravelInformation gradientBounce;
            MinimumTravelInformation boundaryBounce;

            if (PARALLEL) {

                boundaryBounce = getNextBoundaryBounce(position, velocity);

//                MinimumTravelInformation test = getNextGradientBounce(action, gradient, momentum, bounceState);
                gradientBounce = getNextGradientBounceParallel(action, gradient, momentum, bounceState);

//                System.err.println(gradientBounce + " -- " + gradientBounce);

            } else {

                boundaryBounce = getNextBoundaryBounce(position, velocity);
                gradientBounce = getNextGradientBounce(action, gradient, momentum, bounceState);
            }

            if (DEBUG) {
                System.err.println("boundary: " + boundaryBounce);
                System.err.println("gradient: " + gradientBounce);
            }

            bounceState = doBounce(bounceState, boundaryBounce, gradientBounce,
                    position, velocity, momentum, gradient, action);

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

        if (DEBUG_SIGN) {
            printSign(position);
        }

        return 0.0;
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

    private MinimumTravelInformation getNextGradientBounce(WrappedVector inAction,
                                                           WrappedVector inGradient,
                                                           WrappedVector inMomentum,
                                                           BounceState bounceState) {

        final double[] action = inAction.getBuffer();
        final double[] gradient = inGradient.getBuffer();
        final double[] momentum = inMomentum.getBuffer();

        return getNextGradientBounce2(0, action.length, action, gradient, momentum, bounceState);
    }

    private MinimumTravelInformation getNextGradientBounce2(final int begin, final int end,
                                                            final double[] action,
                                                            final double[] gradient,
                                                            final double[] momentum,
                                                            final BounceState bounceState) {

        double minimumRoot = Double.POSITIVE_INFINITY;
        int index = -1;

        for (int i = begin; i < end; ++i) {

            double root = findGradientRoot(i, action[i], gradient[i], momentum[i], bounceState);

            if (root < minimumRoot) {
                minimumRoot = root;
                index = i;
            }
        }

        return new MinimumTravelInformation(minimumRoot, index);
    }

    private MinimumTravelInformation getNextGradientBounceParallel(WrappedVector inAction,
                                                                   WrappedVector inGradient,
                                                                   WrappedVector inMomentum,
                                                                   BounceState bounceState) {

        final double[] action = inAction.getBuffer();
        final double[] gradient = inGradient.getBuffer();
        final double[] momentum = inMomentum.getBuffer();

//        MinimumTravelInformation result = null;
//
//        try {
//            result = customThreadPool.submit(
//                    () -> IntStream.range(0, dim).parallel().mapToObj(
//                            (index) -> new MinimumTravelInformation(
//                                    findGradientRoot(index, action[index], gradient[index],
//                                            momentum[index], bounceState
//                                    ), index)
//                    ).reduce(new MinimumTravelInformation(
//                                    Double.POSITIVE_INFINITY, -1),
//                            (lhs, rhs) -> (lhs.time < rhs.time) ? lhs : rhs)
//            ).get();
//
//        } catch (InterruptedException | ExecutionException e) {
//            e.printStackTrace();
//        }

        TaxonTaskPool.RangeCallable<MinimumTravelInformation> map =
                (start, end, thread) -> getNextGradientBounce2(start, end, action, gradient, momentum, bounceState);

        BinaryOperator<MinimumTravelInformation> reduce =
                (lhs, rhs) -> (lhs.time < rhs.time) ? lhs : rhs;

        return taskPool.mapReduce(map, reduce);
    }

    private static double findGradientRoot(int index,
                                           double action,
                                           double gradient,
                                           double momentum,
                                           BounceState bounceState) {

        double root = Double.POSITIVE_INFINITY;

        if (gradient == 0.0) {
            return root;
        }

        if (bounceState.type == Type.GRADIENT && index == bounceState.index) {
            if (gradient * action > 0.0) {
                root = gradient * 2.0 / action;
            }
        } else {
            root = minimumPositiveRoot(-action / 2, gradient, momentum);
        }

        return root;
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
                                 MinimumTravelInformation boundaryBounce,
                                 MinimumTravelInformation gradientBounce,
                                 WrappedVector position, WrappedVector velocity,
                                 WrappedVector momentum,
                                 WrappedVector gradient, WrappedVector action) {

        double remainingTime = initialBounceState.remainingTime;
        double eventTime = Math.min(boundaryBounce.time, gradientBounce.time);

        final BounceState finalBounceState;
        if (remainingTime < eventTime) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            finalBounceState = new BounceState(Type.NONE, -1, 0.0);

        } else {

            updatePosition(position, velocity, eventTime);
            updateMomentum(momentum, gradient, action, eventTime);

            final Type eventType;
            final int eventIndex;
            if (boundaryBounce.time < gradientBounce.time) { // Reflect against the boundary

                eventType = Type.BOUNDARY;
                eventIndex = boundaryBounce.index;

                reflectMomentum(momentum, position, eventIndex);

            } else { // Bounce caused by the gradient
                eventType = Type.GRADIENT;
                eventIndex = gradientBounce.index;
                setZeroMomentum(momentum, eventIndex);

            }

            reflectVelocity(velocity, eventIndex);

            updateGradient(gradient, eventTime, action);
            updateAction(action, velocity, eventIndex);

            finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
        }

        return finalBounceState;
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

//        if (mask == null) {
//            for (int i = 0, len = m.length; i < len; ++i) {
//                m[i] += eventTime * g[i] - halfEventTimeSquared * a[i];
//            }
//        } else {
//            for (int i = 0, len = m.length; i < len; ++i) {
//                m[i] = maskVector[i] * (m[i] + eventTime * g[i] - halfEventTimeSquared * a[i]);
//            }
//        }

        for (int i = 0, len = m.length; i < len; ++i) {
            m[i] += eventTime * g[i] - halfEventTimeSquared * a[i];
        }

        if (mask != null) {
            applyMask(m);
        }
    }

    private final static boolean DEBUG = false;
    private final static boolean DEBUG_SIGN = false;
    private final static boolean PARALLEL = true;

//    private final ForkJoinPool customThreadPool;
    private final TaxonTaskPool taskPool;
//    private final int dim;
}
