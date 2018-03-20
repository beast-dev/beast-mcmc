/*
 * NewHamiltonianMonteCarloOperator.java
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

package dr.inference.operators.hmc;

import dr.evolution.alignment.PatternList;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class ZigZagOperator extends AbstractParticleOperator {

    public ZigZagOperator(GradientWrtParameterProvider gradientProvider,
                          PrecisionMatrixVectorProductProvider multiplicationProvider,
                          PrecisionColumnProvider columnProvider,
                          double weight, Options runtimeOptions, Parameter mask, PatternList patternList) {

        super(gradientProvider, multiplicationProvider, weight, runtimeOptions, mask, patternList);
        this.columnProvider = columnProvider;
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
                ++count;
            }

            MinimumTravelInformation boundaryBounce = getNextBoundaryBounce(position, velocity);
            MinimumTravelInformation gradientBounce = getNextGradientBounce(action, gradient, momentum);

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

    enum Type {
        NONE,
        BOUNDARY,
        GRADIENT
    }

    class BounceState {
        final Type type;
        final int index;
        final double remainingTime;

        BounceState(Type type, int index, double remainingTime) {
            this.type = type;
            this.index = index;
            this.remainingTime = remainingTime;
        }

        BounceState(double remainingTime) {
            this.type = Type.NONE;
            this.index = -1;
            this.remainingTime = remainingTime;
        }

        boolean isTimeRemaining() {
            return remainingTime > 0.0;
        }

        public String toString() {
            return "remainingTime : " + remainingTime + "\n" +
                    "lastBounceType: " + type + " in dim: " + index;
        }
    }

    private MinimumTravelInformation getNextGradientBounce(ReadableVector action,
                                                           ReadableVector gradient,
                                                           ReadableVector momentum) {

        double minimumRoot = Double.POSITIVE_INFINITY;
        int index = -1;

        for (int i = 0, len = action.getDim(); i < len; ++i) {

            double root = minimumPositiveRoot(action.get(i) / 2, -gradient.get(i), -momentum.get(i));
            if (root < minimumRoot) {
                minimumRoot = root;
                index = i;
            }
        }

        return new MinimumTravelInformation(minimumRoot, index);
    }

    private MinimumTravelInformation getNextBoundaryBounce(ReadableVector position,
                                                           ReadableVector velocity) {


        double minimumTime = Double.POSITIVE_INFINITY;
        int index = -1;

        for (int i = 0, len = position.getDim(); i < len; ++i) {

            double x = position.get(i);
            double v = velocity.get(i);

            if (headingTowardsBoundary(x, v)) { // Also ensures x != 0.0
                double time = Math.abs(x / v);
                if (time < minimumTime) {
                    minimumTime = time;
                    index = i;
                }
            }
        }

        return new MinimumTravelInformation(minimumTime, index);
    }

    private static double minimumPositiveRoot(double a,
                                              double b,
                                              double c) {

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

    private ReadableVector getPrecisionColumn(int index) {

        double[] precisionColumn = columnProvider.getColumn(index);

        if (mask != null) {
            applyMask(precisionColumn);
        }

        return new WrappedVector.Raw(precisionColumn);
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

            }

            reflectVelocity(velocity, eventIndex);

            updateGradient(gradient, eventTime, action);
            updateAction(action, velocity, eventIndex);

            finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
        }

        return finalBounceState;
    }

    private void updateAction(WrappedVector action, ReadableVector velocity, int eventIndex) {

        ReadableVector column = getPrecisionColumn(eventIndex);

        double v = velocity.get(eventIndex);
        for (int i = 0, len = action.getDim(); i < len; ++i) {
            action.set(i,
                    action.get(i) + 2 * v * column.get(i)
            );
        }

        if (mask != null) {
            applyMask(action);
        }
    }

    private static void reflectMomentum(WrappedVector momentum,
                                        WrappedVector position,
                                        int eventIndex) {

        momentum.set(eventIndex, -momentum.get(eventIndex));
        position.set(eventIndex, 0.0); // Exactly on boundary to avoid potential round-off error
    }

    private static void reflectVelocity(WrappedVector velocity,
                                        int eventIndex) {

        velocity.set(eventIndex, -velocity.get(eventIndex));
    }

    private void updateMomentum(WrappedVector momentum,
                                ReadableVector gradient,
                                ReadableVector action,
                                double eventTime) {

        for (int i = 0, len = momentum.getDim(); i < len; ++i) {
            momentum.set(i,
                    momentum.get(i) + eventTime * gradient.get(i) - eventTime * eventTime * action.get(i) / 2
            );
        }

        if (mask != null) {
            applyMask(momentum);
        }
    }

    private final PrecisionColumnProvider columnProvider;

    private final static boolean DEBUG = false;
    private final static boolean DEBUG_SIGN = false;
}
