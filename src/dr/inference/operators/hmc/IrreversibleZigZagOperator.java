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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class IrreversibleZigZagOperator extends AbstractZigZagOperator implements Reportable {

    public IrreversibleZigZagOperator(GradientWrtParameterProvider gradientProvider,
                                      PrecisionMatrixVectorProductProvider multiplicationProvider,
                                      PrecisionColumnProvider columnProvider,
                                      double weight, Options runtimeOptions, Parameter mask,
                                      int threadCount) {

        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, mask, threadCount);
    }

    @Override
    WrappedVector drawInitialMomentum() {
        return null;
    }

    @Override
    WrappedVector drawInitialVelocity(WrappedVector momentum) {
        ReadableVector mass = preconditioning.mass;

        double[] velocity = new double[mass.getDim()];

        for (int i = 0, len = mass.getDim(); i < len; ++i) {
            int sign = (MathUtils.nextDouble() > 0.5) ? 1 : -1;

            velocity[i] = sign / Math.sqrt(mass.get(i));
        }
        if (mask != null) {
            applyMask(velocity);
        }
        return new WrappedVector.Raw(velocity);
    }

    // TODO Not needed
    WrappedVector getSignedByVelocity(WrappedVector original, WrappedVector velocity) {
        double[] signed = new double[original.getDim()];
        for (int i = 0; i < original.getDim(); i++) {
            signed[i] = original.get(i) * velocity.get(i);
        }
        return new WrappedVector.Raw(signed);
    }


    double integrateTrajectory(WrappedVector position) {

        WrappedVector velocity = drawInitialVelocity(null);
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);
        WrappedVector signedGradient = getSignedByVelocity(gradient, velocity);  // TODO Not needed
        WrappedVector signedAction = getSignedByVelocity(action, velocity); // TODO If you get rid of this, then single shared function in AbstractZigZag should be fine

        BounceState bounceState = new BounceState(drawTotalTravelTime());

        int count = 0;

        if (TIMING) {
            timer.startTimer("integrateTrajectory");
        }

        while (bounceState.isTimeRemaining()) {

            final MinimumTravelInformation firstBounce;
//            if (FUSE) {

            if (TIMING) {
                timer.startTimer("getNext");
            }

                firstBounce = getNextBounce(position,
                        velocity, signedAction, signedGradient, null); // TODO Note that velocity is already passed, so no need for "signed" work

            if (TIMING) {
                timer.stopTimer("getNext");
            }

//            } else {
//                MinimumTravelInformation boundaryBounce = getNextBoundaryBounce(
//                        position, velocity);
//
//                MinimumTravelInformation gradientBounce = getNextGradientBounceZigzag(signedGradient, signedAction);
//                firstBounce = (boundaryBounce.time < gradientBounce.time) ?
//                        new MinimumTravelInformation(boundaryBounce.time, boundaryBounce.index, Type.BOUNDARY) :
//                        new MinimumTravelInformation(gradientBounce.time, gradientBounce.index, Type.GRADIENT);
//            }

            bounceState = doBounce(bounceState, firstBounce, position, velocity, signedAction, signedGradient, null);

            ++count;
        }

        if (TIMING) {
            timer.stopTimer("integrateTrajectory");
        }

        return 0.0;
    }

    // TODO Same as in super-class?
    private MinimumTravelInformation getNextBounce(WrappedVector position,
                                                   WrappedVector velocity,
                                                   WrappedVector signedAction,
                                                   WrappedVector signedGradient,
                                                   WrappedVector momentum) {

        return getNextBounce(0, position.getDim(),
                position.getBuffer(), velocity.getBuffer(),
                signedAction.getBuffer(), signedGradient.getBuffer(), null);

    }

    private MinimumTravelInformation getNextBounce(final int begin, final int end,
                                                   final double[] position,
                                                   final double[] velocity,
                                                   final double[] signedAction,
                                                   final double[] signedGradient,
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
            double T = MathUtils.nextExponential(1);
            double gradientTime = getSwitchTime(-signedGradient[i], signedAction[i], T);

            if (gradientTime < minimumTime) {
                minimumTime = gradientTime;
                index = i;
                type = Type.GRADIENT;
            }
        }

        return new MinimumTravelInformation(minimumTime, index, type);
    }

//    MinimumTravelInformation getNextGradientBounceZigzag(WrappedVector signedGradient,
//                                                         WrappedVector signedAction) {
//
//        return getNextGradientBounceZigzag(0, signedGradient.getDim(),
//                signedGradient.getBuffer(), signedAction.getBuffer());
//    }

//    private double findC(double a, double b, double T) {
//        if (a > 0 || b < 0) {
//            return T;
//        } else {
//            return b * b / (4 * a) + T;
//        }
//    }

//    private MinimumTravelInformation getNextGradientBounceZigzag(final int begin, final int end,
//                                                                 final double[] signedGradient,
//                                                                 final double[] signedAction) {
//
//        double minimumRoot = Double.POSITIVE_INFINITY;
//        double root;
//        int index = -1;
//        double a;
//        double b;
//        double c;
//
//        for (int i = begin; i < end; ++i) {
//            if (signedGradient[i] == 0) { // for fix dimension
//                root = Double.POSITIVE_INFINITY;
//            } else {
//                double T = MathUtils.nextExponential(1);
//                root = getSwitchTime(-signedGradient[i], signedAction[i], T);
////                a = -0.5 * signedAction[i];
////                b = signedGradient[i];
////                c = findC(a, b, T);
////                root = minimumPositiveRoot(a, b, c); //if c < 0, root (?)= the large positive root
//            }
//            if (root < minimumRoot) {
//                minimumRoot = root;
//                index = i;
//            }
//        }
//        return new MinimumTravelInformation(minimumRoot, index);
//    }

    private double getSwitchTime(double a, double b, double T) {
        // simulate T such that P(T>= t) = exp(-at-bt^2/2), using uniform random input u
        if (b > 0) {
            if (a < 0)
                return -a / b + Math.sqrt(2 * T / b);
            else       // a >= 0
                return -a / b + Math.sqrt(a * a / (b * b) + 2 * T / b);
        } else if (b == 0) {
            if (a > 0)
                return T / a;
            else
                return Double.POSITIVE_INFINITY;
        } else { // b  < 0
            if (a <= 0)
                return Double.POSITIVE_INFINITY;
            else {
                // a > 0
                double t1 = -a / b;
                if (T <= a * t1 + b * t1 * t1 / 2)
                    return -a / b - Math.sqrt(a * a / (b * b) + 2 * T / b);
                else
                    return Double.POSITIVE_INFINITY;
            }
        }
    }

    @SuppressWarnings("Duplicates")
    @Override
    BounceState doBounce(BounceState initialBounceState, MinimumTravelInformation firstBounce,
                         WrappedVector position, WrappedVector velocity,
                         WrappedVector signedAction, WrappedVector signedGradient, WrappedVector momentum) {

        // TODO Probably shares almost all code with doBounce() in ReversibleZigZagOperator, so move shared
        // TODO code into AbstractZigZagOperator


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

            final Type eventType = firstBounce.type;
            final int eventIndex = firstBounce.index;

            WrappedVector signedColumn = getSignedByVelocity(getPrecisionColumn(eventIndex), velocity);

            updateDynamics(position.getBuffer(), velocity.getBuffer(),
                    signedAction.getBuffer(), signedGradient.getBuffer(),
                    signedColumn.getBuffer(), eventTime, eventIndex);

            reflectVelocity(velocity, eventIndex);
            finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
        }

        if (TIMING) {
            timer.stopTimer("doBounce");
        }

        return finalBounceState;
    }

    @SuppressWarnings("Duplicates")
    private void updateDynamics(double[] p,
                                double[] v,
                                double[] a,
                                double[] g,
                                double[] c,
                                double time,
                                int index) {

        final double twoV = 2 * v[index];

        for (int i = 0, len = p.length; i < len; ++i) {
            final double ai = a[i];

            p[i] = p[i] + time * v[i];
            g[i] = g[i] - time * ai;
            a[i] = ai - twoV * c[i];
        }
        g[index] = -g[index]; // TODO Not needed
        a[index] = -a[index]; // TODO Not needed
    }

    @Override
    public String getOperatorName() {
        return "Irreversible zig-zag operator";
    }
//    private final static boolean FUSE = true;
}
