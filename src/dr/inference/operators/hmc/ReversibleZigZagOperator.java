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
    final WrappedVector drawInitialMomentum() {

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

    @Override
    final WrappedVector drawInitialVelocity(WrappedVector momentum) {

        ReadableVector mass = preconditioning.mass;
        double[] velocity = new double[momentum.getDim()];

        for (int i = 0, len = momentum.getDim(); i < len; ++i) {
            velocity[i] = sign(momentum.get(i)) / Math.sqrt(mass.get(i));
        }

        return new WrappedVector.Raw(velocity);
    }

    @Override
    final BounceState doBounce(BounceState initialBounceState,
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

//    private static void reflectMomentum(WrappedVector momentum,
//                                        WrappedVector position,
//                                        int eventIndex) {
//
//        momentum.set(eventIndex, -momentum.get(eventIndex));
//        position.set(eventIndex, 0.0); // Exactly on boundary to avoid potential round-off error
//    }
//
//    private static void setZeroMomentum(WrappedVector momentum,
//                                        int gradientEventIndex) {
//
//        momentum.set(gradientEventIndex, 0.0); // Exactly zero on gradient event to avoid potential round-off error
//    }
//
//    private static void reflectVelocity(WrappedVector velocity,
//                                        int eventIndex) {
//
//        velocity.set(eventIndex, -velocity.get(eventIndex));
//    }
//
//    private void updateMomentum(WrappedVector momentum,
//                                WrappedVector gradient,
//                                WrappedVector action,
//                                double eventTime) {
//
//        final double[] m = momentum.getBuffer();
//        final double[] g = gradient.getBuffer();
//        final double[] a = action.getBuffer();
//
//        final double halfEventTimeSquared = eventTime * eventTime / 2;
//
//        for (int i = 0, len = m.length; i < len; ++i) {
//            m[i] += eventTime * g[i] - halfEventTimeSquared * a[i];
//        }
//
//        if (mask != null) {
//            applyMask(m);
//        }
//    }
}
