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
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.TaskPool;
import dr.xml.Reportable;

import java.util.function.BinaryOperator;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class ReversibleZigZagOperator extends AbstractZigZagOperator implements Reportable, ReversibleHMCProvider {

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

//    double integrateTrajectory(WrappedVector position, WrappedVector momentum) {
//
//        String signString;
//        if (DEBUG_SIGN) {
//            signString = printSign(position);
//            System.err.println(signString);
//        }
//
//        if (TIMING) {
//            timer.startTimer("warmUp");
//        }
//
//        WrappedVector velocity = drawInitialVelocity(momentum);
//        WrappedVector gradient = getInitialGradient();
//        WrappedVector action = getPrecisionProduct(velocity);
//
//        BounceState bounceState = new BounceState(drawTotalTravelTime());
//
//        if (TIMING) {
//            timer.stopTimer("warmUp");
//        }
//
//        int count = 0;
//        initializeNumEvent();
//        double[] p, v, a, g, m;
//        if (TEST_NATIVE_OPERATOR) {
//
//            p = position.getBuffer();
//            v = velocity.getBuffer();
//            a = action.getBuffer();
//            g = gradient.getBuffer();
//            m = momentum.getBuffer();
//
//            nativeZigZag.operate(columnProvider, p, v, a, g, m,
//                    bounceState.remainingTime);
//        } else {
//
//            if (TEST_CRITICAL_REGION) {
//                nativeZigZag.enterCriticalRegion(position.getBuffer(), velocity.getBuffer(),
//                        action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
//            }
//
//            if (TIMING) {
//                timer.startTimer("integrateTrajectory");
//            }
//
//            while (bounceState.isTimeRemaining()) {
//
//                if (DEBUG) {
//                    debugBefore(position, count);
//                }
//
//                final MinimumTravelInformation firstBounce;
//
//                if (taskPool != null) {
//
//                    firstBounce = getNextBounceParallel(position,
//                            velocity, action, gradient, momentum);
//
//                    if (TEST_NATIVE_BOUNCE) {
//                        testNative(firstBounce, position, velocity, action, gradient, momentum);
//                    }
//
//                } else {
//
//                    firstBounce = getNextBounceSerial(position,
//                            velocity, action, gradient, momentum);
//
//                    if (TEST_NATIVE_BOUNCE) {
//                        testNative(firstBounce, position, velocity, action, gradient, momentum);
//                    }
//
//                }
//
//                bounceState = doBounce(bounceState, firstBounce, position, velocity, action, gradient, momentum);
//                recordOneMoreEvent();
//                if (DEBUG) {
//                    debugAfter(bounceState, position);
//                    String newSignString = printSign(position);
//                    System.err.println(newSignString);
//                    if (bounceState.type != Type.BOUNDARY && signString.compareTo(newSignString) != 0) {
//                        System.err.println("Sign error");
//                    }
//                }
//
//                ++count;
//            }
//
//            if (TIMING) {
//                timer.stopTimer("integrateTrajectory");
//            }
//
//            if (TEST_CRITICAL_REGION) {
//                nativeZigZag.exitCriticalRegion();
//            }
//
//            if (TEST_NATIVE_OPERATOR) {
//
//                if (!close(p, position.getBuffer())) {
//
//                    System.err.println("c: " + new WrappedVector.Raw(p, 0, 10));
//                    System.err.println("c: " + new WrappedVector.Raw(position.getBuffer(), 0, 10));
//                } else {
//                    System.err.println("close");
//                }
//            }
//
//            if (DEBUG_SIGN) {
//                printSign(position);
//            }
//        }
//
//        return 0.0;
//    }


    @Override
    MinimumTravelInformation getNextBounce(WrappedVector position,
                                           WrappedVector velocity,
                                           WrappedVector action,
                                           WrappedVector gradient,
                                           WrappedVector momentum) {

        final MinimumTravelInformation firstBounce;

        if (TIMING) {
            timer.startTimer("getNext");
        }

        if (taskPool != null) {

            firstBounce = getNextBounceParallel(position, velocity, action, gradient, momentum);

        } else {

            firstBounce = getNextBounceSerial(position, velocity, action, gradient, momentum);

        }

        if (TIMING) {
            timer.stopTimer("getNext");
        }

        if (TEST_NATIVE_BOUNCE) {
            testNative(firstBounce, position, velocity, action, gradient, momentum);
        }

        return firstBounce;
    }

    private MinimumTravelInformation getNextBounceSerial(WrappedVector position,
                                                         WrappedVector velocity,
                                                         WrappedVector action,
                                                         WrappedVector gradient,
                                                         WrappedVector momentum) {

        return getNextBounceImpl(0, position.getDim(),
                position.getBuffer(), velocity.getBuffer(),
                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());
    }

    @SuppressWarnings("Duplicates")
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
                (start, end, thread) -> getNextBounceImpl(start, end,
                        position, velocity, action, gradient, momentum);

        BinaryOperator<MinimumTravelInformation> reduce =
                (lhs, rhs) -> (lhs.time < rhs.time) ? lhs : rhs;

        return taskPool.mapReduce(map, reduce);
    }


    @SuppressWarnings("Duplicates")
    private MinimumTravelInformation getNextBounceImpl(final int begin, final int end,
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

    private void testNative(MinimumTravelInformation firstBounce,
                            WrappedVector position,
                            WrappedVector velocity,
                            WrappedVector action,
                            WrappedVector gradient,
                            WrappedVector momentum) {

        if (TIMING) {
            timer.startTimer("getNextC++");
        }

        final MinimumTravelInformation mti = nativeZigZag.getNextEvent(position.getBuffer(), velocity.getBuffer(),
                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer());

        if (TIMING) {
            timer.stopTimer("getNextC++");
        }

        if (!firstBounce.equals(mti)) {
            System.err.println(mti + " ?= " + firstBounce + "\n");
            System.exit(-1);
        }
    }

//    @Override
//    final BounceState doBounce(BounceState initialBounceState,
//                               MinimumTravelInformation firstBounce,
//                               WrappedVector position, WrappedVector velocity,
//                               WrappedVector action, WrappedVector gradient, WrappedVector momentum) {
//
//        if (TIMING) {
//            timer.startTimer("doBounce");
//        }
//
//        double remainingTime = initialBounceState.remainingTime;
//        double eventTime = firstBounce.time;
//
//        final Type eventType = firstBounce.type;
//        final int eventIndex = firstBounce.index;
//
//        WrappedVector column = getPrecisionColumn(eventIndex);
//
//        final BounceState finalBounceState;
//        if (remainingTime < eventTime) { // No event during remaining time
//
//            updatePositionMomentum(position.getBuffer(), velocity.getBuffer(), action.getBuffer(),
//                    gradient.getBuffer(), momentum.getBuffer(), remainingTime);
//
//            //todo: for ZZHMC itself (without
//            //todo: NUTS), updating momentum in the end is not necessary.
//
//            finalBounceState = new BounceState(Type.NONE, -1, 0.0);
//
//        } else {
//
//            if (TIMING) {
//                timer.startTimer("notUpdateAction");
//            }
//
//
//            if (TEST_FUSED_DYNAMICS) {
//
//
//                if (!TEST_NATIVE_INNER_BOUNCE) {
//
//                    updateDynamics(position.getBuffer(), velocity.getBuffer(),
//                            action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
//                            column.getBuffer(), eventTime, eventIndex);
//
//                } else {
//
//                    nativeZigZag.updateDynamics(position.getBuffer(), velocity.getBuffer(),
//                            action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
//                            column.getBuffer(), eventTime, eventIndex, eventType.ordinal());
//                }
//
//                if (firstBounce.type == Type.BOUNDARY) { // Reflect against boundary
//
//                    reflectMomentum(momentum, position, eventIndex);
//
//                } else { // Bounce caused by the gradient
//
//                    setZeroMomentum(momentum, eventIndex);
//
//                }
//
//                reflectVelocity(velocity, eventIndex);
//
//            } else {
//
//                if (!TEST_NATIVE_INNER_BOUNCE) {
//                    updatePosition(position, velocity, eventTime);
//                    updateMomentum(momentum, gradient, action, eventTime);
//
//                    if (firstBounce.type == Type.BOUNDARY) { // Reflect against boundary
//
//                        reflectMomentum(momentum, position, eventIndex);
//
//                    } else { // Bounce caused by the gradient
//
//                        setZeroMomentum(momentum, eventIndex);
//                    }
//
//                    reflectVelocity(velocity, eventIndex);
//                    updateGradient(gradient, eventTime, action);
//
//                } else {
//
//                    if (TEST_CRITICAL_REGION) {
//                        nativeZigZag.innerBounceCriticalRegion(eventTime, eventIndex, eventType.ordinal());
//                    } else {
//                        nativeZigZag.innerBounce(position.getBuffer(), velocity.getBuffer(),
//                                action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
//                                eventTime, eventIndex, eventType.ordinal());
//                    }
//                }
//
//                if (TIMING) {
//                    timer.stopTimer("notUpdateAction");
//                }
//
//                updateAction(action, velocity, eventIndex);
//            }
//
//            finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
//        }
//
//        if (TIMING) {
//            timer.stopTimer("doBounce");
//        }
//
//        return finalBounceState;
//    }


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
            final double gi = g[i];
            final double ai = a[i];

            p[i] = p[i] + time * v[i];
            m[i] = m[i] + time * gi - halfTimeSquared * ai;
            g[i] = gi - time * ai;
            a[i] = ai - twoV * c[i];
        }
    }

    @Override
    void updateDynamics(WrappedVector position,
                        WrappedVector velocity,
                        WrappedVector action,
                        WrappedVector gradient,
                        WrappedVector momentum,
                        WrappedVector column,
                        double eventTime,
                        int eventIndex,
                        Type eventType) {

        if (!TEST_NATIVE_INNER_BOUNCE) {

            updateDynamics(position.getBuffer(), velocity.getBuffer(),
                    action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
                    column.getBuffer(), eventTime, eventIndex);

        } else {

            nativeZigZag.updateDynamics(position.getBuffer(), velocity.getBuffer(),
                    action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(),
                    column.getBuffer(), eventTime, eventIndex, eventType.ordinal());
        }

        if (eventType == Type.BOUNDARY) { // Reflect against boundary

            reflectMomentum(momentum, position, eventIndex);

        } else { // Bounce caused by the gradient

            setZeroMomentum(momentum, eventIndex);

        }
    }

    @Override
    void updatePositionAndMomentum(WrappedVector position,
                                   WrappedVector velocity,
                                   WrappedVector action,
                                   WrappedVector gradient,
                                   WrappedVector momentum,
                                   double time) {

        updatePosition(position.getBuffer(), velocity.getBuffer(), time);
        updateMomentum(action.getBuffer(), gradient.getBuffer(), momentum.getBuffer(), time); // TODO This step is only necessary with NUTS
    }

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum,
                                                 int direction, double time) {

        preconditioning.totalTravelTime = time;
        if (direction == -1) {
            // negate momentum
            negateVector(momentum);
        }
        // integrate
        integrateTrajectory(position, momentum);

        if (direction == -1) {
            //negate momentum again
            negateVector(momentum);
        }
        ReadableVector.Utils.setParameter(position, parameter);
    }

    @Override
    public double[] getInitialPosition() {
        return parameter.getParameterValues();
    }

    @Override
    public double getParameterLogJacobian() { // transform is not allowed yet.
        return 0;
    }

    @Override
    public void setParameter(double[] position) {
        ReadableVector.Utils.setParameter(position, parameter);
    }

    @Override
    public WrappedVector drawMomentum() {
        return drawInitialMomentum();
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
        return gradientProvider.getLikelihood().getLogLikelihood() - getKineticEnergy(momentum) - getParameterLogJacobian();
    }

    @Override
    public double getKineticEnergy(ReadableVector momentum) {

        final int dim = momentum.getDim();
        double energy = 0.0;

        for (int i = 0; i < dim; i++) {
            energy += Math.abs(momentum.get(i));
        }
        return energy;
    }

    @Override
    public double getStepSize() {
        return preconditioning.totalTravelTime;
    }

    private void negateVector(WrappedVector vector) {
        for (int i = 0; i < vector.getDim(); i++) {
            vector.set(i, -vector.get(i));
        }
    }
}
