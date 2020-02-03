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

    double integrateTrajectory(WrappedVector position) {

        WrappedVector velocity = drawInitialVelocity(null);
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);

        BounceState bounceState = new BounceState(drawTotalTravelTime());

        int count = 0;

        while (bounceState.isTimeRemaining()) {

            final MinimumTravelInformation firstBounce;

            MinimumTravelInformation boundaryBounce = getNextBoundaryBounce(
                    position, velocity);

            MinimumTravelInformation gradientBounce = getNextGradientBounceZigzag(velocity, gradient);

            firstBounce = (boundaryBounce.time < gradientBounce.time) ?
                    new MinimumTravelInformation(boundaryBounce.time, boundaryBounce.index, Type.BOUNDARY) :
                    new MinimumTravelInformation(gradientBounce.time, gradientBounce.index, Type.GRADIENT);

            bounceState = doBounce(bounceState, firstBounce, position, velocity, action, gradient, null);

            ++count;
        }
        return 0.0;
    }

    MinimumTravelInformation getNextGradientBounceZigzag(WrappedVector velocity,
                                                         WrappedVector gradient) {

        return getNextGradientBounceZigzag(0, velocity.getDim(),
                velocity.getBuffer(), gradient.getBuffer());
    }

    private MinimumTravelInformation getNextGradientBounceZigzag(final int begin, final int end,
                                                                 final double[] velocity,
                                                                 final double[] gradient) {

        double minimumRoot = Double.POSITIVE_INFINITY;
        double root;
        int index = -1;
        double expRate;

        for (int i = begin; i < end; ++i) {
            if (gradient[i] == 0) { // for fix dimension
                root = Double.POSITIVE_INFINITY;
            } else {
                expRate = -velocity[i] * gradient[i];
                if (expRate < 0) {
                    root = Double.POSITIVE_INFINITY;
                } else {
                    root = MathUtils.nextExponential(expRate);
                }
            }

            if (root < minimumRoot) {
                minimumRoot = root;
                index = i;
            }
        }
        return new MinimumTravelInformation(minimumRoot, index);
    }

    @SuppressWarnings("Duplicates")
    @Override
    BounceState doBounce(BounceState initialBounceState, MinimumTravelInformation firstBounce,
                         WrappedVector position, WrappedVector velocity,
                         WrappedVector action, WrappedVector gradient, WrappedVector momentum) {

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

            WrappedVector column = getPrecisionColumn(eventIndex);

            updateDynamics(position.getBuffer(), velocity.getBuffer(),
                    action.getBuffer(), gradient.getBuffer(), null,
                    column.getBuffer(), eventTime, eventIndex);

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
                                double[] m,
                                double[] c,
                                double time,
                                int index) {

        final double twoV = 2 * v[index];

        for (int i = 0, len = p.length; i < len; ++i) {
            final double gi = g[i];
            final double ai = a[i];

            p[i] = p[i] + time * v[i];
            g[i] = gi - time * ai;
            a[i] = ai - twoV * c[i];
        }
    }

    @Override
    public String getOperatorName() {
        return "Irreversible zig-zag operator";
    }
}
