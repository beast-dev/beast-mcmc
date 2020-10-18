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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperator extends AbstractParticleOperator implements Loggable { //todo: temporarily loggable

    public BouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                  PrecisionMatrixVectorProductProvider multiplicationProvider,
                                  PrecisionColumnProvider columnProvider,
                                  double weight, Options runtimeOptions, NativeCodeOptions nativeOptions, Parameter mask) {
        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, nativeOptions, mask);
    }

    @Override
    public String getOperatorName() {
        return "Bouncy particle operator";
    }

    @Override
    double integrateTrajectory(WrappedVector position, WrappedVector momentum) {

        WrappedVector velocity = drawInitialVelocity();
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);

        BounceState bounceState = new BounceState(drawTotalTravelTime());

        initializeNumEvent();

        while (bounceState.remainingTime > 0) {

            if (bounceState.type == Type.BOUNDARY) {
                updateAction(action, velocity, bounceState.index);
            } else {
                action = getPrecisionProduct(velocity);
            }

            double v_Phi_x = -innerProduct(velocity, gradient);
            double v_Phi_v = innerProduct(velocity, action);

            double tMin = Math.max(0.0, -v_Phi_x / v_Phi_v);
            double U_min = tMin * tMin / 2 * v_Phi_v + tMin * v_Phi_x;

            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, U_min);
            MinimumTravelInformation travelInfo = getTimeToBoundary(position, velocity);
            double refreshTime = getRefreshTime();

            bounceState = doBounce(
                    bounceState.remainingTime, bounceTime, travelInfo, refreshTime,
                    position, velocity, gradient, action
            );

            recordOneMoreEvent();
        }
        storedVelocity = velocity;
        return 0.0;
    }

    private BounceState doBounce(double remainingTime, double bounceTime,
                                 MinimumTravelInformation boundaryInfo,
                                 double refreshTime,
                                 WrappedVector position, WrappedVector velocity,
                                 WrappedVector gradient, WrappedVector action) {

        double timeToBoundary = boundaryInfo.time;
        int boundaryIndex = boundaryInfo.index;
        final BounceState finalBounceState;
        final Type eventType;
        int eventIndex;

        if (remainingTime < Math.min(timeToBoundary, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            finalBounceState = new BounceState(Type.NONE, -1, 0.0);
        } else {
            if (refreshTime < Math.min(timeToBoundary, bounceTime)) { //refreshment event
                eventType = Type.REFRESHMENT;
                eventIndex = -1;
                updatePosition(position, velocity, refreshTime);
                updateGradient(gradient, refreshTime, action);

                refreshVelocity(velocity);
            } else if (timeToBoundary < bounceTime) { // Reflect against the boundary
                eventType = Type.BOUNDARY;
                eventIndex = boundaryIndex;

                updatePosition(position, velocity, timeToBoundary);
                updateGradient(gradient, timeToBoundary, action);

                position.set(boundaryIndex, 0.0);
                velocity.set(boundaryIndex, -1 * velocity.get(boundaryIndex));

                remainingTime -= timeToBoundary;

            } else { // Bounce caused by the gradient
                eventType = Type.GRADIENT;
                eventIndex = -1;

                updatePosition(position, velocity, bounceTime);
                updateGradient(gradient, bounceTime, action);
                updateVelocity(velocity, gradient, preconditioning.mass);
                remainingTime -= bounceTime;
            }
            finalBounceState = new BounceState(eventType, eventIndex, remainingTime);
        }
        return finalBounceState;
    }

    private WrappedVector drawInitialVelocity() {

        ReadableVector mass = preconditioning.mass;
        double[] velocity = new double[mass.getDim()];

        for (int i = 0, len = velocity.length; i < len; i++) {
            velocity[i] = MathUtils.nextGaussian() / Math.sqrt(mass.get(i));
        }

        if (mask != null) {
            applyMask(velocity);
        }

        return new WrappedVector.Raw(velocity);
    }

    private MinimumTravelInformation getTimeToBoundary(ReadableVector position, ReadableVector velocity) {

        assert (position.getDim() == velocity.getDim());

        int index = -1;
        double minTime = Double.MAX_VALUE;

        for (int i = 0, len = position.getDim(); i < len; ++i) {

            // TODO Here is where we check that x_j > x_i for categorical dimensions

            // TODO I believe we can simply the condition below (for fixed boundaries) with:
            // double travelTime = -position.get(i) / velocity.get(i); // This is only true for boundaries at 0
            // if (travelTime > 0.0 && missingDataMask[positionIndex] == 0.0)

            double travelTime = Math.abs(position.get(i) / velocity.get(i));
            if (travelTime > 0.0 && headingTowardsBoundary(position.get(i), velocity.get(i), i)) {

                if (travelTime < minTime) {
                    index = i;
                    minTime = travelTime;
                }
            }
        }

        return new MinimumTravelInformation(minTime, index);
    }

    private double getRefreshTime() {
        return refreshmentRate > 0 ? MathUtils.nextExponential(1) / refreshmentRate : Double.POSITIVE_INFINITY;
    }

    @SuppressWarnings("all")
    private double getBounceTime(double v_phi_v, double v_phi_x, double u_min) {
        double a = v_phi_v / 2;
        double b = v_phi_x;
        double c = - u_min - MathUtils.nextExponential(1);
        return (-b + Math.sqrt(b * b - 4 * a * c)) / 2 / a;
    }

    private static void updateVelocity(WrappedVector velocity, WrappedVector gradient, ReadableVector mass) {

        ReadableVector gDivM = new ReadableVector.Quotient(gradient, mass);

        double vg = innerProduct(velocity, gradient);
        double ggDivM = innerProduct(gradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i, velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
        }
    }

    private void refreshVelocity(WrappedVector velocity) {
        ReadableVector mass = preconditioning.mass;

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i, MathUtils.nextGaussian() / Math.sqrt(mass.get(i)));
        }

        if (mask != null) {
            applyMask(velocity);
        }
    }

    private WrappedVector storedVelocity;

    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[preconditioning.mass.getDim()];
        for (int i = 0; i < preconditioning.mass.getDim(); ++i) {
            final int index = i;
            columns[i] = new NumberColumn("v" + index) {
                @Override
                public double getDoubleValue() {
                    if (storedVelocity != null) {
                        return storedVelocity.get(index);
                    } else {
                        return 0.0;
                    }
                }
            };
        }

        return columns;
    }

    private final double refreshmentRate = 0.0;
}
