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
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperator extends AbstractParticleOperator {

    public BouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                  PrecisionMatrixVectorProductProvider multiplicationProvider,
                                  double weight, Options runtimeOptions, Parameter mask) {
        super(gradientProvider, multiplicationProvider, weight, runtimeOptions, mask);
    }
    
    @Override
    public String getOperatorName() {
        return "Bouncy particle operator";
    }

    @Override
    double integrateTrajectory(WrappedVector position) {

        WrappedVector velocity = drawInitialVelocity();
        WrappedVector gradient = getInitialGradient();

        double remainingTime = drawTotalTravelTime();
        while (remainingTime > 0) {

            ReadableVector Phi_v = getPrecisionProduct(velocity);

            double v_Phi_x = - innerProduct(velocity, gradient);
            double v_Phi_v = innerProduct(velocity, Phi_v);

            double tMin = Math.max(0.0, - v_Phi_x / v_Phi_v);
            double U_min = tMin * tMin / 2 * v_Phi_v + tMin * v_Phi_x;

            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, U_min);
            MinimumTravelInformation travelInfo = getTimeToBoundary(position, velocity);

            remainingTime = doBounce(
                    remainingTime, bounceTime, travelInfo,
                    position, velocity, gradient, Phi_v
            );
        }

        return 0.0;
    }

    private double doBounce(double remainingTime, double bounceTime,
                            MinimumTravelInformation travelInfo,
                            WrappedVector position, WrappedVector velocity,
                            WrappedVector gradient, ReadableVector Phi_v) {

        double timeToBoundary = travelInfo.minTime;
        int boundaryIndex = travelInfo.minIndex;

        if (remainingTime < Math.min(timeToBoundary, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            remainingTime = 0.0;

        } else if (timeToBoundary < bounceTime) { // Reflect against the boundary

            updatePosition(position, velocity, timeToBoundary);
            updateGradient(gradient, timeToBoundary, Phi_v);

            position.set(boundaryIndex, 0.0);
            velocity.set(boundaryIndex, -1 * velocity.get(boundaryIndex));

            remainingTime -= timeToBoundary;

        } else { // Bounce caused by the gradient

            updatePosition(position, velocity, bounceTime);
            updateGradient(gradient, bounceTime, Phi_v);
            updateVelocity(velocity, gradient, preconditioning.mass);

            remainingTime -= bounceTime;

        }

        return remainingTime;
    }
}
