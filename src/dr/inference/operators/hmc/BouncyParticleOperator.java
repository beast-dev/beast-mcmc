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
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperator extends SimpleMCMCOperator implements GibbsOperator {

    public BouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                  PrecisionMatrixVectorProductProvider multiplicationProvider,
                                  double weight,
                                  double randomTimeWidth) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.parameter = gradientProvider.getParameter();
        this.drawDistribution = new NormalDistribution(0, 1);

        setWeight(weight);
        checkParameterBounds(parameter);

        masses = setupPreconditionedMatrix();

        this.totalTravelTime = 0.05; // TODO Determine totalTravelTime
        this.randomTimeWidth = randomTimeWidth;
    }

    @Override
    public double doOperation() {

        WrappedVector position = getInitialPosition();
        WrappedVector velocity = drawInitialVelocity();
        WrappedVector negativeGradient = getInitialNegativeGradient();

        double remainingTime = drawTotalTravelTime(); //totalTravelTime;
        while (remainingTime > 0) {

            ReadableVector Phi_v = getPrecisionProduct(velocity);

            double v_Phi_x = innerProduct(velocity, negativeGradient);
            double v_Phi_v = innerProduct(velocity, Phi_v);

            double tMin = Math.max(0.0, - v_Phi_x / v_Phi_v);
            double U_min = tMin * tMin / 2 * v_Phi_v + tMin * v_Phi_x;

            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, U_min);
            MinimumTravelInformation travelInfo = getTimeToBoundary(position, velocity);

            remainingTime = doBounce(
                    remainingTime, bounceTime, travelInfo,
                    position, velocity, negativeGradient, Phi_v
            );
        }

        setParameter(position);

        return 0.0;
    }

    private double drawTotalTravelTime() {
        double randomFraction = 1.0 + randomTimeWidth * (MathUtils.nextDouble() - 0.5);
        return totalTravelTime * randomFraction;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy Particle operator";
    }

    private double doBounce(double remainingTime, double bounceTime,
                            MinimumTravelInformation travelInfo,
                            WrappedVector position, WrappedVector velocity,
                            WrappedVector negativeGradient, ReadableVector Phi_v) {

        double timeToBoundary = travelInfo.minTime;
        int boundaryIndex = travelInfo.minIndex;

        if (remainingTime < Math.min(timeToBoundary, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            remainingTime = 0.0;

        } else if (timeToBoundary < bounceTime) { // Bounce against the boundary

            updatePosition(position, velocity, timeToBoundary);
            updateNegativeGradient(negativeGradient, timeToBoundary, Phi_v);

            position.set(boundaryIndex, 0.0);
            velocity.set(boundaryIndex, -1 * velocity.get(boundaryIndex));

            remainingTime -= timeToBoundary;

        } else { // Bounce caused by the gradient

            updatePosition(position, velocity, bounceTime);
            updateNegativeGradient(negativeGradient, bounceTime, Phi_v);
            updateVelocity(velocity, negativeGradient, masses);

            remainingTime -= bounceTime;

        }

        return remainingTime;
    }

    private void setParameter(ReadableVector position) {
        for (int j = 0, dim = position.getDim(); j < dim; ++j) {
            parameter.setParameterValueQuietly(j, position.get(j));
        }
        parameter.fireParameterChangedEvent();
    }

    private void updateVelocity(WrappedVector velocity, WrappedVector negativeGradient, ReadableVector masses) {

        ReadableVector gDivM = new ReadableVector.Quotient(negativeGradient, masses);

        double vg = innerProduct(velocity, negativeGradient); // TODO Isn't this already computed
        double ggDivM = innerProduct(negativeGradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i,
                    velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
        }
    }

    private void updateNegativeGradient(WrappedVector negativeGradient, double time, ReadableVector Phi_v) {
        for (int i = 0, len = negativeGradient.getDim(); i < len; ++i) {
            negativeGradient.set(i, negativeGradient.get(i) + time * Phi_v.get(i));
        }
    }

    private void updatePosition(WrappedVector position, WrappedVector velocity, double time) {
        for (int i = 0, len = position.getDim(); i < len; ++i) {
            position.set(i, position.get(i) + time * velocity.get(i));
        }
    }

    @SuppressWarnings("all")
    private double getBounceTime(double v_phi_v, double v_phi_x, double u_min) {
        double a = v_phi_v / 2;
        double b = v_phi_x;
        double c = u_min - MathUtils.nextExponential(1);
        return (-b + Math.sqrt(b * b - 4 * a * c)) / 2 / a;
    }

    private WrappedVector getInitialNegativeGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();
        for (int i = 0, len = gradient.length; i < len; ++i) {
            gradient[i] = -1 * gradient[i];
        }

        return new WrappedVector.Raw(gradient);
    }

    private double innerProduct(ReadableVector x, ReadableVector y) {

        assert (x.getDim() == y.getDim());

        double sum = 0;
        for (int i = 0, dim = x.getDim(); i < dim; ++i) {
            sum += x.get(i) * y.get(i);
        }

        return sum;
    }

    private ReadableVector getPrecisionProduct(ReadableVector velocity) {

        setParameter(velocity);

        double[] product = productProvider.getProduct(parameter);

        return new WrappedVector.Raw(product);
    }

    private WrappedVector drawInitialVelocity() {

        double[] velocity = new double[masses.getDim()];

        for (int i = 0, len = velocity.length; i < len; i++) {
            velocity[i] = (Double) drawDistribution.nextRandom() / Math.sqrt(masses.get(i));
        }
        return new WrappedVector.Raw(velocity);
    }

    private MinimumTravelInformation getTimeToBoundary(ReadableVector position, ReadableVector velocity) {

        assert (position.getDim() == velocity.getDim());

        int index = -1;
        double minTime = Double.MAX_VALUE;

        for (int i = 0, len = position.getDim(); i < len; ++i) {

            double travelTime = Math.abs(position.get(i) / velocity.get(i));

            if (travelTime > 0.0 && headingAwayFromBoundary(position.get(i), velocity.get(i))) {

                if (travelTime < minTime) {
                    index = i;
                    minTime = travelTime;
                }
            }
        }

        return new MinimumTravelInformation(minTime, index);
    }

    private boolean headingAwayFromBoundary(double position, double velocity) {
        return position * velocity < 0.0;
    }

    private class MinimumTravelInformation {

        double minTime;
        int minIndex;

        private MinimumTravelInformation(double minTime, int minIndex) {
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

    private WrappedVector getInitialPosition() {
        return new WrappedVector.Raw(parameter.getParameterValues());
    }

    private void checkParameterBounds(Parameter parameter) {

        for (int i = 0, len = parameter.getDimension(); i < len; ++i) {
            double value = parameter.getParameterValue(i);
            if (value < parameter.getBounds().getLowerLimit(i) ||
                    value > parameter.getBounds().getUpperLimit(i)) {
                throw new IllegalArgumentException("Parameter '" + parameter.getId() + "' is out-of-bounds");
            }
        }
    }

    private WrappedVector setupPreconditionedMatrix() {
        double[] masses = new double[parameter.getDimension()];
        Arrays.fill(masses, 1.0); // TODO
        return new WrappedVector.Raw(masses);
    }

    private final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    private final Parameter parameter;
    private final NormalDistribution drawDistribution;

    private double totalTravelTime;
    private double randomTimeWidth;

    private final WrappedVector masses;
}
