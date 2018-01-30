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
                                  double weight, Options runtimeOptions, Parameter mask) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.parameter = gradientProvider.getParameter();
        this.drawDistribution = new NormalDistribution(0, 1);

        this.runtimeOptions = runtimeOptions;
        this.preconditioning = setupPreconditioning();

        this.mask = mask;

        setWeight(weight);
        checkParameterBounds(parameter);
    }

    @Override
    public double doOperation() {

        if (shouldUpdatePreconditioning()) {
            preconditioning = setupPreconditioning();
        }

        WrappedVector position = getInitialPosition();
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

        setParameter(position);

        return 0.0;
    }

    private double drawTotalTravelTime() {
        double randomFraction = 1.0 + runtimeOptions.randomTimeWidth * (MathUtils.nextDouble() - 0.5);
        return preconditioning.totalTravelTime * randomFraction;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy particle operator";
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

    private void setParameter(ReadableVector position) {
        for (int j = 0, dim = position.getDim(); j < dim; ++j) {
            parameter.setParameterValueQuietly(j, position.get(j));
        }
        parameter.fireParameterChangedEvent();
    }

    private void updateVelocity(WrappedVector velocity, WrappedVector gradient, ReadableVector mass) {

        ReadableVector gDivM = new ReadableVector.Quotient(gradient, mass);

        double vg = innerProduct(velocity, gradient);
        double ggDivM = innerProduct(gradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            if (mask.getParameterValue(i) == 0.0) return; // leave the velocity to be zero if masked out
            else {
                velocity.set(i, velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
            }
        }
    }

    private void updateGradient(WrappedVector gradient, double time, ReadableVector Phi_v) {
        for (int i = 0, len = gradient.getDim(); i < len; ++i) {
            gradient.set(i, gradient.get(i) - time * Phi_v.get(i));
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
    
    private WrappedVector getInitialGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();
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

        ReadableVector mass = preconditioning.mass;
        double[] velocity = new double[mass.getDim()];

        for (int i = 0, len = velocity.length; i < len; i++) {
            if (mask != null){
                velocity[i] = mask.getParameterValue(i)*(Double) drawDistribution.nextRandom() / Math.sqrt(mass.get(i));
            } else {
                velocity[i] = (Double) drawDistribution.nextRandom() / Math.sqrt(mass.get(i));
            }
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

    private Preconditioning setupPreconditioning() {

        double[] mass = new double[parameter.getDimension()];
        Arrays.fill(mass, 1.0);

        // TODO Should use:
        productProvider.getMassVector();
        productProvider.getTimeScale();

        return new Preconditioning(
                new WrappedVector.Raw(mass),
                0.05
        );
    }

    private boolean shouldUpdatePreconditioning() {
        return runtimeOptions.preconditioningUpdateFrequency > 0
                && (getCount() % runtimeOptions.preconditioningUpdateFrequency == 0);
    }

    private class MinimumTravelInformation {

        double minTime;
        int minIndex;

        private MinimumTravelInformation(double minTime, int minIndex) {
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

    public static class Options {

        double randomTimeWidth;
        int preconditioningUpdateFrequency;

        public Options(double randomTimeWidth, int preconditioningUpdateFrequency) {
            this.randomTimeWidth = randomTimeWidth;
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
        }
    }

    private class Preconditioning {
        WrappedVector mass;
        double totalTravelTime;

        private Preconditioning(WrappedVector mass, double totalTravelTime) {
            this.mass = mass;
            this.totalTravelTime = totalTravelTime;
        }
    }

    private final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    private final Parameter parameter;
    private final NormalDistribution drawDistribution;
    private final Options runtimeOptions;
    private Preconditioning preconditioning;

    private final Parameter mask;
}
