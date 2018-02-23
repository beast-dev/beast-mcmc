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
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;
import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public abstract class AbstractParticleOperator extends SimpleMCMCOperator implements GibbsOperator {

    AbstractParticleOperator(GradientWrtParameterProvider gradientProvider,
                                    PrecisionMatrixVectorProductProvider multiplicationProvider,
                                    double weight, Options runtimeOptions, Parameter mask) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.parameter = gradientProvider.getParameter();
        this.mask = mask;

        this.runtimeOptions = runtimeOptions;
        this.preconditioning = setupPreconditioning();

        setWeight(weight);
        checkParameterBounds(parameter);
    }

    @Override
    public double doOperation() {

        if (shouldUpdatePreconditioning()) {
            preconditioning = setupPreconditioning();
        }

        WrappedVector position = getInitialPosition();

        double hastingsRatio = integrateTrajectory(position);

        setParameter(position, parameter);

        return hastingsRatio;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    abstract double integrateTrajectory(WrappedVector position);

    double drawTotalTravelTime() {
        double randomFraction = 1.0 + runtimeOptions.randomTimeWidth * (MathUtils.nextDouble() - 0.5);
        return preconditioning.totalTravelTime * randomFraction;
    }

    void updateVelocity(WrappedVector velocity, WrappedVector gradient, ReadableVector mass) {

        ReadableVector gDivM = new ReadableVector.Quotient(gradient, mass);

        double vg = innerProduct(velocity, gradient);
        double ggDivM = innerProduct(gradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i, velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
        }
    }

    void updateGradient(WrappedVector gradient, double time, ReadableVector Phi_v) {
        for (int i = 0, len = gradient.getDim(); i < len; ++i) {
            gradient.set(i, gradient.get(i) - time * Phi_v.get(i));
        }
    }

    void updatePosition(WrappedVector position, WrappedVector velocity, double time) {
        for (int i = 0, len = position.getDim(); i < len; ++i) {
            position.set(i, position.get(i) + time * velocity.get(i));
        }
    }

    @SuppressWarnings("all")
    protected double getBounceTime(double v_phi_v, double v_phi_x, double u_min) {
        double a = v_phi_v / 2;
        double b = v_phi_x;
        double c = u_min - MathUtils.nextExponential(1);
        return (-b + Math.sqrt(b * b - 4 * a * c)) / 2 / a;
    }
    
    WrappedVector getInitialGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();

        if (mask != null) {
            applyMask(gradient);
        }

        return new WrappedVector.Raw(gradient);
    }

    void applyMask(double[] vector) {

        assert (vector.length == mask.getDimension());

        for (int i = 0, len = vector.length; i < len; ++i) {
            vector[i] *= mask.getParameterValue(i);
        }
    }

    WrappedVector getPrecisionProduct(ReadableVector velocity) {

        setParameter(velocity, parameter);

        double[] product = productProvider.getProduct(parameter);

        if (mask != null) {
            applyMask(product);
        }

        return new WrappedVector.Raw(product);
    }

    WrappedVector drawInitialVelocity() {

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

    MinimumTravelInformation getTimeToBoundary(ReadableVector position, ReadableVector velocity) {

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

    boolean headingAwayFromBoundary(double position, double velocity) {
        return position * velocity < 0.0;
    }

    WrappedVector getInitialPosition() {
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

    Preconditioning setupPreconditioning() {

        double[] mass = new double[parameter.getDimension()];
        Arrays.fill(mass, 1.0);

        // TODO Should use:
        productProvider.getMassVector();
        double time = productProvider.getTimeScale();

        return new Preconditioning(
                new WrappedVector.Raw(mass),
                time
        );
    }

    boolean shouldUpdatePreconditioning() {
        return runtimeOptions.preconditioningUpdateFrequency > 0
                && (getCount() % runtimeOptions.preconditioningUpdateFrequency == 0);
    }

    protected class MinimumTravelInformation {

        final double minTime;
        final int minIndex;

        private MinimumTravelInformation(double minTime, int minIndex) {
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

    public static class Options {

        final double randomTimeWidth;
        final int preconditioningUpdateFrequency;

        public Options(double randomTimeWidth, int preconditioningUpdateFrequency) {
            this.randomTimeWidth = randomTimeWidth;
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
        }
    }

    protected class Preconditioning {

        final WrappedVector mass;
        final double totalTravelTime;

        private Preconditioning(WrappedVector mass, double totalTravelTime) {
            this.mass = mass;
            this.totalTravelTime = totalTravelTime;
        }
    }

    private final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    final Parameter parameter;
    private final Options runtimeOptions;
    final Parameter mask;

    Preconditioning preconditioning;
}
