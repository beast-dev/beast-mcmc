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
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public abstract class AbstractParticleOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final boolean CHECK_MATRIX_ILL_CONDITIONED = false;

    AbstractParticleOperator(GradientWrtParameterProvider gradientProvider,
                             PrecisionMatrixVectorProductProvider multiplicationProvider,
                             PrecisionColumnProvider columnProvider,
                             double weight, Options runtimeOptions, Parameter mask) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.columnProvider = columnProvider;
        this.parameter = gradientProvider.getParameter();
        this.mask = mask;
        this.maskVector = mask.getParameterValues();

        this.runtimeOptions = runtimeOptions;
        this.preconditioning = setupPreconditioning();

        setWeight(weight);
        setMissingDataMask();
        checkParameterBounds(parameter);
    }

    private void setMissingDataMask() {

        int dim = parameter.getDimension();
        missingDataMask = new boolean[dim];
        assert (dim == parameter.getBounds().getBoundsDimension());

        for (int i = 0; i < dim; ++i) {

            missingDataMask[i] = (parameter.getBounds().getUpperLimit(i) == Double.POSITIVE_INFINITY &&
                    parameter.getBounds().getLowerLimit(i) == Double.NEGATIVE_INFINITY);
        }
    }

    @Override
    public double doOperation() {

        if (shouldUpdatePreconditioning()) {
            preconditioning = setupPreconditioning();
        }

        WrappedVector position = getInitialPosition();

        double hastingsRatio = integrateTrajectory(position);

        setParameter(position, parameter);

        if (CHECK_MATRIX_ILL_CONDITIONED & getCount() % 100 == 0) {
            productProvider.getTimeScaleEigen();
        }

        return hastingsRatio;
    }

//    @Override
//    public String getPerformanceSuggestion() {
//        return null;
//    }

    abstract double integrateTrajectory(WrappedVector position);

    double drawTotalTravelTime() {
        double randomFraction = 1.0 + runtimeOptions.randomTimeWidth * (MathUtils.nextDouble() - 0.5);
        return preconditioning.totalTravelTime * randomFraction;
    }

    static void updateGradient(WrappedVector gradient, double time, WrappedVector action) {

        final double[] g = gradient.getBuffer();
        final double[] a = action.getBuffer();

        for (int i = 0, len = g.length; i < len; ++i) {
            g[i] -= time * a[i];
        }
    }

    static void updatePosition(WrappedVector position, WrappedVector velocity, double time) {

        final double[] p = position.getBuffer();
        final double[] v = velocity.getBuffer();

        for (int i = 0, len = p.length; i < len; ++i) {
            p[i] += time * v[i];
        }
    }

    WrappedVector getInitialGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();

        if (mask != null) {
            applyMask(gradient);
        }

        return new WrappedVector.Raw(gradient);
    }

    @SuppressWarnings("unused")
    void applyMask(WrappedVector vector) {
        applyMask(vector.getBuffer());
    }

    void applyMask(double[] vector) {

        assert (vector.length == mask.getDimension());

        for (int i = 0, len = vector.length; i < len; ++i) {
            vector[i] *= maskVector[i];
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

    private WrappedVector getPrecisionColumn(int index) {

        double[] precisionColumn = columnProvider.getColumn(index);

        if (mask != null) {
            applyMask(precisionColumn);
        }

        return new WrappedVector.Raw(precisionColumn);
    }

    void updateAction(WrappedVector action, ReadableVector velocity, int eventIndex) {

        WrappedVector column = getPrecisionColumn(eventIndex);

        final double[] a = action.getBuffer();
        final double[] c = column.getBuffer();

        final double twoV = 2 * velocity.get(eventIndex);

        for (int i = 0, len = a.length; i < len; ++i) {
            a[i] += twoV * c[i];
        }

        if (mask != null) {
            applyMask(a);
        }
    }

    boolean headingTowardsBoundary(double position, double velocity, int positionIndex) {

        if (missingDataMask[positionIndex]) {
            return false;
        } else {
            return position * velocity < 0.0;
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

    private Preconditioning setupPreconditioning() {

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

    private boolean shouldUpdatePreconditioning() {
        return runtimeOptions.preconditioningUpdateFrequency > 0
                && (getCount() % runtimeOptions.preconditioningUpdateFrequency == 0);
    }

    protected class MinimumTravelInformation {

        final double time;
        final int index;

        MinimumTravelInformation(double minTime, int minIndex) {
            this.time = minTime;
            this.index = minIndex;
        }

        public boolean equals(Object obj) {

            if (obj instanceof MinimumTravelInformation) {
                MinimumTravelInformation rhs = (MinimumTravelInformation) obj;
                return this.time == rhs.time && this.index == rhs.index;
            }
            return false;
        }

        public String toString() {
            return "time = " + time + " @ " + index;
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

    enum Type {
        NONE,
        BOUNDARY,
        GRADIENT
    }


    private final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    private final PrecisionColumnProvider columnProvider;
    private final Parameter parameter;
    private final Options runtimeOptions;
    final Parameter mask;
    private final double[] maskVector;

    Preconditioning preconditioning;
    private boolean[] missingDataMask;
}
