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

import dr.evolution.alignment.PatternList;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.FastMatrixParameter;
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

    AbstractParticleOperator(GradientWrtParameterProvider gradientProvider,
                             PrecisionMatrixVectorProductProvider multiplicationProvider,
                             double weight, Options runtimeOptions, Parameter mask, PatternList patternList) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.parameter = gradientProvider.getParameter();
        this.mask = mask;

        this.runtimeOptions = runtimeOptions;
        this.preconditioning = setupPreconditioning();
        this.patternList = patternList;

        setWeight(weight);
        setMissingDataMask();
        checkParameterBounds(parameter);
    }

    private void setMissingDataMask() {

        missingDataMask = new double[parameter.getDimension()];

        int traitDim = ((FastMatrixParameter) parameter).getRowDimension();//todo: better way to get these two dimension? "row/col dimension" isn't explicit.
        int taxaDim = ((FastMatrixParameter) parameter).getColumnDimension();

        for (int i = 0; i < taxaDim; ++i) {

            String taxonName = ((FastMatrixParameter) parameter).getParameter(i).getId();
            int indexInPatternlist = patternList.getTaxonIndex(taxonName);

            for (int j = 0; j < traitDim; ++j) {
                int t = patternList.getPattern(j)[indexInPatternlist];
                missingDataMask[i * traitDim + j] = (t > 1) ? 1 : 0; //now value = 0 in the mask means missing observation;
            }
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

    static void updateGradient(WrappedVector gradient, double time, ReadableVector action) {
        for (int i = 0, len = gradient.getDim(); i < len; ++i) {
            gradient.set(i, gradient.get(i) - time * action.get(i));
        }
    }

    static void updatePosition(WrappedVector position, WrappedVector velocity, double time) {
        for (int i = 0, len = position.getDim(); i < len; ++i) {
            position.set(i, position.get(i) + time * velocity.get(i));
        }
    }

    WrappedVector getInitialGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();

        if (mask != null) {
            applyMask(gradient);
        }

        return new WrappedVector.Raw(gradient);
    }

    void applyMask(WrappedVector vector) {

        assert (vector.getDim() == mask.getDimension());

        for (int i = 0, dim = vector.getDim(); i < dim; ++i) {
            vector.set(i, vector.get(i) * mask.getParameterValue(i));
        }
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

    boolean headingTowardsBoundary(double position, double velocity, int positionIndex) {

        if (missingDataMask[positionIndex] == 1.0) {
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

    private final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    private final Parameter parameter;
    private final Options runtimeOptions;
    final Parameter mask;

    Preconditioning preconditioning;
    PatternList patternList;
    private double[] missingDataMask;
}
