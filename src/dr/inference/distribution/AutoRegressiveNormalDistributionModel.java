/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;

import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.AutoRegressiveNormalDistribution;
import dr.math.distributions.GaussianProcessRandomGenerator;

/**
 * A class that acts as a model for auto-regressive normally distributed data.
 *
 * @author Marc Suchard
 */

public class AutoRegressiveNormalDistributionModel extends AbstractModel implements
        ParametricMultivariateDistributionModel, GaussianProcessRandomGenerator, GradientProvider, HessianProvider {

    public AutoRegressiveNormalDistributionModel(int dim,
                                                 Parameter marginalParameter,
                                                 Parameter decayParameter) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;

        this.marginal = marginalParameter;
        addVariable(marginalParameter);

        this.decay = decayParameter;
        addVariable(decayParameter);

        distribution = createNewDistribution();
        distributionKnown = true;
    }

    public Parameter getMeanParameter() {
        return null;
    }

    // *****************************************************************
    // Interface MultivariateDistribution
    // *****************************************************************

    private void checkDistribution() {
        if (!distributionKnown) {
            distribution = createNewDistribution();
            distributionKnown = true;
        }
    }

    public double logPdf(double[] x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getPrecisionColumn(int index) {
        checkDistribution();
        return distribution.getPrecisionColumn(index);
    }

    public double[] getMean() {
        return new double[dim];
    }

    public String getType() {
        return distribution.getType();
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    public Likelihood getLikelihood() {
        return null;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        distributionKnown = false;
    }

    protected void storeState() {
        storedDistribution = distribution;
        storedDistributionKnown = distributionKnown;
    }

    protected void restoreState() {
        distributionKnown = storedDistributionKnown;
        distribution = storedDistribution;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    private AutoRegressiveNormalDistribution createNewDistribution() {
        return new AutoRegressiveNormalDistribution(getDimension(),
                marginal.getParameterValue(0), decay.getParameterValue(0));
    }

    private final int dim;
    private final Parameter marginal;
    private final Parameter decay;

    private AutoRegressiveNormalDistribution distribution;
    private AutoRegressiveNormalDistribution storedDistribution;

    private boolean distributionKnown;
    private boolean storedDistributionKnown;

    // RandomGenerator interface
    public double[] nextRandom() {
        checkDistribution();
        throw new RuntimeException("Not yet implemented");
    }

    public double logPdf(Object x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    // GradientWrtParameterProvider interface
    @Override
    public double[] getGradientLogDensity(Object x) {
        checkDistribution();
        return distribution.getGradientLogDensity(x);
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        checkDistribution();
        return distribution.getDiagonalHessianLogDensity(x);
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        checkDistribution();
        return distribution.getHessianLogDensity(x);
    }

    public double[] getPrecisionVectorProduct(double[] x) {
        checkDistribution();
        return distribution.getPrecisionVectorProduct(x);
    }

    public double[] getDiagonal() {
        checkDistribution();
        return distribution.getDiagonal();
    }
}
