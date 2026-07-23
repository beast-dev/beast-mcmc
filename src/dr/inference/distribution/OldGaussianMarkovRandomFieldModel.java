/*
 * OldGaussianMarkovRandomFieldModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.distribution;

import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.OldGaussianMarkovRandomField;
import dr.math.distributions.GaussianProcessRandomGenerator;


/**
 * A class that acts as a model for gaussian random walk
 *
 * @author Marc Suchard
 * Pratyusa Datta
 */

public class OldGaussianMarkovRandomFieldModel extends AbstractModel implements
        ParametricMultivariateDistributionModel, GaussianProcessRandomGenerator, GradientProvider, HessianProvider {

    public OldGaussianMarkovRandomFieldModel(int dim,
                                             Parameter incrementPrecisionParameter,
                                             Parameter startParameter) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;

        this.incrementPrecision = incrementPrecisionParameter;
        addVariable(incrementPrecisionParameter);

        this.start = startParameter;
        addVariable(startParameter);

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
        return distribution.getScaleMatrix();
    }


    public double[] getMean() {
        return distribution.getMean();
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
        return distribution.getPrecisionMatrix();
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    public Parameter getIncrementPrecision() { return incrementPrecision; }

    public Parameter getStart() { return start; }
    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    private OldGaussianMarkovRandomField createNewDistribution() {
        return new OldGaussianMarkovRandomField(getDimension(),
                incrementPrecision.getParameterValue(0), start.getParameterValue(0));
    }


    private final int dim;
    private final Parameter incrementPrecision;
    private final Parameter start;

    private OldGaussianMarkovRandomField distribution;
    private OldGaussianMarkovRandomField storedDistribution;

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


}
