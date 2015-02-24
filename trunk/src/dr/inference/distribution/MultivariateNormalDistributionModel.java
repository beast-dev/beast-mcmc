/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;

/**
 * A class that acts as a model for multivariate normally distributed data.
 *
 * @author Marc Suchard
 * @author Max Tolkoff
 */

public class MultivariateNormalDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel, GaussianProcessRandomGenerator {

    public MultivariateNormalDistributionModel(Parameter meanParameter, MatrixParameter precParameter) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
        this.mean = meanParameter;
        addVariable(meanParameter);
        meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                meanParameter.getDimension()));

        this.precision = precParameter;
        addVariable(precParameter);

        distribution = createNewDistribution();
        distributionKnown = true;
    }

    public MatrixParameter getPrecisionMatrixParameter() {
        return precision;
    }

    public Parameter getMeanParameter() {
        return mean;
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
        return precision.getParameterAsMatrix();
    }

    public double[] getMean() {
        return mean.getParameterValues();
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

    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    private MultivariateNormalDistribution createNewDistribution() {
        return new MultivariateNormalDistribution(getMean(), getScaleMatrix());
    }

    private final Parameter mean;
    private final MatrixParameter precision;
    private MultivariateNormalDistribution distribution;
    private MultivariateNormalDistribution storedDistribution;

    private boolean distributionKnown;
    private boolean storedDistributionKnown;

    // RandomGenerator interface
    public double[] nextRandom() {
        checkDistribution();
        return distribution.nextMultivariateNormal();
    }

    public double logPdf(Object x) {
        checkDistribution();
        return distribution.logPdf(x);
    }
}
