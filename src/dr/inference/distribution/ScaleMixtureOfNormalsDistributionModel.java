/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.distributions.NormalDistribution;

/**
 * A model for scaled-mixture-of-normals distributed data.
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */

public class ScaleMixtureOfNormalsDistributionModel extends AbstractModel implements
        ParametricMultivariateDistributionModel, GaussianProcessRandomGenerator,
        GradientProvider, HessianProvider {

    public ScaleMixtureOfNormalsDistributionModel(Parameter globalScale, Parameter localScale) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.globalScale = globalScale;
        this.localScale = localScale;
        this.dim = localScale.getDimension();
        this.distribution = new NormalDistribution(0, 1);

        addVariable(globalScale);
        addVariable(localScale);
    }

    @Override
    public double logPdf(double[] x) {

        assert (x.length == dim);

        double sum = 0.0;
        for (int i = 0; i < x.length; ++i) {
            sum += distribution.logPdf(x[i] / getStandardDeviation(i));
        }

        return sum;
    }

    private double getStandardDeviation(int index) {
        return globalScale.getParameterValue(0) * localScale.getParameterValue(index);
    }

    public double[][] getScaleMatrix() {
        double[][] matrix = new double[dim][dim];

        for (int i = 0; i < dim; ++i) {
            matrix[i][i] = getStandardDeviation(i);
        }

        return matrix;
    }

    @Override
    public double[] getMean() {
        return new double[dim];
    }

    @Override
    public String getType() {
        return "Scaled-mixture-of-normals";
    }

    @Override
    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
       // no intermediates need to be recalculated...
    }

    @Override
    protected void storeState() {
       // Do nothing
    }

    @Override
    protected void restoreState() {
        // Do nothing
    }

    @Override
    protected void acceptState() {
    } // no additional state needs accepting

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[][] getPrecisionMatrix() {

        double[][] matrix = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            double sd = getStandardDeviation(i);
            matrix[i][i] = 1.0 / (sd * sd);
        }

        return matrix;
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    private final Parameter globalScale;
    private final Parameter localScale;
    private final int dim;
    private final NormalDistribution distribution;

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(Object x) {
        double[] vector = (double []) x;
        assert (vector.length == dim);

        return logPdf(vector);
    }

    @Override
    public double[] getGradientLogDensity(Object x) {

        double[] vector = (double[]) x;
        assert (vector.length == dim);

        double[] gradient = new double[dim];

        for (int i = 0; i < dim; ++i) {
            gradient[i] = NormalDistribution.gradLogPdf(vector[i], 0.0, getStandardDeviation(i));
        }

        return gradient;
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }
}
