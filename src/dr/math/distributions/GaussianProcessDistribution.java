/*
 * GaussianMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;

import dr.inference.distribution.LogGaussianProcessModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessDistribution extends RandomFieldDistribution {

    public static final String TYPE = "GaussianProcess";

    private final int dim;
    private final Parameter mean;
    private final LogGaussianProcessModel.GaussianProcessKernel kernel;
    private final RandomField.WeightProvider weightProvider;

    public GaussianProcessDistribution(String name,
                                       int dim,
                                       Parameter mean,
                                       LogGaussianProcessModel.GaussianProcessKernel kernel,
                                       RandomField.WeightProvider weightProvider) {
        super(name);

        this.dim = dim;
        this.mean = mean;
        this.kernel = kernel;
        this.weightProvider = weightProvider;

        addVariable(mean);
        addModel(kernel);

        if (weightProvider != null) {
            addModel(weightProvider);
        }
    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getIncrement(int i, Parameter field) {
        double[] mean = getMean();
        return (field.getParameterValue(i) - mean[i]) - (field.getParameterValue(i + 1) - mean[i + 1]);
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        throw new RuntimeException("Unknown parameter");
    }

    public String getType() {
        return TYPE;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(double[] x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() { }
}
