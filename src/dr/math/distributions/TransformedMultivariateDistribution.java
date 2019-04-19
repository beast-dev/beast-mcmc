/*
 * TransformedMultivariateDistribution.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.GradientProvider;
import dr.util.Transform;

public class TransformedMultivariateDistribution implements MultivariateDistribution, GradientProvider {

    private MultivariateDistribution distribution;
    private Transform.MultivariateTransform transform;
    private Transform.MultivariateTransform inverseTransform;

    public TransformedMultivariateDistribution(MultivariateDistribution distribution, Transform.MultivariateTransform transform) {

        assert (distribution instanceof GradientProvider) : "The transformed distribution should be a gradient provider";

        this.distribution = distribution; // On un-constrained space
        this.transform = transform; // From constrained to un-constrained
        this.inverseTransform = new Transform.InverseMultivariate(transform);

    }

    @Override
    public double logPdf(double[] x) {
        return distribution.logPdf(transform.transform(x)) + transform.getLogJacobian(x);
    }

    @Override
    public String getType() {
        return "Transformed." + distribution.getType();
    }

    @Override
    public int getDimension() {
        return ((GradientProvider) distribution).getDimension();
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    private double[] gradLogPdf(double[] x) {
        double[] transformedValue = transform.transform(x);
        double[] gradient = ((GradientProvider) distribution).getGradientLogDensity(transformedValue);
        return updateGradientLogDensity(gradient, transformedValue);
    }

    private double[] updateGradientLogDensity(double[] gradient, double[] transformedValue) {
        return inverseTransform.updateGradientLogDensity(gradient, transformedValue, 0, transformedValue.length);
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }
}
