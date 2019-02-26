/*
 * NormalDistributionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.UnivariateFunction;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class MarginalizedAlphaStableDistributionModel extends AbstractModel
        implements ParametricDistributionModel, GradientProvider {

    private final Parameter scale;
    private final Parameter alpha;

    public MarginalizedAlphaStableDistributionModel(String name,
                                                    Parameter scale,
                                                    Parameter alpha) {
        super(name);
        this.scale = scale;
        this.alpha = alpha;

        addVariable(scale);
        addVariable(alpha);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    @Override
    public double logPdf(double[] x) {
        double logPdf = 0.0;
        for (double value : x) {
            logPdf += logPdf(value);
        }
        return logPdf;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    @Override
    public double logPdf(double x) {
        return MarginalizedAlphaStableDistribution.logPdf(x,
                scale.getParameterValue(0), alpha.getParameterValue(0));
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x;
        if (obj instanceof double[]) {
            x = (double[]) obj;
        } else {
            x = new double[1];
            x[0] = (Double) obj;
        }

        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = MarginalizedAlphaStableDistribution.gradLogPdf(x[i],
                    scale.getParameterValue(0), alpha.getParameterValue(0));
        }
        return result;
    }
}
