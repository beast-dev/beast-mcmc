/*
 * NormalDistributionModel.java
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
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.math.UnivariateFunction;

/**
 * @author Marc A Suchard
 */

public class CauchyDistributionModel extends AbstractModel implements ParametricDistributionModel, GradientProvider {

    private final Variable<Double> median;
    private final Variable<Double> scale;

    public CauchyDistributionModel(Variable<Double> median, Variable<Double> scale) {

        super(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.median = median;
        this.scale = scale;
        addVariable(median);
        median.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    private double getScale() {
        return scale.getValue(0);
    }

    private double getMedian() {
        return median.getValue(0);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return median;
    }

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);

        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = CauchyDistribution.gradLogPdf(x[i], getMedian(), getScale());
        }
        return result;
    }

    @Override
    public double pdf(double x) {
        return CauchyDistribution.pdf(x, getMedian(), getScale());
    }

    @Override
    public double logPdf(double x) {
        return CauchyDistribution.logPdf(x, getMedian(), getScale());
    }

    @Override
    public double cdf(double x) {
        return CauchyDistribution.cdf(x, getMedian(), getScale());
    }

    @Override
    public double quantile(double y) {
        return CauchyDistribution.quantile(y, getMedian(), getScale());
    }

    @Override
    public double mean() {
        return Double.NaN;
    }

    @Override
    public double variance() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        CauchyDistribution distribution = new CauchyDistribution(getMedian(), getScale());
        return distribution.getProbabilityDensityFunction();
    }
}
