/*
 * GumbelDistributionModel.java
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
import dr.inferencexml.distribution.GumbelDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.GumbelDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for (minimum-oriented) Gumbel type I distributed data.
 * See dr.math.distributions.GumbelDistribution for the exact orientation implemented:
 * if R ~ Exponential(mean = theta), then log(R) ~ GumbelDistributionModel(location = log(theta), scale = 1).
 *
 * @author Filippo Monti
 */

public class GumbelDistributionModel extends AbstractModel implements
        ParametricDistributionModel, GradientProvider, HessianProvider {

    public GumbelDistributionModel(Parameter locationParameter, Parameter scaleParameter) {

        super(GumbelDistributionModelParser.GUMBEL_DISTRIBUTION_MODEL);

        this.locationParameter = locationParameter;
        addVariable(this.locationParameter);
        this.locationParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        this.scaleParameter = scaleParameter;
        addVariable(this.scaleParameter);
        this.scaleParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public Parameter getLocationParameter() {
        return locationParameter;
    }

    public Parameter getScaleParameter() {
        return scaleParameter;
    }

    public final double getLocation() {
        return locationParameter.getValue(0);
    }

    public final double getScale() {
        return scaleParameter.getValue(0);
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return GumbelDistribution.pdf(x, getLocation(), getScale());
    }

    public double logPdf(double x) {
        return GumbelDistribution.logPdf(x, getLocation(), getScale());
    }

    public double cdf(double x) {
        return GumbelDistribution.cdf(x, getLocation(), getScale());
    }

    public double quantile(double y) {
        return GumbelDistribution.quantile(y, getLocation(), getScale());
    }

    public double mean() {
        return GumbelDistribution.mean(getLocation(), getScale());
    }

    public double variance() {
        return GumbelDistribution.variance(getLocation(), getScale());
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    @Override
    public int getDimension() { return 1; }

    @Override
    public double[] getDiagonalHessianLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = GumbelDistribution.hessianLogPdf(x[i], getLocation(), getScale());
        }
        return result;
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        double[] diagonalHessian = getDiagonalHessianLogDensity(x);
        double[][] result = new double[diagonalHessian.length][diagonalHessian.length];
        for (int i = 0; i < diagonalHessian.length; i++) {
            result[i][i] = diagonalHessian[i];
        }
        return result;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = GumbelDistribution.gradLogPdf(x[i], getLocation(), getScale());
        }
        return result;
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return locationParameter;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Parameter locationParameter;
    private final Parameter scaleParameter;
}
