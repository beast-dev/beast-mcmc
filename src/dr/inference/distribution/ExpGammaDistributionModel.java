/*
 * ExpGammaDistributionModel.java
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
import dr.math.distributions.GammaDistribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.apache.commons.math.special.Gamma;

/**
 * If X ~ Gamma(shape,scale), this is the distribution on Y = ln(X)
 *
 * @author Andy Magee
 * @author Yucai Shao
 */

public class ExpGammaDistributionModel extends AbstractModel
        implements ParametricDistributionModel, GradientProvider, HessianProvider, PriorPreconditioningProvider {

    public static final String EXP_GAMMA_DISTRIBUTION_MODEL = "ExpGammaDistributionModel";

    /**
     * Construct a gamma distribution model.
     */
    public ExpGammaDistributionModel(Variable<Double> shape, Variable<Double> scale) {

        super(EXP_GAMMA_DISTRIBUTION_MODEL);

        this.shape = shape;
        addVariable(shape);
        shape.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.scale = scale;
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return Math.exp(x) * GammaDistribution.pdf(Math.exp(x), getShape(), getScale());
    }

    public double logPdf(double x) {
        return x + GammaDistribution.logPdf(Math.exp(x), getShape(), getScale());
    }

    public double cdf(double x) {
        return GammaDistribution.cdf(Math.exp(x), getShape(), getScale());
    }

    public double quantile(double y) {
        try {
            return (new GammaDistributionImpl(getShape(), getScale())).inverseCumulativeProbability(Math.exp(y));
        } catch (MathException e) {
            return Double.NaN;
        }
    }

    public double mean() {
        return Gamma.digamma(getShape()) + Math.log(getScale());
    }

    public double variance() {
        return Gamma.trigamma(getShape());
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

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting


    @Override
    public double getStandardDeviation(int index) {
        return Math.sqrt(variance());
    }

    @Override
    public int getDimension() {
        return 1;
    }

    public double gradGammaPdf(double x, double shape, double scale) {
        return (Math.pow(x, shape - 2.0) * Math.exp(-x / scale) * (shape - 1) * scale - x) / scale;
    }

    public double hessianGammaPdf(double x, double shape, double scale) {
        double numerator = Math.pow(x, shape - 3.0) * Math.exp(-x / scale);
        numerator *= ((shape * shape - 3.0 * shape + 2.0) * scale * scale - 2.0 * (shape - 1.0) * scale * x  + x * x);
        return numerator / (scale * scale);
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {

        double[] x = GradientProvider.toDoubleArray(obj);

        double[] result = new double[x.length];
        double shape = getShape();
        double scale = getScale();
        for (int i = 0; i < x.length; ++i) {
            double expX = Math.exp(x[i]);
            result[i] = expX * gradGammaPdf(expX, shape, scale) / GammaDistribution.pdf(expX, shape, scale);
        }
        return result;
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object obj) {

        double[] x = GradientProvider.toDoubleArray(obj);

        double[] result = new double[x.length];
        double shape = getShape();
        double scale = getScale();

        for (int i = 0; i < x.length; ++i) {
            double expX = Math.exp(x[i]);
            double fExpX = GammaDistribution.pdf(expX, shape, scale);
            double fPrimeExpX = gradGammaPdf(expX, shape, scale);
            double fPrimePrimeExpX = hessianGammaPdf(expX, shape, scale);
            result[i] = expX * (fExpX * (expX * fPrimePrimeExpX + fPrimeExpX) - expX * fPrimeExpX * fPrimeExpX) / (fExpX * fExpX);
        }
        return result;
    }

    @Override
    public double[][] getHessianLogDensity(Object obj) {
        return HessianProvider.expandDiagonals(getDiagonalHessianLogDensity(obj));
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    public double getShape() {
        return shape.getValue(0);
    }

    public double getScale() { return scale.getValue(0); }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Variable<Double> shape;
    private final Variable<Double> scale;

}

