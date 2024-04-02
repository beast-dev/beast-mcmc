/*
 * BetaDistributionModel.java
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
import dr.math.distributions.BetaDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for beta distributed data.
 *
 * @author Marc A. Suchard
 * @author Andy Magee
 */

public class BetaDistributionModel extends AbstractModel implements ParametricDistributionModel, GradientProvider {

    public static final String BETA_DISTRIBUTION_MODEL = "betaDistributionModel";

    public BetaDistributionModel(Variable<Double> alpha, Variable<Double> beta) {
        this(alpha, beta, 0.0, 1.0);
    }


    /**
     * Constructor.
     */
    public BetaDistributionModel(Variable<Double> alpha, Variable<Double> beta, double offset, double length) {

        super(BETA_DISTRIBUTION_MODEL);

        this.alpha = alpha;
        this.beta = beta;
        this.length = length;
        this.offset = offset;

        addVariable(alpha);
        alpha.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        addVariable(beta);
        beta.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        recomputeBetaDistribution();
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        double xScaled = getXScaled(x);
        if (xScaled < 0.0 || xScaled > 1.0) return 0.0;

        return betaDistribution.pdf(xScaled);
    }

    public double logPdf(double x) {
        double xScaled = getXScaled(x);
        if (xScaled < 0.0 || xScaled > 1.0) return Double.NEGATIVE_INFINITY;

        return betaDistribution.logPdf(xScaled);
    }

    public double cdf(double x) {
        if (x < offset) return 0.0;
        return betaDistribution.cdf(getXScaled(x));
    }

    public double quantile(double y) {
        return betaDistribution.quantile(getXScaled(y)) * length + offset;
    }

    public double mean() {
        return betaDistribution.mean() * length + offset;
    }

    public double variance() {
        return betaDistribution.variance() * length * length;
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            double xScale = (x - offset) / length;
            return pdf(xScale);
        }

        public final double getLowerBound() {
            return offset;
        }

        public final double getUpperBound() {
            return offset + length;
        }
    };

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

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        recomputeBetaDistribution();
    }

    protected void storeState() {
        storedBetaDistribution = betaDistribution;
    }

    protected void restoreState() {
        betaDistribution = storedBetaDistribution;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // Private methods
    // **************************************************************

    private void recomputeBetaDistribution() {
        betaDistribution = new BetaDistribution(alpha.getValue(0), beta.getValue(0));
    }

    private double getXScaled(double x) {
        return (x - offset) / length;
    }


    // **************************************************************
    // GradientProvider implementation
    // **************************************************************

    public static double gradLogPdf(double x, double a, double b) {
        return (a - 1.0) / x - (b - 1.0) / (1.0 - x);
    }

    public static double scaledGradLogPdf(double x, double a, double b, double o, double l) {
        return (-(a + b - 2.0) * (o - x) - a * l  + l) / ((o - x) * (l + o - x));
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = scaledGradLogPdf(x[i], alpha.getValue(0), beta.getValue(0), offset, length);
        }
        return result;
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Variable<Double> alpha = null;
    private Variable<Double> beta = null;
    private double offset = 0.0;
    private double length = 0.0;

    private BetaDistribution betaDistribution = null;
    private BetaDistribution storedBetaDistribution = null;

}

