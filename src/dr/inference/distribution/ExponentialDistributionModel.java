/*
 * ExponentialDistributionModel.java
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
import dr.math.distributions.ExponentialDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for exponentially distributed data.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ExponentialDistributionModel.java,v 1.12 2005/05/24 20:25:59 rambaut Exp $
 */

public class ExponentialDistributionModel extends AbstractModel implements
        ParametricDistributionModel, GradientProvider, HessianProvider {

    public static final String EXPONENTIAL_DISTRIBUTION_MODEL = "exponentialDistributionModel";

    /**
     * Constructor.
     */
    public ExponentialDistributionModel(Variable<Double> mean) {

        this(mean, 0.0);
    }


    /**
     * Constructor.
     */
    public ExponentialDistributionModel(Variable<Double> mean, double offset) {

        super(EXPONENTIAL_DISTRIBUTION_MODEL);

        this.mean = mean;
        this.offset = offset;

        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x < offset) return 0.0;
        return ExponentialDistribution.pdf(x - offset, 1.0 / getMean());
    }

    public double logPdf(double x) {
        if (x < offset) return Double.NEGATIVE_INFINITY;
        return ExponentialDistribution.logPdf(x - offset, 1.0 / getMean());
    }

    public double cdf(double x) {
        if (x < offset) return 0.0;
        return ExponentialDistribution.cdf(x - offset, 1.0 / getMean());
    }

    public double quantile(double y) {
        return ExponentialDistribution.quantile(y, 1.0 / getMean()) + offset;
    }

    public double mean() {
        return ExponentialDistribution.mean(1.0 / getMean()) + offset;
    }

    public double variance() {
        return ExponentialDistribution.variance(1.0 / getMean());
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return offset;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
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
        return mean;
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
    // Private methods
    // **************************************************************

    private double getMean() {
        return mean.getValue(0);
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Variable<Double> mean = null;
    private double offset = 0.0;

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object obj) {
        return getDerivativeLogDensity(obj, DerivativeType.DIAGONAL_HESSIAN);
    }

    @Override
    public double[][] getHessianLogDensity(Object obj) {

        double[] diagonalHessian = getDiagonalHessianLogDensity(obj);
        double[][] result = new double[diagonalHessian.length][diagonalHessian.length];
        for (int i = 0; i < diagonalHessian.length; i++) {
            result[i][i] = diagonalHessian[i];
        }
        return result;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        return getDerivativeLogDensity(obj, DerivativeType.GRADIENT);
    }

    private double[] getDerivativeLogDensity(Object obj, DerivativeType derivativeType) {

        double[] x = GradientProvider.toDoubleArray(obj);

        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = derivativeType.getDerivativeLogPdf(x[i] - offset, 1.0 / getMean());
        }
        return result;
    }

    private enum DerivativeType {
        GRADIENT("gradient") {
            @Override
            public double getDerivativeLogPdf(double x, double lambda) {
                return ExponentialDistribution.gradLogPdf(x, lambda);
            }
        },
        DIAGONAL_HESSIAN("diagonalHessian") {
            @Override
            public double getDerivativeLogPdf(double x, double lambda) {
                return ExponentialDistribution.hessianLogPdf(x, lambda);
            }
        };

        private String type;

        DerivativeType(String type) {
            this.type = type;
        }

        public abstract double getDerivativeLogPdf(double x, double lambda);
    }
}

