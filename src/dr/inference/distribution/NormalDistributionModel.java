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
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.NormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for normally distributed data.
 *
 * @author Alexei Drummond
 * @version $Id: NormalDistributionModel.java,v 1.6 2005/05/24 20:25:59 rambaut Exp $
 */

public class NormalDistributionModel extends AbstractModel implements ParametricDistributionModel,
        GaussianProcessRandomGenerator, GradientProvider {
    /**
     * Constructor.
     */
    public NormalDistributionModel(Variable<Double> mean, Variable<Double> scale) {

        super(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.mean = mean;
        this.scale = scale;
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public NormalDistributionModel(Parameter meanParameter, Parameter scale, boolean isPrecision) {
        super(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
        this.hasPrecision = isPrecision;
        this.mean = meanParameter;
        addVariable(meanParameter);
        meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        if (isPrecision) {
            this.precision = scale;
            this.scale = null;  // todo why not keep the name scale to avoid confusion?? For whom?
        } else {
            this.scale = scale;
        }
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getStdev() {
        return getScale();
    }

    public double getScale() {
        if (hasPrecision)
            return 1.0 / Math.sqrt(precision.getValue(0));
        return scale.getValue(0);
    }

    public Variable<Double> getMean() {
        return mean;
    }

    public Variable<Double> getPrecision() {
        if (hasPrecision)
            return precision;
        return null;
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return NormalDistribution.pdf(x, mean(), getScale());
    }

    public double logPdf(double x) {
        return NormalDistribution.logPdf(x, mean(), getScale());
    }

    public double cdf(double x) {
        return NormalDistribution.cdf(x, mean(), getScale());
    }

    public double quantile(double y) {
        return NormalDistribution.quantile(y, mean(), getScale());
    }

    public double mean() {
        return mean.getValue(0);
    }

    public double variance() {
        if (hasPrecision)
            return 1.0 / precision.getValue(0);
        double sd = scale.getValue(0);
        return sd * sd;
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

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Variable<Double> mean;
    private final Variable<Double> scale;
    private Variable<Double> precision;
    private boolean hasPrecision = false;

    public Object nextRandom() {
        double eps = MathUtils.nextGaussian();
        eps *= getScale();
        eps += mean();
        return eps;
    }

    public double logPdf(Object x) {
        double v = (Double) x;
        return logPdf(v);
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() { return 1; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        double[] result = new double[1];
        result[0] = NormalDistribution.gradLogPdf((Double) x, mean(), getScale());
        return result;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        double p = hasPrecision ?
                precision.getValue(0) :
                scale.getValue(0) * scale.getValue(0);
        return new double[][]{{p}};
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
        return mean;
    }

}
