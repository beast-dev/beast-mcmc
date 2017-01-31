/*
 * GammaDistributionModel.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.UnivariateFunction;
import dr.math.distributions.GammaDistribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for gamma distributed data.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GammaDistributionModel.java,v 1.6 2005/05/24 20:25:59 rambaut Exp $
 */

public class GammaDistributionModel extends AbstractModel implements ParametricDistributionModel {

    public enum GammaParameterizationType {
        ShapeScale,
        ShapeRate,
        ShapeMean,
        OneParameter
    }

    public static final String GAMMA_DISTRIBUTION_MODEL = "gammaDistributionModel";
    public static final String ONE_P_GAMMA_DISTRIBUTION_MODEL = "onePGammaDistributionModel";

    /**
     * Construct a gamma distribution model with a default shape scale parameterization.
     */
    public GammaDistributionModel(Variable<Double> shape, Variable<Double> scale) {
        this(GammaParameterizationType.ShapeScale, shape, scale, 0.0);
    }

    /**
     * Construct a one parameter gamma distribution model.
     */
    public GammaDistributionModel(Variable<Double> shape) {
        this(GammaParameterizationType.OneParameter, shape, null, 0.0);
    }


    /**
     * Construct a gamma distribution model.
     */
    public GammaDistributionModel(GammaParameterizationType parameterization, Variable<Double> shape, Variable<Double> parameter2, double offset) {

        super(GAMMA_DISTRIBUTION_MODEL);

        this.offset = offset;

        this.parameterization = parameterization;
        this.shape = shape;
        addVariable(shape);
        shape.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        switch (parameterization) {
            case ShapeScale:
                this.scale = parameter2;
                addVariable(scale);
                scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
                rate = null;
                mean = null;
                break;
            case ShapeRate:
                this.rate = parameter2;
                addVariable(rate);
                rate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
                scale = null;
                mean = null;
                break;
            case ShapeMean:
                this.mean = parameter2;
                addVariable(mean);
                mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
                scale = null;
                rate = null;
                break;
            case OneParameter:
                scale = null;
                rate = null;
                mean = null;
                break;
            default:
                throw new IllegalArgumentException("Unknown parameterization type");
        }
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x < offset) return 0.0;
        return GammaDistribution.pdf(x - offset, getShape(), getScale());
    }

    public double logPdf(double x) {
        if (x < offset) return Double.NEGATIVE_INFINITY;
        return GammaDistribution.logPdf(x - offset, getShape(), getScale());
    }

    public double cdf(double x) {
        if (x < offset) return 0.0;
        return GammaDistribution.cdf(x - offset, getShape(), getScale());
    }

    public double quantile(double y) {
        try {
            return (new GammaDistributionImpl(getShape(), getScale())).inverseCumulativeProbability(y) + offset;
        } catch (MathException e) {
            return Double.NaN;
        }
    }

    public double mean() {
        return GammaDistribution.mean(getShape(), getScale()) + offset;
    }

    public double variance() {
        return GammaDistribution.variance(getShape(), getScale());
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

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    public double getShape() {
        return shape.getValue(0);
    }

    public double getScale() {
        switch (parameterization) {
            case ShapeScale:
                return scale.getValue(0);
            case ShapeRate:
                return (1.0 / rate.getValue(0));
            case ShapeMean:
                return (mean.getValue(0) / getShape());
            case OneParameter:
                return (1.0 / getShape());
            default:
                throw new IllegalArgumentException("Unknown parameterization type");
        }
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
        throw new UnsupportedOperationException("Not implemented");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final GammaParameterizationType parameterization;

    private final Variable<Double> shape;
    private final Variable<Double> scale;
    private final Variable<Double> rate;
    private final Variable<Double> mean;

    private final double offset;

}

