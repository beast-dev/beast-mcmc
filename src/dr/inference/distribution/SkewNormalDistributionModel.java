/*
 * SkewNormalDistributionModel.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.UnivariateFunction;
import dr.math.distributions.NormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for skew-normally distributed data.
 *
 * @author Marc A. Suchard
 */

public class SkewNormalDistributionModel extends AbstractModel implements ParametricDistributionModel {

    public static final String SKEW_NORMAL_DISTRIBUTION_MODEL = "skewNormalDistributionModel";

    /**
     * Constructor.
     */
    public SkewNormalDistributionModel(Variable<Double> location, Variable<Double> scale, Variable<Double> shape) {

        super(SKEW_NORMAL_DISTRIBUTION_MODEL);

        this.location = location;
        this.scale = scale;
        this.shape = shape;

        addVariable(location);
        location.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        addVariable(shape);
        shape.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

    }

    private double shift(double x) {
        return (x - location.getValue(0)) / scale.getValue(0);
    }


    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        double delta = shift(x);
        return 2.0 / scale.getValue(0) * NormalDistribution.pdf(delta, 0, 1) * NormalDistribution.cdf(shape.getValue(0) * delta, 0, 1);
    }

    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    public double cdf(double x) {
        throw new IllegalArgumentException("Not yet implement");
//        return NormalDistribution.cdf(x, mean(), getStdev());
    }

    public double quantile(double y) {
        throw new IllegalArgumentException("Not yet implement");
    }

    public double mean() {
        throw new IllegalArgumentException("Not yet implement");
    }

    public double variance() {
        throw new IllegalArgumentException("Not yet implement");
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

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return location;
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Variable<Double> location;
    private final Variable<Double> scale;
    private final Variable<Double> shape;
}
