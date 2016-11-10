/*
 * TDistributionModel.java
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
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.RandomGenerator;
import dr.math.distributions.TDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc A. Suchard
 * @author Robert E. Weiss
 */

public class TDistributionModel extends AbstractModel implements ParametricDistributionModel, RandomGenerator {


    public TDistributionModel(Parameter location, Parameter scale, Parameter df) {
        super(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.location = location;
        addVariable(location);
        location.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        this.scale = scale;
        addVariable(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.df = df;
        addVariable(df);
        df.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getLocation() {
        return location.getValue(0);
    }

    public double getScale() {
        return scale.getValue(0);
    }

    public double getDF() {
        return df.getValue(0);
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double logPdf(double x) {
        return TDistribution.logPDF(x, getLocation(), getScale(), getDF());
    }

    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double mean() {
        throw new RuntimeException("Not yet implemented.");
    }

    public double variance() {
        throw new RuntimeException("Not yet implemented.");
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

    private final Variable<Double> location;
    private final Variable<Double> scale;
    private final Variable<Double> df;

    public Object nextRandom() {
        throw new RuntimeException("Not yet implemented.");
    }

    public double logPdf(Object x) {
        double v = (Double) x;
        return logPdf(v);
    }
}
