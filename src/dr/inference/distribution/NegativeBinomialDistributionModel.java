/*
 * NegativeBinomialDistributionModel.java
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
import dr.inferencexml.distribution.NegativeBinomialDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.NegativeBinomialDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for Negative Binomially distributed data.
 *
 * @author Andrew Rambaut
 */

public class NegativeBinomialDistributionModel extends AbstractModel implements ParametricDistributionModel {

    /**
     * Constructor.
     */
    public NegativeBinomialDistributionModel(Variable<Double> mean, Variable<Double> alpha) {

        super(NegativeBinomialDistributionModelParser.NEGATIVE_BINOMIAL_DISTRIBUTION_MODEL);

        this.mean = mean;
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.alpha = alpha;
        addVariable(alpha);
        alpha.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return NegativeBinomialDistribution.pdf(x, mean(), alpha());
    }

    public double logPdf(double x) {
        return NegativeBinomialDistribution.logPdf(x, mean(), alpha());
    }

    public double cdf(double x) {
            return NegativeBinomialDistribution.cdf(x, mean(), alpha());
    }

    public double quantile(double y) {
        throw new RuntimeException("Not implemented.");
    }

    public double mean() {
        return mean.getValue(0);
    }

    public double alpha() {
        return alpha.getValue(0);
    }

    public double variance() {
        throw new RuntimeException("Not implemented!");
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
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

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
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
    private final Variable<Double> alpha;

}
