/*
 * PoissonDistributionModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inferencexml.distribution.PoissonDistributionModelParser;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.PoissonDistribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for Poisson distributed data.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */

public class PoissonDistributionModel extends AbstractModel implements ParametricDistributionModel {
    private final PoissonDistributionImpl distribution;

    /**
     * Constructor.
     */
    public PoissonDistributionModel(Variable<Double> mean) {

        super(PoissonDistributionModelParser.POISSON_DISTRIBUTION_MODEL);

        distribution = new PoissonDistributionImpl(mean.getValue(0));

        this.mean = mean;
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return distribution.probability(x);
    }

    public double logPdf(double x) {
        return Math.log(distribution.probability(x));
    }

    public double cdf(double x) {
        try {
            return distribution.cumulativeProbability(x);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double quantile(double y) {
        try {
            return distribution.inverseCumulativeProbability(y);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double mean() {
        return mean.getValue(0);
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

    @Override
    public Variable<Double> getScaleVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // using a depricated method or else we would be reallocating this every call...
        distribution.setMean(mean());
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

}
