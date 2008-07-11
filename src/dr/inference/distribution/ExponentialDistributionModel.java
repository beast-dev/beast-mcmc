/*
 * ExponentialDistributionModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.math.ExponentialDistribution;
import dr.math.UnivariateFunction;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for exponentially distributed data.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ExponentialDistributionModel.java,v 1.12 2005/05/24 20:25:59 rambaut Exp $
 */

public class ExponentialDistributionModel extends AbstractModel implements ParametricDistributionModel {

    public static final String EXPONENTIAL_DISTRIBUTION_MODEL = "exponentialDistributionModel";
    public static final String MEAN = "mean";
    public static final String OFFSET = "offset";

    /**
     * Constructor.
     */
    public ExponentialDistributionModel(Parameter meanParameter) {

        this(meanParameter, 0.0);
    }


    /**
     * Constructor.
     */
    public ExponentialDistributionModel(Parameter meanParameter, double offset) {

        super(EXPONENTIAL_DISTRIBUTION_MODEL);

        this.meanParameter = meanParameter;
        this.offset = offset;

        addParameter(meanParameter);
        meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
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

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
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

    public void handleParameterChangedEvent(Parameter parameter, int index) {
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

    /**
     * Reads an exponential distribution model from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EXPONENTIAL_DISTRIBUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject object = (XMLObject) xo.getChild(MEAN);

            double offset = xo.getAttribute(OFFSET, 0.0);

            Parameter meanParameter;
            try {
                double mean = object.getDoubleChild(0);
                meanParameter = new Parameter.Default(mean);
            } catch (Exception pe) {
                meanParameter = (Parameter) object.getChild(0);
            }
            return new ExponentialDistributionModel(meanParameter, offset);
        }

        public String getParserDescription() {
            return "A model of an exponential distribution.";
        }

        public Class getReturnType() {
            return ExponentialDistributionModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(OFFSET, true),
                new XORRule(
                        new ElementRule(MEAN, Double.class),
                        new ElementRule(MEAN, Parameter.class)
                )
        };
    };

    // **************************************************************
    // Private methods
    // **************************************************************

    private double getMean() {
        return meanParameter.getParameterValue(0);
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Parameter meanParameter = null;
    private double offset = 0.0;

}

