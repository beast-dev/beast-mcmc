/*
 * GammaDistributionModel.java
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
import dr.math.UnivariateFunction;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for gamma distributed data.
 *
 * @author Alexei Drummond
 * @version $Id: GammaDistributionModel.java,v 1.6 2005/05/24 20:25:59 rambaut Exp $
 */

public class GammaDistributionModel extends AbstractModel implements ParametricDistributionModel {

    public static final String GAMMA_DISTRIBUTION_MODEL = "gammaDistributionModel";
    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";


    /**
     * Construct a constant mutation rate model.
     */
    public GammaDistributionModel(Parameter shapeParameter, Parameter scaleParameter) {

        super(GAMMA_DISTRIBUTION_MODEL);

        this.shapeParameter = shapeParameter;
        this.scaleParameter = scaleParameter;
        addParameter(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        addParameter(scaleParameter);
        scaleParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    /**
     * Construct a constant mutation rate model.
     */
    public GammaDistributionModel(Parameter shapeParameter) {

        super(GAMMA_DISTRIBUTION_MODEL);

        this.shapeParameter = shapeParameter;
        addParameter(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return GammaDistribution.pdf(x, getShape(), getScale());
    }

    public double logPdf(double x) {
        return GammaDistribution.logPdf(x, getShape(), getScale());
    }

    public double cdf(double x) {
        return GammaDistribution.cdf(x, getShape(), getScale());
    }

    public double quantile(double y) {
        return GammaDistribution.quantile(y, getShape(), getScale());
    }

    public double mean() {
        return GammaDistribution.mean(getShape(), getScale());
    }

    public double variance() {
        return GammaDistribution.variance(getShape(), getScale());
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
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
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
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
     * Reads a gamma distribution model from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GAMMA_DISTRIBUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter shapeParam;
            Parameter scaleParam;

            XMLObject cxo = (XMLObject) xo.getChild(SHAPE);
            if (cxo.getChild(0) instanceof Parameter) {
                shapeParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                shapeParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            cxo = (XMLObject) xo.getChild(SCALE);
            if (cxo.getChild(0) instanceof Parameter) {
                scaleParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                scaleParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new GammaDistributionModel(shapeParam, scaleParam);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(SHAPE,
                        new XMLSyntaxRule[]{
                                new XORRule(
                                        new ElementRule(Parameter.class),
                                        new ElementRule(Double.class)
                                )}
                ),
                new ElementRule(SCALE,
                        new XMLSyntaxRule[]{
                                new XORRule(
                                        new ElementRule(Parameter.class),
                                        new ElementRule(Double.class)
                                )}
                )
        };

        public String getParserDescription() {
            return "Describes a gamma distribution with a given shape and scale " +
                    "that can be used in a distributionLikelihood element.";
        }

        public Class getReturnType() {
            return GammaDistributionModel.class;
        }
    };

    // **************************************************************
    // Private methods
    // **************************************************************

    private double getShape() {
        return shapeParameter.getParameterValue(0);
    }

    private double getScale() {
        if (scaleParameter == null) return (1.0 / getShape());
        return scaleParameter.getParameterValue(0);
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Parameter shapeParameter = null;
    private Parameter scaleParameter = null;

}

