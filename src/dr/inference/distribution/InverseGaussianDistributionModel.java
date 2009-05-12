package dr.inference.distribution;

import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;
import dr.inference.model.Model;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.InverseGaussianDistribution;
import dr.math.UnivariateFunction;
import dr.xml.*;
import org.w3c.dom.Element;
import org.w3c.dom.Document;

/*
 * InverseGaussianDistributionModel.java
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



/**
 *
 *
 * @author Wai Lok Sibon Li
 * @author Alexei Drummond
 * @version $Id: InverseGaussianDistributionModel.java,v 1.8 2009/03/30 20:25:59 rambaut Exp $
 */

public class InverseGaussianDistributionModel extends AbstractModel implements ParametricDistributionModel {


    public static final String INVERSEGAUSSIAN_DISTRIBUTION_MODEL = "InverseGaussianDistributionModel";
    public static final String MEAN = "mean";
    public static final String SHAPE = "shape";
    public static final String OFFSET = "offset";


    /**
     * Constructor.
     */
    public InverseGaussianDistributionModel(Parameter meanParameter, Parameter shapeParameter, double offset) {

        super(INVERSEGAUSSIAN_DISTRIBUTION_MODEL);

        this.meanParameter = meanParameter;
        this.shapeParameter = shapeParameter;
        this.offset = offset;
        addParameter(meanParameter);
        addParameter(shapeParameter);
        this.meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public final double getShape() {
        return shapeParameter.getParameterValue(0);
    }

    public final void setShape(double S) {
        shapeParameter.setParameterValue(0, S);
    }

    public final Parameter getShapeParameter() {
        return shapeParameter;
    }

    /* Unused method */
    //private double getStDev() {
        //return Math.sqrt(InverseGaussianDistribution.variance(getM(), getShape()));//Math.sqrt((getM()*getM()*getM())/getShape());
    //}

    /**
     * @return the mean
     */
    public final double getM() {
        return meanParameter.getParameterValue(0);
    }

    public final void setM(double M) {
        meanParameter.setParameterValue(0, M);
    }

    public final Parameter getMParameter() {
        return meanParameter;
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return InverseGaussianDistribution.pdf(x - offset, getM(), getShape());
    }

    public double logPdf(double x) {
        if (x - offset <= 0.0) return Double.NEGATIVE_INFINITY;
        return InverseGaussianDistribution.logPdf(x - offset, getM(), getShape());
    }

    public double cdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return InverseGaussianDistribution.cdf(x - offset, getM(), getShape());
    }

    public double quantile(double y) {
        return InverseGaussianDistribution.quantile(y, getM(), getShape()) + offset;
    }

    /**
     * @return the mean of the distribution
     */
    public double mean() {
        return InverseGaussianDistribution.mean(getM(), getShape()) + offset;
    }

    /**
     * @return the variance of the distribution.
     */
    public double variance() {
        return InverseGaussianDistribution.variance(getM(), getShape());
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            System.out.println("just checking if this ever gets used anyways... probably have to change the getLowerBound in LogNormalDistributionModel if it does");
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
            //return Double.NEGATIVE_INFINITY;
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

    public void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
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
     * Reads an inverse gaussian distribution model from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INVERSEGAUSSIAN_DISTRIBUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter meanParam;
            Parameter shapeParam;
            double offset = xo.getAttribute(OFFSET, 0.0);

            XMLObject cxo = (XMLObject) xo.getChild(MEAN);
            if (cxo.getChild(0) instanceof Parameter) {
                meanParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                meanParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            cxo = (XMLObject) xo.getChild(SHAPE);
            if (cxo.getChild(0) instanceof Parameter) {
                shapeParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                shapeParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new InverseGaussianDistributionModel(meanParam, shapeParam, offset);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(OFFSET, true),
                new ElementRule(MEAN,
                        new XMLSyntaxRule[]{
                                new XORRule(
                                        new ElementRule(Parameter.class),
                                        new ElementRule(Double.class)
                                )}
                ),
                new ElementRule(SHAPE,
                        new XMLSyntaxRule[]{
                                new XORRule(
                                        new ElementRule(Parameter.class),
                                        new ElementRule(Double.class)
                                )}
                )
        };

        public String getParserDescription() {
            return "Describes a inverse gaussian distribution with a given mean and shape " +
                    "that can be used in a distributionLikelihood element";
        }

        public Class getReturnType() {
            return InverseGaussianDistributionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Parameter meanParameter;
    private final Parameter shapeParameter;
    private final double offset;

}