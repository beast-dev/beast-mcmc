/*
 * NormalDistributionModel.java
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
import dr.math.distributions.NormalDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for normally distributed data.
 *
 * @author Alexei Drummond
 * @version $Id: NormalDistributionModel.java,v 1.6 2005/05/24 20:25:59 rambaut Exp $
 */

public class NormalDistributionModel extends AbstractModel implements ParametricDistributionModel {


    public static final String NORMAL_DISTRIBUTION_MODEL = "normalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String PREC = "precision";

    /**
     * Constructor.
     */
    public NormalDistributionModel(Parameter meanParameter, Parameter stdevParameter) {

        super(NORMAL_DISTRIBUTION_MODEL);

        this.meanParameter = meanParameter;
        this.stdevParameter = stdevParameter;
        addParameter(meanParameter);
        meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addParameter(stdevParameter);
        stdevParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public NormalDistributionModel(Parameter meanParameter, Parameter scale, boolean isPrecision) {
        super(NORMAL_DISTRIBUTION_MODEL);
        this.hasPrecision = isPrecision;
        this.meanParameter = meanParameter;
        addParameter(meanParameter);
        meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        if (isPrecision) {
            this.precisionParameter = scale;
        } else {
            this.stdevParameter = scale;
        }
        addParameter(scale);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getStdev() {
        if (hasPrecision)
            return 1.0 / Math.sqrt(precisionParameter.getParameterValue(0));
        return stdevParameter.getParameterValue(0);
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return NormalDistribution.pdf(x, mean(), getStdev());
    }

    public double logPdf(double x) {
        return NormalDistribution.logPdf(x, mean(), getStdev());
    }

    public double cdf(double x) {
        return NormalDistribution.cdf(x, mean(), getStdev());
    }

    public double quantile(double y) {
        return NormalDistribution.quantile(y, mean(), getStdev());
    }

    public double mean() {
        return meanParameter.getParameterValue(0);
    }

    public double variance() {
        return NormalDistribution.variance(mean(), getStdev());
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
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

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
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

    /**
     * Reads a normal distribution model from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NORMAL_DISTRIBUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter meanParam;
            Parameter stdevParam;
            Parameter precParam;

            XMLObject cxo = (XMLObject) xo.getChild(MEAN);
            if (cxo.getChild(0) instanceof Parameter) {
                meanParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                meanParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            if (xo.getChild(STDEV) != null) {

                cxo = (XMLObject) xo.getChild(STDEV);
                if (cxo.getChild(0) instanceof Parameter) {
                    stdevParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
                }

                return new NormalDistributionModel(meanParam, stdevParam);
            }

            cxo = (XMLObject) xo.getChild(PREC);
            if (cxo.getChild(0) instanceof Parameter) {
                precParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                precParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
            return new NormalDistributionModel(meanParam, precParam, true);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MEAN,
                        new XMLSyntaxRule[]{
                                new XORRule(
                                        new ElementRule(Parameter.class),
                                        new ElementRule(Double.class)
                                )}
                ),
                new XORRule(
                        new ElementRule(STDEV,
                                new XMLSyntaxRule[]{
                                        new XORRule(
                                                new ElementRule(Parameter.class),
                                                new ElementRule(Double.class)
                                        )}
                        ),
                        new ElementRule(PREC,
                                new XMLSyntaxRule[]{
                                        new XORRule(
                                                new ElementRule(Parameter.class),
                                                new ElementRule(Double.class)
                                        )}
                        )
                )
        };

        public String getParserDescription() {
            return "Describes a normal distribution with a given mean and standard deviation " +
                    "that can be used in a distributionLikelihood element";
        }

        public Class getReturnType() {
            return NormalDistributionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Parameter meanParameter;
    private Parameter stdevParameter;
    private Parameter precisionParameter;
    private boolean hasPrecision = false;

}
