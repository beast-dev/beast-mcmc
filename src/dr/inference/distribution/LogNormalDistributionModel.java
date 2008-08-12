/*
 * LogNormalDistributionModel.java
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
 * A class that acts as a model for log-normally distributed data.
 *
 * @author Alexei Drummond
 * @version $Id: LogNormalDistributionModel.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */

public class LogNormalDistributionModel extends AbstractModel implements ParametricDistributionModel {


    public static final String NORMAL_DISTRIBUTION_MODEL = "logNormalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String PRECISION = "prec";
    public static final String OFFSET = "offset";
    public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";

    boolean isMeanInRealSpace = true;
    boolean usesStDev = true;

    /**
     * Constructor.
     */
    public LogNormalDistributionModel(Parameter meanParameter, Parameter stdevParameter, double offset, boolean meanInRealSpace) {

        super(NORMAL_DISTRIBUTION_MODEL);

        isMeanInRealSpace = meanInRealSpace;

        this.meanParameter = meanParameter;
        this.scaleParameter = stdevParameter;
        this.offset = offset;
        addParameter(meanParameter);
        if (isMeanInRealSpace) {
            meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        } else {
            meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        }
        addParameter(stdevParameter);
        stdevParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public LogNormalDistributionModel(Parameter meanParameter, Parameter scaleParameter,
                                      double offset, boolean meanInRealSpace, boolean usesStDev) {

        super(NORMAL_DISTRIBUTION_MODEL);

        isMeanInRealSpace = meanInRealSpace;
        this.usesStDev = usesStDev;

        this.meanParameter = meanParameter;
        this.scaleParameter = scaleParameter;
        this.offset = offset;
        addParameter(meanParameter);
        if (isMeanInRealSpace) {
            meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        } else {
            meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        }
        addParameter(this.scaleParameter);
        this.scaleParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }


    public final double getS() {
        return scaleParameter.getParameterValue(0);
    }

    public final void setS(double S) {
        scaleParameter.setParameterValue(0, S);
    }

    public final Parameter getSParameter() {
        return scaleParameter;
    }


    private double getStDev() {
        if (usesStDev)
            return getS();
        else
            return Math.sqrt(1.0 / getS());
    }

    /**
     * @return the mean (always in log space)
     */
    public final double getM() {
        if (isMeanInRealSpace) {
            double stDev = getStDev();
            double M = Math.log(meanParameter.getParameterValue(0)) - (0.5 * stDev * stDev);
            return M;
        } else {
            return meanParameter.getParameterValue(0);

        }
    }

    public final void setM(double M) {
        if (isMeanInRealSpace) {
            double stDev = getStDev();
            meanParameter.setParameterValue(0, Math.exp(M + (0.5 * stDev * stDev)));
        } else {
            meanParameter.setParameterValue(0, M);
        }
    }

    public final Parameter getMnParameter() {
        return meanParameter;
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return NormalDistribution.pdf(Math.log(x - offset), getM(), getStDev()) / (x - offset);
    }

    public double logPdf(double x) {
        if (x - offset <= 0.0) return Double.NEGATIVE_INFINITY;
        return NormalDistribution.logPdf(Math.log(x - offset), getM(), getStDev()) - Math.log(x - offset);
    }

    public double cdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return NormalDistribution.cdf(Math.log(x - offset), getM(), getStDev());
    }

    public double quantile(double y) {
        return Math.exp(NormalDistribution.quantile(y, getM(), getStDev())) + offset;
    }

    /**
     * @return the mean of the distribution
     */
    public double mean() {
        return Math.exp(getM() + (getStDev() * getStDev() / 2)) + offset;
    }

    /**
     * @return the variance of the log normal distribution.
     */
    public double variance() {
        double S2 = getStDev();
        S2 *= S2;

        return Math.exp(S2 + 2 * getM()) * (Math.exp(S2) - 1);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(Math.log(x));
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
            double offset = xo.getAttribute(OFFSET, 0.0);

            final boolean meanInRealSpace = xo.getBooleanAttribute(MEAN_IN_REAL_SPACE);

            XMLObject cxo = (XMLObject) xo.getChild(MEAN);
            if (cxo.getChild(0) instanceof Parameter) {
                meanParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                meanParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            cxo = (XMLObject) xo.getChild(PRECISION);
            if (cxo != null) {
                if (cxo.getChild(0) instanceof Parameter) {
                    precParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    precParam = new Parameter.Default(cxo.getDoubleChild(0));
                }

                return new LogNormalDistributionModel(meanParam, precParam, offset, false);
            } else {

                cxo = (XMLObject) xo.getChild(STDEV);
                if (cxo.getChild(0) instanceof Parameter) {
                    stdevParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
                }

                return new LogNormalDistributionModel(meanParam, stdevParam, offset, meanInRealSpace);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(MEAN_IN_REAL_SPACE),
                AttributeRule.newDoubleRule(OFFSET, true),
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
                        new ElementRule(PRECISION,
                                new XMLSyntaxRule[]{
                                        new XORRule(
                                                new ElementRule(Parameter.class),
                                                new ElementRule(Double.class)
                                        )}

                        ))
        };

        public String getParserDescription() {
            return "Describes a normal distribution with a given mean and standard deviation " +
                    "that can be used in a distributionLikelihood element";
        }

        public Class getReturnType() {
            return LogNormalDistributionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Parameter meanParameter;
    private Parameter scaleParameter;
    //	private Parameter precParameter;
    private final double offset;

}
