/*
 * LogNormalDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class LogNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String LOGNORMAL_DISTRIBUTION_MODEL = "logNormalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String MU = "mu";
    public static final String SIGMA = "sigma";
    public static final String PRECISION = "precision";
    public static final String OFFSET = "offset";
    public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";
    public static final String STDEV_IN_REAL_SPACE = "stdevInRealSpace";

    public String getParserName() {
        return LOGNORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double offset = xo.getAttribute(OFFSET, 0.0);

        if (xo.hasAttribute(MEAN_IN_REAL_SPACE) || xo.hasAttribute(STDEV_IN_REAL_SPACE)) {
            // If either of these set then use the old logic for parsing. This means that mean
            // means mu if meanInRealSpace is false and stdev means sigma if stdevInRealSpace
            // is false

            if (!xo.hasChildNamed(MEAN) || !xo.hasChildNamed(STDEV)) {
                throw new XMLParseException("If meanInRealSpace attribute is given then the lognormal model must be parameterized with mean and stdev.");
            }

            final boolean meanInRealSpace = xo.getAttribute(MEAN_IN_REAL_SPACE, false);
            if (xo.hasAttribute(STDEV_IN_REAL_SPACE)) {
                boolean stdevInRealSpace = xo.getAttribute(STDEV_IN_REAL_SPACE, false);
                if (!meanInRealSpace && stdevInRealSpace) {
                    throw new XMLParseException("Cannot parameterise Lognormal model with mu and stdev");
                }
                if (meanInRealSpace && !stdevInRealSpace) {
                    throw new XMLParseException("Cannot parameterise Lognormal model with mean and sigma");
                }
            }

            XMLObject cxo = xo.getChild(MEAN);
            Parameter meanParam;
            if (cxo.getChild(0) instanceof Parameter) {
                meanParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                meanParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            cxo = xo.getChild(STDEV);
            Parameter stdevParam;
            if (cxo.getChild(0) instanceof Parameter) {
                stdevParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new LogNormalDistributionModel(meanParam, stdevParam, offset, meanInRealSpace);
        } else {
            // otherwise we decide the parameterization by which parameters are specified.

            if (xo.hasChildNamed(MEAN) && xo.hasChildNamed(STDEV)) {
                XMLObject cxo = xo.getChild(MEAN);
                Parameter meanParam;
                if (cxo.getChild(0) instanceof Parameter) {
                    meanParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    meanParam = new Parameter.Default(cxo.getDoubleChild(0));
                }

                cxo = xo.getChild(STDEV);
                Parameter stdevParam;
                if (cxo.getChild(0) instanceof Parameter) {
                    stdevParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
                }
                return new LogNormalDistributionModel(LogNormalDistributionModel.Parameterization.MEAN_STDEV, meanParam, stdevParam, offset);
            } else if (xo.hasChildNamed(MU)) {
                XMLObject cxo = xo.getChild(MU);
                Parameter muParam;
                if (cxo.getChild(0) instanceof Parameter) {
                    muParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    muParam = new Parameter.Default(cxo.getDoubleChild(0));
                }

                if (xo.hasChildNamed(SIGMA)) {
                    cxo = xo.getChild(SIGMA);
                    Parameter sigmaParam;
                    if (cxo.getChild(0) instanceof Parameter) {
                        sigmaParam = (Parameter) cxo.getChild(Parameter.class);
                    } else {
                        sigmaParam = new Parameter.Default(cxo.getDoubleChild(0));
                    }
                    return new LogNormalDistributionModel(LogNormalDistributionModel.Parameterization.MU_SIGMA, muParam, sigmaParam, offset);
                } else if (xo.hasChildNamed(PRECISION)) {
                    cxo = xo.getChild(PRECISION);
                    Parameter precParam;
                    if (cxo.getChild(0) instanceof Parameter) {
                        precParam = (Parameter) cxo.getChild(Parameter.class);
                    } else {
                        precParam = new Parameter.Default(cxo.getDoubleChild(0));
                    }
                    return new LogNormalDistributionModel(LogNormalDistributionModel.Parameterization.MU_PRECISION, muParam, precParam, offset);

                } else {
                    throw new XMLParseException("Lognormal model must be parameterized as [mean, stdev], [mu, sigma], or [mu, precision]");
                }
            } else {
                throw new XMLParseException("Lognormal model must be parameterized as [mean, stdev], [mu, sigma], or [mu, precision]");
            }

        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(MEAN_IN_REAL_SPACE, true),
            AttributeRule.newBooleanRule(STDEV_IN_REAL_SPACE, true),
            AttributeRule.newDoubleRule(OFFSET, true),
            new XORRule(new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
                    new ElementRule(MU,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                    ))
            ,
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
}
