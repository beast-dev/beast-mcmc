/*
 * NormalDistributionModelParser.java
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

import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class NormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "normalDistributionModel";
    public static final String MEAN = "mean";
    public static final String MU = "mu";
    public static final String SCALE = "scale";
    public static final String STDEV = "stdev";
    public static final String PREC = "precision";

    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter meanParam;
        Parameter scaleParam;
        Parameter precParam;

        XMLObject cxo = xo.getChild(MEAN);
        if (cxo.getChild(0) instanceof Parameter) {
            meanParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            meanParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(STDEV);
        if (cxo == null) {
            cxo = xo.getChild(SCALE);
        }

        if (cxo != null) {
            if (cxo.getChild(0) instanceof Parameter) {
                scaleParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                scaleParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new NormalDistributionModel(meanParam, scaleParam);
        } else {
            cxo = xo.getChild(PREC);
            if (cxo.getChild(0) instanceof Parameter) {
                precParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                precParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
        }
        return new NormalDistributionModel(meanParam, precParam, true);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new XORRule(
                    new XORRule(
                            new ElementRule(SCALE,
                                    new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Parameter.class),
                                                    new ElementRule(Double.class)
                                            )}
                            ),
                            new ElementRule(STDEV,
                                    new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Parameter.class),
                                                    new ElementRule(Double.class)
                                            )}
                            )
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

}
