/*
 * TwoPieceLocationScaleDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.TwoPieceLocationScaleDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class TwoPieceLocationScaleDistributionModelParser extends AbstractXMLObjectParser {

    public static final String DISTRIBUTION_MODEL = "twoPieceLocationScaleDistributionModel";
    public static final String LOCATION = "location";
    public static final String SIGMA = "sigma";
    public static final String GAMMA = "gamma";
    public static final String PARAMETERIZATION = "parameterization";

    public String getParserName() {
        return DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter locationParam;
        Parameter sigmaParam;
        Parameter gammaParam;

        XMLObject cxo = xo.getChild(LOCATION);
        if (cxo.getChild(0) instanceof Parameter) {
            locationParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            locationParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        String parameterizationLabel = (String) xo.getAttribute(PARAMETERIZATION);
        TwoPieceLocationScaleDistributionModel.Parameterization parameterization =
                TwoPieceLocationScaleDistributionModel.Parameterization.parseFromString(parameterizationLabel);
        if (parameterization == null) {
            throw new XMLParseException("Unrecognized parameterization '" + parameterizationLabel + "'");
        }

        cxo = xo.getChild(SIGMA);
        if (cxo.getChild(0) instanceof Parameter) {
            sigmaParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            sigmaParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(GAMMA);
        if (cxo.getChild(0) instanceof Parameter) {
            gammaParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            gammaParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        Distribution distribution = (Distribution) xo.getChild(Distribution.class);

        return new TwoPieceLocationScaleDistributionModel(locationParam, distribution,
                sigmaParam, gammaParam, parameterization);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringArrayRule(PARAMETERIZATION),
            new ElementRule(LOCATION,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new ElementRule(SIGMA,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}, true
            ),
            new ElementRule(GAMMA,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}, true
            ),
            new ElementRule(Distribution.class),
    };

    public String getParserDescription() {
        return "Describes a two-piece location-scale distribution with a given location and two scales " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return TwoPieceLocationScaleDistributionModelParser.class;
    }

}
